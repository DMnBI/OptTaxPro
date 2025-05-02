#!/usr/bin/env python

import os
import sys
import argparse
import logging
import glob
import tempfile

import pandas as pd
import subprocess as sp
import multiprocessing as mp
from tqdm import tqdm
from Bio import SeqIO

from commons import (
	VSEARCH_BASE,
	auto_run, 
	run_cmd
)

def parse_args(argv = sys.argv[1:]):
	parser = argparse.ArgumentParser()
	# optional arguments
	parser.add_argument("-q", "--quiet",
		dest="quiet",
		action="store_true",
		help="quiet; do not print any message on the screen")
	parser.add_argument("--log",
		dest="log",
		metavar="FILE",
		help="log file name; default stderr")
	parser.add_argument("--t_dir",
		dest="t_dir",
		metavar="PATH",
		default="./",
		help="temporal directory to saving intermediate files; default ./")
	parser.add_argument("-p", "--processes",
		dest="p",
		metavar="INT",
		type=int,
		default=10,
		help="the number of processes to run concurrently; default 10")
	parser.add_argument("--suffix",
		dest="suffix",
		metavar="STR",
		default='.fasta',
		help="suffix of input files; default .fasta")

	# VSEARCH arguments
	vse_group = parser.add_argument_group("VSEARCH arguments")
	vse_group.add_argument("-t", "--threads",
		dest='t',
		metavar="INT",
		type=int,
		default=5,
		help="the number of threads for each process; default 5")
	vse_group.add_argument("--cutoffs",
		dest='cutoffs',
		metavar="FLOAT",
		nargs='+',
		type=float,
		default=[0.94, 0.86, 0.75],
		help="homology search cutoffs for iterative search; default [0.94, 0.86, 0.75]")
	vse_group.add_argument("--maxaccepts",
		dest='maxaccepts',
		metavar="INT",
		type=int,
		default=5000,
		help="the number of maximum accepted hits; default 5000")

	# required arguments
	req_group = parser.add_argument_group("required arguments")
	req_group.add_argument("-i", "--i_dir",
		dest="i_dir",
		metavar="PATH",
		required=True,
		help="directory involving raw fastq files; required")
	req_group.add_argument("-o", "--o_dir",
		dest="o_dir",
		metavar="PATH",
		required=True,
		help="directory for saving preprocessed samples; required")
	req_group.add_argument("-d", "--db",
		dest="db",
		metavar="FILE",
		required=True,
		help="16S database file; required")

	return parser.parse_args(argv)

def get_unmapped(args):
	seq_file, seq_type, m6_file, out_file = args

	m6 = pd.read_csv(m6_file, sep='\t', header=None, usecols=[0, 1, 2])
	m6.columns = ['query', 'db', 'identity']
	mapped = {sid:True for sid in m6['query'].unique()}

	seq_index = SeqIO.index(seq_file, seq_type)
	unmapped = [sid for sid in seq_index if sid not in mapped]
	if len(unmapped) <= 0:
		return

	with open(out_file, 'w') as ostream:
		for sid in unmapped:
			print(seq_index[sid].format("fasta"), file=ostream)

def cat(args):
	m6_list, merged = args
	with open(merged, 'w') as ostream:
		_ = run_cmd(['cat'] + m6_list, stdout=ostream)

def main(argv = sys.argv[1:]):
	args = parse_args(argv)
	logging.basicConfig(
		filename=args.log,
		format='[%(asctime)s] %(message)s', 
		datefmt='%Y/%m/%d %I:%M:%S', 
		level=logging.INFO if not args.quiet else 100,
	)

	pname = os.path.basename(sys.argv[0])
	logging.info(f"[{pname}] {sys.executable} {' '.join(sys.argv)}")
	rpath, _ = os.path.split(os.path.realpath(__file__))

	assert os.path.isfile(args.db), "DB file is missing"

	if not os.path.isdir(args.o_dir):
		os.makedirs(args.o_dir, exist_ok=True)
	t_dir = f"{args.t_dir}/intermediate"
	os.makedirs(t_dir, exist_ok=True)

	targets = glob.glob(f"{args.i_dir}/*{args.suffix}")
	suffix = args.suffix
	stype = 'fastq' if args.suffix == ".fastq" else 'fasta'

	# RUN homology search iteratively
	for id_cutoff in args.cutoffs:
		# SEARCH aginast DB
		logging.info(f"[{pname}] RUN VSEARCH (cutoff={id_cutoff:.2f})")
		cmds = []
		m6_list = []
		for seq_file in targets:
			m6_file = f"{t_dir}/{os.path.basename(seq_file).replace(suffix, f'.{id_cutoff*100:.0f}.m6')}"
			m6_list.append(m6_file)

			cmd = VSEARCH_BASE.copy()
			cmd += ['--id', id_cutoff,
				'--maxaccepts', args.maxaccepts,
				'--threads', args.t,
				'--usearch_global', seq_file,
				'--db', args.db, '--blast6out', m6_file
			]
			cmds.append(cmd)
		auto_run(run_cmd, cmds, args.p, quiet=args.quiet, desc=f'search (cutoff={id_cutoff:.2f})')
		logging.info(f"[{pname}] Done. (VSEARCH; cutoff={id_cutoff:.2f})")

		# GET unmapped sequences
		logging.info(f"[{pname}] GET unmapped seqs (cutoff={id_cutoff:.2f})")
		cmds = []
		unmapped = []
		for seq_file, m6_file in zip(targets, m6_list):
			unmapped_file = m6_file.replace(".m6", ".fna")
			unmapped.append(unmapped_file)

			cmd = (seq_file, stype, m6_file, unmapped_file)
			cmds.append(cmd)
		auto_run(get_unmapped, cmds, args.p // args.t, quiet=args.quiet, desc=f"get unmapped (cutoff={id_cutoff:.2f})")
		logging.info(f"[{pname}] Done. (get unmapped; cutoff={id_cutoff:.2f})")

		targets = unmapped.copy()
		stype = 'fasta'
		suffix = f".{id_cutoff*100:.0f}.fna"

	# merge m6 files
	logging.info(f"[{pname}] merge m6 files")
	cmds = []
	first_suffix = f".{args.cutoffs[0]*100:.0f}.m6"
	targets = glob.glob(f"{t_dir}/*{first_suffix}")
	for target in targets:
		m6_list = [target] + [target.replace(first_suffix, f".{id_cutoff*100:.0f}.m6") for id_cutoff in args.cutoffs[1:]]
		merged = target.replace(first_suffix, ".merged.m6").replace(t_dir, args.o_dir)
		cmds.append((m6_list, merged))
	auto_run(cat, cmds, args.p, quiet=args.quiet, desc=f"cat m6 files")

	logging.info(f"[{pname}] Done. (iterative search)")
	return 0

if __name__ == "__main__":
	exit(main())
#!/usr/bin/env python

import os
import sys
import argparse
import logging
import glob
import tempfile

import subprocess as sp
import multiprocessing as mp
from tqdm import tqdm

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
	vse_group.add_argument("-c", "--cutoff",
		dest='cutoff',
		metavar="FLOAT",
		type=float,
		default=0.99,
		help="clustering cutoff for picking OTU; default 0.99")
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

	return parser.parse_args(argv)


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

	if not os.path.isdir(args.o_dir):
		os.makedirs(args.o_dir, exist_ok=True)

	targets = glob.glob(f"{args.i_dir}/*{args.suffix}")

	# RUN VSEARCH cluster_fast
	logging.info(f"[{pname}] RUN VSEARCH (cluster_fast)")
	cmds = []
	for fastq in targets:
		centroid = f"{args.o_dir}/{os.path.basename(fastq).replace(args.suffix, '.rep.fna')}"
		uc = centroid.replace(".rep.fna", ".uc")
		cmd = VSEARCH_BASE.copy()
		cmd += ['--id', args.cutoff,
			'--maxaccepts', args.maxaccepts,
			'--threads', args.t,
			'--cluster_fast', fastq,
			'--centroids', centroid, '--uc', uc
		]
		cmds.append(cmd)
	auto_run(run_cmd, cmds, args.p, quiet=args.quiet, desc='cluster_fast')
	logging.info(f"[{pname}] Done. (cluster_fast)")

	return 0

if __name__ == "__main__":
	exit(main())
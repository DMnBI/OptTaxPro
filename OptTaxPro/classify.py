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

from utils.commons import run_cmd
from utils.taxonomy_assignment import (
	ACC2TAXID, 
	HSG_FILE, 
	RANKS, 
	CUTOFFS
)

DPATH = f"{os.path.split(os.path.realpath(__file__))[0]}/data"
DB = f"{DPATH}/NRDB.fna"

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
		default=".rep.fna",
		help="suffix of input files; default .rep.fna")
	parser.add_argument("--t_dir",
		dest="t_dir",
		metavar="PATH",
		default="./",
		help="temporal directory to saving intermediate files; default ./")

	# VSEARCH arguments
	vse_group = parser.add_argument_group("VSEARCH arguments")
	vse_group.add_argument("-d", "--db",
		dest="db",
		metavar="FILE",
		default=DB,
		help=f"16S database file; default {DB}")
	vse_group.add_argument("-t", "--threads",
		dest='t',
		metavar="INT",
		type=int,
		default=5,
		help="the number of threads for each process; default 5")
	vse_group.add_argument("--search-cutoffs",
		dest='search_cutoffs',
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

	# assignment arguments
	asg_group = parser.add_argument_group("assignment arguments")
	asg_group.add_argument("-r", "--ranks",
		dest='ranks',
		metavar="RANK",
		nargs="+",
		default=RANKS,
		help=f"target ranks to assign; default {' '.join(RANKS)}")
	asg_group.add_argument("--assign-cutoffs",
		dest='assign_cutoffs',
		metavar="FLOAT",
		nargs="+",
		default=CUTOFFS,
		help=f"identity cutoffs for each rank; default {' '.join([str(c) for c in CUTOFFS])}")
	asg_group.add_argument("--acc2taxid",
		dest='acc2taxid',
		metavar="FILE",
		default=ACC2TAXID,
		help=f"accession to taxid mapping file; default {ACC2TAXID}")
	asg_group.add_argument("--hsg",
		dest="hsg",
		metavar="FILE",
		default=HSG_FILE,
		help=f"HSG definition file; default {HSG_FILE}")
	asg_group.add_argument("--remove-self",
		dest="rm_self",
		action="store_true",
		help="remove self hit according to the assembly_accession")
	asg_group.add_argument("--u_dir",
		dest="u_dir",
		metavar="PATH",
		help="directory involving uclust results")
	asg_group.add_argument("--uc-suffix",
		dest="uc_suffix",
		metavar="STR",
		default=".uc",
		help="suffix of uclust result files; default .uc")
	asg_group.add_argument("--add-name",
		dest="add_name",
		action="store_true",
		help="add scientific_name column")

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
	t_dir = tempfile.mkdtemp(dir=args.t_dir, prefix="99_classify_", suffix="_tmp")
	targets = glob.glob(f"{args.i_dir}/*{args.suffix}")

	assert len(args.ranks) == len(args.assign_cutoffs), "Length mismatch between ranks and cutoffs"
	assert os.path.isfile(args.acc2taxid), "acc2taxid file is missing"
	assert os.path.isfile(args.hsg), "hsg file is missing"
	assert os.path.isfile(args.db), "DB file is missing"

	# RUN iterative search
	logging.info(f"[{pname}] RUN iterative search")
	cmd = ["python", f"{rpath}/utils/iterative_search.py",
		'-i', args.i_dir, '-o', t_dir,
		'--db', args.db, 
		'-t', args.t, '-p', args.p,
		'--t_dir', t_dir,
		'--maxaccepts', args.maxaccepts,
		'--suffix', args.suffix
	]
	cmd += ['--cutoffs'] + args.search_cutoffs
	if args.quiet:
		cmd += ['quiet']
	if args.log is not None:
		cmd += ['--log', args.log]
	_ = run_cmd(cmd)

	# RUN taxonomy assignment
	logging.info(f"[{pname}] RUN taxonomy assignment")
	cmd = ["python", f"{rpath}/utils/taxonomy_assignment.py",
		'-i', t_dir, '-o', args.o_dir,
		'-p', args.p,
		'--suffix', ".merged.m6",
		'--acc2taxid', args.acc2taxid,
		'--hsg', args.hsg,
	]
	cmd += ['--cutoffs'] + args.assign_cutoffs
	cmd += ['--ranks'] + args.ranks
	if args.quiet:
		cmd += ['quiet']
	if args.log is not None:
		cmd += ['--log', args.log]
	if args.rm_self:
		cmd += ['--remove-self']
	if args.u_dir is not None:
		cmd += ['--u_dir', args.u_dir]
	if args.uc_suffix is not None:
		cmd += ['--uc-suffix', args.uc_suffix]
	if args.add_name is not None:
		cmd += ['--add-name']
	_ = run_cmd(cmd)

	logging.info(f"[{pname}] Done. (clustering)")
	return 0

if __name__ == "__main__":
	exit(main())
	
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
		default=".fastq",
		help="suffix of input files; default .fastq")

	# VSEARCH arguments
	otu_group = parser.add_argument_group("VSEARCH arguments")
	otu_group.add_argument("--cutoff",
		dest='cutoff',
		metavar="FLOAT",
		type=float,
		default=0.99,
		help="OTU clustering cutoff; default 0.97")
	otu_group.add_argument("--threads",
		dest="threads",
		metavar="INT",
		type=int,
		default=5,
		help="the number of threads for VSEARCH; default 5")

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

	# RUN cluster_fast
	logging.info(f"[{pname}] RUN cluster_fast")
	cmd = ["python", f"{rpath}/utils/cluster_fast.py",
		'-i', args.i_dir, '-o', args.o_dir,
		'-t', args.threads, '-p', args.p,
		'-c', args.cutoff, 
		'--suffix', args.suffix
	]
	if args.quiet:
		cmd += ['quiet']
	if args.log is not None:
		cmd += ['--log', args.log]
	_ = run_cmd(cmd)

	logging.info(f"[{pname}] Done. (clustering)")
	return 0

if __name__ == "__main__":
	exit(main())
	
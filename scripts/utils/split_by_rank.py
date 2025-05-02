#!/usr/bin/env python

import sys
import argparse
import os.path
import logging

import pandas as pd
import numpy as np
from Bio import SeqIO

def parse_args(argv = sys.argv[1:]):
	parser = argparse.ArgumentParser()

	# optional arguments
	parser.add_argument("-q", "--quiet",
		dest="quiet",
		action="store_true",
		help="quiet; do not print any message")
	parser.add_argument("--log",
		dest="log",
		metavar="STR",
		help="log file name; default stderr")
	parser.add_argument("--rank",
		dest='rank',
		choices=('species', 'genus'),
		default='species',
		help="target rank to be split; default species")
	parser.add_argument("--skip-single",
		dest='skip_single',
		action='store_true',
		help="skip exporting singleton file")

	# required arguments
	req_group = parser.add_argument_group("required arguments")
	req_group.add_argument("-r", "--reference",
		dest="ref",
		metavar="FASTA",
		required=True,
		help="reference sequences to be split; required")
	req_group.add_argument("-o", "--out-dir",
		dest="o_dir",
		metavar="PATH",
		required=True,
		help="output file path; required")
	req_group.add_argument("-a", "--acc2taxid",
		dest="acc2taxid",
		metavar="TXT",
		required=True,
		help="accession to taxid mapping file; required")

	return parser.parse_args(argv)

def group_by_rank(sid_list, acc2taxid, rank):
	groups = {}
	for sid in sid_list:
		acc = sid.split("__")[1].split('.')[0]
		taxid = acc2taxid.loc[acc, f"{rank}_taxid"]

		# Handling Escherichia-Shigella
		if taxid == 620:
			taxid = 561

		if taxid not in groups:
			groups[taxid] = []

		groups[taxid].append(sid)

	return groups

def print_summary(groups, rank):
	size_list = np.array([len(groups[species]) for species in groups])

	print(f"Among {size_list.sum()} sequences,")
	print(f"Total {len(groups)} {rank} were detected.")
	print(f"Min: {size_list.min()}, Med: {np.median(size_list)}, Max: {size_list.max()}")

def push_to_hdd(groups, parser, o_dir, skip_single=False):
	skipped = 0
	for taxid in groups:
		if skip_single and len(groups[taxid]) == 1:
			skipped += 1
			continue

		ostream = open(f"{o_dir}/{taxid}.fna", 'w')
		for sid in groups[taxid]:
			print(parser[sid].format('fasta'), file=ostream)
		ostream.close()

	if skip_single:
		print(f"{skipped} singletons were skipped")

def main(argv = sys.argv[1:]):
	args = parse_args(argv)
	logging.basicConfig(
		filename=args.log,
		format='[%(asctime)s] %(message)s', 
		datefmt='%Y/%m/%d %I:%M:%S', 
		level=logging.INFO if not args.quiet else 100,
	)

	logging.info(f"{sys.executable} {' '.join(sys.argv)}")
	pname = os.path.basename(sys.argv[0])

	acc2taxid = pd.read_csv(args.acc2taxid, sep='\t', index_col="accession")
	parser = SeqIO.index(args.ref, 'fasta')

	groups = group_by_rank([sid for sid in parser], acc2taxid, args.rank)
	print_summary(groups, args.rank)
	push_to_hdd(groups, parser, args.o_dir, args.skip_single)

	logging.info(f"Done. ({pname})")
	return 0

# main
if __name__ == "__main__":
	exit(main())

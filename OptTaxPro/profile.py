#!/usr/bin/env python

import os
import sys
import argparse
import logging
import glob
import tempfile

import pandas as pd
import multiprocessing as mp
from tqdm import tqdm

from utils.taxonomy_assignment import load_uc, expand_to_uc

DEFAULT_RANKS = ['species', 'HSG', 'genus', 'family']

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
		default=".ta.csv",
		help="suffix of input files; default .ta.csv")
	parser.add_argument("--output-prefix",
		dest="prefix",
		metavar="STR",
		default="profile",
		help="prefix of output tables; profile")
	parser.add_argument("--t_dir",
		dest="t_dir",
		metavar="PATH",
		default="./",
		help="temporal directory to saving intermediate files; default ./")

	# profile arguments
	prof_group = parser.add_argument_group("profile arguments")
	prof_group.add_argument("--u_dir",
		dest="u_dir",
		metavar="FILE",
		help="directory involving uclust results")
	prof_group.add_argument("--uc-suffix",
		dest="u_suffix",
		metavar="STR",
		default=".uc",
		help="suffix of uclust file; default .uc")
	prof_group.add_argument("--ranks",
		dest="ranks",
		metavar="RANK",
		nargs="+",
		default=DEFAULT_RANKS,
		help=f"target ranks to generate profile; default {' '.join(DEFAULT_RANKS)}")
	prof_group.add_argument("--base-col",
		dest="base_col",
		choices=('assigned', 'scientific_name'),
		default="assigned",
		help="base column to be row of the profile results; default assigned")
	prof_group.add_argument("--filtering-cutoffs",
		dest="filter_cutoffs",
		metavar="FLOAT",
		nargs='*',
		type=float,
		default=None,
		help="filtering cutoff for low abundance")
	prof_group.add_argument("--filtering-pivot",
		dest="filter_pivot",
		choices=("max", "min", "median", "avg"),
		default='max',
		help="method for obtain filtering pivots; default max")
	

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


def get_profile(fargs):
	def calc_for_rank(df, rank, base):
		#####
		# subdf for each rank
		# 
		# [subdf]
		# seqid | rank    | taxid | identity | assigned
		# ----------------------------------------------
		# seq1  | species | 561   | 99       | 561
		# seq1  | species | 562   | 99       | 562
		# seq2  | species | 1281  | 100      | 1281
		# seq3  | species | 561   | 99.5     | 561
		# ...
		#
		subdf = df.query("rank == @rank")

		#####
		# split read count when multiple besthits occur
		#
		# [w]
		# seqid | weight
		# ---------------
		# seq1  | 0.5
		# seq2  | 1
		# seq3  | 1
		# ...
		#
		w = subdf.groupby('seqid')[base].apply(lambda x: 1 / x.shape[0])
		w = pd.DataFrame({'weight': w})

		#####
		# merged subdf
		# 
		# [subdf]
		# seqid | rank    | taxid | identity | assigned | weight
		# -------------------------------------------------------
		# seq1  | species | 561   | 99       | 561      | 0.5
		# seq1  | species | 562   | 99       | 562      | 0.5
		# seq2  | species | 1281  | 100      | 1281     | 1
		# seq3  | species | 561   | 99.5     | 561      | 1
		# ...
		#
		subdf = pd.merge(left=subdf, right=w, left_on='seqid', right_index=True, how='left')

		#####
		# aggregating read counts according to thier weight
		#
		# [counts]
		# assigned | count
		# -----------------
		# 561      | 1.5
		# 562      | 0.5
		# 1281     | 1
		# ...
		#
		# [total] = 3
		counts = subdf.groupby(base)['weight'].sum()
		total = counts.sum()

		#####
		# return relative proportion
		#
		# assigned | profile
		# 561      | 0.5
		# 562      | 0.166667
		# 1281     | 0.333333
		# ...
		return counts / total

	ta_file, args = fargs

	df = pd.read_csv(ta_file, sep=',')
	sample_name = os.path.basename(ta_file)[:-len(args.suffix)]

	if args.u_dir is not None:
		uc_file = f"{args.u_dir}/{sample_name}{args.uc_suffix}"
		df = expand_to_uc(df, args.uc)

	profiled = {}
	for rank in args.ranks:
		res = calc_for_rank(df, rank, args.base_col)
		profiled[rank] = {sample_name:res}

	return profiled

def filtering_table(table, cutoff, method, others_key=-2):
	if cutoff == 0.0:
		return table.copy()

	table['pivot'] = (
		table.max(axis=1) if method == 'max' else \
		table.min(axis=1) if method == 'min' else \
		table.mean(axis=1) if method == 'avg' else \
		table.median(axis=1) 
	)

	others = table.query('pivot < @cutoff')
	table = table.query('pivot >= @cutoff').copy()
	table.loc[others_key] = others.sum(axis=0)

	return table.drop('pivot', axis=1)

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

	logging.info(f"[{pname}] Measure profiles (ranks={' '.join(args.ranks)})")
	fargs = []
	for ta_file in targets:
		fargs.append((ta_file, args))

	rank_df = {rank:{} for rank in args.ranks}
	pool = mp.Pool(args.p)
	for ret in tqdm(pool.imap_unordered(get_profile, fargs), total=len(fargs), desc='profiling'):
		for rank in args.ranks:
			rank_df[rank].update(ret[rank])

	tables = {rank: pd.DataFrame(rank_df[rank]).fillna(0.0) for rank in rank_df}
	if args.filter_cutoffs is not None:
		logging.info(f"[{pname}] filtering low abundances (CUTOFFs={' '.join([f'{c:.2f}' for c in args.filter_cutoffs])})")
		others_key = -2 if args.base_col == 'assigned' else 'others'

		filtered = {cutoff:{} for cutoff in args.filter_cutoffs}
		for cutoff in args.filter_cutoffs:
			for rank in args.ranks:
				filtered[cutoff][rank] = filtering_table(tables[rank], cutoff, args.filter_pivot, others_key)

	logging.info(f"[{pname}] Saving profiled tables")
	for rank in args.ranks:
		base_name = f"{args.o_dir}/{args.prefix}.{rank}.csv"
		if args.filter_cutoffs is not None:
			for cutoff in args.filter_cutoffs:
				out_name = base_name.replace(".csv", f".cut{cutoff*100:.0f}.csv")
				filtered[cutoff][rank].to_csv(out_name, sep=',')
		else:
			tables[rank].to_csv(base_name, sep=',')

	logging.info(f"[{pname}] Done. (clustering)")
	return 0

if __name__ == "__main__":
	exit(main())

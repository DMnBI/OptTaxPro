#!/usr/bin/env python

import os
import sys
import argparse
import logging
import glob
import tempfile

import pandas as pd
import numpy as np
import subprocess as sp
import multiprocessing as mp
from tqdm import tqdm
from Bio import SeqIO
from ete3 import NCBITaxa

if __name__ == "__main__":
	from commons import auto_run

DPATH = os.path.split(os.path.realpath(__file__))[0].replace('/utils', '/data')
ACC2TAXID = f"{DPATH}/acc2taxid.txt"
HSG_FILE = f"{DPATH}/HSG.csv"

RANKS = ['species', 'HSG', 'genus', 'family', 'order', 'class', 'phylum']
CUTOFFS = [97, 97, 94, 86, 82, 78, 75]
ncbi = NCBITaxa()

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
#	parser.add_argument("--t_dir",
#		dest="t_dir",
#		metavar="PATH",
#		default="./",
#		help="temporal directory to saving intermediate files; default ./")
	parser.add_argument("-p", "--processes",
		dest="p",
		metavar="INT",
		type=int,
		default=10,
		help="the number of processes to run concurrently; default 10")
	parser.add_argument("--suffix",
		dest="suffix",
		metavar="STR",
		default='.merged.m6',
		help="suffix of input files; default .merged.m6")

	# assignment arguments
	asg_group = parser.add_argument_group("assignment arguments")
	asg_group.add_argument("-r", "--ranks",
		dest='ranks',
		metavar="RANK",
		nargs="+",
		default=RANKS,
		help=f"target ranks to assign; default {' '.join(RANKS)}")
	asg_group.add_argument("-c", "--cutoffs",
		dest='cutoffs',
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

def load_m6(file_name, rm_self=False):
	m6 = pd.read_csv(file_name, sep='\t', usecols=[0, 1, 2])
	m6.columns = ['seqid', 'db', 'identity']

	if rm_self:
		m6['itself'] = [row['seqid'].split('__')[0] == row['db'].split('__')[0] for _, row in m6.iterrows()]
		m6 = m6.query('~itself').drop('itself', axis=1)

	return m6

def load_uc(file_name):
	uc = pd.read_csv(file_name, sep='\t', header=None)
	uc.columns = ['record_type', 'cluster_num', 'size', 'identity', 'strand', 'col6', 'col7', 'CIGAR', 'query', 'rep']
	uc['rep'] = [row['query'] if row['rep'] == '*' else row['rep'] for _, row in uc.iterrows()]
	uc = uc.query("record_type != 'C'")

	return uc[['query', 'rep']]

#####
# Assignment results
#
# seqid | rank    | taxid | identity | assigned
# ----------------------------------------------
#  seq1 | species | 561   | 96.1     | -1
#  seq1 | species | 562   | 96.1     | -1
#  seq1 | HSG     | 560:0 | 96.1     | -1
#  seq1 | genus   | 560   | 96.1     | 560
#  ...
#
def make_assignment_table(m6_file, acc2taxid, hsg_list, ranks, cutoffs, rm_self=False):
	def wrap_to_df(seqid, taxids, hsgs, identity):
		tmp = {'rank': [], 'taxid': []}
		ranks = ncbi.get_rank(taxids)
		for taxid, rank in ranks.items():
			tmp['rank'].append(rank)
			tmp['taxid'].append(taxid)
		for hsg_id in hsgs:
			tmp['rank'].append('HSG')
			tmp['taxid'].append(hsg_id)

		df = pd.DataFrame(tmp)
		df['seqid'] = seqid
		df['identity'] = identity

		return df[['seqid', 'rank', 'taxid', 'identity']]

	rank2cutoff = {rank:float(cutoff) for rank, cutoff in zip(ranks, cutoffs)}

	m6 = load_m6(m6_file, rm_self=rm_self)
	tmp = []
	for query in m6['seqid'].unique():
		# find best hits
		subdf = m6.query("seqid == @query")
		maxid = subdf['identity'].max()
		besthits = subdf.query('identity == @maxid')

		# find ancestors
		acc_list = [db.split('__')[0] for db in besthits['db']]
		species_list = np.unique([acc2taxid.loc[acc, 'species_taxid'] for acc in acc_list]).tolist()
		taxids = []
		for species in species_list:
			lineage = ncbi.get_lineage(species)
			taxids += lineage
		taxids = np.unique(taxids).tolist()
		hsgs = np.unique([hsg_list.loc[species_taxid, 'HSG'] for species_taxid in species_list])

		# wrap into dataframe / add HSGs
		df = wrap_to_df(query, taxids, hsgs, maxid)
		tmp.append(df)

	# filter unintended ranks
	assigned = pd.concat(tmp)
	assigned = assigned.query('rank in @rank2cutoff').reset_index(drop=True)
	assigned['seqid'] = assigned['seqid'].astype('category')
	assigned['rank'] = pd.Categorical(assigned['rank'], categories=ranks, ordered=True)

	# filter low identity
	assigned['cutoff'] = [rank2cutoff[rank] for rank in assigned['rank']]
	assigned['confident'] = assigned['identity'] >= assigned['cutoff']
	assigned['assigned'] = [row['taxid'] if row['confident'] else -1 for _, row in assigned.iterrows()]
	assigned = assigned.drop(['cutoff', 'confident'], axis=1)

	# sort assignment table
	assigned = assigned.sort_values(['seqid', 'rank'])

	return assigned

def expand_to_uc(assigned, uc_file):
	uc = load_uc(uc_file)
	expanded = pd.merge(left=assigned, right=uc, left_on='seqid', right_on='rep', how='outer')
	expanded = expanded.drop(['seqid', 'rep'], axis=1)
	expanded = expanded.sort_values(['query', 'rank'])

	expanded = expanded[['query', 'rank', 'taxid', 'identity', 'assigned']]
	expanded.columns = ['seqid', 'rank', 'taxid', 'identity', 'assigned']

	return expanded

def taxid_to_name(taxid):
	if isinstance(taxid, str) and ":" in taxid:
		gid, hsg_n = taxid.split(':')
		name = ncbi.get_taxid_translator([gid])[int(gid)]
		name = f"{name}:{hsg_n}"
	else:
		taxid = int(taxid)
		name = ncbi.get_taxid_translator([taxid])
		name = name[taxid] if taxid in name else "unclassified"

	return name

def assignment(fargs):
	m6_file, args = fargs

	# LOAD data
	acc2taxid = pd.read_csv(args.acc2taxid, sep='\t')
	acc2taxid = acc2taxid.drop_duplicates("assembly_accession").set_index("assembly_accession")
	hsg_list = pd.read_csv(args.hsg, sep=',', index_col='species')

	# make assignments
	assigned = make_assignment_table(m6_file, acc2taxid, hsg_list, args.ranks, args.cutoffs, args.rm_self)
	if args.u_dir is not None:
		uc_file = f"{args.u_dir}/{os.path.basename(m6_file).replace(args.suffix, '')}{args.uc_suffix}"
		assigned = expand_to_uc(assigned, uc_file)

	# add scientific name
	if args.add_name:
		assigned['scientific_name'] = [taxid_to_name(taxid) for taxid in assigned['assigned']]

	# save results
	out_name = f"{args.o_dir}/{os.path.basename(m6_file).replace(args.suffix, '.ta.csv')}"
	assigned.to_csv(out_name, sep=',', index=False)

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

	assert len(args.ranks) == len(args.cutoffs), "Length mismatch between ranks and cutoffs"
	assert os.path.isfile(args.acc2taxid), "acc2taxid file is missing"
	assert os.path.isfile(args.hsg), "hsg file is missing"

	logging.info(f"[{pname}] Make taxonomy assignments")
	targets = glob.glob(f"{args.i_dir}/*{args.suffix}")
	cmds = []
	for m6_file in targets:
		cmds.append((m6_file, args))
	auto_run(assignment, cmds, args.p, quiet=args.quiet, desc=f'taxonomy assignment')

	logging.info(f"[{pname}] Done. (taxonomy assignment)")
	return 0

if __name__ == "__main__":
	exit(main())

#!/usr/bin/env python

import sys
import argparse
import os.path
import logging

import itertools as it
import pandas as pd
from tqdm import tqdm

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
	parser.add_argument("-o", "--output",
		dest="output",
		metavar="STR",
		help="output file name; default stdout")
	parser.add_argument("-m", "--method",
		dest="method",
		choices=('median', 'max', 'shared_ratio'),
		default='median',
		help="aggregating method; default median")
	parser.add_argument("--by-strain",
		dest="by_strain",
		action="store_true",
		help="aggregate into the strains")

	# required arguments
	req_group = parser.add_argument_group("required arguments")
	req_group.add_argument("-i", "--m6",
		dest='m6',
		metavar="M6",
		required=True,
		help="allpair alignment results; required")
	req_group.add_argument("-a", "--acc2taxid",
		dest="acc2taxid",
		metavar="TXT",
		required=True,
		help="accession to taxid mapping file; required")

	return parser.parse_args(argv)

def load_m6(file_name):
	try:
		m6 = pd.read_csv(file_name, sep='\t', header=None)
		m6.columns = ['query', 'subject', 'identity', 'length', 'mismatch', 'gapopen', 'qeury_s', 'query_e', 'subject_s', 'subject_e', 'evalue', 'bitscore']

		return m6
	except:
		return None

def get_acc(seqid):
	return seqid.split('__')[1].split('.')[0]

def make_raw_table(m6, disable=False):
	tmp = {}
	for _, row in tqdm(m6.iterrows(), total=m6.shape[0], desc='make raw table', disable=disable):
		query, subject = row['query'], row['subject']
		if query not in tmp:
			tmp[query] = {query: 100.0}
		if subject not in tmp:
			tmp[subject] = {subject: 100.0}

		tmp[query][subject] = tmp[subject][query] = row['identity']
	
	return pd.DataFrame(tmp).fillna(0.0)

def agg_by_species(df, acc2taxid, method, disable=False):
	df['species'] = [acc2taxid.loc[get_acc(sid), 'species_taxid'] for sid in tqdm(df.index, desc='convert species', disable=disable)]

	return df.groupby('species').median() if method == 'median' else df.groupby('species').max()

def agg_by_strain(df, acc2taxid, method, disable=False):
	df['strain'] = [sid.split('__')[0] for sid in df.index]

	return df.groupby('strain').median() if method == 'median' else df.groupby('strain').max()

def make_shared_ratio_table(m6, acc2taxid, logging):
	def shared_ratio(m6, sa, sb, _self=False):
		def div(a, b, _self=False):
			return (a / b) if b != 0 else 100.0 if _self else 0.0

		subdf1 = m6.query("query_tax == @sa and subject_tax == @sb")
		subdf2 = m6.query("query_tax == @sb and subject_tax == @sa")
		subdf = pd.concat([subdf1, subdf2])

		return div(subdf.query("identity >= 99.0").shape[0], subdf.shape[0], _self=_self)

	m6['query_tax'] = m6['query'].apply(lambda x: acc2taxid.loc[get_acc(x), 'species_taxid'])
	m6['subject_tax'] = m6['subject'].apply(lambda x: acc2taxid.loc[get_acc(x), 'species_taxid'])

	species = m6['query_tax'].to_list() + m6['subject_tax'].to_list()
	species = list(set(species))

	tmp = {}
	for sp in species:
		tmp[sp] = {sp: shared_ratio(m6, sp, sp, _self=True)}

	for sa, sb in it.combinations(species, 2):
		tmp[sa][sb] = tmp[sb][sa] = shared_ratio(m6, sa, sb, _self=False)

	return pd.DataFrame(tmp)


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

	acc2taxid = pd.read_csv(args.acc2taxid, sep='\t', index_col='accession')
	m6 = load_m6(args.m6)
	if m6 is None:
		print(f"{args.m6} is empty")
		logging.info(f"Done. ({pname})")
		return -1

	if args.method in ['median', 'max']:
		raw = make_raw_table(m6, disable=args.quiet)
		logging.info(f"raw_table: {raw.shape}")
		
		func = agg_by_strain if args.by_strain else agg_by_species
		df = func(raw, acc2taxid, args.method, disable=args.quiet)
		df = func(df.T, acc2taxid, args.method, disable=args.quiet)
		logging.info(f"{'strain' if args.by_strain else 'species'} table: {df.shape}")
	else:
		df = make_shared_ratio_table(m6, acc2taxid, logging)

	df.to_csv(args.output if args.output else sys.stdout, sep=',')

	logging.info(f"Done. ({pname})")
	return 0

# main
if __name__ == "__main__":
	exit(main())

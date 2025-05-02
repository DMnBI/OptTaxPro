#!/usr/bin/env python

import sys
import glob
import argparse
import os.path
import logging
from tqdm import tqdm

import pandas as pd
from ete3 import NCBITaxa

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

	# required arguments
	req_group = parser.add_argument_group("required arguments")
	req_group.add_argument("-i", "--in-dir",
		dest="i_dir",
		metavar="PATH",
		required=True,
		help="input directory involving target files; required")
	req_group.add_argument("-a", "--acc2taxid",
		dest="acc2taxid",
		metavar="TXT",
		required=True,
		help="accession to taxid mapping file; required")

	return parser.parse_args(argv)

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
	rpath, _ = os.path.split(os.path.realpath(__file__))

	items = glob.glob(f"{args.i_dir}/*")

	tmp = []
	for item in items:
		genus = item.split('/')[-1].split('.')[0]
		df = pd.read_csv(item, sep=',')
		
		vc = df['single'].value_counts()
		groups = vc.loc[vc > 1].index

	#	df['HSG'] = [f"{genus}:{row['single']}" if row['single'] in groups else row['species'] for _, row in df.iterrows()]
		df['HSG'] = [f"{genus}:{row['single']}" for _, row in df.iterrows()]
		df['singleton'] = [row['single'] not in groups for _, row in df.iterrows()]
		
		tmp.append(df.drop(['single'], axis=1))

	df = pd.concat(tmp)
	df = df.set_index('species')

	ncbi = NCBITaxa()
	species_list = pd.read_csv(args.acc2taxid, sep='\t', index_col='accession')['species_taxid'].unique()
	for species in species_list:
		if species not in df.index:
			lineage = ncbi.get_lineage(species)
			rank = ncbi.get_rank(lineage)
			r2t = {r:t for t, r in rank.items()}
			genus = r2t['genus']
			df.loc[species] = [ncbi.get_taxid_translator([species])[species], f"{genus}:0", True]

	def convert_HSG_name(HSG, ncbi):
		def get_scientific_name(tid, ncbi):
			tid = int(tid)
			return ncbi.get_taxid_translator([tid])[tid]

		if isinstance(HSG, str) and ':' in HSG:
			genus, hsg = HSG.split(':')
			return f"{get_scientific_name(genus, ncbi)}:HSG{hsg}"
		else:
			return get_scientific_name(HSG, ncbi)

	df['HSG name'] = [convert_HSG_name(hsg, ncbi) for hsg in df['HSG']]

	df.to_csv(args.output, sep=',', index=True)

	logging.info(f"Done. ({pname})")
	return 0

# main
if __name__ == "__main__":
	exit(main())

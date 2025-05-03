#!/usr/bin/env python

import sys
import argparse
import os.path
import logging

import pandas as pd
import numpy as np

from sklearn_extra.cluster import KMedoids
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score
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
	parser.add_argument("-m", "--method",
		dest="method",
		choices=("identity", "shared_ratio"),
		default="identity",
		help="metric of distance between two species; default identity")

	# required arguments
	req_group = parser.add_argument_group("required arguments")
	req_group.add_argument("-i", "--idmat",
		dest="idmat",
		metavar="CSV",
		required=True,
		help="species by species identity matrix; required")

	return parser.parse_args(argv)

def load_mat(file_name, method):
	mat = pd.read_csv(file_name, sep=',', index_col=0)
	mat = (100 if method == "identity" else 1) - mat

	diag0 = mat.copy()
	np.fill_diagonal(diag0.to_numpy(), 0)

	return mat, diag0

def PAM_clustering(mat):
	tmp = {"k": [], "Silhouette": []}
	for k in range(2, mat.shape[0]):
		PAM = KMedoids(n_clusters=k, metric='precomputed', method='pam', init='build', max_iter=10000)
		bins = PAM.fit(mat)

		silhouette = silhouette_score(mat, bins.labels_, metric='precomputed')

		tmp['k'].append(k)
		tmp['Silhouette'].append(silhouette)

	df = pd.DataFrame(tmp)
	opt_k = tmp['k'][df['Silhouette'].argmax()]

	PAM = KMedoids(n_clusters=opt_k, metric='precomputed', method='pam', init='build', max_iter=10000)
	bins = PAM.fit(mat)

	df = pd.DataFrame({"PAM": bins.labels_}, index=mat.index)
	df.index.name = 'species'

	return df

def single_linkage(mat, method):
	single = AgglomerativeClustering(
			n_clusters=None,
			metric='precomputed',
			distance_threshold=(1.0 if method == 'identity' else 0.7) + 1e-4,
			linkage='single'
		)

	bins = single.fit(mat)

	df = pd.DataFrame({"single": bins.labels_}, index=mat.index)
	df.index.name='species'

	return df

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

	ncbi = NCBITaxa()
	mat, diag0 = load_mat(args.idmat, method=args.method)

	try:
#		pam = PAM_clustering(diag0)
		single = single_linkage(mat, args.method)
	except:
		return -1

	df = single
#	df = pd.merge(left=pam, right=single, left_index=True, right_index=True)
	df['name'] = [ncbi.get_taxid_translator([x])[x] for x in df.index]

	df.to_csv(args.output if args.output else sys.stdout, sep=',')

	logging.info(f"Done. ({pname})")
	return 0

# main
if __name__ == "__main__":
	exit(main())

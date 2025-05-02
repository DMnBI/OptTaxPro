#!/usr/bin/env python

import sys
import argparse
import os.path
import logging

import glob
import tempfile

import subprocess as sp
import multiprocessing as mp
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
	parser.add_argument("-t", "--tmp-dir",
		dest="t_dir",
		metavar="PATH",
		default='./',
		help="a path for temporary directory; default ./")
	parser.add_argument("-c", "--concurrent",
		dest="c",
		metavar="INT",
		type=int,
		default=1,
		help="a number of processes which run concurrently; default 1")

	# required arguments
	req_group = parser.add_argument_group("required arguments")
	req_group.add_argument("-s", "--sequence",
		dest="sequence",
		metavar="PATH",
		required=True,
		help="a directory involving rna files; required")
	req_group.add_argument("-a", "--acc2taxid",
		dest="acc2taxid",
		metavar="TXT",
		required=True,
		help="accession to taxid mapping file; required")
	req_group.add_argument("-o", "--output",
		dest="output",
		metavar="STR",
		required=True,
		help="output file name; required")

	return parser.parse_args(argv)

def auto_run(cmds, c, desc=None):
	pool = mp.Pool(c)
	for _ in tqdm(pool.imap_unordered(sp.run, cmds), desc=desc, total=len(cmds)):
		pass
	pool.close()
	pool.join()

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

	t_dir = tempfile.mkdtemp(dir=args.t_dir, prefix='4_build_HSG_')
	EACH_GENUS = f"{t_dir}/each_genus"
	M6_DIR = f"{t_dir}/m6"
	IDMAT_DIR = f"{t_dir}/idmat"
	HSG_DIR = f"{t_dir}/HSG"

	# All pair-wise alignment for each genus
	os.makedirs(EACH_GENUS, exist_ok=True)
	cmd = [f"{rpath}/utils/split_by_rank.py",
		'-r', args.sequence,
		'-a', args.acc2taxid,
		'-o', EACH_GENUS,
		'--rank', 'genus'
	]
	_ = sp.run(cmd)

	os.makedirs(M6_DIR, exist_ok=True)
	targets = glob.glob(f"{EACH_GENUS}/*.fna")
	cmds = []

	for target in targets:
		outname = os.path.split(target)[-1].replace(".fna", '.m6')
		cmd = ['vsearch', '--quiet',
			'--allpairs_global', target,
			'--blast6out', f"{M6_DIR}/{outname}",
			'--acceptall',
			'--match', '5', '--mismatch', '-5',
			'--gapopen', '5I/0E', '--gapext', '5I/0E',
			'--iddef', '1',
			'--threads', '20',
		]
		cmds.append(cmd)

	auto_run(cmds, 4, desc="Intra-genus search")

	# Build identity matrix for each genus
	os.makedirs(IDMAT_DIR, exist_ok=True)
	m6_list = glob.glob(f"{M6_DIR}/*.m6")
	cmds = []
	for m6 in m6_list:
		taxid = os.path.split(m6)[-1].split('.')[0]
		cmd = ['python', f"{rpath}/utils/make_species_idmat.py",
			'-i', m6,
			'-a', args.acc2taxid,
			'-m', 'max',
			'-o', f"{IDMAT_DIR}/{taxid}.max.idmat.csv",
			'-q',
		]
		cmds.append(cmd)

	auto_run(cmds, args.c, desc="Make species id matrice")

	# Define homologous species group for each genus
	os.makedirs(HSG_DIR, exist_ok=True)
	idmats = glob.glob(f"{IDMAT_DIR}/*.csv")
	cmds = []
	for idmat in idmats:
		taxid = os.path.split(idmat)[-1].split('.')[0]
		cmd = ['python', f"{rpath}/utils/make_HSG.py",
			'-i', idmat,
			'-m', 'identity',
			'-o', f"{HSG_DIR}/{taxid}.HSG.csv",
			'-q'
		]
		cmds.append(cmd)

	auto_run(cmds, args.c, desc="Build HGS for each genus")

	# Concat HSGs into a single file
	cmd = [f'{rpath}/utils/merge_HSG.py',
		'-i', HSG_DIR,
		'-a', args.acc2taxid,
		'-o', args.output,
	]
	_ = sp.run(cmd)

	logging.info(f"Done. ({pname})")
	return 0

# main
if __name__ == "__main__":
	exit(main())

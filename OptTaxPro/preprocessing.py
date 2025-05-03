#!/usr/bin/env python

import os
import sys
import argparse
import logging
import glob
import tempfile
import json

import pandas as pd
import numpy as np
import subprocess as sp
from tqdm import tqdm
from Bio import SeqIO

from utils.commons import auto_run, run_cmd

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
		default=".fastq",
		help="suffix of input files; default .fastq")
	parser.add_argument("--t_dir",
		dest="t_dir",
		metavar="PATH",
		default="./",
		help="temporal directory to saving intermediate files; default ./")

	# cutadapt arguments
	cut_group = parser.add_argument_group("cutadapt arguments")
	cut_group.add_argument("-m", "--min-len",
		dest='m',
		metavar="INT",
		type=int,
		default=1300,
		help="minimum length to retain; default 1300")
	cut_group.add_argument("-M", "--max-len",
		dest='M',
		metavar="INT",
		type=int,
		default=1850,
		help="maximum length to retain; default 1850")
	cut_group.add_argument("-O",
		dest='O',
		metavar='INT',
		type=int,
		default=15,
		help="minimum overlap length between read and adapter; deault 15")
	cut_group.add_argument("-e",
		dest='e',
		metavar="FLOAT",
		type=float,
		default=0.15,
		help="maximum allowed error rate; default 0.15")

	# quality control arguments
	qc_group = parser.add_argument_group("quality control arguments")
	qc_group.add_argument("--qc-cutoff",
		dest='qc_cutoff',
		metavar="FLOAT",
		type=float,
		default=0.99,
		help="minimum average read quality; default 0.99")

	# singleton arguments
	single_group = parser.add_argument_group("singleton removal arguments")
	single_group.add_argument("--skip-singleton-removal",
		dest="skip_single",
		action="store_true",
		help='Skip singleton removal process')
	single_group.add_argument("--otu-cutoff",
		dest='otu_cutoff',
		metavar="FLOAT",
		type=float,
		default=0.97,
		help="OTU clustering cutoff; default 0.97")
	single_group.add_argument("--threads",
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

def quality_control(args):
	def calc_accuracy(phred_score):
		P = 10 ** (-phred_score / 10)
		return 1 - P

	trimmed, qc, cutoff = args
	quality_info = open(trimmed.replace(".fastq", ".quality.txt"), 'w')
	print("\t".join(['seqid', 'mean_quality']), file=quality_info)
	passed = []
	for record in SeqIO.parse(trimmed, "fastq"):
		phred_scores = record.letter_annotations["phred_quality"]
		read_quality = np.mean([calc_accuracy(phred_score) for phred_score in phred_scores])
		if read_quality >= cutoff:
			passed.append(record)

		print(f"{record.id}\t{read_quality}", file=quality_info)
	quality_info.close()

	SeqIO.write(passed, qc, "fastq")

def remove_singleton(args):
	fastq, uc, output = args

	header = ['record_type', 'cluster_no', 'slen', 'identity', 'strand', 'dump1', 'dump2', 'CIGAR', 'query', 'target']
	uc = pd.read_csv(uc, sep='\t', header=None, names=header)
	singletons = uc.query("record_type == 'C' and slen == 1")['query'].tolist()
	singletons = {sid:True for sid in singletons}

	records = SeqIO.parse(fastq, "fastq")
	removed = [record for record in records if record.id not in singletons]
	SeqIO.write(removed, output, 'fastq')

def summary_preprocessing(o_dir, t_dir, qc_cutoff, summary):
	def load_json(json_file, sample_name):
		data = json.load(open(json_file))
		
		df = pd.json_normalize(data)
		cols = ['read_counts.input'] + [f"read_counts.filtered.{col}" for col in ['too_short', 'too_long', 'too_many_n', 'discard_untrimmed']]
		df = df[cols]
		df.columns = ['raw', 'short', 'long', 'N', 'no primer']
		df['sample'] = sample_name

		return df

	def load_quality(quality_file, qc_cutoff):
		quality = pd.read_csv(quality_file, sep='\t')
		low = quality.query("mean_quality < @qc_cutoff")
		return {'low_quality': low.shape[0], 'high_quality': quality.shape[0] - low.shape[0]}

	tmp = []
	for json_file in glob.glob(f"{t_dir}/*.json"):
		sample_name = os.path.basename(json_file).replace(".trimmed.json", "")
		subdf = load_json(json_file, sample_name)

		quality_file = f"{t_dir}/{sample_name}.trimmed.quality.txt"
		qc = load_quality(quality_file, qc_cutoff)
		subdf['low_quality'] = qc['low_quality']

		final_fastq = f"{o_dir}/{sample_name}.clean.fastq"
		nseqs = len(SeqIO.index(final_fastq, 'fastq'))
		subdf['singeton'] = qc['high_quality'] - nseqs
		subdf['clean'] = nseqs

		tmp.append(subdf)
	df = pd.concat(tmp)

	final_cols = ['sample', 'raw', 'short', 'long', 'N', 'no primer', 'low_quality', 'singeton', 'clean']
	df = df[final_cols].sort_values('sample')
	df.to_csv(summary, index=False)

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
	t_dir = tempfile.mkdtemp(dir=args.t_dir, prefix="99_preprocessing_", suffix="_tmp")
	targets = glob.glob(f"{args.i_dir}/*{args.suffix}")
	

	# run cutadapt
	logging.info(f"[{pname}] RUN cutadapt")
	cmds = []
	trimmed_list = []
	for fastq in targets:
		trimmed = f"{t_dir}/{os.path.basename(fastq).replace(args.suffix, '.trimmed.fastq')}"
		trimmed_list.append(trimmed)

		json_file = trimmed.replace(".fastq", ".json")
		cmd = ['cutadapt', '--quiet',
			'-g', 'AGRGTTYGATYMTGGCTCAG...AAGTCGTAACAAGGTARCY',
			'--rc', '--trimmed-only', '-O', str(args.O),
			'-m', str(args.m), '-M', str(args.M),
			'--max-n', '0',
			'-j', str(args.p), '-e', str(args.e),
			'-o', trimmed, '--json', json_file,
			fastq
		]
		cmds.append(cmd)
	auto_run(run_cmd, cmds, args.p, quiet=args.quiet, desc='cutadapt')
	logging.info(f"[{pname}] Done. (cutadapt)")

	# run quality control
	logging.info(f"[{pname}] RUN quality control")
	cmds = []
	qc_list = []
	for trimmed in trimmed_list:
		qc = trimmed.replace(".trimmed.fastq", ".qc.fastq")
		qc_list.append(qc)

		cmd = (trimmed, qc, args.qc_cutoff)
		cmds.append(cmd)
	auto_run(quality_control, cmds, args.p, quiet=args.quiet, desc="quality control")
	logging.info(f"[{pname}] Done. (quality control)")

	if not args.skip_single:
		# run singleton removal
		logging.info(f"[{pname}] FIND singletons (cutoff={args.otu_cutoff:.2f}) ")
		cmd = ["python", f"{rpath}/utils/cluster_fast.py",
			'-i', t_dir, '-o', t_dir,
			'-t', str(args.threads), '-p', str(args.p // args.threads),
			'-c', str(args.otu_cutoff), 
			'--suffix', '.qc.fastq'
		]
		if args.quiet:
			cmd += ['quiet']
		if args.log is not None:
			cmd += ['--log', args.log]
		_ = run_cmd(cmd)

		logging.info(f"[{pname}] REMOVE singletons")
		cmds = []
		for qc in qc_list:
			uc = qc.replace('.qc.fastq', '.uc')
			output = f"{args.o_dir}/{os.path.basename(qc).replace('.qc.fastq', '.clean.fastq')}"
			cmds.append((qc, uc, output))
		auto_run(remove_singleton, cmds, args.p, quiet=args.quiet, desc='remove singletons')
		logging.info(f"[{pname}] Done. (singleton removal)")
	else:
		for qc in qc_list:
			output = f"{args.o_dir}/{os.path.basename(qc).replace('.qc.fastq', '.clean.fastq')}"
			sp.run(['ln', '-s', os.path.realpath(qc), output])

	# make preprocessing summary
	summary = f"{args.o_dir}/preprocessing_summary.csv"
	summary_preprocessing(args.o_dir, t_dir, args.qc_cutoff, summary)

	logging.info(f"[{pname}] Done. (preprocessing)")
	return 0

# main
if __name__ == "__main__":
	exit(main())
	
import subprocess as sp
import multiprocessing as mp
from tqdm import tqdm

VSEARCH_BASE = ['vsearch', '--quiet',
	'--match', '5', '--mismatch', '-5',
	'--gapopen', '5I/0E', '--gapext', '5I/0E',
	'--iddef', '1'
]

def run_cmd(cmd, stdout=None):
	return sp.run([str(x) for x in cmd], stdout=stdout)

def auto_run(func, cmds, c, quiet=False, desc=""):
	pool = mp.Pool(c)
	for _ in tqdm(pool.imap_unordered(func, cmds), total=len(cmds), disable=quiet, desc=desc):
		pass
	pool.close()
	pool.join()

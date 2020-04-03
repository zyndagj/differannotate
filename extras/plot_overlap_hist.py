#!/usr/bin/env python
#
###############################################################################
# Author: Greg Zynda
# Last Modified: 03/19/2020
###############################################################################

import sys, os, re
if len(sys.argv) == 2:
	import matplotlib
	matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import multiprocessing as mp

def main():
	dir_path = sys.argv[1]
	all_beds = glob(os.path.join(dir_path, "*bed"))
	cat_re = re.compile(r"/(\w+)_[B\-\+]_")
	categories = set(map(lambda x: cat_re.search(x).group(1), all_beds))
	strands = '+-B'
	sam_re = re.compile(r"[B\-\+]_(\w+)_(\w+)_")
	#samples = set(map(lambda x: sam_re.search(x).group(1), all_beds)+map(lambda x: sam_re.search(x).group(2), all_beds))
	for cat in categories:
		for strand in strands:
			goodGlob = glob(os.path.join(dir_path, "%s_%c_*+B+A*bed"%(cat, strand)))
			samlist = map(lambda x: sam_re.search(x).group(1), goodGlob)
			plt.figure(dpi=200, figsize=(10,4))
			if len(samlist) > 1:
				for sam in sorted(samlist):
					goodF = glob(os.path.join(dir_path, "%s_%c_%s*+B+A.bed"%(cat, strand, sam)))[0]
					badF = glob(os.path.join(dir_path, "%s_%c_%s*+B-a.bed"%(cat, strand, sam)))[0]
					V = np.concatenate([np.loadtxt(f, delimiter='\t', usecols=4) for f in (goodF, badF)])
					if len(V):
						plt.hist(V, bins=40, log=True, alpha=0.6, label=sam, histtype='step', lw=2)
					else:
						print V
				plt.legend(loc=9)
			else:
				goodF = glob(os.path.join(dir_path, "%s_%c_*+B+A*bed"%(cat, strand)))[0]
				badF = glob(os.path.join(dir_path, "%s_%c_*+B-a*bed"%(cat, strand)))[0]
				V = np.concatenate([np.loadtxt(f, delimiter='\t', usecols=4) for f in (goodF, badF)])
				plt.hist(V, bins=40, log=True)
			plt.axvline(75, color='r')
			plt.ylabel("Feature Count")
			plt.xlabel("Reciprocal Overlap")
			plt.title("Reciprocal Overlap of %s on %s"%(cat, '+/-' if strand == 'B' else strand))
			plt.xlim((0,100))
			plt.tight_layout()
			if len(sys.argv) == 3:
				plt.show()
				sys.exit()
			plt.savefig(os.path.join(dir_path, "%s_%c_overlap.png"%(cat, strand)))
			plt.close()

if __name__ == "__main__":
	main()


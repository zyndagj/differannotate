#!/usr/bin/env python
#
###############################################################################
# Author: Greg Zynda
# Last Modified: 04/03/2020
###############################################################################

import numpy as np
import matplotlib
import sys, os, re
from glob import glob
import multiprocessing as mp
from pysam import FastaFile
from Meth5py import Meth5py
from itertools import imap
import multiprocessing as mp
from operator import itemgetter

serial = len(sys.argv) == 3
if not serial: matplotlib.use("Agg")
import matplotlib.pyplot as plt

plot = False

def main():
	dir_path = sys.argv[1]
	all_beds = glob(os.path.join(dir_path, "*bed"))
	cat_re = re.compile(r"/(\w+)_[B\-\+]_")
	categories = set(map(lambda x: cat_re.search(x).group(1), all_beds))
	strands = '+-B'
	if serial: global M5
	M5 = Meth5py('th.bam.mr','arabidopsis_thaliana.fa')
	cons = M5.contexts
	#sam_re = re.compile(r"[B\-\+]_(\w+)_(\w+)_")
	#samples = set(map(lambda x: sam_re.search(x).group(1), all_beds)+map(lambda x: sam_re.search(x).group(2), all_beds))
	print "Feature,Strand,Context,MSE(TP,FP),MSE(TP,FN),MSE(FP,FN)"
	if serial:
		global sd
		global meth
		meth = {c:M5.fetch(c) for c in M5.sorted_chroms}
		sd = {"+":"-", "-":"+"}
		ret = map(gen_fig, [(dir_path, c, s, con, 20) for c in categories for con in cons for s in strands])
	else:
		del M5
		pool = mp.Pool(4, worker_init)
		ret = pool.map(gen_fig, [(dir_path, c, s, con, 20) for c in categories for s in strands for con in cons])
		pool.close()
		pool.join()
	transcript = set(('exon','gene','cds','mrna'))
	for con in cons:
		plt.figure(dpi=200, figsize=(6,4))
		for s in '+-':
			#fL = filter(lambda x: x[0] in transcript and x[1] == s and x[2] == con, ret)
			#fL = filter(lambda x: x[1] == s and x[2] == con, ret)
			fL = filter(lambda x: x[0] in transcript and x[1] == s, ret)
			fX, fY = map(itemgetter(3), fL), map(itemgetter(4), fL)
			plt.scatter(np.sqrt(fX), np.sqrt(fY), label=s)
			#plt.scatter(fX, fY, label=s)
		plt.legend()
		plt.title("%s Methylation"%(con))
		plt.xlabel("median(|TP-FP|)")
		plt.ylabel("median(|TP-FN|)")
		plt.savefig(con+".png")
		plt.close()

def worker_init():
	global M5
	global sd
	global meth
	sd = {"+":"-", "-":"+"}
	M5 = Meth5py('th.bam.mr','arabidopsis_thaliana.fa')
	meth = {c:M5.fetch(c) for c in M5.sorted_chroms}
	print "Opened M5"

def gen_fig(arg_list):
	#print arg_list
	dir_path, cat, strand, con, bins = arg_list
	a11 = glob(os.path.join(dir_path, "%s_%c_*+A-b*bed"%(cat, strand)))[0]
	# Separate with and without overlap
	rnn = glob(os.path.join(dir_path, "%s_%c_*+B-a*bed"%(cat, strand)))[0]
	both = glob(os.path.join(dir_path, "%s_%c_*+B+A*bed"%(cat, strand)))[0]
	plt.figure(dpi=200, figsize=(6,4))
	av = gen_vals(a11, con, strand, M5, bins=20, method=np.mean)
	rv = gen_vals(rnn, con, strand, M5, bins=20, method=np.mean)
	bv = gen_vals(both, con, strand, M5, bins=20, method=np.mean)
	plt.plot(av, label="Araport11")
	plt.plot(rv, label="RNNotate")
	plt.plot(bv, label="Both")
	aa, ra, ba = map(np.array, (av, rv, bv))
	middle = np.mean(av+rv+bv)
	tpfp = np.mean(np.square(ba-ra)).round(6)
	tpfn = np.mean(np.square(ba-aa)).round(6)
	fpfn = np.mean(np.square(aa-ra)).round(6)
	print ','.join(map(str,(cat, strand, con, tpfp, tpfn, fpfn)))
	plt.xlim((0,bins*3-1))
	if strand in sd:
		oav = gen_vals(a11.replace("_%c_"%(strand),"_%c_"%(sd[strand]),1), con, sd[strand], M5, bins=20, method=np.mean)
		orv = gen_vals(rnn.replace("_%c_"%(strand),"_%c_"%(sd[strand]),1), con, sd[strand], M5, bins=20, method=np.mean)
		obv = gen_vals(both.replace("_%c_"%(strand),"_%c_"%(sd[strand]),1), con, sd[strand], M5, bins=20, method=np.mean)
		plt.ylim((0,max(av+rv+bv+oav+orv+obv)*1.05))
	plt.axvline(bins-0.5, color='0.75')
	plt.axvline((bins*2)-0.5, color='0.75')
	plt.xticks(np.arange(bins/2, bins*3, bins), ["1kb Upstream", "Feature", "1kb Downstream"])
	plt.title("Mean %s Methylation of %s on the %s Strand"%(con, cat, strand if strand != 'B' else '+/-'))
	plt.ylabel("Methylation Frequency")
	#plt.xlabel("Relative Location")
	plt.legend(loc=1)
	plt.tight_layout()
	if serial:
		plt.show()
	else:
		plt.savefig(os.path.join(dir_path, "%s_%c_%s_methylation.png"%(cat, strand, con)))
	plt.close()
	return (cat, strand, con, tpfp, tpfn)

def gen_vals(infile, con, strand, M5, bins=10, min_overlap=0, max_overlap=100, method=np.mean):
	IF = open(infile,'r')
	before_list = [[] for i in range(bins)]
	feat_list = [[] for i in range(bins)]
	after_list = [[] for i in range(bins)]
	for rec in imap(lambda x: x.rstrip('\n').split('\t'), filter(lambda x: x[0] != "#", IF)):
		c,s,e,st = rec[0], int(rec[1]), int(rec[2]), rec[3]
		if st == '+':
			#append_list(filter_meth(M5.fetch(c,s-999,s), con, strand, M5, bins, method), before_list)
			#append_list(filter_meth(M5.fetch(c,s+1,e), con, strand, M5, bins, method), feat_list)
			#append_list(filter_meth(M5.fetch(c,e+1,e+1000), con, strand, M5, bins, method), after_list)
			append_list(filter_meth(meth[c][s-1000:s], con, strand, M5, bins, method), before_list)
			append_list(filter_meth(meth[c][s:e], con, strand, M5, bins, method), feat_list)
			append_list(filter_meth(meth[c][e:e+1000], con, strand, M5, bins, method), after_list)
		else:
			#append_list(filter_meth(M5.fetch(c,s-999,s), con, strand, M5, bins, method)[::-1], after_list)
			#append_list(filter_meth(M5.fetch(c,s+1,e), con, strand, M5, bins, method)[::-1], feat_list)
			#append_list(filter_meth(M5.fetch(c,e+1,e+1000), con, strand, M5, bins, method)[::-1], before_list)
			append_list(filter_meth(meth[c][s-1000:s], con, strand, M5, bins, method)[::-1], after_list)
			append_list(filter_meth(meth[c][s:e], con, strand, M5, bins, method)[::-1], feat_list)
			append_list(filter_meth(meth[c][e:e+1000], con, strand, M5, bins, method)[::-1], before_list)
	before_vals = map(method, before_list)
	feat_vals = map(method, feat_list)
	after_vals = map(method, after_list)
	#for l in (before_list, feat_list, after_list):
	#	print map(len, l)
	#print np.round(before_vals+feat_vals+after_vals,3)
	IF.close()
	return before_vals+feat_vals+after_vals

def filter_meth(meth_list, context, strand, M5, bins=20, method=np.mean):
	ci = M5.contexts.index(context)
	if strand != 'B': si = M5.strands.index(strand)
	meth_a = np.array(meth_list,dtype=np.float)
	split_a = np.array_split(meth_a, bins)
	ret_a = []
	for a in split_a:
		if strand == 'B':
			view = a[a[:,0] == ci]
		else:
			view = a[np.logical_and(a[:,0] == ci, a[:,1] == si)]
		ratios = view[:,2]/view[:,3]
		ret_a.append(method(ratios) if len(ratios) else -1)
	return ret_a

def append_list(ret, accum):
	for i,v in enumerate(ret):
		if v != -1:
			accum[i].append(v)

if __name__ == "__main__":
	main()


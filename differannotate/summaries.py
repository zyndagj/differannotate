#!/usr/bin/env python
#
###############################################################################
# Author: Greg Zynda
# Last Modified: 01/26/2020
###############################################################################
# BSD 3-Clause License
# 
# Copyright (c) 2019, Greg Zynda
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
# 
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# 
# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
###############################################################################

import logging, sys
from .constants import FORMAT, BaseIndex
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.WARN, format=FORMAT)

from . import comparisons
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import numpy as np
from time import time

def tabular_region(GI, p=95, fig_ext='png', temd=False):
	chrom_set = GI.get_chrom_set()	# intersecting set chroms from all files
	max_chrom_len = max(map(len, chrom_set)+[len("Chrom")])
	max_elem_len = max(map(len, list(GI.element_dict)+list(GI.order_dict)+list(GI.sufam_dict)))
	max_name_len = max(map(len, list(GI.gff3_names)))
	#feature_set = # features only in reference
	non_te_elements = set(GI.element_dict) - GI.te_names
	te_elements = set(GI.element_dict) & GI.te_names
	for chrom in chrom_set:
		_print_table_region(GI, chrom, GI.element_dict, 1, max_chrom_len, \
			max_elem_len, max_name_len, p, fig_ext)
		if temd:
			_print_table_region(GI, chrom, GI.order_dict, 2, max_chrom_len, \
				max_elem_len, max_name_len, p, fig_ext)
			_print_table_region(GI, chrom, GI.sufam_dict, 3, max_chrom_len, \
				max_elem_len, max_name_len, p, fig_ext)
def _print_table_region(GI, chrom, elem_list, col, mcl, mel, mnl, p=95, fig_ext='png',fmt='csv'):
	target = ("", "Element", "TE_Order", "TE_Superfamily")
	mel = max(map(len, target)+[mel])
	header = ("Chrom","S",target[col],"Sample","TP", "FP", "FN", "SENS", "PREC")
	if fmt == 'csv':
		template = ','.join(["{}"]*9)
	elif fmt == 'tsv':
		template = "{:<{mcl}} {:^3} {:<{mel}} {:<{mn}} "+' '.join(["{:>5}"]*5)
	else:
		sys.exit("Unhandled format: %s"%(fmt))
	print(template.format(*header, mcl=mcl, mn=mnl, mel=mel))
	cname = GI.gff3_names[0]
	for elem in elem_list:
		eid = elem_list[elem]
		#print col, eid
		for sstr, sval in zip(('+/-','+','-'), (False, '+', '-')):
			#if not np.sum(A): continue
			length_array_dict = {}
			proportion_array_dict = {}
			tp_list, fp_list = [], []
			for i, name in enumerate(GI.gff3_names):
				Ab, aB, AB = GI.calc_intersect_2(chrom, cname, name, eid, col, p, strand=sval)
				#if col == 2:
				#	sAb, saB, sAB = GI.calc_intersect_2(chrom, cname, name, eid, col, p, strand=sval, ret_set=True)
				#	print "%s-%s"%(cname.upper(), name), sAb
				#	print "%s-%s"%(cname, name.upper()), saB
				#	print "%s-%s"%(cname.upper(), name.upper()), sAB
				tp, fp, fn, sen, pre = _calc_stats_region(Ab, aB, AB)
				if i != 0:
					tp_list.append(tp)
					fp_list.append(fp)
				else:
					if not tp and not fp: break
				if not i:
					print(template.format(chrom, sstr, elem, name, tp,fp,fn,sen,pre, mcl=mcl, mn=mnl, mel=mel))
				else:
					print(template.format('', '', '', name, tp,fp,fn,sen,pre, mcl=mcl, mn=mnl, mel=mel))
				# Store array of interval lengths
				length_array_dict[name] = GI.get_length_array(chrom, name, eid, col, sval)
				if GI.FA:
					proportion_array_dict[name] = GI.get_proportion_arrays(chrom, name, eid, col, sval)
			if not fig_ext or not (sum(tp_list) or sum(fp_list)): continue
			sstrand = 'B' if sstr == '+/-' else sstr
			# Generate length boxplot
			fig_name = "length_bp_%s_%s_%s.%s"%(chrom, sstrand, elem, fig_ext)
			plt.figure(dpi=200)
			assert(np.all(length_array_dict[length_array_dict.keys()[0]] == length_array_dict.values()[0]))
			len_list = [length_array_dict[name] for name in GI.gff3_names]
			plt.boxplot(map(np.sqrt, len_list), labels=GI.gff3_names)
			plt.ylabel('sqrt(length)')
			plt.title('%s %s %s length distribution'%(chrom, sstr, elem))
			plt.savefig(fig_name)
			del len_list
			plt.close()
			# Generate proportion boxplot
			if GI.FA:
				fig_name = "proportion_%s_%s_%s.%s"%(chrom, sstrand, elem, fig_ext)
				plt.figure(dpi=200)
				colors = {'G':'gold', 'T':'tomato', 'A':'seagreen', 'C':'royalblue'}
				space, width = 2.0/33, 4.0/33
				pos_stop = space/2.0+width/2.0
				pos = np.array([-3.0*pos_stop, -pos_stop, pos_stop, 3*pos_stop])
				bpl = []
				for i,name in enumerate(GI.gff3_names):
					for j in range(4):
						bp = plt.boxplot(proportion_array_dict[name][j], positions=[pos[j]+(i+1.0)], patch_artist=True, widths=[width], showfliers=False)
						bpl.append(bp)
						plt.setp(bpl[-1]["boxes"], facecolor=colors[BaseIndex[j]])
				plt.xticks(np.arange(len(GI.gff3_names))+1, GI.gff3_names)
				max_pro = 0.5
				#for pa in proportion_array_dict.values():
				#	for a in pa:
				#		if len(a) > 0: max_pro = max(max_pro, max(a))
				for bp in bpl:
					max_pro = max(max_pro, max([max(w.get_ydata()) for w in bp['whiskers']]))
				plt.legend([bp["boxes"][0] for bp in bpl[:4]], [BaseIndex[i] for i in range(4)], loc='upper right')
				plt.ylim(0, min(max_pro*1.2, 1))
				plt.xlim(0.5, len(GI.gff3_names)+0.5)
				plt.ylabel('Proportion')
				plt.title('%s %s %s Nucleotide Proportion'%(chrom, sstr, elem))
				plt.savefig(fig_name)
				plt.close()
			# Generate venn figure
			if len(GI.gff3_names) not in (2,3): continue
			fig_name = "region_%s_%s_%s.%s"%(chrom, sstrand, elem, fig_ext)
			logger.debug("Generating %s"%(fig_name))
			if len(GI.gff3_names) in (2,3):
				plt.figure(figsize=(4,4), dpi=200)
				plt.title("%s %s %s"%(chrom, sstr, elem))
			if len(GI.gff3_names) == 2: # (Ab, aB, AB)
				n1, n2 = GI.gff3_names
				ret = GI.calc_intersect_2(chrom, n1, n2, eid, col, p, strand=sval)
				if Ab or aB or AB:
					venn2(subsets=ret, set_labels=GI.gff3_names)
				else:
					logger.warn("Empty plot for %s_%s_%s"%(chrom, sstrand, elem))
			elif len(GI.gff3_names) == 3: # (Abc, aBc, ABc, abC, AbC, aBC, ABC)
				n1, n2, n3 = GI.gff3_names
				ret = GI.calc_intersect_3(chrom, n1, n2, n3, eid, col, p, strand=sval)
				if np.sum(ret):
					venn3(subsets=ret, set_labels=GI.gff3_names)
				else:
					strand = 'B' if sstr == '+/-' else sstr
					logger.warn("Empty plot for %s_%s_%s"%(chrom, strand, elem))
			if len(GI.gff3_names) in (2,3):
				plt.savefig(fig_name)
				plt.close()
	print("")
	logger.info("Finished region table")

def tabular(GI, strand=True, fig_ext='png', temd=False):
	chrom_set = GI.get_chrom_set()	# intersecting set chroms from all files
	max_chrom_len = max(map(len, chrom_set)+[len("Chrom")])
	max_elem_len = max(map(len, list(GI.element_dict)+list(GI.order_dict)+list(GI.sufam_dict)))
	max_name_len = max(map(len, list(GI.gff3_names)))
	#feature_set = # features only in reference
	non_te_elements = set(GI.element_dict)-GI.te_names
	te_elements = set(GI.element_dict) & GI.te_names
	for chrom in chrom_set:
		_print_table(GI, chrom, GI.element_dict, 1, max_chrom_len, \
			max_elem_len, max_name_len, fig_ext)
		if temd:
			_print_table(GI, chrom, GI.order_dict, 2, max_chrom_len, \
				max_elem_len, max_name_len, fig_ext)
			_print_table(GI, chrom, GI.sufam_dict, 3, max_chrom_len, \
				max_elem_len, max_name_len, fig_ext)
				
def _print_table(GI, chrom, elem_list, col, mcl, mel, mnl, fig_ext='png',fmt='csv'):
	target = ("", "Element", "TE_Order", "TE_Superfamily")
	mel = max(map(len, target)+[mel])
	header = ("Chrom","S",target[col],"Sample","TP", "FP", "TN", "FN", "SENS", "SPEC", "PREC")
	if fmt == 'csv':
		template = ','.join(["{}"]*11)
	elif fmt == 'tsv':
		template = "{:<{mcl}} {:^3} {:<{mel}} {:<{mn}} "+' '.join(["{:>8}"]*7)
	else:
		sys.exit("Unhandled format: %s"%(fmt))
	print(template.format(*header, mcl=mcl, mn=mnl, mel=mel))
	for elem in elem_list:
		fa, ra, ba = _gen_arrays(GI, chrom, elem, col)
		for s, A in zip(('+/-','+','-'), (ba,fa,ra)):
			#if not np.sum(A): continue
			tp, fp, tn, fn, sen, spe, pre = _calc_stats(A)
			if not tp[0]+fp[0]: continue
			for i, name in enumerate(GI.gff3_names):
				if not i:
					print(template.format(chrom, s, elem, name, tp[i],fp[i],tn[i],fn[i],sen[i],spe[i],pre[i], mcl=mcl, mn=mnl, mel=mel))
				else:
					print(template.format('','','', name, tp[i],fp[i],tn[i],fn[i],sen[i],spe[i],pre[i], mcl=mcl, mn=mnl, mel=mel))
			if not fig_ext or A.shape[0] not in (2,3) or not (tp[1:].sum() or fp[1:].sum()): continue
			# Generate figure
			strand = 'B' if s == '+/-' else s
			fig_name = "base_%s_%s_%s.%s"%(chrom, strand, elem, fig_ext)
			logger.debug("Generating %s"%(fig_name))
			plt.figure(figsize=(4,4), dpi=200)
			plt.title("%s %s %s"%(chrom, s, elem))
			if A.shape[0] == 2: # (Ab, aB, AB)
				Ab = _venn2_helper(A, 1, 0)
				aB = _venn2_helper(A, 0, 1)
				AB = _venn2_helper(A, 1, 1)
				venn_sets = (Ab, aB, AB)
				assert(venn_sets == (fn[1], fp[1], tp[1]))
				if A.sum():
					venn2(subsets=venn_sets, set_labels=GI.gff3_names)
				else:
					logger.warn("Empty plot for %s_%s_%s"%(chrom, strand, elem))
			elif A.shape[0] == 3: # (Abc, aBc, ABc, abC, AbC, aBC, ABC)
				Abc = _venn3_helper(A, 1, 0, 0)
				aBc = _venn3_helper(A, 0, 1, 0)
				ABc = _venn3_helper(A, 1, 1, 0)
				abC = _venn3_helper(A, 0, 0, 1)
				AbC = _venn3_helper(A, 1, 0, 1)
				aBC = _venn3_helper(A, 0, 1, 1)
				ABC = _venn3_helper(A, 1, 1, 1)
				if A.sum():
					venn3(subsets=(Abc, aBc, ABc, abC, AbC, aBC, ABC), set_labels=GI.gff3_names)
				else:
					logger.warn("Empty plot for %s_%s_%s"%(chrom, strand, elem))
			plt.savefig(fig_name)
			plt.close()
	print("")

def _venn3_helper(array, rv0=1, rv1=0, rv2=0):
	start = time()
	m0 = array[0,:] == rv0
	m1 = array[1,:] == rv1
	m2 = array[2,:] == rv2
	assert(len(m0) == len(m1) and len(m1) == len(m2))
	m0a1 = np.logical_and(m0, m1)
	m0a1a2 = np.logical_and(m0a1, m2)
	ret = np.sum(m0a1a2)
	logger.debug("%.3f seconds"%(time()-start))
	return ret
def _venn2_helper(array, rv0=1, rv1=0):
	start = time()
	m0 = array[0,:] == rv0
	m1 = array[1,:] == rv1
	assert(len(m0) == len(m1))
	m0a1 = np.logical_and(m0, m1)
	ret = np.sum(m0a1)
	logger.debug("%.3f seconds"%(time()-start))
	return ret

def _gen_arrays(GI, chrom, elem, col):
	start = time()
	if col == 1:
		elem_id = GI.element_dict[elem]
	elif col == 2:
		elem_id = GI.order_dict[elem]
	elif col == 3:
		elem_id = GI.sufam_dict[elem]
	ba, da = GI.elem_array(chrom, elem_id, col, False)
	fa, ra = GI.elem_array(chrom, elem_id, col, True)
	assert(not da)
	logger.debug("%.3f seconds"%(time()-start))
	return fa, ra, ba
sd = 3
def _calc_stats(A):
	tp = comparisons.tp(A)
	fp = comparisons.fp(A)
	tn = comparisons.tn(A)
	fn = comparisons.fn(A)
	sen = comparisons.sensitivity(A)
	spe = comparisons.specificity(A)
	pre = comparisons.precision(A)
	return tp, fp, tn, fn, np.round(sen,sd), np.round(spe,sd), np.round(pre,sd)
def _calc_stats_region(Ab, aB, AB):
	tp = AB
	fp = aB
	fn = Ab
	sen = tp/float(tp+fn) if tp+fn else np.nan
	pre = tp/float(tp+fp) if tp+fp else np.nan
	return tp, fp, fn, np.round(sen,sd), np.round(pre,sd)

def main():
	pass

if __name__ == "__main__":
	main()

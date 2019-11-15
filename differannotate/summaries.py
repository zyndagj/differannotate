#!/usr/bin/env python
#
###############################################################################
# Author: Greg Zynda
# Last Modified: 11/15/2019
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

import logging
from .constants import FORMAT
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.WARN, format=FORMAT)

from . import comparisons
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import numpy as np

def tabular_region(GI, p=95):
	chrom_set = GI.get_chrom_set()	# intersecting set chroms from all files
	max_chrom_len = max(map(len, chrom_set)+[len("Chrom")])
	max_elem_len = max(map(len, list(GI.element_dict)+list(GI.order_dict)+list(GI.sufam_dict)))
	max_name_len = max(map(len, list(GI.gff3_names)))
	#feature_set = # features only in reference
	non_te_elements = set(GI.element_dict) - GI.te_names
	te_elements = set(GI.element_dict) & GI.te_names
	for chrom in chrom_set:
		_print_table_region(GI, chrom, GI.element_dict, 1, max_chrom_len, max_elem_len, max_name_len, p)
		_print_table_region(GI, chrom, GI.order_dict, 2, max_chrom_len, max_elem_len, max_name_len, p)
		_print_table_region(GI, chrom, GI.sufam_dict, 3, max_chrom_len, max_elem_len, max_name_len, p)
def _print_table_region(GI, chrom, elem_list, col, mcl, mel, mnl, p=95):
	target = ("", "Element", "TE_Order", "TE_Superfamily")
	mel = max(map(len, target)+[mel])
	header = ("Chrom","S",target[col],"Sample","TP", "FP", "FN", "SENS", "PREC")
	template = "{:<{mcl}} {:^3} {:<{mel}} {:<{mn}} "+' '.join(["{:>5}"]*5)
	print(template.format(*header, mcl=mcl, mn=mnl, mel=mel))
	cname = GI.gff3_names[0]
	for elem in elem_list:
		eid = elem_list[elem]
		#print col, eid
		for sstr, sval in zip(('+/-','+','-'), (False, '+', '-')):
			#if not np.sum(A): continue
			for i, name in enumerate(GI.gff3_names):
				Ab, aB, AB = GI.calc_intersect_2(chrom, cname, name, eid, col, p, strand=sval)
				#if col == 2:
				#	sAb, saB, sAB = GI.calc_intersect_2(chrom, cname, name, eid, col, p, strand=sval, ret_set=True)
				#	print "%s-%s"%(cname.upper(), name), sAb
				#	print "%s-%s"%(cname, name.upper()), saB
				#	print "%s-%s"%(cname.upper(), name.upper()), sAB
				tp, fp, fn, sen, pre = _calc_stats_region(Ab, aB, AB)
				if not i:
					print(template.format(chrom, sstr, elem, name, tp,fp,fn,sen,pre, mcl=mcl, mn=mnl, mel=mel))
				else:
					print(template.format('', '', '', name, tp,fp,fn,sen,pre, mcl=mcl, mn=mnl, mel=mel))
			if len(GI.gff3_names) in (2,3):
				plt.figure()
				plt.title("%s %s %s"%(chrom, sstr, elem))
			if len(GI.gff3_names) == 2: # (Ab, aB, AB)
				n1, n2 = GI.gff3_names
				ret = GI.calc_intersect_2(chrom ,n1, n2, eid, col, p, strand=sval)
				if Ab or aB or AB:
					venn2(subsets=ret, set_labels=GI.gff3_names)
				else:
					strand = 'B' if sstr == '+/-' else sstr
					logger.warn("Empty plot for %s_%s_%s"%(chrom, strand, elem))
			elif len(GI.gff3_names) == 3: # (Abc, aBc, ABc, abC, AbC, aBC, ABC)
				n1, n2, n3 = GI.gff3_names
				ret = GI.calc_intersect_2(chrom ,n1, n2, n3, eid, col, p, strand=sval)
				if sum(ret):
					venn3(subsets=ret, set_labels=GI.gff3_names)
				else:
					strand = 'B' if sstr == '+/-' else sstr
					logger.warn("Empty plot for %s_%s_%s"%(chrom, strand, elem))
			if len(GI.gff3_names) in (2,3):
				plt.savefig("region_%s_%s_%s.png"%(chrom,'B' if sstr == '+/-' else sstr,elem))
				plt.close()
	print("")

def tabular(GI, strand=True):
	chrom_set = GI.get_chrom_set()	# intersecting set chroms from all files
	max_chrom_len = max(map(len, chrom_set)+[len("Chrom")])
	max_elem_len = max(map(len, list(GI.element_dict)+list(GI.order_dict)+list(GI.sufam_dict)))
	max_name_len = max(map(len, list(GI.gff3_names)))
	#feature_set = # features only in reference
	non_te_elements = set(GI.element_dict)-GI.te_names
	te_elements = set(GI.element_dict) & GI.te_names
	for chrom in chrom_set:
		_print_table(GI, chrom, GI.element_dict, 1, max_chrom_len, max_elem_len, max_name_len)
		_print_table(GI, chrom, GI.order_dict, 2, max_chrom_len, max_elem_len, max_name_len)
		_print_table(GI, chrom, GI.sufam_dict, 3, max_chrom_len, max_elem_len, max_name_len)
				
def _print_table(GI, chrom, elem_list, col, mcl, mel, mnl):
	target = ("", "Element", "TE_Order", "TE_Superfamily")
	mel = max(map(len, target)+[mel])
	header = ("Chrom","S",target[col],"Sample","TP", "FP", "TN", "FN", "SENS", "SPEC", "PREC")
	template = "{:<{mcl}} {:^3} {:<{mel}} {:<{mn}} "+' '.join(["{:>5}"]*7)
	print(template.format(*header, mcl=mcl, mn=mnl, mel=mel))
	for elem in elem_list:
		fa, ra, ba = _gen_arrays(GI, chrom, elem, col)
		for s, A in zip(('+/-','+','-'), (ba,fa,ra)):
			#if not np.sum(A): continue
			tp, fp, tn, fn, sen, spe, pre = _calc_stats(A)
			for i, name in enumerate(GI.gff3_names):
				if not i:
					print(template.format(chrom, s, elem, name, tp[i],fp[i],tn[i],fn[i],sen[i],spe[i],pre[i], mcl=mcl, mn=mnl, mel=mel))
				else:
					print(template.format('','','', name, tp[i],fp[i],tn[i],fn[i],sen[i],spe[i],pre[i], mcl=mcl, mn=mnl, mel=mel))
			if A.shape[0] in (2,3):
				plt.figure()
				plt.title("%s %s %s"%(chrom, s, elem))
			if A.shape[0] == 2: # (Ab, aB, AB)
				Ab = _venn2_helper(A, 1, 0)
				aB = _venn2_helper(A, 0, 1)
				AB = _venn2_helper(A, 1, 1)
				assert(Ab == fn[1])
				assert(aB == fp[1])
				assert(AB == tp[1])
				if A.sum():
					venn2(subsets=(fn[1], fp[1], tp[1]), set_labels=GI.gff3_names)
				else:
					strand = 'B' if s == '+/-' else s
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
					strand = 'B' if s == '+/-' else s
					logger.warn("Empty plot for %s_%s_%s"%(chrom, strand, elem))
			if A.shape[0] in (2,3):
				plt.savefig("base_%s_%s_%s.png"%(chrom,'B' if s == '+/-' else s,elem))
				plt.close()
	print("")

def _venn3_helper(array, rv0=1, rv1=0, rv2=0):
	m0 = array[0,:] == rv0
	m1 = array[1,:] == rv1
	m2 = array[2,:] == rv2
	assert(len(m0) == len(m1) and len(m1) == len(m2))
	m0a1 = np.logical_and(m0, m1)
	m0a1a2 = np.logical_and(m0a1, m2)
	return sum(m0a1a2)
def _venn2_helper(array, rv0=1, rv1=0):
	m0 = array[0,:] == rv0
	m1 = array[1,:] == rv1
	assert(len(m0) == len(m1))
	m0a1 = np.logical_and(m0, m1)
	return sum(m0a1)

def _gen_arrays(GI, chrom, elem, col):
	if col == 1:
		elem_id = GI.element_dict[elem]
	elif col == 2:
		elem_id = GI.order_dict[elem]
	elif col == 3:
		elem_id = GI.sufam_dict[elem]
	ba, da = GI.elem_array(chrom, elem_id, col, False)
	fa, ra = GI.elem_array(chrom, elem_id, col, True)
	assert(not da)
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

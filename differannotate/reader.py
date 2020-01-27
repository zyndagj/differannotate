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

import logging, re, os, sys
import numpy as np
import multiprocessing as mp
from functools import partial
from pysam import FastaFile
from collections import Counter
from collections import defaultdict as dd
from differannotate.constants import FORMAT, BaseIndex

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.WARN, format=FORMAT)

from differannotate.datastructures import *
from differannotate.comparisons import overlap_r, _overlap_r_tup

class gff3_interval:
	def __init__(self, gff3, name='control', fasta=None, include_chrom=False, force=False, \
			chrom_names=['chromosome','contig','supercontig'], \
			te_names=['transposable_element', 'transposable_element_gene', 'transposon_fragment']):
		self._order_re = re.compile('[Oo]rder=(?P<order>[^;/]+)')
		self._sufam_re = re.compile('[Ss]uperfamily=(?P<sufam>[^;]+)')
		self.element_dict = dict_index()
		self.order_dict = dict_index()
		self.sufam_dict = dict_index()
		self.strand_dict = {'+':0, '-':1, 0:'+', 1:'-'}
		self.chrom_names = set(chrom_names)
		self.te_names = set(te_names)
		self.include_chrom = include_chrom
		self.chrom_lens = None
		self.FA = False
		self.pool = False
		if fasta and os.path.exists(fasta+'.fai'):
			self.FA = FastaFile(fasta)
			self.chrom_lens = self._parse_fai(fasta+'.fai')
			self.pool = mp.Pool(4, worker_init, (fasta,))
		# create the initial interval tree
		self.gff3_trees = {name:self._2tree(gff3)}
		self.gff3_names = [name]
	def __del__(self):
		if self.pool:
			self.pool.close()
			self.pool.join()
	def _parse_fai(self, fai_file):
		'''
		Parses a fa.fai into a python dictionary
		
		# Parameters
		fai_file (str): Path to fai file
		
		# Returns
		dict: {chr:size, chr:size... }
		'''
		with open(fai_file,'r') as FAI:
			return dict(map(lambda y: (y[0], int(y[1])), map(lambda y: y.split('\t'), FAI.readlines())))
	def add_gff3(self, gff3, name):
		self.gff3_trees[name] = self._2tree(gff3, control=False)
		self.gff3_names.append(name)
	def _2tree(self, gff3, control=True):
		#Chr1    TAIR10  transposable_element_gene       433031  433819  .       -       .       ID=AT1G02228;Note=transposable_element_gene;Name=AT1G02228;Derives_from=AT1TE01405
		exclude = set(self.chrom_names) if self.include_chrom else set([])
		interval_tree = dd(iterit)
		with open(gff3,'r') as IF:
			for line in filter(lambda x: x[0] != "#", IF):
				tmp = line.rstrip('\n').split('\t')
				try:
					chrom, strand, element, attributes = tmp[0], tmp[6], tmp[2].lower(), tmp[8]
				except:
					print(tmp)
					sys.exit(1)
				if element not in exclude and (control or element in self.element_dict):
					element_id = self.element_dict[element]
					strand_id = self.strand_dict[strand]
					start, end = map(int, tmp[3:5])
					if element in self.te_names:
						te_order, te_sufam = self._extract_order_sufam(attributes)
						try:
							te_order_id = self.order_dict[te_order]
							te_sufam_id = self.sufam_dict[te_sufam]
						except:
							te_order_id, te_sufam_id = 0,0
						# Interval tree is not inclusive on the upper limit
						interval_tree[chrom].add(start-1, end, (strand_id, element_id, te_order_id, te_sufam_id))
					else:
						interval_tree[chrom].add(start-1, end, (strand_id, element_id))
		return interval_tree
	def _extract_order_sufam(self, attribute_string):
		order_match = self._order_re.search(attribute_string)
		sufam_match = self._sufam_re.search(attribute_string)
		order_str = order_match.group('order') if order_match else ''
		sufam_str = sufam_match.group('sufam') if sufam_match else ''
		return (order_str, sufam_str)
	def _get_max(self, chrom):
		if self.chrom_lens:
			return self.chrom_lens[chrom]
		return max([iit[chrom].max for iit in self.gff3_trees.values()])
	def get_chrom_set(self):
		ret_set = set(self.gff3_trees[self.gff3_names[0]])
		for n in self.gff3_names[1:]:
			ret_set &= set(self.gff3_trees[n])
		return ret_set
	def elem_array(self, chrom, eid, col=1, strand=True):
		'''
		Creates a binary numpy array to represent the presence of
		a specific element id

		# Parameters
		chrom (str): Target chromosome
		eid (int): Element id
		col (int): Can target {1:element, 2:te_order, 3:te_sufam}
		strand (bool): Return stranded results

		# Returns
		np.ndarray: Forward (or both strands)
		np.ndarray: Reverse strand
		'''
		max_size = self._get_max(chrom)
		num_rows = len(self.gff3_names)
		p_array = np.zeros((num_rows, max_size), dtype=np.bool)
		n_array = np.zeros((num_rows, max_size), dtype=np.bool)
		for i,name in enumerate(self.gff3_names):
			iit = self.gff3_trees[name][chrom]
			for interval in iit.search(0,max_size):
				D = interval.data
				if len(D) < col+1 or D[col] != eid:
					continue
				strand_id = D[0]
				s,e = interval.start, interval.end
				if strand_id == 0 or not strand:
					p_array[i, s:e] = 1
				elif strand_id == 1 and strand:
					n_array[i, s:e] = 1
		if strand:
			return p_array, n_array
		else:
			return p_array, []
	def calc_intersect_2(self, chrom, name1, name2, elem, col, p=95, strand=False, ret_set=False):
		eid = self._get_eid(elem)
		# (Ab, aB, AB)
		for n in (name1, name2): assert(chrom in self.gff3_trees[n])
		n1_tree = self.gff3_trees[name1][chrom]
		n2_tree = self.gff3_trees[name2][chrom]
		# Leftover
		n1_set = n1_tree.to_set(eid, col, strand)	#Ab
		n2_set = n2_tree.to_set(eid, col, strand)	#aB
		# Used
		n1_int_set, n2_int_set = set(), set()	#AB
		for interval_tup in map(interval2tuple, n1_tree.iifilter(eid, col, strand)):
			if interval_tup in n1_int_set:
				continue
			n1_s, n1_e = interval_tup[0], interval_tup[1]
			for n2int_tup in map(interval2tuple, n2_tree.searchfilter(n1_s, n1_e, eid, col, strand)):
				if n2int_tup in n2_set and _overlap_r_tup(interval_tup, n2int_tup, p):
					n1_int_set.add(interval_tup)
					n1_set.remove(interval_tup)
					n2_int_set.add(n2int_tup)
					n2_set.remove(n2int_tup)
					break
		assert len(n1_int_set) == len(n2_int_set)
		if ret_set:
			return n1_set, n2_set, n1_int_set
		return len(n1_set), len(n2_set), len(n1_int_set)
	def calc_intersect_3(self, chrom, name1, name2, name3, elem, col, p=95, strand=False, ret_set=False):
		# (Abc, aBc, ABc, abC, AbC, aBC, ABC)
		eid = self._get_eid(elem)
		for n in (name1, name2, name3): assert(chrom in self.gff3_trees[n])
		n1_set = self.gff3_trees[name1][chrom].to_set(eid, col, strand)
		func = self.calc_intersect_2
		Ab12, aB12, AB12 = func(chrom, name1, name2, elem, col, p, strand, ret_set=True)
		#print "Calc3",elem,col,p
		#print 'Ab12', Ab12, 'aB12', aB12, 'AB12', AB12
		Ab13, aB13, AB13 = func(chrom, name1, name3, elem, col, p, strand, ret_set=True)
		#print 'Ab13', Ab13, 'aB13', aB13, 'AB13', AB13
		Ab23, aB23, AB23 = func(chrom, name2, name3, elem, col, p, strand, ret_set=True)
		#print 'Ab23', Ab23, 'aB23', aB23, 'AB23', AB23
		Abc123 = Ab12 & Ab13
		aBc123 = _set_int(aB12, Ab23, p)
		ABc123 = AB12 - (AB13 | _set_mutate(n1_set, AB23,p))
		abC123 = _set_int(aB13, aB23, p)
		AbC123 = AB13 - (AB12 | _set_mutate(n1_set, AB23,p))
		aBC123 = _set_mutate(n1_set, AB23, p) - (AB12 | AB13)
		ABC123 = _set_int(_set_int(AB12, AB13, p), AB23, p)
		ret = (Abc123 , aBc123 , ABc123 , abC123 , AbC123 , aBC123 , ABC123)
		#print 'Abc123',Abc123 , 'aBc123',aBc123 , 'ABc123',ABc123 , 'abC123',abC123 , 'AbC123',AbC123 , 'aBC123',aBC123 , 'ABC123',ABC123
		if ret_set:
			return ret
		return tuple(map(len, ret))
	def get_length_array(self, chrom, name, elem, col, strand=False):
		eid = self._get_eid(elem)
		return map(_tuple_size, self.gff3_trees[name][chrom].to_set(eid, col, strand))
	def region_analysis(self, p=95):
		pass
		# TODO
	def fetch(self, chrom, start, end):
		outA = np.zeros((end-start, len(gff3_f2i)+2), dtype=np.uint8)
		for interval in self.interval_tree[chrom].search(start,end):
			s = max(interval.start, start)-start
			e = min(interval.end, end)-start
			element_id, te_order_id, te_sufam_id = interval.data
			outA[s:e,element_id] = 1
			outA[s:e,-2] = te_order_id
			outA[s:e,-1] = te_sufam_id
		return outA
	def _get_eid(self, elem):
		try:
			eid = int(elem)
		except ValueError:
			eid = self.element_dict[elem]
		return eid
	def get_proportion_arrays(self, chrom, name, elem, col, strand=False):
		eid = self._get_eid(elem)
		interval_set = self.gff3_trees[name][chrom].to_set(eid, col, strand)
		if not interval_set: return [[],[],[],[]]
		if self.pool:
			partial_wtp = partial(worker_tuple_proportion, chrom=chrom)
			proportion_arrays = zip(*self.pool.imap(partial_wtp, interval_set, chunksize=10))
		else:
			proportion_arrays = zip(*map(lambda x: self._tuple_proportion(chrom, x), interval_set))
		assert(len(interval_set) == len(proportion_arrays[0]))
		return proportion_arrays
	def _tuple_proportion(self, chrom, interval_tuple):
		start, end = interval_tuple[0], interval_tuple[1]
		seq = self.FA.fetch(chrom, start, end).upper()
		assert(len(seq) == end - start)
		count_dict = Counter(seq)
		total = sum(count_dict.values())
		return tuple((float(count_dict[base])/total for base in ('A','T','G','C')))
	#A	G	C	T	Total	A	G	C	T	Total
	#2	6	5	3		0.125	0.375	0.3125	0.1875	1
	#30	10	11	26		0.38961	0.12987	0.14285	0.33766	1
	#0.3440	0.17204	0.17204	0.31182	1	0.25730	0.25243	0.22767	0.26258	1 Proportion Calculation

def worker_init(fasta):
	global FA
	FA = FastaFile(fasta)
def worker_tuple_proportion(interval_tuple, chrom):
	start, end = interval_tuple[0], interval_tuple[1]
	seq = FA.fetch(chrom, start, end).upper()
	assert(len(seq) == end - start)
	count_dict = Counter(seq)
	total = sum(count_dict.values())
	return tuple((float(count_dict[base])/total for base in ('A','T','G','C')))

def _tuple_size(interval_tuple):
        return interval_tuple[1] - interval_tuple[0]

def _map_size(interval_set):
        return map(_tuple_size, interval_set)

def _set_int(prior, second, p=95):
	base = prior & second
	outBase = prior & second
	for tupP in prior - base:
		for tupS in second - outBase:
			if _overlap_r_tup(tupP, tupS, p):
				outBase.add(tupP)
				break
	return outBase
def _set_mutate(prior, second, p=95):
	base = prior & second
	outBase = prior & second
	for tupS in second - base:
		added = False
		for tupP in prior - outBase:
			if _overlap_r_tup(tupP, tupS, p):
				outBase.add(tupP)
				added = True
				break
		if not added:
			outBase.add(tupS)
	return outBase

if __name__ == "__main__":
	import doctest
	doctest.testmod()

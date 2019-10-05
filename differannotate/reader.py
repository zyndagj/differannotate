#!/usr/bin/env python
#
###############################################################################
# Author: Greg Zynda
# Last Modified: 10/04/2019
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

import logging, re, os
from quicksect import IntervalTree
import numpy as np
from collections import defaultdict as dd
from .constants import FORMAT

logging.basicConfig(level=logging.WARN, format=FORMAT)

class dict_index(dict):
	'''
	A modified dictionary class meant to store and increment unique values
	IDs as values are retrieved from keys.

	# Usage
	>>> DI = dict_index()
	>>> DI['cat']
	0
	>>> DI['bear']
	1
	>>> DI['cat']
	0
	>>> DI['cat'] = 10
	>>> DI['cat']
	0
	>>> DI.getkey(0)
	'cat'
	>>> DI.getkey(1)
	'bear'
	>>> DI.getkey(2)
	Traceback (most recent call last):
	...
	KeyError: 2
	>>> DI.getkey('dog')
	Traceback (most recent call last):
	...
	TypeError: dog
	'''
	def __init__(self):
		super(dict_index,self).__init__()
		self.cur = 0
	def __getitem__(self, key):
		try:
			return super(dict_index,self).__getitem__(key)
		except:
			super(dict_index,self).__setitem__(key, self.cur)
			self.cur += 1
			return super(dict_index,self).__getitem__(key)
	def __setitem__(self, key, value):
		pass
	def getkey(self, val):
		'''
		# Parameters
		val (int): Should be < len(dict_index)

		# Raises
		TypeError: if val is not an integer
		KeyError: if val does not exist as a value in the dict_index
		'''
		if not isinstance(val, int):
			raise TypeError(val)
		if val >= self.cur:
			raise KeyError(val)
		keys = super(dict_index,self).keys()
		vals = super(dict_index,self).values()
		return keys[vals.index(val)]

class iterit(IntervalTree):
	def __init__(self):
		super(iterit,self).__init__()
		self.min = None
		self.max = None
	def add(self, start, end, other=None):
		if self.min == None:
			self.min = start
			self.max = end
		else:
			if start < self.min:
				self.min = start
			if end > self.max:
				self.max = end
		super(iterit,self).add(start, end, other)
	def iterintervals(self):
		return super(iterit,self).search(self.min, self.max)

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
		if fasta and os.path.exists(fasta+'.fai'):
			self.chrom_lens = self._parse_fai(fasta+'.fai')
		# create the initial interval tree
		self.gff3_trees = {name:self._2tree(gff3)}
		self.gff3_names = [name]
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
		self.gff3_trees[name] = self._2tree(gff3)
		self.gff3_names.append(name)
	def _2tree(self, gff3):
		#Chr1    TAIR10  transposable_element_gene       433031  433819  .       -       .       ID=AT1G02228;Note=transposable_element_gene;Name=AT1G02228;Derives_from=AT1TE01405
		exclude = set(self.chrom_names) if self.include_chrom else set([])
		interval_tree = dd(iterit)
		with open(gff3,'r') as IF:
			for line in filter(lambda x: x[0] != "#", IF):
				tmp = line.rstrip('\n').split('\t')
				chrom, strand, element, attributes = tmp[0], tmp[6], tmp[2].lower(), tmp[8]
				if element not in exclude:
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
				if interval.data[col] != eid:
					continue
				strand_id = interval.data[0]
				s,e = interval.start, interval.end
				if strand_id == 0 or not strand:
					p_array[i, s:e] = 1
				elif strand_id == 1 and strand:
					n_array[i, s:e] = 1
		if strand:
			return p_array, n_array
		else:
			return p_array, []
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

if __name__ == "__main__":
	import doctest
	doctest.testmod()

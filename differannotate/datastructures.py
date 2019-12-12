#!/usr/bin/env python
#
###############################################################################
# Author: Greg Zynda
# Last Modified: 12/11/2019
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

from quicksect import IntervalTree
import logging
from differannotate.constants import FORMAT

logger = logging.getLogger(__name__)
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
		self.set_cache = {}
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
	def iifilter(self, eid, col, strand=False):
		'''
		>>> IT = iterit()
		>>> IT.add(0, 10, (0, 0))
		>>> IT.add(5, 15, (1, 1))
		>>> IT.add(10, 20, (1, 2))
		>>> ret = IT.iifilter(1,0)
		>>> len(ret)
		2
		>>> for i in map(interval2tuple, ret): print i
		(5, 15, 1, 1)
		(10, 20, 1, 2)
		'''
		assert(col >= 1)
		if _strand(strand):
			sid = _get_strand(strand)
			return list(filter(lambda x: len(x.data) > col and x.data[col] == eid and x.data[0] == sid, self.iterintervals()))
		else:
			return list(filter(lambda x: len(x.data) > col and x.data[col] == eid, self.iterintervals()))
	def searchfilter(self, start, end, eid, col, strand=False):
		assert(col >= 1)
		if _strand(strand):
			sid = _get_strand(strand)
			return list(filter(lambda x: len(x.data) > col and x.data[col] == eid and x.data[0] == sid, super(iterit,self).search(start, end)))
		else:
			return list(filter(lambda x: len(x.data) > col and x.data[col] == eid, super(iterit,self).search(start, end)))
	def to_set(self,eid=False, col=False, strand=False):
		cache_name = (eid, col, strand)
		if cache_name in self.set_cache:
			return self.set_cache[cache_name].copy()
		if eid or col or strand:
			ret = set(map(interval2tuple, self.iifilter(eid, col, strand)))
		else:
			ret = set(map(interval2tuple, self.iterintervals()))
		self.set_cache[cache_name] = ret
		return ret.copy()

def _strand(strand):
	return not isinstance(strand, bool)
strand_dict = {'+':0, '-':1, 0:'+', 1:'-'}
def _get_strand(strand):
	if isinstance(strand, int):
		return strand
	elif isinstance(strand, str):
		return strand_dict[strand]
	else:
		raise ValueError(strand)

def interval2tuple(interval):
	'''
	Converts an interval to a tuple

	# Usage
	>>> IT = iterit()
	>>> IT.add(0, 10, (0, 0))
	>>> IT.add(5, 15, (1, 1))
	>>> for i in map(interval2tuple, IT.iterintervals()): print i
	(0, 10, 0, 0)
	(5, 15, 1, 1)
	'''
	if interval.data:
		return (interval.start, interval.end)+tuple(interval.data)
	else:
		return (interval.start, interval.end)


if __name__ == "__main__":
	import doctest
	doctest.testmod()

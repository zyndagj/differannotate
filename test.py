#!/usr/bin/env python
#
###############################################################################
# Author: Greg Zynda
# Last Modified: 09/21/2019
###############################################################################
# BSD 3-Clause License
# 
# Copyright (c) 2017, Texas Advanced Computing Center
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

class nd(dict):
	def __init__(self):
		super(nd,self).__init__()
		self.cur = 0
	def __getitem__(self, key):
		try:
			return super(nd,self).__getitem__(key)
		except:
			super(nd,self).__setitem__(key, self.cur)
			self.cur += 1
			return super(nd,self).__getitem__(key)

class ititer(IntervalTree):
	def __init__(self):
		super(ititer,self).__init__()
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
		super(ititer,self).add(start, end, other)
	def iterintervals(self):
		return super(ititer,self).search(self.min, self.max)

def main():
	ND = nd()
	print(ND['cats'])
	print(ND['bears'])
	print(ND['cats'])
	print(list(ND.keys()))
	print(list(ND.values()))
	
	IT = ititer()
	IT.add(2,5)
	IT.add(5,10)
	print(IT.search(0,5))
	print(IT.min)
	print(IT.max)
	print(IT.iterintervals())

if __name__ == "__main__":
	main()


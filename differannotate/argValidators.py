#!/usr/bin/env python
#
###############################################################################
# Author: Greg Zynda
# Last Modified: 11/16/2019
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

import argparse, os

class argChecker:
	'''
	Checks to ensure an argument is contained in a whitelist
	
	Arguments
	-----------
	options -- iterator of valid options
	afterValid -- argument type

	Returns
	-----------
	matching option

	Usage
	-----------
	parser.add_argument("-t", metavar='STR', help=\
		"method (value|auto|percent) (Default: %(default)s)", default="auto",\
		type=argChecker(('value','auto','percent'),'replication method').check)
	'''
	def __init__(self, options, afterValid):
		self.options = options
		self.av = afterValid
	def check(self, x):
		if x in self.options:
			return x
		else:
			raise argparse.ArgumentTypeError("%s not a valid %s"%(x, self.av))
class fileCheck:
	'''
	Checks to make sure a file given as an argument is of the correct type
	
	Returns
	-----------
	input -- when correct

	Usage
	-----------
	fCheck = fileCheck()
	parser.add_argument('-R', metavar='FASTA', help='Reference for alignment', required=True, type=fCheck.fasta)
	'''
	def _check(self, file, exts):
		ext = os.path.splitext(file)[1][1:]
		fName = os.path.split(file)[1]
		if not ext in exts:
			raise argparse.ArgumentTypeError("%s not a %s"%(fName, exts[0]))
		if not os.path.exists(file):
			raise argparse.ArgumentTypeError("%s does not exist"%(file))
	def fastq(self, file):
		self._check(file, ['fastq','fq'])
		return file
	def fasta(self, file):
		self._check(file, ['fasta','fa'])
		return file
	def csv(self, file):
		self._check(file, ['csv'])
		return file
	def gff3(self, file):
		self._check(file, ['gff3','gff','gtf'])
		return file

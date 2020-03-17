#!/usr/bin/env python
#
###############################################################################
# Author: Greg Zynda
# Last Modified: 03/16/2020
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
import numpy as np
from quicksect import Interval
from time import time
from .constants import FORMAT

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.WARN, format=FORMAT)

from differannotate.datastructures import interval2tuple

def _super(array, control_row=0, control_target=1, treat_target=1):
	start = time()
	mask = array[control_row,:] == control_target
	ret = np.sum(array[:,mask] == treat_target, axis=1)
	logger.debug("%.4f seconds"%(time()-start))
	return ret
def tp(array, control_row=0):
	return _super(array, control_row, 1, 1)
def fp(array, control_row=0):
	return _super(array, control_row, 0, 1)
def tn(array, control_row=0):
	return _super(array, control_row, 0, 0)
def fn(array, control_row=0):
	return _super(array, control_row, 1, 0)
def _all(array, control_row=0):
	tpv = tp(array, control_row)
	fpv = fp(array, control_row)
	tnv = tn(array, control_row)
	fnv = fn(array, control_row)
	return tpv, fpv, tnv, fnv
def sensitivity(array, control_row=0):
	tpv, fpv, tnv, fnv = _all(array, control_row)
	return tpv/(tpv+fnv).astype(np.float)
def specificity(array, control_row=0):
	tpv, fpv, tnv, fnv = _all(array, control_row)
	return tnv/(tnv+fpv).astype(np.float)
def precision(array, control_row=0):
	tpv, fpv, tnv, fnv = _all(array, control_row)
	return tpv/(tpv+fpv).astype(np.float)

def _overlap_b(A, B):
	'''
	Calculates overlap in bases of A and B. Not inclusive
	>>> _overlap_b(Interval(0,10), Interval(5,15)
	5
	>>> _overlap_b(Interval(0,10), Interval(10,15)
	0
	'''
	aS, aE = _get_se(A)
	bS, bE = _get_se(B)
	return max(0, min(aE, bE) - max(aS, bS))
def _get_se(obj):
	if isinstance(obj, tuple):
		return obj[0], obj[1]
	elif isinstance(obj, Interval):
		return obj.start, obj.end
	else:
		raise TypeError(obj)

def _overlap_p(A, B):
	'''
	Calculates overlap percentage of A by B
	>>> _overlap_b(Interval(0,10), Interval(5,15)
	50.0
	>>> _overlap_b(Interval(0,10), Interval(10,15)
	0.0
	'''
	bases_overlap = _overlap_b(A,B)
	if isinstance(A, tuple):
		size_A = A[1] - A[0]
	elif isinstance(A, Interval):
		size_A = A.end - A.start
	else:
		raise TypeError(A)
	#logger.info("TEST")
	return float(bases_overlap)/float(size_A)*100.0
def _overlap_p_tup(A, B):
	'''
	Calculates overlap percentage of A by B
	'''
	bases_overlap = max(0, min(A[1], B[1]) - max(A[0], B[0]))
	size_A = A[1] - A[0]
	return float(bases_overlap)/float(size_A)*100.0

def overlap(A, B, overlap_p=95):
	'''
	>= 95% of A is covered by B
	'''
	if overlap_p > 100:
		logger.warn("Input overlap percentage > 100 [%i] falling back to 100"%(overlap_p))
		overlap_p = 100
	elif overlap_p < 1:
		logger.error("Input overlap should be expressed as a percentage and not a fraction")
		raise ValueError(overlap_p)
	dA, dB = map(_decode, [A,B])
	return _overlap_p(dA, dB) >= overlap_p
def _overlap_tup(A, B, overlap_p=95):
	return _overlap_p_tup(A, B) >= overlap_p
	
def overlap_r(A, B, overlap_p=95):
	'''
	>= 95% reciprocal overlap between A and B
	'''
	if overlap_p > 100:
		logger.warn("Input overlap percentage > 100 [%i] falling back to 100"%(overlap_p))
		overlap_p = 100
	elif overlap_p < 1:
		logger.error("Input overlap should be expressed as a percentage and not a fraction")
		raise ValueError(overlap_p)
	dA, dB = map(_decode, [A,B])
	oab = _overlap_tup(dA, dB, overlap_p)
	oba = _overlap_tup(dB, dA, overlap_p)
	return oab and oba
def _overlap_r_tup_val(A, B):
	bases_overlap = max(0, min(A[1], B[1]) - max(A[0], B[0]))
	if not bases_overlap: return False
	size_A = A[1] - A[0]
	size_B = B[1] - B[0]
	oab = float(bases_overlap)/float(size_A)*100.0
	oba = float(bases_overlap)/float(size_B)*100.0
	return min(oab, oba)
def _overlap_r_tup(A, B, overlap_p=95):
	'''
	>= 95% reciprocal overlap between A and B
	'''
	return  _overlap_r_tup_val(A, B) >= overlap_p
def _decode(obj):
	if isinstance(obj, tuple):
		return obj
	elif isinstance(obj, Interval):
		return interval2tuple(obj)
	else:
		raise TypeError(obj)

if __name__ == "__main__":
	import doctest
	doctest.testmod()

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

import logging
import numpy as np
from .constants import FORMAT

logging.basicConfig(level=logging.WARN, format=FORMAT)

def _super(array, control_row=0, control_target=1, treat_target=1):
	mask = array[control_row,:] == control_target
	return np.sum(array[:,mask] == treat_target, axis=1)
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

if __name__ == "__main__":
	import doctest
	doctest.testmod()


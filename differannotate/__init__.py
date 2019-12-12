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

import argparse, sys, logging, os
from differannotate.constants import FORMAT
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format=FORMAT)
from differannotate.argValidators import fileCheck, argChecker
from differannotate import reader, summaries

def main():
	fCheck = fileCheck() #class for checking parameters
	parser = argparse.ArgumentParser(description="A tool for comparing GFF3 annotations")
	parser.add_argument('-C', '--control', metavar='GFF3', help='Control GFF3', required=True, type=fCheck.gff3)
	parser.add_argument('-R', '--reference', metavar='FASTA', \
		help='Control Reference (accurate base metrics)', type=fCheck.fasta)
	parser.add_argument('--cname', metavar='STR', help='Name of control GFF3', default='control', type=str)
	parser.add_argument('-T', '--treat', metavar='GFF3', \
		help='Space separated list of GFF3 files for comparison against the control', \
		type=fCheck.gff3, required=True, nargs='+')
	parser.add_argument('-N', '--names', metavar='STR', \
		help='Space separated list of names for treatment GFF3 files (must be same order)', \
		type=str, required=True, nargs='+')
	parser.add_argument('-p', '--percent', metavar='INT', \
		help='Reciprocal percent overlap threshold', type=int, default=90)
	parser.add_argument('--plot', action="store_true", help="Plot venn diagrams of results")
	parser.add_argument('-e', '--ext', metavar='EXT', \
		help='Figure extension [%(default)s]', default='png', \
		type=argChecker(('pdf','png','eps'),'figure extension').check)
	parser.add_argument('-v', '--verbose', action="store_true", help='Enable verbose logging')
	parser.add_argument('--temd', action="store_true", help='Analyze TE metadata')
	args = parser.parse_args()
	################################
	# Configure logging
	################################
	if args.verbose:
		logger.setLevel(logging.DEBUG)
		logger.debug("DEBUG logging enabled")
	else:
		logger.setLevel(logging.INFO)
	################################
	# Check arguments
	################################
	if len(args.names) != len(args.treat):
		logger.error("treat(%i) != names(%i)"%(len(args.treat), len(args.names)))
		raise ValueError
	################################
	# Create GFF3 intervals
	################################
	GI = reader.gff3_interval(args.control, fasta=args.reference)
	for f, n in zip(args.treat, args.names):
		GI.add_gff3(f, n)
	################################
	# Generate results
	################################
	fig_ext = args.ext if args.plot else False
	if args.reference:
		logger.info("Basepair resolution results")
		summaries.tabular(GI, fig_ext=fig_ext, temd=args.temd)
	logger.info("Interval results")
	summaries.tabular_region(GI, p=args.percent, fig_ext=fig_ext, temd=args.temd)
	logger.info("Done")

if __name__ == "__main__":
	main()

#!/usr/bin/env python

import unittest, sys, os, logging
from glob import glob
from time import time
from shutil import rmtree
import numpy as np
try:
	from StringIO import StringIO
	from unittest.mock import patch
except:
	from io import StringIO
	from mock import patch
# Create root logger before main import
from differannotate.constants import FORMAT
logStream = StringIO()	# buffer for capturing log info
logging.basicConfig(stream=logStream, level=logging.DEBUG, format=FORMAT)

import differannotate
from differannotate import reader, comparisons

class TestReader(unittest.TestCase):
	def setUp(self):
		tpath = os.path.dirname(__file__)
		self.fa = os.path.join(tpath, 'test.fa')
		self.fai = os.path.join(tpath, 'test.fa.fai')
		self.gff3_1 = os.path.join(tpath, 'test_1.gff3')
		self.gff3_2 = os.path.join(tpath, 'test_2.gff3')
		self.gff3_max = {'Chr1':1500, 'Chr2':15}
		self.fa_max = {'Chr1':2000, 'Chr2':20}
		#
		self.elems_sum1 = {'transposon_fragment':(100*3,100*2,100*5), \
			'cds':(0,0,0), \
			'exon':(100,100,200), \
			'lnc_rna':(100,0,100), \
			'antisense_lncrna':(0,100,100), \
			'transposable_element':(100+500,100*2,100+100+500), \
			'protein':(0,100,100), \
			'gene':(100,100,200)}
		self.order_sum1 = {'LTR':(100,0,100), \
			'DNA':(0,100,100), \
			'RC':(500,100,500+100)}
		self.sufam_sum1 = {'Copia':(100,0,100), \
			'HAT':(0,100,100), \
			'Helitron':(500,100,500+100)}
		#
		self.elems_sum2 = {'transposon_fragment':(107+100*2,104+100,107+104+100*3), \
			'cds':(0,0,0), \
			'exon':(100,100,200), \
			'lnc_rna':(100,0,100), \
			'antisense_lncrna':(0,100,100), \
			'transposable_element':(94+500,96+100,94+96+500), \
			'protein':(0,0,0), \
			'gene':(100,100,200)}
		self.order_sum2 = {'LTR':(107,0,107), \
			'DNA':(0,100,100), \
			'RC':(500,104,500+104)}
		self.sufam_sum2 = {'Copia':(107,0,107), \
			'HAT':(0,100,100), \
			'Helitron':(500,104,500+104)}
	def tearDown(self): ## Runs after every test function ##
		logStream.truncate(0) # Wipe log
		#if os.path.exists('mse_tmp'): rmtree('mse_tmp')
	def test_dict_index(self):
		DI = reader.dict_index()
		self.assertEqual(DI['cat'], 0)
		self.assertEqual(DI['bear'], 1)
		self.assertEqual(DI['cat'], 0)
		DI['cat'] = 10
		self.assertEqual(DI['cat'], 0)
		self.assertEqual(DI.getkey(0), 'cat')
		self.assertEqual(DI.getkey(1), 'bear')
		with self.assertRaises(KeyError):
			DI.getkey(2)
		with self.assertRaises(TypeError):
			DI.getkey('dog')
	def test_iterit(self):
		IIT = reader.iterit()
		IIT.add(10, 100)
		self.assertEqual(IIT.max, 100)
		self.assertEqual(IIT.min, 10)
		IIT.add(4, 30)
		self.assertEqual(IIT.max, 100)
		self.assertEqual(IIT.min, 4)
		self.assertEqual(len(list(IIT.iterintervals())), 2)
	def test_gff3_1_nofa(self):
		GI = reader.gff3_interval(self.gff3_1)
		self._test_elem_sets(GI)
		self.assertEqual(GI.gff3_names, ['control'])
		for chrom in ('Chr1','Chr2'):
			self.assertEqual(GI._get_max(chrom), self.gff3_max[chrom])
		# Check sums
		for elem in GI.element_dict.keys():
			fa, ra, ba = self._gen_arrays(GI, 'Chr1', elem, 1)
			for i, a in enumerate((fa, ra, ba)):
				self.assertEqual(a.shape, (1, self.gff3_max['Chr1']))
				self.assertEqual(np.sum(a), self.elems_sum1[elem][i])
	def test_gff3_1_fa(self):
		GI = reader.gff3_interval(self.gff3_1, fasta=self.fa)
		self._test_elem_sets(GI)
		self.assertEqual(GI.gff3_names, ['control'])
		for chrom in ('Chr1','Chr2'):
			self.assertEqual(GI._get_max(chrom), self.fa_max[chrom])
		# Check sums
		for elem in GI.element_dict.keys():
			fa, ra, ba = self._gen_arrays(GI, 'Chr1', elem, 1)
			for i, a in enumerate((fa, ra, ba)):
				self.assertEqual(a.shape, (1, self.fa_max['Chr1']))
				self.assertEqual(np.sum(a), self.elems_sum1[elem][i])
	def test_gff3_12_nofa(self):
		GI = reader.gff3_interval(self.gff3_1)
		GI.add_gff3(self.gff3_2, 'treat')
		self._test_elem_sets(GI)
		self.assertEqual(GI.gff3_names, ['control','treat'])
		for chrom in ('Chr1','Chr2'):
			self.assertEqual(GI._get_max(chrom), self.gff3_max[chrom])
		# Check sums
		for elem in GI.element_dict.keys():
			fa, ra, ba = self._gen_arrays(GI, 'Chr1', elem, 1)
			for i, a in enumerate((fa, ra, ba)):
				self.assertEqual(a.shape, (2, self.gff3_max['Chr1']))
				self.assertEqual(np.sum(a[0,:]), self.elems_sum1[elem][i])
				self.assertEqual(np.sum(a[1,:]), self.elems_sum2[elem][i])
	def test_gff3_12_nofa(self):
		GI = reader.gff3_interval(self.gff3_1, fasta=self.fa)
		GI.add_gff3(self.gff3_2, 'treat')
		self._test_elem_sets(GI)
		self.assertEqual(GI.gff3_names, ['control','treat'])
		for chrom in ('Chr1','Chr2'):
			self.assertEqual(GI._get_max(chrom), self.fa_max[chrom])
		# Check sums
		for elem in GI.element_dict.keys():
			fa, ra, ba = self._gen_arrays(GI, 'Chr1', elem, 1)
			for i, a in enumerate((fa, ra, ba)):
				self.assertEqual(a.shape, (2, self.fa_max['Chr1']))
				self.assertEqual(np.sum(a[0,:]), self.elems_sum1[elem][i])
				self.assertEqual(np.sum(a[1,:]), self.elems_sum2[elem][i])
	def _test_elem_sets(self, GI):
		self.assertEqual(set(GI.element_dict.keys()), set(self.elems_sum1.keys()))
		self.assertEqual(set(GI.order_dict.keys()), set(self.order_sum1.keys()))
		self.assertEqual(set(GI.sufam_dict.keys()), set(self.sufam_sum1.keys()))
	def _gen_arrays(self, GI, chrom, elem, col):
		ba, da = GI.elem_array(chrom, GI.element_dict[elem], col, False)
		fa, ra = GI.elem_array(chrom, GI.element_dict[elem], col, True)
		self.assertFalse(da)
		return fa, ra, ba
class TestComparisons(unittest.TestCase):
	def setUp(self):
		tpath = os.path.dirname(__file__)
		self.fa = os.path.join(tpath, 'test.fa')
		self.fai = os.path.join(tpath, 'test.fa.fai')
		self.gff3_1 = os.path.join(tpath, 'test_1.gff3')
		self.gff3_2 = os.path.join(tpath, 'test_2.gff3')
		self.A = np.array([[1,1,0,0],[1,1,1,1],[0,1,1,0]])
	def tearDown(self): ## Runs after every test function ##
		logStream.truncate(0) # Wipe log
		#if os.path.exists('mse_tmp'): rmtree('mse_tmp')
	def test_tp(self):
		self.assertTrue(np.array_equal(comparisons.tp(self.A), [2,2,1]))
	def test_fp(self):
		self.assertTrue(np.array_equal(comparisons.fp(self.A), [0,2,1]))
	def test_tn(self):
		self.assertTrue(np.array_equal(comparisons.tn(self.A), [2,0,1]))
	def test_fn(self):
		self.assertTrue(np.array_equal(comparisons.fn(self.A), [0,0,1]))
	def test_sensitivity(self):
		self.assertTrue(np.array_equal(comparisons.sensitivity(self.A), \
			[2.0/(2+0),2.0/(2+0),1.0/(1+1)]))
	def test_specificity(self):
		self.assertTrue(np.array_equal(comparisons.specificity(self.A), \
			[2.0/(2+0),0,1.0/(1+1)]))
	def test_precision(self):
		self.assertTrue(np.array_equal(comparisons.precision(self.A), \
			[2.0/(2+0), 2.0/(2+2), 1.0/(1+1)]))
#	def test_train_cli_01(self):
#		if not self.test_model: return
#		testArgs = ['teamRNN', \
#			'-b', '-f']
#		with patch('sys.argv', testArgs):
#			teamRNN.main()
#		output = logStream.getvalue()
#		splitOut = output.split('\n')
#		self.assertTrue('Done' in splitOut[-2])
#		#for f in glob('test_cli/*'): print f
#		self.assertTrue(os.path.exists('test_cli/plain_s15x10_o68_1xbilstm60_merge-concat_statefulF_learn%s_drop0.h5'%(str(self.learning_rate))))
#		self.assertTrue(os.path.exists('test_cli/config.pkl'))

if __name__ == "__main__":
	unittest.main()

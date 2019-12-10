#!/usr/bin/env python

import unittest, sys, os, logging
from operator import itemgetter
from glob import glob
from time import time
from shutil import rmtree
try:
	print "Detected python 2"
	from StringIO import StringIO
	from mock import patch
except:
	print "Detected python 3"
	from io import StringIO
	from unittest.mock import patch
# Create root logger before main import
logStream = StringIO()	# buffer for capturing log info
logger = logging.getLogger(__name__)
FORMAT = "[%(levelname)s - %(filename)s:%(lineno)s - %(funcName)15s] %(message)s"
logging.basicConfig(stream=logStream, level=logging.INFO, format=FORMAT)
import numpy as np
from quicksect import Interval
import differannotate
from differannotate import reader, comparisons, summaries, datastructures

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
		DI = datastructures.dict_index()
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
		IIT = datastructures.iterit()
		IIT.add(10, 100)
		self.assertEqual(IIT.max, 100)
		self.assertEqual(IIT.min, 10)
		IIT.add(4, 30)
		self.assertEqual(IIT.max, 100)
		self.assertEqual(IIT.min, 4)
		self.assertEqual(len(list(IIT.iterintervals())), 2)
	def test_iifilter(self):
		IIT = datastructures.iterit()
		IIT.add(0, 10, (0, 0))
		IIT.add(5, 15, (1, 1))
		IIT.add(10, 20, (1, 2))
		IIT.add(10, 20, (1, 2, 1))
		ret = map(datastructures.interval2tuple, IIT.iifilter(1, 1))
		self.assertEqual(len(ret), 1)
		self.assertEqual(ret[0], (5, 15, 1, 1))
		ret = map(datastructures.interval2tuple, IIT.iifilter(2, 1, strand=0))
		self.assertEqual(len(ret), 0)
		ret = map(datastructures.interval2tuple, IIT.iifilter(2, 1, strand=1))
		self.assertEqual(len(ret), 2)
		ret = map(datastructures.interval2tuple, IIT.iifilter(1, 2, strand=1))
		self.assertEqual(len(ret), 1)
		self.assertEqual(ret[0], (10, 20, 1, 2, 1))
	def test_tuple_size(self):
		self.assertEqual(reader._tuple_size((0, 10, 0, 0)), 10)
		self.assertEqual(reader._tuple_size((5, 15, 1, 1)), 10)
		self.assertEqual(reader._tuple_size((10, 20, 1, 2)), 10)
		self.assertEqual(reader._tuple_size((10, 20, 1, 2, 1)), 10)
	def test_tuple_proportion(self):
		#TATTAGGCTGTGATGTGCTT
		#01234567890123456789
		GI = reader.gff3_interval(self.gff3_1, fasta=self.fa)
		self.assertEqual(GI._tuple_proportion('Chr2', (0, 10, 0, 0)), (0.2,0.4,0.3,0.1))
		self.assertEqual(GI._tuple_proportion('Chr2', (5, 15, 1, 1)), (0.1,0.3,0.5,0.1))
	def test_map_size(self):
		IIT = datastructures.iterit()
		IIT.add(0, 10, (0, 0))
		IIT.add(5, 15, (1, 1))
		IIT.add(10, 20, (1, 2))
		IIT.add(10, 20, (1, 2, 1))
		self.assertEqual(reader._map_size(IIT.to_set(1,1)), [10])
		self.assertEqual(reader._map_size(IIT.to_set(2,1,strand=0)), [])
		self.assertEqual(reader._map_size(IIT.to_set(2,1,strand=1)), [10,10])
		self.assertEqual(reader._map_size(IIT.to_set(1,2,strand=1)), [10])
	def test_get_proportion_arrays(self):
		#TATTAGGCTGTGATGTGCTT
		#01234567890123456789
		#  ----- ------
		GI = reader.gff3_interval(self.gff3_1, fasta=self.fa)
		prop_array = GI.get_proportion_arrays('Chr2', 'control', 'exon', 1, strand=False)
		self.assertEqual(sorted(prop_array[0]), sorted([float(1)/5, float(1)/6])) #A
		self.assertEqual(sorted(prop_array[1]), sorted([float(2)/5, float(3)/6])) #T
		self.assertEqual(sorted(prop_array[2]), sorted([float(2)/5, float(2)/6])) #G
		self.assertEqual(sorted(prop_array[3]), sorted([float(0)/5, float(0)/6])) #C
	def test_searchfilter(self):
		IIT = datastructures.iterit()
		IIT.add(0, 10, (0, 0))
		IIT.add(5, 15, (1, 1))
		IIT.add(10, 20, (2, 1))
		IIT.add(10, 20, (2, 1, 1))
		ret = map(reader.interval2tuple, IIT.searchfilter(6,6,1,1))
		self.assertEqual(len(ret), 1)
		self.assertEqual(ret[0], (5, 15, 1, 1))
		ret = map(reader.interval2tuple, IIT.searchfilter(10,11,1,2))
		self.assertEqual(len(ret), 1)
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
	def test_gff3_12_fa(self):
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
	def test_gff3_12_fa_region_2_gene(self):
		GI = reader.gff3_interval(self.gff3_1, fasta=self.fa)
		GI.add_gff3(self.gff3_2, 'treat')
		ret = GI.calc_intersect_2('Chr1',  'control', 'treat', GI.element_dict['gene'], 1, p=99)
		self.assertEquals(ret, (1,1,1))
	def test_gff3_12_fa_region_2_gene_strand(self):
		GI = reader.gff3_interval(self.gff3_1, fasta=self.fa)
		GI.add_gff3(self.gff3_2, 'treat')
		self.assertEquals(GI.calc_intersect_2('Chr1',  'control', 'treat', 'gene', 1, p=99, strand='+'), (0,0,1))
		self.assertEquals(GI.calc_intersect_2('Chr1',  'control', 'treat', 'gene', 1, p=99, strand='-'), (1,1,0))
		self.assertEquals(GI.calc_intersect_2('Chr1',  'control', 'treat', 'gene', 1, p=99, strand=0), (0,0,1))
		self.assertEquals(GI.calc_intersect_2('Chr1',  'control', 'treat', 'gene', 1, p=99, strand=1), (1,1,0))
	def test_gff3_12_fa_region_2_error(self):
		GI = reader.gff3_interval(self.gff3_1, fasta=self.fa)
		GI.add_gff3(self.gff3_2, 'treat')
		with self.assertRaises(AssertionError):
			GI.calc_intersect_2('Chr2',  'control', 'treat', GI.element_dict['gene'], 1, p=99)
	def test_gff3_12_fa_region_2(self):
		GI = reader.gff3_interval(self.gff3_1, fasta=self.fa)
		GI.add_gff3(self.gff3_2, 'treat')
		func = GI.calc_intersect_2
		self.assertEquals(func('Chr1',  'control', 'treat', 'gene', 1, p=99), (1,1,1))
		self.assertEquals(func('Chr1',  'control', 'treat', 'cds', 1, p=50), (0,0,0))
		self.assertEquals(func('Chr1',  'control', 'treat', 'lnc_rna', 1, p=99), (0,0,1))
		self.assertEquals(func('Chr1',  'control', 'treat', 'antisense_lncrna', 1, p=99), (0,0,1))
		self.assertEquals(func('Chr1',  'control', 'treat', 'exon', 1, p=97), (1,1,1))
		self.assertEquals(func('Chr1',  'control', 'treat', 'exon', 1, p=96), (0,0,2))
		self.assertEquals(func('Chr1',  'control', 'treat', 'protein', 1, p=50), (1,0,0))
	def test_gff3_122_fa_region_3(self):
		GI = reader.gff3_interval(self.gff3_1, fasta=self.fa)
		GI.add_gff3(self.gff3_2, 'treat1')
		GI.add_gff3(self.gff3_2, 'treat2')
		func = GI.calc_intersect_3
		args = ('Chr1', 'control', 'treat1', 'treat2')
		self.assertEquals(func(*args, elem='gene', col=1, p=99), (1,0,0,0,0,1,1))
		self.assertEquals(func('Chr1',  'control','treat1','treat2', 'cds', 1, p=50), \
			(0,0,0,0,0,0,0))
		self.assertEquals(func('Chr1',  'control','treat1','treat2', 'lnc_rna', 1, p=99), \
			(0,0,0,0,0,0,1))
		self.assertEquals(func('Chr1',  'control','treat1','treat2', 'antisense_lncrna', 1, p=99), \
			(0,0,0,0,0,0,1))
		self.assertEquals(func('Chr1',  'control','treat1','treat2', 'exon', 1, p=97), \
			(1,0,0,0,0,1,1))
		self.assertEquals(func('Chr1',  'control','treat1','treat2', 'exon', 1, p=96), \
			(0,0,0,0,0,0,2))
		self.assertEquals(func('Chr1',  'control','treat1','treat2', 'protein', 1, p=50), \
			(1,0,0,0,0,0,0))
	def _test_elem_sets(self, GI):
		self.assertEqual(set(GI.element_dict.keys()), set(self.elems_sum1.keys()))
		self.assertEqual(set(GI.order_dict.keys()), set(self.order_sum1.keys()))
		self.assertEqual(set(GI.sufam_dict.keys()), set(self.sufam_sum1.keys()))
	def _gen_arrays(self, GI, chrom, elem, col):
		ba, da = GI.elem_array(chrom, GI.element_dict[elem], col, False)
		fa, ra = GI.elem_array(chrom, GI.element_dict[elem], col, True)
		self.assertFalse(da)
		return fa, ra, ba
	def test_gff3_12_get_chrom_set(self):
		GI = reader.gff3_interval(self.gff3_1)
		self.assertEqual(GI.get_chrom_set(), set(('Chr1','Chr2')))
		GI.add_gff3(self.gff3_2, 'treat')
		self.assertEqual(GI.get_chrom_set(), set(('Chr1',)))
class TestComparisons(unittest.TestCase):
	def setUp(self):
		tpath = os.path.dirname(__file__)
		self.fa = os.path.join(tpath, 'test.fa')
		self.fai = os.path.join(tpath, 'test.fa.fai')
		self.gff3_1 = os.path.join(tpath, 'test_1.gff3')
		self.gff3_2 = os.path.join(tpath, 'test_2.gff3')
		self.A = np.array([[1,1,0,0],[1,1,1,1],[0,1,1,0]])
		self.i5 = (Interval(0,10), Interval(5,15))
		self.i0 = (Interval(0,10), Interval(10,15))
		self.i3 = (Interval(0,9), Interval(0,3))
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
	def test_overlap_b(self):
		self.assertEqual(comparisons._overlap_b(*self.i5), 5)
		self.assertEqual(comparisons._overlap_b(*self.i0), 0)
		self.assertEqual(comparisons._overlap_b(*self.i3), 3)
	def test_overlap_p(self):
		self.assertEqual(comparisons._overlap_p(*self.i5), 50)
		self.assertEqual(comparisons._overlap_p(*self.i5[::-1]), 50)
		self.assertEqual(comparisons._overlap_p(*self.i0), 0)
		self.assertEqual(comparisons._overlap_p(*self.i3), float(1.0/3)*100.0)
		self.assertEqual(comparisons._overlap_p(*self.i3[::-1]), 100)
	def test_overlap(self):
		self.assertEqual(comparisons.overlap(*self.i5, overlap_p=60), False)
		self.assertEqual(comparisons.overlap(*self.i5, overlap_p=50), True)
		self.assertEqual(comparisons.overlap(*self.i5[::-1], overlap_p=49), True)
		self.assertEqual(comparisons.overlap(*self.i5[::-1], overlap_p=51), False)
		self.assertEqual(comparisons.overlap(*self.i0, overlap_p=50), False)
		self.assertEqual(comparisons.overlap(*self.i3, overlap_p=50), False)
		self.assertEqual(comparisons.overlap(*self.i3, overlap_p=30), True)
		self.assertEqual(comparisons.overlap(*self.i3[::-1], overlap_p=95), True)
		with self.assertRaises(ValueError):
			comparisons.overlap(*self.i3[::-1], overlap_p=0.5)
	def test_overlap_p(self):
		self.assertEqual(comparisons.overlap_r(*self.i5, overlap_p=40), True)
		self.assertEqual(comparisons.overlap_r(*self.i5[::-1], overlap_p=40), True)
		self.assertEqual(comparisons.overlap_r(*self.i0), False)
		self.assertEqual(comparisons.overlap_r(*self.i0[::-1]), False)
		self.assertEqual(comparisons.overlap_r(*self.i3, overlap_p=40), False)
		self.assertEqual(comparisons.overlap_r(*self.i3, overlap_p=30),True)
class TestSummaries(unittest.TestCase):
	def setUp(self):
		tpath = os.path.dirname(__file__)
		self.fa = os.path.join(tpath, 'test.fa')
		self.fai = os.path.join(tpath, 'test.fa.fai')
		self.gff3_1 = os.path.join(tpath, 'test_1.gff3')
		self.gff3_2 = os.path.join(tpath, 'test_2.gff3')
	def test_gff3_12_tabular(self):
		GI = reader.gff3_interval(self.gff3_1)
		GI.add_gff3(self.gff3_2, 'treat')
		print("")
		summaries.tabular(GI)
		for chrom in GI.get_chrom_set():
			for d in (GI.element_dict, GI.order_dict, GI.sufam_dict):
				for elem in d:
					for strand in '+-B':
						image = 'base_%s_%s_%s.png'%(chrom, strand, elem)
						self.assertTrue(os.path.exists(image))
						os.remove(image)
	def test_gff3_12_tabular_region(self):
		GI = reader.gff3_interval(self.gff3_1, fasta=self.fa)
		GI.add_gff3(self.gff3_2, 'treat')
		print("")
		summaries.tabular_region(GI, p=94)
		for chrom in GI.get_chrom_set():
			for d in (GI.element_dict, GI.order_dict, GI.sufam_dict):
				for elem in d:
					for strand in '+-B':
						for prefix in ('region','length_bp','proportion'):
							image = '%s_%s_%s_%s.png'%(prefix, chrom, strand, elem)
							self.assertTrue(os.path.exists(image))
							os.remove(image)
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

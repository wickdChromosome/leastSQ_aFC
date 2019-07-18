#!/usr/bin/env python3

import unittest
import argparse
import csv
import multiprocessing
import pdb


#functions to be tested
from parser import extract_haplotypes
from parser import get_haplotypes
from parser import get_expressions

from calcs import nonlin_solve
from calcs import linear_estimate




def read_list(filen):
	
	with open(filen, 'r') as f:

		reader = csv.reader(f)
		read_list = list(map(lambda x : x[0] ,list(reader)))
		return read_list		




class testSpheus(unittest.TestCase):
	
	def test_expressions_parsing(self):
		"""
		Using test files, test :
		-If parsed expressions have the intended amount of columns
		-If dataset is a dictionary

		"""
		#input test vars
		variant_id_list = read_list("test_data/expressions_test.variant_id_list")
		expr_file = "test_data/expressions_test.gz" 
		gene_list = read_list("test_data/expressions_test.gene_list")

		out = get_expressions(variant_id_list, gene_list, expr_file)

		#check if number of parsed samples equals that expected by format
		self.assertEqual(len(out.keys()), 299)
		#check if it returns a dictionary
		self.assertEqual(type(out), dict, "get_expressions does not return a dictionary")


	def test_haplotype_extraction(self):
		"""
		Using test files, test :
		-If inner dictionary list dimensions correct 

		"""
		
		variant_list = read_list("test_data/haplotypes_test.variant_id_list")
		line = "chr1\t903285\tchr1_903285_T_C_b38\tT\tC\t.\tPASS\tAF=0.0286396;AN=1676;AC=48\tGT:PG:PB:PI:PW:PC:PM\t0|0:0/0:.:.:0|0:.:.\t0|0:0/0:.:.:0|0:.:."
		haplotype_dict = {}

		#if it returns 2, the test ran (2 is the number of idividual haplotypes in this example)
		self.assertEqual(extract_haplotypes(line, haplotype_dict, variant_list), 2 )





	def test_haplotype_parsing(self):
		"""
		Using test files, test :
		-If parsed haplotypes have the intended amount of columns
		-How many values were nan (when haplotype couldn't be extracted)
		-If dataset contains anything other than 0,1 or nan
		-If dataset is a dictionary with an embedded list

		"""

		#input test vars
		variant_id_list = read_list("test_data/haplotypes_test.variant_id_list")
		vcf_file = "test_data/haplotypes_test.sort.gz" 

		#get output of function
		out = get_haplotypes(variant_id_list, vcf_file)
		#check if number of parsed samples equals that expected by format
		self.assertEqual(len(out.keys()), 3000, "get_expressions does not return all samples in dataset")
		#check if it returns a dictionary
		self.assertEqual(type(out), dict, "get_expressions does not return a dictionary")


#	def test_simple_afc():
		"""
		Using test files, test :
		-If returned slope is within reasonable limits
		-If returned expr is within reasonable limits
		-If returned error is within reasonable limits
		-If output vectors dimensions are as expected

		"""	


#	def test_optimized_afc():
		"""
		Using test files, test :
		-If returned slope is within reasonable limits
		-If returned expr is within reasonable limits
		-If returned error is within reasonable limits
		-If output vector dimensions are as expected

		"""	

		

if __name__ == '__main__':

	unittest.main()

## DATA PARSING FUNCTIONS
import gzip
from helpers import printProgressBar as progbar
import datetime
import pdb
import numpy as np
import pysam
import logging
import time
import sys
import os

def extract_haplotypes( line, haplotype_dict, useful_variants, len_ref):



	if isinstance(line, bytes) and not isinstance(line, str):

		line = line.decode()


	#if this line begins with # ignore it
	if line[0] == '#':

		return

	#pdb.set_trace()

	var_content = line.rstrip().split("\t")

	#if this is the first run
	if len_ref[0] == True:
            
		len_ref[1] = len(var_content)
		len_ref[0] = False
	else:
           
		#inconsistency in columns in file 
		if len_ref[1] != len(var_content):

			print("Haplotype file has inconsistent column lengths")
			sys.exit(0)


	hap1 = []
	hap2 = []

	if var_content[2] in useful_variants:

		#split all genotypes into two matrices for this variant

		colcount = 0
		for x in var_content[9:]:

			hap1_item = np.nan
			hap2_item = np.nan

			if "/" not in x.split(":")[0]:

				[hap1_item, hap2_item] = x.split(":")[0].split('|')

			else:

				[hap1_item, hap2_item] = [np.nan, np.nan]

			hap1.append(hap1_item)
			hap2.append(hap2_item)

		#put variant id haplotypes with that variant_id
		varname = var_content[2]
		haplotype_dict[varname] = [hap1, hap2]


	return len(hap1)

def get_haplotypes( useful_variants, vcf_file):
	"""
	Extract haplotypes from VCF file for each individual
	@params:

		args.vcf	- Global    : input arguments to script

	@returns:

		haplo_dict	- { [variant_id] : [hap1, hap2] }

	"""
	#get all the haplotypes from vcf file using a tabix index
	tabix_haplotypes = pysam.Tabixfile(vcf_file,"r")

	#number of gene we have processed so far
	variant_count = 0

	#the length of the first read record for inconsistency detection, if [0] true then it is the 
	#first row we are processing
	len_ref = [True, 0]

	haplo_dict = {}

	print("getting header")
	
	#get the columns in the file
	command = "zcat " + str(vcf_file)  + " | grep -m 1 '#CHROM' > haplo_header" 
	print(command) 
	os.system(command)
	with open("haplo_header") as f:

		haplo_header = f.readline().replace("\n","").split('\t')	


	#extracted haplotype container
	for variant in useful_variants:

		variant_count += 1
		progbar(variant_count, len(useful_variants), prefix = 'Reading EQTLS:', suffix = 'variants done', length = 80)

		#get coordinates
		[chrom, pos] = variant.split("_")[:2]

		#query vcf using tabix
		try:
			records = tabix_haplotypes.fetch(chrom, int(pos) - 1, int(pos))

		except:

			logging.info("Variant not found: " + str(chrom) + ":" + str(pos) + " in the vcf file")
			continue


		for record in records:

			extract_haplotypes(record,haplo_dict, useful_variants, len_ref)


	return [haplo_dict, haplo_header[4:]]




def get_expressions( useful_genes, expr_file):
	"""
	Extract expressions for each sample from expressions file
	@params:

		args.expr	- Global    : input arguments to script (expressions file)

	Returns:

		expressions_dict  -	{ [chr#_pos#_geneid] = [expressions] }

	"""

	#open expressions file
	expression_stream = gzip.open(expr_file, "r")
    
	#reset line number
	linenum = 0

	expressions_dict = {}

	expressions_header = [] 

	#initialize progress bar
	for line in expression_stream:

		linenum += 1
        
		#skip first line, as those are the labels


		if isinstance(line, bytes) and not isinstance(line, str):

					line = line.decode()
		if line[0] != "#":

			#parse line
			line_content = line.rstrip().split(",")
			#if variant pos and gene match some value
			if line_content[0].split(".")[0] in useful_genes :

				#save the expression data for all the samples

				var_expr = line_content[1:]
				expressions_dict[line_content[0].split(".")[0]] = var_expr
				#processed another variant




			elif line.split(',')[0] == 'Name':
  
				#this is our header
				expressions_header = line.replace("\n","").split(',')

	return [expressions_dict, expressions_header]


import pysam
import pandas as pd
import datetime
import numpy as np
import os
import logging
import gzip
import argparse
import warnings

#sub-libs
from helpers import printProgressBar as progbar
from helpers import intro
from parser import get_haplotypes
from parser import get_expressions
from calcs import nonlin_solve

def get_skipindex(args):


        """
        Calculate for which people there is expressions data
        """

        #get columns from expressions file
        expr_columns = ""

        with gzip.open(args.expr, "r") as f:

                for line in f:
                        

                        line = line.decode()            

                        if "#Chr" in line.split('\t'):

                                expr_columns = line.split('\t')[4:]
                                print('found it')
                                break

        sample_ids = expr_columns

        #get columns from vcf
        columns = ""

        with gzip.open(args.vcf, "r") as f:

                for line in f:
                
                        line = line.decode()
                        if "#CHROM" == line.split('\t')[0]:

                                print('found it')
                                columns = line.split('\t')
                                break

        person_ids = columns


        #now only use columns if there is an expressions entry for it
        #for now just overwrite
        person_ids = list(pd.read_csv("test_data/testdata.vcf.headers",sep='\t').columns)[9:]
        columns2notuse = set(person_ids).symmetric_difference(set(sample_ids))

        skipindex = [[],[]]

        for item in columns2notuse:

                try:
                     skipindex[0].append(person_ids.index(item))
                except:
                     pass

        for item in columns2notuse:

                try:
                     skipindex[1].append(sample_ids.index(item))
                except:
                     pass


        #use this list of indexes later to strip the haplo list 

        return skipindex

def get_expr_genes(args):
    
        """
        Get gene names in the expressions file, to compare to eqtl gene names
        """

        expr_file = args.expr
        os.system("zcat " + expr_file + " | cut -d, -f1 | tail -n +2 > expr_genes.txt") 
        #pdb.set_trace()
        with open("expr_genes.txt") as f:

                lines = f.read().splitlines()   
                return lines

def main():
    
        """
        The main function which gets compares the genes and variants in the three required input files
        and does the least squares optimization on them
        """

        #set up logger
        logging.basicConfig(filename='logs/spheusinfo.log',level=logging.INFO)
        logging.basicConfig(filename='logs/spheusinfo.debug.log',level=logging.DEBUG)
        
        logging.info("Run starting at " + str(datetime.datetime.now))

        #parse input arguments
        parser = argparse.ArgumentParser()
        
        parser.add_argument("--vcf", required=True, help="Genotype VCF")
        parser.add_argument("--expr", required=True, help="Phenotype file")
        parser.add_argument("--eqtl", required=True, help="File containing QTL to calculate allelic fold change for. Should contain tab separated columns 'pid' with phenotype (gene) IDs and 'sid' with SNP IDs. Optionally can include the columns 'sid_chr' and 'sid_pos', which will facilitate tabix retrieval of genotypes, greatly reducing runtime.")
        parser.add_argument("--geno", required=False, default="GT", help="Which field in VCF to use as the genotype. By default 'GT' = genotype. Setting to 'DS' will use dosage rounded to the nearest integer (IE 1.75 = 2 = 1|1).")
        parser.add_argument("--output", "--o", required=True, help="Output file")
        parser.add_argument("--cov", "-c", required=False, help="Covariate matrix")
        parser.add_argument("--isnorm", "-n", required=True, help="1 - matrix has been normalized using DESeq2(or another method) | 0 - matrix hasn't been normalized yet")
        parser.add_argument("--islog", "-l", required=True, help="1 - matrix has been taken the log2 of | 0 - matrix hasn't been log transformed yet")


        #disable warnings
        warnings.filterwarnings("ignore")

        #make input args accessible by all sub-functions
        global args

        args = parser.parse_args()

        #roll intro
        intro(args)

        """
        Find matching genes and variants in expressions, haplotype, and eqtl files, and individual IDs
        """
        #put eqtl data into dataframe
        eqtl_dataf = pd.read_csv(args.eqtl, sep="\t", index_col=False)
        eqtl_dataf['gene_id'] = [x.split(".")[0] for x in eqtl_dataf.gene_id.values]

        
        #initialize the optional args 
        cov_dataf = args.cov
        neednorm = args.isnorm
        needlog = args.islog

        outname = args.output
        is_cov = 1
        
        #if there are covariates provided, read them in
        if cov_dataf != None:
        
            #check if there is a covariate matrix
            cov_dataf = pd.read_csv(cov_dataf, sep="\t", index_col=False)
            
        else:
            
            is_cov = 0

        
        #get list of variants and genes in eqtl file
        eqtl_genes = eqtl_dataf['gene_id'].values
        eqtl_genes = [x.split(".")[0] for x in eqtl_genes]
        eqtl_variants = eqtl_dataf['variant_id'].values

        #now the same for gene expressions
        expr_genes = get_expr_genes(args)
        expr_genes = [x.split(".")[0] for x in expr_genes]

        #compare the sets, but useful variants is not a full list
        useful_genes = set(eqtl_genes) & set(expr_genes)
        useful_variants = set(eqtl_variants)

        #print out number of common genes found in the files
        print("Found useful genes: " + str(len(useful_genes)))
               
        #get individual IDs that should be kept due to overlap in haplo and expr files
        skipindex = get_skipindex(args) 


        """
        Get variants of interest from eqtl file 
        """
        #iterate progress bar
        progbar(0, 100, prefix = 'Reading EQTLS:', suffix = '% done', length = 80)

        #log dataframe specs
        logging.debug("Eqtl dataframe has shape of: " + str(eqtl_dataf.shape))

        #iterate progress bar
        progbar(50, 100, prefix = 'Reading EQTLS:', suffix = '% done', length = 80)

        #get list of genes in eqtl file
        genes = eqtl_dataf["gene_id"]

        #get list of variant ids from eqtl file
        eqtl_variant_ids = eqtl_dataf['variant_id']

        #progress bar full, task finished
        progbar(100, 100, prefix = 'Reading EQTLS:', suffix = '% done', length = 80)
        
        logging.info("Done with reading eqtl data")
        

        """
        Get haplotype vectors        
        """
        #haplotype_vector :      { variant_id : [hap1, hap2] }
        [haplotype_vector, haplotype_header] = get_haplotypes( useful_variants, vcf_file = args.vcf )

        logging.debug("Haplotype vector has " + str(len(haplotype_vector.keys())) + " variants in it")

        progbar( len(eqtl_variant_ids.values), len(eqtl_variant_ids.values), prefix = 'Parsing haplotypes:', suffix = ' variants parsed', length = 75)

        
        """
        Expression vector function desc
        """
        #expressions_vector : {gene_name : [expression1, expression2..]}
        [expressions_vector, expressions_header] = get_expressions(useful_genes, expr_file = args.expr)
        #now convert the haplotypes and expressions into dataframes that are easy to understand
        expressions_df = pd.DataFrame(expressions_vector)
        haplotype0_df = pd.DataFrame(haplotype_vector).transpose()[0].apply(pd.Series)
        haplotype1_df = pd.DataFrame(haplotype_vector).transpose()[1].apply(pd.Series)

        #try splitting the gene names from the version ids - try different delimiters between them
        try:
            expressions_df.index = list(map(lambda x : x.split('-')[0] + '-' + x.split('-')[1], expressions_header[1:]))
        except:
            expressions_df.index = list(map(lambda x : x.split('.')[0] + '-' + x.split('.')[1], expressions_header[1:]))

        
        haplotype0_df.columns = haplotype_header[5:] 
        haplotype1_df.columns = haplotype_header[5:] 

        #drop individuals that are not in both expressions_df and haplotype_df (intercept of two sets)
        expressions_inds = set(expressions_df.index.values) 
        haplotype_inds = set(haplotype0_df.columns.values) 
        
        #make sure that there were no additional apostrophes
        expressions_df.index = [x.replace('"',"").replace("'","") for x in expressions_df.index]
        
        #select only useful individuals that are in all datasets
        useful_inds = list(set.intersection(haplotype_inds, expressions_inds))
        haplotype0_df = haplotype0_df[useful_inds]
        haplotype1_df = haplotype1_df[useful_inds]

        expressions_df = expressions_df[expressions_df.index.isin(useful_inds)]        
        """
        Calculate aFCs, then dump the results into a file
        """
        nonlin_solve(haplotype0_df, haplotype1_df, eqtl_dataf, expressions_df ,useful_genes, cov_dataf, neednorm, outname, is_cov, needlog)
                        



if __name__ == "__main__":
        main()

# Least-squares aFC
## Overview
This script uses a non-least squares fit to calculate allelic fold change. It also calculates the linear approximation for the aFCs for comparison.


## Usage
```
spheus --flag1 --flag2
```

### Use flags

#### Required
```
--vcf The VCF file with phasing data for each variant ID
--expr The expressions file with expressions data for each individual for each gene
--eqtl The eqtl file tie 
--isnorm -n Are the expressions normalized already? 1 for yes, 0 for no
--islog -n Are the expressions log transformed using log2 already? 1 for yes, 0 for no
--output -o
```

#### Optional
```
--cov -c The covariate matrix
```

## Inputs

## Gene counts

The input gene counts should be in the format:
```
Name, sample_id1, sample_id2..
```
Where Name is a column that has the gene ID (such as ENSG00000224533.11), which 
should be in the same format as the gene IDs in the EQTL file

### Normalize counts using DESEQ2
```R
#import library
library(DESeq2)
#read in table 
cts <- rix(read.csv2("example_table.tsv",sep='\t'))
#Set all extra columns other than the counts to NULL
cts$gene_id <- NULL
cts$Names <- NULL
#Delete original column names
colnames(cts) <- NULL
#Read dataset in
dds <- DESeqDataSetFromMatrix(countData = data.matrix(cts), colData = rep.int(1,ncol(cts)) , design = data.matrix(rep.int(1,ncol(cts))))
#Do the transform
dds <- DESeq(dds)
dds <- estimateSizeFactors(dds)
#Get the normalized counts
normalized_counts <- counts(dds, normalized=TRUE)
#Dump the result
write.csv(normalized_counts, "result.csv", row.names=FALSE)
```

## Haplotypes (phasing data)

This should be a VCF file in the format:
```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO  sample1  sample2...
```
where #CHROM is the chr # POS is the position, ID is the variant ID in the format chr1_13550_G_A_b38, REF is the
reference allele, ALT is the alternative allele. 

 
## Eqtl list
This file should contain gene IDs, variant IDs that match - it can also contain other stuff, but it needs to contain at least these two columns:
```
gene_id	variant_id other_stuff1 other_stuff2...
``` 
## Covariate matrix (optional)
This matrix should contain a feature (such as gender) for each individual. The matrix should be in the format num_individuals x num_features
```
ID individual_1 individual_2 individual_3...
```





# Predicting TWAS features (FeaturePred)

FeaturePred is a tool designed to simplify the process of predicting features (e.g. gene expression) in a target sample. It is designed to work with FUSION formated SNP-weights files and PLINK formatted genotype data. It uses a FUSION released script to convert the weights files into a PLINK .score file, which is then used to predict the feature in a target sample using PLINK.

Access to weights files and more information on FUSION can be found [here](http://gusevlab.org/projects/fusion/).



## Getting started

### Prerequisites

* Install the following R packages:
  * data.table
  * optparse
  * foreach
  * doMC

* Target sample genetic data should be in binary PLINK format (.bed/.bim/.fam) with RSIDs that match the FUSION released 1000 Genomes reference data 
* Impute gene expression levels in the target sample:
  * Instructions on how to impute gene expression levels are [here](http://gusevlab.org/projects/fusion/).



### Input files

##### --twas_results

The output of [**FUSION.assoc_test.R** ](https://github.com/gusevlab/fusion_twas/blob/master/FUSION.assoc_test.R) or a file containing the following columns FILE, P0, P1, TWAS.Z, TWAS.P.  Per chromosome files should be combined into a single file. An example is available [here](http://gitlab.psycm.cf.ac.uk/mpmop/gene-expression-risk-scoring/blob/master/ukbiobank-2017-1160-prePRS-fusion.tsv.GW). The file can whitespace or comma delimited. If the file name ends .gz, the file will be assumed to gzipped.

##### --target_gene_exp

A file containing gene expression values for each individual in your target sample. The first two columns are should FID and IID, then each column should contain gene expression data. An example is available [here](http://gitlab.psycm.cf.ac.uk/mpmop/gene-expression-risk-scoring/blob/master/CMC.BRAIN.RNASEQ_GeneX_all_MINI.csv). The gene expression column names must match the values in the FILE column in the --twas_results file. IFRisk ignores the substring before the last '/' and the '.wgt.RDat' string when matching. For example, the column name for the gene expression corresponding to the first value of the [example TWAS results](http://gitlab.psycm.cf.ac.uk/mpmop/gene-expression-risk-scoring/blob/master/ukbiobank-2017-1160-prePRS-fusion.tsv.GW) should be 'CMC.LOC643837'. The file must be readable by the fread function in R. The file can whitespace or comma delimited. If the file name ends .gz, the file will be assumed to gzipped.



### Output files

##### '-GeRS.csv' 

This comma delimited file will contain the gene expression risk scores in the target sample specified. The first two columns are FID and IID, and then each column will contain scores based on the different p-value thresholds specified.

##### '-NGene_Table.csv'

This comma delimited file will contain information on the number of genes surpassing the different p-value threshold specified before and after pruning.

##### '.log'

This is a log file containing general information on the time taken, any errors, the number of genes at different stages and more.



### Optional parameters

##### --prune_thresh

R-squared threshold for pruning genes. 

Default value = 0.9

##### --cor_window

Window for deriving pruning blocks in bases. 

Default value = 5e6

##### --pTs

The p-value thresholds used to derive the risk scores. There must not be spaces between the values.

Default value = '5e-1,1e-1,5e-2,1e-2,1e-3,1e-4,1e-5,1e-6'

##### --prune_mhc

Option to retain only the most significant gene within the MHC region. 

Default value = T



## Examples

##### When using default settings:

```R
Rscript IRIS.V1.0.R \
	--twas_results ukbiobank-2017-1160-prePRS-fusion.tsv.GW \
	--target_gene_exp CMC.BRAIN.RNASEQ_GeneX_all_MINI.csv \
	--output demo
```

##### When using specific p-value thresholds

```R
Rscript IRIS.V1.0.R \
	--twas_results ukbiobank-2017-1160-prePRS-fusion.tsv.GW \
	--target_gene_exp CMC.BRAIN.RNASEQ_GeneX_all_MINI.csv \
	--pTs 1e-5,0.01,0.5 \
	--output demo
```



## Help

This script was written by Dr Oliver Pain under the supervision of Dr Richard Anney whilst at the MRC Centre for Neuropsychiatric Genetics and Genomics, Cardiff University.

If you have any questions or comments please email Ollie (paino@cardiff.ac.uk).








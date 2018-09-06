# Predicting TWAS features (FeaturePred)

FeaturePred is a tool designed to simplify the process of predicting features (e.g. gene expression) in a target sample. It is designed to work with FUSION formated SNP-weights files and PLINK formatted genotype data. It uses a FUSION released script to convert the weights files into a PLINK .SCORE file, which is then used to predict the feature in a target sample using PLINK.

Access to weights files and more information on FUSION can be found [here](http://gusevlab.org/projects/fusion/).



## Getting started

### Prerequisites

* Install the following R packages:
  * data.table

  * optparse

  * foreach

  * doMC

    

* Target sample genetic data

  * Binary PLINK format (.bed/.bim/.fam)

  * RSIDs should match the FUSION LD reference data, which can be downloaded [here](https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2).

    

* FUSION formatted SNP-weights

  * Can be downloaded from the [FUSION website](http://gusevlab.org/projects/fusion/)
  * To make your own FUSION format SNP-weights, see the [FUSION website](http://gusevlab.org/projects/fusion/) and an easy to use pipeline that I wrote [here](http://gitlab.psycm.cf.ac.uk/mpmop/Calculating-FUSION-TWAS-weights-pipeline)



### Provided software and data

* fusion_twas-master - This contains all the fusion scripts and the 1KG reference released by FUSION.

* plink2 - PLINK v1.90b5.4

* test_data/CMC_weights_mini - 10 CMC-based TWAS weights files from the FUSION website

* test_data/CMC_weights_mini.pos - A .pos file listing the weights files in CMC_weights_mini

  

### Input files

##### --PLINK_prefix

Path to genome-wide PLINK binaries (.bed/.bim/.fam)

##### --PLINK_prefix_chr

Path to per chromosome PLINK binaries (.bed/.bim/.fam)

##### --weights

Path for .pos file describing features

##### --weights_dir

Directory containing the weights corresponding to the features in the .pos file

##### --ref_ld_chr

Path to FUSION 1KG reference

##### --make_score_script

Path 'make_score.R

##### --n_cores

Specify the number of cores available for parallel computing.

##### --memory

RAM available in MB.

##### --plink

Path to PLINK software

##### --output

Name of output directory



### Output files

In the specified output directory, the following files will be produced:

##### --FeaturePredictions.csv

Comma delimited file containing FID, IID, and the predicted values for each feature.

##### --FeaturePredictions.log

A log file with information on the analysis.



### Optional parameters

##### --n_cores

Specify the number of cores available for parallel computing.

Default = 1

##### --memory

RAM available in MB.

Default = 2000

##### --save_score

Specify as T if temporary .SCORE files should kept.

Default = F

##### --save_profile

Specify as T if temporary .profile files should kept.

Default = F



## Examples

##### When using default settings:

```R
Rscript FeaturePred.V1.0.R \
	--PLINK_prefix_chr fusion_twas-master/LDREF/1000G.EUR. \
	--weights test_data/CMC_weights_mini.pos \
	--weights_dir test_data/CMC_weights_mini \
	--ref_ld_chr fusion_twas-master/LDREF/1000G.EUR. \
	--make_score_script fusion_twas-master/utils/make_score.R \
	--plink ./plink2 \
	--output demo
```

##### Running in parallel on cluster:

```r
qsub -cwd -b y -l h_vmem=50G,mem_free=50G -e /dev/null -o /dev/null Rscript FeaturePred.V1.0.R \
	--PLINK_prefix_chr fusion_twas-master/LDREF/1000G.EUR. \
	--weights test_data/CMC_weights_mini.pos \
	--weights_dir test_data/CMC_weights_mini \
	--ref_ld_chr fusion_twas-master/LDREF/1000G.EUR. \
	--make_score_script fusion_twas-master/utils/make_score.R \
	--plink ./plink2 \
	--n_cores 5 \
	--memory 50000 \
	--output demo
```



## Help

This script was written by Dr Oliver Pain under the supervision of Dr Richard Anney whilst at the MRC Centre for Neuropsychiatric Genetics and Genomics, Cardiff University.

If you have any questions or comments please email Ollie (paino@cardiff.ac.uk).








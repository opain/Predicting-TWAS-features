# Predicting TWAS features (FeaturePred)

FeaturePred is a tool designed to simplify the process of predicting features (e.g. gene expression) in a target sample. It is designed to work with FUSION formated SNP-weights files and PLINK formatted genotype data. It uses a FUSION released script to convert the weights files into a PLINK .SCORE file, which is then used to predict the feature in a target sample using PLINK.

Access to weights files and more information on FUSION can be found [here](http://gusevlab.org/projects/fusion/).



## Getting started

### Prerequisites

* R and the required packages:
  ```R
  install.packages(c('data.table','optparse','foreach','doMC'))
  ```

* FUSION software:

  ```
  git clone https://github.com/gusevlab/fusion_twas.git
  ```

* FUSION LD reference data ([download](https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2))

* [PLINK 1.9 software](https://www.cog-genomics.org/plink2)

* Target sample genetic data

  * Binary PLINK format (.bed/.bim/.fam) 
  * RSIDs should match the FUSION LD reference data (1000 Genomes phase 3)

* FUSION formatted SNP-weights

  * Can be downloaded from the [FUSION website](http://gusevlab.org/projects/fusion/)
  * To make your own FUSION format SNP-weights, see the [FUSION website](http://gusevlab.org/projects/fusion/) and an easy to use pipeline that I wrote [here](http://gitlab.psycm.cf.ac.uk/mpmop/Calculating-FUSION-TWAS-weights-pipeline)

### Parameters

| Flag                | Description                                                  | Default |
| :------------------ | ------------------------------------------------------------ | :-----: |
| --PLINK_prefix      | Path to genome-wide PLINK binaries (.bed/.bim/.fam)          |   NA    |
| --PLINK_prefix_chr  | Path to per chromosome PLINK binaries (.bed/.bim/.fam)       |   NA    |
| --weights           | Path for .pos file describing features                       |   NA    |
| --weights_dir       | Directory containing the weights listed in the .pos file     |   NA    |
| --ref_ld_chr        | Path to FUSION 1KG reference                                 |   NA    |
| --make_score_script | Path 'make_score.R' script (within FUSION software)          |   NA    |
| --n_cores           | Specify the number of cores available for parallel computing. |    1    |
| --memory            | RAM available in MB.                                         |  2000   |
| --plink             | Path to PLINK software                                       |   NA    |
| --save_score        | Specify as T if temporary .SCORE files should kept.          |  FALSE  |
| --save_profile      | Specify as T if temporary .profile files should kept.        |  FALSE  |
| --output            | Name of output directory                                     |   NA    |



### Output files

In the specified output directory, the following files will be produced:

| Name                   | Description                                                  |
| ---------------------- | ------------------------------------------------------------ |
| FeaturePredictions.csv | Comma delimited file containing FID, IID, and the predicted values for each feature. |
| FeaturePredictions.log | Log file.                                                    |
| SCORE_failed.txt       | Text file listing weights which couldn't be converted to a .SCORE file. |
| Prediction_failed.txt  | Text file listing features that could not be predicted with the reason why. |



## Examples

These examples use the weights and .pos file provided [here](http://gitlab.psycm.cf.ac.uk/mpmop/Predicting-TWAS-features/tree/master/test_data).

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
qsub -cwd -b y -l h_vmem=50G,mem_free=50G Rscript FeaturePred.V1.0.R \
	--PLINK_prefix_chr fusion_twas-master/LDREF/1000G.EUR. \
	--weights test_data/CMC_weights_mini.pos \
	--weights_dir test_data/CMC_weights_mini \
	--ref_ld_chr fusion_twas-master/LDREF/1000G.EUR. \
	--make_score_script fusion_twas-master/utils/make_score.R \
	--plink ./plink2 \
	--n_cores 6 \
	--memory 50000 \
	--output demo
```



## Help

This script was written by Dr Oliver Pain under the supervision of Dr Richard Anney whilst at the MRC Centre for Neuropsychiatric Genetics and Genomics, Cardiff University.

If you have any questions or comments please email Ollie (paino@cardiff.ac.uk).








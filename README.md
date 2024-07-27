# Predicting TWAS features (FeaturePred)

FeaturePred is a tool designed to simplify the process of predicting features (e.g. gene expression) in a target sample. It is designed to work with FUSION formated SNP-weights files and PLINK formatted genotype data. It uses a FUSION released script to convert the weights files into a PLINK .SCORE file, which is then used to predict the feature in a target sample using PLINK. The script harmonises the target sample to the reference data automatically and efficiently handles very large target sample datasets.

Access to weights files and more information on FUSION can be found [here](http://gusevlab.org/projects/fusion/).

## Getting started

### Prerequisites

* R and the required packages:

```R
install.packages(c('data.table','optparse','foreach','doMC'))
```

* FUSION software:

```R
git clone https://github.com/gusevlab/fusion_twas.git
```

* FUSION LD reference data ([download](https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2))
* [PLINK 1.9 software](https://www.cog-genomics.org/plink2)
* [pigz software](https://zlib.net/pigz/) for parallel gz compression
* Target sample genetic data

  * Binary PLINK format (.bed/.bim/.fam) 
  * RSIDs should match the FUSION LD reference data (1000 Genomes phase 3)
* FUSION formatted SNP-weights

  * Can be downloaded from the [FUSION website](http://gusevlab.org/projects/fusion/)
  * To make your own FUSION format SNP-weights, see the [FUSION website](http://gusevlab.org/projects/fusion/) and a pipeline that I wrote [here](http://gitlab.psycm.cf.ac.uk/mpmop/Calculating-FUSION-TWAS-weights-pipeline)



### Parameters

| Flag               | Description                                                  | Default |
| :----------------- | ------------------------------------------------------------ | :-----: |
| --PLINK_prefix     | Path to genome-wide PLINK binaries (.bed/.bim/.fam) [required] |   NA    |
| --PLINK_prefix_chr | Path to per chromosome PLINK binaries (.bed/.bim/.fam) [required] |   NA    |
| --weights          | Path for .pos file describing features [required]            |   NA    |
| --weights_dir      | Directory containing the weights listed in the .pos file [required] |   NA    |
| --ref_ld_chr       | Path to FUSION 1KG reference [required]                      |   NA    |
| --score_files      | Path to SCORE files corresponding to weights [optional]      |   NA    |
| --ref_expr         | Path to reference expression data [optional]                 |   NA    |
| --n_cores          | Specify the number of cores available for parallel computing [optional] |    1    |
| --memory           | RAM available in MB [optional]                               |  2000   |
| --plink            | Path to PLINK software [required]                            |   NA    |
| --save_score       | Specify as T if temporary .SCORE files should kept [optional] |  TRUE   |
| --save_ref_expr    | Save reference expression data [optional]                    |  TRUE   |
| --output           | Name of output directory                                     |   NA    |
| --pigz             | Path to pigz binary [required]                               |   NA    |
| --targ_pred        | Set to FALSE to create SCORE file and expression reference only [optional] |  TRUE   |
| --ref_maf          | Path to per chromosome PLINK freq files [required]           |   NA    |
| --chr              | Specify chromosome number [optional]                         |   NA    |

### Output files

In the specified output directory, the following files will be produced:

| Name                                                  | Description                                                  |
| ----------------------------------------------------- | ------------------------------------------------------------ |
| FeaturePredictions_\<PANEL\>\_chr\<chr\>\.txt.gz .csv | Space delimited file containing FID, IID, and the predicted values for each feature. |
| FeaturePredictions.log                                | Log file.                                                    |
| SCORE_failed.txt                                      | Text file listing weights which couldn't be converted to a .SCORE file (if any). |
| Prediction_failed.txt                                 | Text file listing features that couldn't be predicted (if any). |



## Examples

These examples use the weights and .pos file provided [here](https://github.com/opain/Predicting-TWAS-features/tree/master/test_data).

##### When using default settings:

```shell
Rscript FeaturePred.V2.0.R \
	--PLINK_prefix_chr FUSION/LDREF/1000G.EUR. \
	--weights test_data/CMC.BRAIN.RNASEQ/CMC.BRAIN.RNASEQ.pos \
	--weights_dir test_data/CMC.BRAIN.RNASEQ \
	--ref_ld_chr FUSION/LDREF/1000G.EUR. \
	--plink plink \
	--ref_maf FUSION/LDREF/1000G.EUR. \
	--pigz pigz \
	--chr 1 \
	--output demo
```

##### Running in parallel on cluster:

```shell
sbatch -p shared -n 6 --mem 50G Rscript FeaturePred.V2.0.R \
	--PLINK_prefix_chr FUSION/LDREF/1000G.EUR. \
	--weights test_data/CMC.BRAIN.RNASEQ/CMC.BRAIN.RNASEQ.pos \
	--weights_dir test_data/CMC.BRAIN.RNASEQ \
	--ref_ld_chr FUSION/LDREF/1000G.EUR. \
	--plink plink \
	--ref_maf FUSION/LDREF/1000G.EUR. \
	--pigz pigz \
	--chr 1 \
	--output demo \
	--n_cores 6 \
	--memory 50000
```



## Help

This script was written by Dr Oliver Pain (oliver.pain@kcl.ac.uk)

If you have any questions or comments use the [google group](https://groups.google.com/forum/#!forum/twas-related-r-scripts).








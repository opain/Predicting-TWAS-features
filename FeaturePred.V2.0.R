#!/usr/bin/Rscript
# This script was written by Oliver Pain.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--PLINK_prefix_chr", action="store", default=NA, type='character',
		help="Path to per chromosome PLINK binaries [required]"),
make_option("--weights", action="store", default=NA, type='character',
		help="Path for .pos file describing features [required]"),
make_option("--weights_dir", action="store", default=NA, type='character',
		help="Directory containing the weights corresponding to the features in the .pos file [required]"),
make_option("--ref_ld_chr", action="store", default=NA, type='character',
		help="Path to FUSION 1KG reference [required]"),
make_option("--n_cores", action="store", default=1, type='numeric',
		help="Specify the number of cores available [required]"),
make_option("--score_files", action="store", default=NA, type='character',
		help="Path to SCORE files corresponding to weights [optional]"),
make_option("--ref_expr", action="store", default=NA, type='character',
		help="Path to reference expression data [optional]"),
make_option("--memory", action="store", default=2000, type='numeric',
		help="RAM available in MB [required]"),
make_option("--plink", action="store", default=NA, type='character',
		help="Path to PLINK software [required]"),
make_option("--save_score", action="store", default='T', type='logical',
		help="Save SCORE files [optional]"),
make_option("--save_ref_expr", action="store", default='T', type='logical',
		help="Save reference expression data [optional]"),
make_option("--output", action="store", default=NA, type='character',
    help="Name of output directory [required]"),
make_option("--all_mod", action="store", default=F, type='logical',
    help="Set to T for predictions based on all TWAS models for each gene, rather than the best model alone  [required]"),
make_option("--pigz", action="store", default=NA, type='character',
		help="Path to pigz binary [required]"),
make_option("--targ_pred", action="store", default=T, type='logical',
		help="Set to F to create SCORE file and expression reference only [optional]"),
make_option("--ref_maf", action="store", default=NA, type='character',
		help="Path to per chromosome PLINK .frq files [required]"),
make_option("--chr", action="store", default=NA, type='numeric',
		help="Specify chromosome number [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

system(paste('mkdir -p ',opt$output))

if(is.na(opt$chr)){
	CHROMS<-1:22
} else {
	CHROMS<-opt$chr
}

sink(file = paste(opt$output,'/FeaturePredictions.log',sep=''), append = F)
cat(
'#################################################################
# FeaturePred V2.0
# V1.0 20/02/2019
#################################################################

Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

if(opt$targ_pred == T){
	if(is.na(opt$PLINK_prefix_chr)){
		stop('PLINK_prefix_chr has not been specified.\n')
	}
  missing_target_files <- paste0(opt$PLINK_prefix_chr, CHROMS, '.bed')[!file.exists(paste0(opt$PLINK_prefix_chr, CHROMS, '.bed'))]
  if (length(missing_target_files) > 0) {
    cat(paste0(missing_target_files, collapse = "\n"), "\n")
    stop("The following --PLINK_prefix_chr files are missing:", paste0(missing_target_files, collapse= ', '), "\n")
  }
}

if(is.na(opt$ref_expr)){
  if(is.na(opt$ref_ld_chr)){
    stop('--ref_ld_chr must be specified.')
  }
  missing_ref_files <- paste0(opt$ref_ld_chr, CHROMS, '.bed')[!file.exists(paste0(opt$ref_ld_chr, CHROMS, '.bed'))]
  if (length(missing_ref_files) > 0) {
    cat(paste0(missing_ref_files, collapse = "\n"), "\n")
    stop("The following --ref_ld_chr files are missing:", paste0(missing_ref_files, collapse= ', '), "\n")
  }
}

if(is.na(opt$weights)){
  stop('--weights must be specified.')
}
if(!file.exists(opt$weights)){
  stop('--weights files does not exist.')
}

if(is.na(opt$plink)){
	stop('--plink must be specified.\n')
}

plink_error<-system(paste0(opt$plink),ignore.stdout=T, ignore.stderr=T)
if(plink_error == 127){
  stop('--plink cannot be found. Check the path for plink software.\n')
} else {
	plink_log<-system(paste0(opt$plink,' --help --noweb'),intern=T)
	if(length(grep('PLINK v1.9',plink_log[1])) == 0) {
		cat('\nWarning: Check you are using PLINK v1.9!\n\n')
		system(paste0('rm plink.log'))
	}
}

# Check for required frq files
if(opt$targ_pred == T | is.na(opt$ref_expr)){
  missing_frq_files <- paste0(opt$ref_maf, CHROMS, '.frq')[!file.exists(paste0(opt$ref_maf, CHROMS, '.frq'))]
  if (length(missing_frq_files) > 0) {
    cat(paste0(missing_frq_files, collapse = "\n"), "\n")
    stop("The following --ref_maf files are missing:", paste0(missing_frq_files, collapse= ', '), "\n")
  }
}

if(is.na(opt$pigz)){
  stop('--pigz must be specified.\n')
}

pigz_error<-system(paste0(opt$pigz, ' -h'),ignore.stdout=T, ignore.stderr=T)
if(pigz_error == 127){
  stop('--pigz cannot be found. Check the path for pigz software.\n')
}

library(data.table)
library(foreach)
library(doMC)
registerDoMC(opt$n_cores)

###################################
# Check SNP overlap between FUSION reference and target dataset
###################################

if(opt$targ_pred == T){

	sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
	cat('Harmonsising reference and target data...',sep='')
	sink()

	# Read in the SNPs in reference sample
	Ref<-NULL
	for(i in 1:22){
	  Ref<-rbind(Ref, fread(paste0(opt$ref_ld_chr,i,'.bim')))
	}

	if(sum(duplicated(Ref$V2)) > 0){
		sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
		cat('Duplicates are present in the reference data.\n')
		sink()
		q()
	}

	Target<-NULL
	for(i in 1:22){
	  Target<-rbind(Target, fread(paste0(opt$PLINK_prefix_chr,i,'.bim')))
	}

	if(sum(duplicated(Ref$V2)) > 0){
		sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
		cat('Duplicates are present in the target data.\n')
		sink()
		q()
	}

	# Idenitfy SNPs that are in both datasets and have matching allele codes
	Ref$IUPAC<-NA
	Ref$IUPAC[Ref$V5 == 'A' & Ref$V6 =='T' | Ref$V5 == 'T' & Ref$V6 =='A']<-'W'
	Ref$IUPAC[Ref$V5 == 'C' & Ref$V6 =='G' | Ref$V5 == 'G' & Ref$V6 =='C']<-'S'
	Ref$IUPAC[Ref$V5 == 'A' & Ref$V6 =='G' | Ref$V5 == 'G' & Ref$V6 =='A']<-'R'
	Ref$IUPAC[Ref$V5 == 'C' & Ref$V6 =='T' | Ref$V5 == 'T' & Ref$V6 =='C']<-'Y'
	Ref$IUPAC[Ref$V5 == 'G' & Ref$V6 =='T' | Ref$V5 == 'T' & Ref$V6 =='G']<-'K'
	Ref$IUPAC[Ref$V5 == 'A' & Ref$V6 =='C' | Ref$V5 == 'C' & Ref$V6 =='A']<-'M'

	if(sum(is.na(Ref$IUPAC)) > 0){
	  sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
	  cat(sum(is.na(Ref$IUPAC)),"non-SNP variants removed from reference data (Ref_non_SNPs.txt).\n")
	  sink()

	  write.table(Ref[is.na(Ref$IUPAC),], paste0(opt$output,'/Ref_non_SNPs.txt'), col.names=T, row.names=F, quote=F)

	  Ref<-Ref[!is.na(Ref$IUPAC),]
	}

	Target$IUPAC<-NA
	Target$IUPAC[Target$V5 == 'A' & Target$V6 =='T' | Target$V5 == 'T' & Target$V6 =='A']<-'W'
	Target$IUPAC[Target$V5 == 'C' & Target$V6 =='G' | Target$V5 == 'G' & Target$V6 =='C']<-'S'
	Target$IUPAC[Target$V5 == 'A' & Target$V6 =='G' | Target$V5 == 'G' & Target$V6 =='A']<-'R'
	Target$IUPAC[Target$V5 == 'C' & Target$V6 =='T' | Target$V5 == 'T' & Target$V6 =='C']<-'Y'
	Target$IUPAC[Target$V5 == 'G' & Target$V6 =='T' | Target$V5 == 'T' & Target$V6 =='G']<-'K'
	Target$IUPAC[Target$V5 == 'A' & Target$V6 =='C' | Target$V5 == 'C' & Target$V6 =='A']<-'M'

	if(sum(is.na(Target$IUPAC)) > 0){
	  sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
	  cat(sum(is.na(Target$IUPAC)),'non-SNP variants removed from target data (Targ_non_SNPs.txt).\n')
	  sink()

	  write.table(Target[is.na(Target$IUPAC),], paste0(opt$output,'/Targ_non_SNPs.txt'), col.names=T, row.names=F, quote=F)

	  Target<-Target[!is.na(Target$IUPAC),]
	}

	Ref_Target<-merge(Ref, Target, by='V2')

	n_flip<-sum(Ref_Target$IUPAC.x == 'R' & Ref_Target$IUPAC.y == 'Y' |
							Ref_Target$IUPAC.x == 'Y' & Ref_Target$IUPAC.y == 'R' |
							Ref_Target$IUPAC.x == 'K' & Ref_Target$IUPAC.y == 'M' |
							Ref_Target$IUPAC.x == 'M' & Ref_Target$IUPAC.y == 'K' )

	if(n_flip > 0){
		sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
		cat(n_flip,'SNPs will be flipped!\n')
		sink()
		flip_list<-		Ref_Target$V2[Ref_Target$IUPAC.x == 'R' & Ref_Target$IUPAC.y == 'Y' |
									Ref_Target$IUPAC.x == 'Y' & Ref_Target$IUPAC.y == 'R' |
									Ref_Target$IUPAC.x == 'K' & Ref_Target$IUPAC.y == 'M' |
									Ref_Target$IUPAC.x == 'M' & Ref_Target$IUPAC.y == 'K' ]
		write.table(flip_list, paste0(opt$output,'/flip.snplist'), col.names=F, row.names=F, quote=F)
	}

	# Idenitfy variants that match
	incl<-Ref_Target$V2[	Ref_Target$IUPAC.x == 'R' & Ref_Target$IUPAC.y == 'R' |
												Ref_Target$IUPAC.x == 'Y' & Ref_Target$IUPAC.y == 'Y' |
												Ref_Target$IUPAC.x == 'K' & Ref_Target$IUPAC.y == 'K' |
												Ref_Target$IUPAC.x == 'M' & Ref_Target$IUPAC.y == 'M' ]

	sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
	cat(paste0('Done!\n'))
	cat(paste0('Number of SNPs in FUSION LD Reference = ',length(Ref$V2),'.\n'))
	cat(paste0('Number of SNPs in target PLINK files = ',length(Target$V2),'.\n'))
	cat(paste0('Number of SNPs in both = ',length(incl),'.\n'))
	cat(paste0('Percentage of SNPs in FUSION LD Reference that are in the target PLINK files = ',round(length(incl)/length(Ref$V2)*100,2),'%.\n'))
	sink()

	# Save list of SNPs that intersect the target and ref
	fwrite(list(incl),paste0(opt$output,'/intersect.snplist'),col.names=F, na='NA', quote=F)

	rm(Ref)
	rm(Target)
	rm(Ref_Target)
	rm(incl)
	gc()
}

###################################
# Convert FUSION weights into SCORE files
###################################

# Read in the .pos file
pos<-data.frame(fread(opt$weights, nThread=opt$n_cores))
if(all(names(pos) != 'PANEL')){
  stop("--weight file does not contain PANEL column.\n")
}
pos<-pos[pos$CHR %in% CHROMS,]

sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
cat('The .pos file contains ',dim(pos)[1],' features.\n',sep='')
sink()

# Attach weights directory to WGT values in pos file
pos$FILE<-paste0(opt$weights_dir,'/',pos$PANEL,'/',sub('.*/','',pos$WGT))

# Remove .wgt.RDat from the WGT values
pos$WGT<-gsub('.wgt.RDat','',pos$WGT)
pos$WGT<-gsub('.*/','',pos$WGT)

if(is.na(opt$score_files)){

  # Read in the SNPs in reference sample
  Ref<-NULL
  for(i in 1:22){
    Ref<-rbind(Ref, fread(paste0(opt$ref_ld_chr,i,'.bim')))
  }

  Ref$IUPAC<-NA
  Ref$IUPAC[Ref$V5 == 'A' & Ref$V6 =='T' | Ref$V5 == 'T' & Ref$V6 =='A']<-'W'
  Ref$IUPAC[Ref$V5 == 'C' & Ref$V6 =='G' | Ref$V5 == 'G' & Ref$V6 =='C']<-'S'
  Ref$IUPAC[Ref$V5 == 'A' & Ref$V6 =='G' | Ref$V5 == 'G' & Ref$V6 =='A']<-'R'
  Ref$IUPAC[Ref$V5 == 'C' & Ref$V6 =='T' | Ref$V5 == 'T' & Ref$V6 =='C']<-'Y'
  Ref$IUPAC[Ref$V5 == 'G' & Ref$V6 =='T' | Ref$V5 == 'T' & Ref$V6 =='G']<-'K'
  Ref$IUPAC[Ref$V5 == 'A' & Ref$V6 =='C' | Ref$V5 == 'C' & Ref$V6 =='A']<-'M'

	opt$score_files<-paste0(opt$output,'/SCORE_FILES')

	system(paste0('mkdir -p ',opt$score_files))

	# Create SCORE file for each set of weights using FUSION make_score.R
	sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
	cat('Converting weights files into PLINK SCORE files...',sep='')
	sink()

	tmp<-foreach(i=1:length(pos$FILE), .combine=c) %dopar% {
	  load(pos$FILE[i])

	  # Resolve strand flips with reference
	  snps$IUPAC<-NA
	  snps$IUPAC[snps$V5 == 'A' & snps$V6 =='T' | snps$V5 == 'T' & snps$V6 =='A']<-'W'
	  snps$IUPAC[snps$V5 == 'C' & snps$V6 =='G' | snps$V5 == 'G' & snps$V6 =='C']<-'S'
	  snps$IUPAC[snps$V5 == 'A' & snps$V6 =='G' | snps$V5 == 'G' & snps$V6 =='A']<-'R'
	  snps$IUPAC[snps$V5 == 'C' & snps$V6 =='T' | snps$V5 == 'T' & snps$V6 =='C']<-'Y'
	  snps$IUPAC[snps$V5 == 'G' & snps$V6 =='T' | snps$V5 == 'T' & snps$V6 =='G']<-'K'
	  snps$IUPAC[snps$V5 == 'A' & snps$V6 =='C' | snps$V5 == 'C' & snps$V6 =='A']<-'M'

	  snps_ref<-merge(snps, Ref[,c('V2','V5','V6','IUPAC')], by='V2')

	  flip_list<-		snps_ref$V2[snps_ref$IUPAC.x == 'R' & snps_ref$IUPAC.y == 'Y' |
	                              snps_ref$IUPAC.x == 'Y' & snps_ref$IUPAC.y == 'R' |
	                              snps_ref$IUPAC.x == 'K' & snps_ref$IUPAC.y == 'M' |
	                              snps_ref$IUPAC.x == 'M' & snps_ref$IUPAC.y == 'K' ]

	  snp_allele_comp<-function(x=NA){
	    x_new<-x
	    x_new[x == 'A']<-'T'
	    x_new[x == 'T']<-'A'
	    x_new[x == 'G']<-'C'
	    x_new[x == 'C']<-'G'
	    x_new[!(x %in% c('A','T','G','C'))]<-NA
	    return(x_new)
	  }

	  snps$IUPAC<-NULL
	  snps$V5[snps$V2 %in% flip_list]<-snp_allele_comp(snps$V5[snps$V2 %in% flip_list])
	  snps$V6[snps$V2 %in% flip_list]<-snp_allele_comp(snps$V6[snps$V2 %in% flip_list])

	  # Remove models with invalid (all NA or invariant weights)
	  valid_snp_columns <- apply(wgt.matrix, 2, function(col) !all(is.na(col)))
	  non_invariant_columns <- apply(wgt.matrix, 2, function(col) length(unique(col[!is.na(col)])) > 1)
	  valid_models <- valid_snp_columns & non_invariant_columns
	  valid_models <- names(valid_models[valid_models])
	  
	  cv.performance<-cv.performance[,colnames(cv.performance) %in% valid_models]
	  wgt.matrix<-wgt.matrix[,colnames(wgt.matrix) %in% valid_models]
	  
	  if(opt$all_mod == F){
  	  best = which.min(cv.performance[2,])

      if ( names(best) == "top1" ) {
    		keep = which.max(wgt.matrix[,best]^2)
  	  } else {
    		keep = 1:nrow(wgt.matrix)
  	  }

  	  write.table(format(cbind( (snps[,c(2,5,6)]) , wgt.matrix[,best])[keep,], digits=3), paste0(opt$score_files,'/',pos$WGT[i],'.SCORE') , quote=F , row.names=F , col.names=F , sep='\t' )

  	  # Write out a snplist for each SCORE file to reduce compution time for scoring
  		system(paste0('cut -f 1 ',opt$score_files,'/',pos$WGT[i],'.SCORE > ',opt$score_files,'/',pos$WGT[i],'.snplist'))
	  } else {
	    for(mod in colnames(wgt.matrix)){
	      if (mod == "top1" ) {
	        wgt.matrix[-which.max(wgt.matrix[,mod]^2),mod]<-0
	      }

  	    write.table(format(cbind((snps[,c(2,5,6)]) , wgt.matrix[,mod]),digits=3), paste0(opt$score_files,'/',pos$WGT[i],'.',mod,'.SCORE') , quote=F , row.names=F , col.names=F , sep='\t' )
	    }
	    # Write out a snplist for each SCORE file to reduce compution time for scoring
	    write.table(snps$V2, paste0(opt$score_files,'/',pos$WGT[i],'.snplist') , quote=F , row.names=F , col.names=F , sep='\t' )
	  }
	}

	sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
	cat('Done!\n',sep='')
	sink()
} else {
	sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
	cat('Using precomputed score files.\n',sep='')
	sink()
}

###################################
# Predict features in the reference sample
###################################

# Read in reference fam file
ref_fam<-fread(paste0(opt$ref_ld_chr,'1.fam'), nThread=opt$n_cores)
ref_fam<-ref_fam[,1:2]
names(ref_fam)<-c('FID','IID')

if(is.na(opt$ref_expr)){

	# Create directory for the PROFILE files
	system(paste0('mkdir ',opt$output,'/REF_PROFILE_FILES'))

	# Calculate profile scores (i.e. feature predictions)
	sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
	cat('Predicting features in reference sample...\n',sep='')
	sink()

	# Create column IDs to be combined to the feature predictions.
	write.table(ref_fam, paste0(opt$output,'/REF_PROFILE_FILES/REF.IDs'), col.names=T, row.names=F, quote=F, sep='\t')

	error_table_all<-NULL

	error_table<-foreach(i=1:length(pos$FILE), .combine=rbind) %dopar% {

	  if(opt$all_mod == F){

  		# Calculate feature predictions
  		error<-system(paste0(opt$plink,' --bfile ',opt$ref_ld_chr,pos$CHR[i],' --extract ',opt$score_files,'/',pos$WGT[i],'.snplist --allow-no-sex --read-freq ',opt$ref_maf,pos$CHR[i],'.frq --score ',opt$score_files,'/',pos$WGT[i],'.SCORE 1 2 4 --out ',opt$output,'/REF_PROFILE_FILES/',pos$WGT[i],' --memory ', floor((opt$memory*0.4)/opt$n_cores)),ignore.stdout=T, ignore.stderr=T)
  		# Delete temporary files and extract feature prediction column to reduce disk space
  		system(paste0("awk '{print $6}' ",opt$output,'/REF_PROFILE_FILES/',pos$WGT[i],'.profile | tail -n +2 > ',opt$output,'/REF_PROFILE_FILES/',pos$WGT[i],'.profile_mini'),intern=T)
  		system(paste0("echo ",pos$PANEL,'.',pos$WGT[i]," | cat - ",opt$output,'/REF_PROFILE_FILES/',pos$WGT[i],'.profile_mini > ',opt$output,'/REF_PROFILE_FILES/',pos$WGT[i],'.profile_mini_tmp && mv ',opt$output,'/REF_PROFILE_FILES/',pos$WGT[i],'.profile_mini_tmp ',opt$output,'/REF_PROFILE_FILES/',pos$WGT[i],'.profile_mini'),intern=T)
  		system(paste0('rm ',opt$output,'/REF_PROFILE_FILES/',pos$WGT[i],'.profile'),ignore.stdout=T, ignore.stderr=T)
  		system(paste0('rm ',opt$output,'/REF_PROFILE_FILES/',pos$WGT[i],'.nosex'),ignore.stdout=T, ignore.stderr=T)
  		data.frame(N=i,Error=error)

	  } else {

	    load(pos$FILE[i])
	    error_tmp<-NULL
	    for(mod in colnames(wgt.matrix)){
  	    # Calculate feature predictions
  	    error<-system(paste0(opt$plink,' --bfile ',opt$ref_ld_chr,pos$CHR[i],' --extract ',opt$score_files,'/',pos$WGT[i],'.snplist --allow-no-sex --read-freq ',opt$ref_maf,pos$CHR[i],'.frq --score ',opt$score_files,'/',pos$WGT[i],'.',mod,'.SCORE 1 2 4 --out ',opt$output,'/REF_PROFILE_FILES/',pos$WGT[i],'.',mod,' --memory ', floor((opt$memory*0.4)/opt$n_cores)),ignore.stdout=T, ignore.stderr=T)
  	    if(file.exists(paste0(opt$output,'/REF_PROFILE_FILES/',pos$WGT[i],'.',mod,'.profile'))){
    	    # Delete temporary files and extract feature prediction column to reduce disk space
    	    system(paste0("awk '{print $6}' ",opt$output,'/REF_PROFILE_FILES/',pos$WGT[i],'.',mod,'.profile | tail -n +2 > ',opt$output,'/REF_PROFILE_FILES/',pos$WGT[i],'.',mod,'.profile_mini'),intern=T)
    	    system(paste0("echo ",pos$PANEL,'.',pos$WGT[i],'.',mod," | cat - ",opt$output,'/REF_PROFILE_FILES/',pos$WGT[i],'.',mod,'.profile_mini > ',opt$output,'/REF_PROFILE_FILES/',pos$WGT[i],'.',mod,'.profile_mini_tmp && mv ',opt$output,'/REF_PROFILE_FILES/',pos$WGT[i],'.',mod,'.profile_mini_tmp ',opt$output,'/REF_PROFILE_FILES/',pos$WGT[i],'.',mod,'.profile_mini'),intern=T)
    	    system(paste0('rm ',opt$output,'/REF_PROFILE_FILES/',pos$WGT[i],'.',mod,'.profile'),ignore.stdout=T, ignore.stderr=T)
    	    system(paste0('rm ',opt$output,'/REF_PROFILE_FILES/',pos$WGT[i],'.',mod,'.nosex'),ignore.stdout=T, ignore.stderr=T)
    	    error_tmp<-rbind(error_tmp, data.frame(N=i,Model=mod,Error=error))
  	    }
	    }
	    error_tmp
	  }
	}

	# Split feature predictions into list of <1000 files, then paste each list of files in batches, and then past all batches.
	system(paste0("find ", opt$output,"/REF_PROFILE_FILES -type f -name '*.profile_mini' | split -l 500 -a 4 -d - ",opt$output,"/profile_mini_lists"))
	system(paste0("echo ", opt$output, "/REF_PROFILE_FILES/REF.IDs | cat - ",opt$output,"/profile_mini_lists0000 > ",opt$output,"/profile_mini_lists0000_temp && mv ",opt$output,"/profile_mini_lists0000_temp ",opt$output,"/profile_mini_lists0000"))
	tmp<-foreach(k=list.files(path=opt$output, pattern="^profile_mini_lists*"), .combine=c) %dopar% {
		exit_status<-system(paste0("paste $(cat ",opt$output,"/",k,") > ", opt$output,"/merge_",k))
		if(exit_status != 0){
		  stop()
		}
	}

	exit_status<-system(paste0("paste ", opt$output,"/merge_profile_mini_lists* > ", opt$output,"/REF_expr.txt"))
	if(exit_status != 0){
	  stop()
	}
	
	# Delete temporary files
	system(paste0("rm ",opt$output,'/profile_mini_lists*'))
	system(paste0("rm ",opt$output,'/merge_profile_mini_lists*'))
	system(paste0("rm ",opt$output,"/REF_PROFILE_FILES/*.profile_mini"))

	# Output file containing list of features that couldn't be predicted.
	error_table$Error[error_table$Error > 0]<-1
	error_table<-error_table[error_table$Error > 0,]
	if(dim(error_table)[1] > 0){
		for(i in 1:dim(error_table)[1]){
			k<-error_table$N[i]
			error_table$ID[i]<-pos$WGT[k]
			tmp1<-read.table(paste0(opt$output,'/REF_PROFILE_FILES/',pos$WGT[k],'.log'),sep='*')
			NoValid<-sum(grepl('Error: No valid entries in --score file.',tmp1$V1))
			error_table$Reason[i]<-'Error: No valid entries in --score file.'
		}
		error_table$Error<-NULL
		error_table<-error_table[c('ID','Reason')]
		error_table_all<-rbind(error_table_all,error_table)

		write.table(error_table_all, paste0(opt$output,'/Prediction_failed.txt'), col.names=T, row.names=F, quote=T)
		sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
			cat(paste0(dim(error_table_all)[1],' feature/s cannot not be predicted (',opt$output,'/Prediction_failed.txt)\n'))
		sink()
		rm(error_table)
		rm(error_table_all)
		gc()
	}
	system(paste0("rm ",opt$output,"/REF_PROFILE_FILES/*.log"))

	sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
		cat('Done!\n')
	sink()
} else {
	sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
	cat('Using precomputed reference expression data.\n',sep='')
	sink()
}

###################################
# Calculate the mean and SD of feature predictions in the reference
###################################

if(is.na(opt$ref_expr)){
	REF_expr<-fread(paste0(opt$output,"/REF_expr.txt"), nThread=opt$n_cores)

	ref_scale<-data.frame(	ID=names(REF_expr[,-1:-2]),
							Mean=sapply(REF_expr[,-1:-2], function(x) mean(x)),
							SD=sapply(REF_expr[,-1:-2], function(x) sd(x)))

	if(opt$save_ref_expr == T){
		system(paste0('mkdir ',opt$output,'/Reference_Expression'))
		fwrite(ref_scale, paste0(opt$output,'/Reference_Expression/Scale.txt'), sep=' ', nThread=opt$n_cores, na='NA')
		REF_expr_scaled<-REF_expr
		REF_expr_scaled[ , (names(REF_expr_scaled)[-1:-2]) := lapply(.SD, function(x) round(as.numeric(scale(x)),3)), .SDcols = (names(REF_expr_scaled)[-1:-2])]
		for(panel in unique(pos$PANEL)){
			REF_expr_scaled_pan<-REF_expr_scaled[,c(T,T,grepl(panel, names(REF_expr_scaled)[-1:-2])), with=FALSE]
			fwrite(REF_expr_scaled_pan, paste0(opt$output,'/Reference_Expression/Reference_Expression_',panel,'.txt'), nThread = opt$n_cores, sep=' ', na='NA')
			system(paste0(opt$pigz,' ',opt$output,'/Reference_Expression/Reference_Expression_',panel,'.txt'))
		}
	}

	system(paste0('rm ',opt$output,'/REF_expr.txt'))
	rm(REF_expr)
	rm(REF_expr_scaled)
	gc()
} else {
	ref_scale<-data.frame(fread(paste0(opt$ref_expr,'/Scale.txt'), nThread=opt$n_cores, fill=T))
}

if(opt$targ_pred == T){

	###################################
	# Extract single individual from reference to be merged with the target to insert missing SNPs.
	###################################

	write.table(ref_fam[1,], paste0(opt$output,'/ref_keep'), col.names=F, row.names=F, quote=F)
	rm(ref_fam)
	gc()

	foreach(chr=CHROMS, .combine=c) %dopar% {
		system(paste0(opt$plink,' --bfile ',opt$ref_ld_chr,chr,' --keep ',opt$output,'/ref_keep --make-bed --out ',opt$output,'/ref_indiv_chr',chr,' --memory ', floor((opt$memory*0.4)/opt$n_cores)),ignore.stdout=T, ignore.stderr=T)

		# Update ref indiv ID to be distinct from target samples
		ref_fam<-fread(paste0(opt$output,'/ref_indiv_chr',chr,'.fam'))
		ref_fam$V1<-paste0('REF_',ref_fam$V1)
		ref_fam$V2<-paste0('REF_',ref_fam$V2)
		write.table(ref_fam, paste0(opt$output,'/ref_indiv_chr',chr,'.fam'), col.names=F, row.names=F, quote=F)
		rm(ref_fam)
		gc()
	}

	###################################
	# Calculate feature predictions in the target sample
	###################################

	# Crate a directory to store the predicted expression
	system(paste0('mkdir ',opt$output,'/TARG_PROFILE_FILES'))

	# Create an file containing IDs to be merged with predicted expression
	targ_fam<-fread(paste0(opt$PLINK_prefix_chr,'1.fam'), nThread=opt$n_cores)
	targ_fam_id<-targ_fam[,1:2]
	names(targ_fam_id)<-c('FID','IID')
	write.table(targ_fam_id, paste0(opt$output,'/TARG_PROFILE_FILES/TARG.IDs'), col.names=T, row.names=F, quote=F, sep=' ')

	rm(targ_fam)
	rm(targ_fam_id)
	gc()

	# Calculate and process gene expression by chromosome to avoid large memory requirement.
	for(chr in CHROMS[CHROMS %in% unique(pos$CHR)]){

		# Create a directory for predicted expression of genes on chromosome
		system(paste0('mkdir ',opt$output,'/TARG_PROFILE_FILES/chr',chr))

	  # Remove any duplicate or 3+ allele variants
	  system(paste0(opt$plink,' --bfile ',opt$PLINK_prefix_chr,chr,' --bmerge ',opt$PLINK_prefix_chr,chr,' --merge-mode 6 --make-bed --out ',opt$output,'/targ_chr',chr,'_noDup'))

	  if(file.exists(paste0(opt$output,'/targ_chr',chr,'_noDup.missnp'))){
	    missnp<-fread(paste0(opt$output,'/targ_chr',chr,'_noDup.missnp'), header=F)$V1
	    log_chr<-readLines(paste0(opt$output,'/targ_chr',chr,'_noDup.log'))
	    log_chr<-log_chr[grepl("Warning: Multiple positions seen for variant", log_chr)]
	    log_chr<-unique(gsub("'\\.","",gsub(".*variant '","",log_chr)))
	    missnp<-c(missnp, log_chr)
	    write.table(missnp, paste0(opt$output,'/targ_chr',chr,'_noDup.missnp'), col.names=F, row.names=F, quote=F)

	    # Extract SNPs in reference and flip SNPs if necessary.
	    # And remove any duplicate 3+ allele variants
	    system(paste0(opt$plink,' --bfile ',opt$PLINK_prefix_chr,chr,' --make-bed --exclude ',opt$output,'/targ_chr',chr,'_noDup.missnp --out ',opt$output,'/targ_chr',chr,' --memory ',floor((opt$memory*0.4))))

	    if(n_flip > 0){
	      system(paste0(opt$plink,' --bfile ',opt$output,'/targ_chr',chr,' --make-bed --extract ',opt$output,'/intersect.snplist --flip ',opt$output,'/flip.snplist --out ',opt$output,'/targ_chr',chr,' --memory ',floor((opt$memory*0.4))))
	    } else {
	      system(paste0(opt$plink,' --bfile ',opt$output,'/targ_chr',chr,' --make-bed --extract ',opt$output,'/intersect.snplist --out ',opt$output,'/targ_chr',chr,' --memory ',floor((opt$memory*0.4))))
	    }
	  } else {
	    # Extract SNPs in reference and flip SNPs if necessary.
	    if(n_flip > 0){
	      system(paste0(opt$plink,' --bfile ',opt$PLINK_prefix_chr,chr,' --make-bed --extract ',opt$output,'/intersect.snplist --flip ',opt$output,'/flip.snplist --out ',opt$output,'/targ_chr',chr,' --memory ',floor((opt$memory*0.4))))
	    } else {
	      system(paste0(opt$plink,' --bfile ',opt$PLINK_prefix_chr,chr,' --make-bed --extract ',opt$output,'/intersect.snplist --out ',opt$output,'/targ_chr',chr,' --memory ',floor((opt$memory*0.4))))
	    }
	  }

		# Merge and remove single individual with the target data to insert missing reference SNPs.
		system(paste0(opt$plink,' --bfile ',opt$output,'/targ_chr',chr,' --bmerge ',opt$output,'/ref_indiv_chr',chr,' --make-bed --out ',opt$output,'/ref_targ_chr',chr,' --memory ',floor((opt$memory*0.4))))
		system(paste0(opt$plink,' --bfile ',opt$output,'/ref_targ_chr',chr,' --remove ',opt$output,'/ref_indiv_chr',chr,'.fam --make-bed --out ',opt$output,'/targ_chr',chr,' --memory ',floor((opt$memory*0.4))))

		# Delete the temporary files
		system(paste0('rm ',opt$output,'/ref_targ_chr*'))

		sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
		cat('Predicting features in chromosome ',chr,'...',sep='')
		sink()

		# Subset features on the chromosome
		pos_chr<-pos[pos$CHR == chr,]

		# Predict each panel seperately to avoid large files when many panels in one run.
		for(panel in unique(pos_chr$PANEL)){
			pos_chr_panel<-pos_chr[pos_chr$PANEL == panel,]
			TARG_expr<-foreach(i=1:length(pos_chr_panel$FILE), .combine=cbind) %dopar% {
			  if(opt$all_mod == F){

  				# Calculate feature predictions
  				tmp<-system(paste0(opt$plink,' --bfile ',opt$output,'/targ_chr',chr,' --extract ',opt$score_files,'/',pos_chr_panel$WGT[i],'.snplist --allow-no-sex --read-freq ',opt$ref_maf,chr,'.frq --score ',opt$score_files,'/',pos_chr_panel$WGT[i],'.SCORE 1 2 4 --out ',opt$output,'/TARG_PROFILE_FILES/chr',chr,'/',panel,'_',pos_chr_panel$WGT[i],' --memory ', floor((opt$memory*0.4)/opt$n_cores)),ignore.stdout=T, ignore.stderr=T)

  				if(tmp != 0){
  						# Delete temporary files
  						system(paste0('rm ',opt$output,'/TARG_PROFILE_FILES/chr',chr,'/',panel,'_',pos_chr_panel$WGT[i],'.*'),ignore.stdout=T, ignore.stderr=T)
  						return(NULL)
  				}

  				# Read in the predictions, extract SCORE column and change header.
  				feature<-fread(paste0(opt$output,'/TARG_PROFILE_FILES/chr',chr,'/',panel,'_',pos_chr_panel$WGT[i],'.profile'), nThread=1)
  				feature<-feature[,6]
  				names(feature)<-paste0(panel,'.',pos_chr_panel$WGT[i])

  				# Delete temporary files
  				system(paste0('rm ',opt$output,'/TARG_PROFILE_FILES/chr',chr,'/',panel,'_',pos_chr_panel$WGT[i],'.*'),ignore.stdout=T, ignore.stderr=T)

  				# Scale expression to the reference and round
  				feature<-feature-ref_scale$Mean[ref_scale$ID == names(feature)]
  				feature<-feature/ref_scale$SD[ref_scale$ID == names(feature)]
  				feature<-round(feature,3)

  				feature
			  } else {
			    load(pos_chr_panel$FILE[i])
			    feature_tmp<-NULL
			    for(mod in colnames(wgt.matrix)){
			      # Calculate feature predictions
			      tmp<-system(paste0(opt$plink,' --bfile ',opt$output,'/targ_chr',chr,' --extract ',opt$score_files,'/',pos_chr_panel$WGT[i],'.snplist --allow-no-sex --read-freq ',opt$ref_maf,chr,'.frq --score ',opt$score_files,'/',pos_chr_panel$WGT[i],'.',mod,'.SCORE 1 2 4 --out ',opt$output,'/TARG_PROFILE_FILES/chr',chr,'/',panel,'_',pos_chr_panel$WGT[i],'.',mod,' --memory ', floor((opt$memory*0.4)/opt$n_cores)),ignore.stdout=T, ignore.stderr=T)

			      if(tmp != 0){
			        # Delete temporary files
			        system(paste0('rm ',opt$output,'/TARG_PROFILE_FILES/chr',chr,'/',panel,'_',pos_chr_panel$WGT[i],'.',mod,'.*'),ignore.stdout=T, ignore.stderr=T)
			        return(NULL)
			      }

			      # Read in the predictions, extract SCORE column and change header.
			      feature<-fread(paste0(opt$output,'/TARG_PROFILE_FILES/chr',chr,'/',panel,'_',pos_chr_panel$WGT[i],'.',mod,'.profile'), nThread=1)
			      feature<-feature[,6]
			      names(feature)<-paste0(panel,'.',pos_chr_panel$WGT[i],'.',mod)

			      # Delete temporary files
			      system(paste0('rm ',opt$output,'/TARG_PROFILE_FILES/chr',chr,'/',panel,'_',pos_chr_panel$WGT[i],'.',mod,'.*'),ignore.stdout=T, ignore.stderr=T)

			      # Scale expression to the reference and round
			      feature<-feature-ref_scale$Mean[ref_scale$ID == names(feature)]
			      feature<-feature/ref_scale$SD[ref_scale$ID == names(feature)]
			      feature<-round(feature,3)
			      feature_tmp<-cbind(feature_tmp, feature)
			    }
			  feature_tmp
			  }
			}

			# Remove any columns containing NA
			TARG_expr<-TARG_expr[,which(unlist(lapply(TARG_expr, function(x)!is.na(x[1])))),with=F]

			# Save TARG_expr and compress
			fwrite(TARG_expr, paste0(opt$output,'/TARG_PROFILE_FILES/chr',chr,'/',panel,'_expression.txt'), nThread = opt$n_cores, sep=' ', na='NA')
			system(paste0(opt$pigz,' ',opt$output,'/TARG_PROFILE_FILES/chr',chr,'/',panel,'_expression.txt'))

			# Delete TARG_expr and run garbage collection.
			# Garbage collection turns out to be very important here as memory can spike in subsequent loop.
			rm(TARG_expr)
			rm(pos_chr_panel)
			gc()
		}

		system(paste0("rm ",opt$output,"/targ_chr",chr,"*"),ignore.stdout=T, ignore.stderr=T)

		sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
		cat('Done!\n')
		sink()

		rm(pos_chr)
		gc()
	}

	# Combine per chromosome predicted expression values and insert FID and IID columns
	sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
	cat('Combining per chromomsome files...')
	sink()
	for(panel in unique(pos$PANEL)){
		pos_panel<-pos[pos$PANEL == panel,]
		for(chr in unique(pos_panel$CHR)){
			system(paste0(opt$pigz,' -d ',opt$output,'/TARG_PROFILE_FILES/chr',chr,'/',panel,'_expression.txt.gz'))
		}
		system(paste0("paste -d ' ' $(echo ",opt$output,'/TARG_PROFILE_FILES/TARG.IDs $(ls ',opt$output,'/TARG_PROFILE_FILES/chr*/',panel,'_expression.txt)) | ',opt$pigz,' > ',opt$output,'/FeaturePredictions_',panel,'.txt.gz'))
		system(paste0('rm ',opt$output,'/TARG_PROFILE_FILES/chr*/',panel,'_expression.txt'))
	}
	sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
	cat('Done!\n')
	sink()

  system(paste0("rm -r ",opt$output,"/TARG_PROFILE_FILES"))
}

# Delete temporary files
system(paste0("rm ",opt$output,"/ref_indiv*"))
system(paste0("rm ",opt$output,"/ref_keep"))
system(paste0("rm ",opt$output,"/intersect.snplist"))

if(opt$score_files == paste0(opt$output,"/SCORE_FILES") & opt$save_score == F){
	system(paste0("rm -r ",opt$output,"/SCORE_FILES"))
}

if(is.na(opt$ref_expr)){
	system(paste0("rm -r ",opt$output,"/REF_PROFILE_FILES"))
}

end.time <- Sys.time()
time.taken <- end.time - start.time

sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()

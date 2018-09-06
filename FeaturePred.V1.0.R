#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at Cardiff University under the supervision of Richard Anney.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--PLINK_prefix", action="store", default=NA, type='character',
		help="Path to PLINK binaries [required]"),
make_option("--PLINK_prefix_chr", action="store", default=NA, type='character',
		help="Path to per chromosome PLINK binaries [required]"),
make_option("--weights", action="store", default=NA, type='character',
		help="Path for .pos file describing features [required]"),
make_option("--weights_dir", action="store", default=NA, type='character',
		help="Directory containing the weights corresponding to the features in the .pos file [required]"),
make_option("--ref_ld_chr", action="store", default=NA, type='character',
		help="Path to FUSION 1KG reference [required]"),
make_option("--make_score_script", action="store", default=NA, type='character',
		help="Path 'make_score.R' script from FUSION [required]"),
make_option("--n_cores", action="store", default=1, type='numeric',
		help="Specify the number of cores available [required]"),
make_option("--memory", action="store", default=2000, type='numeric',
		help="RAM available in MB [required]"),
make_option("--plink", action="store", default='NA', type='character',
		help="Path to PLINK software [required]"),
make_option("--output", action="store", default=NA, type='character',
		help="Name of output directory [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

if(file.exists(paste(opt$output,'/FeaturePredictions.csv',sep=''))){
	cat('Error: A file named',paste0(opt$output,'/FeaturePredictions.csv'),'already exists.\n')
	q()
}

system(paste('mkdir -p ',opt$output))

sink(file = paste(opt$output,'/FeaturePredictions.log',sep=''), append = F)
cat(
'#################################################################
# FeaturePred V1.0
# V1.0 05/09/2018
#################################################################

Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')

if(is.na(opt$PLINK_prefix) & is.na(opt$PLINK_prefix_chr)){
	cat('Error: No target sample PLINK files have been specified\n')
	q()
}
if(!is.na(opt$PLINK_prefix) & !is.na(opt$PLINK_prefix_chr)){
	cat('Error: Both PLINK_prefix and PLINK_prefix_chr have been specified.\n')
	q()
}

sink()

suppressMessages(library(data.table))
suppressMessages(library(foreach))
suppressMessages(library(doMC))
registerDoMC(opt$n_cores)

###################################
# Check SNP overlap between FUSION reference and target dataset
###################################

ref_ld_chr_files<-sub('.*/', '', opt$ref_ld_chr)
ref_ld_chr_dir<-sub(ref_ld_chr_files, '', opt$ref_ld_chr)
temp = list.files(path=ref_ld_chr_dir, pattern=paste0(ref_ld_chr_files,'*.bim'))

Ref<-do.call(rbind, lapply(paste0(ref_ld_chr_dir,temp), function(x) data.frame(fread(x))))

# Read in the SNPs in CLOZUK
if(!is.na(opt$PLINK_prefix)){
	Target<-data.frame(fread(paste(opt$PLINK_prefix,'.bim',sep='')))
}

if(!is.na(opt$PLINK_prefix_chr)){
	PLINK_prefix_chr_files<-sub('.*/', '', opt$PLINK_prefix_chr)
	PLINK_prefix_chr_dir<-sub(PLINK_prefix_chr_files, '', opt$PLINK_prefix_chr)
	temp = list.files(path=PLINK_prefix_chr_dir, pattern=paste0(PLINK_prefix_chr_files,'*.bim'))
	Target<-do.call(rbind, lapply(paste0(ref_ld_chr_dir,temp), function(x) data.frame(fread(x))))
}

# Get intersect of the two based on RSID
Overlap<-intersect(Ref$V2, Target$V2)

sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
cat('Number of SNPs in FUSION LD Reference = ',length(Ref$V2),'
Number of SNPs in target PLINK files = ',length(Target$V2),'
Number of SNPs in both = ',length(Overlap),'
Percentage of SNPs in FUSION LD Reference that are in the target PLINK files = ',length(Overlap)/length(Ref$V2)*100,'%\n',sep='')
sink()

# Save list of SNPs that intersect the target and ref
fwrite(list(Overlap),paste0(opt$output,'/intersect.snplist'),col.names=F)

sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
cat('Extracting intersect from target sample...',sep='')
sink()

# Run PLINK to extract intersecting SNPs
system(paste0(opt$plink,' --bfile ', opt$PLINK_prefix,' --make-bed --extract ',opt$output,'/intersect.snplist --out ',opt$output,'/intersect_target --memory ', floor(opt$memory*.9)),ignore.stdout=T, ignore.stderr=T)

sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
cat('Done!\n',sep='')
sink()

# Create new opt$PLINK_prefix_new variable for the intersect files.
opt$PLINK_prefix_new<-paste0(opt$output,'/intersect_target')

######################################
# Flip variants in target to match the FUSION 1KG reference
######################################

# Extract intersecting SNPs and put in the same order.
Ref<-Ref[(Ref$V2 %in% Overlap),]
Target<-Target[(Target$V2 %in% Overlap),]
Target<-Target[match(Ref$V2, Target$V2),]

# Find SNPs with mismatching allele codes.
mismatch<-Target[which(Target$V5 != Ref$V5 & Target$V5 != Ref$V6),]

sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
cat('Number of SNPs that have mismatched allele codes = ',dim(mismatch)[1],'\n',sep='')
sink()

if(dim(mismatch)[1] != 0){
	# Save list of SNPs that need to flipped
	fwrite(mismatch['V2'],paste0(opt$output,'/mismatch.snplist'),col.names=F)
	
	sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
	cat('Flipping mismatch allele in target sample...',sep='')
	sink()
	
	# Run PLINK to flip mismatch alleles.
	system(paste0(opt$plink,' --bfile ', opt$PLINK_prefix_new,' --make-bed --flip ',opt$output,'/mismatch.snplist --out ',opt$output,'/intersect_flipped_target --memory ', floor(opt$memory*.9)),ignore.stdout=T, ignore.stderr=T)
	
	sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
	cat('Done!\n',sep='')
	sink()
	
	# Delete intersect file
	system(paste0('rm ',opt$PLINK_prefix_new,'.bed'))
	system(paste0('rm ',opt$PLINK_prefix_new,'.bim'))
	system(paste0('rm ',opt$PLINK_prefix_new,'.fam'))
	
	# Reassign opt$PLINK_prefix_new to the flipped files.
	opt$PLINK_prefix_new<-paste0(opt$output,'/intersect_flipped_target')
	
	# Check whether flipping SNPs reduced the number of SNPs with mismatched allele codes
	if(!is.na(opt$PLINK_prefix_new)){
		Target<-data.frame(fread(paste(opt$PLINK_prefix_new,'.bim',sep='')))
	}
	
	if(!is.na(opt$PLINK_prefix_chr)){
		PLINK_prefix_chr_files<-sub('.*/', '', opt$PLINK_prefix_chr)
		PLINK_prefix_chr_dir<-sub(PLINK_prefix_chr_files, '', opt$PLINK_prefix_chr)
		temp = list.files(path=PLINK_prefix_chr_dir, pattern=paste0(PLINK_prefix_chr_files,'*.bim'))
		Target<-do.call(rbind, lapply(paste0(ref_ld_chr_dir,temp), function(x) data.frame(fread(x))))
	}
	
	Target<-Target[match(Ref$V2, Target$V2),]
	mismatch2<-Target[which(Target$V5 != Ref$V5 & Target$V5 != Ref$V6),]
	
	sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
	cat('Number of SNPs that have mismatched allele codes after flipping= ',length(mismatch2),'\n',sep='')
	sink()
}

######################################
# Convert weights into SCORE files for PLINK
######################################

# Read in the .pos file
pos<-data.frame(fread(opt$weights))

sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
cat('The .pos file contains ',dim(pos)[1],' features.\n',sep='')
sink()

# Attach weights directory to WGT values in pos file
pos$FILE<-paste0(opt$weights_dir,'/',sub('.*/','',pos$WGT))

# Remove .wgt.RDat from the WGT values
pos$WGT<-gsub('.wgt.RDat','',pos$WGT)
pos$WGT<-gsub('.*/','',pos$WGT)

# Create directory for the SCORE files
system(paste0('mkdir ',opt$output,'/SCORE_FILES'))

# Create SCORE file for each set of weights using FUSION make_score.R
sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
cat('Converting weights files into PLINK SCORE files...',sep='')
sink()

error<-foreach(i=1:length(pos$FILE), .combine=c) %dopar% {
	system(paste0('Rscript ',opt$make_score_script,' ',pos$FILE[i],' > ', opt$output,'/SCORE_FILES/',pos$WGT[i],'.SCORE'))
}

sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
cat('Done!\n',sep='')
sink()

sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
cat(sum(error),' errors were encountered when creating SCORE files\n',sep='')
sink()

if(sum(error) != 0){
}

#######################################
# Predict features in target sample using PLINK
#######################################

# Create directory for the SCORE files
system(paste0('mkdir ',opt$output,'/PROFILE_FILES'))

# Calculate profile scores (i.e. feature predictions)
sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
cat('Predicting features in target sample...',sep='')
sink()

error<-foreach(i=1:length(pos$FILE), .combine=c) %dopar% {
	system(paste0(opt$plink,' --bfile ',opt$PLINK_prefix_new,' --score ',opt$output,'/SCORE_FILES/',pos$WGT[i],'.SCORE 1 2 4 --out ',opt$output,'/PROFILE_FILES/',pos$WGT[i],' --memory ', floor((opt$memory*0.8)/opt$n_cores+2)),ignore.stdout=T, ignore.stderr=T)
}

sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
cat('Done!\n',sep='')
sink()

sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
cat(sum(error),' errors were encountered during feature prediction\n',sep='')
sink()

if(sum(error) != 0){
}

#######################################
# Combine .profile files into a single file and write out
#######################################

# List all profile files
files<-list.files(path=paste0(opt$output,'/PROFILE_FILES/'),pattern='*profile')
files<-gsub('.profile', '', files)

fam<-read.table(paste0(opt$PLINK_prefix_new,'.fam'), header=F, stringsAsFactors=F)
IDs<-fam[c(1,2)]
names(IDs)<-c('FID','IID')

AllFeat<-data.frame(IDs,foreach(i=1:length(files), .combine=cbind) %dopar% {
	Feat<-read.table(paste0(opt$output,'/PROFILE_FILES/',files[i],'.profile'), header=T, stringsAsFactors=F)
	names(Feat)[6]<-files[i]
	Feat[6]
})

fwrite(AllFeat, paste0(opt$output,'/FeaturePredictions.csv'), sep=',')

sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
cat((dim(AllFeat)[2]-2),' features have been predicted in the target sample\n',sep='')
cat('Feature predictions have been saved as ',opt$output,'/FeaturePredictions.csv\n',sep='')
sink()

#########################################
# Delete temporary files
#########################################

system(paste0('rm ', opt$output,'/intersect*'))
system(paste0('rm -r ', opt$output,'/SCORE_FILES'))
system(paste0('rm -r ', opt$output,'/PROFILE_FILES'))

#########################################
# Write a summary to the log file
#########################################
end.time <- Sys.time()
time.taken <- end.time - start.time

sink(file = paste0(opt$output,'/FeaturePredictions.log'), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()

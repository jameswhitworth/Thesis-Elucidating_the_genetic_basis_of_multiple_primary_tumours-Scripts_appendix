###############################################################################################################################
##Scripts to count single nucleotide variant and indel calls affecting elements of interest and compare frequency vs controls##
###############################################################################################################################

########################################################################################################################################################################################################################
##Scripts contain functions that were not utilised in interpretation of results as follows:
#1.Counting of variants per gene/element of interest where the numerator was taken as the total number of variants observed in that gene in cases /controls and the denominator was taken as the number of cases/controls multiplied the number of variant sites observed
#2.Counts of individual variants or variants per gene/element where heterozygous and homozygous variants where aggregated (homozygous counting as 2 and heterozgous as 1) to consider the total frequency of variant alleles (denoted "hets_homs_agg") in cases or controls rather than the frequency of individual variants or individuals with variants in a given gene/element
#3.Case control comparison between cases and 1958 birth cohort controls, though these were retained in outputs in case supporting evidence of results of interest might have been gleaned
########################################################################################################################################################################################################################

##############################################################################################################################################################
##Filtering of merged VCF files compiled by NIHR BioResource Rare Diseases project and stored on University of Cambridge high performance computeing cluster##
##############################################################################################################################################################

##Produce list of per chromosome files for filtering
ls /scratch/WGS10K/data/release/20170614-A/merged-vcf/no_hgmd/*.vcf.gz > 1_chromosome_files.txt
sed -r 's/.{7}$//' 1_chromosome_files.txt > 2_chromosome_files.txt
sed -r 's/.{59}//' 2_chromosome_files.txt > chromosome_files.txt
rm 1_chromosome_files.txt
rm 2_chromosome_files.txt

##Filter per chromosome VCF files using BED file corresponding to coordinates of interest. BED file containing all coordinates for gene lists
for i in `cat chromosome_files.txt`; do
bcftools view -R TELgenes_exons.bed -e 'FILTER!="PASS"' -o ${i}_filtered.vcf -O v /scratch/WGS10K/data/release/20170614-A/merged-vcf/no_hgmd/${i}.vcf.gz
done

##Merge filtered per chromosome VCF files and retain lines corresponding to samples of interest
ls *_filtered.vcf > filtered_samples_for_bgzip_and_tabix.txt

for i in `cat filtered_samples_for_bgzip_and_tabix.txt`; do
	vcf-sort -c ${i} > sorted_${i}
	bgzip sorted_${i}
	tabix sorted_${i}.gz
done

rm filtered_samples_for_bgzip_and_tabix.txt

ls *.vcf.gz | tr "\n" " " > files_to_merge.txt

bcftools concat -Oz `cat files_to_merge.txt` -o merged.vcf.gz

##Filter merged VCF based on quality parameters (insert missing genotypes where criteria not fulfilled)
bcftools filter -e 'FMT/DP<10 || FMT/GQ<30 || (FMT/AD[1]/(FMT/AD[0] + FMT/AD[1]) <0.3) || (FMT/AD[1]/(FMT/AD[0] + FMT/AD[1] + FMT/AD[2]) <0.3) || (FMT/AD[2]/(FMT/AD[0] + FMT/AD[1] + FMT/AD[2]) <0.3)' -S . -o TELgenes_MPT_BRIDGE_merged.vcf.gz -O z merged.vcf.gz


##Transferred to local server

###################################################
##Further data addition, annotation and filtering##
###################################################

##Script launched with following commands

#./VEP_and_filter_for_candidate_search.sh TELgenes ALL XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX

######################################################################################
##Filtering of 1958 birth cohort VCF files compiled by Institute of Cancer Research ##
######################################################################################

##Produce list of files for filtering
ls /home/jww39/WES_1958BC_hg19_Individuals/*.vcf.gz > 1_samples.txt
sed -r 's/.{7}$//' 1_samples.txt > 2_samples.txt
sed -r 's/.{40}//' 2_samples.txt > samples.txt
rm 1_samples.txt
rm 2_samples.txt


## Filter merged VCF
for i in `cat samples.txt`; do
  bcftools view -R ${1}_exons.bed -o ${i}_filtered.vcf -O v /home/jww39/WES_1958BC_hg19_Individuals/${i}.vcf.gz
done

##Split multi-allelic sites
ls *_filtered.vcf > filtered_samples_for_bgzip_and_tabix.txt
for i in `cat filtered_samples_for_bgzip_and_tabix.txt`; do
	vcf-sort -c ${i} > sorted_${i}
	bgzip sorted_${i}
	tabix sorted_${i}.gz
done
rm filtered_samples_for_bgzip_and_tabix.txt
ls sorted_*hg19.QC_filtered.vcf.gz | sed 's/.\{7\}$//' > sorted_files_for_allele_split.txt
for i in `cat sorted_files_for_allele_split.txt`; do
  java -jar /data/Resources/Software/Javas/GenomeAnalysisTK.jar -T LeftAlignAndTrimVariants -R /data/Resources/References/hg19/hg19.fa --variant ${i}.vcf.gz -o ${i}_split.vcf.gz --splitMultiallelics
done

##Merge VCF files
ls sorted_*_hg19.QC_filtered_split.vcf.gz | tr "\n" " " > files_to_merge.txt
bcftools merge -m none `cat files_to_merge.txt` -o ICR_merged.vcf.gz -O z
tabix ICR_merged.vcf.gz

##Filter merged VCF based on quality parameters (insert missing genotypes where criteria not fulfilled)
bcftools filter -e 'FMT/DP<10 || FMT/GQ<30 || (FMT/AD[1]/(FMT/AD[0] + FMT/AD[1]) <0.3) || (FMT/AD[1]/(FMT/AD[0] + FMT/AD[1] + FMT/AD[2]) <0.3) || (FMT/AD[2]/(FMT/AD[0] + FMT/AD[1] + FMT/AD[2]) <0.3)' -S . -o ${1}_ICR_merged_filtered.vcf.gz -O z ICR_merged.vcf.gz
tabix ${1}_ICR_merged_filtered.vcf.gz

################################################################################
##Merge 1958 birth cohort and NIHR BioResource Rare Diseases merged VCF files ##
################################################################################

tabix ${1}_MPT_BRIDGE_merged.vcf.gz
bcftools merge -m none ${1}_ICR_merged_filtered.vcf.gz ${1}_MPT_BRIDGE_merged.vcf.gz -o ${1}_MPT_BRIDGE_ICR_merged.vcf.gz -O z

#############################
##Annotate merged VCF file ##
#############################

##Annotate with Variant Effect Predictor
/home/jww39/ensembl-vep/vep -i ${1}_MPT_BRIDGE_ICR_merged.vcf -offline --assembly GRCh37 -o VEP_out_${1}_MPT_BRIDGE_ICR_merged.vcf --pick --pick_order canonical --everything --vcf --plugin LoF --plugin CADD,/home/jww39/.vep/Plugins/whole_genome_SNVs.tsv.gz,/home/jww39/.vep/Plugins/InDels.tsv.gz --vcf_info_field ANN

##Filter with Variant Effect Predictor filter script based on annotations  
/home/jww39/ensembl-vep/filter_vep -i VEP_out_${1}_MPT_BRIDGE_ICR_merged.vcf -o VEP_out_filtered_${1}_MPT_BRIDGE_ICR_merged.vcf -filter "Consequence in /home/jww39/candidate_gene_searches/consequence.txt" --only_matched
/home/jww39/ensembl-vep/filter_vep -i VEP_out_filtered_${1}_MPT_BRIDGE_ICR_merged.vcf -o VEP_out_double_filtered_${1}_MPT_BRIDGE_ICR_merged.vcf -filter "(EUR_AF < 0.01 or not EUR_AF) and (UK10KWGS_AF < 0.01 or not UK10KWGS_AF) and (WGS10K_AF < 0.05 or not WGS10K_AF)" --only_matched
  
#Prepare resulting VCF for reading into R
sed -i 's/#CHROM/CHROM/g' VEP_out_double_filtered_${1}_MPT_BRIDGE_ICR_merged.vcf
sed '/^#/ d' VEP_out_double_filtered_${1}_MPT_BRIDGE_ICR_merged.vcf > VEP_out_double_filtered_${1}_MPT_BRIDGE_ICR_merged_without_header.vcf


######################################################################################################################
##R script to identify samples with each variant in table, count variants and compare frequency in cases vs controls##
######################################################################################################################

##Following commands used to launch script according to phenotypic subgroup:

#Rscript candidate_gene_output.R TELgenes ALL XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX

args = commandArgs(trailingOnly=TRUE)

##Read the vcf output in
vcf_to_read <- paste("VEP_out_double_filtered",args[1],"MPT_BRIDGE_ICR_merged.vcf",sep = "_") 
df <- read.table(vcf_to_read, header = T, comment.char = "")
df2 <- df[10:ncol(df)]

BRIDGE_list <- read.delim("BRIDGE_euro_controls.txt", header = F)
ICR_list <- read.delim("ICR_list.txt", header = F)

##Preparing a list of cases
if(args[2] == "ALL")case_list <- read.delim("high_resid_MPMT_samples.txt", header = F) ##File containing sample names with high telomere length estimate. Script run again for those with low estimate

##Identifying samples containing each variant
samples_raw_hets <- data.frame()
for (i in 1:nrow(df2)) {
  
  a <- unlist(df2[i,])
  b <- grep("^0/1", a, perl = TRUE)
  c <- as.vector(names(a[b]))
  d <- paste(c, collapse = ';')
  samples_raw_hets <- append(samples_raw_hets, d)
} 
samples_hets <- unlist(samples_raw_hets)

samples_raw_homs <- data.frame()
for (i in 1:nrow(df2)) {
  
  a <- unlist(df2[i,])
  b <- grep("^1/1", a, perl = TRUE)
  c <- as.vector(names(a[b]))
  d <- paste(c, collapse = ';')
  samples_raw_homs <- append(samples_raw_homs, d)
  
} 
samples_homs <- unlist(samples_raw_homs)
 
##Extracting which samples are among the cases (hets)
samples_df_hets <- as.data.frame(samples_hets)

cases_list_hets <- data.frame()
for (i in 1:nrow(samples_df_hets)) {
  
  cases_split_hets <- strsplit(samples_hets[i], ";")
  cases_split_hets_unlisted <- unlist(cases_split_hets)
  case_samples_hets <- cases_split_hets_unlisted[(which(cases_split_hets_unlisted %in% case_list[,1]))] #Gives items that match the list
  case_samples_hets_string <- paste(case_samples_hets, collapse = ";")
  cases_list_hets <- as.vector(append(cases_list_hets, case_samples_hets_string))
  
}
cases_hets <- unlist(cases_list_hets)

##Counting how many cases have the variant (hets)
library(stringr)
cases_df_hets <- as.data.frame(cases_hets)


df_for_cases_variant_count_hets <- data.frame()
for (i in 1:nrow(cases_df_hets)) { 
  count_of_samples_with_variant_hets <- (str_count(cases_df_hets[i,], "[A-Z]0"))
  df_for_cases_variant_count_hets <- as.vector(append(df_for_cases_variant_count_hets, count_of_samples_with_variant_hets))
}
cases_variant_count_hets <- cbind(df_for_cases_variant_count_hets)
cases_variant_count_hets <- unlist(cases_variant_count_hets)

##Extracting BRIDGE samples (hets)

BRIDGES_list_hets <- data.frame()
for (i in 1:nrow(samples_df_hets)) {
  
  BRIDGES_split_hets <- strsplit(samples_hets[i], ";")
  BRIDGES_split_hets_unlisted <- unlist(BRIDGES_split_hets)
  BRIDGE_samples_hets <- BRIDGES_split_hets_unlisted[(which(BRIDGES_split_hets_unlisted %in% BRIDGE_list[,1]))] #Gives items that match the list
  BRIDGE_samples_hets_string <- paste(BRIDGE_samples_hets, collapse = ";")
  BRIDGES_list_hets <- as.vector(append(BRIDGES_list_hets, BRIDGE_samples_hets_string))
  
}
BRIDGES_hets <- unlist(BRIDGES_list_hets)

##Counting how many BRIDGE samples have the variant (hets)


library(stringr)
BRIDGES_df_hets <- as.data.frame(BRIDGES_hets)
df_for_BRIDGES_variant_count_hets <- data.frame()
for (i in 1:nrow(BRIDGES_df_hets)) { 
  count_of_samples_with_variant_hets <- (str_count(BRIDGES_df_hets[i,], "[A-Z]0"))
  df_for_BRIDGES_variant_count_hets <- as.vector(append(df_for_BRIDGES_variant_count_hets, count_of_samples_with_variant_hets))
}
BRIDGES_variant_count_hets <- cbind(df_for_BRIDGES_variant_count_hets)
BRIDGES_variant_count_hets <- unlist(BRIDGES_variant_count_hets)

##Extracting 1958BC samples (hets)


ICRS_list_hets <- data.frame()
for (i in 1:nrow(samples_df_hets)) {
  
  ICRS_split_hets <- strsplit(samples_hets[i], ";")
  ICRS_split_hets_unlisted <- unlist(ICRS_split_hets)
  ICR_samples_hets <- ICRS_split_hets_unlisted[(which(ICRS_split_hets_unlisted %in% ICR_list[,1]))] #Gives items that match the list
  ICR_samples_hets_string <- paste(ICR_samples_hets, collapse = ";")
  ICRS_list_hets <- as.vector(append(ICRS_list_hets, ICR_samples_hets_string))
  
}
ICRS_hets <- unlist(ICRS_list_hets)

##Counting how many 1958BC samples have the variant (hets)
library(stringr)
ICRS_df_hets <- as.data.frame(ICRS_hets)
df_for_ICRS_variant_count_hets <- data.frame()
for (i in 1:nrow(ICRS_df_hets)) { 
  count_of_samples_with_variant_hets <- (str_count(ICRS_df_hets[i,], "EGAR0")) 
  df_for_ICRS_variant_count_hets <- as.vector(append(df_for_ICRS_variant_count_hets, count_of_samples_with_variant_hets))
}

ICRS_variant_count_hets <- cbind(df_for_ICRS_variant_count_hets)
ICRS_variant_count_hets <- unlist(ICRS_variant_count_hets)

##Extracting which samples are among the cases (homs)
samples_df_homs <- as.data.frame(samples_homs)
cases_list_homs <- data.frame()
for (i in 1:nrow(samples_df_homs)) {
  
  cases_split_homs <- strsplit(samples_homs[i], ";")
  cases_split_homs_unlisted <- unlist(cases_split_homs)
  case_samples_homs <- cases_split_homs_unlisted[(which(cases_split_homs_unlisted %in% case_list[,1]))] #Gives items that match the list
  case_samples_homs_string <- paste(case_samples_homs, collapse = ";")
  cases_list_homs <- as.vector(append(cases_list_homs, case_samples_homs_string))
  
}

cases_homs <- unlist(cases_list_homs)

##Counting how many cases have the variant (homs)
library(stringr)
cases_df_homs <- as.data.frame(cases_homs)
df_for_cases_variant_count_homs <- data.frame()
for (i in 1:nrow(cases_df_homs)) { 
  count_of_samples_with_variant_homs <- (str_count(cases_df_homs[i,], "[A-Z]0"))
  df_for_cases_variant_count_homs <- as.vector(append(df_for_cases_variant_count_homs, count_of_samples_with_variant_homs))
}

cases_variant_count_homs <- cbind(df_for_cases_variant_count_homs)
cases_variant_count_homs <- unlist(cases_variant_count_homs)

##Extracting BRIDGE samples (homs)
BRIDGES_list_homs <- data.frame()
for (i in 1:nrow(samples_df_homs)) {
  
  BRIDGES_split_homs <- strsplit(samples_homs[i], ";")
  BRIDGES_split_homs_unlisted <- unlist(BRIDGES_split_homs)
  BRIDGE_samples_homs <- BRIDGES_split_homs_unlisted[(which(BRIDGES_split_homs_unlisted %in% BRIDGE_list[,1]))] #Gives items that match the list
  BRIDGE_samples_homs_string <- paste(BRIDGE_samples_homs, collapse = ";")
  BRIDGES_list_homs <- as.vector(append(BRIDGES_list_homs, BRIDGE_samples_homs_string))
  
}

BRIDGES_homs <- unlist(BRIDGES_list_homs)

##Counting how many BRIDGE samples have the variant (homs)
library(stringr)
BRIDGES_df_homs <- as.data.frame(BRIDGES_homs)
df_for_BRIDGES_variant_count_homs <- data.frame()
for (i in 1:nrow(BRIDGES_df_homs)) { 
  count_of_samples_with_variant_homs <- (str_count(BRIDGES_df_homs[i,], "[A-Z]0"))
  df_for_BRIDGES_variant_count_homs <- as.vector(append(df_for_BRIDGES_variant_count_homs, count_of_samples_with_variant_homs))
}

BRIDGES_variant_count_homs <- cbind(df_for_BRIDGES_variant_count_homs)
BRIDGES_variant_count_homs <- unlist(BRIDGES_variant_count_homs)

##Extracting 1958BC samples (homs)
ICRS_list_homs <- data.frame()
for (i in 1:nrow(samples_df_homs)) {
  
  ICRS_split_homs <- strsplit(samples_homs[i], ";")
  ICRS_split_homs_unlisted <- unlist(ICRS_split_homs)
  ICR_samples_homs <- ICRS_split_homs_unlisted[(which(ICRS_split_homs_unlisted %in% ICR_list[,1]))] #Gives items that match the list
  ICR_samples_homs_string <- paste(ICR_samples_homs, collapse = ";")
  ICRS_list_homs <- as.vector(append(ICRS_list_homs, ICR_samples_homs_string))
  
}

ICRS_homs <- unlist(ICRS_list_homs)

##Counting how many 1958BC samples have the variant (homs)
library(stringr)
ICRS_df_homs <- as.data.frame(ICRS_homs)
df_for_ICRS_variant_count_homs <- data.frame()
for (i in 1:nrow(ICRS_df_homs)) { 
  count_of_samples_with_variant_homs <- (str_count(ICRS_df_homs[i,], "EGAR0")) 
  df_for_ICRS_variant_count_homs <- as.vector(append(df_for_ICRS_variant_count_homs, count_of_samples_with_variant_homs))
}

ICRS_variant_count_homs <- cbind(df_for_ICRS_variant_count_homs)
ICRS_variant_count_homs <- unlist(ICRS_variant_count_homs)

##Produce columns containing variant counts 
cases_variant_count_hets_homs <- cases_variant_count_hets + cases_variant_count_homs
BRIDGES_variant_count_hets_homs <- BRIDGES_variant_count_hets + BRIDGES_variant_count_homs
ICRS_variant_count_hets_homs <- ICRS_variant_count_hets + ICRS_variant_count_homs

cases_variant_count_hets_homs_agg <- cases_variant_count_hets + (2 * cases_variant_count_homs)
BRIDGES_variant_count_hets_homs_agg <- BRIDGES_variant_count_hets + (2 * BRIDGES_variant_count_homs)
ICRS_variant_count_hets_homs_agg <- ICRS_variant_count_hets + (2 * ICRS_variant_count_homs)


##Produce a table with variant information and counts
INFO_col_character <- as.character(df$INFO)
library(stringr)
gene_col <- str_match(INFO_col_character , "ENSG\\d{11}")

MPT_BRIDGE_ICR_counts <- cbind(
  df[,1:9],
  cases_variant_count_hets, #10
  cases_variant_count_homs, #11
  cases_variant_count_hets_homs, #12
  cases_variant_count_hets_homs_agg, #13
  BRIDGES_variant_count_hets, #14
  BRIDGES_variant_count_homs, #15
  BRIDGES_variant_count_hets_homs, #16
  BRIDGES_variant_count_hets_homs_agg, #17
  ICRS_variant_count_hets, #18
  ICRS_variant_count_homs, #19
  ICRS_variant_count_hets_homs, #20
  ICRS_variant_count_hets_homs_agg, #21
  gene_col, #22
  cases_hets,
  BRIDGES_hets,
  ICRS_hets,
  cases_homs,
  BRIDGES_homs,
  ICRS_homs,
  df[,10:ncol(df)] #23 onwards
)

############################
##Case control comparisons##
############################

##Fishers exact with fdr adjustment for each individual variant observed (hets cases vs BRIDGE)
rownum <- nrow(MPT_BRIDGE_ICR_counts)
BRIDGE_pvalue_hets <- as.vector(c(), mode = "any")
for(i in 1:rownum){
  
  fishtable = matrix(c(MPT_BRIDGE_ICR_counts[i,10],(nrow(case_list) - MPT_BRIDGE_ICR_counts[i,10]),MPT_BRIDGE_ICR_counts[i,14],(nrow(BRIDGE_list) - MPT_BRIDGE_ICR_counts[i,14])), 
                     nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue_hets <- append(BRIDGE_pvalue_hets, fisherout[[1]])
  
}
BRIDGE_qvalue_hets <- p.adjust(BRIDGE_pvalue_hets, method = "fdr", n = rownum)

##Fishers exact with fdr adjustment for each individual variant observed (hets cases vs ICR)
rownum <- nrow(MPT_BRIDGE_ICR_counts)
ICR_pvalue_hets <- as.vector(c(), mode = "any")
for(i in 1:rownum){
  
  fishtable = matrix(c(MPT_BRIDGE_ICR_counts[i,10],(nrow(case_list) - MPT_BRIDGE_ICR_counts[i,10]),MPT_BRIDGE_ICR_counts[i,18],(nrow(BRIDGE_list) - MPT_BRIDGE_ICR_counts[i,18])), 
                     nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  ICR_pvalue_hets <- append(ICR_pvalue_hets, fisherout[[1]])
  
}
ICR_qvalue_hets <- p.adjust(ICR_pvalue_hets, method = "fdr", n = rownum)

##Fishers exact with fdr adjustment for each individual variant observed (homs cases vs BRIDGE)
rownum <- nrow(MPT_BRIDGE_ICR_counts)
BRIDGE_pvalue_homs <- as.vector(c(), mode = "any")
for(i in 1:rownum){
  
  fishtable = matrix(c(MPT_BRIDGE_ICR_counts[i,11],(nrow(case_list) - MPT_BRIDGE_ICR_counts[i,11]),MPT_BRIDGE_ICR_counts[i,15],(nrow(BRIDGE_list) - MPT_BRIDGE_ICR_counts[i,15])), 
                     nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue_homs <- append(BRIDGE_pvalue_homs, fisherout[[1]])
  
}
BRIDGE_qvalue_homs <- p.adjust(BRIDGE_pvalue_homs, method = "fdr", n = rownum)

##Fishers exact with fdr adjustment for each individual variant observed (homs cases vs ICR)
rownum <- nrow(MPT_BRIDGE_ICR_counts)
ICR_pvalue_homs <- as.vector(c(), mode = "any")
for(i in 1:rownum){
  
  fishtable = matrix(c(MPT_BRIDGE_ICR_counts[i,11],(nrow(case_list) - MPT_BRIDGE_ICR_counts[i,11]),MPT_BRIDGE_ICR_counts[i,19],(nrow(BRIDGE_list) - MPT_BRIDGE_ICR_counts[i,19])), 
                     nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  ICR_pvalue_homs <- append(ICR_pvalue_homs, fisherout[[1]])
  
}
ICR_qvalue_homs <- p.adjust(ICR_pvalue_homs, method = "fdr", n = rownum)

##Fishers exact with fdr adjustment for each individual variant observed (hets_homs cases vs BRIDGE)
rownum <- nrow(MPT_BRIDGE_ICR_counts)
BRIDGE_pvalue_hets_homs <- as.vector(c(), mode = "any")
for(i in 1:rownum){
  
  fishtable = matrix(c(MPT_BRIDGE_ICR_counts[i,12],(nrow(case_list) - MPT_BRIDGE_ICR_counts[i,12]),MPT_BRIDGE_ICR_counts[i,16],(nrow(BRIDGE_list) - MPT_BRIDGE_ICR_counts[i,16])), 
                     nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue_hets_homs <- append(BRIDGE_pvalue_hets_homs, fisherout[[1]])
  
}
BRIDGE_qvalue_hets_homs <- p.adjust(BRIDGE_pvalue_hets_homs, method = "fdr", n = rownum)

##Fishers exact with fdr adjustment for each individual variant observed (hets_homs cases vs ICR)
rownum <- nrow(MPT_BRIDGE_ICR_counts)
ICR_pvalue_hets_homs <- as.vector(c(), mode = "any")
for(i in 1:rownum){
  
  fishtable = matrix(c(MPT_BRIDGE_ICR_counts[i,12],(nrow(case_list) - MPT_BRIDGE_ICR_counts[i,12]),MPT_BRIDGE_ICR_counts[i,20],(nrow(BRIDGE_list) - MPT_BRIDGE_ICR_counts[i,20])), 
                     nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  ICR_pvalue_hets_homs <- append(ICR_pvalue_hets_homs, fisherout[[1]])
  
}
ICR_qvalue_hets_homs <- p.adjust(ICR_pvalue_hets_homs, method = "fdr", n = rownum)

##Fishers exact with fdr adjustment for each individual variant observed (hets_homs_agg cases vs BRIDGE)
rownum <- nrow(MPT_BRIDGE_ICR_counts)
BRIDGE_pvalue_hets_homs_agg <- as.vector(c(), mode = "any")
for(i in 1:rownum){
  
  fishtable = matrix(c(MPT_BRIDGE_ICR_counts[i,13],((nrow(case_list)*2) - MPT_BRIDGE_ICR_counts[i,13]),MPT_BRIDGE_ICR_counts[i,17],((nrow(BRIDGE_list)*2) - MPT_BRIDGE_ICR_counts[i,17])), 
                     nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue_hets_homs_agg <- append(BRIDGE_pvalue_hets_homs_agg, fisherout[[1]])
  
}
BRIDGE_qvalue_hets_homs_agg <- p.adjust(BRIDGE_pvalue_hets_homs_agg, method = "fdr", n = rownum)

##Fishers exact with fdr adjustment for each individual variant observed (hets_homs_agg cases vs ICR)
rownum <- nrow(MPT_BRIDGE_ICR_counts)
ICR_pvalue_hets_homs_agg <- as.vector(c(), mode = "any")
for(i in 1:rownum){
  
  fishtable = matrix(c(MPT_BRIDGE_ICR_counts[i,13],((nrow(case_list)*2) - MPT_BRIDGE_ICR_counts[i,13]),MPT_BRIDGE_ICR_counts[i,21],((nrow(ICR_list)*2) - MPT_BRIDGE_ICR_counts[i,21])), 
                     nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  ICR_pvalue_hets_homs_agg <- append(ICR_pvalue_hets_homs_agg, fisherout[[1]])
  
}
ICR_qvalue_hets_homs_agg <- p.adjust(ICR_pvalue_hets_homs_agg, method = "fdr", n = rownum)

##Collation of results for individual variants
MPT_BRIDGE_ICR_final <- cbind(BRIDGE_pvalue_hets,
                              BRIDGE_qvalue_hets,
                              ICR_pvalue_hets,
                              ICR_qvalue_hets,
                              BRIDGE_pvalue_homs,
                              BRIDGE_qvalue_homs,
                              ICR_pvalue_homs,
                              ICR_qvalue_homs,
                              BRIDGE_pvalue_hets_homs,
                              BRIDGE_qvalue_hets_homs,
                              ICR_pvalue_hets_homs,
                              ICR_qvalue_hets_homs,
                              BRIDGE_pvalue_hets_homs_agg,
                              BRIDGE_qvalue_hets_homs_agg,
                              ICR_pvalue_hets_homs_agg,
                              ICR_qvalue_hets_homs_agg,
                              MPT_BRIDGE_ICR_counts)

##Output potentially significant results
MPT_BRIDGE_ICR_final$BRIDGE_qvalue_hets <- as.numeric(as.character(MPT_BRIDGE_ICR_final$BRIDGE_qvalue_hets))
MPT_BRIDGE_ICR_final$BRIDGE_qvalue_homs <- as.numeric(as.character(MPT_BRIDGE_ICR_final$BRIDGE_qvalue_homs))
MPT_BRIDGE_ICR_final$BRIDGE_qvalue_hets_homs <- as.numeric(as.character(MPT_BRIDGE_ICR_final$BRIDGE_qvalue_hets_homs))
MPT_BRIDGE_ICR_final$BRIDGE_qvalue_hets_homs_agg <- as.numeric(as.character(MPT_BRIDGE_ICR_final$BRIDGE_qvalue_hets_homs_agg))

MPT_BRIDGE_ICR_final$ICR_qvalue_hets <- as.numeric(as.character(MPT_BRIDGE_ICR_final$ICR_qvalue_hets))
MPT_BRIDGE_ICR_final$ICR_qvalue_homs <- as.numeric(as.character(MPT_BRIDGE_ICR_final$ICR_qvalue_homs))
MPT_BRIDGE_ICR_final$ICR_qvalue_hets_homs <- as.numeric(as.character(MPT_BRIDGE_ICR_final$ICR_qvalue_hets_homs))
MPT_BRIDGE_ICR_final$ICR_qvalue_hets_homs_agg <- as.numeric(as.character(MPT_BRIDGE_ICR_final$ICR_qvalue_hets_homs_agg))

sig_index <- which(
  (
    (BRIDGE_qvalue_hets < 0.05) 
    |
      (BRIDGE_qvalue_homs < 0.05) 
    |
      (BRIDGE_qvalue_hets_homs < 0.05) 
    |
      (BRIDGE_qvalue_hets_homs_agg < 0.05)  
  )
)
sig_index <- unique(sig_index)

variants_table_sig <- MPT_BRIDGE_ICR_final[sig_index,]

##Output results
write.csv(variants_table_sig, paste(args[1], "variant_table.csv", sep = "_"))
write.csv(MPT_BRIDGE_ICR_final, paste(args[1], "variant_table_all.csv", sep = "_"))

##Per gene analysis
genes <- unique(MPT_BRIDGE_ICR_final$gene_col)
genes <- genes[!is.na(genes)]

##Counts of variants per gene for cases, BRIDGE and ICR (hets)
variant_per_gene_count_cases_hets <- data.frame()

for(i in 1:length(genes)){
  a <- sum(MPT_BRIDGE_ICR_counts[,10][which(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i])]) 
  variant_per_gene_count_cases_hets  <- append(variant_per_gene_count_cases_hets,a)
}

variant_per_gene_count_cases_hets <- t(data.frame(variant_per_gene_count_cases_hets))
gene_names <- as.vector(genes)
rownames(variant_per_gene_count_cases_hets) <- gene_names

variant_per_gene_count_BRIDGES_hets <- data.frame()

for(i in 1:length(genes)){
  a <- sum(MPT_BRIDGE_ICR_counts[,14][which(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i])]) 
  variant_per_gene_count_BRIDGES_hets  <- append(variant_per_gene_count_BRIDGES_hets,a)
}

variant_per_gene_count_BRIDGES_hets <- t(data.frame(variant_per_gene_count_BRIDGES_hets))
gene_names <- as.vector(genes)
rownames(variant_per_gene_count_BRIDGES_hets) <- gene_names

variant_per_gene_count_ICRS_hets <- data.frame()

for(i in 1:length(genes)){
  a <- sum(MPT_BRIDGE_ICR_counts[,18][which(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i])]) 
  variant_per_gene_count_ICRS_hets  <- append(variant_per_gene_count_ICRS_hets,a)
}

variant_per_gene_count_ICRS_hets <- t(data.frame(variant_per_gene_count_ICRS_hets))
gene_names <- as.vector(genes)
rownames(variant_per_gene_count_ICRS_hets) <- gene_names

##Fishers exact with fdr adjustment for counts of variants per gene (hets)
cases_trials_per_gene <- as.vector(c(), mode = "any")
BRIDGE_pvalue_per_gene_hets <- as.vector(c(), mode = "any")
BRIDGE_trials_per_gene <- as.vector(c(), mode = "any")

for(i in 1:nrow(variant_per_gene_count_cases_hets))
{
  fishtable = matrix(
    c(variant_per_gene_count_cases_hets[i,1],
      sum(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i]) * nrow(case_list) - variant_per_gene_count_cases_hets[i,1],
      variant_per_gene_count_BRIDGES_hets[i,1],
      sum(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i]) * nrow(BRIDGE_list) - variant_per_gene_count_BRIDGES_hets[i,1]),
    
    nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue_per_gene_hets <- append(BRIDGE_pvalue_per_gene_hets, fisherout[[1]])
  cases_trials_per_gene <- append(cases_trials_per_gene, sum(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i]) * nrow(case_list))
  BRIDGE_trials_per_gene <- append(BRIDGE_trials_per_gene, sum(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i]) * nrow(BRIDGE_list))
}

BRIDGE_qvalue_per_gene_hets <- p.adjust(BRIDGE_pvalue_per_gene_hets, method = "fdr", n = 137)
cases_per_gene_prop_hets <- variant_per_gene_count_cases_hets / cases_trials_per_gene
BRIDGES_per_gene_prop_hets <- variant_per_gene_count_BRIDGES_hets / BRIDGE_trials_per_gene

ICR_pvalue_per_gene_hets <- as.vector(c(), mode = "any")
ICR_trials_per_gene <- as.vector(c(), mode = "any")

for(i in 1:nrow(variant_per_gene_count_cases_hets))
{
  fishtable = matrix(
    c(variant_per_gene_count_cases_hets[i,1],
      sum(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i]) * nrow(case_list) - variant_per_gene_count_cases_hets[i,1],
      variant_per_gene_count_BRIDGES_hets[i,1],
      sum(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i]) * nrow(ICR_list) - variant_per_gene_count_ICRS_hets[i,1]),
    
    nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  ICR_pvalue_per_gene_hets <- append(ICR_pvalue_per_gene_hets, fisherout[[1]])
  ICR_trials_per_gene <- append(ICR_trials_per_gene, sum(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i]) * nrow(ICR_list))
}
ICR_qvalue_per_gene_hets <- p.adjust(ICR_pvalue_per_gene_hets, method = "fdr", n = 137)
ICRS_per_gene_prop_hets <- variant_per_gene_count_ICRS_hets / ICR_trials_per_gene

##Counts of variants per gene for cases, BRIDGE and ICR (homs)
variant_per_gene_count_cases_homs <- data.frame()

for(i in 1:length(genes)){
  a <- sum(MPT_BRIDGE_ICR_counts[,11][which(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i])]) 
  variant_per_gene_count_cases_homs  <- append(variant_per_gene_count_cases_homs,a)
}

variant_per_gene_count_cases_homs <- t(data.frame(variant_per_gene_count_cases_homs))
gene_names <- as.vector(genes)
rownames(variant_per_gene_count_cases_homs) <- gene_names

variant_per_gene_count_BRIDGES_homs <- data.frame()

for(i in 1:length(genes)){
  a <- sum(MPT_BRIDGE_ICR_counts[,15][which(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i])]) 
  variant_per_gene_count_BRIDGES_homs  <- append(variant_per_gene_count_BRIDGES_homs,a)
}

variant_per_gene_count_BRIDGES_homs <- t(data.frame(variant_per_gene_count_BRIDGES_homs))
gene_names <- as.vector(genes)
rownames(variant_per_gene_count_BRIDGES_homs) <- gene_names

variant_per_gene_count_ICRS_homs <- data.frame()

for(i in 1:length(genes)){
  a <- sum(MPT_BRIDGE_ICR_counts[,19][which(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i])]) 
  variant_per_gene_count_ICRS_homs  <- append(variant_per_gene_count_ICRS_homs,a)
}

variant_per_gene_count_ICRS_homs <- t(data.frame(variant_per_gene_count_ICRS_homs))
gene_names <- as.vector(genes)
rownames(variant_per_gene_count_ICRS_homs) <- gene_names

##Fishers exact with fdr adjustment for counts of variants per gene (homs)
BRIDGE_pvalue_per_gene_homs <- as.vector(c(), mode = "any")
for(i in 1:nrow(variant_per_gene_count_cases_homs))
{
  fishtable = matrix(
    c(variant_per_gene_count_cases_homs[i,1],
      sum(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i]) * nrow(case_list) - variant_per_gene_count_cases_homs[i,1],
      variant_per_gene_count_BRIDGES_homs[i,1],
      sum(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i]) * nrow(BRIDGE_list) - variant_per_gene_count_BRIDGES_homs[i,1]),
    
    nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue_per_gene_homs <- append(BRIDGE_pvalue_per_gene_homs, fisherout[[1]])
}
BRIDGE_qvalue_per_gene_homs <- p.adjust(BRIDGE_pvalue_per_gene_homs, method = "fdr", n = 137)
cases_per_gene_prop_homs <- variant_per_gene_count_cases_homs / cases_trials_per_gene
BRIDGES_per_gene_prop_homs <- variant_per_gene_count_BRIDGES_homs / BRIDGE_trials_per_gene

ICR_pvalue_per_gene_homs <- as.vector(c(), mode = "any")
for(i in 1:nrow(variant_per_gene_count_cases_homs))
{
  fishtable = matrix(
    c(variant_per_gene_count_cases_homs[i,1],
      sum(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i]) * nrow(case_list) - variant_per_gene_count_cases_homs[i,1],
      variant_per_gene_count_BRIDGES_homs[i,1],
      sum(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i]) * nrow(ICR_list) - variant_per_gene_count_ICRS_homs[i,1]),
    
    nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  ICR_pvalue_per_gene_homs <- append(ICR_pvalue_per_gene_homs, fisherout[[1]])
}
ICR_qvalue_per_gene_homs <- p.adjust(ICR_pvalue_per_gene_homs, method = "fdr", n = 137)
ICRS_per_gene_prop_homs <- variant_per_gene_count_ICRS_homs / ICR_trials_per_gene

##Counts of variants per gene for cases, BRIDGE and ICR (hets_homs)
variant_per_gene_count_cases_hets_homs <- data.frame()

for(i in 1:length(genes)){
  a <- sum(MPT_BRIDGE_ICR_counts[,12][which(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i])]) 
  variant_per_gene_count_cases_hets_homs  <- append(variant_per_gene_count_cases_hets_homs,a)
}

variant_per_gene_count_cases_hets_homs <- t(data.frame(variant_per_gene_count_cases_hets_homs))
gene_names <- as.vector(genes)
rownames(variant_per_gene_count_cases_hets_homs) <- gene_names

variant_per_gene_count_BRIDGES_hets_homs <- data.frame()

for(i in 1:length(genes)){
  a <- sum(MPT_BRIDGE_ICR_counts[,16][which(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i])]) 
  variant_per_gene_count_BRIDGES_hets_homs  <- append(variant_per_gene_count_BRIDGES_hets_homs,a)
}

variant_per_gene_count_BRIDGES_hets_homs <- t(data.frame(variant_per_gene_count_BRIDGES_hets_homs))
gene_names <- as.vector(genes)
rownames(variant_per_gene_count_BRIDGES_hets_homs) <- gene_names

variant_per_gene_count_ICRS_hets_homs <- data.frame()

for(i in 1:length(genes)){
  a <- sum(MPT_BRIDGE_ICR_counts[,20][which(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i])]) 
  variant_per_gene_count_ICRS_hets_homs  <- append(variant_per_gene_count_ICRS_hets_homs,a)
}

variant_per_gene_count_ICRS_hets_homs <- t(data.frame(variant_per_gene_count_ICRS_hets_homs))
gene_names <- as.vector(genes)
rownames(variant_per_gene_count_ICRS_hets_homs) <- gene_names

##Fishers exact with fdr adjustment for counts of variants per gene (hets_homs)
BRIDGE_pvalue_per_gene_hets_homs <- as.vector(c(), mode = "any")

for(i in 1:nrow(variant_per_gene_count_cases_hets_homs))
{
  fishtable = matrix(
    c(variant_per_gene_count_cases_hets_homs[i,1],
      sum(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i]) * nrow(case_list) - variant_per_gene_count_cases_hets_homs[i,1],
      variant_per_gene_count_BRIDGES_hets_homs[i,1],
      sum(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i]) * nrow(BRIDGE_list) - variant_per_gene_count_BRIDGES_hets_homs[i,1]),
    
    nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue_per_gene_hets_homs <- append(BRIDGE_pvalue_per_gene_hets_homs, fisherout[[1]])
}
BRIDGE_qvalue_per_gene_hets_homs <- p.adjust(BRIDGE_pvalue_per_gene_hets_homs, method = "fdr", n = 137)
cases_per_gene_prop_hets_homs <- variant_per_gene_count_cases_hets_homs / cases_trials_per_gene
BRIDGES_per_gene_prop_hets_homs <- variant_per_gene_count_BRIDGES_hets_homs / BRIDGE_trials_per_gene

ICR_pvalue_per_gene_hets_homs <- as.vector(c(), mode = "any")

for(i in 1:nrow(variant_per_gene_count_cases_hets_homs))
{
  fishtable = matrix(
    c(variant_per_gene_count_cases_hets_homs[i,1],
      sum(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i]) * nrow(case_list) - variant_per_gene_count_cases_hets_homs[i,1],
      variant_per_gene_count_BRIDGES_hets_homs[i,1],
      sum(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i]) * nrow(ICR_list) - variant_per_gene_count_ICRS_hets_homs[i,1]),
    
    nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  ICR_pvalue_per_gene_hets_homs <- append(ICR_pvalue_per_gene_hets_homs, fisherout[[1]])
}
ICR_qvalue_per_gene_hets_homs <- p.adjust(ICR_pvalue_per_gene_hets_homs, method = "fdr", n = 137)
ICRS_per_gene_prop_hets_homs <- variant_per_gene_count_ICRS_hets_homs / ICR_trials_per_gene

##Counts of variants per gene for cases, BRIDGE and ICR (hets_homs_agg)
variant_per_gene_count_cases_hets_homs_agg <- data.frame()

for(i in 1:length(genes)){
  a <- sum(MPT_BRIDGE_ICR_counts[,13][which(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i])]) 
  variant_per_gene_count_cases_hets_homs_agg  <- append(variant_per_gene_count_cases_hets_homs_agg,a)
}

variant_per_gene_count_cases_hets_homs_agg <- t(data.frame(variant_per_gene_count_cases_hets_homs_agg))
gene_names <- as.vector(genes)
rownames(variant_per_gene_count_cases_hets_homs_agg) <- gene_names

##

variant_per_gene_count_BRIDGES_hets_homs_agg <- data.frame()

for(i in 1:length(genes)){
  a <- sum(MPT_BRIDGE_ICR_counts[,17][which(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i])]) 
  variant_per_gene_count_BRIDGES_hets_homs_agg  <- append(variant_per_gene_count_BRIDGES_hets_homs_agg,a)
}

variant_per_gene_count_BRIDGES_hets_homs_agg <- t(data.frame(variant_per_gene_count_BRIDGES_hets_homs_agg))
gene_names <- as.vector(genes)
rownames(variant_per_gene_count_BRIDGES_hets_homs_agg) <- gene_names

variant_per_gene_count_ICRS_hets_homs_agg <- data.frame()

for(i in 1:length(genes)){
  a <- sum(MPT_BRIDGE_ICR_counts[,21][which(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i])]) 
  variant_per_gene_count_ICRS_hets_homs_agg  <- append(variant_per_gene_count_ICRS_hets_homs_agg,a)
}

variant_per_gene_count_ICRS_hets_homs_agg <- t(data.frame(variant_per_gene_count_ICRS_hets_homs_agg))
gene_names <- as.vector(genes)
rownames(variant_per_gene_count_ICRS_hets_homs_agg) <- gene_names

##Fishers exact with fdr adjustment for counts of variants per gene (hets_homs_agg)
BRIDGE_pvalue_per_gene_hets_homs_agg <- as.vector(c(), mode = "any")

for(i in 1:nrow(variant_per_gene_count_cases_hets_homs_agg))
{
  fishtable = matrix(
    c(variant_per_gene_count_cases_hets_homs_agg[i,1],
      sum(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i]) * (nrow(case_list)*2) - variant_per_gene_count_cases_hets_homs_agg[i,1],
      variant_per_gene_count_BRIDGES_hets_homs_agg[i,1],
      sum(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i]) * (nrow(BRIDGE_list)*2) - variant_per_gene_count_BRIDGES_hets_homs_agg[i,1]),
    
    nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue_per_gene_hets_homs_agg <- append(BRIDGE_pvalue_per_gene_hets_homs_agg, fisherout[[1]])
}
BRIDGE_qvalue_per_gene_hets_homs_agg <- p.adjust(BRIDGE_pvalue_per_gene_hets_homs_agg, method = "fdr", n = 137)
cases_per_gene_prop_hets_homs_agg <- variant_per_gene_count_cases_hets_homs_agg / (cases_trials_per_gene *2)
BRIDGES_per_gene_prop_hets_homs_agg <- variant_per_gene_count_BRIDGES_hets_homs_agg / (BRIDGE_trials_per_gene *2)

ICR_pvalue_per_gene_hets_homs_agg <- as.vector(c(), mode = "any")

for(i in 1:nrow(variant_per_gene_count_cases_hets_homs_agg))
{
  fishtable = matrix(
    c(variant_per_gene_count_cases_hets_homs_agg[i,1],
      sum(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i]) * nrow(case_list) - variant_per_gene_count_cases_hets_homs_agg[i,1],
      variant_per_gene_count_BRIDGES_hets_homs_agg[i,1],
      sum(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i]) * nrow(ICR_list) - variant_per_gene_count_ICRS_hets_homs_agg[i,1]),
    
    nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  ICR_pvalue_per_gene_hets_homs_agg <- append(ICR_pvalue_per_gene_hets_homs_agg, fisherout[[1]])
}
ICR_qvalue_per_gene_hets_homs_agg <- p.adjust(ICR_pvalue_per_gene_hets_homs_agg, method = "fdr", n = 137)
ICRS_per_gene_prop_hets_homs_agg <- variant_per_gene_count_ICRS_hets_homs_agg / (ICR_trials_per_gene *2)

##Counts of individuals with variants per gene for cases and controls (hets)
trimmed_cases_hets <- data.frame(MPT_BRIDGE_ICR_counts$cases_hets)

indv_with_variant_per_gene_cases_hets <- data.frame()

for(i in 1:length(genes)){
  a <- trimmed_cases_hets[which(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i]),] #gets the case IDs from any line where the gene name occurs. CHANGE TO PICK FROM TRIMMED VARIANT TABLE
  b <- a[which(a != "")]
  c <- paste(b, collapse = ";")
  indv_with_variant_per_gene_cases_hets <- append(indv_with_variant_per_gene_cases_hets, c)
}

indv_with_variant_per_gene_cases_hets <- unlist(indv_with_variant_per_gene_cases_hets)
indv_with_variant_per_gene_cases_hets <- strsplit(indv_with_variant_per_gene_cases_hets, ";", fixed = T)
unique_indv_with_variant_per_gene_cases_hets <- lapply(indv_with_variant_per_gene_cases_hets, unique)

indv_per_gene_count_cases_hets <- data.frame()

for(i in 1:length(unique_indv_with_variant_per_gene_cases_hets)){
  a <- length(which(unique_indv_with_variant_per_gene_cases_hets[[i]] != ""))
  indv_per_gene_count_cases_hets <- append(indv_per_gene_count_cases_hets, a)
}

indv_per_gene_count_cases_hets <- unlist(indv_per_gene_count_cases_hets)

trimmed_BRIDGES_hets <- data.frame(MPT_BRIDGE_ICR_counts$BRIDGES_hets)

indv_with_variant_per_gene_BRIDGES_hets <- data.frame()

for(i in 1:length(genes)){
  a <- trimmed_BRIDGES_hets[which(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i]),]
  b <- a[which(a != "")]
  c <- paste(b, collapse = ";")
  indv_with_variant_per_gene_BRIDGES_hets <- append(indv_with_variant_per_gene_BRIDGES_hets, c)
}

indv_with_variant_per_gene_BRIDGES_hets <- unlist(indv_with_variant_per_gene_BRIDGES_hets)
indv_with_variant_per_gene_BRIDGES_hets <- strsplit(indv_with_variant_per_gene_BRIDGES_hets, ";", fixed = T)
unique_indv_with_variant_per_gene_BRIDGES_hets <- lapply(indv_with_variant_per_gene_BRIDGES_hets, unique)


indv_per_gene_count_BRIDGES_hets <- data.frame()

for(i in 1:length(unique_indv_with_variant_per_gene_BRIDGES_hets)){
  a <- length(which(unique_indv_with_variant_per_gene_BRIDGES_hets[[i]] != ""))
  indv_per_gene_count_BRIDGES_hets <- append(indv_per_gene_count_BRIDGES_hets, a)
}

indv_per_gene_count_BRIDGES_hets <- unlist(indv_per_gene_count_BRIDGES_hets)

trimmed_ICRS_hets <- data.frame(MPT_BRIDGE_ICR_counts$ICRS_hets)

indv_with_variant_per_gene_ICRS_hets <- data.frame()

for(i in 1:length(genes)){
  a <- trimmed_ICRS_hets[which(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i]),]
  b <- a[which(a != "")]
  c <- paste(b, collapse = ";")
  indv_with_variant_per_gene_ICRS_hets <- append(indv_with_variant_per_gene_ICRS_hets, c)
}

indv_with_variant_per_gene_ICRS_hets <- unlist(indv_with_variant_per_gene_ICRS_hets)
indv_with_variant_per_gene_ICRS_hets <- strsplit(indv_with_variant_per_gene_ICRS_hets, ";", fixed = T)
unique_indv_with_variant_per_gene_ICRS_hets <- lapply(indv_with_variant_per_gene_ICRS_hets, unique)


indv_per_gene_count_ICRS_hets <- data.frame()

for(i in 1:length(unique_indv_with_variant_per_gene_ICRS_hets)){
  a <- length(which(unique_indv_with_variant_per_gene_ICRS_hets[[i]] != ""))
  indv_per_gene_count_ICRS_hets <- append(indv_per_gene_count_ICRS_hets, a)
}

indv_per_gene_count_ICRS_hets <- unlist(indv_per_gene_count_ICRS_hets)

##Fishers exact with fdr adjustment for counts of individuals with variants per gene (hets)
BRIDGE_pvalue_per_indv_per_gene_hets <- as.vector(c(), mode = "any")
cases_indv_per_gene_hets <- as.vector(c(), mode = "any")
BRIDGES_indv_per_gene_hets <- as.vector(c(), mode = "any")

for(i in 1:length(indv_per_gene_count_cases_hets))
{
  fishtable = matrix(
    c(indv_per_gene_count_cases_hets[i],
      nrow(case_list) - indv_per_gene_count_cases_hets[i],
      indv_per_gene_count_BRIDGES_hets[i],
      nrow(BRIDGE_list) - indv_per_gene_count_BRIDGES_hets[i]), 
    nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue_per_indv_per_gene_hets <- append(BRIDGE_pvalue_per_indv_per_gene_hets, fisherout[[1]])
}
BRIDGE_qvalue_per_indv_per_gene_hets <- p.adjust(BRIDGE_pvalue_per_indv_per_gene_hets, method = "fdr", n = 137)
cases_per_indv_per_gene_prop_hets <- indv_per_gene_count_cases_hets / nrow(case_list)
BRIDGES_per_indv_per_gene_prop_hets <- indv_per_gene_count_BRIDGES_hets / nrow(BRIDGE_list)

ICR_pvalue_per_indv_per_gene_hets <- as.vector(c(), mode = "any")
cases_indv_per_gene_hets <- as.vector(c(), mode = "any")
ICRS_indv_per_gene_hets <- as.vector(c(), mode = "any")

for(i in 1:length(indv_per_gene_count_cases_hets))
{
  fishtable = matrix(
    c(indv_per_gene_count_cases_hets[i],
      nrow(case_list) - indv_per_gene_count_cases_hets[i],
      indv_per_gene_count_ICRS_hets[i],
      nrow(BRIDGE_list) - indv_per_gene_count_ICRS_hets[i]), 
    nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  ICR_pvalue_per_indv_per_gene_hets <- append(ICR_pvalue_per_indv_per_gene_hets, fisherout[[1]])
}
ICR_qvalue_per_indv_per_gene_hets <- p.adjust(ICR_pvalue_per_indv_per_gene_hets, method = "fdr", n = 137)
ICRS_per_indv_per_gene_prop_hets <- indv_per_gene_count_ICRS_hets / nrow(ICR_list)

##Counts of individuals with variants per gene for cases and controls (homs)
trimmed_cases_homs <- data.frame(MPT_BRIDGE_ICR_counts$cases_homs)

indv_with_variant_per_gene_cases_homs <- data.frame()

for(i in 1:length(genes)){
  a <- trimmed_cases_homs[which(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i]),] #gets the case IDs from any line where the gene name occurs. CHANGE TO PICK FROM TRIMMED VARIANT TABLE
  b <- a[which(a != "")]
  c <- paste(b, collapse = ";")
  indv_with_variant_per_gene_cases_homs <- append(indv_with_variant_per_gene_cases_homs, c)
}

indv_with_variant_per_gene_cases_homs <- unlist(indv_with_variant_per_gene_cases_homs)
indv_with_variant_per_gene_cases_homs <- strsplit(indv_with_variant_per_gene_cases_homs, ";", fixed = T)
unique_indv_with_variant_per_gene_cases_homs <- lapply(indv_with_variant_per_gene_cases_homs, unique)

indv_per_gene_count_cases_homs <- data.frame()

for(i in 1:length(unique_indv_with_variant_per_gene_cases_homs)){
  a <- length(which(unique_indv_with_variant_per_gene_cases_homs[[i]] != ""))
  indv_per_gene_count_cases_homs <- append(indv_per_gene_count_cases_homs, a)
}

indv_per_gene_count_cases_homs <- unlist(indv_per_gene_count_cases_homs)

trimmed_BRIDGES_homs <- data.frame(MPT_BRIDGE_ICR_counts$BRIDGES_homs)

indv_with_variant_per_gene_BRIDGES_homs <- data.frame()

for(i in 1:length(genes)){
  a <- trimmed_BRIDGES_homs[which(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i]),]
  b <- a[which(a != "")]
  c <- paste(b, collapse = ";")
  indv_with_variant_per_gene_BRIDGES_homs <- append(indv_with_variant_per_gene_BRIDGES_homs, c)
}

indv_with_variant_per_gene_BRIDGES_homs <- unlist(indv_with_variant_per_gene_BRIDGES_homs)
indv_with_variant_per_gene_BRIDGES_homs <- strsplit(indv_with_variant_per_gene_BRIDGES_homs, ";", fixed = T)
unique_indv_with_variant_per_gene_BRIDGES_homs <- lapply(indv_with_variant_per_gene_BRIDGES_homs, unique)


indv_per_gene_count_BRIDGES_homs <- data.frame()

for(i in 1:length(unique_indv_with_variant_per_gene_BRIDGES_homs)){
  a <- length(which(unique_indv_with_variant_per_gene_BRIDGES_homs[[i]] != ""))
  indv_per_gene_count_BRIDGES_homs <- append(indv_per_gene_count_BRIDGES_homs, a)
}

indv_per_gene_count_BRIDGES_homs <- unlist(indv_per_gene_count_BRIDGES_homs)

trimmed_ICRS_homs <- data.frame(MPT_BRIDGE_ICR_counts$ICRS_homs)

indv_with_variant_per_gene_ICRS_homs <- data.frame()

for(i in 1:length(genes)){
  a <- trimmed_ICRS_homs[which(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i]),]
  b <- a[which(a != "")]
  c <- paste(b, collapse = ";")
  indv_with_variant_per_gene_ICRS_homs <- append(indv_with_variant_per_gene_ICRS_homs, c)
}

indv_with_variant_per_gene_ICRS_homs <- unlist(indv_with_variant_per_gene_ICRS_homs)
indv_with_variant_per_gene_ICRS_homs <- strsplit(indv_with_variant_per_gene_ICRS_homs, ";", fixed = T)
unique_indv_with_variant_per_gene_ICRS_homs <- lapply(indv_with_variant_per_gene_ICRS_homs, unique)


indv_per_gene_count_ICRS_homs <- data.frame()

for(i in 1:length(unique_indv_with_variant_per_gene_ICRS_homs)){
  a <- length(which(unique_indv_with_variant_per_gene_ICRS_homs[[i]] != ""))
  indv_per_gene_count_ICRS_homs <- append(indv_per_gene_count_ICRS_homs, a)
}

indv_per_gene_count_ICRS_homs <- unlist(indv_per_gene_count_ICRS_homs)

##Fishers exact with fdr adjustment for counts of individuals with variants per gene (homs)
BRIDGE_pvalue_per_indv_per_gene_homs <- as.vector(c(), mode = "any")
cases_indv_per_gene_homs <- as.vector(c(), mode = "any")
BRIDGES_indv_per_gene_homs <- as.vector(c(), mode = "any")

for(i in 1:length(indv_per_gene_count_cases_homs))
{
  fishtable = matrix(
    c(indv_per_gene_count_cases_homs[i],
      nrow(case_list) - indv_per_gene_count_cases_homs[i],
      indv_per_gene_count_BRIDGES_homs[i],
      nrow(BRIDGE_list) - indv_per_gene_count_BRIDGES_homs[i]), 
    nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue_per_indv_per_gene_homs <- append(BRIDGE_pvalue_per_indv_per_gene_homs, fisherout[[1]])
}
BRIDGE_qvalue_per_indv_per_gene_homs <- p.adjust(BRIDGE_pvalue_per_indv_per_gene_homs, method = "fdr", n = 137)
cases_per_indv_per_gene_prop_homs <- indv_per_gene_count_cases_homs / nrow(case_list)
BRIDGES_per_indv_per_gene_prop_homs <- indv_per_gene_count_BRIDGES_homs / nrow(BRIDGE_list)

ICR_pvalue_per_indv_per_gene_homs <- as.vector(c(), mode = "any")
cases_indv_per_gene_homs <- as.vector(c(), mode = "any")
ICRS_indv_per_gene_homs <- as.vector(c(), mode = "any")

for(i in 1:length(indv_per_gene_count_cases_homs))
{
  fishtable = matrix(
    c(indv_per_gene_count_cases_homs[i],
      nrow(case_list) - indv_per_gene_count_cases_homs[i],
      indv_per_gene_count_ICRS_homs[i],
      nrow(BRIDGE_list) - indv_per_gene_count_ICRS_homs[i]), 
    nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  ICR_pvalue_per_indv_per_gene_homs <- append(ICR_pvalue_per_indv_per_gene_homs, fisherout[[1]])
}
ICR_qvalue_per_indv_per_gene_homs <- p.adjust(ICR_pvalue_per_indv_per_gene_homs, method = "fdr", n = 137)
ICRS_per_indv_per_gene_prop_homs <- indv_per_gene_count_ICRS_homs / nrow(ICR_list)

##Counts of individuals with variants per gene for cases (hets_homs)
trimmed_cases_hets_homs <- paste(trimmed_cases_hets$MPT_BRIDGE_ICR_counts.cases_hets, trimmed_cases_homs$MPT_BRIDGE_ICR_counts.cases_homs, sep = ";")
trimmed_cases_hets_homs <- data.frame(trimmed_cases_hets_homs)

indv_with_variant_per_gene_cases_hets_homs <- data.frame()

for(i in 1:length(genes)){
  a <- trimmed_cases_hets_homs[which(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i]),] 
  b <- a[which(a != ";")]
  c <- paste(b, collapse = ";")
  indv_with_variant_per_gene_cases_hets_homs <- append(indv_with_variant_per_gene_cases_hets_homs, c)
}

indv_with_variant_per_gene_cases_hets_homs <- unlist(indv_with_variant_per_gene_cases_hets_homs)
indv_with_variant_per_gene_cases_hets_homs <- strsplit(indv_with_variant_per_gene_cases_hets_homs, ";", fixed = T)
unique_indv_with_variant_per_gene_cases_hets_homs <- lapply(indv_with_variant_per_gene_cases_hets_homs, unique)

indv_per_gene_count_cases_hets_homs <- data.frame()

for(i in 1:length(unique_indv_with_variant_per_gene_cases_hets_homs)){
  a <- length(which(unique_indv_with_variant_per_gene_cases_hets_homs[[i]] != ""))
  indv_per_gene_count_cases_hets_homs <- append(indv_per_gene_count_cases_hets_homs, a)
}

indv_per_gene_count_cases_hets_homs <- unlist(indv_per_gene_count_cases_hets_homs)

##Counts of individuals with variants per gene for BRIDGE controls (hets_homs)
trimmed_BRIDGES_hets_homs <- paste(trimmed_BRIDGES_hets$MPT_BRIDGE_ICR_counts.BRIDGES_hets, trimmed_BRIDGES_homs$MPT_BRIDGE_ICR_counts.BRIDGES_homs, sep = ";")
trimmed_BRIDGES_hets_homs <- data.frame(trimmed_BRIDGES_hets_homs)

indv_with_variant_per_gene_BRIDGES_hets_homs <- data.frame()

for(i in 1:length(genes)){
  a <- trimmed_BRIDGES_hets_homs[which(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i]),]
  b <- a[which(a != ";")]
  c <- paste(b, collapse = ";")
  indv_with_variant_per_gene_BRIDGES_hets_homs <- append(indv_with_variant_per_gene_BRIDGES_hets_homs, c)
}

indv_with_variant_per_gene_BRIDGES_hets_homs <- unlist(indv_with_variant_per_gene_BRIDGES_hets_homs)
indv_with_variant_per_gene_BRIDGES_hets_homs <- strsplit(indv_with_variant_per_gene_BRIDGES_hets_homs, ";", fixed = T)
unique_indv_with_variant_per_gene_BRIDGES_hets_homs <- lapply(indv_with_variant_per_gene_BRIDGES_hets_homs, unique)

indv_per_gene_count_BRIDGES_hets_homs <- data.frame()

for(i in 1:length(unique_indv_with_variant_per_gene_BRIDGES_hets_homs)){
  a <- length(which(unique_indv_with_variant_per_gene_BRIDGES_hets_homs[[i]] != ""))
  indv_per_gene_count_BRIDGES_hets_homs <- append(indv_per_gene_count_BRIDGES_hets_homs, a)
}

indv_per_gene_count_BRIDGES_hets_homs <- unlist(indv_per_gene_count_BRIDGES_hets_homs)

##Counts of individuals with variants per gene for ICR controls (hets_homs)
trimmed_ICRS_hets_homs <- paste(trimmed_ICRS_hets$MPT_BRIDGE_ICR_counts.ICRS_hets, trimmed_ICRS_homs$MPT_BRIDGE_ICR_counts.ICRS_homs, sep = ";")
trimmed_ICRS_hets_homs <- data.frame(trimmed_ICRS_hets_homs)

indv_with_variant_per_gene_ICRS_hets_homs <- data.frame()

for(i in 1:length(genes)){
  a <- trimmed_ICRS_hets_homs[which(MPT_BRIDGE_ICR_counts$gene_col %in% genes[i]),]
  b <- a[which(a != ";")]
  c <- paste(b, collapse = ";")
  indv_with_variant_per_gene_ICRS_hets_homs <- append(indv_with_variant_per_gene_ICRS_hets_homs, c)
}

indv_with_variant_per_gene_ICRS_hets_homs <- unlist(indv_with_variant_per_gene_ICRS_hets_homs)
indv_with_variant_per_gene_ICRS_hets_homs <- strsplit(indv_with_variant_per_gene_ICRS_hets_homs, ";", fixed = T)
unique_indv_with_variant_per_gene_ICRS_hets_homs <- lapply(indv_with_variant_per_gene_ICRS_hets_homs, unique)

indv_per_gene_count_ICRS_hets_homs <- data.frame()

for(i in 1:length(unique_indv_with_variant_per_gene_ICRS_hets_homs)){
  a <- length(which(unique_indv_with_variant_per_gene_ICRS_hets_homs[[i]] != ""))
  indv_per_gene_count_ICRS_hets_homs <- append(indv_per_gene_count_ICRS_hets_homs, a)
}

indv_per_gene_count_ICRS_hets_homs <- unlist(indv_per_gene_count_ICRS_hets_homs)

##Fishers exact with fdr adjustment for counts of individuals with variants per gene (hets_homs)
BRIDGE_pvalue_per_indv_per_gene_hets_homs <- as.vector(c(), mode = "any")
cases_indv_per_gene_hets_homs <- as.vector(c(), mode = "any")
BRIDGES_indv_per_gene_hets_homs <- as.vector(c(), mode = "any")

for(i in 1:length(indv_per_gene_count_cases_hets_homs))
{
  fishtable = matrix(
    c(indv_per_gene_count_cases_hets_homs[i],
      nrow(case_list) - indv_per_gene_count_cases_hets_homs[i],
      indv_per_gene_count_BRIDGES_hets_homs[i],
      nrow(BRIDGE_list) - indv_per_gene_count_BRIDGES_hets_homs[i]), 
    nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue_per_indv_per_gene_hets_homs <- append(BRIDGE_pvalue_per_indv_per_gene_hets_homs, fisherout[[1]])
}
BRIDGE_qvalue_per_indv_per_gene_hets_homs <- p.adjust(BRIDGE_pvalue_per_indv_per_gene_hets_homs, method = "fdr", n = 137)
cases_per_indv_per_gene_prop_hets_homs <- indv_per_gene_count_cases_hets_homs / nrow(case_list)
BRIDGES_per_indv_per_gene_prop_hets_homs <- indv_per_gene_count_BRIDGES_homs / nrow(BRIDGE_list)

ICR_pvalue_per_indv_per_gene_hets_homs <- as.vector(c(), mode = "any")
cases_indv_per_gene_hets_homs <- as.vector(c(), mode = "any")
ICRS_indv_per_gene_hets_homs <- as.vector(c(), mode = "any")

for(i in 1:length(indv_per_gene_count_cases_hets_homs))
{
  fishtable = matrix(
    c(indv_per_gene_count_cases_hets_homs[i],
      nrow(case_list) - indv_per_gene_count_cases_hets_homs[i],
      indv_per_gene_count_ICRS_hets_homs[i],
      nrow(BRIDGE_list) - indv_per_gene_count_ICRS_hets_homs[i]), 
    nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  ICR_pvalue_per_indv_per_gene_hets_homs <- append(ICR_pvalue_per_indv_per_gene_hets_homs, fisherout[[1]])
}
ICR_qvalue_per_indv_per_gene_hets_homs <- p.adjust(ICR_pvalue_per_indv_per_gene_hets_homs, method = "fdr", n = 137) 
ICRS_per_indv_per_gene_prop_hets_homs <- indv_per_gene_count_ICRS_hets_homs / nrow(ICR_list)

##Collation of per gene results
variants_per_gene_table <- as.data.frame(cbind(rownames(variant_per_gene_count_cases_hets),
                                               
                                               variant_per_gene_count_cases_hets,
                                               variant_per_gene_count_BRIDGES_hets,
                                               variant_per_gene_count_ICRS_hets,
                                               BRIDGE_pvalue_per_gene_hets,
                                               BRIDGE_qvalue_per_gene_hets,
                                               ICR_pvalue_per_gene_hets,
                                               ICR_qvalue_per_gene_hets,
                                               cases_per_gene_prop_hets,
                                               BRIDGES_per_gene_prop_hets,
                                               ICRS_per_gene_prop_hets,
                                               
                                               variant_per_gene_count_cases_homs,
                                               variant_per_gene_count_BRIDGES_homs,
                                               variant_per_gene_count_ICRS_homs,
                                               BRIDGE_pvalue_per_gene_homs,
                                               BRIDGE_qvalue_per_gene_homs,
                                               ICR_pvalue_per_gene_homs,
                                               ICR_qvalue_per_gene_homs,
                                               cases_per_gene_prop_homs,
                                               BRIDGES_per_gene_prop_homs,
                                               ICRS_per_gene_prop_homs,
                                               
                                               variant_per_gene_count_cases_hets_homs,
                                               variant_per_gene_count_BRIDGES_hets_homs,
                                               variant_per_gene_count_ICRS_hets_homs,
                                               BRIDGE_pvalue_per_gene_hets_homs,
                                               BRIDGE_qvalue_per_gene_hets_homs,
                                               ICR_pvalue_per_gene_hets_homs,
                                               ICR_qvalue_per_gene_hets_homs,
                                               cases_per_gene_prop_hets_homs,
                                               BRIDGES_per_gene_prop_hets_homs,
                                               ICRS_per_gene_prop_hets_homs,
                                               
                                               variant_per_gene_count_cases_hets_homs_agg,
                                               variant_per_gene_count_BRIDGES_hets_homs_agg,
                                               variant_per_gene_count_ICRS_hets_homs_agg,
                                               BRIDGE_pvalue_per_gene_hets_homs_agg,
                                               BRIDGE_qvalue_per_gene_hets_homs_agg,
                                               ICR_pvalue_per_gene_hets_homs_agg,
                                               ICR_qvalue_per_gene_hets_homs_agg,
                                               cases_per_gene_prop_hets_homs_agg,
                                               BRIDGES_per_gene_prop_hets_homs_agg,
                                               ICRS_per_gene_prop_hets_homs_agg,
                                               
                                               indv_per_gene_count_cases_hets,
                                               indv_per_gene_count_BRIDGES_hets,
                                               indv_per_gene_count_ICRS_hets,
                                               BRIDGE_pvalue_per_indv_per_gene_hets,
                                               BRIDGE_qvalue_per_indv_per_gene_hets,
                                               cases_per_indv_per_gene_prop_hets,
                                               BRIDGES_per_indv_per_gene_prop_hets,
                                               ICR_pvalue_per_indv_per_gene_hets,
                                               ICR_qvalue_per_indv_per_gene_hets,
                                               ICRS_per_indv_per_gene_prop_hets,
                                               
                                               indv_per_gene_count_cases_homs,
                                               indv_per_gene_count_BRIDGES_homs,
                                               indv_per_gene_count_ICRS_homs,
                                               BRIDGE_pvalue_per_indv_per_gene_homs,
                                               BRIDGE_qvalue_per_indv_per_gene_homs,
                                               cases_per_indv_per_gene_prop_homs,
                                               BRIDGES_per_indv_per_gene_prop_homs,
                                               ICR_pvalue_per_indv_per_gene_homs,
                                               ICR_qvalue_per_indv_per_gene_homs,
                                               ICRS_per_indv_per_gene_prop_homs,
                                               
                                               indv_per_gene_count_cases_hets_homs,
                                               indv_per_gene_count_BRIDGES_hets_homs,
                                               indv_per_gene_count_ICRS_hets_homs,
                                               BRIDGE_pvalue_per_indv_per_gene_hets_homs,
                                               BRIDGE_qvalue_per_indv_per_gene_hets_homs,
                                               cases_per_indv_per_gene_prop_hets_homs,
                                               BRIDGES_per_indv_per_gene_prop_hets_homs,
                                               ICR_pvalue_per_indv_per_gene_hets_homs,
                                               ICR_qvalue_per_indv_per_gene_hets_homs,
                                               ICRS_per_indv_per_gene_prop_hets_homs))

colnames(variants_per_gene_table) <- c("gene", 
                                       
                                       "variant_per_gene_count_cases_hets",
                                       "variant_per_gene_count_BRIDGES_hets",
                                       "variant_per_gene_count_ICRS_hets",
                                       "BRIDGE_pvalue_per_gene_hets",
                                       "BRIDGE_qvalue_per_gene_hets",
                                       "ICR_pvalue_per_gene_hets",
                                       "ICR_qvalue_per_gene_hets",
                                       "cases_per_gene_prop_hets",
                                       "BRIDGES_per_gene_prop_hets",
                                       "ICRS_per_gene_prop_hets",
                                       
                                       "variant_per_gene_count_cases_homs",
                                       "variant_per_gene_count_BRIDGES_homs",
                                       "variant_per_gene_count_ICRS_homs",
                                       "BRIDGE_pvalue_per_gene_homs",
                                       "BRIDGE_qvalue_per_gene_homs",
                                       "ICR_pvalue_per_gene_homs",
                                       "ICR_qvalue_per_gene_homs",
                                       "cases_per_gene_prop_homs",
                                       "BRIDGES_per_gene_prop_homs",
                                       "ICRS_per_gene_prop_homs",
                                       
                                       "variant_per_gene_count_cases_hets_homs",
                                       "variant_per_gene_count_BRIDGES_hets_homs",
                                       "variant_per_gene_count_ICRS_hets_homs",
                                       "BRIDGE_pvalue_per_gene_hets_homs",
                                       "BRIDGE_qvalue_per_gene_hets_homs",
                                       "ICR_pvalue_per_gene_hets_homs",
                                       "ICR_qvalue_per_gene_hets_homs",
                                       "cases_per_gene_prop_hets_homs",
                                       "BRIDGES_per_gene_prop_hets_homs",
                                       "ICRS_per_gene_prop_hets_homs",
                                       
                                       "variant_per_gene_count_cases_hets_homs_agg",
                                       "variant_per_gene_count_BRIDGES_hets_homs_agg",
                                       "variant_per_gene_count_ICRS_hets_homs_agg",
                                       "BRIDGE_pvalue_per_gene_hets_homs_agg",
                                       "BRIDGE_qvalue_per_gene_hets_homs_agg",
                                       "ICR_pvalue_per_gene_hets_homs_agg",
                                       "ICR_qvalue_per_gene_hets_homs_agg",
                                       "cases_per_gene_prop_hets_homs_agg",
                                       "BRIDGES_per_gene_prop_hets_homs_agg",
                                       "ICRS_per_gene_prop_hets_homs_agg",
                                       
                                       "indv_per_gene_count_cases_hets",
                                       "indv_per_gene_count_BRIDGES_hets",
                                       "indv_per_gene_count_ICRS_hets",
                                       "BRIDGE_pvalue_per_indv_per_gene_hets",
                                       "BRIDGE_qvalue_per_indv_per_gene_hets",
                                       "cases_per_indv_per_gene_prop_hets",
                                       "BRIDGES_per_indv_per_gene_prop_hets",
                                       "ICR_pvalue_per_indv_per_gene_hets",
                                       "ICR_qvalue_per_indv_per_gene_hets",
                                       "ICRS_per_indv_per_gene_prop_hets",
                                       
                                       "indv_per_gene_count_cases_homs",
                                       "indv_per_gene_count_BRIDGES_homs",
                                       "indv_per_gene_count_ICRS_homs",
                                       "BRIDGE_pvalue_per_indv_per_gene_homs",
                                       "BRIDGE_qvalue_per_indv_per_gene_homs",
                                       "cases_per_indv_per_gene_prop_homs",
                                       "BRIDGES_per_indv_per_gene_prop_homs",
                                       "ICR_pvalue_per_indv_per_gene_homs",
                                       "ICR_qvalue_per_indv_per_gene_homs",
                                       "ICRS_per_indv_per_gene_prop_homs",
                                       
                                       "indv_per_gene_count_cases_hets_homs",
                                       "indv_per_gene_count_BRIDGES_hets_homs",
                                       "indv_per_gene_count_ICRS_hets_homs",
                                       "BRIDGE_pvalue_per_indv_per_gene_hets_homs",
                                       "BRIDGE_qvalue_per_indv_per_gene_hets_homs",
                                       "cases_per_indv_per_gene_prop_hets_homs",
                                       "BRIDGES_per_indv_per_gene_prop_hets_homs",
                                       "ICR_pvalue_per_indv_per_gene_hets_homs",
                                       "ICR_qvalue_per_indv_per_gene_hets_homs",
                                       "ICRS_per_indv_per_gene_prop_hets_homs")

variants_per_gene_table$BRIDGE_qvalue_per_gene_hets <- as.numeric(as.character(variants_per_gene_table$BRIDGE_qvalue_per_gene_hets))
variants_per_gene_table$BRIDGE_qvalue_per_gene_homs <- as.numeric(as.character(variants_per_gene_table$BRIDGE_qvalue_per_gene_homs))
variants_per_gene_table$BRIDGE_qvalue_per_gene_hets_homs <- as.numeric(as.character(variants_per_gene_table$BRIDGE_qvalue_per_gene_hets_homs))
variants_per_gene_table$BRIDGE_qvalue_per_gene_hets_homs_agg <- as.numeric(as.character(variants_per_gene_table$BRIDGE_qvalue_per_gene_hets_homs_agg))

variants_per_gene_table$BRIDGE_qvalue_per_indv_per_gene_hets <- as.numeric(as.character(variants_per_gene_table$BRIDGE_qvalue_per_indv_per_gene_hets))
variants_per_gene_table$BRIDGE_qvalue_per_indv_per_gene_homs <- as.numeric(as.character(variants_per_gene_table$BRIDGE_qvalue_per_indv_per_gene_homs))
variants_per_gene_table$BRIDGE_qvalue_per_indv_per_gene_hets_homs <- as.numeric(as.character(variants_per_gene_table$BRIDGE_qvalue_per_indv_per_gene_hets_homs))

sig_index <- which(
(
  
((BRIDGE_qvalue_per_gene_hets < 0.05) && (cases_per_gene_prop_hets > BRIDGES_per_gene_prop_hets))
|
((BRIDGE_qvalue_per_gene_homs < 0.05) && (cases_per_gene_prop_homs > BRIDGES_per_gene_prop_homs))
|
((BRIDGE_qvalue_per_gene_hets_homs < 0.05) && (cases_per_gene_prop_hets_homs > BRIDGES_per_gene_prop_hets_homs))
|
((BRIDGE_qvalue_per_gene_hets_homs_agg < 0.05) && (cases_per_gene_prop_hets_homs_agg > BRIDGES_per_gene_prop_hets_homs_agg))
|

((BRIDGE_qvalue_per_indv_per_gene_hets < 0.05) && (cases_per_indv_per_gene_prop_hets > BRIDGES_per_indv_per_gene_prop_hets))
|
((BRIDGE_qvalue_per_indv_per_gene_homs < 0.05) && (cases_per_indv_per_gene_prop_homs > BRIDGES_per_indv_per_gene_prop_homs))
|
((BRIDGE_qvalue_per_indv_per_gene_hets_homs < 0.05) && (cases_per_indv_per_gene_prop_hets_homs > BRIDGES_per_indv_per_gene_prop_hets_homs))
)
)
sig_index <- unique(sig_index)

variants_per_gene_table_sig <- variants_per_gene_table[sig_index,]


##Output results
write.csv(variants_per_gene_table_sig, paste(args[1],"_variants_per_gene_table_sig.csv", sep = "")) 
write.csv(variants_per_gene_table, paste(args[1],"_variants_per_gene_table_all.csv", sep = "")) 

###############################################################################################################################
##Scripts to count single nucleotide variant and indel calls affecting elements of interest and compare frequency vs controls##
###############################################################################################################################

########################################################################################################################################################################################################################
##Scripts contain functions that were not utilised in interpretation of results as follows:
#1.Counting of variants per gene/element of interest where the numerator was taken as the total number of variants observed in that gene in cases /controls and the denominator was taken as the number of cases/controls multiplied the number of variant sites observed
#2.Counts of individual variants or variants per gene/element where heterozygous and homozygous variants where aggregated (homozygous counting as 2 and heterozgous as 1) to consider the total frequency of variant alleles (denoted "hets_homs_agg") in cases or controls rather than the frequency of individual variants or individuals with variants in a given gene/element
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

##Filter per chromosome VCF files using BED file corresponding to coordinates of interest
for i in `cat chromosome_files.txt`; do
bcftools view -R gtex_positions.txt -e 'FILTER!="PASS"' -o ${i}_filtered.vcf -O v /scratch/WGS10K/data/release/20170614-A/merged-vcf/no_hgmd/${i}.vcf.gz
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
tabix merged.vcf.gz

##Filter merged VCF based on quality parameters (insert missing genotypes where criteria not fulfilled)
bcftools filter -e 'FMT/DP<10 || FMT/GQ<30 || (FMT/AD[1]/(FMT/AD[0] + FMT/AD[1]) <0.3) || (FMT/AD[1]/(FMT/AD[0] + FMT/AD[1] + FMT/AD[2]) <0.3) || (FMT/AD[2]/(FMT/AD[0] + FMT/AD[1] + FMT/AD[2]) <0.3)' -S . -o gtex_merged_filtered.vcf.gz -O z merged.vcf.gz

##Transferred to local server

######################
##Further  filtering##
######################

##Execute VEP filter script 
/home/jww39/ensembl-vep/filter_vep -i VEP_out_gtex_merged_filtered.vcf -o VEP_out_filtered_MPT_BRIDGE_merged.recode.vcf -filter "(EUR_AF < 0.01 or not EUR_AF) and (UK10KWGS_AF < 0.01 or not UK10KWGS_AF) and (WGS10K_AF < 0.05 or not WGS10K_AF)" --only_matched

#Prepare resulting VCF for reading into R
sed -i 's/#CHROM/CHROM/g' VEP_out_filtered_MPT_BRIDGE_merged.recode.vcf
sed '/^#/ d' VEP_out_filtered_MPT_BRIDGE_merged.recode.vcf > VEP_out_filtered_MPT_BRIDGE_merged.recode_without_header.vcf
mv VEP_out_filtered_MPT_BRIDGE_merged.recode.vcf VEP_out_filtered_MPT_BRIDGE_merged.recode_with_header.vcf
mv VEP_out_filtered_MPT_BRIDGE_merged.recode_without_header.vcf  VEP_out_filtered_MPT_BRIDGE_merged.recode.vcf

######################################################################################################################
##R script to identify samples with each variant in table, count variants and compare frequency in cases vs controls##
######################################################################################################################

#########################
##Heterozygous variants##
#########################

###########################################################################################################
##Following commands used to launch script according to phenotypic subgroup:
#Rscript variant_gtex_vlookup_hets.R All XXX XXX XXX
#Rscript variant_gtex_vlookup_hets.R Phaeochromoctoma ACC XXX XXX
#Rscript variant_gtex_vlookup_hets.R CNS CNSNervesheath XXX XXX
#Rscript variant_gtex_vlookup_hets.R Breast XXX XXX XXX
#Rscript variant_gtex_vlookup_hets.R Haemlymphoid XXX XXX XXX
#Rscript variant_gtex_vlookup_hets.R Colorectal XXX XXX XXX
#Rscript variant_gtex_vlookup_hets.R Oesophagus XXX XXX XXX
#Rscript variant_gtex_vlookup_hets.R Cardiacmyxoma XXX XXX XXX
#Rscript variant_gtex_vlookup_hets.R Lung XXX XXX XXX
#Rscript variant_gtex_vlookup_hets.R Salivarygland XXX XXX XXX
#Rscript variant_gtex_vlookup_hets.R Connectivetissuesofttissuesarcoma XXX XXX XXX
#Rscript variant_gtex_vlookup_hets.R PNSNervesheathbenign PNSNervesheath Nervesheathbenign XXX
#Rscript variant_gtex_vlookup_hets.R Ovary XXX XXX XXX
#Rscript variant_gtex_vlookup_hets.R Pancreas XXX XXX XXX
#Rscript variant_gtex_vlookup_hets.R Pituitary XXX XXX XXX
#Rscript variant_gtex_vlookup_hets.R Prostate XXX XXX XXX
#Rscript variant_gtex_vlookup_hets.R NMSC Melanoma Skinbenign XXX
#Rscript variant_gtex_vlookup_hets.R Smallbowel GINET XXX XXX
#Rscript variant_gtex_vlookup_hets.R Gastric XXX XXX XXX
#Rscript variant_gtex_vlookup_hets.R Testicular XXX XXX XXX
#Rscript variant_gtex_vlookup_hets.R Thyroid XXX XXX XXX
#Rscript variant_gtex_vlookup_hets.R Endometrial Uterineleiomyoma Uterinesarcoma XXX
#Rscript variant_gtex_vlookup_hets.R Haemlymphoid Haemmyeloid Haempolycythaemia Haemthrombocythaemia
############################################################################################################

args = commandArgs(trailingOnly=TRUE)

##Read in necessary files
gtex <- read.csv("/home/jww39/non_coding/gtex/gtex_list.csv", header = T)
var <- read.delim("VEP_out_filtered_MPT_BRIDGE_merged.recode.vcf", header = T)

##Identify samples with phenotype of interest
MPT_table <- read.csv("anonymised_MPT_datasheet_22_11_17.csv")
rownames(MPT_table) <- MPT_table$WGS.ID.1
 MPT_table$Age.1 <- as.numeric(as.character(MPT_table$Age.1)) 
MPT_table$Age.2 <- as.numeric(as.character(MPT_table$Age.2))
MPT_table$Age.3 <- as.numeric(as.character(MPT_table$Age.3))
MPT_table$Age.4 <- as.numeric(as.character(MPT_table$Age.4))
MPT_table$Age.5 <- as.numeric(as.character(MPT_table$Age.5))
MPT_table$Age.6 <- as.numeric(as.character(MPT_table$Age.6))
MPT_table$Age.7 <- as.numeric(as.character(MPT_table$Age.7))

t1 <- args[1]
t2 <- args[2]
t3 <- args[3]
t4 <- args[4]

ting <- rownames(MPT_table[
  
  (
    
    ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
     | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
     | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
     | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
     | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
     | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
     | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
    
    |
      
      ((MPT_table$Single.word.1 == t2 & MPT_table$Age.1 < 70) 
       | (MPT_table$Single.word.2 == t2 & MPT_table$Age.2 < 70) 
       | (MPT_table$Single.word.3 == t2 & MPT_table$Age.3 < 70)
       | (MPT_table$Single.word.4 == t2 & MPT_table$Age.4 < 70)
       | (MPT_table$Single.word.5 == t2 & MPT_table$Age.5 < 70)
       | (MPT_table$Single.word.6 == t2 & MPT_table$Age.6 < 70)
       | (MPT_table$Single.word.7 == t2 & MPT_table$Age.7 < 70))
    
    |
      
      ((MPT_table$Single.word.1 == t3 & MPT_table$Age.1 < 70) 
       | (MPT_table$Single.word.2 == t3 & MPT_table$Age.2 < 70) 
       | (MPT_table$Single.word.3 == t3 & MPT_table$Age.3 < 70)
       | (MPT_table$Single.word.4 == t3 & MPT_table$Age.4 < 70)
       | (MPT_table$Single.word.5 == t3 & MPT_table$Age.5 < 70)
       | (MPT_table$Single.word.6 == t3 & MPT_table$Age.6 < 70)
       | (MPT_table$Single.word.7 == t3 & MPT_table$Age.7 < 70))
    
    |
      
      ((MPT_table$Single.word.1 == t4 & MPT_table$Age.1 < 70) 
       | (MPT_table$Single.word.2 == t4 & MPT_table$Age.2 < 70) 
       | (MPT_table$Single.word.3 == t4 & MPT_table$Age.3 < 70)
       | (MPT_table$Single.word.4 == t4 & MPT_table$Age.4 < 70)
       | (MPT_table$Single.word.5 == t4 & MPT_table$Age.5 < 70)
       | (MPT_table$Single.word.6 == t4 & MPT_table$Age.6 < 70)
       | (MPT_table$Single.word.7 == t4 & MPT_table$Age.7 < 70))
    
  )
  
  & 
    
    (MPT_table$proband == "YES")
  
  &
    
    (MPT_table$exclusion != "EXCLUDE")
  
  ,]) 

indv <- grep("R0*", ting, value = T)
indv <- gsub("_A", "", indv)
indv <- data.frame(indv)


#Put column in variants sheet equivalent to gtex id e.g. 12_58268149_G_T_b37
var_gtex_ID_equiv <- paste(var$CHROM, var$POS, var$REF, var$ALT, "b37", sep = "_")
var <- cbind(var_gtex_ID_equiv, var)


##Matching variants to gtex eqtl and extracting gtex information
gtex_variant_id <- data.frame()

for(i in 1:nrow(var)){
  a <- gtex$variant_id[match(var$var_gtex_ID_equiv[i], gtex$variant_id)]
  a <- as.character(a)
  gtex_variant_id <- append(gtex_variant_id,a)
}
gtex_variant_id <- unlist(gtex_variant_id)

gtex_gene_id <- data.frame()

for(i in 1:nrow(var)){
  a <- gtex$gene_id[match(var$var_gtex_ID_equiv[i], gtex$variant_id)]
  a <- as.character(a)
  gtex_gene_id <- append(gtex_gene_id,a)
}
gtex_gene_id <- unlist(gtex_gene_id)

tss_distance <- data.frame()

for(i in 1:nrow(var)){
  a <- gtex$tss_distance[match(var$var_gtex_ID_equiv[i], gtex$variant_id)]
  a <- as.numeric(a)
  tss_distance <- append(tss_distance,a)
}
tss_distance <- unlist(tss_distance)

tissue <- data.frame()

for(i in 1:nrow(var)){
  a <- which(gtex$variant_id %in% var$var_gtex_ID_equiv[i])
  b <- gtex$tissue[a]
  c <- toString(b)
  tissue <- append(tissue,c)
}

tissue_matrix <- t(data.frame(tissue))
rownames(tissue_matrix) <- 1:nrow(tissue_matrix)
tissue_matrix[tissue_matrix==""] <- NA
library(splitstackshape)
tissue_cols <- data.frame(cSplit(tissue_matrix, 'V1', sep=", ", type.convert=FALSE))
colnames(tissue_cols) <- paste("tissue", 1:ncol(tissue_cols), sep = "_")

##Count variants

#Make list of case samples and control samples
case_list <- indv
if(t1 == "All")case_list <- read.csv("MPT_probands_euro.txt", header = F)

control_list <- read.csv("BRIDGE_euro_controls.txt", header = F)

##Identifying samples (cases and controls) containing each variant
var_gtex_match_genotypes_only <- var[,11:ncol(var)]

samples_raw <- data.frame()
for (i in 1:nrow(var_gtex_match_genotypes_only)) {
  
  a <- unlist(var_gtex_match_genotypes_only[i,])
  b <- grep("^0\\/1", a, perl = TRUE)
  c <- as.vector(names(a[b]))
  d <- paste(c, collapse = ';')
  samples_raw <- append(samples_raw, d)
  
} 
samples <- unlist(samples_raw)

##Extract which samples are among the cases
samples_df <- as.data.frame(samples)
cases_list <- data.frame()
for (i in 1:nrow(samples_df)) {
  
  cases_split <- strsplit(samples[i], ";") # Works. Is a character
  cases_split_unlisted <- unlist(cases_split)
  case_samples <- cases_split_unlisted[(which(cases_split_unlisted %in% case_list[,1]))] #Gives items that match the list
  case_samples_string <- paste(case_samples, collapse = ";")
  cases_list <- as.vector(append(cases_list, case_samples_string))
  
}
cases <- unlist(cases_list)


##Count how many cases have the variant
library(stringr)
cases_df <- as.data.frame(cases)
df_for_cases_variant_count <- data.frame()
for (i in 1:nrow(cases_df)) { 
  count_of_samples_with_variant <- (str_count(cases_df[i,], "[A-Z]0"))
  df_for_cases_variant_count <- as.vector(append(df_for_cases_variant_count, count_of_samples_with_variant))
}
cases_variant_count <- cbind(df_for_cases_variant_count)
cases_variant_count <- unlist(cases_variant_count)

##Extract which samples are among the controls
controls_list <- data.frame()
for (i in 1:nrow(samples_df)) {
  
  controls_split <- strsplit(samples[i], ";") # Works. Is a character
  controls_split_unlisted <- unlist(controls_split)
  controls_samples <- controls_split_unlisted[(which(controls_split_unlisted %in% control_list[,1]))] #Gives items that match the list
  controls_samples_string <- paste(controls_samples, collapse = ";")
  controls_list <- as.vector(append(controls_list, controls_samples_string))
  
}
controls <- unlist(controls_list)

##Count how many controls have the variant
controls_df <- as.data.frame(controls)
df_for_controls_variant_count <- data.frame()
for (i in 1:nrow(controls_df)) { 
  count_of_samples_with_variant <- (str_count(controls_df[i,], "[A-Z][0-9]"))
  df_for_controls_variant_count <- as.vector(append(df_for_controls_variant_count, count_of_samples_with_variant))
}
controls_variant_count <- cbind(df_for_controls_variant_count)
controls_variant_count <- unlist(controls_variant_count)

##Produce a table with variant information and counts
var_gtex_match <- cbind(
  var[,1:9],
  gtex_variant_id, #10
  gtex_gene_id, #11
  tss_distance, #12
  cases, #13
  cases_variant_count, #14
  controls, #15
  controls_variant_count, #16
  tissue_cols[1:ncol(tissue_cols)], #17
  var[,10:ncol(var)]
)

colnames(var_gtex_match)[14] <- paste("cases_variant_count","_n=",(nrow(case_list)),sep = "")
colnames(var_gtex_match)[16] <- paste("controls_variant_count_variant_count","_n=",(nrow(control_list)),sep = "")


##Reduce table down depending on what the tissue/samples of interest are
tissue_equiv <- read.csv("gtex_tissues_with_single_word_equivalents.csv", header = T)

##Get lines with single word of interest in reference table and convert to gtex tissue labels
gtex_tis_1 <- which(tissue_equiv$single_word_equivalent_1 %in% t1)
gtex_tis_2 <- which(tissue_equiv$single_word_equivalent_1 %in% t2)
gtex_tis_3 <- which(tissue_equiv$single_word_equivalent_1 %in% t3)
gtex_tis_4 <- which(tissue_equiv$single_word_equivalent_1 %in% t4)
gtex_tis_all <- c(gtex_tis_1,gtex_tis_2,gtex_tis_3,gtex_tis_4)
gtex_tis_all <- tissue_equiv$gtex.tissue.filename[gtex_tis_all]

##Find which lines in the gtex match table contain the gtex tissue name/s of interest
tissue_lines_of_interest <- data.frame()

for (i in 1:ncol(tissue_cols)) {
  
  a <- which(tissue_cols[,i] %in% gtex_tis_all)
  tissue_lines_of_interest <- append(tissue_lines_of_interest, a)
  
}

tissue_lines_of_interest <- unlist(tissue_lines_of_interest)

##mkae new gtex match table only contining lines which contain the tissue of interest
if(t1 != "All")var_gtex_match <- var_gtex_match[tissue_lines_of_interest,]

##Hypothesis testing

##Fishers exact with fdr adjustment for each individual variant observed
rownum <- nrow(var_gtex_match)
BRIDGE_pvalue <- as.vector(c(), mode = "any")
for(i in 1:rownum){
  
  fishtable = matrix(c(var_gtex_match[i,14],(nrow(case_list) - var_gtex_match[i,14]),var_gtex_match[i,16],(nrow(control_list) - var_gtex_match[i,16])), 
                     nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue <- append(BRIDGE_pvalue, fisherout[[1]])
  
}

BRIDGE_qvalue <- p.adjust(BRIDGE_pvalue, method = "fdr", n = rownum)

##Collate results
var_gtex_match_with_pq <- cbind(BRIDGE_pvalue, BRIDGE_qvalue, var_gtex_match)

var_gtex_match_with_pq$BRIDGE_pvalue <- as.numeric(as.character(var_gtex_match_with_pq$BRIDGE_pvalue))
var_gtex_match_with_pq$BRIDGE_qvalue <- as.numeric(as.character(var_gtex_match_with_pq$BRIDGE_qvalue))
var_gtex_match_sig_index <- which(var_gtex_match_with_pq$BRIDGE_qvalue < 0.1)
var_gtex_match_sig <- var_gtex_match_with_pq[var_gtex_match_sig_index,]

options <- paste(args[1],args[2],args[3],args[4],sep = "_")

##Output results
write.csv(var_gtex_match_sig, paste("/home/jww39/non_coding/gtex/with_internal_af_filter/hets/", options,"_var_gtex_match.csv", sep = ""))
write.csv(var_gtex_match_with_pq, paste("/home/jww39/non_coding/gtex/with_internal_af_filter/hets/", options,"_var_gtex_match_all.csv", sep = ""))

##Per gene analysis

##Counts of variants per gene in cases and controls
genes <- unique(var_gtex_match$gtex_gene_id)
genes <- genes[!is.na(genes)]

variant_per_gene_count_cases <- data.frame()

for(i in 1:length(genes)){
  a <- sum(var_gtex_match[,14][which(var_gtex_match$gtex_gene_id %in% genes[i])]) 
  variant_per_gene_count_cases  <- append(variant_per_gene_count_cases,a)
}

variant_per_gene_count_cases <- t(data.frame(variant_per_gene_count_cases))
gene_names <- as.vector(genes)
rownames(variant_per_gene_count_cases) <- gene_names

variant_per_gene_count_controls <- data.frame()

for(i in 1:length(genes)){
  a <- sum(var_gtex_match[,16][which(var_gtex_match$gtex_gene_id %in% genes[i])])
  variant_per_gene_count_controls  <- append(variant_per_gene_count_controls,a)
}

variant_per_gene_count_controls <- t(data.frame(variant_per_gene_count_controls))
gene_names <- as.vector(genes)
rownames(variant_per_gene_count_controls) <- gene_names

##Fishers exact with fdr adjustment for counts of variants per gene
BRIDGE_pvalue_per_gene <- as.vector(c(), mode = "any")
cases_trials_per_gene <- as.vector(c(), mode = "any")
controls_trials_per_gene <- as.vector(c(), mode = "any")

for(i in 1:nrow(variant_per_gene_count_cases))
{
  fishtable = matrix(
    c(variant_per_gene_count_cases[i,1],
      sum(var_gtex_match$gtex_gene_id %in% genes[i]) * nrow(case_list) - variant_per_gene_count_cases[i,1],
      variant_per_gene_count_controls[i,1],
      sum(var_gtex_match$gtex_gene_id %in% genes[i]) * nrow(control_list) - variant_per_gene_count_controls[i,1]),
    
    nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue_per_gene <- append(BRIDGE_pvalue_per_gene, fisherout[[1]])
  cases_trials_per_gene <- append(cases_trials_per_gene, sum(var_gtex_match$gtex_gene_id %in% genes[i]) * nrow(case_list))
  controls_trials_per_gene <- append(controls_trials_per_gene, sum(var_gtex_match$gtex_gene_id %in% genes[i]) * nrow(control_list))
}

BRIDGE_qvalue_per_gene <- p.adjust(BRIDGE_pvalue_per_gene, method = "fdr", n = 83)

cases_per_gene_prop <- variant_per_gene_count_cases / cases_trials_per_gene
controls_per_gene_prop <- variant_per_gene_count_controls / controls_trials_per_gene

##Counts for individuals with variants per gene for cases and controls
trimmed_cases <- data.frame(var_gtex_match$cases)

indv_with_variant_per_gene_cases <- data.frame()

for(i in 1:length(genes)){
  a <- trimmed_cases[which(var_gtex_match$gtex_gene_id %in% genes[i]),]
  b <- a[which(a != "")]
  c <- paste(b, collapse = ";")
  indv_with_variant_per_gene_cases <- append(indv_with_variant_per_gene_cases, c)
}

indv_with_variant_per_gene_cases <- unlist(indv_with_variant_per_gene_cases)
indv_with_variant_per_gene_cases <- strsplit(indv_with_variant_per_gene_cases, ";", fixed = T)
unique_indv_with_variant_per_gene_cases <- lapply(indv_with_variant_per_gene_cases, unique)

indv_per_gene_count_cases <- data.frame()

for(i in 1:length(unique_indv_with_variant_per_gene_cases)){
  a <- length(which(unique_indv_with_variant_per_gene_cases[[i]] != ""))
  indv_per_gene_count_cases <- append(indv_per_gene_count_cases, a)
}

indv_per_gene_count_cases <- unlist(indv_per_gene_count_cases)

trimmed_controls <- data.frame(var_gtex_match$controls)

indv_with_variant_per_gene_controls <- data.frame()

for(i in 1:length(genes)){
  a <- trimmed_controls[which(var_gtex_match$gtex_gene_id %in% genes[i]),]
  b <- a[which(a != "")]
  c <- paste(b, collapse = ";")
  indv_with_variant_per_gene_controls <- append(indv_with_variant_per_gene_controls, c)
}

indv_with_variant_per_gene_controls <- unlist(indv_with_variant_per_gene_controls)
indv_with_variant_per_gene_controls <- strsplit(indv_with_variant_per_gene_controls, ";", fixed = T)
unique_indv_with_variant_per_gene_controls <- lapply(indv_with_variant_per_gene_controls, unique)

indv_per_gene_count_controls <- data.frame()

for(i in 1:length(unique_indv_with_variant_per_gene_controls)){
  a <- length(which(unique_indv_with_variant_per_gene_controls[[i]] != ""))
  indv_per_gene_count_controls <- append(indv_per_gene_count_controls, a)
}

indv_per_gene_count_controls <- unlist(indv_per_gene_count_controls)

##Fishers exact with fdr adjustment for counts of individuals with variants per gene
BRIDGE_pvalue_per_indv_per_gene <- as.vector(c(), mode = "any")
cases_indv_per_gene <- as.vector(c(), mode = "any")
controls_indv_per_gene <- as.vector(c(), mode = "any")

for(i in 1:length(indv_per_gene_count_cases))
{
  fishtable = matrix(
    c(indv_per_gene_count_cases[i],
      nrow(case_list) - indv_per_gene_count_cases[i],
      indv_per_gene_count_controls[i],
      nrow(control_list) - indv_per_gene_count_controls[i]), 
    nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue_per_indv_per_gene <- append(BRIDGE_pvalue_per_indv_per_gene, fisherout[[1]])
}

BRIDGE_qvalue_per_indv_per_gene <- p.adjust(BRIDGE_pvalue_per_indv_per_gene, method = "fdr", n = 83) 

cases_per_indv_per_gene_prop <- indv_per_gene_count_cases / nrow(case_list)
controls_per_indv_per_gene_prop <- indv_per_gene_count_controls / nrow(control_list)

##Collate results
variants_per_gene_table <- as.data.frame(cbind(rownames(variant_per_gene_count_cases),
                                               
                                               variant_per_gene_count_cases,
                                               cases_trials_per_gene,
                                               cases_per_gene_prop,
                                               variant_per_gene_count_controls,
                                               controls_trials_per_gene,
                                               controls_per_gene_prop,
                                               BRIDGE_pvalue_per_gene,
                                               BRIDGE_qvalue_per_gene,
                                               
                                               indv_per_gene_count_cases,
                                               cases_per_indv_per_gene_prop, 
                                               indv_per_gene_count_controls,
                                               controls_per_indv_per_gene_prop,
                                               BRIDGE_pvalue_per_indv_per_gene,
                                               BRIDGE_qvalue_per_indv_per_gene))

colnames(variants_per_gene_table) <- c("gene", 
                                       
                                       "variant_per_gene_count_cases", 
                                       "cases_trials_per_gene", 
                                       "cases_per_gene_prop",
                                       "variant_per_gene_count_controls",
                                       "controls_trials_per_gene",
                                       "controls_per_gene_prop",
                                       "BRIDGE_pvalue_per_gene",
                                       "BRIDGE_qvalue_per_gene",
                                       
                                       "indv_per_gene_count_cases",
                                       "cases_per_indv_per_gene_prop", 
                                       "indv_per_gene_count_controls",
                                       "controls_per_indv_per_gene_prop",
                                       "BRIDGE_pvalue_per_indv_per_gene",
                                       "BRIDGE_qvalue_per_indv_per_gene")


variants_per_gene_table$BRIDGE_qvalue_per_gene <- as.numeric(as.character(variants_per_gene_table$BRIDGE_qvalue_per_gene))
variants_per_gene_table$BRIDGE_pvalue_per_gene <- as.numeric(as.character(variants_per_gene_table$BRIDGE_pvalue_per_gene))
variants_per_gene_table$BRIDGE_qvalue_per_indv_per_gene <- as.numeric(as.character(variants_per_gene_table$BRIDGE_qvalue_per_indv_per_gene))
variants_per_gene_table$BRIDGE_pvalue_per_indv_per_gene <- as.numeric(as.character(variants_per_gene_table$BRIDGE_pvalue_per_indv_per_gene))

per_gene_sig_index <- which(variants_per_gene_table$BRIDGE_qvalue_per_gene < 0.1)
per_indv_per_gene_sig_index <- which(variants_per_gene_table$BRIDGE_qvalue_per_indv_per_gene < 0.1) 
sig_index <- c(per_gene_sig_index, per_indv_per_gene_sig_index)
sig_index <- unique(sig_index)

variants_per_gene_table_sig <- variants_per_gene_table[sig_index,]

##Output results
options <- paste(args[1],args[2],args[3],args[4],sep = "_")
write.csv(variants_per_gene_table_sig, paste("/home/jww39/non_coding/gtex/with_internal_af_filter/hets/", options,"_variants_per_gene_table.csv", sep = "")) 
write.csv(variants_per_gene_table, paste("/home/jww39/non_coding/gtex/with_internal_af_filter/hets/", options,"_variants_per_gene_table_all.csv", sep = "")) 

#######################
##Homozygous variants##
#######################

####################################################################################################################
##Following commands used to launch script according to phenotypic subgroup:
#Rscript variant_gtex_vlookup_homs.R All XXX XXX XXX
#Rscript variant_gtex_vlookup_homs.R Phaeochromoctoma ACC XXX XXX
#Rscript variant_gtex_vlookup_homs.R CNS CNSNervesheath XXX XXX
#Rscript variant_gtex_vlookup_homs.R Breast XXX XXX XXX
#Rscript variant_gtex_vlookup_homs.R Haemlymphoid XXX XXX XXX
#Rscript variant_gtex_vlookup_homs.R Colorectal XXX XXX XXX
#Rscript variant_gtex_vlookup_homs.R Oesophagus XXX XXX XXX
#Rscript variant_gtex_vlookup_homs.R Cardiacmyxoma XXX XXX XXX
#Rscript variant_gtex_vlookup_homs.R Lung XXX XXX XXX
#Rscript variant_gtex_vlookup_homs.R Salivarygland XXX XXX XXX
#Rscript variant_gtex_vlookup_homs.R Connectivetissuesofttissuesarcoma XXX XXX XXX
#Rscript variant_gtex_vlookup_homs.R PNSNervesheathbenign PNSNervesheath Nervesheathbenign XXX
#Rscript variant_gtex_vlookup_homs.R Ovary XXX XXX XXX
#Rscript variant_gtex_vlookup_homs.R Pancreas XXX XXX XXX
#Rscript variant_gtex_vlookup_homs.R Pituitary XXX XXX XXX
#Rscript variant_gtex_vlookup_homs.R Prostate XXX XXX XXX
#Rscript variant_gtex_vlookup_homs.R NMSC Melanoma Skinbenign XXX
#Rscript variant_gtex_vlookup_homs.R Smallbowel GINET XXX XXX
#Rscript variant_gtex_vlookup_homs.R Gastric XXX XXX XXX
#Rscript variant_gtex_vlookup_homs.R Testicular XXX XXX XXX
#Rscript variant_gtex_vlookup_homs.R Thyroid XXX XXX XXX
#Rscript variant_gtex_vlookup_homs.R Endometrial Uterineleiomyoma Uterinesarcoma XXX
#Rscript variant_gtex_vlookup_homs.R Haemlymphoid Haemmyeloid Haempolycythaemia Haemthrombocythaemia
###############################################################################################################

##Read in necessary files
gtex <- read.csv("/home/jww39/non_coding/gtex/gtex_list.csv", header = T)
var <- read.delim("VEP_out_filtered_MPT_BRIDGE_merged.recode.vcf", header = T)

##Identify samples with phenotype of interest
MPT_table <- read.csv("anonymised_MPT_datasheet_22_11_17.csv")
rownames(MPT_table) <- MPT_table$WGS.ID.1
 MPT_table$Age.1 <- as.numeric(as.character(MPT_table$Age.1)) 
MPT_table$Age.2 <- as.numeric(as.character(MPT_table$Age.2))
MPT_table$Age.3 <- as.numeric(as.character(MPT_table$Age.3))
MPT_table$Age.4 <- as.numeric(as.character(MPT_table$Age.4))
MPT_table$Age.5 <- as.numeric(as.character(MPT_table$Age.5))
MPT_table$Age.6 <- as.numeric(as.character(MPT_table$Age.6))
MPT_table$Age.7 <- as.numeric(as.character(MPT_table$Age.7))

t1 <- args[1]
t2 <- args[2]
t3 <- args[3]
t4 <- args[4]

ting <- rownames(MPT_table[
  
  (
    
    ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
     | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
     | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
     | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
     | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
     | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
     | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
    
    |
      
      ((MPT_table$Single.word.1 == t2 & MPT_table$Age.1 < 70) 
       | (MPT_table$Single.word.2 == t2 & MPT_table$Age.2 < 70) 
       | (MPT_table$Single.word.3 == t2 & MPT_table$Age.3 < 70)
       | (MPT_table$Single.word.4 == t2 & MPT_table$Age.4 < 70)
       | (MPT_table$Single.word.5 == t2 & MPT_table$Age.5 < 70)
       | (MPT_table$Single.word.6 == t2 & MPT_table$Age.6 < 70)
       | (MPT_table$Single.word.7 == t2 & MPT_table$Age.7 < 70))
    
    |
      
      ((MPT_table$Single.word.1 == t3 & MPT_table$Age.1 < 70) 
       | (MPT_table$Single.word.2 == t3 & MPT_table$Age.2 < 70) 
       | (MPT_table$Single.word.3 == t3 & MPT_table$Age.3 < 70)
       | (MPT_table$Single.word.4 == t3 & MPT_table$Age.4 < 70)
       | (MPT_table$Single.word.5 == t3 & MPT_table$Age.5 < 70)
       | (MPT_table$Single.word.6 == t3 & MPT_table$Age.6 < 70)
       | (MPT_table$Single.word.7 == t3 & MPT_table$Age.7 < 70))
    
    |
      
      ((MPT_table$Single.word.1 == t4 & MPT_table$Age.1 < 70) 
       | (MPT_table$Single.word.2 == t4 & MPT_table$Age.2 < 70) 
       | (MPT_table$Single.word.3 == t4 & MPT_table$Age.3 < 70)
       | (MPT_table$Single.word.4 == t4 & MPT_table$Age.4 < 70)
       | (MPT_table$Single.word.5 == t4 & MPT_table$Age.5 < 70)
       | (MPT_table$Single.word.6 == t4 & MPT_table$Age.6 < 70)
       | (MPT_table$Single.word.7 == t4 & MPT_table$Age.7 < 70))
    
  )
  
  & 
    
    (MPT_table$proband == "YES")
  
  &
    
    (MPT_table$exclusion != "EXCLUDE")
  
  ,]) 

indv <- grep("R0*", ting, value = T)
indv <- gsub("_A", "", indv)
indv <- data.frame(indv)


#Put column in variants sheet equivalent to gtex id e.g. 12_58268149_G_T_b37
var_gtex_ID_equiv <- paste(var$CHROM, var$POS, var$REF, var$ALT, "b37", sep = "_")
var <- cbind(var_gtex_ID_equiv, var)


##Matching variants to gtex eqtl and extracting gtex information
gtex_variant_id <- data.frame()

for(i in 1:nrow(var)){
  a <- gtex$variant_id[match(var$var_gtex_ID_equiv[i], gtex$variant_id)]
  a <- as.character(a)
  gtex_variant_id <- append(gtex_variant_id,a)
}
gtex_variant_id <- unlist(gtex_variant_id)

gtex_gene_id <- data.frame()

for(i in 1:nrow(var)){
  a <- gtex$gene_id[match(var$var_gtex_ID_equiv[i], gtex$variant_id)]
  a <- as.character(a)
  gtex_gene_id <- append(gtex_gene_id,a)
}
gtex_gene_id <- unlist(gtex_gene_id)

tss_distance <- data.frame()

for(i in 1:nrow(var)){
  a <- gtex$tss_distance[match(var$var_gtex_ID_equiv[i], gtex$variant_id)]
  a <- as.numeric(a)
  tss_distance <- append(tss_distance,a)
}
tss_distance <- unlist(tss_distance)

tissue <- data.frame()

for(i in 1:nrow(var)){
  a <- which(gtex$variant_id %in% var$var_gtex_ID_equiv[i])
  b <- gtex$tissue[a]
  c <- toString(b)
  tissue <- append(tissue,c)
}

tissue_matrix <- t(data.frame(tissue))
rownames(tissue_matrix) <- 1:nrow(tissue_matrix)
tissue_matrix[tissue_matrix==""] <- NA
library(splitstackshape)
tissue_cols <- data.frame(cSplit(tissue_matrix, 'V1', sep=", ", type.convert=FALSE))
colnames(tissue_cols) <- paste("tissue", 1:ncol(tissue_cols), sep = "_")

##Count variants

#Make list of case samples and control samples
case_list <- indv
if(t1 == "All")case_list <- read.csv("MPT_probands_euro.txt", header = F)

control_list <- read.csv("BRIDGE_euro_controls.txt", header = F)

##Identifying samples (cases and controls) containing each variant
var_gtex_match_genotypes_only <- var[,11:ncol(var)]

samples_raw <- data.frame()
for (i in 1:nrow(var_gtex_match_genotypes_only)) {
  
  a <- unlist(var_gtex_match_genotypes_only[i,])
  b <- grep("^1\\/1", a, perl = TRUE)
  c <- as.vector(names(a[b]))
  d <- paste(c, collapse = ';')
  samples_raw <- append(samples_raw, d)
}  
samples <- unlist(samples_raw)

##Extract which samples are among the cases
samples_df <- as.data.frame(samples)
cases_list <- data.frame()
for (i in 1:nrow(samples_df)) {
  
  cases_split <- strsplit(samples[i], ";") # Works. Is a character
  cases_split_unlisted <- unlist(cases_split)
  case_samples <- cases_split_unlisted[(which(cases_split_unlisted %in% case_list[,1]))] #Gives items that match the list
  case_samples_string <- paste(case_samples, collapse = ";")
  cases_list <- as.vector(append(cases_list, case_samples_string))
  
}
cases <- unlist(cases_list)


##Count how many cases have the variant
library(stringr)
cases_df <- as.data.frame(cases)
df_for_cases_variant_count <- data.frame()
for (i in 1:nrow(cases_df)) { 
  count_of_samples_with_variant <- (str_count(cases_df[i,], "[A-Z]0"))
  df_for_cases_variant_count <- as.vector(append(df_for_cases_variant_count, count_of_samples_with_variant))
}
cases_variant_count <- cbind(df_for_cases_variant_count)
cases_variant_count <- unlist(cases_variant_count)

##Extract which samples are among the controls
controls_list <- data.frame()
for (i in 1:nrow(samples_df)) {
  
  controls_split <- strsplit(samples[i], ";") # Works. Is a character
  controls_split_unlisted <- unlist(controls_split)
  controls_samples <- controls_split_unlisted[(which(controls_split_unlisted %in% control_list[,1]))] #Gives items that match the list
  controls_samples_string <- paste(controls_samples, collapse = ";")
  controls_list <- as.vector(append(controls_list, controls_samples_string))
  
}
controls <- unlist(controls_list)

##Count how many controls have the variant
controls_df <- as.data.frame(controls)
df_for_controls_variant_count <- data.frame()
for (i in 1:nrow(controls_df)) { 
  count_of_samples_with_variant <- (str_count(controls_df[i,], "[A-Z][0-9]"))
  df_for_controls_variant_count <- as.vector(append(df_for_controls_variant_count, count_of_samples_with_variant))
}
controls_variant_count <- cbind(df_for_controls_variant_count)
controls_variant_count <- unlist(controls_variant_count)

##Produce a table with variant information and counts
var_gtex_match <- cbind(
  var[,1:9],
  gtex_variant_id, #10
  gtex_gene_id, #11
  tss_distance, #12
  cases, #13
  cases_variant_count, #14
  controls, #15
  controls_variant_count, #16
  tissue_cols[1:ncol(tissue_cols)], #17
  var[,10:ncol(var)]
)

colnames(var_gtex_match)[14] <- paste("cases_variant_count","_n=",(nrow(case_list)),sep = "")
colnames(var_gtex_match)[16] <- paste("controls_variant_count_variant_count","_n=",(nrow(control_list)),sep = "")


##Reduce table down depending on what the tissue/samples of interest are
tissue_equiv <- read.csv("gtex_tissues_with_single_word_equivalents.csv", header = T)

##Get lines with single word of interest in reference table and convert to gtex tissue labels
gtex_tis_1 <- which(tissue_equiv$single_word_equivalent_1 %in% t1)
gtex_tis_2 <- which(tissue_equiv$single_word_equivalent_1 %in% t2)
gtex_tis_3 <- which(tissue_equiv$single_word_equivalent_1 %in% t3)
gtex_tis_4 <- which(tissue_equiv$single_word_equivalent_1 %in% t4)
gtex_tis_all <- c(gtex_tis_1,gtex_tis_2,gtex_tis_3,gtex_tis_4)
gtex_tis_all <- tissue_equiv$gtex.tissue.filename[gtex_tis_all]

##Find which lines in the gtex match table contain the gtex tissue name/s of interest
tissue_lines_of_interest <- data.frame()

for (i in 1:ncol(tissue_cols)) {
  
  a <- which(tissue_cols[,i] %in% gtex_tis_all)
  tissue_lines_of_interest <- append(tissue_lines_of_interest, a)
  
}

tissue_lines_of_interest <- unlist(tissue_lines_of_interest)

##mkae new gtex match table only contining lines which contain the tissue of interest
if(t1 != "All")var_gtex_match <- var_gtex_match[tissue_lines_of_interest,]

##Hypothesis testing

##Fishers exact with fdr adjustment for each individual variant observed
rownum <- nrow(var_gtex_match)
BRIDGE_pvalue <- as.vector(c(), mode = "any")
for(i in 1:rownum){
  
  fishtable = matrix(c(var_gtex_match[i,14],(nrow(case_list) - var_gtex_match[i,14]),var_gtex_match[i,16],(nrow(control_list) - var_gtex_match[i,16])), 
                     nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue <- append(BRIDGE_pvalue, fisherout[[1]])
  
}

BRIDGE_qvalue <- p.adjust(BRIDGE_pvalue, method = "fdr", n = rownum)

##Collate results
var_gtex_match_with_pq <- cbind(BRIDGE_pvalue, BRIDGE_qvalue, var_gtex_match)

var_gtex_match_with_pq$BRIDGE_pvalue <- as.numeric(as.character(var_gtex_match_with_pq$BRIDGE_pvalue))
var_gtex_match_with_pq$BRIDGE_qvalue <- as.numeric(as.character(var_gtex_match_with_pq$BRIDGE_qvalue))
var_gtex_match_sig_index <- which(var_gtex_match_with_pq$BRIDGE_qvalue < 0.1)
var_gtex_match_sig <- var_gtex_match_with_pq[var_gtex_match_sig_index,]

options <- paste(args[1],args[2],args[3],args[4],sep = "_")

##Output results
write.csv(var_gtex_match_sig, paste("/home/jww39/non_coding/gtex/with_internal_af_filter/hets/", options,"_var_gtex_match.csv", sep = ""))
write.csv(var_gtex_match_with_pq, paste("/home/jww39/non_coding/gtex/with_internal_af_filter/hets/", options,"_var_gtex_match_all.csv", sep = ""))

##Per gene analysis

##Counts of variants per gene in cases and controls
genes <- unique(var_gtex_match$gtex_gene_id)
genes <- genes[!is.na(genes)]

variant_per_gene_count_cases <- data.frame()

for(i in 1:length(genes)){
  a <- sum(var_gtex_match[,14][which(var_gtex_match$gtex_gene_id %in% genes[i])]) 
  variant_per_gene_count_cases  <- append(variant_per_gene_count_cases,a)
}

variant_per_gene_count_cases <- t(data.frame(variant_per_gene_count_cases))
gene_names <- as.vector(genes)
rownames(variant_per_gene_count_cases) <- gene_names

variant_per_gene_count_controls <- data.frame()

for(i in 1:length(genes)){
  a <- sum(var_gtex_match[,16][which(var_gtex_match$gtex_gene_id %in% genes[i])])
  variant_per_gene_count_controls  <- append(variant_per_gene_count_controls,a)
}

variant_per_gene_count_controls <- t(data.frame(variant_per_gene_count_controls))
gene_names <- as.vector(genes)
rownames(variant_per_gene_count_controls) <- gene_names

##Fishers exact with fdr adjustment for counts of variants per gene
BRIDGE_pvalue_per_gene <- as.vector(c(), mode = "any")
cases_trials_per_gene <- as.vector(c(), mode = "any")
controls_trials_per_gene <- as.vector(c(), mode = "any")

for(i in 1:nrow(variant_per_gene_count_cases))
{
  fishtable = matrix(
    c(variant_per_gene_count_cases[i,1],
      sum(var_gtex_match$gtex_gene_id %in% genes[i]) * nrow(case_list) - variant_per_gene_count_cases[i,1],
      variant_per_gene_count_controls[i,1],
      sum(var_gtex_match$gtex_gene_id %in% genes[i]) * nrow(control_list) - variant_per_gene_count_controls[i,1]),
    
    nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue_per_gene <- append(BRIDGE_pvalue_per_gene, fisherout[[1]])
  cases_trials_per_gene <- append(cases_trials_per_gene, sum(var_gtex_match$gtex_gene_id %in% genes[i]) * nrow(case_list))
  controls_trials_per_gene <- append(controls_trials_per_gene, sum(var_gtex_match$gtex_gene_id %in% genes[i]) * nrow(control_list))
}

BRIDGE_qvalue_per_gene <- p.adjust(BRIDGE_pvalue_per_gene, method = "fdr", n = 83)

cases_per_gene_prop <- variant_per_gene_count_cases / cases_trials_per_gene
controls_per_gene_prop <- variant_per_gene_count_controls / controls_trials_per_gene

##Counts for individuals with variants per gene for cases and controls
trimmed_cases <- data.frame(var_gtex_match$cases)

indv_with_variant_per_gene_cases <- data.frame()

for(i in 1:length(genes)){
  a <- trimmed_cases[which(var_gtex_match$gtex_gene_id %in% genes[i]),]
  b <- a[which(a != "")]
  c <- paste(b, collapse = ";")
  indv_with_variant_per_gene_cases <- append(indv_with_variant_per_gene_cases, c)
}

indv_with_variant_per_gene_cases <- unlist(indv_with_variant_per_gene_cases)
indv_with_variant_per_gene_cases <- strsplit(indv_with_variant_per_gene_cases, ";", fixed = T)
unique_indv_with_variant_per_gene_cases <- lapply(indv_with_variant_per_gene_cases, unique)

indv_per_gene_count_cases <- data.frame()

for(i in 1:length(unique_indv_with_variant_per_gene_cases)){
  a <- length(which(unique_indv_with_variant_per_gene_cases[[i]] != ""))
  indv_per_gene_count_cases <- append(indv_per_gene_count_cases, a)
}

indv_per_gene_count_cases <- unlist(indv_per_gene_count_cases)

trimmed_controls <- data.frame(var_gtex_match$controls)

indv_with_variant_per_gene_controls <- data.frame()

for(i in 1:length(genes)){
  a <- trimmed_controls[which(var_gtex_match$gtex_gene_id %in% genes[i]),]
  b <- a[which(a != "")]
  c <- paste(b, collapse = ";")
  indv_with_variant_per_gene_controls <- append(indv_with_variant_per_gene_controls, c)
}

indv_with_variant_per_gene_controls <- unlist(indv_with_variant_per_gene_controls)
indv_with_variant_per_gene_controls <- strsplit(indv_with_variant_per_gene_controls, ";", fixed = T)
unique_indv_with_variant_per_gene_controls <- lapply(indv_with_variant_per_gene_controls, unique)

indv_per_gene_count_controls <- data.frame()

for(i in 1:length(unique_indv_with_variant_per_gene_controls)){
  a <- length(which(unique_indv_with_variant_per_gene_controls[[i]] != ""))
  indv_per_gene_count_controls <- append(indv_per_gene_count_controls, a)
}

indv_per_gene_count_controls <- unlist(indv_per_gene_count_controls)

##Fishers exact with fdr adjustment for counts of individuals with variants per gene
BRIDGE_pvalue_per_indv_per_gene <- as.vector(c(), mode = "any")
cases_indv_per_gene <- as.vector(c(), mode = "any")
controls_indv_per_gene <- as.vector(c(), mode = "any")

for(i in 1:length(indv_per_gene_count_cases))
{
  fishtable = matrix(
    c(indv_per_gene_count_cases[i],
      nrow(case_list) - indv_per_gene_count_cases[i],
      indv_per_gene_count_controls[i],
      nrow(control_list) - indv_per_gene_count_controls[i]), 
    nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue_per_indv_per_gene <- append(BRIDGE_pvalue_per_indv_per_gene, fisherout[[1]])
}

BRIDGE_qvalue_per_indv_per_gene <- p.adjust(BRIDGE_pvalue_per_indv_per_gene, method = "fdr", n = 83) 

cases_per_indv_per_gene_prop <- indv_per_gene_count_cases / nrow(case_list)
controls_per_indv_per_gene_prop <- indv_per_gene_count_controls / nrow(control_list)

##Collate results
variants_per_gene_table <- as.data.frame(cbind(rownames(variant_per_gene_count_cases),
                                               
                                               variant_per_gene_count_cases,
                                               cases_trials_per_gene,
                                               cases_per_gene_prop,
                                               variant_per_gene_count_controls,
                                               controls_trials_per_gene,
                                               controls_per_gene_prop,
                                               BRIDGE_pvalue_per_gene,
                                               BRIDGE_qvalue_per_gene,
                                               
                                               indv_per_gene_count_cases,
                                               cases_per_indv_per_gene_prop, 
                                               indv_per_gene_count_controls,
                                               controls_per_indv_per_gene_prop,
                                               BRIDGE_pvalue_per_indv_per_gene,
                                               BRIDGE_qvalue_per_indv_per_gene))

colnames(variants_per_gene_table) <- c("gene", 
                                       
                                       "variant_per_gene_count_cases", 
                                       "cases_trials_per_gene", 
                                       "cases_per_gene_prop",
                                       "variant_per_gene_count_controls",
                                       "controls_trials_per_gene",
                                       "controls_per_gene_prop",
                                       "BRIDGE_pvalue_per_gene",
                                       "BRIDGE_qvalue_per_gene",
                                       
                                       "indv_per_gene_count_cases",
                                       "cases_per_indv_per_gene_prop", 
                                       "indv_per_gene_count_controls",
                                       "controls_per_indv_per_gene_prop",
                                       "BRIDGE_pvalue_per_indv_per_gene",
                                       "BRIDGE_qvalue_per_indv_per_gene")


variants_per_gene_table$BRIDGE_qvalue_per_gene <- as.numeric(as.character(variants_per_gene_table$BRIDGE_qvalue_per_gene))
variants_per_gene_table$BRIDGE_pvalue_per_gene <- as.numeric(as.character(variants_per_gene_table$BRIDGE_pvalue_per_gene))
variants_per_gene_table$BRIDGE_qvalue_per_indv_per_gene <- as.numeric(as.character(variants_per_gene_table$BRIDGE_qvalue_per_indv_per_gene))
variants_per_gene_table$BRIDGE_pvalue_per_indv_per_gene <- as.numeric(as.character(variants_per_gene_table$BRIDGE_pvalue_per_indv_per_gene))

per_gene_sig_index <- which(variants_per_gene_table$BRIDGE_qvalue_per_gene < 0.1)
per_indv_per_gene_sig_index <- which(variants_per_gene_table$BRIDGE_qvalue_per_indv_per_gene < 0.1) 
sig_index <- c(per_gene_sig_index, per_indv_per_gene_sig_index)
sig_index <- unique(sig_index)

variants_per_gene_table_sig <- variants_per_gene_table[sig_index,]

##Output results
options <- paste(args[1],args[2],args[3],args[4],sep = "_")
write.csv(variants_per_gene_table_sig, paste("/home/jww39/non_coding/gtex/with_internal_af_filter/homs/", options,"_variants_per_gene_table.csv", sep = "")) 
write.csv(variants_per_gene_table, paste("/home/jww39/non_coding/gtex/with_internal_af_filter/homs/", options,"_variants_per_gene_table_all.csv", sep = "")) 

#########################################################
##Heterozygous variants summed with homozygous variants##
#########################################################

####################################################################################################################
##Following commands used to launch script according to phenotypic subgroup:
#Rscript variant_gtex_vlookup_hets_homs.R All XXX XXX XXX
#Rscript variant_gtex_vlookup_hets_homs.R Phaeochromoctoma ACC XXX XXX
#Rscript variant_gtex_vlookup_hets_homs.R CNS CNSNervesheath XXX XXX
#Rscript variant_gtex_vlookup_hets_homs.R Breast XXX XXX XXX
#Rscript variant_gtex_vlookup_hets_homs.R Haemlymphoid XXX XXX XXX
#Rscript variant_gtex_vlookup_hets_homs.R Colorectal XXX XXX XXX
#Rscript variant_gtex_vlookup_hets_homs.R Oesophagus XXX XXX XXX
#Rscript variant_gtex_vlookup_hets_homs.R Cardiacmyxoma XXX XXX XXX
#Rscript variant_gtex_vlookup_hets_homs.R Lung XXX XXX XXX
#Rscript variant_gtex_vlookup_hets_homs.R Salivarygland XXX XXX XXX
#Rscript variant_gtex_vlookup_hets_homs.R Connectivetissuesofttissuesarcoma XXX XXX XXX
#Rscript variant_gtex_vlookup_hets_homs.R PNSNervesheathbenign PNSNervesheath Nervesheathbenign XXX
#Rscript variant_gtex_vlookup_hets_homs.R Ovary XXX XXX XXX
#Rscript variant_gtex_vlookup_hets_homs.R Pancreas XXX XXX XXX
#Rscript variant_gtex_vlookup_hets_homs.R Pituitary XXX XXX XXX
#Rscript variant_gtex_vlookup_hets_homs.R Prostate XXX XXX XXX
#Rscript variant_gtex_vlookup_hets_homs.R NMSC Melanoma Skinbenign XXX
#Rscript variant_gtex_vlookup_hets_homs.R Smallbowel GINET XXX XXX
#Rscript variant_gtex_vlookup_hets_homs.R Gastric XXX XXX XXX
#Rscript variant_gtex_vlookup_hets_homs.R Testicular XXX XXX XXX
#Rscript variant_gtex_vlookup_hets_homs.R Thyroid XXX XXX XXX
#Rscript variant_gtex_vlookup_hets_homs.R Endometrial Uterineleiomyoma Uterinesarcoma XXX
#Rscript variant_gtex_vlookup_hets_homs.R Haemlymphoid Haemmyeloid Haempolycythaemia Haemthrombocythaemia
###############################################################################################################

##Read in necessary files
gtex <- read.csv("/home/jww39/non_coding/gtex/gtex_list.csv", header = T)
var <- read.delim("VEP_out_filtered_MPT_BRIDGE_merged.recode.vcf", header = T)

##Identify samples with phenotype of interest
MPT_table <- read.csv("anonymised_MPT_datasheet_22_11_17.csv")
rownames(MPT_table) <- MPT_table$WGS.ID.1
 MPT_table$Age.1 <- as.numeric(as.character(MPT_table$Age.1)) 
MPT_table$Age.2 <- as.numeric(as.character(MPT_table$Age.2))
MPT_table$Age.3 <- as.numeric(as.character(MPT_table$Age.3))
MPT_table$Age.4 <- as.numeric(as.character(MPT_table$Age.4))
MPT_table$Age.5 <- as.numeric(as.character(MPT_table$Age.5))
MPT_table$Age.6 <- as.numeric(as.character(MPT_table$Age.6))
MPT_table$Age.7 <- as.numeric(as.character(MPT_table$Age.7))

t1 <- args[1]
t2 <- args[2]
t3 <- args[3]
t4 <- args[4]

ting <- rownames(MPT_table[
  
  (
    
    ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
     | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
     | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
     | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
     | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
     | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
     | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
    
    |
      
      ((MPT_table$Single.word.1 == t2 & MPT_table$Age.1 < 70) 
       | (MPT_table$Single.word.2 == t2 & MPT_table$Age.2 < 70) 
       | (MPT_table$Single.word.3 == t2 & MPT_table$Age.3 < 70)
       | (MPT_table$Single.word.4 == t2 & MPT_table$Age.4 < 70)
       | (MPT_table$Single.word.5 == t2 & MPT_table$Age.5 < 70)
       | (MPT_table$Single.word.6 == t2 & MPT_table$Age.6 < 70)
       | (MPT_table$Single.word.7 == t2 & MPT_table$Age.7 < 70))
    
    |
      
      ((MPT_table$Single.word.1 == t3 & MPT_table$Age.1 < 70) 
       | (MPT_table$Single.word.2 == t3 & MPT_table$Age.2 < 70) 
       | (MPT_table$Single.word.3 == t3 & MPT_table$Age.3 < 70)
       | (MPT_table$Single.word.4 == t3 & MPT_table$Age.4 < 70)
       | (MPT_table$Single.word.5 == t3 & MPT_table$Age.5 < 70)
       | (MPT_table$Single.word.6 == t3 & MPT_table$Age.6 < 70)
       | (MPT_table$Single.word.7 == t3 & MPT_table$Age.7 < 70))
    
    |
      
      ((MPT_table$Single.word.1 == t4 & MPT_table$Age.1 < 70) 
       | (MPT_table$Single.word.2 == t4 & MPT_table$Age.2 < 70) 
       | (MPT_table$Single.word.3 == t4 & MPT_table$Age.3 < 70)
       | (MPT_table$Single.word.4 == t4 & MPT_table$Age.4 < 70)
       | (MPT_table$Single.word.5 == t4 & MPT_table$Age.5 < 70)
       | (MPT_table$Single.word.6 == t4 & MPT_table$Age.6 < 70)
       | (MPT_table$Single.word.7 == t4 & MPT_table$Age.7 < 70))
    
  )
  
  & 
    
    (MPT_table$proband == "YES")
  
  &
    
    (MPT_table$exclusion != "EXCLUDE")
  
  ,]) 

indv <- grep("R0*", ting, value = T)
indv <- gsub("_A", "", indv)
indv <- data.frame(indv)


#Put column in variants sheet equivalent to gtex id e.g. 12_58268149_G_T_b37
var_gtex_ID_equiv <- paste(var$CHROM, var$POS, var$REF, var$ALT, "b37", sep = "_")
var <- cbind(var_gtex_ID_equiv, var)


##Matching variants to gtex eqtl and extracting gtex information
gtex_variant_id <- data.frame()

for(i in 1:nrow(var)){
  a <- gtex$variant_id[match(var$var_gtex_ID_equiv[i], gtex$variant_id)]
  a <- as.character(a)
  gtex_variant_id <- append(gtex_variant_id,a)
}
gtex_variant_id <- unlist(gtex_variant_id)

gtex_gene_id <- data.frame()

for(i in 1:nrow(var)){
  a <- gtex$gene_id[match(var$var_gtex_ID_equiv[i], gtex$variant_id)]
  a <- as.character(a)
  gtex_gene_id <- append(gtex_gene_id,a)
}
gtex_gene_id <- unlist(gtex_gene_id)

tss_distance <- data.frame()

for(i in 1:nrow(var)){
  a <- gtex$tss_distance[match(var$var_gtex_ID_equiv[i], gtex$variant_id)]
  a <- as.numeric(a)
  tss_distance <- append(tss_distance,a)
}
tss_distance <- unlist(tss_distance)

tissue <- data.frame()

for(i in 1:nrow(var)){
  a <- which(gtex$variant_id %in% var$var_gtex_ID_equiv[i])
  b <- gtex$tissue[a]
  c <- toString(b)
  tissue <- append(tissue,c)
}

tissue_matrix <- t(data.frame(tissue))
rownames(tissue_matrix) <- 1:nrow(tissue_matrix)
tissue_matrix[tissue_matrix==""] <- NA
library(splitstackshape)
tissue_cols <- data.frame(cSplit(tissue_matrix, 'V1', sep=", ", type.convert=FALSE))
colnames(tissue_cols) <- paste("tissue", 1:ncol(tissue_cols), sep = "_")

##Count variants

#Make list of case samples and control samples
case_list <- indv
if(t1 == "All")case_list <- read.csv("MPT_probands_euro.txt", header = F)

control_list <- read.csv("BRIDGE_euro_controls.txt", header = F)

##Identifying samples (cases and controls) containing each variant
var_gtex_match_genotypes_only <- var[,11:ncol(var)]

samples_raw <- data.frame()
for (i in 1:nrow(var_gtex_match_genotypes_only)) {
  
  a <- unlist(var_gtex_match_genotypes_only[i,])
  b <- grep("^[0-9]\\/[0-9]", a, perl = TRUE)
  c <- as.vector(names(a[b]))
  d <- paste(c, collapse = ';')
  samples_raw <- append(samples_raw, d)
} 
samples <- unlist(samples_raw)

##Extract which samples are among the cases
samples_df <- as.data.frame(samples)
cases_list <- data.frame()
for (i in 1:nrow(samples_df)) {
  
  cases_split <- strsplit(samples[i], ";") # Works. Is a character
  cases_split_unlisted <- unlist(cases_split)
  case_samples <- cases_split_unlisted[(which(cases_split_unlisted %in% case_list[,1]))] #Gives items that match the list
  case_samples_string <- paste(case_samples, collapse = ";")
  cases_list <- as.vector(append(cases_list, case_samples_string))
  
}
cases <- unlist(cases_list)


##Count how many cases have the variant
library(stringr)
cases_df <- as.data.frame(cases)
df_for_cases_variant_count <- data.frame()
for (i in 1:nrow(cases_df)) { 
  count_of_samples_with_variant <- (str_count(cases_df[i,], "[A-Z]0"))
  df_for_cases_variant_count <- as.vector(append(df_for_cases_variant_count, count_of_samples_with_variant))
}
cases_variant_count <- cbind(df_for_cases_variant_count)
cases_variant_count <- unlist(cases_variant_count)

##Extract which samples are among the controls
controls_list <- data.frame()
for (i in 1:nrow(samples_df)) {
  
  controls_split <- strsplit(samples[i], ";") # Works. Is a character
  controls_split_unlisted <- unlist(controls_split)
  controls_samples <- controls_split_unlisted[(which(controls_split_unlisted %in% control_list[,1]))] #Gives items that match the list
  controls_samples_string <- paste(controls_samples, collapse = ";")
  controls_list <- as.vector(append(controls_list, controls_samples_string))
  
}
controls <- unlist(controls_list)

##Count how many controls have the variant
controls_df <- as.data.frame(controls)
df_for_controls_variant_count <- data.frame()
for (i in 1:nrow(controls_df)) { 
  count_of_samples_with_variant <- (str_count(controls_df[i,], "[A-Z][0-9]"))
  df_for_controls_variant_count <- as.vector(append(df_for_controls_variant_count, count_of_samples_with_variant))
}
controls_variant_count <- cbind(df_for_controls_variant_count)
controls_variant_count <- unlist(controls_variant_count)

##Produce a table with variant information and counts
var_gtex_match <- cbind(
  var[,1:9],
  gtex_variant_id, #10
  gtex_gene_id, #11
  tss_distance, #12
  cases, #13
  cases_variant_count, #14
  controls, #15
  controls_variant_count, #16
  tissue_cols[1:ncol(tissue_cols)], #17
  var[,10:ncol(var)]
)

colnames(var_gtex_match)[14] <- paste("cases_variant_count","_n=",(nrow(case_list)),sep = "")
colnames(var_gtex_match)[16] <- paste("controls_variant_count_variant_count","_n=",(nrow(control_list)),sep = "")


##Reduce table down depending on what the tissue/samples of interest are
tissue_equiv <- read.csv("gtex_tissues_with_single_word_equivalents.csv", header = T)

##Get lines with single word of interest in reference table and convert to gtex tissue labels
gtex_tis_1 <- which(tissue_equiv$single_word_equivalent_1 %in% t1)
gtex_tis_2 <- which(tissue_equiv$single_word_equivalent_1 %in% t2)
gtex_tis_3 <- which(tissue_equiv$single_word_equivalent_1 %in% t3)
gtex_tis_4 <- which(tissue_equiv$single_word_equivalent_1 %in% t4)
gtex_tis_all <- c(gtex_tis_1,gtex_tis_2,gtex_tis_3,gtex_tis_4)
gtex_tis_all <- tissue_equiv$gtex.tissue.filename[gtex_tis_all]

##Find which lines in the gtex match table contain the gtex tissue name/s of interest
tissue_lines_of_interest <- data.frame()

for (i in 1:ncol(tissue_cols)) {
  
  a <- which(tissue_cols[,i] %in% gtex_tis_all)
  tissue_lines_of_interest <- append(tissue_lines_of_interest, a)
  
}

tissue_lines_of_interest <- unlist(tissue_lines_of_interest)

##mkae new gtex match table only contining lines which contain the tissue of interest
if(t1 != "All")var_gtex_match <- var_gtex_match[tissue_lines_of_interest,]

##Hypothesis testing

##Fishers exact with fdr adjustment for each individual variant observed
rownum <- nrow(var_gtex_match)
BRIDGE_pvalue <- as.vector(c(), mode = "any")
for(i in 1:rownum){
  
  fishtable = matrix(c(var_gtex_match[i,14],(nrow(case_list) - var_gtex_match[i,14]),var_gtex_match[i,16],(nrow(control_list) - var_gtex_match[i,16])), 
                     nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue <- append(BRIDGE_pvalue, fisherout[[1]])
  
}

BRIDGE_qvalue <- p.adjust(BRIDGE_pvalue, method = "fdr", n = rownum)

##Collate results
var_gtex_match_with_pq <- cbind(BRIDGE_pvalue, BRIDGE_qvalue, var_gtex_match)

var_gtex_match_with_pq$BRIDGE_pvalue <- as.numeric(as.character(var_gtex_match_with_pq$BRIDGE_pvalue))
var_gtex_match_with_pq$BRIDGE_qvalue <- as.numeric(as.character(var_gtex_match_with_pq$BRIDGE_qvalue))
var_gtex_match_sig_index <- which(var_gtex_match_with_pq$BRIDGE_qvalue < 0.1)
var_gtex_match_sig <- var_gtex_match_with_pq[var_gtex_match_sig_index,]

options <- paste(args[1],args[2],args[3],args[4],sep = "_")

##Output results
write.csv(var_gtex_match_sig, paste("/home/jww39/non_coding/gtex/with_internal_af_filter/hets/", options,"_var_gtex_match.csv", sep = ""))
write.csv(var_gtex_match_with_pq, paste("/home/jww39/non_coding/gtex/with_internal_af_filter/hets/", options,"_var_gtex_match_all.csv", sep = ""))

##Per gene analysis

##Counts of variants per gene in cases and controls
genes <- unique(var_gtex_match$gtex_gene_id)
genes <- genes[!is.na(genes)]

variant_per_gene_count_cases <- data.frame()

for(i in 1:length(genes)){
  a <- sum(var_gtex_match[,14][which(var_gtex_match$gtex_gene_id %in% genes[i])]) 
  variant_per_gene_count_cases  <- append(variant_per_gene_count_cases,a)
}

variant_per_gene_count_cases <- t(data.frame(variant_per_gene_count_cases))
gene_names <- as.vector(genes)
rownames(variant_per_gene_count_cases) <- gene_names

variant_per_gene_count_controls <- data.frame()

for(i in 1:length(genes)){
  a <- sum(var_gtex_match[,16][which(var_gtex_match$gtex_gene_id %in% genes[i])])
  variant_per_gene_count_controls  <- append(variant_per_gene_count_controls,a)
}

variant_per_gene_count_controls <- t(data.frame(variant_per_gene_count_controls))
gene_names <- as.vector(genes)
rownames(variant_per_gene_count_controls) <- gene_names

##Fishers exact with fdr adjustment for counts of variants per gene
BRIDGE_pvalue_per_gene <- as.vector(c(), mode = "any")
cases_trials_per_gene <- as.vector(c(), mode = "any")
controls_trials_per_gene <- as.vector(c(), mode = "any")

for(i in 1:nrow(variant_per_gene_count_cases))
{
  fishtable = matrix(
    c(variant_per_gene_count_cases[i,1],
      sum(var_gtex_match$gtex_gene_id %in% genes[i]) * nrow(case_list) - variant_per_gene_count_cases[i,1],
      variant_per_gene_count_controls[i,1],
      sum(var_gtex_match$gtex_gene_id %in% genes[i]) * nrow(control_list) - variant_per_gene_count_controls[i,1]),
    
    nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue_per_gene <- append(BRIDGE_pvalue_per_gene, fisherout[[1]])
  cases_trials_per_gene <- append(cases_trials_per_gene, sum(var_gtex_match$gtex_gene_id %in% genes[i]) * nrow(case_list))
  controls_trials_per_gene <- append(controls_trials_per_gene, sum(var_gtex_match$gtex_gene_id %in% genes[i]) * nrow(control_list))
}

BRIDGE_qvalue_per_gene <- p.adjust(BRIDGE_pvalue_per_gene, method = "fdr", n = 83)

cases_per_gene_prop <- variant_per_gene_count_cases / cases_trials_per_gene
controls_per_gene_prop <- variant_per_gene_count_controls / controls_trials_per_gene

##Counts for individuals with variants per gene for cases and controls
trimmed_cases <- data.frame(var_gtex_match$cases)

indv_with_variant_per_gene_cases <- data.frame()

for(i in 1:length(genes)){
  a <- trimmed_cases[which(var_gtex_match$gtex_gene_id %in% genes[i]),]
  b <- a[which(a != "")]
  c <- paste(b, collapse = ";")
  indv_with_variant_per_gene_cases <- append(indv_with_variant_per_gene_cases, c)
}

indv_with_variant_per_gene_cases <- unlist(indv_with_variant_per_gene_cases)
indv_with_variant_per_gene_cases <- strsplit(indv_with_variant_per_gene_cases, ";", fixed = T)
unique_indv_with_variant_per_gene_cases <- lapply(indv_with_variant_per_gene_cases, unique)

indv_per_gene_count_cases <- data.frame()

for(i in 1:length(unique_indv_with_variant_per_gene_cases)){
  a <- length(which(unique_indv_with_variant_per_gene_cases[[i]] != ""))
  indv_per_gene_count_cases <- append(indv_per_gene_count_cases, a)
}

indv_per_gene_count_cases <- unlist(indv_per_gene_count_cases)

trimmed_controls <- data.frame(var_gtex_match$controls)

indv_with_variant_per_gene_controls <- data.frame()

for(i in 1:length(genes)){
  a <- trimmed_controls[which(var_gtex_match$gtex_gene_id %in% genes[i]),]
  b <- a[which(a != "")]
  c <- paste(b, collapse = ";")
  indv_with_variant_per_gene_controls <- append(indv_with_variant_per_gene_controls, c)
}

indv_with_variant_per_gene_controls <- unlist(indv_with_variant_per_gene_controls)
indv_with_variant_per_gene_controls <- strsplit(indv_with_variant_per_gene_controls, ";", fixed = T)
unique_indv_with_variant_per_gene_controls <- lapply(indv_with_variant_per_gene_controls, unique)

indv_per_gene_count_controls <- data.frame()

for(i in 1:length(unique_indv_with_variant_per_gene_controls)){
  a <- length(which(unique_indv_with_variant_per_gene_controls[[i]] != ""))
  indv_per_gene_count_controls <- append(indv_per_gene_count_controls, a)
}

indv_per_gene_count_controls <- unlist(indv_per_gene_count_controls)

##Fishers exact with fdr adjustment for counts of individuals with variants per gene
BRIDGE_pvalue_per_indv_per_gene <- as.vector(c(), mode = "any")
cases_indv_per_gene <- as.vector(c(), mode = "any")
controls_indv_per_gene <- as.vector(c(), mode = "any")

for(i in 1:length(indv_per_gene_count_cases))
{
  fishtable = matrix(
    c(indv_per_gene_count_cases[i],
      nrow(case_list) - indv_per_gene_count_cases[i],
      indv_per_gene_count_controls[i],
      nrow(control_list) - indv_per_gene_count_controls[i]), 
    nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue_per_indv_per_gene <- append(BRIDGE_pvalue_per_indv_per_gene, fisherout[[1]])
}

BRIDGE_qvalue_per_indv_per_gene <- p.adjust(BRIDGE_pvalue_per_indv_per_gene, method = "fdr", n = 83) 

cases_per_indv_per_gene_prop <- indv_per_gene_count_cases / nrow(case_list)
controls_per_indv_per_gene_prop <- indv_per_gene_count_controls / nrow(control_list)

##Collate results
variants_per_gene_table <- as.data.frame(cbind(rownames(variant_per_gene_count_cases),
                                               
                                               variant_per_gene_count_cases,
                                               cases_trials_per_gene,
                                               cases_per_gene_prop,
                                               variant_per_gene_count_controls,
                                               controls_trials_per_gene,
                                               controls_per_gene_prop,
                                               BRIDGE_pvalue_per_gene,
                                               BRIDGE_qvalue_per_gene,
                                               
                                               indv_per_gene_count_cases,
                                               cases_per_indv_per_gene_prop, 
                                               indv_per_gene_count_controls,
                                               controls_per_indv_per_gene_prop,
                                               BRIDGE_pvalue_per_indv_per_gene,
                                               BRIDGE_qvalue_per_indv_per_gene))

colnames(variants_per_gene_table) <- c("gene", 
                                       
                                       "variant_per_gene_count_cases", 
                                       "cases_trials_per_gene", 
                                       "cases_per_gene_prop",
                                       "variant_per_gene_count_controls",
                                       "controls_trials_per_gene",
                                       "controls_per_gene_prop",
                                       "BRIDGE_pvalue_per_gene",
                                       "BRIDGE_qvalue_per_gene",
                                       
                                       "indv_per_gene_count_cases",
                                       "cases_per_indv_per_gene_prop", 
                                       "indv_per_gene_count_controls",
                                       "controls_per_indv_per_gene_prop",
                                       "BRIDGE_pvalue_per_indv_per_gene",
                                       "BRIDGE_qvalue_per_indv_per_gene")


variants_per_gene_table$BRIDGE_qvalue_per_gene <- as.numeric(as.character(variants_per_gene_table$BRIDGE_qvalue_per_gene))
variants_per_gene_table$BRIDGE_pvalue_per_gene <- as.numeric(as.character(variants_per_gene_table$BRIDGE_pvalue_per_gene))
variants_per_gene_table$BRIDGE_qvalue_per_indv_per_gene <- as.numeric(as.character(variants_per_gene_table$BRIDGE_qvalue_per_indv_per_gene))
variants_per_gene_table$BRIDGE_pvalue_per_indv_per_gene <- as.numeric(as.character(variants_per_gene_table$BRIDGE_pvalue_per_indv_per_gene))

per_gene_sig_index <- which(variants_per_gene_table$BRIDGE_qvalue_per_gene < 0.1)
per_indv_per_gene_sig_index <- which(variants_per_gene_table$BRIDGE_qvalue_per_indv_per_gene < 0.1) 
sig_index <- c(per_gene_sig_index, per_indv_per_gene_sig_index)
sig_index <- unique(sig_index)

variants_per_gene_table_sig <- variants_per_gene_table[sig_index,]

##Output results
options <- paste(args[1],args[2],args[3],args[4],sep = "_")
write.csv(variants_per_gene_table_sig, paste("/home/jww39/non_coding/gtex/with_internal_af_filter/hets_homs/", options,"_variants_per_gene_table.csv", sep = "")) 
write.csv(variants_per_gene_table, paste("/home/jww39/non_coding/gtex/with_internal_af_filter/hets_homs/", options,"_variants_per_gene_table_all.csv", sep = "")) 








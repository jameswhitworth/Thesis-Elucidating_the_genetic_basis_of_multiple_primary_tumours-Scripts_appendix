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
ls /rds/project/who1000/rds-who1000-wgs10k/WGS10K/data/release/20170614-A/merged-vcf/no_hgmd/*.vcf.gz > 1_chromosome_files.txt
sed -r 's/.{7}$//' 1_chromosome_files.txt > 2_chromosome_files.txt
sed -r 's/.{90}//' 2_chromosome_files.txt > chromosome_files.txt
rm 1_chromosome_files.txt
rm 2_chromosome_files.txt

##Filter per chromosome VCF files using BED file corresponding to coordinates of interest. BED file containing all coordinates for gene lists
for i in `cat chromosome_files.txt`; do
bcftools view -R enhancers_promoters.bed -e 'FILTER!="PASS"' -o ${i}_filtered.vcf -O v /rds/project/who1000/rds-who1000-wgs10k/WGS10K/data/release/20170614-A/merged-vcf/no_hgmd/${i}.vcf.gz
done

##Merge filtered per chromosome VCF files and retain lines corresponding to samples of interest
ls *_filtered.vcf > filtered_samples_for_bgzip_and_tabix.txt

for i in `cat filtered_samples_for_bgzip_and_tabix.txt`; do
	vcf-sort -c ${i} > sorted_${i}
done

rm filtered_samples_for_bgzip_and_tabix.txt
ls sorted*.vcf | tr "\n" " " > files_to_merge.txt

bcftools concat -Oz `cat files_to_merge.txt` -o merged.vcf.gz

##Filter merged VCF based on quality parameters (insert missing genotypes where criteria not fulfilled)
bcftools filter -e 'FMT/DP<10 || FMT/GQ<30 || (FMT/AD[1]/(FMT/AD[0] + FMT/AD[1]) <0.3) || (FMT/AD[1]/(FMT/AD[0] + FMT/AD[1] + FMT/AD[2]) <0.3) || (FMT/AD[2]/(FMT/AD[0] + FMT/AD[1] + FMT/AD[2]) <0.3)' -S . -o enhancers_promoters_merged_filtered.vcf.gz -O z merged.vcf.gz

##Transferred to local server

###################################################
##Further data addition, annotation and filtering##
###################################################

##Execute VEP filter script
/home/jww39/ensembl-vep/filter_vep -i enhancers_promoters_merged_filtered.vcf -o VEP_filtered_enhancers_promoters_merged.vcf -filter "(
1000G_EUR_AF < 0.01 or not 1000G_EUR_AF) and (UK10KWGS_AF < 0.01 or not UK10KWGS_AF) and (WGS10K_AF < 0.05 or not WGS10K_AF)" --only_matched

#Prepare resulting VCF for reading into R
sed -i 's/#CHROM/CHROM/g' VEP_filtered_enhancers_promoters_merged.vcf
sed '/^#/ d' VEP_filtered_enhancers_promoters_merged.vcf > VEP_filtered_enhancers_promoters_merged_without_header.vcf
mv VEP_filtered_enhancers_promoters_merged.vcf VEP_filtered_enhancers_promoters_merged_with_header.vcf
mv VEP_filtered_enhancers_promoters_merged_without_header.vcf  VEP_filtered_enhancers_promoters_merged.vcf

###########################################################
##R script to identify samples with each variant in table##
###########################################################
args = commandArgs(trailingOnly=TRUE)

MPT_BRIDGE <- read.delim("VEP_filtered_enhancers_promoters_merged.vcf", header = T) #name of input vcf

##Identifying samples (cases and controls) containing each variant
MPT_BRIDGE_gt_only <- MPT_BRIDGE[,10:ncol(MPT_BRIDGE)] 

samples_hets_raw <- data.frame()

for (i in 1:nrow(MPT_BRIDGE_gt_only )) { 
  
  a <- unlist(MPT_BRIDGE_gt_only[i,])
  b <- grep("^0\\/1", a, perl = TRUE) ##Part to change for hets/homs##
  c <- as.vector(names(a[b]))
  d <- paste(c, collapse = ';')
  samples_hets_raw <- append(samples_hets_raw, d)  
} 
samples_hets <- unlist(samples_hets_raw)

samples_homs_raw <- data.frame()

for (i in 1:nrow(MPT_BRIDGE_gt_only )) {
  
  a <- unlist(MPT_BRIDGE_gt_only[i,])
  b <- grep("^1\\/1", a, perl = TRUE) ##Part to change for hets/homs##
  c <- as.vector(names(a[b]))
  d <- paste(c, collapse = ';')
  samples_homs_raw <- append(samples_homs_raw, d)
} 
samples_homs <- unlist(samples_homs_raw)

##Output results
write(samples_hets, "samples_hets.txt")
write(samples_homs, "samples_homs.txt")

#########################################################################
##R script to count variants and compare frequency in cases vs controls##
#########################################################################

##################################################################################################################################
##Following commands used to launch script according to phenotypic subgroup:
#Rscript non_coding_SNV_lookup_subsequent.R ALL XXX XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE ACC XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE Aerodigestivetract XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE Bladder XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE Bonebenign XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE Breast XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE Cervix XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE CNS XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE CNShaemangioblastoma XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE CNSmeningioma XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE CNSNervesheath XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE Colorectal XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE Colorectalpolyps XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE Connectivetissuesofttissuesarcoma XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE Endometrium XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE GINET XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE GIST XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE Haemlymphoid XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE Haemmyeloid XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE Kidney XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE Kidneyoncocytoma XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE Lung XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE Lungcarcinoid XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE Melanoma XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE NMSC XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE Oesophagus XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE Ovary XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE Ovarysexcord-gonadalstromal XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE Pancreas XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE Paraganglioma XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE Parathyroid XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE Phaeochromocytoma XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE Pituitary XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE PNET XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE PNSNervesheathbenign XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE Prostate XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE Salivarygland XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE Smallbowel XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE Testicular XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE Thyroid XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R SINGLE Uvealmelanoma XXX XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R BOTH Breast Colorectal XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R BOTH Breast Breast XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R BOTH Breast NMSC XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R BOTH Breast Endometrium XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R BOTH Breast Ovary XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R BOTH Breast Haemlymphoid XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R BOTH Breast Melanoma XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R BOTH Breast Thyroid XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R BOTH Endometrium Ovary XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R BOTH Breast Kidney XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R BOTH Kidney Kidney XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R BOTH Colorectal NMSC XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R BOTH NMSC NMSC XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R BOTH Breast Lung XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R BOTH Colorectal Colorectal XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R BOTH Breast Connectivetissuesofttissuesarcoma XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R BOTH Colorectal Endometrium XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R BOTH Kidney Thyroid XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R BOTH Breast CNSmeningioma XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R BOTH Paraganglioma Paraganglioma XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R BOTH Breast Cervix XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R BOTH Colorectal Prostate XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R BOTH Colorectal Thyroid XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R BOTH Kidney Lung XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM2 Breast Ovary XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM2 Colorectal Endometrium XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM2 Thyroidmedullary Phaeochromocytoma XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM2 NMSC Melanoma XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM2 Colorectal Gastric XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM2 Parathyroid Bonebenign XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM2 Breast Gastric XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM2 Thyroid Pituitary XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM2 NMSC Haemlymphoid XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM2 Haemlymphoid Haemmyeloid XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM2 Connectivetissuesofttissuesarcoma Bladder XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM2 Breast Pancreas XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM2 NMSC Bonesarcoma XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM2 NMSC Haemlymphoid XXX XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM3 Haemmyeloid Aerodigestivetract Anus XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM3 Melanoma Pancreas CNS XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM3 Kidney Kidneyangiomyolipoma CNS XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM3 CNSmeningioma CNS CNSNervesheath XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM3 Haemlymphoid CNS Connectivetissuesofttissuesarcoma XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM3 Breast Thyroid Endometrium XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM3 Phaeochromocytoma Paraganglioma GIST XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM3 Wilms Connectivetissuesofttissuesarcoma Haemmyeloid XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM3 Cardiacmyxoma Thyroid Ovarysexcord-gonadalstromal XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM3 CNS PNSNervesheath PNSNervesheathbenign XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM4 Kidney Phaeochromocytoma Paraganglioma CNShaemangioblastoma XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM4 CNS CNShaemangioblastoma CNSmeningioma CNSNervesheath XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM4 Kidney Uterineleiomyoma Uterinesarcoma Cutaneousleiomyoma XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM4 Kidney Adrenaloncocytoma Kidneyoncocytoma Fibrofolliculoma XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM4 Haemmyeloid Aerodigestivetract Anus Melanoma XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM4 Breast Aerodigestivetract Lung Ovary XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM4 Colorectal Breast Gastric Ovarysexcord-gonadalstromal XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM4 Colorectal Endometrium Ovary Sebaceous XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM5 Haemmyeloid Aerodigestivetract Oesophagus Cervix Penis XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM5 Connectivetissuesofttissuesarcoma Bonesarcoma GIST Skinsarcoma Uterinesarcoma XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM6 Uvealmelanoma Kidney Melanoma Lung Mesothelioma CNSmeningioma XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM6 GINET Lungcarcinoid Ovaryneuroendocrine Paraganglioma Phaeochromocytoma PNET XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM7 Retinoblastoma Connectivetissuesofttissuesarcoma Bonesarcoma GIST Skinsarcoma Uterinesarcoma Melanoma XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM7 Connectivetissuesofttissuesarcoma Bonesarcoma GIST Skinsarcoma Uterinesarcoma Melanoma Thyroid XXX
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM8 Pituitary Parathyroid ACC GINET Lungcarcinoid Ovaryneuroendocrine Paraganglioma Phaeochromocytoma
#Rscript non_coding_SNV_lookup_subsequent.R 1FROM8 Breast ACC CNS Connectivetissuesofttissuesarcoma Bonesarcoma GIST Skinsarcoma Uterinesarcoma
#Rscript non_coding_SNV_lookup_subsequent.R 2FROM3 Breast Thyroid Endometrium XXX XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 2FROM4 Breast Aerodigestivetract Lung Ovary XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 2FROM4 Colorectal Breast Gastric Ovarysexcord-gonadalstromal XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 2FROM4 Colorectal Endometrium Ovary Sebaceous XXX XXX XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 2FROM6 Uvealmelanoma Kidney Melanoma Lung Mesothelioma CNSmeningioma XXX XXX
#Rscript non_coding_SNV_lookup_subsequent.R 2FROM7 Connectivetissuesofttissuesarcoma Bonesarcoma GIST Skinsarcoma Uterinesarcoma Melanoma Thyroid XXX
#Rscript non_coding_SNV_lookup_subsequent.R 2FROM8 Breast ACC CNS Connectivetissuesofttissuesarcoma Bonesarcoma GIST Skinsarcoma Uterinesarcoma
############################################################################################################################################################

args = commandArgs(trailingOnly=TRUE)

##Read in necessary files
MPT_BRIDGE <- read.delim("VEP_filtered_enhancers_promoters_merged.vcf", header = T) #name of input vcf
elements <- read.csv("enhancers_promoters_elements.csv", header = T)
BRIDGE_list <- read.delim("BRIDGE_euro_controls.txt", header = F)

##Identify samples with phenotype of interest
MPT_table <- read.csv("anonymised_MPT_datasheet_22_11_17.csv") ##Datasheet of research participants and phenotype


rownames(MPT_table) <- MPT_table$WGS.ID.1
MPT_table$Age.1 <- as.numeric(as.character(MPT_table$Age.1)) 
MPT_table$Age.2 <- as.numeric(as.character(MPT_table$Age.2))
MPT_table$Age.3 <- as.numeric(as.character(MPT_table$Age.3))
MPT_table$Age.4 <- as.numeric(as.character(MPT_table$Age.4))
MPT_table$Age.5 <- as.numeric(as.character(MPT_table$Age.5))
MPT_table$Age.6 <- as.numeric(as.character(MPT_table$Age.6))
MPT_table$Age.7 <- as.numeric(as.character(MPT_table$Age.7))

t1 <- args[2]
t2 <- args[3]
t3 <- args[4]
t4 <- args[5]
t5 <- args[6]
t6 <- args[7]
t7 <- args[8]
t8 <- args[9]

##Single
if(args[1] == "SINGLE")ting <- rownames(MPT_table[
  
  ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
   | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
   | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
   | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
   | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
   | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
   | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
  
  &  
    
    (MPT_table$proband == "YES")
  
  &
    
    (MPT_table$exclusion != "EXCLUDE")
  
  ,]) 

##2 from 2
if(args[1] == "BOTH")ting <- rownames(MPT_table[
  
  (
    
    ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
     | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
     | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
     | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
     | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
     | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
     | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
    
    &
      
      ((MPT_table$Single.word.1 == t2 & MPT_table$Age.1 < 70) 
       | (MPT_table$Single.word.2 == t2 & MPT_table$Age.2 < 70) 
       | (MPT_table$Single.word.3 == t2 & MPT_table$Age.3 < 70)
       | (MPT_table$Single.word.4 == t2 & MPT_table$Age.4 < 70)
       | (MPT_table$Single.word.5 == t2 & MPT_table$Age.5 < 70)
       | (MPT_table$Single.word.6 == t2 & MPT_table$Age.6 < 70)
       | (MPT_table$Single.word.7 == t2 & MPT_table$Age.7 < 70))
    
  )
  
  & 
    
    (MPT_table$proband == "YES")
  
  &
    
    (MPT_table$exclusion != "EXCLUDE")
  
  ,]) 

##1 from 2
if(args[1] == "1FROM2")ting <- rownames(MPT_table[
  
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
    
  )
  
  & 
    
    (MPT_table$proband == "YES")
  
  &
    
    (MPT_table$exclusion != "EXCLUDE")
  
  ,])
  
###1 from 3
if(args[1] == "1FROM3")ting <- rownames(MPT_table[
  
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
    
  )
  
  & 
    
    (MPT_table$proband == "YES")
  
  &
    
    (MPT_table$exclusion != "EXCLUDE")
  
  ,])

###1 from 4
if(args[1] == "1FROM4")ting <- rownames(MPT_table[
  
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

##1 from 5
if(args[1] == "1FROM5")ting <- rownames(MPT_table[
  
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
    
    |
      
      ((MPT_table$Single.word.1 == t5 & MPT_table$Age.1 < 70) 
       | (MPT_table$Single.word.2 == t5 & MPT_table$Age.2 < 70) 
       | (MPT_table$Single.word.3 == t5 & MPT_table$Age.3 < 70)
       | (MPT_table$Single.word.4 == t5 & MPT_table$Age.4 < 70)
       | (MPT_table$Single.word.5 == t5 & MPT_table$Age.5 < 70)
       | (MPT_table$Single.word.6 == t5 & MPT_table$Age.6 < 70)
       | (MPT_table$Single.word.7 == t5 & MPT_table$Age.7 < 70))
    
  )
  
  & 
    
    (MPT_table$proband == "YES")
  
  &
    
    (MPT_table$exclusion != "EXCLUDE")
  
  ,])

##1 from 6
if(args[1] == "1FROM6")ting <- rownames(MPT_table[
  
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
    
    |
      
      ((MPT_table$Single.word.1 == t5 & MPT_table$Age.1 < 70) 
       | (MPT_table$Single.word.2 == t5 & MPT_table$Age.2 < 70) 
       | (MPT_table$Single.word.3 == t5 & MPT_table$Age.3 < 70)
       | (MPT_table$Single.word.4 == t5 & MPT_table$Age.4 < 70)
       | (MPT_table$Single.word.5 == t5 & MPT_table$Age.5 < 70)
       | (MPT_table$Single.word.6 == t5 & MPT_table$Age.6 < 70)
       | (MPT_table$Single.word.7 == t5 & MPT_table$Age.7 < 70))
    
    |
      
      ((MPT_table$Single.word.1 == t6 & MPT_table$Age.1 < 70) 
       | (MPT_table$Single.word.2 == t6 & MPT_table$Age.2 < 70) 
       | (MPT_table$Single.word.3 == t6 & MPT_table$Age.3 < 70)
       | (MPT_table$Single.word.4 == t6 & MPT_table$Age.4 < 70)
       | (MPT_table$Single.word.5 == t6 & MPT_table$Age.5 < 70)
       | (MPT_table$Single.word.6 == t6 & MPT_table$Age.6 < 70)
       | (MPT_table$Single.word.7 == t6 & MPT_table$Age.7 < 70))
    
  )
  
  & 
    
    (MPT_table$proband == "YES")
  
  &
    
    (MPT_table$exclusion != "EXCLUDE")
  
  ,])

###1 from 7
if(args[1] == "1FROM7")ting <- rownames(MPT_table[
  
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
    
    |
      
      ((MPT_table$Single.word.1 == t5 & MPT_table$Age.1 < 70) 
       | (MPT_table$Single.word.2 == t5 & MPT_table$Age.2 < 70) 
       | (MPT_table$Single.word.3 == t5 & MPT_table$Age.3 < 70)
       | (MPT_table$Single.word.4 == t5 & MPT_table$Age.4 < 70)
       | (MPT_table$Single.word.5 == t5 & MPT_table$Age.5 < 70)
       | (MPT_table$Single.word.6 == t5 & MPT_table$Age.6 < 70)
       | (MPT_table$Single.word.7 == t5 & MPT_table$Age.7 < 70))
    
    |
      
      ((MPT_table$Single.word.1 == t6 & MPT_table$Age.1 < 70) 
       | (MPT_table$Single.word.2 == t6 & MPT_table$Age.2 < 70) 
       | (MPT_table$Single.word.3 == t6 & MPT_table$Age.3 < 70)
       | (MPT_table$Single.word.4 == t6 & MPT_table$Age.4 < 70)
       | (MPT_table$Single.word.5 == t6 & MPT_table$Age.5 < 70)
       | (MPT_table$Single.word.6 == t6 & MPT_table$Age.6 < 70)
       | (MPT_table$Single.word.7 == t6 & MPT_table$Age.7 < 70))
    
    |
      
      ((MPT_table$Single.word.1 == t7 & MPT_table$Age.1 < 70) 
       | (MPT_table$Single.word.2 == t7 & MPT_table$Age.2 < 70) 
       | (MPT_table$Single.word.3 == t7 & MPT_table$Age.3 < 70)
       | (MPT_table$Single.word.4 == t7 & MPT_table$Age.4 < 70)
       | (MPT_table$Single.word.5 == t7 & MPT_table$Age.5 < 70)
       | (MPT_table$Single.word.6 == t7 & MPT_table$Age.6 < 70)
       | (MPT_table$Single.word.7 == t7 & MPT_table$Age.7 < 70))
    
  )
  
  & 
    
    (MPT_table$proband == "YES")
  
  &
    
    (MPT_table$exclusion != "EXCLUDE")
  
  ,])

##1 from 8
if(args[1] == "1FROM8")ting <- rownames(MPT_table[
  
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
    
    |
      
      ((MPT_table$Single.word.1 == t5 & MPT_table$Age.1 < 70) 
       | (MPT_table$Single.word.2 == t5 & MPT_table$Age.2 < 70) 
       | (MPT_table$Single.word.3 == t5 & MPT_table$Age.3 < 70)
       | (MPT_table$Single.word.4 == t5 & MPT_table$Age.4 < 70)
       | (MPT_table$Single.word.5 == t5 & MPT_table$Age.5 < 70)
       | (MPT_table$Single.word.6 == t5 & MPT_table$Age.6 < 70)
       | (MPT_table$Single.word.7 == t5 & MPT_table$Age.7 < 70))
    
    |
      
      ((MPT_table$Single.word.1 == t6 & MPT_table$Age.1 < 70) 
       | (MPT_table$Single.word.2 == t6 & MPT_table$Age.2 < 70) 
       | (MPT_table$Single.word.3 == t6 & MPT_table$Age.3 < 70)
       | (MPT_table$Single.word.4 == t6 & MPT_table$Age.4 < 70)
       | (MPT_table$Single.word.5 == t6 & MPT_table$Age.5 < 70)
       | (MPT_table$Single.word.6 == t6 & MPT_table$Age.6 < 70)
       | (MPT_table$Single.word.7 == t6 & MPT_table$Age.7 < 70))
    
    |
      
      ((MPT_table$Single.word.1 == t7 & MPT_table$Age.1 < 70) 
       | (MPT_table$Single.word.2 == t7 & MPT_table$Age.2 < 70) 
       | (MPT_table$Single.word.3 == t7 & MPT_table$Age.3 < 70)
       | (MPT_table$Single.word.4 == t7 & MPT_table$Age.4 < 70)
       | (MPT_table$Single.word.5 == t7 & MPT_table$Age.5 < 70)
       | (MPT_table$Single.word.6 == t7 & MPT_table$Age.6 < 70)
       | (MPT_table$Single.word.7 == t7 & MPT_table$Age.7 < 70))
    
    |
      
      ((MPT_table$Single.word.1 == t8 & MPT_table$Age.1 < 70) 
       | (MPT_table$Single.word.2 == t8 & MPT_table$Age.2 < 70) 
       | (MPT_table$Single.word.3 == t8 & MPT_table$Age.3 < 70)
       | (MPT_table$Single.word.4 == t8 & MPT_table$Age.4 < 70)
       | (MPT_table$Single.word.5 == t8 & MPT_table$Age.5 < 70)
       | (MPT_table$Single.word.6 == t8 & MPT_table$Age.6 < 70)
       | (MPT_table$Single.word.7 == t8 & MPT_table$Age.7 < 70))
    
  )
  
  & 
    
    (MPT_table$proband == "YES")
  
  &
    
    (MPT_table$exclusion != "EXCLUDE")
  
  ,])

##2 from 3
if(args[1] == "2FROM3")ting <- rownames(MPT_table[
  
  (
    
    (
      ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
       | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
       | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
       | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
       | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
       | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
       | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
      
      &
        
        ((MPT_table$Single.word.1 == t2 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t2 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t2 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t2 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t2 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t2 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t2 & MPT_table$Age.7 < 70))
    )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t3 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t3 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t3 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t3 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t3 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t3 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t3 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t2 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t2 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t2 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t2 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t2 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t2 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t2 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t3 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t3 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t3 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t3 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t3 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t3 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t3 & MPT_table$Age.7 < 70))
      )
    
  )
  
  & 
    
    (MPT_table$proband == "YES")
  
  &
    
    (MPT_table$exclusion != "EXCLUDE")
  
  ,])

##2 from 4
if(args[1] == "2FROM4")ting <- rownames(MPT_table[
  
  (
    
    (
      ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
       | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
       | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
       | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
       | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
       | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
       | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
      
      &
        
        ((MPT_table$Single.word.1 == t2 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t2 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t2 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t2 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t2 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t2 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t2 & MPT_table$Age.7 < 70))
    )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t3 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t3 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t3 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t3 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t3 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t3 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t3 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t4 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t4 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t4 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t4 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t4 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t4 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t4 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t2 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t2 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t2 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t2 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t2 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t2 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t2 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t3 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t3 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t3 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t3 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t3 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t3 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t3 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t3 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t3 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t3 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t3 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t3 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t3 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t3 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t4 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t4 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t4 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t4 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t4 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t4 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t4 & MPT_table$Age.7 < 70))
      )
    
  )
  
  & 
    
    (MPT_table$proband == "YES")
  
  &
    
    (MPT_table$exclusion != "EXCLUDE")
  
  ,])

##2 from 5
if(args[1] == "2FROM5")ting <- rownames(MPT_table[
  
  (
    
    (
      ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
       | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
       | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
       | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
       | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
       | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
       | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
      
      &
        
        ((MPT_table$Single.word.1 == t2 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t2 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t2 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t2 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t2 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t2 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t2 & MPT_table$Age.7 < 70))
    )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t3 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t3 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t3 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t3 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t3 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t3 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t3 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t4 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t4 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t4 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t4 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t4 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t4 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t4 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t5 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t5 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t5 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t5 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t5 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t5 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t5 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t2 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t2 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t2 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t2 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t2 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t2 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t2 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t3 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t3 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t3 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t3 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t3 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t3 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t3 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t2 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t2 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t2 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t2 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t2 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t2 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t2 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t4 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t4 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t4 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t4 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t4 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t4 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t4 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t2 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t2 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t2 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t2 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t2 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t2 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t2 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t5 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t5 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t5 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t5 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t5 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t5 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t5 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t3 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t3 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t3 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t3 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t3 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t3 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t3 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t4 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t4 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t4 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t4 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t4 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t4 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t4 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t3 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t3 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t3 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t3 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t3 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t3 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t3 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t5 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t5 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t5 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t5 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t5 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t5 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t5 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t4 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t4 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t4 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t4 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t4 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t4 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t4 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t5 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t5 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t5 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t5 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t5 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t5 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t5 & MPT_table$Age.7 < 70))
      )
    
  )
  
  & 
    
    (MPT_table$proband == "YES")
  
  &
    
    (MPT_table$exclusion != "EXCLUDE")
  
  ,])

##2 from 6
if(args[1] == "2FROM6")ting <- rownames(MPT_table[
  
  (
    
    (
      ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
       | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
       | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
       | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
       | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
       | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
       | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
      
      &
        
        ((MPT_table$Single.word.1 == t2 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t2 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t2 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t2 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t2 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t2 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t2 & MPT_table$Age.7 < 70))
    )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t3 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t3 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t3 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t3 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t3 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t3 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t3 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t4 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t4 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t4 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t4 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t4 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t4 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t4 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t5 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t5 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t5 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t5 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t5 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t5 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t5 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t6 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t6 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t6 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t6 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t6 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t6 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t6 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t2 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t2 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t2 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t2 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t2 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t2 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t2 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t3 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t3 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t3 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t3 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t3 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t3 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t3 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t2 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t2 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t2 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t2 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t2 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t2 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t2 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t4 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t4 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t4 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t4 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t4 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t4 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t4 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t2 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t2 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t2 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t2 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t2 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t2 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t2 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t5 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t5 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t5 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t5 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t5 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t5 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t5 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t2 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t2 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t2 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t2 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t2 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t2 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t2 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t6 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t6 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t6 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t6 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t6 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t6 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t6 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t3 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t3 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t3 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t3 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t3 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t3 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t3 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t4 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t4 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t4 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t4 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t4 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t4 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t4 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t3 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t3 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t3 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t3 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t3 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t3 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t3 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t5 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t5 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t5 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t5 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t5 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t5 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t5 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t3 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t3 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t3 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t3 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t3 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t3 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t3 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t6 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t6 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t6 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t6 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t6 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t6 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t6 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t4 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t4 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t4 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t4 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t4 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t4 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t4 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t5 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t5 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t5 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t5 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t5 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t5 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t5 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t4 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t4 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t4 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t4 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t4 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t4 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t4 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t6 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t6 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t6 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t6 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t6 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t6 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t6 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t5 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t5 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t5 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t5 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t5 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t5 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t5 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t6 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t6 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t6 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t6 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t6 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t6 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t6 & MPT_table$Age.7 < 70))
      )
    
  )
  
  & 
    
    (MPT_table$proband == "YES")
  
  &
    
    (MPT_table$exclusion != "EXCLUDE")
  
  ,])

##2 from 7
if(args[1] == "2FROM7")ting <- rownames(MPT_table[
  
  (
    
    (
      ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
       | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
       | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
       | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
       | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
       | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
       | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
      
      &
        
        ((MPT_table$Single.word.1 == t2 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t2 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t2 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t2 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t2 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t2 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t2 & MPT_table$Age.7 < 70))
    )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t3 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t3 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t3 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t3 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t3 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t3 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t3 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t4 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t4 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t4 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t4 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t4 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t4 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t4 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t5 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t5 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t5 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t5 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t5 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t5 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t5 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t6 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t6 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t6 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t6 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t6 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t6 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t6 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t7 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t7 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t7 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t7 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t7 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t7 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t7 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t2 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t2 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t2 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t2 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t2 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t2 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t2 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t3 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t3 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t3 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t3 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t3 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t3 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t3 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t2 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t2 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t2 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t2 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t2 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t2 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t2 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t4 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t4 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t4 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t4 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t4 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t4 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t4 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t2 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t2 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t2 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t2 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t2 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t2 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t2 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t5 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t5 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t5 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t5 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t5 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t5 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t5 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t2 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t2 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t2 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t2 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t2 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t2 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t2 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t6 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t6 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t6 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t6 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t6 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t6 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t6 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t2 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t2 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t2 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t2 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t2 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t2 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t2 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t7 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t7 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t7 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t7 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t7 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t7 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t7 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t3 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t3 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t3 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t3 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t3 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t3 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t3 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t4 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t4 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t4 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t4 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t4 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t4 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t4 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t3 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t3 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t3 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t3 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t3 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t3 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t3 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t5 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t5 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t5 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t5 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t5 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t5 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t5 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t3 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t3 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t3 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t3 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t3 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t3 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t3 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t6 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t6 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t6 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t6 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t6 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t6 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t6 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t3 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t3 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t3 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t3 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t3 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t3 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t3 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t7 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t7 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t7 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t7 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t7 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t7 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t7 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t4 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t4 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t4 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t4 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t4 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t4 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t4 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t5 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t5 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t5 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t5 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t5 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t5 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t5 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t4 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t4 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t4 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t4 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t4 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t4 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t4 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t6 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t6 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t6 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t6 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t6 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t6 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t6 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t4 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t4 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t4 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t4 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t4 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t4 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t4 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t7 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t7 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t7 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t7 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t7 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t7 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t7 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t5 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t5 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t5 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t5 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t5 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t5 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t5 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t6 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t6 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t6 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t6 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t6 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t6 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t6 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t5 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t5 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t5 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t5 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t5 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t5 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t5 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t7 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t7 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t7 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t7 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t7 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t7 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t7 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t6 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t6 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t6 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t6 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t6 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t6 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t6 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t7 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t7 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t7 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t7 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t7 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t7 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t7 & MPT_table$Age.7 < 70))
      )
    
  )
  
  & 
    
    (MPT_table$proband == "YES")
  
  &
    
    (MPT_table$exclusion != "EXCLUDE")
  
  ,])

##2 from 8
if(args[1] == "2FROM8")ting <- rownames(MPT_table[
  
  (
    
    (
      ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
       | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
       | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
       | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
       | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
       | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
       | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
      
      &
        
        ((MPT_table$Single.word.1 == t2 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t2 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t2 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t2 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t2 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t2 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t2 & MPT_table$Age.7 < 70))
    )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t3 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t3 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t3 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t3 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t3 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t3 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t3 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t4 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t4 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t4 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t4 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t4 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t4 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t4 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t5 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t5 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t5 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t5 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t5 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t5 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t5 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t6 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t6 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t6 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t6 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t6 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t6 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t6 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t7 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t7 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t7 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t7 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t7 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t7 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t7 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t1 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t1 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t1 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t1 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t1 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t1 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t1 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t8 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t8 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t8 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t8 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t8 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t8 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t8 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t2 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t2 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t2 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t2 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t2 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t2 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t2 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t3 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t3 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t3 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t3 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t3 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t3 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t3 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t2 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t2 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t2 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t2 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t2 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t2 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t2 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t4 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t4 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t4 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t4 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t4 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t4 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t4 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t2 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t2 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t2 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t2 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t2 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t2 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t2 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t5 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t5 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t5 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t5 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t5 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t5 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t5 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t2 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t2 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t2 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t2 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t2 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t2 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t2 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t6 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t6 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t6 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t6 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t6 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t6 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t6 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t2 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t2 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t2 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t2 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t2 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t2 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t2 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t7 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t7 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t7 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t7 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t7 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t7 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t7 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t2 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t2 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t2 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t2 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t2 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t2 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t2 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t8 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t8 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t8 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t8 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t8 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t8 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t8 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t3 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t3 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t3 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t3 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t3 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t3 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t3 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t4 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t4 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t4 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t4 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t4 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t4 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t4 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t3 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t3 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t3 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t3 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t3 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t3 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t3 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t5 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t5 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t5 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t5 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t5 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t5 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t5 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t3 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t3 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t3 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t3 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t3 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t3 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t3 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t6 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t6 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t6 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t6 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t6 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t6 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t6 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t3 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t3 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t3 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t3 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t3 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t3 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t3 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t7 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t7 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t7 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t7 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t7 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t7 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t7 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t3 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t3 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t3 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t3 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t3 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t3 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t3 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t8 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t8 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t8 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t8 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t8 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t8 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t8 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t4 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t4 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t4 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t4 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t4 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t4 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t4 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t5 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t5 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t5 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t5 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t5 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t5 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t5 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t4 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t4 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t4 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t4 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t4 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t4 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t4 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t6 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t6 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t6 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t6 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t6 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t6 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t6 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t4 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t4 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t4 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t4 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t4 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t4 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t4 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t7 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t7 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t7 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t7 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t7 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t7 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t7 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t4 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t4 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t4 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t4 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t4 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t4 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t4 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t8 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t8 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t8 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t8 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t8 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t8 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t8 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t5 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t5 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t5 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t5 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t5 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t5 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t5 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t6 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t6 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t6 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t6 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t6 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t6 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t6 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t5 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t5 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t5 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t5 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t5 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t5 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t5 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t7 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t7 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t7 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t7 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t7 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t7 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t7 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t5 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t5 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t5 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t5 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t5 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t5 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t5 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t8 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t8 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t8 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t8 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t8 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t8 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t8 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t6 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t6 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t6 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t6 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t6 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t6 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t6 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t7 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t7 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t7 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t7 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t7 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t7 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t7 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t6 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t6 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t6 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t6 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t6 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t6 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t6 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t8 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t8 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t8 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t8 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t8 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t8 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t8 & MPT_table$Age.7 < 70))
      )
    
    |
      
      (
        ((MPT_table$Single.word.1 == t7 & MPT_table$Age.1 < 70) 
         | (MPT_table$Single.word.2 == t7 & MPT_table$Age.2 < 70) 
         | (MPT_table$Single.word.3 == t7 & MPT_table$Age.3 < 70)
         | (MPT_table$Single.word.4 == t7 & MPT_table$Age.4 < 70)
         | (MPT_table$Single.word.5 == t7 & MPT_table$Age.5 < 70)
         | (MPT_table$Single.word.6 == t7 & MPT_table$Age.6 < 70)
         | (MPT_table$Single.word.7 == t7 & MPT_table$Age.7 < 70))
        
        &
          
          ((MPT_table$Single.word.1 == t8 & MPT_table$Age.1 < 70) 
           | (MPT_table$Single.word.2 == t8 & MPT_table$Age.2 < 70) 
           | (MPT_table$Single.word.3 == t8 & MPT_table$Age.3 < 70)
           | (MPT_table$Single.word.4 == t8 & MPT_table$Age.4 < 70)
           | (MPT_table$Single.word.5 == t8 & MPT_table$Age.5 < 70)
           | (MPT_table$Single.word.6 == t8 & MPT_table$Age.6 < 70)
           | (MPT_table$Single.word.7 == t8 & MPT_table$Age.7 < 70))
      )
    
  )
  
  & 
    
    (MPT_table$proband == "YES")
  
  &
    
    (MPT_table$exclusion != "EXCLUDE")
  
  ,])
  
##Preparing a list of cases
if(args[1] != "ALL")indv <- grep("R0*", ting, value = T)
if(args[1] != "ALL")indv <- gsub("_A", "", indv)
if(args[1] != "ALL")indv <- data.frame(indv)
if(args[1] != "ALL")case_list <- indv
if(args[1] == "ALL")case_list <- read.delim("MPT_probands_euro.txt", header = F)

##Read in lists of samples with variants generated from previous R script
samples_hets <- readLines("samples_hets.txt", n= -1L)
samples_homs <- readLines("samples_homs.txt", n= -1L)

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

##Produce columns containing variant counts 
cases_variant_count_hets_homs <- cases_variant_count_hets + cases_variant_count_homs
BRIDGES_variant_count_hets_homs <- BRIDGES_variant_count_hets + BRIDGES_variant_count_homs

cases_variant_count_hets_homs_agg <- cases_variant_count_hets + (2 * cases_variant_count_homs)
BRIDGES_variant_count_hets_homs_agg <- BRIDGES_variant_count_hets + (2 * BRIDGES_variant_count_homs)


##Match variants to element that they are within
MPT_BRIDGE_chrom_pos_only <- data.frame(MPT_BRIDGE$CHROM, MPT_BRIDGE$POS, paste(MPT_BRIDGE$CHROM, MPT_BRIDGE$POS, sep = "_"))
colnames(MPT_BRIDGE_chrom_pos_only) <- c("CHROM", "POS", "CHROM_POS")

library(sqldf)

variant_element_match <- sqldf("select * 
                               from MPT_BRIDGE_chrom_pos_only f1 
                               left join elements f2 on (f1.POS >= f2.start and f1.POS <= f2.end and f1.CHROM==f2.chrom) ")

##Take the first match (only interested in identifier from the interim dataframe and produce column to put into MPT_BRIDGE_counts)
MPT_BRIDGE_variant_element_match_indices <- match(MPT_BRIDGE_chrom_pos_only$CHROM_POS, variant_element_match$CHROM_POS)
variant_element_hits_IDs <- variant_element_match$identifier[MPT_BRIDGE_variant_element_match_indices]

##Produce a table with variant information and counts
MPT_BRIDGE_counts <- cbind(
  MPT_BRIDGE[,1:9],
  cases_variant_count_hets, #10
  cases_variant_count_homs, #11
  cases_variant_count_hets_homs, #12
  cases_variant_count_hets_homs_agg, #13
  BRIDGES_variant_count_hets, #14
  BRIDGES_variant_count_homs, #15
  BRIDGES_variant_count_hets_homs, #16
  BRIDGES_variant_count_hets_homs_agg, #17
  cases_hets,#18
  BRIDGES_hets,#19
  cases_homs,#20
  BRIDGES_homs,#21
  variant_element_hits_IDs, ##22
  MPT_BRIDGE[,10:ncol(MPT_BRIDGE)] #23 onwards
)

############################
##Case control comparisons##
############################

##Fishers exact with fdr adjustment for each individual variant observed (hets cases vs BRIDGE)
rownum <- nrow(MPT_BRIDGE_counts)
BRIDGE_pvalue_hets <- as.vector(c(), mode = "any")
for(i in 1:rownum){
  
  fishtable = matrix(c(MPT_BRIDGE_counts[i,10],
                       (nrow(case_list) - MPT_BRIDGE_counts[i,10]),
                       MPT_BRIDGE_counts[i,14],
                       (nrow(BRIDGE_list) - MPT_BRIDGE_counts[i,14])), 
                     nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue_hets <- append(BRIDGE_pvalue_hets, fisherout[[1]])
  
}
BRIDGE_qvalue_hets <- p.adjust(BRIDGE_pvalue_hets, method = "fdr", n = rownum)

##Fishers exact with fdr adjustment for each individual variant observed (homs cases vs BRIDGE)
rownum <- nrow(MPT_BRIDGE_counts)
BRIDGE_pvalue_homs <- as.vector(c(), mode = "any")
for(i in 1:rownum){
  
  fishtable = matrix(c(MPT_BRIDGE_counts[i,11],(nrow(case_list) - MPT_BRIDGE_counts[i,11]),MPT_BRIDGE_counts[i,15],(nrow(BRIDGE_list) - MPT_BRIDGE_counts[i,15])), 
                     nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue_homs <- append(BRIDGE_pvalue_homs, fisherout[[1]])
  
}
BRIDGE_qvalue_homs <- p.adjust(BRIDGE_pvalue_homs, method = "fdr", n = rownum)

##Fishers exact with fdr adjustment for each individual variant observed (hets_homs cases vs BRIDGE)
rownum <- nrow(MPT_BRIDGE_counts)
BRIDGE_pvalue_hets_homs <- as.vector(c(), mode = "any")
for(i in 1:rownum){
  
  fishtable = matrix(c(MPT_BRIDGE_counts[i,12],(nrow(case_list) - MPT_BRIDGE_counts[i,12]),MPT_BRIDGE_counts[i,16],(nrow(BRIDGE_list) - MPT_BRIDGE_counts[i,16])), 
                     nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue_hets_homs <- append(BRIDGE_pvalue_hets_homs, fisherout[[1]])
  
}
BRIDGE_qvalue_hets_homs <- p.adjust(BRIDGE_pvalue_hets_homs, method = "fdr", n = rownum)

##Fishers exact with fdr adjustment for each individual variant observed (hets_homs_agg cases vs BRIDGE)
rownum <- nrow(MPT_BRIDGE_counts)
BRIDGE_pvalue_hets_homs_agg <- as.vector(c(), mode = "any")
for(i in 1:rownum){
  
  fishtable = matrix(c(MPT_BRIDGE_counts[i,13],((nrow(case_list)*2) - MPT_BRIDGE_counts[i,13]),MPT_BRIDGE_counts[i,17],((nrow(BRIDGE_list)*2) - MPT_BRIDGE_counts[i,17])), 
                     nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue_hets_homs_agg <- append(BRIDGE_pvalue_hets_homs_agg, fisherout[[1]])
  
}
BRIDGE_qvalue_hets_homs_agg <- p.adjust(BRIDGE_pvalue_hets_homs_agg, method = "fdr", n = rownum)

##Collation of results for individual variants
tumour_query <- paste(args[1],args[2],args[3],args[4],args[5],args[6],args[7],args[8],args[9],sep = "")
tumour_col <- rep(tumour_query, nrow(MPT_BRIDGE))

MPT_BRIDGE_final <- cbind(BRIDGE_pvalue_hets,
                          BRIDGE_qvalue_hets,
                          BRIDGE_pvalue_homs,
                          BRIDGE_qvalue_homs,
                          BRIDGE_pvalue_hets_homs,
                          BRIDGE_qvalue_hets_homs,
                          BRIDGE_pvalue_hets_homs_agg,
                          BRIDGE_qvalue_hets_homs_agg,
                          MPT_BRIDGE_counts,
                          tumour_col)

colnames(MPT_BRIDGE_final)[ncol(MPT_BRIDGE_final)-1] <- paste("tumour_col_n=",(nrow(case_list)),sep = "")

##Output potentially significant results
MPT_BRIDGE_final$BRIDGE_qvalue_hets <- as.numeric(as.character(MPT_BRIDGE_final$BRIDGE_qvalue_hets))
MPT_BRIDGE_final$BRIDGE_qvalue_homs <- as.numeric(as.character(MPT_BRIDGE_final$BRIDGE_qvalue_homs))
MPT_BRIDGE_final$BRIDGE_qvalue_hets_homs <- as.numeric(as.character(MPT_BRIDGE_final$BRIDGE_qvalue_hets_homs))
MPT_BRIDGE_final$BRIDGE_qvalue_hets_homs_agg <- as.numeric(as.character(MPT_BRIDGE_final$BRIDGE_qvalue_hets_homs_agg))

BRIDGE_qvalue_hets_sig_index <- which(MPT_BRIDGE_final$BRIDGE_qvalue_hets <0.05)
BRIDGE_qvalue_homs_sig_index <- which(MPT_BRIDGE_final$BRIDGE_qvalue_homs <0.05)
BRIDGE_qvalue_hets_homs_sig_index <- which(MPT_BRIDGE_final$BRIDGE_qvalue_hets_homs <0.05)
BRIDGE_qvalue_hets_homs_agg_sig_index <- which(MPT_BRIDGE_final$BRIDGE_qvalue_hets_homs_agg <0.05)

sig_indexes <- unique(c(BRIDGE_qvalue_hets_sig_index,
                        BRIDGE_qvalue_homs_sig_index,
                        BRIDGE_qvalue_hets_homs_sig_index,
                        BRIDGE_qvalue_hets_homs_agg_sig_index
))

MPT_BRIDGE_final_sig <- MPT_BRIDGE_final[sig_indexes,]

##Per element analysis ("gene" term used in scripts as adapted from previous script)
genes <- unique(MPT_BRIDGE_final$variant_element_hits_IDs)
genes <- genes[!is.na(genes)]

##Counts of variants per element for cases and BRIDGE controls (hets)
variant_per_gene_count_cases_hets <- data.frame()

for(i in 1:length(genes)){
  a <- sum(MPT_BRIDGE_counts[,10][which(MPT_BRIDGE_counts$variant_element_hits_IDs %in% genes[i])]) 
  variant_per_gene_count_cases_hets  <- append(variant_per_gene_count_cases_hets,a)
}

variant_per_gene_count_cases_hets <- t(data.frame(variant_per_gene_count_cases_hets))
gene_names <- as.vector(genes)
rownames(variant_per_gene_count_cases_hets) <- gene_names

variant_per_gene_count_BRIDGES_hets <- data.frame()

for(i in 1:length(genes)){
  a <- sum(MPT_BRIDGE_counts[,14][which(MPT_BRIDGE_counts$variant_element_hits_IDs %in% genes[i])]) 
  variant_per_gene_count_BRIDGES_hets  <- append(variant_per_gene_count_BRIDGES_hets,a)
}

variant_per_gene_count_BRIDGES_hets <- t(data.frame(variant_per_gene_count_BRIDGES_hets))
gene_names <- as.vector(genes)
rownames(variant_per_gene_count_BRIDGES_hets) <- gene_names

##Fishers exact with fdr adjustment for counts of variants per element (hets)
cases_trials_per_gene <- as.vector(c(), mode = "any")
BRIDGE_pvalue_per_gene_hets <- as.vector(c(), mode = "any")
BRIDGE_trials_per_gene <- as.vector(c(), mode = "any")

for(i in 1:nrow(variant_per_gene_count_cases_hets))
{
  fishtable = matrix(
    c(variant_per_gene_count_cases_hets[i,1],
      sum(MPT_BRIDGE_counts$variant_element_hits_IDs %in% genes[i]) * nrow(case_list) - variant_per_gene_count_cases_hets[i,1],
      variant_per_gene_count_BRIDGES_hets[i,1],
      sum(MPT_BRIDGE_counts$variant_element_hits_IDs %in% genes[i]) * nrow(BRIDGE_list) - variant_per_gene_count_BRIDGES_hets[i,1]),
    
    nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue_per_gene_hets <- append(BRIDGE_pvalue_per_gene_hets, fisherout[[1]])
  cases_trials_per_gene <- append(cases_trials_per_gene, sum(MPT_BRIDGE_counts$variant_element_hits_IDs %in% genes[i]) * nrow(case_list))
  BRIDGE_trials_per_gene <- append(BRIDGE_trials_per_gene, sum(MPT_BRIDGE_counts$variant_element_hits_IDs %in% genes[i]) * nrow(BRIDGE_list))
}

BRIDGE_qvalue_per_gene_hets <- p.adjust(BRIDGE_pvalue_per_gene_hets, method = "fdr", n = nrow(elements))

cases_per_gene_prop_hets <- variant_per_gene_count_cases_hets / cases_trials_per_gene
BRIDGES_per_gene_prop_hets <- variant_per_gene_count_BRIDGES_hets / BRIDGE_trials_per_gene

##Counts of variants per element for cases and BRIDGE controls (homs)
variant_per_gene_count_cases_homs <- data.frame()

for(i in 1:length(genes)){
  a <- sum(MPT_BRIDGE_counts[,11][which(MPT_BRIDGE_counts$variant_element_hits_IDs %in% genes[i])]) 
  variant_per_gene_count_cases_homs  <- append(variant_per_gene_count_cases_homs,a)
}

variant_per_gene_count_cases_homs <- t(data.frame(variant_per_gene_count_cases_homs))
gene_names <- as.vector(genes)
rownames(variant_per_gene_count_cases_homs) <- gene_names

variant_per_gene_count_BRIDGES_homs <- data.frame()

for(i in 1:length(genes)){
  a <- sum(MPT_BRIDGE_counts[,15][which(MPT_BRIDGE_counts$variant_element_hits_IDs %in% genes[i])]) 
  variant_per_gene_count_BRIDGES_homs  <- append(variant_per_gene_count_BRIDGES_homs,a)
}

variant_per_gene_count_BRIDGES_homs <- t(data.frame(variant_per_gene_count_BRIDGES_homs))
gene_names <- as.vector(genes)
rownames(variant_per_gene_count_BRIDGES_homs) <- gene_names

##Fishers exact with fdr adjustment for counts of variants per element (homs)
BRIDGE_pvalue_per_gene_homs <- as.vector(c(), mode = "any")

for(i in 1:nrow(variant_per_gene_count_cases_homs))
{
  fishtable = matrix(
    c(variant_per_gene_count_cases_homs[i,1],
      sum(MPT_BRIDGE_counts$variant_element_hits_IDs %in% genes[i]) * nrow(case_list) - variant_per_gene_count_cases_homs[i,1],
      variant_per_gene_count_BRIDGES_homs[i,1],
      sum(MPT_BRIDGE_counts$variant_element_hits_IDs %in% genes[i]) * nrow(BRIDGE_list) - variant_per_gene_count_BRIDGES_homs[i,1]),
    
    nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue_per_gene_homs <- append(BRIDGE_pvalue_per_gene_homs, fisherout[[1]])
}

BRIDGE_qvalue_per_gene_homs <- p.adjust(BRIDGE_pvalue_per_gene_homs, method = "fdr", n = nrow(elements)) 

cases_per_gene_prop_homs <- variant_per_gene_count_cases_homs / cases_trials_per_gene
BRIDGES_per_gene_prop_homs <- variant_per_gene_count_BRIDGES_homs / BRIDGE_trials_per_gene

##Counts of variants per element for cases, BRIDGE and (hets_homs)
variant_per_gene_count_cases_hets_homs <- data.frame()

for(i in 1:length(genes)){
  a <- sum(MPT_BRIDGE_counts[,12][which(MPT_BRIDGE_counts$variant_element_hits_IDs %in% genes[i])]) 
  variant_per_gene_count_cases_hets_homs  <- append(variant_per_gene_count_cases_hets_homs,a)
}

variant_per_gene_count_cases_hets_homs <- t(data.frame(variant_per_gene_count_cases_hets_homs))
gene_names <- as.vector(genes)
rownames(variant_per_gene_count_cases_hets_homs) <- gene_names

variant_per_gene_count_BRIDGES_hets_homs <- data.frame()

for(i in 1:length(genes)){
  a <- sum(MPT_BRIDGE_counts[,16][which(MPT_BRIDGE_counts$variant_element_hits_IDs %in% genes[i])]) 
  variant_per_gene_count_BRIDGES_hets_homs  <- append(variant_per_gene_count_BRIDGES_hets_homs,a)
}

variant_per_gene_count_BRIDGES_hets_homs <- t(data.frame(variant_per_gene_count_BRIDGES_hets_homs))
gene_names <- as.vector(genes)
rownames(variant_per_gene_count_BRIDGES_hets_homs) <- gene_names

##Fishers exact with fdr adjustment for counts of variants per element (hets_homs)
BRIDGE_pvalue_per_gene_hets_homs <- as.vector(c(), mode = "any")

for(i in 1:nrow(variant_per_gene_count_cases_hets_homs))
{
  fishtable = matrix(
    c(variant_per_gene_count_cases_hets_homs[i,1],
      sum(MPT_BRIDGE_counts$variant_element_hits_IDs %in% genes[i]) * nrow(case_list) - variant_per_gene_count_cases_hets_homs[i,1],
      variant_per_gene_count_BRIDGES_hets_homs[i,1],
      sum(MPT_BRIDGE_counts$variant_element_hits_IDs %in% genes[i]) * nrow(BRIDGE_list) - variant_per_gene_count_BRIDGES_hets_homs[i,1]),
    
    nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue_per_gene_hets_homs <- append(BRIDGE_pvalue_per_gene_hets_homs, fisherout[[1]])
}

BRIDGE_qvalue_per_gene_hets_homs <- p.adjust(BRIDGE_pvalue_per_gene_hets_homs, method = "fdr", n = nrow(elements))

cases_per_gene_prop_hets_homs <- variant_per_gene_count_cases_hets_homs / cases_trials_per_gene
BRIDGES_per_gene_prop_hets_homs <- variant_per_gene_count_BRIDGES_hets_homs / BRIDGE_trials_per_gene

##Counts of variants per element for cases, BRIDGE (hets_homs_agg)
variant_per_gene_count_cases_hets_homs_agg <- data.frame()

for(i in 1:length(genes)){
  a <- sum(MPT_BRIDGE_counts[,13][which(MPT_BRIDGE_counts$variant_element_hits_IDs %in% genes[i])]) 
  variant_per_gene_count_cases_hets_homs_agg  <- append(variant_per_gene_count_cases_hets_homs_agg,a)
}

variant_per_gene_count_cases_hets_homs_agg <- t(data.frame(variant_per_gene_count_cases_hets_homs_agg))
gene_names <- as.vector(genes)
rownames(variant_per_gene_count_cases_hets_homs_agg) <- gene_names

variant_per_gene_count_BRIDGES_hets_homs_agg <- data.frame()

for(i in 1:length(genes)){
  a <- sum(MPT_BRIDGE_counts[,17][which(MPT_BRIDGE_counts$variant_element_hits_IDs %in% genes[i])]) 
  variant_per_gene_count_BRIDGES_hets_homs_agg  <- append(variant_per_gene_count_BRIDGES_hets_homs_agg,a)
}

variant_per_gene_count_BRIDGES_hets_homs_agg <- t(data.frame(variant_per_gene_count_BRIDGES_hets_homs_agg))
gene_names <- as.vector(genes)
rownames(variant_per_gene_count_BRIDGES_hets_homs_agg) <- gene_names

##Fishers exact with fdr adjustment for counts of variants per element (hets_homs_agg)
BRIDGE_pvalue_per_gene_hets_homs_agg <- as.vector(c(), mode = "any")


for(i in 1:nrow(variant_per_gene_count_cases_hets_homs_agg))
{
  fishtable = matrix(
    c(variant_per_gene_count_cases_hets_homs_agg[i,1],
      sum(MPT_BRIDGE_counts$variant_element_hits_IDs %in% genes[i]) * (nrow(case_list)*2) - variant_per_gene_count_cases_hets_homs_agg[i,1],
      variant_per_gene_count_BRIDGES_hets_homs_agg[i,1],
      sum(MPT_BRIDGE_counts$variant_element_hits_IDs %in% genes[i]) * (nrow(BRIDGE_list)*2) - variant_per_gene_count_BRIDGES_hets_homs_agg[i,1]),
    
    nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue_per_gene_hets_homs_agg <- append(BRIDGE_pvalue_per_gene_hets_homs_agg, fisherout[[1]])
}

BRIDGE_qvalue_per_gene_hets_homs_agg <- p.adjust(BRIDGE_pvalue_per_gene_hets_homs_agg, method = "fdr", n = nrow(elements)) 

cases_per_gene_prop_hets_homs_agg <- variant_per_gene_count_cases_hets_homs_agg / (cases_trials_per_gene *2)
BRIDGES_per_gene_prop_hets_homs_agg <- variant_per_gene_count_BRIDGES_hets_homs_agg / (BRIDGE_trials_per_gene *2)

##Counts for individuals with variants per element for cases and controls (hets).
trimmed_cases_hets <- data.frame(MPT_BRIDGE_counts$cases_hets)

indv_with_variant_per_gene_cases_hets <- data.frame()

for(i in 1:length(genes)){
  a <- trimmed_cases_hets[which(MPT_BRIDGE_counts$variant_element_hits_IDs %in% genes[i]),]
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

trimmed_BRIDGES_hets <- data.frame(MPT_BRIDGE_counts$BRIDGES_hets)

indv_with_variant_per_gene_BRIDGES_hets <- data.frame()

for(i in 1:length(genes)){
  a <- trimmed_BRIDGES_hets[which(MPT_BRIDGE_counts$variant_element_hits_IDs %in% genes[i]),]
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

##Fishers exact with fdr adjustment for counts for individuals with variants per element (hets)
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

BRIDGE_qvalue_per_indv_per_gene_hets <- p.adjust(BRIDGE_pvalue_per_indv_per_gene_hets, method = "fdr", n = nrow(elements)) 
cases_per_indv_per_gene_prop_hets <- indv_per_gene_count_cases_hets / nrow(case_list)
BRIDGES_per_indv_per_gene_prop_hets <- indv_per_gene_count_BRIDGES_hets / nrow(BRIDGE_list)

##Counts for individuals with variants per element for cases and controls (homs).
trimmed_cases_homs <- data.frame(MPT_BRIDGE_counts$cases_homs)

indv_with_variant_per_gene_cases_homs <- data.frame()

for(i in 1:length(genes)){
  a <- trimmed_cases_homs[which(MPT_BRIDGE_counts$variant_element_hits_IDs %in% genes[i]),]
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

trimmed_BRIDGES_homs <- data.frame(MPT_BRIDGE_counts$BRIDGES_homs)

indv_with_variant_per_gene_BRIDGES_homs <- data.frame()

for(i in 1:length(genes)){
  a <- trimmed_BRIDGES_homs[which(MPT_BRIDGE_counts$variant_element_hits_IDs %in% genes[i]),]
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

##Fishers exact with fdr adjustment for counts for individuals with variants per element (homs)
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

BRIDGE_qvalue_per_indv_per_gene_homs <- p.adjust(BRIDGE_pvalue_per_indv_per_gene_homs, method = "fdr", n = nrow(elements))
cases_per_indv_per_gene_prop_homs <- indv_per_gene_count_cases_homs / nrow(case_list)
BRIDGES_per_indv_per_gene_prop_homs <- indv_per_gene_count_BRIDGES_homs / nrow(BRIDGE_list)

##Counts for individuals with variants per element for cases (hets_homs).
trimmed_cases_hets_homs <- paste(trimmed_cases_hets$MPT_BRIDGE_counts.cases_hets, trimmed_cases_homs$MPT_BRIDGE_counts.cases_homs, sep = ";")
trimmed_cases_hets_homs <- data.frame(trimmed_cases_hets_homs)

indv_with_variant_per_gene_cases_hets_homs <- data.frame()

for(i in 1:length(genes)){
  a <- trimmed_cases_hets_homs[which(MPT_BRIDGE_counts$variant_element_hits_IDs %in% genes[i]),]
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

##Counts for individuals with variants per element for BRIDGE controls (hets_homs).
trimmed_BRIDGES_hets_homs <- paste(trimmed_BRIDGES_hets$MPT_BRIDGE_counts.BRIDGES_hets, trimmed_BRIDGES_homs$MPT_BRIDGE_counts.BRIDGES_homs, sep = ";")
trimmed_BRIDGES_hets_homs <- data.frame(trimmed_BRIDGES_hets_homs)

indv_with_variant_per_gene_BRIDGES_hets_homs <- data.frame()

for(i in 1:length(genes)){
  a <- trimmed_BRIDGES_hets_homs[which(MPT_BRIDGE_counts$variant_element_hits_IDs %in% genes[i]),]
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

##Fishers exact with fdr adjustment for counts for individuals with variants per element (hets_homs)
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

BRIDGE_qvalue_per_indv_per_gene_hets_homs <- p.adjust(BRIDGE_pvalue_per_indv_per_gene_hets_homs, method = "fdr", n = nrow(elements)) 
cases_per_indv_per_gene_prop_hets_homs <- indv_per_gene_count_cases_hets_homs / nrow(case_list)
BRIDGES_per_indv_per_gene_prop_hets_homs <- indv_per_gene_count_BRIDGES_hets_homs / nrow(BRIDGE_list)

##Collate results
tumour_query <- paste(args[1],args[2],args[3],args[4],args[5],args[6],args[7],args[8],args[9],sep = "")
tumour_col <- rep(tumour_query, nrow(variant_per_gene_count_cases_hets))
variants_per_gene_table <- as.data.frame(cbind(rownames(variant_per_gene_count_cases_hets),
                                               
                                               tumour_col,
                                               
                                               variant_per_gene_count_cases_hets,
                                               variant_per_gene_count_BRIDGES_hets,
                                               BRIDGE_pvalue_per_gene_hets,
                                               BRIDGE_qvalue_per_gene_hets,
                                               cases_per_gene_prop_hets,
                                               BRIDGES_per_gene_prop_hets,
                                               
                                               variant_per_gene_count_cases_homs,
                                               variant_per_gene_count_BRIDGES_homs,
                                               BRIDGE_pvalue_per_gene_homs,
                                               BRIDGE_qvalue_per_gene_homs,
                                               cases_per_gene_prop_homs,
                                               BRIDGES_per_gene_prop_homs,
                                               
                                               variant_per_gene_count_cases_hets_homs,
                                               variant_per_gene_count_BRIDGES_hets_homs,
                                               BRIDGE_pvalue_per_gene_hets_homs,
                                               BRIDGE_qvalue_per_gene_hets_homs,
                                               cases_per_gene_prop_hets_homs,
                                               BRIDGES_per_gene_prop_hets_homs,
                                               
                                               variant_per_gene_count_cases_hets_homs_agg,
                                               variant_per_gene_count_BRIDGES_hets_homs_agg,
                                               BRIDGE_pvalue_per_gene_hets_homs_agg,
                                               BRIDGE_qvalue_per_gene_hets_homs_agg,
                                               cases_per_gene_prop_hets_homs_agg,
                                               BRIDGES_per_gene_prop_hets_homs_agg,
                                               
                                               indv_per_gene_count_cases_hets,
                                               indv_per_gene_count_BRIDGES_hets,
                                               BRIDGE_pvalue_per_indv_per_gene_hets,
                                               BRIDGE_qvalue_per_indv_per_gene_hets,
                                               cases_per_indv_per_gene_prop_hets,
                                               BRIDGES_per_indv_per_gene_prop_hets,
                                               
                                               indv_per_gene_count_cases_homs,
                                               indv_per_gene_count_BRIDGES_homs,
                                               BRIDGE_pvalue_per_indv_per_gene_homs,
                                               BRIDGE_qvalue_per_indv_per_gene_homs,
                                               cases_per_indv_per_gene_prop_homs,
                                               BRIDGES_per_indv_per_gene_prop_homs,
                                               
                                               indv_per_gene_count_cases_hets_homs,
                                               indv_per_gene_count_BRIDGES_hets_homs,
                                               BRIDGE_pvalue_per_indv_per_gene_hets_homs,
                                               BRIDGE_qvalue_per_indv_per_gene_hets_homs,
                                               cases_per_indv_per_gene_prop_hets_homs,
                                               BRIDGES_per_indv_per_gene_prop_hets_homs))

colnames(variants_per_gene_table) <- c("gene", 
                                       
                                       paste("tumour_col_n=",(nrow(case_list)),sep = ""),
                                       
                                       "variant_per_gene_count_cases_hets",
                                       "variant_per_gene_count_BRIDGES_hets",
                                       "BRIDGE_pvalue_per_gene_hets",
                                       "BRIDGE_qvalue_per_gene_hets",
                                       "cases_per_gene_prop_hets",
                                       "BRIDGES_per_gene_prop_hets",
                                       
                                       "variant_per_gene_count_cases_homs",
                                       "variant_per_gene_count_BRIDGES_homs",
                                       "BRIDGE_pvalue_per_gene_homs",
                                       "BRIDGE_qvalue_per_gene_homs",
                                       "cases_per_gene_prop_homs",
                                       "BRIDGES_per_gene_prop_homs",
                                       
                                       "variant_per_gene_count_cases_hets_homs",
                                       "variant_per_gene_count_BRIDGES_hets_homs",
                                       "BRIDGE_pvalue_per_gene_hets_homs",
                                       "BRIDGE_qvalue_per_gene_hets_homs",
                                       "cases_per_gene_prop_hets_homs",
                                       "BRIDGES_per_gene_prop_hets_homs",
                                       
                                       "variant_per_gene_count_cases_hets_homs_agg",
                                       "variant_per_gene_count_BRIDGES_hets_homs_agg",
                                       "BRIDGE_pvalue_per_gene_hets_homs_agg",
                                       "BRIDGE_qvalue_per_gene_hets_homs_agg",
                                       "cases_per_gene_prop_hets_homs_agg",
                                       "BRIDGES_per_gene_prop_hets_homs_agg",
                                       
                                       "indv_per_gene_count_cases_hets",
                                       "indv_per_gene_count_BRIDGES_hets",
                                       "BRIDGE_pvalue_per_indv_per_gene_hets",
                                       "BRIDGE_qvalue_per_indv_per_gene_hets",
                                       "cases_per_indv_per_gene_prop_hets",
                                       "BRIDGES_per_indv_per_gene_prop_hets",
                                       
                                       "indv_per_gene_count_cases_homs",
                                       "indv_per_gene_count_BRIDGES_homs",
                                       "BRIDGE_pvalue_per_indv_per_gene_homs",
                                       "BRIDGE_qvalue_per_indv_per_gene_homs",
                                       "cases_per_indv_per_gene_prop_homs",
                                       "BRIDGES_per_indv_per_gene_prop_homs",
                                       
                                       "indv_per_gene_count_cases_hets_homs",
                                       "indv_per_gene_count_BRIDGES_hets_homs",
                                       "BRIDGE_pvalue_per_indv_per_gene_hets_homs",
                                       "BRIDGE_qvalue_per_indv_per_gene_hets_homs",
                                       "cases_per_indv_per_gene_prop_hets_homs",
                                       "BRIDGES_per_indv_per_gene_prop_hets_homs")

variants_per_gene_table$BRIDGE_qvalue_per_gene_hets <- as.numeric(as.character (variants_per_gene_table$BRIDGE_qvalue_per_gene_hets))
variants_per_gene_table$BRIDGE_qvalue_per_gene_homs <- as.numeric(as.character (variants_per_gene_table$BRIDGE_qvalue_per_gene_homs))
variants_per_gene_table$BRIDGE_qvalue_per_gene_hets_homs <- as.numeric(as.character (variants_per_gene_table$BRIDGE_qvalue_per_gene_hets_homs))
variants_per_gene_table$BRIDGE_qvalue_per_gene_hets_homs_agg <- as.numeric(as.character (variants_per_gene_table$BRIDGE_qvalue_per_gene_hets_homs_agg))
variants_per_gene_table$BRIDGE_qvalue_per_indv_per_gene_hets <- as.numeric(as.character (variants_per_gene_table$BRIDGE_qvalue_per_indv_per_gene_hets))
variants_per_gene_table$BRIDGE_qvalue_per_indv_per_gene_homs <- as.numeric(as.character (variants_per_gene_table$BRIDGE_qvalue_per_indv_per_gene_homs))
variants_per_gene_table$BRIDGE_qvalue_per_indv_per_gene_hets_homs <- as.numeric(as.character (variants_per_gene_table$BRIDGE_qvalue_per_indv_per_gene_hets_homs))


per_gene_sig_index_hets <- which(variants_per_gene_table$BRIDGE_qvalue_per_gene_hets < 0.05)
per_gene_sig_index_homs <- which(variants_per_gene_table$BRIDGE_qvalue_per_gene_homs < 0.05)
per_gene_sig_index_hets_homs <- which(variants_per_gene_table$BRIDGE_qvalue_per_gene_hets_homs < 0.05)
per_gene_sig_index_hets_homs_agg <- which(variants_per_gene_table$BRIDGE_qvalue_per_gene_hets_homs_agg < 0.05)

per_indv_per_gene_sig_index_hets <- which(variants_per_gene_table$BRIDGE_qvalue_per_indv_per_gene_hets < 0.05) 
per_indv_per_gene_sig_index_homs <- which(variants_per_gene_table$BRIDGE_qvalue_per_indv_per_gene_homs < 0.05) 
per_indv_per_gene_sig_index_hets_homs <- which(variants_per_gene_table$BRIDGE_qvalue_per_indv_per_gene_hets_homs < 0.05) 

sig_index <- c(per_gene_sig_index_hets, 
               per_gene_sig_index_homs, 
               per_gene_sig_index_hets_homs, 
               per_gene_sig_index_hets_homs_agg, 
               per_indv_per_gene_sig_index_hets,
               per_indv_per_gene_sig_index_homs,
               per_indv_per_gene_sig_index_hets_homs)

sig_index <- unique(sig_index)
variants_per_gene_table_sig<- variants_per_gene_table[sig_index,]

##Output results
write.csv(MPT_BRIDGE_final_sig, paste("/home/jww39/non_coding/UCNE_enhancers_promoters/with_internal_af_filter/enhancers_promoters_only/", tumour_query,".csv", sep = "")) 
write.csv(variants_per_gene_table_sig, paste("/home/jww39/non_coding/UCNE_enhancers_promoters/with_internal_af_filter/enhancers_promoters_only/", tumour_query,"_variants_per_gene_table.csv", sep = "")) 

write.csv(MPT_BRIDGE_final, paste("/home/jww39/non_coding/UCNE_enhancers_promoters/with_internal_af_filter/enhancers_promoters_only/", tumour_query,"_all.csv", sep = "")) 
write.csv(variants_per_gene_table, paste("/home/jww39/non_coding/UCNE_enhancers_promoters/with_internal_af_filter/enhancers_promoters_only/", tumour_query,"_variants_per_gene_table_all.csv", sep = "")) 


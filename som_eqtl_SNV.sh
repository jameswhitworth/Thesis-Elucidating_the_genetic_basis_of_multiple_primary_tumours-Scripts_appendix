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

##Filter per chromosome VCF files using BED file corresponding to coordinates of interest. BED file containing all coordinates for gene lists
for i in `cat chromosome_files.txt`; do
bcftools view -R somatic_eQTL.bed -e 'FILTER!="PASS"' -o ${i}_filtered.vcf -O v /scratch/WGS10K/data/release/20170614-A/merged-vcf/no_hgmd/${i}.vcf.gz
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
bcftools filter -e 'FMT/DP<10 || FMT/GQ<30 || (FMT/AD[1]/(FMT/AD[0] + FMT/AD[1]) <0.3) || (FMT/AD[1]/(FMT/AD[0] + FMT/AD[1] + FMT/AD[2]) <0.3) || (FMT/AD[2]/(FMT/AD[0] + FMT/AD[1] + FMT/AD[2]) <0.3)' -S . -o somatic_eQTL_merged_filtered.vcf.gz -O z merged.vcf.gz

##Transferred to local server

######################
##Further  filtering##
######################

##Execute VEP filter script 
/home/jww39/ensembl-vep/filter_vep -i VEP_out_somatic_eQTL_merged.vcf -o VEP_out_filtered_somatic_eQTL_merged.vcf -filter "(EUR_AF < 0.01 or not EUR_AF) and (UK10KWGS_AF < 0.01 or not UK10KWGS_AF) and (WGS10K_AF < 0.05 or not WGS10K_AF)" --only_matched

#Prepare resulting VCF for reading into R
sed -i 's/#CHROM/CHROM/g' VEP_out_filtered_somatic_eQTL_merged.vcf
sed '/^#/ d' VEP_out_filtered_somatic_eQTL_merged.vcf > VEP_out_filtered_somatic_eQTL_merged_without_header.vcf
mv VEP_out_filtered_somatic_eQTL_merged.vcf VEP_out_filtered_somatic_eQTL_merged_with_header.vcf
mv VEP_out_filtered_somatic_eQTL_merged_without_header.vcf  VEP_out_filtered_somatic_eQTL_merged.vcf

######################################################################################################################
##R script to identify samples with each variant in table, count variants and compare frequency in cases vs controls##
######################################################################################################################

#########################
##Heterozygous variants##
#########################

############################################################################################################################
##Following commands used to launch script according to phenotypic subgroup:
#Rscript som_eQTL_lookup_hets.R ALL XXX XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE ACC XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE Aerodigestivetract XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE Bladder XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE Bonebenign XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE Breast XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE Cervix XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE CNS XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE CNShaemangioblastoma XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE CNSmeningioma XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE CNSNervesheath XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE Colorectal XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE Colorectalpolyps XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE Connectivetissuesofttissuesarcoma XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE Endometrium XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE GINET XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE GIST XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE Haemlymphoid XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE Haemmyeloid XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE Kidney XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE Kidneyoncocytoma XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE Lung XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE Lungcarcinoid XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE Melanoma XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE NMSC XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE Oesophagus XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE Ovary XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE Ovarysexcord-gonadalstromal XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE Pancreas XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE Paraganglioma XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE Parathyroid XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE Phaeochromocytoma XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE Pituitary XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE PNET XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE PNSNervesheathbenign XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE Prostate XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE Salivarygland XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE Smallbowel XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE Testicular XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE Thyroid XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R SINGLE Uvealmelanoma XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R BOTH Breast Colorectal XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R BOTH Breast Breast XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R BOTH Breast NMSC XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R BOTH Breast Endometrium XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R BOTH Breast Ovary XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R BOTH Breast Haemlymphoid XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R BOTH Breast Melanoma XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R BOTH Breast Thyroid XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R BOTH Endometrium Ovary XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R BOTH Breast Kidney XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R BOTH Kidney Kidney XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R BOTH Colorectal NMSC XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R BOTH NMSC NMSC XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R BOTH Breast Lung XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R BOTH Colorectal Colorectal XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R BOTH Breast Connectivetissuesofttissuesarcoma XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R BOTH Colorectal Endometrium XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R BOTH Kidney Thyroid XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R BOTH Breast CNSmeningioma XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R BOTH Paraganglioma Paraganglioma XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R BOTH Breast Cervix XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R BOTH Colorectal Prostate XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R BOTH Colorectal Thyroid XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R BOTH Kidney Lung XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM2 Breast Ovary XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM2 Colorectal Endometrium XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM2 Thyroidmedullary Phaeochromocytoma XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM2 NMSC Melanoma XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM2 Colorectal Gastric XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM2 Parathyroid Bonebenign XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM2 Breast Gastric XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM2 Thyroid Pituitary XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM2 NMSC Haemlymphoid XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM2 Haemlymphoid Haemmyeloid XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM2 Connectivetissuesofttissuesarcoma Bladder XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM2 Breast Pancreas XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM2 NMSC Bonesarcoma XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM2 NMSC Haemlymphoid XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM3 Haemmyeloid Aerodigestivetract Anus XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM3 Melanoma Pancreas CNS XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM3 Kidney Kidneyangiomyolipoma CNS XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM3 CNSmeningioma CNS CNSNervesheath XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM3 Haemlymphoid CNS Connectivetissuesofttissuesarcoma XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM3 Breast Thyroid Endometrium XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM3 Phaeochromocytoma Paraganglioma GIST XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM3 Wilms Connectivetissuesofttissuesarcoma Haemmyeloid XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM3 Cardiacmyxoma Thyroid Ovarysexcord-gonadalstromal XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM3 CNS PNSNervesheath PNSNervesheathbenign XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM4 Kidney Phaeochromocytoma Paraganglioma CNShaemangioblastoma XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM4 CNS CNShaemangioblastoma CNSmeningioma CNSNervesheath XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM4 Kidney Uterineleiomyoma Uterinesarcoma Cutaneousleiomyoma XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM4 Kidney Adrenaloncocytoma Kidneyoncocytoma Fibrofolliculoma XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM4 Haemmyeloid Aerodigestivetract Anus Melanoma XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM4 Breast Aerodigestivetract Lung Ovary XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM4 Colorectal Breast Gastric Ovarysexcord-gonadalstromal XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM4 Colorectal Endometrium Ovary Sebaceous XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM5 Haemmyeloid Aerodigestivetract Oesophagus Cervix Penis XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM5 Connectivetissuesofttissuesarcoma Bonesarcoma GIST Skinsarcoma Uterinesarcoma XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM6 Uvealmelanoma Kidney Melanoma Lung Mesothelioma CNSmeningioma XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM6 GINET Lungcarcinoid Ovaryneuroendocrine Paraganglioma Phaeochromocytoma PNET XXX XXX
#Rscript som_eQTL_lookup_hets.R 1FROM7 Retinoblastoma Connectivetissuesofttissuesarcoma Bonesarcoma GIST Skinsarcoma Uterinesarcoma Melanoma XXX
#Rscript som_eQTL_lookup_hets.R 1FROM7 Connectivetissuesofttissuesarcoma Bonesarcoma GIST Skinsarcoma Uterinesarcoma Melanoma Thyroid XXX
#Rscript som_eQTL_lookup_hets.R 1FROM8 Pituitary Parathyroid ACC GINET Lungcarcinoid Ovaryneuroendocrine Paraganglioma Phaeochromocytoma
#Rscript som_eQTL_lookup_hets.R 1FROM8 Breast ACC CNS Connectivetissuesofttissuesarcoma Bonesarcoma GIST Skinsarcoma Uterinesarcoma
#Rscript som_eQTL_lookup_hets.R 2FROM3 Breast Thyroid Endometrium XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 2FROM4 Breast Aerodigestivetract Lung Ovary XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 2FROM4 Colorectal Breast Gastric Ovarysexcord-gonadalstromal XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 2FROM4 Colorectal Endometrium Ovary Sebaceous XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets.R 2FROM6 Uvealmelanoma Kidney Melanoma Lung Mesothelioma CNSmeningioma XXX XXX
#Rscript som_eQTL_lookup_hets.R 2FROM7 Connectivetissuesofttissuesarcoma Bonesarcoma GIST Skinsarcoma Uterinesarcoma Melanoma Thyroid XXX
#Rscript som_eQTL_lookup_hets.R 2FROM8 Breast ACC CNS Connectivetissuesofttissuesarcoma Bonesarcoma GIST Skinsarcoma Uterinesarcoma
#############################################################################################################################################

args = commandArgs(trailingOnly=TRUE)

##Read in necessary files
som_eQTL <- read.csv("/home/jww39/non_coding/somatic_eQTL/somatic_eQTLS.csv", header = T)
var <- read.delim("VEP_out_filtered_somatic_eQTL_merged.vcf", header = T)

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

##Matching variants to element that they are within
var_max_POS <- var$POS + nchar(as.character(var$REF))
var_chrom_pos_only <- data.frame(var$CHROM, var$POS, var_max_POS, paste(var$CHROM, var$POS, sep = "_"))

colnames(var_chrom_pos_only) <- c("CHROM", "POS", "max_POS", "CHROM_POS")

colnames(som_eQTL)[3] <- "chrom"
colnames(som_eQTL)[4] <- "start"
colnames(som_eQTL)[6] <- "end"

library(sqldf)

variant_element_match <- sqldf("select * from var_chrom_pos_only f1 
                               left join som_eQTL f2 on 
                               (
(f1.POS >= f2.start and f1.POS <= f2.end and f1.CHROM==f2.chrom) or 
(f1.max_POS >= f2.start and f1.max_POS <= f2.end and f1.CHROM==f2.chrom) or
(f1.POS <= f2.start and f1.max_POS >= f2.end and f1.CHROM==f2.chrom)
)")

##Variant counting

#Make list of control samples

control_list <- read.csv("BRIDGE_euro_controls.txt", header = F)

##Identifying samples (cases and controls) containing each variant
var_genotypes_only <- var[,10:ncol(var)]


samples_raw <- data.frame()
for (i in 1:nrow(var_genotypes_only)) {
  
  a <- unlist(var_genotypes_only[i,])
  b <- grep("^0\\/1", a, perl = TRUE)
  c <- as.vector(names(a[b]))
  d <- paste(c, collapse = ';')
  samples_raw <- append(samples_raw, d)
  
} 
samples <- unlist(samples_raw)

##Extracting which samples are among the cases

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

##Counting how many cases have the variant
library(stringr)

cases_df <- as.data.frame(cases)
df_for_cases_variant_count <- data.frame()
for (i in 1:nrow(cases_df)) { 
  count_of_samples_with_variant <- (str_count(cases_df[i,], "[A-Z]0"))
  df_for_cases_variant_count <- as.vector(append(df_for_cases_variant_count, count_of_samples_with_variant))
}
cases_variant_count <- cbind(df_for_cases_variant_count)
cases_variant_count <- unlist(cases_variant_count)

##Extracting which samples are among the controls
controls_list <- data.frame()
for (i in 1:nrow(samples_df)) {
  
  controls_split <- strsplit(samples[i], ";") # Works. Is a character
  controls_split_unlisted <- unlist(controls_split)
  controls_samples <- controls_split_unlisted[(which(controls_split_unlisted %in% control_list[,1]))] #Gives items that match the list
  controls_samples_string <- paste(controls_samples, collapse = ";")
  controls_list <- as.vector(append(controls_list, controls_samples_string))
  
}
controls <- unlist(controls_list)

##Counting how many controls have the variant
controls_df <- as.data.frame(controls)
df_for_controls_variant_count <- data.frame()
for (i in 1:nrow(controls_df)) { 
  count_of_samples_with_variant <- (str_count(controls_df[i,], "[A-Z][0-9]"))
  df_for_controls_variant_count <- as.vector(append(df_for_controls_variant_count, count_of_samples_with_variant))
} 
controls_variant_count <- cbind(df_for_controls_variant_count)
controls_variant_count <- unlist(controls_variant_count)

##Produce a table with variant information and counts
variant_element_match_with_gt <- cbind(
  var[,1:9],
  variant_element_match,
  cases, 
  cases_variant_count, 
  controls, 
  controls_variant_count,
  var[,10:ncol(var)]
)
colnames(variant_element_match_with_gt)[29] <- paste("cases_variant_count","_n=",(nrow(case_list)),sep = "")
colnames(variant_element_match_with_gt)[31] <- paste("controls_variant_count_variant_count","_n=",(nrow(control_list)),sep = "")

##Fishers exact with fdr adjustment for each individual variant observed
rownum <- nrow(variant_element_match_with_gt)
BRIDGE_pvalue <- as.vector(c(), mode = "any")
for(i in 1:rownum){
  
  fishtable = matrix(c(variant_element_match_with_gt[i,29],(nrow(case_list) - variant_element_match_with_gt[i,29]),variant_element_match_with_gt[i,31],(nrow(control_list) - variant_element_match_with_gt[i,31])), 
                     nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue <- append(BRIDGE_pvalue, fisherout[[1]])
  
}
BRIDGE_qvalue <- p.adjust(BRIDGE_pvalue, method = "fdr", n = rownum)

##Collate results
variant_element_match_with_gt_and_pq <- cbind(BRIDGE_pvalue, BRIDGE_qvalue, variant_element_match_with_gt)
variant_element_match_with_gt_and_pq$BRIDGE_pvalue <- as.numeric(as.character(variant_element_match_with_gt_and_pq$BRIDGE_pvalue))
variant_element_match_with_gt_and_pq$BRIDGE_qvalue <- as.numeric(as.character(variant_element_match_with_gt_and_pq$BRIDGE_qvalue))
sig_index <- which(variant_element_match_with_gt_and_pq$BRIDGE_qvalue < 0.1)
variant_element_match_sig <- variant_element_match_with_gt_and_pq[sig_index,]

##Output results
options <- paste(args[1],args[2],args[3],args[4],args[5],args[6],args[7],args[8],args[9],sep = "")
write.csv(variant_element_match_sig, paste("/home/jww39/non_coding/somatic_eQTL/with_internal_af_filter/hets/", options,"_variant_element_match.csv", sep = ""))
write.csv(variant_element_match_with_gt_and_pq, paste("/home/jww39/non_coding/somatic_eQTL/with_internal_af_filter/hets/", options,"_variant_element_match_all.csv", sep = ""))


##Per gene analysis

##Counts of variants per gene for cases and controls. Export to table and do fishers exact on counts vs controls. 
genes <- unique(variant_element_match_with_gt$Gene)
genes <- genes[!is.na(genes)]

variant_per_gene_count_cases <- data.frame()

for(i in 1:length(genes)){
  a <- sum(variant_element_match_with_gt[,29][which(variant_element_match_with_gt$Gene %in% genes[i])]) 
  variant_per_gene_count_cases  <- append(variant_per_gene_count_cases,a)
}

variant_per_gene_count_cases <- t(data.frame(variant_per_gene_count_cases))
gene_names <- as.vector(genes)
rownames(variant_per_gene_count_cases) <- gene_names

variant_per_gene_count_controls <- data.frame()

for(i in 1:length(genes)){
  a <- sum(variant_element_match_with_gt[,31][which(variant_element_match_with_gt$Gene %in% genes[i])])
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
      sum(variant_element_match_with_gt$Gene %in% genes[i]) * nrow(case_list) - variant_per_gene_count_cases[i,1],
      variant_per_gene_count_controls[i,1],
      sum(variant_element_match_with_gt$Gene %in% genes[i]) * nrow(control_list) - variant_per_gene_count_controls[i,1]),
    
    nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue_per_gene <- append(BRIDGE_pvalue_per_gene, fisherout[[1]])
  cases_trials_per_gene <- append(cases_trials_per_gene, sum(variant_element_match_with_gt$Gene %in% genes[i]) * nrow(case_list))
  controls_trials_per_gene <- append(controls_trials_per_gene, sum(variant_element_match_with_gt$Gene %in% genes[i]) * nrow(control_list))
}

BRIDGE_qvalue_per_gene <- p.adjust(BRIDGE_pvalue_per_gene, method = "fdr", n = 27)

cases_per_gene_prop <- variant_per_gene_count_cases / cases_trials_per_gene
controls_per_gene_prop <- variant_per_gene_count_controls / controls_trials_per_gene

##Counts for individuals with variants per gene for cases and controls 
trimmed_cases <- data.frame(variant_element_match_with_gt$cases)

indv_with_variant_per_gene_cases <- data.frame()

for(i in 1:length(genes)){
  a <- trimmed_cases[which(variant_element_match_with_gt$Gene %in% genes[i]),]
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


trimmed_controls <- data.frame(variant_element_match_with_gt$controls)

indv_with_variant_per_gene_controls <- data.frame()

for(i in 1:length(genes)){
  a <- trimmed_controls[which(variant_element_match_with_gt$Gene %in% genes[i]),]
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

BRIDGE_qvalue_per_indv_per_gene <- p.adjust(BRIDGE_pvalue_per_indv_per_gene, method = "fdr", n = 27) 

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
options <- paste(args[1],args[2],args[3],args[4],args[5],args[6],args[7],args[8],args[9],sep = "")
write.csv(variants_per_gene_table_sig, paste("/home/jww39/non_coding/somatic_eQTL/with_internal_af_filter/hets/", options,"_variants_per_gene_table.csv", sep = "")) 
write.csv(variants_per_gene_table, paste("/home/jww39/non_coding/somatic_eQTL/with_internal_af_filter/hets/", options,"_variants_per_gene_table_all.csv", sep = "")) 


#######################
##Homozygous variants##
#######################

############################################################################################################################
##Following commands used to launch script according to phenotypic subgroup:
#Rscript som_eQTL_lookup_homs.R ALL XXX XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE ACC XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE Aerodigestivetract XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE Bladder XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE Bonebenign XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE Breast XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE Cervix XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE CNS XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE CNShaemangioblastoma XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE CNSmeningioma XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE CNSNervesheath XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE Colorectal XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE Colorectalpolyps XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE Connectivetissuesofttissuesarcoma XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE Endometrium XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE GINET XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE GIST XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE Haemlymphoid XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE Haemmyeloid XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE Kidney XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE Kidneyoncocytoma XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE Lung XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE Lungcarcinoid XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE Melanoma XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE NMSC XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE Oesophagus XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE Ovary XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE Ovarysexcord-gonadalstromal XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE Pancreas XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE Paraganglioma XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE Parathyroid XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE Phaeochromocytoma XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE Pituitary XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE PNET XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE PNSNervesheathbenign XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE Prostate XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE Salivarygland XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE Smallbowel XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE Testicular XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE Thyroid XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R SINGLE Uvealmelanoma XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R BOTH Breast Colorectal XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R BOTH Breast Breast XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R BOTH Breast NMSC XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R BOTH Breast Endometrium XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R BOTH Breast Ovary XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R BOTH Breast Haemlymphoid XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R BOTH Breast Melanoma XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R BOTH Breast Thyroid XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R BOTH Endometrium Ovary XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R BOTH Breast Kidney XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R BOTH Kidney Kidney XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R BOTH Colorectal NMSC XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R BOTH NMSC NMSC XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R BOTH Breast Lung XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R BOTH Colorectal Colorectal XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R BOTH Breast Connectivetissuesofttissuesarcoma XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R BOTH Colorectal Endometrium XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R BOTH Kidney Thyroid XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R BOTH Breast CNSmeningioma XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R BOTH Paraganglioma Paraganglioma XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R BOTH Breast Cervix XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R BOTH Colorectal Prostate XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R BOTH Colorectal Thyroid XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R BOTH Kidney Lung XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM2 Breast Ovary XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM2 Colorectal Endometrium XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM2 Thyroidmedullary Phaeochromocytoma XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM2 NMSC Melanoma XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM2 Colorectal Gastric XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM2 Parathyroid Bonebenign XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM2 Breast Gastric XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM2 Thyroid Pituitary XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM2 NMSC Haemlymphoid XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM2 Haemlymphoid Haemmyeloid XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM2 Connectivetissuesofttissuesarcoma Bladder XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM2 Breast Pancreas XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM2 NMSC Bonesarcoma XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM2 NMSC Haemlymphoid XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM3 Haemmyeloid Aerodigestivetract Anus XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM3 Melanoma Pancreas CNS XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM3 Kidney Kidneyangiomyolipoma CNS XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM3 CNSmeningioma CNS CNSNervesheath XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM3 Haemlymphoid CNS Connectivetissuesofttissuesarcoma XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM3 Breast Thyroid Endometrium XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM3 Phaeochromocytoma Paraganglioma GIST XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM3 Wilms Connectivetissuesofttissuesarcoma Haemmyeloid XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM3 Cardiacmyxoma Thyroid Ovarysexcord-gonadalstromal XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM3 CNS PNSNervesheath PNSNervesheathbenign XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM4 Kidney Phaeochromocytoma Paraganglioma CNShaemangioblastoma XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM4 CNS CNShaemangioblastoma CNSmeningioma CNSNervesheath XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM4 Kidney Uterineleiomyoma Uterinesarcoma Cutaneousleiomyoma XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM4 Kidney Adrenaloncocytoma Kidneyoncocytoma Fibrofolliculoma XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM4 Haemmyeloid Aerodigestivetract Anus Melanoma XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM4 Breast Aerodigestivetract Lung Ovary XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM4 Colorectal Breast Gastric Ovarysexcord-gonadalstromal XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM4 Colorectal Endometrium Ovary Sebaceous XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM5 Haemmyeloid Aerodigestivetract Oesophagus Cervix Penis XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM5 Connectivetissuesofttissuesarcoma Bonesarcoma GIST Skinsarcoma Uterinesarcoma XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM6 Uvealmelanoma Kidney Melanoma Lung Mesothelioma CNSmeningioma XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM6 GINET Lungcarcinoid Ovaryneuroendocrine Paraganglioma Phaeochromocytoma PNET XXX XXX
#Rscript som_eQTL_lookup_homs.R 1FROM7 Retinoblastoma Connectivetissuesofttissuesarcoma Bonesarcoma GIST Skinsarcoma Uterinesarcoma Melanoma XXX
#Rscript som_eQTL_lookup_homs.R 1FROM7 Connectivetissuesofttissuesarcoma Bonesarcoma GIST Skinsarcoma Uterinesarcoma Melanoma Thyroid XXX
#Rscript som_eQTL_lookup_homs.R 1FROM8 Pituitary Parathyroid ACC GINET Lungcarcinoid Ovaryneuroendocrine Paraganglioma Phaeochromocytoma
#Rscript som_eQTL_lookup_homs.R 1FROM8 Breast ACC CNS Connectivetissuesofttissuesarcoma Bonesarcoma GIST Skinsarcoma Uterinesarcoma
#Rscript som_eQTL_lookup_homs.R 2FROM3 Breast Thyroid Endometrium XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 2FROM4 Breast Aerodigestivetract Lung Ovary XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 2FROM4 Colorectal Breast Gastric Ovarysexcord-gonadalstromal XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 2FROM4 Colorectal Endometrium Ovary Sebaceous XXX XXX XXX XXX
#Rscript som_eQTL_lookup_homs.R 2FROM6 Uvealmelanoma Kidney Melanoma Lung Mesothelioma CNSmeningioma XXX XXX
#Rscript som_eQTL_lookup_homs.R 2FROM7 Connectivetissuesofttissuesarcoma Bonesarcoma GIST Skinsarcoma Uterinesarcoma Melanoma Thyroid XXX
#Rscript som_eQTL_lookup_homs.R 2FROM8 Breast ACC CNS Connectivetissuesofttissuesarcoma Bonesarcoma GIST Skinsarcoma Uterinesarcoma
#############################################################################################################################################

args = commandArgs(trailingOnly=TRUE)

##Read in necessary files
som_eQTL <- read.csv("/home/jww39/non_coding/somatic_eQTL/somatic_eQTLS.csv", header = T)
var <- read.delim("VEP_out_filtered_somatic_eQTL_merged.vcf", header = T)

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

##Matching variants to element that they are within
var_max_POS <- var$POS + nchar(as.character(var$REF))
var_chrom_pos_only <- data.frame(var$CHROM, var$POS, var_max_POS, paste(var$CHROM, var$POS, sep = "_"))

colnames(var_chrom_pos_only) <- c("CHROM", "POS", "max_POS", "CHROM_POS")

colnames(som_eQTL)[3] <- "chrom"
colnames(som_eQTL)[4] <- "start"
colnames(som_eQTL)[6] <- "end"

library(sqldf)

variant_element_match <- sqldf("select * from var_chrom_pos_only f1 
                               left join som_eQTL f2 on 
                               (
(f1.POS >= f2.start and f1.POS <= f2.end and f1.CHROM==f2.chrom) or 
(f1.max_POS >= f2.start and f1.max_POS <= f2.end and f1.CHROM==f2.chrom) or
(f1.POS <= f2.start and f1.max_POS >= f2.end and f1.CHROM==f2.chrom)
)")

##Variant counting

#Make list of control samples

control_list <- read.csv("BRIDGE_euro_controls.txt", header = F)

##Identifying samples (cases and controls) containing each variant
var_genotypes_only <- var[,10:ncol(var)]


samples_raw <- data.frame()
for (i in 1:nrow(var_genotypes_only)) {
  
  a <- unlist(var_genotypes_only[i,])
  b <- grep("^1\\/1", a, perl = TRUE)
  c <- as.vector(names(a[b]))
  d <- paste(c, collapse = ';')
  samples_raw <- append(samples_raw, d)
  
} 
samples <- unlist(samples_raw)

##Extracting which samples are among the cases

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

##Counting how many cases have the variant
library(stringr)

cases_df <- as.data.frame(cases)
df_for_cases_variant_count <- data.frame()
for (i in 1:nrow(cases_df)) { 
  count_of_samples_with_variant <- (str_count(cases_df[i,], "[A-Z]0"))
  df_for_cases_variant_count <- as.vector(append(df_for_cases_variant_count, count_of_samples_with_variant))
}
cases_variant_count <- cbind(df_for_cases_variant_count)
cases_variant_count <- unlist(cases_variant_count)

##Extracting which samples are among the controls
controls_list <- data.frame()
for (i in 1:nrow(samples_df)) {
  
  controls_split <- strsplit(samples[i], ";") # Works. Is a character
  controls_split_unlisted <- unlist(controls_split)
  controls_samples <- controls_split_unlisted[(which(controls_split_unlisted %in% control_list[,1]))] #Gives items that match the list
  controls_samples_string <- paste(controls_samples, collapse = ";")
  controls_list <- as.vector(append(controls_list, controls_samples_string))
  
}
controls <- unlist(controls_list)

##Counting how many controls have the variant
controls_df <- as.data.frame(controls)
df_for_controls_variant_count <- data.frame()
for (i in 1:nrow(controls_df)) { 
  count_of_samples_with_variant <- (str_count(controls_df[i,], "[A-Z][0-9]"))
  df_for_controls_variant_count <- as.vector(append(df_for_controls_variant_count, count_of_samples_with_variant))
} 
controls_variant_count <- cbind(df_for_controls_variant_count)
controls_variant_count <- unlist(controls_variant_count)

##Produce a table with variant information and counts
variant_element_match_with_gt <- cbind(
  var[,1:9],
  variant_element_match,
  cases, 
  cases_variant_count, 
  controls, 
  controls_variant_count,
  var[,10:ncol(var)]
)
colnames(variant_element_match_with_gt)[29] <- paste("cases_variant_count","_n=",(nrow(case_list)),sep = "")
colnames(variant_element_match_with_gt)[31] <- paste("controls_variant_count_variant_count","_n=",(nrow(control_list)),sep = "")

##Fishers exact with fdr adjustment for each individual variant observed
rownum <- nrow(variant_element_match_with_gt)
BRIDGE_pvalue <- as.vector(c(), mode = "any")
for(i in 1:rownum){
  
  fishtable = matrix(c(variant_element_match_with_gt[i,29],(nrow(case_list) - variant_element_match_with_gt[i,29]),variant_element_match_with_gt[i,31],(nrow(control_list) - variant_element_match_with_gt[i,31])), 
                     nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue <- append(BRIDGE_pvalue, fisherout[[1]])
  
}
BRIDGE_qvalue <- p.adjust(BRIDGE_pvalue, method = "fdr", n = rownum)

##Collate results
variant_element_match_with_gt_and_pq <- cbind(BRIDGE_pvalue, BRIDGE_qvalue, variant_element_match_with_gt)
variant_element_match_with_gt_and_pq$BRIDGE_pvalue <- as.numeric(as.character(variant_element_match_with_gt_and_pq$BRIDGE_pvalue))
variant_element_match_with_gt_and_pq$BRIDGE_qvalue <- as.numeric(as.character(variant_element_match_with_gt_and_pq$BRIDGE_qvalue))
sig_index <- which(variant_element_match_with_gt_and_pq$BRIDGE_qvalue < 0.1)
variant_element_match_sig <- variant_element_match_with_gt_and_pq[sig_index,]

##Output results
options <- paste(args[1],args[2],args[3],args[4],args[5],args[6],args[7],args[8],args[9],sep = "")
write.csv(variant_element_match_sig, paste("/home/jww39/non_coding/somatic_eQTL/with_internal_af_filter/hets/", options,"_variant_element_match.csv", sep = ""))
write.csv(variant_element_match_with_gt_and_pq, paste("/home/jww39/non_coding/somatic_eQTL/with_internal_af_filter/hets/", options,"_variant_element_match_all.csv", sep = ""))


##Per gene analysis

##Counts of variants per gene for cases and controls. Export to table and do fishers exact on counts vs controls. 
genes <- unique(variant_element_match_with_gt$Gene)
genes <- genes[!is.na(genes)]

variant_per_gene_count_cases <- data.frame()

for(i in 1:length(genes)){
  a <- sum(variant_element_match_with_gt[,29][which(variant_element_match_with_gt$Gene %in% genes[i])]) 
  variant_per_gene_count_cases  <- append(variant_per_gene_count_cases,a)
}

variant_per_gene_count_cases <- t(data.frame(variant_per_gene_count_cases))
gene_names <- as.vector(genes)
rownames(variant_per_gene_count_cases) <- gene_names

variant_per_gene_count_controls <- data.frame()

for(i in 1:length(genes)){
  a <- sum(variant_element_match_with_gt[,31][which(variant_element_match_with_gt$Gene %in% genes[i])])
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
      sum(variant_element_match_with_gt$Gene %in% genes[i]) * nrow(case_list) - variant_per_gene_count_cases[i,1],
      variant_per_gene_count_controls[i,1],
      sum(variant_element_match_with_gt$Gene %in% genes[i]) * nrow(control_list) - variant_per_gene_count_controls[i,1]),
    
    nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue_per_gene <- append(BRIDGE_pvalue_per_gene, fisherout[[1]])
  cases_trials_per_gene <- append(cases_trials_per_gene, sum(variant_element_match_with_gt$Gene %in% genes[i]) * nrow(case_list))
  controls_trials_per_gene <- append(controls_trials_per_gene, sum(variant_element_match_with_gt$Gene %in% genes[i]) * nrow(control_list))
}

BRIDGE_qvalue_per_gene <- p.adjust(BRIDGE_pvalue_per_gene, method = "fdr", n = 27)

cases_per_gene_prop <- variant_per_gene_count_cases / cases_trials_per_gene
controls_per_gene_prop <- variant_per_gene_count_controls / controls_trials_per_gene

##Counts for individuals with variants per gene for cases and controls 
trimmed_cases <- data.frame(variant_element_match_with_gt$cases)

indv_with_variant_per_gene_cases <- data.frame()

for(i in 1:length(genes)){
  a <- trimmed_cases[which(variant_element_match_with_gt$Gene %in% genes[i]),]
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


trimmed_controls <- data.frame(variant_element_match_with_gt$controls)

indv_with_variant_per_gene_controls <- data.frame()

for(i in 1:length(genes)){
  a <- trimmed_controls[which(variant_element_match_with_gt$Gene %in% genes[i]),]
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

BRIDGE_qvalue_per_indv_per_gene <- p.adjust(BRIDGE_pvalue_per_indv_per_gene, method = "fdr", n = 27) 

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
options <- paste(args[1],args[2],args[3],args[4],args[5],args[6],args[7],args[8],args[9],sep = "")
write.csv(variants_per_gene_table_sig, paste("/home/jww39/non_coding/somatic_eQTL/with_internal_af_filter/homs/", options,"_variants_per_gene_table.csv", sep = "")) 
write.csv(variants_per_gene_table, paste("/home/jww39/non_coding/somatic_eQTL/with_internal_af_filter/homs/", options,"_variants_per_gene_table_all.csv", sep = "")) 



#########################################################
##Heterozygous variants summed with homozygous variants##
#########################################################

############################################################################################################################
##Following commands used to launch script according to phenotypic subgroup:
#Rscript som_eQTL_lookup_hets_homs.R ALL XXX XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE ACC XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE Aerodigestivetract XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE Bladder XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE Bonebenign XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE Breast XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE Cervix XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE CNS XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE CNShaemangioblastoma XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE CNSmeningioma XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE CNSNervesheath XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE Colorectal XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE Colorectalpolyps XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE Connectivetissuesofttissuesarcoma XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE Endometrium XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE GINET XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE GIST XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE Haemlymphoid XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE Haemmyeloid XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE Kidney XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE Kidneyoncocytoma XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE Lung XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE Lungcarcinoid XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE Melanoma XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE NMSC XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE Oesophagus XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE Ovary XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE Ovarysexcord-gonadalstromal XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE Pancreas XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE Paraganglioma XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE Parathyroid XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE Phaeochromocytoma XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE Pituitary XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE PNET XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE PNSNervesheathbenign XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE Prostate XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE Salivarygland XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE Smallbowel XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE Testicular XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE Thyroid XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R SINGLE Uvealmelanoma XXX XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R BOTH Breast Colorectal XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R BOTH Breast Breast XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R BOTH Breast NMSC XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R BOTH Breast Endometrium XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R BOTH Breast Ovary XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R BOTH Breast Haemlymphoid XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R BOTH Breast Melanoma XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R BOTH Breast Thyroid XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R BOTH Endometrium Ovary XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R BOTH Breast Kidney XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R BOTH Kidney Kidney XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R BOTH Colorectal NMSC XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R BOTH NMSC NMSC XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R BOTH Breast Lung XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R BOTH Colorectal Colorectal XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R BOTH Breast Connectivetissuesofttissuesarcoma XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R BOTH Colorectal Endometrium XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R BOTH Kidney Thyroid XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R BOTH Breast CNSmeningioma XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R BOTH Paraganglioma Paraganglioma XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R BOTH Breast Cervix XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R BOTH Colorectal Prostate XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R BOTH Colorectal Thyroid XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R BOTH Kidney Lung XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM2 Breast Ovary XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM2 Colorectal Endometrium XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM2 Thyroidmedullary Phaeochromocytoma XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM2 NMSC Melanoma XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM2 Colorectal Gastric XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM2 Parathyroid Bonebenign XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM2 Breast Gastric XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM2 Thyroid Pituitary XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM2 NMSC Haemlymphoid XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM2 Haemlymphoid Haemmyeloid XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM2 Connectivetissuesofttissuesarcoma Bladder XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM2 Breast Pancreas XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM2 NMSC Bonesarcoma XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM2 NMSC Haemlymphoid XXX XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM3 Haemmyeloid Aerodigestivetract Anus XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM3 Melanoma Pancreas CNS XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM3 Kidney Kidneyangiomyolipoma CNS XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM3 CNSmeningioma CNS CNSNervesheath XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM3 Haemlymphoid CNS Connectivetissuesofttissuesarcoma XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM3 Breast Thyroid Endometrium XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM3 Phaeochromocytoma Paraganglioma GIST XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM3 Wilms Connectivetissuesofttissuesarcoma Haemmyeloid XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM3 Cardiacmyxoma Thyroid Ovarysexcord-gonadalstromal XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM3 CNS PNSNervesheath PNSNervesheathbenign XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM4 Kidney Phaeochromocytoma Paraganglioma CNShaemangioblastoma XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM4 CNS CNShaemangioblastoma CNSmeningioma CNSNervesheath XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM4 Kidney Uterineleiomyoma Uterinesarcoma Cutaneousleiomyoma XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM4 Kidney Adrenaloncocytoma Kidneyoncocytoma Fibrofolliculoma XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM4 Haemmyeloid Aerodigestivetract Anus Melanoma XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM4 Breast Aerodigestivetract Lung Ovary XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM4 Colorectal Breast Gastric Ovarysexcord-gonadalstromal XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM4 Colorectal Endometrium Ovary Sebaceous XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM5 Haemmyeloid Aerodigestivetract Oesophagus Cervix Penis XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM5 Connectivetissuesofttissuesarcoma Bonesarcoma GIST Skinsarcoma Uterinesarcoma XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM6 Uvealmelanoma Kidney Melanoma Lung Mesothelioma CNSmeningioma XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM6 GINET Lungcarcinoid Ovaryneuroendocrine Paraganglioma Phaeochromocytoma PNET XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM7 Retinoblastoma Connectivetissuesofttissuesarcoma Bonesarcoma GIST Skinsarcoma Uterinesarcoma Melanoma XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM7 Connectivetissuesofttissuesarcoma Bonesarcoma GIST Skinsarcoma Uterinesarcoma Melanoma Thyroid XXX
#Rscript som_eQTL_lookup_hets_homs.R 1FROM8 Pituitary Parathyroid ACC GINET Lungcarcinoid Ovaryneuroendocrine Paraganglioma Phaeochromocytoma
#Rscript som_eQTL_lookup_hets_homs.R 1FROM8 Breast ACC CNS Connectivetissuesofttissuesarcoma Bonesarcoma GIST Skinsarcoma Uterinesarcoma
#Rscript som_eQTL_lookup_hets_homs.R 2FROM3 Breast Thyroid Endometrium XXX XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 2FROM4 Breast Aerodigestivetract Lung Ovary XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 2FROM4 Colorectal Breast Gastric Ovarysexcord-gonadalstromal XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 2FROM4 Colorectal Endometrium Ovary Sebaceous XXX XXX XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 2FROM6 Uvealmelanoma Kidney Melanoma Lung Mesothelioma CNSmeningioma XXX XXX
#Rscript som_eQTL_lookup_hets_homs.R 2FROM7 Connectivetissuesofttissuesarcoma Bonesarcoma GIST Skinsarcoma Uterinesarcoma Melanoma Thyroid XXX
#Rscript som_eQTL_lookup_hets_homs.R 2FROM8 Breast ACC CNS Connectivetissuesofttissuesarcoma Bonesarcoma GIST Skinsarcoma Uterinesarcoma
#############################################################################################################################################

args = commandArgs(trailingOnly=TRUE)

##Read in necessary files
som_eQTL <- read.csv("/home/jww39/non_coding/somatic_eQTL/somatic_eQTLS.csv", header = T)
var <- read.delim("VEP_out_filtered_somatic_eQTL_merged.vcf", header = T)

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

##Matching variants to element that they are within
var_max_POS <- var$POS + nchar(as.character(var$REF))
var_chrom_pos_only <- data.frame(var$CHROM, var$POS, var_max_POS, paste(var$CHROM, var$POS, sep = "_"))

colnames(var_chrom_pos_only) <- c("CHROM", "POS", "max_POS", "CHROM_POS")

colnames(som_eQTL)[3] <- "chrom"
colnames(som_eQTL)[4] <- "start"
colnames(som_eQTL)[6] <- "end"

library(sqldf)

variant_element_match <- sqldf("select * from var_chrom_pos_only f1 
                               left join som_eQTL f2 on 
                               (
(f1.POS >= f2.start and f1.POS <= f2.end and f1.CHROM==f2.chrom) or 
(f1.max_POS >= f2.start and f1.max_POS <= f2.end and f1.CHROM==f2.chrom) or
(f1.POS <= f2.start and f1.max_POS >= f2.end and f1.CHROM==f2.chrom)
)")

##Variant counting

#Make list of control samples

control_list <- read.csv("BRIDGE_euro_controls.txt", header = F)

##Identifying samples (cases and controls) containing each variant
var_genotypes_only <- var[,10:ncol(var)]


samples_raw <- data.frame()
for (i in 1:nrow(var_genotypes_only)) {
  
  a <- unlist(var_genotypes_only[i,])
  b <- grep("^[0-9]\\/[0-9]", a, perl = TRUE)
  c <- as.vector(names(a[b]))
  d <- paste(c, collapse = ';')
  samples_raw <- append(samples_raw, d)
  
} 
samples <- unlist(samples_raw)

##Extracting which samples are among the cases

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

##Counting how many cases have the variant
library(stringr)

cases_df <- as.data.frame(cases)
df_for_cases_variant_count <- data.frame()
for (i in 1:nrow(cases_df)) { 
  count_of_samples_with_variant <- (str_count(cases_df[i,], "[A-Z]0"))
  df_for_cases_variant_count <- as.vector(append(df_for_cases_variant_count, count_of_samples_with_variant))
}
cases_variant_count <- cbind(df_for_cases_variant_count)
cases_variant_count <- unlist(cases_variant_count)

##Extracting which samples are among the controls
controls_list <- data.frame()
for (i in 1:nrow(samples_df)) {
  
  controls_split <- strsplit(samples[i], ";") # Works. Is a character
  controls_split_unlisted <- unlist(controls_split)
  controls_samples <- controls_split_unlisted[(which(controls_split_unlisted %in% control_list[,1]))] #Gives items that match the list
  controls_samples_string <- paste(controls_samples, collapse = ";")
  controls_list <- as.vector(append(controls_list, controls_samples_string))
  
}
controls <- unlist(controls_list)

##Counting how many controls have the variant
controls_df <- as.data.frame(controls)
df_for_controls_variant_count <- data.frame()
for (i in 1:nrow(controls_df)) { 
  count_of_samples_with_variant <- (str_count(controls_df[i,], "[A-Z][0-9]"))
  df_for_controls_variant_count <- as.vector(append(df_for_controls_variant_count, count_of_samples_with_variant))
} 
controls_variant_count <- cbind(df_for_controls_variant_count)
controls_variant_count <- unlist(controls_variant_count)

##Produce a table with variant information and counts
variant_element_match_with_gt <- cbind(
  var[,1:9],
  variant_element_match,
  cases, 
  cases_variant_count, 
  controls, 
  controls_variant_count,
  var[,10:ncol(var)]
)
colnames(variant_element_match_with_gt)[29] <- paste("cases_variant_count","_n=",(nrow(case_list)),sep = "")
colnames(variant_element_match_with_gt)[31] <- paste("controls_variant_count_variant_count","_n=",(nrow(control_list)),sep = "")

##Fishers exact with fdr adjustment for each individual variant observed
rownum <- nrow(variant_element_match_with_gt)
BRIDGE_pvalue <- as.vector(c(), mode = "any")
for(i in 1:rownum){
  
  fishtable = matrix(c(variant_element_match_with_gt[i,29],(nrow(case_list) - variant_element_match_with_gt[i,29]),variant_element_match_with_gt[i,31],(nrow(control_list) - variant_element_match_with_gt[i,31])), 
                     nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue <- append(BRIDGE_pvalue, fisherout[[1]])
  
}
BRIDGE_qvalue <- p.adjust(BRIDGE_pvalue, method = "fdr", n = rownum)

##Collate results
variant_element_match_with_gt_and_pq <- cbind(BRIDGE_pvalue, BRIDGE_qvalue, variant_element_match_with_gt)
variant_element_match_with_gt_and_pq$BRIDGE_pvalue <- as.numeric(as.character(variant_element_match_with_gt_and_pq$BRIDGE_pvalue))
variant_element_match_with_gt_and_pq$BRIDGE_qvalue <- as.numeric(as.character(variant_element_match_with_gt_and_pq$BRIDGE_qvalue))
sig_index <- which(variant_element_match_with_gt_and_pq$BRIDGE_qvalue < 0.1)
variant_element_match_sig <- variant_element_match_with_gt_and_pq[sig_index,]

##Output results
options <- paste(args[1],args[2],args[3],args[4],args[5],args[6],args[7],args[8],args[9],sep = "")
write.csv(variant_element_match_sig, paste("/home/jww39/non_coding/somatic_eQTL/with_internal_af_filter/hets/", options,"_variant_element_match.csv", sep = ""))
write.csv(variant_element_match_with_gt_and_pq, paste("/home/jww39/non_coding/somatic_eQTL/with_internal_af_filter/hets/", options,"_variant_element_match_all.csv", sep = ""))


##Per gene analysis

##Counts of variants per gene for cases and controls. Export to table and do fishers exact on counts vs controls. 
genes <- unique(variant_element_match_with_gt$Gene)
genes <- genes[!is.na(genes)]

variant_per_gene_count_cases <- data.frame()

for(i in 1:length(genes)){
  a <- sum(variant_element_match_with_gt[,29][which(variant_element_match_with_gt$Gene %in% genes[i])]) 
  variant_per_gene_count_cases  <- append(variant_per_gene_count_cases,a)
}

variant_per_gene_count_cases <- t(data.frame(variant_per_gene_count_cases))
gene_names <- as.vector(genes)
rownames(variant_per_gene_count_cases) <- gene_names

variant_per_gene_count_controls <- data.frame()

for(i in 1:length(genes)){
  a <- sum(variant_element_match_with_gt[,31][which(variant_element_match_with_gt$Gene %in% genes[i])])
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
      sum(variant_element_match_with_gt$Gene %in% genes[i]) * nrow(case_list) - variant_per_gene_count_cases[i,1],
      variant_per_gene_count_controls[i,1],
      sum(variant_element_match_with_gt$Gene %in% genes[i]) * nrow(control_list) - variant_per_gene_count_controls[i,1]),
    
    nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  BRIDGE_pvalue_per_gene <- append(BRIDGE_pvalue_per_gene, fisherout[[1]])
  cases_trials_per_gene <- append(cases_trials_per_gene, sum(variant_element_match_with_gt$Gene %in% genes[i]) * nrow(case_list))
  controls_trials_per_gene <- append(controls_trials_per_gene, sum(variant_element_match_with_gt$Gene %in% genes[i]) * nrow(control_list))
}

BRIDGE_qvalue_per_gene <- p.adjust(BRIDGE_pvalue_per_gene, method = "fdr", n = 27)

cases_per_gene_prop <- variant_per_gene_count_cases / cases_trials_per_gene
controls_per_gene_prop <- variant_per_gene_count_controls / controls_trials_per_gene

##Counts for individuals with variants per gene for cases and controls 
trimmed_cases <- data.frame(variant_element_match_with_gt$cases)

indv_with_variant_per_gene_cases <- data.frame()

for(i in 1:length(genes)){
  a <- trimmed_cases[which(variant_element_match_with_gt$Gene %in% genes[i]),]
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


trimmed_controls <- data.frame(variant_element_match_with_gt$controls)

indv_with_variant_per_gene_controls <- data.frame()

for(i in 1:length(genes)){
  a <- trimmed_controls[which(variant_element_match_with_gt$Gene %in% genes[i]),]
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

BRIDGE_qvalue_per_indv_per_gene <- p.adjust(BRIDGE_pvalue_per_indv_per_gene, method = "fdr", n = 27) 

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
options <- paste(args[1],args[2],args[3],args[4],args[5],args[6],args[7],args[8],args[9],sep = "")
write.csv(variants_per_gene_table_sig, paste("/home/jww39/non_coding/somatic_eQTL/with_internal_af_filter/hets_homs/", options,"_variants_per_gene_table.csv", sep = "")) 
write.csv(variants_per_gene_table, paste("/home/jww39/non_coding/somatic_eQTL/with_internal_af_filter/hets_homs/", options,"_variants_per_gene_table_all.csv", sep = "")) 



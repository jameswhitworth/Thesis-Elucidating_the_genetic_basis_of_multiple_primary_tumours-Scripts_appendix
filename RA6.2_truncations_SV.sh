############################################################################################################
##R script to count structural variant calls affecting genes of interest and compare frequency vs controls##
############################################################################################################

############################################################################################################
##Following commands used to launch script according to phenotypic subgroup:
#Rscript trunc_SV_detection.R ALL XXX XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE ACC XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE Aerodigestivetract XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE Bladder XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE Bonebenign XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE Breast XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE Cervix XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE CNS XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE CNShaemangioblastoma XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE CNSmeningioma XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE CNSNervesheath XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE Colorectal XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE Colorectalpolyps XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE Connectivetissuesofttissuesarcoma XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE Endometrium XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE GINET XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE GIST XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE Haemlymphoid XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE Haemmyeloid XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE Kidney XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE Kidneyoncocytoma XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE Lung XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE Lungcarcinoid XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE Melanoma XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE NMSC XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE Oesophagus XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE Ovary XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE Ovarysexcord-gonadalstromal XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE Pancreas XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE Paraganglioma XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE Parathyroid XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE Phaeochromocytoma XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE Pituitary XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE PNET XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE PNSNervesheathbenign XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE Prostate XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE Salivarygland XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE Smallbowel XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE Testicular XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE Thyroid XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R SINGLE Uvealmelanoma XXX XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R BOTH Breast Colorectal XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R BOTH Breast NMSC XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R BOTH Breast Endometrium XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R BOTH Breast Ovary XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R BOTH Breast Haemlymphoid XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R BOTH Breast Melanoma XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R BOTH Breast Thyroid XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R BOTH Endometrium Ovary XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R BOTH Breast Kidney XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R BOTH Kidney Kidney XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R BOTH Colorectal NMSC XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R BOTH NMSC NMSC XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R BOTH Breast Lung XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R BOTH Breast Connectivetissuesofttissuesarcoma XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R BOTH Colorectal Endometrium XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R BOTH Kidney Thyroid XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R BOTH Breast CNSmeningioma XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R BOTH Paraganglioma Paraganglioma XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R BOTH Breast Cervix XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R BOTH Colorectal Prostate XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R BOTH Colorectal Thyroid XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R BOTH Kidney Lung XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM2 Breast Ovary XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM2 Colorectal Endometrium XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM2 Thyroidmedullary Phaeochromocytoma XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM2 NMSC Melanoma XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM2 Colorectal Gastric XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM2 Parathyroid Bonebenign XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM2 Breast Gastric XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM2 Thyroid Pituitary XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM2 NMSC Haemlymphoid XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM2 Haemlymphoid Haemmyeloid XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM2 Connectivetissuesofttissuesarcoma Bladder XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM2 Breast Pancreas XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM2 NMSC Bonesarcoma XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM2 NMSC Haemlymphoid XXX XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM3 Haemmyeloid Aerodigestivetract Anus XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM3 Melanoma Pancreas CNS XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM3 Kidney Kidneyangiomyolipoma CNS XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM3 CNSmeningioma CNS CNSNervesheath XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM3 Haemlymphoid CNS Connectivetissuesofttissuesarcoma XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM3 Breast Thyroid Endometrium XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM3 Phaeochromocytoma Paraganglioma GIST XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM3 Wilms Connectivetissuesofttissuesarcoma Haemmyeloid XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM3 Cardiacmyxoma Thyroid Ovarysexcord-gonadalstromal XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM3 CNS PNSNervesheath PNSNervesheathbenign XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM4 Kidney Phaeochromocytoma Paraganglioma CNShaemangioblastoma XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM4 CNS CNShaemangioblastoma CNSmeningioma CNSNervesheath XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM4 Kidney Uterineleiomyoma Uterinesarcoma Cutaneousleiomyoma XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM4 Kidney Adrenaloncocytoma Kidneyoncocytoma Fibrofolliculoma XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM4 Haemmyeloid Aerodigestivetract Anus Melanoma XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM4 Breast Aerodigestivetract Lung Ovary XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM4 Colorectal Breast Gastric Ovarysexcord-gonadalstromal XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM4 Colorectal Endometrium Ovary Sebaceous XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM5 Haemmyeloid Aerodigestivetract Oesophagus Cervix Penis XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM5 Connectivetissuesofttissuesarcoma Bonesarcoma GIST Skinsarcoma Uterinesarcoma XXX XXX XXX
#Rscript trunc_SV_detection.R 1FROM6 Uvealmelanoma Kidney Melanoma Lung Mesothelioma CNSmeningioma XXX XXX
#Rscript trunc_SV_detection.R 1FROM6 GINET Lungcarcinoid Ovaryneuroendocrine Paraganglioma Phaeochromocytoma PNET XXX XXX
#Rscript trunc_SV_detection.R 1FROM7 Retinoblastoma Connectivetissuesofttissuesarcoma Bonesarcoma GIST Skinsarcoma Uterinesarcoma Melanoma XXX
#Rscript trunc_SV_detection.R 1FROM7 Connectivetissuesofttissuesarcoma Bonesarcoma GIST Skinsarcoma Uterinesarcoma Melanoma Thyroid XXX
#Rscript trunc_SV_detection.R 1FROM8 Pituitary Parathyroid ACC GINET Lungcarcinoid Ovaryneuroendocrine Paraganglioma Phaeochromocytoma
#Rscript trunc_SV_detection.R 1FROM8 Breast ACC CNS Connectivetissuesofttissuesarcoma Bonesarcoma GIST Skinsarcoma Uterinesarcoma
#Rscript trunc_SV_detection.R 2FROM3 Breast Thyroid Endometrium XXX XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 2FROM4 Breast Aerodigestivetract Lung Ovary XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 2FROM4 Colorectal Breast Gastric Ovarysexcord-gonadalstromal XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 2FROM4 Colorectal Endometrium Ovary Sebaceous XXX XXX XXX XXX
#Rscript trunc_SV_detection.R 2FROM6 Uvealmelanoma Kidney Melanoma Lung Mesothelioma CNSmeningioma XXX XXX
#Rscript trunc_SV_detection.R 2FROM7 Connectivetissuesofttissuesarcoma Bonesarcoma GIST Skinsarcoma Uterinesarcoma Melanoma Thyroid XXX
#Rscript trunc_SV_detection.R 2FROM8 Breast ACC CNS Connectivetissuesofttissuesarcoma Bonesarcoma GIST Skinsarcoma Uterinesarcoma
#############################################################################################################################################################

args = commandArgs(trailingOnly=TRUE)

##Read in necessary files
elements <- read.csv("trunc_gene_list_regions_exons.csv", header = T)
elements <- elements[which(!(is.na(elements$start))),]

MPT_20170104 <- read.delim("MPT_probands_euro_20170104.txt", header = F)
BRIDGE_list <- read.delim("BRIDGE_euro_controls_20170104.txt", header = F)

##Deletions called by Manta
variants <- read.delim("mantacalls_allsamples_del.ann.filt.0.01_subset.txt",header = T) ##Filtered variant call file provided by NIHR BioResource Rare Disease project
variants <- variants[which(variants$MANTA_GQ > 29),]

##Adapting the input tables for use
elements$start <- as.numeric(elements$start)
elements$end <- as.numeric(elements$end)
elements$chrom <- as.character(elements$chrom)

variants$MANTA_CHROM <- as.character(variants$MANTA_CHROM)
variants$MANTA_START <- as.numeric(variants$MANTA_START)
variants$MANTA_END <- as.numeric(variants$MANTA_END)

library(stringr)
split_MANTA_CIPOS <- str_split_fixed(variants$MANTA_CIPOS, ",", 2)
split_MANTA_CIPOS <- data.frame(split_MANTA_CIPOS)
colnames(split_MANTA_CIPOS) <- c("MANTA_CIPOS_MIN_POS", "MANTA_CIPOS_MAX_POS")
split_MANTA_CIPOS$MANTA_CIPOS_MIN_POS <- as.numeric(as.character(split_MANTA_CIPOS$MANTA_CIPOS_MIN_POS))
split_MANTA_CIPOS$MANTA_CIPOS_MAX_POS <- as.numeric(as.character(split_MANTA_CIPOS$MANTA_CIPOS_MAX_POS))
split_MANTA_CIPOS <- replace(split_MANTA_CIPOS, is.na(split_MANTA_CIPOS), 0)
split_MANTA_CIEND <- str_split_fixed(variants$MANTA_CIEND, ",", 2)
split_MANTA_CIEND <- data.frame(split_MANTA_CIEND)
colnames(split_MANTA_CIEND) <- c("MANTA_CIEND_MIN_END", "MANTA_CIEND_MAX_END")
split_MANTA_CIEND$MANTA_CIEND_MIN_END <- as.numeric(as.character(split_MANTA_CIEND$MANTA_CIEND_MIN_END))
split_MANTA_CIEND$MANTA_CIEND_MAX_END <- as.numeric(as.character(split_MANTA_CIEND$MANTA_CIEND_MAX_END))
split_MANTA_CIEND <- replace(split_MANTA_CIEND, is.na(split_MANTA_CIEND), 0)
MANTA_MIN_START <- variants$MANTA_START + split_MANTA_CIPOS$MANTA_CIPOS_MIN_POS
MANTA_MAX_START <-variants$MANTA_START + split_MANTA_CIPOS$MANTA_CIPOS_MAX_POS
MANTA_MIN_END <- variants$MANTA_END + split_MANTA_CIEND$MANTA_CIEND_MIN_END
MANTA_MAX_END <- variants$MANTA_END + split_MANTA_CIEND$MANTA_CIEND_MAX_END

variants <- cbind(variants,MANTA_MIN_START,MANTA_MAX_START,split_MANTA_CIPOS,MANTA_MIN_END,MANTA_MAX_END,split_MANTA_CIEND)

##Match variants to elements
library(sqldf)

##Deletion of entire gene
variant_element_match_full_del <- sqldf("select * from variants f1 
                               left join elements f2 on 
                               (
(f1.MANTA_CHROM==f2.chrom) and 

(f1.MANTA_MAX_START < f2.start and f1.MANTA_MIN_END > f2.end) 
)")

mut_type_col <- vector(mode = "character", nrow(variant_element_match_full_del))
mut_type_col[which(!is.na(variant_element_match_full_del$identifier))] <- "full_del"
variant_element_match_full_del <- cbind(variant_element_match_full_del, mut_type_col)

##Deletion within gene
variant_element_match_intra_del <- sqldf("select * from variants f1 
                               left join elements f2 on 
                               (
(f1.MANTA_CHROM==f2.chrom) and 

(f1.MANTA_MIN_START > f2.start and f1.MANTA_MAX_END < f2.end) 
)")

mut_type_col <- vector(mode = "character", nrow(variant_element_match_intra_del))
mut_type_col[which(!is.na(variant_element_match_intra_del$identifier))] <- "intra_del"
variant_element_match_intra_del <- cbind(variant_element_match_intra_del, mut_type_col)

##Deletion involving one end of gene but not other
variant_element_match_start_end_del <- sqldf("select * from variants f1 
                               left join elements f2 on 
                               (
(f1.MANTA_CHROM==f2.chrom) and 

(

((f1.MANTA_MAX_START < f2.start) and (f1.MANTA_MIN_END > f2.start and f1.MANTA_MAX_END < f2.end)) or

((f1.MANTA_MIN_END > f2.end) and (f1.MANTA_MAX_START < f2.end and f1.MANTA_MIN_START > f2.start))

)

)")

mut_type_col <- vector(mode = "character", nrow(variant_element_match_start_end_del))
mut_type_col[which(!is.na(variant_element_match_start_end_del$identifier))] <- "start_end_del"
variant_element_match_start_end_del <- cbind(variant_element_match_start_end_del, mut_type_col)

##Collate results
variant_element_match_manta_del <- 
rbind(
variant_element_match_full_del[which(!is.na(variant_element_match_full_del$identifier)),],
variant_element_match_intra_del[which(!is.na(variant_element_match_intra_del$identifier)),],
variant_element_match_start_end_del[which(!is.na(variant_element_match_start_end_del$identifier)),]
)

##Translocations called by Manta
variants <-  read.delim("mantacalls_allsamples_bnd.ann.filt.0.01_subset.txt",header = T) ##Filtered variant call file provided by NIHR BioResource Rare Disease project
variants <- variants[which(variants$MANTA_GQ > 29),]
elements <- read.csv("trunc_gene_list_regions_genes.csv", header = T) 

##Adapting the input tables for use
elements$start <- as.numeric(elements$start)
elements$end <- as.numeric(elements$end)
elements$chrom <- as.character(elements$chrom)

variants$MANTA_CHROM <- as.character(variants$MANTA_CHROM)
variants$MANTA_START <- as.numeric(variants$MANTA_START)
variants$MANTA_END <- as.numeric(variants$MANTA_END)

library(stringr)
split_MANTA_CIPOS <- str_split_fixed(variants$MANTA_CIPOS, ",", 2)
split_MANTA_CIPOS <- data.frame(split_MANTA_CIPOS)
colnames(split_MANTA_CIPOS) <- c("MANTA_CIPOS_MIN_POS", "MANTA_CIPOS_MAX_POS")
split_MANTA_CIPOS$MANTA_CIPOS_MIN_POS <- as.numeric(as.character(split_MANTA_CIPOS$MANTA_CIPOS_MIN_POS))
split_MANTA_CIPOS$MANTA_CIPOS_MAX_POS <- as.numeric(as.character(split_MANTA_CIPOS$MANTA_CIPOS_MAX_POS))
split_MANTA_CIPOS <- replace(split_MANTA_CIPOS, is.na(split_MANTA_CIPOS), 0)
split_MANTA_CIEND <- str_split_fixed(variants$MANTA_CIEND, ",", 2)
split_MANTA_CIEND <- data.frame(split_MANTA_CIEND)
colnames(split_MANTA_CIEND) <- c("MANTA_CIEND_MIN_END", "MANTA_CIEND_MAX_END")
split_MANTA_CIEND$MANTA_CIEND_MIN_END <- as.numeric(as.character(split_MANTA_CIEND$MANTA_CIEND_MIN_END))
split_MANTA_CIEND$MANTA_CIEND_MAX_END <- as.numeric(as.character(split_MANTA_CIEND$MANTA_CIEND_MAX_END))
split_MANTA_CIEND <- replace(split_MANTA_CIEND, is.na(split_MANTA_CIEND), 0)
MANTA_MIN_START <- variants$MANTA_START + split_MANTA_CIPOS$MANTA_CIPOS_MIN_POS
MANTA_MAX_START <-variants$MANTA_START + split_MANTA_CIPOS$MANTA_CIPOS_MAX_POS
MANTA_MIN_END <- variants$MANTA_END + split_MANTA_CIEND$MANTA_CIEND_MIN_END
MANTA_MAX_END <- variants$MANTA_END + split_MANTA_CIEND$MANTA_CIEND_MAX_END

variants <- cbind(variants,MANTA_MIN_START,MANTA_MAX_START,split_MANTA_CIPOS,MANTA_MIN_END,MANTA_MAX_END,split_MANTA_CIEND)

##Matching variants to elements
library(sqldf)

##Translocation with breakpoint within gene
variant_element_match_bnd <- sqldf("select * from variants f1 
                               left join elements f2 on 
                               (
(f1.MANTA_CHROM==f2.chrom) and 

(
(f1.MANTA_MIN_START > f2.start and f1.MANTA_MAX_START < f2.end) or (f1.MANTA_MIN_END > f2.start and f1.MANTA_MAX_END < f2.end)
)

)")

##Collate results
mut_type_col <- vector(mode = "character", nrow(variant_element_match_bnd))
mut_type_col[which(!is.na(variant_element_match_bnd$identifier))] <- "bnd"
variant_element_match_bnd <- cbind(variant_element_match_bnd, mut_type_col)
variant_element_match_bnd <- variant_element_match_bnd[which(!is.na(variant_element_match_bnd$identifier)),]

##Inversions called by Manta
variants <-  read.delim("mantacalls_allsamples_inv.ann.filt.0.01_subset.txt", header = T) ##Filtered variant call file provided by NIHR BioResource Rare Disease project
variants <- variants[which(variants$MANTA_GQ > 29),]
elements <- read.csv("trunc_gene_list_regions_genes.csv", header = T) 

##Adapting the input tables for use
elements$start <- as.numeric(elements$start)
elements$end <- as.numeric(elements$end)
elements$chrom <- as.character(elements$chrom)

variants$MANTA_CHROM <- as.character(variants$MANTA_CHROM)
variants$MANTA_START <- as.numeric(variants$MANTA_START)
variants$MANTA_END <- as.numeric(variants$MANTA_END)

library(stringr)
split_MANTA_CIPOS <- str_split_fixed(variants$MANTA_CIPOS, ",", 2)
split_MANTA_CIPOS <- data.frame(split_MANTA_CIPOS)
colnames(split_MANTA_CIPOS) <- c("MANTA_CIPOS_MIN_POS", "MANTA_CIPOS_MAX_POS")
split_MANTA_CIPOS$MANTA_CIPOS_MIN_POS <- as.numeric(as.character(split_MANTA_CIPOS$MANTA_CIPOS_MIN_POS))
split_MANTA_CIPOS$MANTA_CIPOS_MAX_POS <- as.numeric(as.character(split_MANTA_CIPOS$MANTA_CIPOS_MAX_POS))
split_MANTA_CIPOS <- replace(split_MANTA_CIPOS, is.na(split_MANTA_CIPOS), 0)
split_MANTA_CIEND <- str_split_fixed(variants$MANTA_CIEND, ",", 2)
split_MANTA_CIEND <- data.frame(split_MANTA_CIEND)
colnames(split_MANTA_CIEND) <- c("MANTA_CIEND_MIN_END", "MANTA_CIEND_MAX_END")
split_MANTA_CIEND$MANTA_CIEND_MIN_END <- as.numeric(as.character(split_MANTA_CIEND$MANTA_CIEND_MIN_END))
split_MANTA_CIEND$MANTA_CIEND_MAX_END <- as.numeric(as.character(split_MANTA_CIEND$MANTA_CIEND_MAX_END))
split_MANTA_CIEND <- replace(split_MANTA_CIEND, is.na(split_MANTA_CIEND), 0)
MANTA_MIN_START <- variants$MANTA_START + split_MANTA_CIPOS$MANTA_CIPOS_MIN_POS
MANTA_MAX_START <-variants$MANTA_START + split_MANTA_CIPOS$MANTA_CIPOS_MAX_POS
MANTA_MIN_END <- variants$MANTA_END + split_MANTA_CIEND$MANTA_CIEND_MIN_END
MANTA_MAX_END <- variants$MANTA_END + split_MANTA_CIEND$MANTA_CIEND_MAX_END

variants <- cbind(variants,MANTA_MIN_START,MANTA_MAX_START,split_MANTA_CIPOS,MANTA_MIN_END,MANTA_MAX_END,split_MANTA_CIEND)


##Matching variants to elements
library(sqldf)

##Inversion with breakpoint within gene
variant_element_match_inv <- sqldf("select * from variants f1 
                               left join elements f2 on 
                               (
(f1.MANTA_CHROM==f2.chrom) and 

(
(f1.MANTA_MIN_START > f2.start and f1.MANTA_MAX_START < f2.end) or (f1.MANTA_MIN_END > f2.start and f1.MANTA_MAX_END < f2.end)
)

)")

##Collate results
mut_type_col <- vector(mode = "character", nrow(variant_element_match_inv))
mut_type_col[which(!is.na(variant_element_match_inv$identifier))] <- "inv"
variant_element_match_inv <- cbind(variant_element_match_inv, mut_type_col)
variant_element_match_inv <- variant_element_match_inv[which(!is.na(variant_element_match_inv$identifier)),]

##Insertions called by Manta
variants <-  read.delim("mantacalls_allsamples_ins.ann.filt.0.01_subset.txt",header = T) ##Filtered variant call file provided by NIHR BioResource Rare Disease project
variants <- variants[which(variants$MANTA_GQ > 29),]
elements <- read.csv("trunc_gene_list_regions_genes.csv", header = T) 

##Adapting the input tables for use
elements$start <- as.numeric(elements$start)
elements$end <- as.numeric(elements$end)
elements$chrom <- as.character(elements$chrom)

variants$MANTA_CHROM <- as.character(variants$MANTA_CHROM)
variants$MANTA_START <- as.numeric(variants$MANTA_START)
variants$MANTA_END <- as.numeric(variants$MANTA_END)

library(stringr)
split_MANTA_CIPOS <- str_split_fixed(variants$MANTA_CIPOS, ",", 2)
split_MANTA_CIPOS <- data.frame(split_MANTA_CIPOS)
colnames(split_MANTA_CIPOS) <- c("MANTA_CIPOS_MIN_POS", "MANTA_CIPOS_MAX_POS")
split_MANTA_CIPOS$MANTA_CIPOS_MIN_POS <- as.numeric(as.character(split_MANTA_CIPOS$MANTA_CIPOS_MIN_POS))
split_MANTA_CIPOS$MANTA_CIPOS_MAX_POS <- as.numeric(as.character(split_MANTA_CIPOS$MANTA_CIPOS_MAX_POS))
split_MANTA_CIPOS <- replace(split_MANTA_CIPOS, is.na(split_MANTA_CIPOS), 0)
split_MANTA_CIEND <- str_split_fixed(variants$MANTA_CIEND, ",", 2)
split_MANTA_CIEND <- data.frame(split_MANTA_CIEND)
colnames(split_MANTA_CIEND) <- c("MANTA_CIEND_MIN_END", "MANTA_CIEND_MAX_END")
split_MANTA_CIEND$MANTA_CIEND_MIN_END <- as.numeric(as.character(split_MANTA_CIEND$MANTA_CIEND_MIN_END))
split_MANTA_CIEND$MANTA_CIEND_MAX_END <- as.numeric(as.character(split_MANTA_CIEND$MANTA_CIEND_MAX_END))
split_MANTA_CIEND <- replace(split_MANTA_CIEND, is.na(split_MANTA_CIEND), 0)
MANTA_MIN_START <- variants$MANTA_START + split_MANTA_CIPOS$MANTA_CIPOS_MIN_POS
MANTA_MAX_START <-variants$MANTA_START + split_MANTA_CIPOS$MANTA_CIPOS_MAX_POS
MANTA_MIN_END <- variants$MANTA_END + split_MANTA_CIEND$MANTA_CIEND_MIN_END
MANTA_MAX_END <- variants$MANTA_END + split_MANTA_CIEND$MANTA_CIEND_MAX_END

variants <- cbind(variants,MANTA_MIN_START,MANTA_MAX_START,split_MANTA_CIPOS,MANTA_MIN_END,MANTA_MAX_END,split_MANTA_CIEND)


##Matching variants to elements
library(sqldf)

##Insertion with breakpoint within gene
variant_element_match_ins <- sqldf("select * from variants f1 
                               left join elements f2 on 
                               (
(f1.MANTA_CHROM==f2.chrom) and 

(
(f1.MANTA_MIN_START > f2.start and f1.MANTA_MAX_START < f2.end) or (f1.MANTA_MIN_END > f2.start and f1.MANTA_MAX_END < f2.end)
)

)")

##Collate results
mut_type_col <- vector(mode = "character", nrow(variant_element_match_ins))
mut_type_col[which(!is.na(variant_element_match_ins$identifier))] <- "ins"
variant_element_match_ins <- cbind(variant_element_match_ins, mut_type_col)
variant_element_match_ins <- variant_element_match_ins[which(!is.na(variant_element_match_ins$identifier)),]

##Copy number loss called by Canvas
variants <-  read.delim("canvascalls_allsamples_loss.ann.filt.0.01_subset.txt",header = T) ##Filtered variant call file provided by NIHR BioResource Rare Disease project
variants <- variants[which(variants$CANVAS_QUAL > 29),]
elements <- read.csv("trunc_gene_list_regions_exons.csv", header = T) 

##Adapting the input tables for use
elements$start <- as.numeric(elements$start)
elements$end <- as.numeric(elements$end)
elements$chrom <- as.character(elements$chrom)

variants$CANVAS_CHROM <- as.character(variants$CANVAS_CHROM)
variants$CANVAS_START <- as.numeric(variants$CANVAS_START)
variants$CANVAS_END <- as.numeric(variants$CANVAS_END)

##Matching variants to elements
library(sqldf)

##Deletion of entire gene
variant_element_match_full_del <- sqldf("select * from variants f1 
                                        left join elements f2 on 
                                        (
                                        (f1.CANVAS_CHROM==f2.chrom) and 
                                        
                                        (f1.CANVAS_START < f2.start and f1.CANVAS_END > f2.end) 
                                        )")

mut_type_col <- vector(mode = "character", nrow(variant_element_match_full_del))
mut_type_col[which(!is.na(variant_element_match_full_del$identifier))] <- "full_del"
variant_element_match_full_del <- cbind(variant_element_match_full_del, mut_type_col)

##Deletion within gene
variant_element_match_intra_del <- sqldf("select * from variants f1 
                                         left join elements f2 on 
                                         (
                                         (f1.CANVAS_CHROM==f2.chrom) and 
                                         
                                         (f1.CANVAS_START > f2.start and f1.CANVAS_END < f2.end) 
                                         )")

mut_type_col <- vector(mode = "character", nrow(variant_element_match_intra_del))
mut_type_col[which(!is.na(variant_element_match_intra_del$identifier))] <- "intra_del"
variant_element_match_intra_del <- cbind(variant_element_match_intra_del, mut_type_col)

##Deletion involving one end of gene but not other
variant_element_match_start_end_del <- sqldf("select * from variants f1 
                                             left join elements f2 on 
                                             (
                                             (f1.CANVAS_CHROM==f2.chrom) and 
                                             
                                             (
                                             
                                             ((f1.CANVAS_START < f2.start) and (f1.CANVAS_END > f2.start and f1.CANVAS_END < f2.end)) or
                                             
                                             ((f1.CANVAS_END > f2.end) and (f1.CANVAS_START < f2.end and f1.CANVAS_START > f2.start))
                                             
                                             )
                                             
                                             )")

mut_type_col <- vector(mode = "character", nrow(variant_element_match_start_end_del))
mut_type_col[which(!is.na(variant_element_match_start_end_del$identifier))] <- "start_end_del"
variant_element_match_start_end_del <- cbind(variant_element_match_start_end_del, mut_type_col)

##Collate results
variant_element_match_canvas_del <- 
rbind(
variant_element_match_full_del[which(!is.na(variant_element_match_full_del$identifier)),],
variant_element_match_intra_del[which(!is.na(variant_element_match_intra_del$identifier)),],
variant_element_match_start_end_del[which(!is.na(variant_element_match_start_end_del$identifier)),]
)

##Merge annotated variant tables with different SV modalities

library(data.table)
SV_table <- rbindlist(list(variant_element_match_manta_del, variant_element_match_bnd, variant_element_match_inv, variant_element_match_ins, variant_element_match_canvas_del), fill=TRUE)

##Remove lines with no variant element match
SV_table <- SV_table[which(!is.na(SV_table$mut_type_col)),]

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
if(args[1] != "ALL")case_list <- data.frame(case_list[which(case_list$indv %in% MPT_20170104$V1),])   
if(args[1] != "ALL")colnames(case_list) <- "V1"
if(args[1] == "ALL")case_list <- read.delim("MPT_probands_euro_20170104.txt", header = F)

#######################################################################
##Counts of variants and comparison of frequency in cases vs controls##
#######################################################################

##Read in gene lists
full_gene_list <- read.delim("trunc_gene_list.txt", header = F)
webgestalt_gene_list <- read.delim("trunc_gene_list_webgestalt.txt", header = F)
loftool_gene_list <- read.delim("trunc_gene_list_loftool.txt", header = F)
cgp_gene_list <- read.delim("trunc_gene_list_cgp.txt", header = F)
mania_gene_list <- read.delim("trunc_gene_list_mania.txt", header = F)
repair_gene_list <- read.delim("trunc_gene_list_repair.txt", header = F)


#############################
##Full gene list#############
#############################

##Reduce table to relevant gene list
SV_table_full <- SV_table[which(SV_table$identifier %in% full_gene_list$V1),]


##Subset table into case and control group
SV_table_full_cases <- SV_table_full[which(SV_table_full$BRIDGE_ID %in% case_list$V1),]
SV_table_full_controls <- SV_table_full[which(SV_table_full$BRIDGE_ID %in% BRIDGE_list$V1),]

##Count variants affecting each gene and put into table
cases_count <- lapply(full_gene_list$V1, function(x) length(unique(SV_table_full_cases$BRIDGE_ID[which(SV_table_full_cases$identifier == x)])))
controls_count <- lapply(full_gene_list$V1, function(x) length(unique(SV_table_full_controls$BRIDGE_ID[which(SV_table_full_controls$identifier == x)])))

cases_count_hets <- lapply(full_gene_list$V1, function(x) length(unique(SV_table_full_cases$BRIDGE_ID[which(SV_table_full_cases$identifier == x & (SV_table_full_cases$MANTA_GT == "0/1" | SV_table_full_cases$CANVAS_GT == "0/1"))])))
controls_count_hets <- lapply(full_gene_list$V1, function(x) length(unique(SV_table_full_controls$BRIDGE_ID[which(SV_table_full_controls$identifier == x & (SV_table_full_controls$MANTA_GT == "0/1" | SV_table_full_controls$CANVAS_GT == "0/1"))])))

cases_count_homs <- lapply(full_gene_list$V1, function(x) length(unique(SV_table_full_cases$BRIDGE_ID[which(SV_table_full_cases$identifier == x & (SV_table_full_cases$MANTA_GT == "1/1" | SV_table_full_cases$CANVAS_GT == "1/1"))])))
controls_count_homs <- lapply(full_gene_list$V1, function(x) length(unique(SV_table_full_controls$BRIDGE_ID[which(SV_table_full_controls$identifier == x & (SV_table_full_controls$MANTA_GT == "1/1" | SV_table_full_controls$CANVAS_GT == "1/1"))])))

cases_controls_count <- data.frame(full_gene_list, unlist(cases_count),unlist(controls_count),unlist(cases_count_hets),unlist(controls_count_hets),unlist(cases_count_homs),unlist(controls_count_homs))
colnames(cases_controls_count) <- c("gene", "cases_count", "controls_count","cases_count_hets", "controls_count_hets","cases_count_homs", "controls_count_homs")

##Hypothesis testing and compiling into table
p_val_hets_homs <- unlist(lapply(lapply(cases_controls_count$gene, function(x)fisher.test(matrix(c(cases_controls_count$cases_count[x],nrow(case_list) - cases_controls_count$cases_count[x],cases_controls_count$controls_count[x],nrow(BRIDGE_list) - cases_controls_count$controls_count[x]),nrow = 2))),"[[", 1))
q_val_hets_homs <- p.adjust(p_val_hets_homs, method = "fdr", n = nrow(full_gene_list))
cases_prop_hets_homs <- as.numeric(cases_count) / nrow(case_list)
controls_prop_hets_homs <- as.numeric(controls_count) / nrow(BRIDGE_list)

p_val_hets <- unlist(lapply(lapply(cases_controls_count$gene, function(x)fisher.test(matrix(c(cases_controls_count$cases_count_hets[x],nrow(case_list) - cases_controls_count$cases_count_hets[x],cases_controls_count$controls_count_hets[x],nrow(BRIDGE_list) - cases_controls_count$controls_count_hets[x]),nrow = 2))),"[[", 1))
q_val_hets <- p.adjust(p_val_hets, method = "fdr", n = nrow(full_gene_list))
cases_prop_hets <- as.numeric(cases_count_hets) / nrow(case_list)
controls_prop_hets <- as.numeric(controls_count_hets) / nrow(BRIDGE_list)

p_val_homs <- unlist(lapply(lapply(cases_controls_count$gene, function(x)fisher.test(matrix(c(cases_controls_count$cases_count_homs[x],nrow(case_list) - cases_controls_count$cases_count_homs[x],cases_controls_count$controls_count_homs[x],nrow(BRIDGE_list) - cases_controls_count$controls_count_homs[x]),nrow = 2))),"[[", 1))
q_val_homs <- p.adjust(p_val_homs, method = "fdr", n = nrow(full_gene_list))
cases_prop_homs <- as.numeric(cases_count_homs) / nrow(case_list)
controls_prop_homs <- as.numeric(controls_count_homs) / nrow(BRIDGE_list)

gene_list_col <- rep("full", nrow(full_gene_list))
cases_controls_count_full <- cbind(cases_controls_count,cases_prop_hets,controls_prop_hets,cases_prop_homs,controls_prop_homs,cases_prop_hets_homs,controls_prop_hets_homs,p_val_hets, q_val_hets,p_val_homs,q_val_homs,p_val_hets_homs,q_val_hets_homs,gene_list_col)

########################
##Webgestalt gene list##
########################

##Reduce table to relevant gene list
SV_table_webgestalt <- SV_table[which(SV_table$identifier %in% webgestalt_gene_list$V1),]

##Subset table into case and control group
SV_table_webgestalt_cases <- SV_table_webgestalt[which(SV_table_webgestalt$BRIDGE_ID %in% case_list$V1),]
SV_table_webgestalt_controls <- SV_table_webgestalt[which(SV_table_webgestalt$BRIDGE_ID %in% BRIDGE_list$V1),]

##Count variants affecting each gene and put into table
cases_count <- lapply(webgestalt_gene_list$V1, function(x) length(unique(SV_table_webgestalt_cases$BRIDGE_ID[which(SV_table_webgestalt_cases$identifier == x)])))
controls_count <- lapply(webgestalt_gene_list$V1, function(x) length(unique(SV_table_webgestalt_controls$BRIDGE_ID[which(SV_table_webgestalt_controls$identifier == x)])))

cases_count_hets <- lapply(webgestalt_gene_list$V1, function(x) length(unique(SV_table_webgestalt_cases$BRIDGE_ID[which(SV_table_webgestalt_cases$identifier == x & (SV_table_webgestalt_cases$MANTA_GT == "0/1" | SV_table_webgestalt_cases$CANVAS_GT == "0/1"))])))
controls_count_hets <- lapply(webgestalt_gene_list$V1, function(x) length(unique(SV_table_webgestalt_controls$BRIDGE_ID[which(SV_table_webgestalt_controls$identifier == x & (SV_table_webgestalt_controls$MANTA_GT == "0/1" | SV_table_webgestalt_controls$CANVAS_GT == "0/1"))])))

cases_count_homs <- lapply(webgestalt_gene_list$V1, function(x) length(unique(SV_table_webgestalt_cases$BRIDGE_ID[which(SV_table_webgestalt_cases$identifier == x & (SV_table_webgestalt_cases$MANTA_GT == "1/1" | SV_table_webgestalt_cases$CANVAS_GT == "1/1"))])))
controls_count_homs <- lapply(webgestalt_gene_list$V1, function(x) length(unique(SV_table_webgestalt_controls$BRIDGE_ID[which(SV_table_webgestalt_controls$identifier == x & (SV_table_webgestalt_controls$MANTA_GT == "1/1" | SV_table_webgestalt_controls$CANVAS_GT == "1/1"))])))

cases_controls_count <- data.frame(webgestalt_gene_list, unlist(cases_count),unlist(controls_count),unlist(cases_count_hets),unlist(controls_count_hets),unlist(cases_count_homs),unlist(controls_count_homs))
colnames(cases_controls_count) <- c("gene", "cases_count", "controls_count","cases_count_hets", "controls_count_hets","cases_count_homs", "controls_count_homs")

##Hypothesis testing and compiling into table
p_val_hets_homs <- unlist(lapply(lapply(cases_controls_count$gene, function(x)fisher.test(matrix(c(cases_controls_count$cases_count[x],nrow(case_list) - cases_controls_count$cases_count[x],cases_controls_count$controls_count[x],nrow(BRIDGE_list) - cases_controls_count$controls_count[x]),nrow = 2))),"[[", 1))
q_val_hets_homs <- p.adjust(p_val_hets_homs, method = "fdr", n = nrow(webgestalt_gene_list))
cases_prop_hets_homs <- as.numeric(cases_count) / nrow(case_list)
controls_prop_hets_homs <- as.numeric(controls_count) / nrow(BRIDGE_list)

p_val_hets <- unlist(lapply(lapply(cases_controls_count$gene, function(x)fisher.test(matrix(c(cases_controls_count$cases_count_hets[x],nrow(case_list) - cases_controls_count$cases_count_hets[x],cases_controls_count$controls_count_hets[x],nrow(BRIDGE_list) - cases_controls_count$controls_count_hets[x]),nrow = 2))),"[[", 1))
q_val_hets <- p.adjust(p_val_hets, method = "fdr", n = nrow(webgestalt_gene_list))
cases_prop_hets <- as.numeric(cases_count_hets) / nrow(case_list)
controls_prop_hets <- as.numeric(controls_count_hets) / nrow(BRIDGE_list)

p_val_homs <- unlist(lapply(lapply(cases_controls_count$gene, function(x)fisher.test(matrix(c(cases_controls_count$cases_count_homs[x],nrow(case_list) - cases_controls_count$cases_count_homs[x],cases_controls_count$controls_count_homs[x],nrow(BRIDGE_list) - cases_controls_count$controls_count_homs[x]),nrow = 2))),"[[", 1))
q_val_homs <- p.adjust(p_val_homs, method = "fdr", n = nrow(webgestalt_gene_list))
cases_prop_homs <- as.numeric(cases_count_homs) / nrow(case_list)
controls_prop_homs <- as.numeric(controls_count_homs) / nrow(BRIDGE_list)

gene_list_col <- rep("webgestalt", nrow(webgestalt_gene_list))
cases_controls_count_webgestalt <- cbind(cases_controls_count,cases_prop_hets,controls_prop_hets,cases_prop_homs,controls_prop_homs,cases_prop_hets_homs,controls_prop_hets_homs,p_val_hets, q_val_hets,p_val_homs,q_val_homs,p_val_hets_homs,q_val_hets_homs,gene_list_col)

#####################
##Loftool gene list##
#####################

##Reduce table to relevant gene list
SV_table_loftool <- SV_table[which(SV_table$identifier %in% loftool_gene_list$V1),]

##Subset table into case and control group
SV_table_loftool_cases <- SV_table_loftool[which(SV_table_loftool$BRIDGE_ID %in% case_list$V1),]
SV_table_loftool_controls <- SV_table_loftool[which(SV_table_loftool$BRIDGE_ID %in% BRIDGE_list$V1),]

##Count variants affecting each gene and put into table
cases_count <- lapply(loftool_gene_list$V1, function(x) length(unique(SV_table_loftool_cases$BRIDGE_ID[which(SV_table_loftool_cases$identifier == x)])))
controls_count <- lapply(loftool_gene_list$V1, function(x) length(unique(SV_table_loftool_controls$BRIDGE_ID[which(SV_table_loftool_controls$identifier == x)])))

cases_count_hets <- lapply(loftool_gene_list$V1, function(x) length(unique(SV_table_loftool_cases$BRIDGE_ID[which(SV_table_loftool_cases$identifier == x & (SV_table_loftool_cases$MANTA_GT == "0/1" | SV_table_loftool_cases$CANVAS_GT == "0/1"))])))
controls_count_hets <- lapply(loftool_gene_list$V1, function(x) length(unique(SV_table_loftool_controls$BRIDGE_ID[which(SV_table_loftool_controls$identifier == x & (SV_table_loftool_controls$MANTA_GT == "0/1" | SV_table_loftool_controls$CANVAS_GT == "0/1"))])))

cases_count_homs <- lapply(loftool_gene_list$V1, function(x) length(unique(SV_table_loftool_cases$BRIDGE_ID[which(SV_table_loftool_cases$identifier == x & (SV_table_loftool_cases$MANTA_GT == "1/1" | SV_table_loftool_cases$CANVAS_GT == "1/1"))])))
controls_count_homs <- lapply(loftool_gene_list$V1, function(x) length(unique(SV_table_loftool_controls$BRIDGE_ID[which(SV_table_loftool_controls$identifier == x & (SV_table_loftool_controls$MANTA_GT == "1/1" | SV_table_loftool_controls$CANVAS_GT == "1/1"))])))

cases_controls_count <- data.frame(loftool_gene_list, unlist(cases_count),unlist(controls_count),unlist(cases_count_hets),unlist(controls_count_hets),unlist(cases_count_homs),unlist(controls_count_homs))
colnames(cases_controls_count) <- c("gene", "cases_count", "controls_count","cases_count_hets", "controls_count_hets","cases_count_homs", "controls_count_homs")

##Hypothesis testing and compiling into table
p_val_hets_homs <- unlist(lapply(lapply(cases_controls_count$gene, function(x)fisher.test(matrix(c(cases_controls_count$cases_count[x],nrow(case_list) - cases_controls_count$cases_count[x],cases_controls_count$controls_count[x],nrow(BRIDGE_list) - cases_controls_count$controls_count[x]),nrow = 2))),"[[", 1))
q_val_hets_homs <- p.adjust(p_val_hets_homs, method = "fdr", n = nrow(loftool_gene_list))
cases_prop_hets_homs <- as.numeric(cases_count) / nrow(case_list)
controls_prop_hets_homs <- as.numeric(controls_count) / nrow(BRIDGE_list)

p_val_hets <- unlist(lapply(lapply(cases_controls_count$gene, function(x)fisher.test(matrix(c(cases_controls_count$cases_count_hets[x],nrow(case_list) - cases_controls_count$cases_count_hets[x],cases_controls_count$controls_count_hets[x],nrow(BRIDGE_list) - cases_controls_count$controls_count_hets[x]),nrow = 2))),"[[", 1))
q_val_hets <- p.adjust(p_val_hets, method = "fdr", n = nrow(loftool_gene_list))
cases_prop_hets <- as.numeric(cases_count_hets) / nrow(case_list)
controls_prop_hets <- as.numeric(controls_count_hets) / nrow(BRIDGE_list)

p_val_homs <- unlist(lapply(lapply(cases_controls_count$gene, function(x)fisher.test(matrix(c(cases_controls_count$cases_count_homs[x],nrow(case_list) - cases_controls_count$cases_count_homs[x],cases_controls_count$controls_count_homs[x],nrow(BRIDGE_list) - cases_controls_count$controls_count_homs[x]),nrow = 2))),"[[", 1))
q_val_homs <- p.adjust(p_val_homs, method = "fdr", n = nrow(loftool_gene_list))
cases_prop_homs <- as.numeric(cases_count_homs) / nrow(case_list)
controls_prop_homs <- as.numeric(controls_count_homs) / nrow(BRIDGE_list)

gene_list_col <- rep("loftool", nrow(loftool_gene_list))
cases_controls_count_loftool <- cbind(cases_controls_count,cases_prop_hets,controls_prop_hets,cases_prop_homs,controls_prop_homs,cases_prop_hets_homs,controls_prop_hets_homs,p_val_hets, q_val_hets,p_val_homs,q_val_homs,p_val_hets_homs,q_val_hets_homs,gene_list_col)

#################
##Cgp gene list##
#################

##Reduce table to relevant gene list
SV_table_cgp <- SV_table[which(SV_table$identifier %in% cgp_gene_list$V1),]

##Subset table into case and control group
SV_table_cgp_cases <- SV_table_cgp[which(SV_table_cgp$BRIDGE_ID %in% case_list$V1),]
SV_table_cgp_controls <- SV_table_cgp[which(SV_table_cgp$BRIDGE_ID %in% BRIDGE_list$V1),]

##Count variants affecting each gene and put into table
cases_count <- lapply(cgp_gene_list$V1, function(x) length(unique(SV_table_cgp_cases$BRIDGE_ID[which(SV_table_cgp_cases$identifier == x)])))
controls_count <- lapply(cgp_gene_list$V1, function(x) length(unique(SV_table_cgp_controls$BRIDGE_ID[which(SV_table_cgp_controls$identifier == x)])))

cases_count_hets <- lapply(cgp_gene_list$V1, function(x) length(unique(SV_table_cgp_cases$BRIDGE_ID[which(SV_table_cgp_cases$identifier == x & (SV_table_cgp_cases$MANTA_GT == "0/1" | SV_table_cgp_cases$CANVAS_GT == "0/1"))])))
controls_count_hets <- lapply(cgp_gene_list$V1, function(x) length(unique(SV_table_cgp_controls$BRIDGE_ID[which(SV_table_cgp_controls$identifier == x & (SV_table_cgp_controls$MANTA_GT == "0/1" | SV_table_cgp_controls$CANVAS_GT == "0/1"))])))

cases_count_homs <- lapply(cgp_gene_list$V1, function(x) length(unique(SV_table_cgp_cases$BRIDGE_ID[which(SV_table_cgp_cases$identifier == x & (SV_table_cgp_cases$MANTA_GT == "1/1" | SV_table_cgp_cases$CANVAS_GT == "1/1"))])))
controls_count_homs <- lapply(cgp_gene_list$V1, function(x) length(unique(SV_table_cgp_controls$BRIDGE_ID[which(SV_table_cgp_controls$identifier == x & (SV_table_cgp_controls$MANTA_GT == "1/1" | SV_table_cgp_controls$CANVAS_GT == "1/1"))])))

cases_controls_count <- data.frame(cgp_gene_list, unlist(cases_count),unlist(controls_count),unlist(cases_count_hets),unlist(controls_count_hets),unlist(cases_count_homs),unlist(controls_count_homs))
colnames(cases_controls_count) <- c("gene", "cases_count", "controls_count","cases_count_hets", "controls_count_hets","cases_count_homs", "controls_count_homs")

##Hypothesis testing and compiling into table
p_val_hets_homs <- unlist(lapply(lapply(cases_controls_count$gene, function(x)fisher.test(matrix(c(cases_controls_count$cases_count[x],nrow(case_list) - cases_controls_count$cases_count[x],cases_controls_count$controls_count[x],nrow(BRIDGE_list) - cases_controls_count$controls_count[x]),nrow = 2))),"[[", 1))
q_val_hets_homs <- p.adjust(p_val_hets_homs, method = "fdr", n = nrow(cgp_gene_list))
cases_prop_hets_homs <- as.numeric(cases_count) / nrow(case_list)
controls_prop_hets_homs <- as.numeric(controls_count) / nrow(BRIDGE_list)

p_val_hets <- unlist(lapply(lapply(cases_controls_count$gene, function(x)fisher.test(matrix(c(cases_controls_count$cases_count_hets[x],nrow(case_list) - cases_controls_count$cases_count_hets[x],cases_controls_count$controls_count_hets[x],nrow(BRIDGE_list) - cases_controls_count$controls_count_hets[x]),nrow = 2))),"[[", 1))
q_val_hets <- p.adjust(p_val_hets, method = "fdr", n = nrow(cgp_gene_list))
cases_prop_hets <- as.numeric(cases_count_hets) / nrow(case_list)
controls_prop_hets <- as.numeric(controls_count_hets) / nrow(BRIDGE_list)

p_val_homs <- unlist(lapply(lapply(cases_controls_count$gene, function(x)fisher.test(matrix(c(cases_controls_count$cases_count_homs[x],nrow(case_list) - cases_controls_count$cases_count_homs[x],cases_controls_count$controls_count_homs[x],nrow(BRIDGE_list) - cases_controls_count$controls_count_homs[x]),nrow = 2))),"[[", 1))
q_val_homs <- p.adjust(p_val_homs, method = "fdr", n = nrow(cgp_gene_list))
cases_prop_homs <- as.numeric(cases_count_homs) / nrow(case_list)
controls_prop_homs <- as.numeric(controls_count_homs) / nrow(BRIDGE_list)

gene_list_col <- rep("cgp", nrow(cgp_gene_list))
cases_controls_count_cgp <- cbind(cases_controls_count,cases_prop_hets,controls_prop_hets,cases_prop_homs,controls_prop_homs,cases_prop_hets_homs,controls_prop_hets_homs,p_val_hets, q_val_hets,p_val_homs,q_val_homs,p_val_hets_homs,q_val_hets_homs,gene_list_col)

###################
##Mania gene list##
###################


##Reduce table to relevant gene list
SV_table_mania <- SV_table[which(SV_table$identifier %in% mania_gene_list$V1),]

##Subset table into case and control group
SV_table_mania_cases <- SV_table_mania[which(SV_table_mania$BRIDGE_ID %in% case_list$V1),]
SV_table_mania_controls <- SV_table_mania[which(SV_table_mania$BRIDGE_ID %in% BRIDGE_list$V1),]

##Count variants affecting each gene and put into table
cases_count <- lapply(mania_gene_list$V1, function(x) length(unique(SV_table_mania_cases$BRIDGE_ID[which(SV_table_mania_cases$identifier == x)])))
controls_count <- lapply(mania_gene_list$V1, function(x) length(unique(SV_table_mania_controls$BRIDGE_ID[which(SV_table_mania_controls$identifier == x)])))

cases_count_hets <- lapply(mania_gene_list$V1, function(x) length(unique(SV_table_mania_cases$BRIDGE_ID[which(SV_table_mania_cases$identifier == x & (SV_table_mania_cases$MANTA_GT == "0/1" | SV_table_mania_cases$CANVAS_GT == "0/1"))])))
controls_count_hets <- lapply(mania_gene_list$V1, function(x) length(unique(SV_table_mania_controls$BRIDGE_ID[which(SV_table_mania_controls$identifier == x & (SV_table_mania_controls$MANTA_GT == "0/1" | SV_table_mania_controls$CANVAS_GT == "0/1"))])))

cases_count_homs <- lapply(mania_gene_list$V1, function(x) length(unique(SV_table_mania_cases$BRIDGE_ID[which(SV_table_mania_cases$identifier == x & (SV_table_mania_cases$MANTA_GT == "1/1" | SV_table_mania_cases$CANVAS_GT == "1/1"))])))
controls_count_homs <- lapply(mania_gene_list$V1, function(x) length(unique(SV_table_mania_controls$BRIDGE_ID[which(SV_table_mania_controls$identifier == x & (SV_table_mania_controls$MANTA_GT == "1/1" | SV_table_mania_controls$CANVAS_GT == "1/1"))])))

cases_controls_count <- data.frame(mania_gene_list, unlist(cases_count),unlist(controls_count),unlist(cases_count_hets),unlist(controls_count_hets),unlist(cases_count_homs),unlist(controls_count_homs))
colnames(cases_controls_count) <- c("gene", "cases_count", "controls_count","cases_count_hets", "controls_count_hets","cases_count_homs", "controls_count_homs")

##Hypothesis testing and compiling into table
p_val_hets_homs <- unlist(lapply(lapply(cases_controls_count$gene, function(x)fisher.test(matrix(c(cases_controls_count$cases_count[x],nrow(case_list) - cases_controls_count$cases_count[x],cases_controls_count$controls_count[x],nrow(BRIDGE_list) - cases_controls_count$controls_count[x]),nrow = 2))),"[[", 1))
q_val_hets_homs <- p.adjust(p_val_hets_homs, method = "fdr", n = nrow(mania_gene_list))
cases_prop_hets_homs <- as.numeric(cases_count) / nrow(case_list)
controls_prop_hets_homs <- as.numeric(controls_count) / nrow(BRIDGE_list)

p_val_hets <- unlist(lapply(lapply(cases_controls_count$gene, function(x)fisher.test(matrix(c(cases_controls_count$cases_count_hets[x],nrow(case_list) - cases_controls_count$cases_count_hets[x],cases_controls_count$controls_count_hets[x],nrow(BRIDGE_list) - cases_controls_count$controls_count_hets[x]),nrow = 2))),"[[", 1))
q_val_hets <- p.adjust(p_val_hets, method = "fdr", n = nrow(mania_gene_list))
cases_prop_hets <- as.numeric(cases_count_hets) / nrow(case_list)
controls_prop_hets <- as.numeric(controls_count_hets) / nrow(BRIDGE_list)

p_val_homs <- unlist(lapply(lapply(cases_controls_count$gene, function(x)fisher.test(matrix(c(cases_controls_count$cases_count_homs[x],nrow(case_list) - cases_controls_count$cases_count_homs[x],cases_controls_count$controls_count_homs[x],nrow(BRIDGE_list) - cases_controls_count$controls_count_homs[x]),nrow = 2))),"[[", 1))
q_val_homs <- p.adjust(p_val_homs, method = "fdr", n = nrow(mania_gene_list))
cases_prop_homs <- as.numeric(cases_count_homs) / nrow(case_list)
controls_prop_homs <- as.numeric(controls_count_homs) / nrow(BRIDGE_list)

gene_list_col <- rep("mania", nrow(mania_gene_list))
cases_controls_count_mania <- cbind(cases_controls_count,cases_prop_hets,controls_prop_hets,cases_prop_homs,controls_prop_homs,cases_prop_hets_homs,controls_prop_hets_homs,p_val_hets, q_val_hets,p_val_homs,q_val_homs,p_val_hets_homs,q_val_hets_homs,gene_list_col)

####################
##Repair gene list##
####################

##Reduce table to relevant gene list
SV_table_repair <- SV_table[which(SV_table$identifier %in% repair_gene_list$V1),]

##Subset table into case and control group
SV_table_repair_cases <- SV_table_repair[which(SV_table_repair$BRIDGE_ID %in% case_list$V1),]
SV_table_repair_controls <- SV_table_repair[which(SV_table_repair$BRIDGE_ID %in% BRIDGE_list$V1),]

##Count variants affecting each gene and put into table
cases_count <- lapply(repair_gene_list$V1, function(x) length(unique(SV_table_repair_cases$BRIDGE_ID[which(SV_table_repair_cases$identifier == x)])))
controls_count <- lapply(repair_gene_list$V1, function(x) length(unique(SV_table_repair_controls$BRIDGE_ID[which(SV_table_repair_controls$identifier == x)])))

cases_count_hets <- lapply(repair_gene_list$V1, function(x) length(unique(SV_table_repair_cases$BRIDGE_ID[which(SV_table_repair_cases$identifier == x & (SV_table_repair_cases$MANTA_GT == "0/1" | SV_table_repair_cases$CANVAS_GT == "0/1"))])))
controls_count_hets <- lapply(repair_gene_list$V1, function(x) length(unique(SV_table_repair_controls$BRIDGE_ID[which(SV_table_repair_controls$identifier == x & (SV_table_repair_controls$MANTA_GT == "0/1" | SV_table_repair_controls$CANVAS_GT == "0/1"))])))

cases_count_homs <- lapply(repair_gene_list$V1, function(x) length(unique(SV_table_repair_cases$BRIDGE_ID[which(SV_table_repair_cases$identifier == x & (SV_table_repair_cases$MANTA_GT == "1/1" | SV_table_repair_cases$CANVAS_GT == "1/1"))])))
controls_count_homs <- lapply(repair_gene_list$V1, function(x) length(unique(SV_table_repair_controls$BRIDGE_ID[which(SV_table_repair_controls$identifier == x & (SV_table_repair_controls$MANTA_GT == "1/1" | SV_table_repair_controls$CANVAS_GT == "1/1"))])))

cases_controls_count <- data.frame(repair_gene_list, unlist(cases_count),unlist(controls_count),unlist(cases_count_hets),unlist(controls_count_hets),unlist(cases_count_homs),unlist(controls_count_homs))
colnames(cases_controls_count) <- c("gene", "cases_count", "controls_count","cases_count_hets", "controls_count_hets","cases_count_homs", "controls_count_homs")

##Hypothesis testing and compiling into table
p_val_hets_homs <- unlist(lapply(lapply(cases_controls_count$gene, function(x)fisher.test(matrix(c(cases_controls_count$cases_count[x],nrow(case_list) - cases_controls_count$cases_count[x],cases_controls_count$controls_count[x],nrow(BRIDGE_list) - cases_controls_count$controls_count[x]),nrow = 2))),"[[", 1))
q_val_hets_homs <- p.adjust(p_val_hets_homs, method = "fdr", n = nrow(repair_gene_list))
cases_prop_hets_homs <- as.numeric(cases_count) / nrow(case_list)
controls_prop_hets_homs <- as.numeric(controls_count) / nrow(BRIDGE_list)

p_val_hets <- unlist(lapply(lapply(cases_controls_count$gene, function(x)fisher.test(matrix(c(cases_controls_count$cases_count_hets[x],nrow(case_list) - cases_controls_count$cases_count_hets[x],cases_controls_count$controls_count_hets[x],nrow(BRIDGE_list) - cases_controls_count$controls_count_hets[x]),nrow = 2))),"[[", 1))
q_val_hets <- p.adjust(p_val_hets, method = "fdr", n = nrow(repair_gene_list))
cases_prop_hets <- as.numeric(cases_count_hets) / nrow(case_list)
controls_prop_hets <- as.numeric(controls_count_hets) / nrow(BRIDGE_list)

p_val_homs <- unlist(lapply(lapply(cases_controls_count$gene, function(x)fisher.test(matrix(c(cases_controls_count$cases_count_homs[x],nrow(case_list) - cases_controls_count$cases_count_homs[x],cases_controls_count$controls_count_homs[x],nrow(BRIDGE_list) - cases_controls_count$controls_count_homs[x]),nrow = 2))),"[[", 1))
q_val_homs <- p.adjust(p_val_homs, method = "fdr", n = nrow(repair_gene_list))
cases_prop_homs <- as.numeric(cases_count_homs) / nrow(case_list)
controls_prop_homs <- as.numeric(controls_count_homs) / nrow(BRIDGE_list)

gene_list_col <- rep("repair", nrow(repair_gene_list))
cases_controls_count_repair <- cbind(cases_controls_count,
                                     cases_prop_hets,
                                     controls_prop_hets,
                                     cases_prop_homs,
                                     controls_prop_homs,
                                     cases_prop_hets_homs,
                                     controls_prop_hets_homs,
                                     p_val_hets, q_val_hets,p_val_homs,q_val_homs,p_val_hets_homs,q_val_hets_homs,gene_list_col)

#########################################################
##Compile results from all gene lists into single table##
#########################################################
compiled_SV_table <- rbind(cases_controls_count_full,cases_controls_count_webgestalt,cases_controls_count_loftool,cases_controls_count_cgp,cases_controls_count_mania,cases_controls_count_repair)

tumour_query <- paste(args[1],args[2],args[3],args[4],args[5],args[6],args[7],args[8],args[9],sep = "")
tumour_col <- rep(tumour_query, nrow(compiled_SV_table))

compiled_SV_table <- cbind(compiled_SV_table,tumour_col)

compiled_SV_table_sig <- compiled_SV_table[which((compiled_SV_table$q_val_hets < 0.05 & compiled_SV_table$cases_prop_hets > compiled_SV_table$controls_prop_hets) | (compiled_SV_table$q_val_homs < 0.05 & compiled_SV_table$cases_prop_homs > compiled_SV_table$controls_prop_homs) | (compiled_SV_table$q_val_hets_homs < 0.05 & compiled_SV_table$cases_prop_hets_homs > compiled_SV_table$controls_prop_hets_homs)),]

##Output results
write.csv(compiled_SV_table, paste("/home/jww39/new_truncations/with_internal_af_filter/SV_search/", tumour_query,"_SV_table_all.csv", sep = "")) 
write.csv(compiled_SV_table_sig, paste("/home/jww39/new_truncations/with_internal_af_filter/SV_search/", tumour_query,"_SV_table_sig.csv", sep = "")) 

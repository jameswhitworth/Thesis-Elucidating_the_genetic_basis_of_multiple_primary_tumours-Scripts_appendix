###############################################################################################################
##R script to count structural variant calls affecting elements of interest and compare frequency vs controls##
###############################################################################################################

############################################################################################################
##Following commands used to launch script according to phenotypic subgroup:
#Rscript gtex_SV_detection.R ALL XXX XXX XXX XXX
#Rscript gtex_SV_detection.R 1FROM2 Phaeochromoctoma ACC XXX XXX
#Rscript gtex_SV_detection.R 1FROM2 CNS CNSNervesheath XXX XXX
#Rscript gtex_SV_detection.R SINGLE Breast XXX XXX XXX
#Rscript gtex_SV_detection.R SINGLE Haemlymphoid XXX XXX XXX
#Rscript gtex_SV_detection.R SINGLE Colorectal XXX XXX XXX
#Rscript gtex_SV_detection.R SINGLE Oesophagus XXX XXX XXX
#Rscript gtex_SV_detection.R SINGLE Cardiacmyxoma XXX XXX XXX
#Rscript gtex_SV_detection.R SINGLE Lung XXX XXX XXX
#Rscript gtex_SV_detection.R SINGLE Salivarygland XXX XXX XXX
#Rscript gtex_SV_detection.R SINGLE Connectivetissuesofttissuesarcoma XXX XXX XXX
#Rscript gtex_SV_detection.R 1FROM3 PNSNervesheathbenign PNSNervesheath Nervesheathbenign XXX
#Rscript gtex_SV_detection.R SINGLE Ovary XXX XXX XXX
#Rscript gtex_SV_detection.R SINGLE Pancreas XXX XXX XXX
#Rscript gtex_SV_detection.R SINGLE Pituitary XXX XXX XXX
#Rscript gtex_SV_detection.R SINGLE Prostate XXX XXX XXX
#Rscript gtex_SV_detection.R 1FROM3 NMSC Melanoma Skinbenign XXX
#Rscript gtex_SV_detection.R 1FROM2 Smallbowel GINET XXX XXX
#Rscript gtex_SV_detection.R SINGLE Gastric XXX XXX XXX
#Rscript gtex_SV_detection.R SINGLE Testicular XXX XXX XXX
#Rscript gtex_SV_detection.R SINGLE Thyroid XXX XXX XXX
#Rscript gtex_SV_detection.R 1FROM3 Endometrial Uterineleiomyoma Uterinesarcoma XXX
#Rscript gtex_SV_detection.R 1FROM4 Haemlymphoid Haemmyeloid Haempolycythaemia Haemthrombocythaemia
############################################################################################################
args = commandArgs(trailingOnly=TRUE)

##Read in necessary files
elements <- read.csv("gtex_list.csv", header = T) 
colnames(elements)[2] <- "chrom"
colnames(elements)[3] <- "start"
colnames(elements)[7] <- "identifier"
elements <- elements[elements$GOF_TS == "TS",]

MPT_20170104 <- read.delim("MPT_probands_euro_20170104.txt", header = F)
BRIDGE_list <- read.delim("BRIDGE_euro_controls_20170104.txt", header = F)

##Deletions called by Manta
variants <- read.delim("mantacalls_allsamples_del.ann.filt.0.01_subset.txt",header = T)
variants <- variants[which(variants$MANTA_GQ > 29),]

##Adapt input tables for use
elements$start <- as.numeric(elements$start)
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

##Deletion of entire element
variant_element_match_full_del <- sqldf("select * from variants f1 
                                        left join elements f2 on 
                                        (
                                        (f1.MANTA_CHROM==f2.chrom) and 
                                        
                                        (f1.MANTA_MAX_START < f2.start and f1.MANTA_MIN_END > f2.start) 
                                        )")

mut_type_col <- vector(mode = "character", nrow(variant_element_match_full_del))
mut_type_col[which(!is.na(variant_element_match_full_del$identifier))] <- "full_del"
variant_element_match_full_del <- cbind(variant_element_match_full_del, mut_type_col)

##Collate results
variant_element_match_manta_del <- 
  rbind(
    variant_element_match_full_del[which(!is.na(variant_element_match_full_del$identifier)),]
  )

##Copy number loss called by Canvas
variants <-  read.delim("canvascalls_allsamples_loss.ann.filt.0.01_subset.txt",header = T)
variants <- variants[which(variants$CANVAS_QUAL > 29),]

variants$CANVAS_CHROM <- as.character(variants$CANVAS_CHROM)
variants$CANVAS_START <- as.numeric(variants$CANVAS_START)
variants$CANVAS_END <- as.numeric(variants$CANVAS_END)

##Match variants to elements
library(sqldf)

##Deletion of entire element
variant_element_match_full_del <- sqldf("select * from variants f1 
                                        left join elements f2 on 
                                        (
                                        (f1.CANVAS_CHROM==f2.chrom) and 
                                        
                                        (f1.CANVAS_START < f2.start and f1.CANVAS_END > f2.start) 
                                        )")

mut_type_col <- vector(mode = "character", nrow(variant_element_match_full_del))
mut_type_col[which(!is.na(variant_element_match_full_del$identifier))] <- "full_del"
variant_element_match_full_del <- cbind(variant_element_match_full_del, mut_type_col)

##Collate results
variant_element_match_canvas_del <- 
  rbind(
    variant_element_match_full_del[which(!is.na(variant_element_match_full_del$identifier)),]
  )


##Merge results tables 
library(data.table)
SV_table <- rbindlist(list(variant_element_match_manta_del, variant_element_match_canvas_del), fill=TRUE)
SV_table <- SV_table[which(!is.na(SV_table$mut_type_col)),]

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
if(args[1] != "ALL")colnames(case_list) <- "indv"
if(args[1] == "ALL")case_list <- read.delim("MPT_probands_euro_20170104.txt", header = F)

##Reduce the table down depending on what the tissue/samples of interest are
tissue_equiv <- read.csv("gtex_tissues_with_single_word_equivalents.csv", header = T)

##Get lines with single word of interest in reference table and convert to gtex tissue labels
if(args[1] != "ALL")gtex_tis_1 <- which(tissue_equiv$single_word_equivalent_1 %in% t1)
if(args[1] != "ALL")gtex_tis_2 <- which(tissue_equiv$single_word_equivalent_1 %in% t2)
if(args[1] != "ALL")gtex_tis_3 <- which(tissue_equiv$single_word_equivalent_1 %in% t3)
if(args[1] != "ALL")gtex_tis_4 <- which(tissue_equiv$single_word_equivalent_1 %in% t4)
if(args[1] != "ALL")gtex_tis_all <- c(gtex_tis_1,gtex_tis_2,gtex_tis_3,gtex_tis_4)
if(args[1] != "ALL")gtex_tis_all <- tissue_equiv$gtex.tissue.filename[gtex_tis_all]

##Find which lines in the gtex match table contain the gtex tissue name/s of interest and subset
if(args[1] != "ALL")SV_table <- SV_table[which(SV_table$tissue %in% gtex_tis_all),]

##Read in gene list
full_gene_list <- data.frame(unique(elements$identifier))
colnames(full_gene_list) <- "V1"

#Reduce table to relevant gene list
SV_table_full <- SV_table[which(SV_table$identifier %in% as.character(full_gene_list$V1)),]


##Subset table into case and control group
colnames(case_list) <- "V1"
SV_table_full_cases <- SV_table_full[which(SV_table_full$BRIDGE_ID %in% case_list$V1),]
SV_table_full_controls <- SV_table_full[which(SV_table_full$BRIDGE_ID %in% BRIDGE_list$V1),]

##Count individuals with variants affecting each gene and put into table
cases_count <- lapply(full_gene_list$V1, function(x) length(unique(SV_table_full_cases$BRIDGE_ID[which(SV_table_full_cases$identifier == x)])))
controls_count <- lapply(full_gene_list$V1, function(x) length(unique(SV_table_full_controls$BRIDGE_ID[which(SV_table_full_controls$identifier == x)])))

cases_count_hets <- lapply(full_gene_list$V1, function(x) length(unique(SV_table_full_cases$BRIDGE_ID[which(SV_table_full_cases$identifier == x & (SV_table_full_cases$MANTA_GT == "0/1" | SV_table_full_cases$CANVAS_GT == "0/1"))])))
controls_count_hets <- lapply(full_gene_list$V1, function(x) length(unique(SV_table_full_controls$BRIDGE_ID[which(SV_table_full_controls$identifier == x & (SV_table_full_controls$MANTA_GT == "0/1" | SV_table_full_controls$CANVAS_GT == "0/1"))])))

cases_count_homs <- lapply(full_gene_list$V1, function(x) length(unique(SV_table_full_cases$BRIDGE_ID[which(SV_table_full_cases$identifier == x & (SV_table_full_cases$MANTA_GT == "1/1" | SV_table_full_cases$CANVAS_GT == "1/1"))])))
controls_count_homs <- lapply(full_gene_list$V1, function(x) length(unique(SV_table_full_controls$BRIDGE_ID[which(SV_table_full_controls$identifier == x & (SV_table_full_controls$MANTA_GT == "1/1" | SV_table_full_controls$CANVAS_GT == "1/1"))])))

cases_controls_count <- data.frame(full_gene_list, unlist(cases_count),unlist(controls_count),unlist(cases_count_hets),unlist(controls_count_hets),unlist(cases_count_homs),unlist(controls_count_homs))
colnames(cases_controls_count) <- c("gene", "cases_count", "controls_count","cases_count_hets", "controls_count_hets","cases_count_homs", "controls_count_homs")


##Hypothesis testing and compilation into table
p_val_hets_homs <- unlist(lapply(lapply(1:length(cases_controls_count$gene), function(x) fisher.test(matrix(c(cases_controls_count$cases_count[x],nrow(case_list) - cases_controls_count$cases_count[x], cases_controls_count$controls_count[x], nrow(BRIDGE_list) - cases_controls_count$controls_count[x]),nrow = 2))),"[[", 1))
q_val_hets_homs <- p.adjust(p_val_hets_homs, method = "fdr", n = nrow(full_gene_list))
cases_prop_hets_homs <- as.numeric(cases_count) / nrow(case_list)
controls_prop_hets_homs <- as.numeric(controls_count) / nrow(BRIDGE_list)

p_val_hets <- unlist(lapply(lapply(1:length(cases_controls_count$gene), function(x) fisher.test(matrix(c(cases_controls_count$cases_count_hets[x],nrow(case_list) - cases_controls_count$cases_count_hets[x],cases_controls_count$controls_count_hets[x],nrow(BRIDGE_list) - cases_controls_count$controls_count_hets[x]),nrow = 2))),"[[", 1))
q_val_hets <- p.adjust(p_val_hets, method = "fdr", n = nrow(full_gene_list))
cases_prop_hets <- as.numeric(cases_count_hets) / nrow(case_list)
controls_prop_hets <- as.numeric(controls_count_hets) / nrow(BRIDGE_list)

p_val_homs <- unlist(lapply(lapply(1:length(cases_controls_count$gene), function(x)fisher.test(matrix(c(cases_controls_count$cases_count_homs[x],nrow(case_list) - cases_controls_count$cases_count_homs[x],cases_controls_count$controls_count_homs[x],nrow(BRIDGE_list) - cases_controls_count$controls_count_homs[x]),nrow = 2))),"[[", 1))
q_val_homs <- p.adjust(p_val_homs, method = "fdr", n = nrow(full_gene_list))
cases_prop_homs <- as.numeric(cases_count_homs) / nrow(case_list)
controls_prop_homs <- as.numeric(controls_count_homs) / nrow(BRIDGE_list)

cases_controls_count_full <- cbind(cases_controls_count,cases_prop_hets,controls_prop_hets,cases_prop_homs,controls_prop_homs,cases_prop_hets_homs,controls_prop_hets_homs,p_val_hets, q_val_hets,p_val_homs,q_val_homs,p_val_hets_homs,q_val_hets_homs)

##Compile results from all gene lists into single table. Add reference to phenotypic subgroup
compiled_SV_table <- rbind(cases_controls_count_full)

tumour_query <- paste(args[1],args[2],args[3],args[4],args[5],sep = "")
tumour_col <- rep(tumour_query, nrow(compiled_SV_table))

compiled_SV_table <- cbind(compiled_SV_table,tumour_col)

compiled_SV_table_sig <- compiled_SV_table[which((compiled_SV_table$q_val_hets < 0.05 & compiled_SV_table$cases_prop_hets > compiled_SV_table$controls_prop_hets) | (compiled_SV_table$q_val_homs < 0.05 & compiled_SV_table$cases_prop_homs > compiled_SV_table$controls_prop_homs) | (compiled_SV_table$q_val_hets_homs < 0.05 & compiled_SV_table$cases_prop_hets_homs > compiled_SV_table$controls_prop_hets_homs)),]

##Output results
write.csv(compiled_SV_table, paste("/home/jww39/non_coding/gtex/with_internal_af_filter/gtex_SV_search/", tumour_query,"_SV_table_all.csv", sep = "")) 
write.csv(compiled_SV_table_sig, paste("/home/jww39/non_coding/gtex/with_internal_af_filter/gtex_SV_search/", tumour_query,"_SV_table_sig.csv", sep = "")) 

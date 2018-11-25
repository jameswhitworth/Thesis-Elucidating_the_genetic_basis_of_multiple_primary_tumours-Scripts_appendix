############################################################################################################
##R script to count structural variant calls affecting genes of interest and compare frequency vs controls##
############################################################################################################

############################################################################################################
##Following commands used to launch script according to phenotypic subgroup:
#Rscript TELgenes_SV_detection.R LOW
#Rscript TELgenes_SV_detection.R HIGH
############################################################################################################
args = commandArgs(trailingOnly=TRUE)

##Read in necessary files
elements <- read.csv("TELgenes_gene_list_regions_exons.csv", header = T) 
elements <- elements[which(!(is.na(elements$start))),]

MPT_20170104 <- read.delim("MPT_probands_euro_20170104.txt", header = F)
BRIDGE_list <- read.delim("BRIDGE_euro_controls_20170104.txt", header = F)


##Deletions called by Manta
variants <- read.delim("mantacalls_allsamples_del.ann.filt.0.01_subset.txt",header = T)
variants <- variants[which(variants$MANTA_GQ > 29),]

##Adapt input tables for use
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

##Deletion with breakpoints within gene
variant_element_match_intra_del <- sqldf("select * from variants f1 
                               left join elements f2 on 
                               (
(f1.MANTA_CHROM==f2.chrom) and 

(f1.MANTA_MIN_START > f2.start and f1.MANTA_MAX_END < f2.end) 
)")

mut_type_col <- vector(mode = "character", nrow(variant_element_match_intra_del))
mut_type_col[which(!is.na(variant_element_match_intra_del$identifier))] <- "intra_del"
variant_element_match_intra_del <- cbind(variant_element_match_intra_del, mut_type_col)

##Deletion with breakpoints incorporating one, but not the other, end of gene
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
variants <-  read.delim("mantacalls_allsamples_bnd.ann.filt.0.01_subset.txt",header = T)
variants <- variants[which(variants$MANTA_GQ > 29),]

elements <- read.csv("TELgenes_gene_list_regions_genes.csv", header = T) 

##Adapt input tables for use
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

##Translocation with breakpint within gene
variant_element_match_bnd <- sqldf("select * from variants f1 
                               left join elements f2 on 
                               (
(f1.MANTA_CHROM==f2.chrom) and 

(
(f1.MANTA_MIN_START > f2.start and f1.MANTA_MAX_START < f2.end) or (f1.MANTA_MIN_END > f2.start and f1.MANTA_MAX_END < f2.end)
)

)")

mut_type_col <- vector(mode = "character", nrow(variant_element_match_bnd))
mut_type_col[which(!is.na(variant_element_match_bnd$identifier))] <- "bnd"
variant_element_match_bnd <- cbind(variant_element_match_bnd, mut_type_col)
variant_element_match_bnd <- variant_element_match_bnd[which(!is.na(variant_element_match_bnd$identifier)),]

##Inversions called by Manta
variants <-  read.delim("mantacalls_allsamples_inv.ann.filt.0.01_subset.txt", header = T)
variants <- variants[which(variants$MANTA_GQ > 29),]

elements <- read.csv("TELgenes_gene_list_regions_genes.csv", header = T) 

##Adapt input tables for use
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

##Inversion with breakpoint within gene
variant_element_match_inv <- sqldf("select * from variants f1 
                               left join elements f2 on 
                               (
(f1.MANTA_CHROM==f2.chrom) and 

(
(f1.MANTA_MIN_START > f2.start and f1.MANTA_MAX_START < f2.end) or (f1.MANTA_MIN_END > f2.start and f1.MANTA_MAX_END < f2.end)
)

)")

mut_type_col <- vector(mode = "character", nrow(variant_element_match_inv))
mut_type_col[which(!is.na(variant_element_match_inv$identifier))] <- "inv"
variant_element_match_inv <- cbind(variant_element_match_inv, mut_type_col)
variant_element_match_inv <- variant_element_match_inv[which(!is.na(variant_element_match_inv$identifier)),]


##Insertions called by Manta
variants <-  read.delim("mantacalls_allsamples_ins.ann.filt.0.01_subset.txt",header = T)
variants <- variants[which(variants$MANTA_GQ > 29),]

elements <- read.csv("TELgenes_gene_list_regions_genes.csv", header = T) 

##Adapt input tables for use
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

##Insertions with breakpoint within gene
variant_element_match_ins <- sqldf("select * from variants f1 
                               left join elements f2 on 
                               (
(f1.MANTA_CHROM==f2.chrom) and 

(
(f1.MANTA_MIN_START > f2.start and f1.MANTA_MAX_START < f2.end) or (f1.MANTA_MIN_END > f2.start and f1.MANTA_MAX_END < f2.end)
)
)")

mut_type_col <- vector(mode = "character", nrow(variant_element_match_ins))
mut_type_col[which(!is.na(variant_element_match_ins$identifier))] <- "ins"
variant_element_match_ins <- cbind(variant_element_match_ins, mut_type_col)
variant_element_match_ins <- variant_element_match_ins[which(!is.na(variant_element_match_ins$identifier)),]


##Copy number loss called by Canvas
variants <-  read.delim("canvascalls_allsamples_loss.ann.filt.0.01_subset.txt",header = T)
variants <- variants[which(variants$CANVAS_QUAL > 29),]

elements <- read.csv("TELgenes_gene_list_regions_exons.csv", header = T) 

elements$start <- as.numeric(elements$start)
elements$end <- as.numeric(elements$end)
elements$chrom <- as.character(elements$chrom)

variants$CANVAS_CHROM <- as.character(variants$CANVAS_CHROM)
variants$CANVAS_START <- as.numeric(variants$CANVAS_START)
variants$CANVAS_END <- as.numeric(variants$CANVAS_END)

##Match variants to elements
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



##Deletion with breakpoints within gene
variant_element_match_intra_del <- sqldf("select * from variants f1 
                                         left join elements f2 on 
                                         (
                                         (f1.CANVAS_CHROM==f2.chrom) and 
                                         
                                         (f1.CANVAS_START > f2.start and f1.CANVAS_END < f2.end) 
                                         )")

mut_type_col <- vector(mode = "character", nrow(variant_element_match_intra_del))
mut_type_col[which(!is.na(variant_element_match_intra_del$identifier))] <- "intra_del"
variant_element_match_intra_del <- cbind(variant_element_match_intra_del, mut_type_col)

##Deletion with breakpoints incorporating one, but not the other, end of gene
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

##Merge results tables
library(data.table)

SV_table <- rbindlist(list(variant_element_match_manta_del, variant_element_match_bnd, variant_element_match_inv, variant_element_match_ins, variant_element_match_canvas_del), fill=TRUE)

#Remove lines with no variant element match
SV_table <- SV_table[which(!is.na(SV_table$mut_type_col)),]

##Create list of cases based on previous telomere length estimates
if(args[1] == "LOW")case_list <- read.delim("low_resid_MPMT_samples.txt", header = F)
if(args[1] == "HIGH")case_list <- read.delim("high_resid_MPMT_samples.txt", header = F)
case_list <- data.frame(case_list[which(case_list$V1 %in% MPT_20170104$V1),])
colnames(case_list) <- "V1"

##Read in gene list
gene_list <- read.delim("TELgenes_gene_list.txt", header = F)


##Reduce table to relevant gene list
SV_table <- SV_table[which(SV_table$identifier %in% gene_list$V1),]


##Subset table into case and control group
SV_table_cases <- SV_table[which(SV_table$BRIDGE_ID %in% case_list$V1),]
SV_table_controls <- SV_table[which(SV_table$BRIDGE_ID %in% BRIDGE_list$V1),]


##Count variants affecting each gene and put into table
cases_count <- lapply(gene_list$V1, function(x) length(unique(SV_table_cases$BRIDGE_ID[which(SV_table_cases$identifier == x)])))
controls_count <- lapply(gene_list$V1, function(x) length(unique(SV_table_controls$BRIDGE_ID[which(SV_table_controls$identifier == x)])))

cases_count_hets <- lapply(gene_list$V1, function(x) length(unique(SV_table_cases$BRIDGE_ID[which(SV_table_cases$identifier == x & (SV_table_cases$MANTA_GT == "0/1" | SV_table_cases$CANVAS_GT == "0/1"))])))
controls_count_hets <- lapply(gene_list$V1, function(x) length(unique(SV_table_controls$BRIDGE_ID[which(SV_table_controls$identifier == x & (SV_table_controls$MANTA_GT == "0/1" | SV_table_controls$CANVAS_GT == "0/1"))])))

cases_count_homs <- lapply(gene_list$V1, function(x) length(unique(SV_table_cases$BRIDGE_ID[which(SV_table_cases$identifier == x & (SV_table_cases$MANTA_GT == "1/1" | SV_table_cases$CANVAS_GT == "1/1"))])))
controls_count_homs <- lapply(gene_list$V1, function(x) length(unique(SV_table_controls$BRIDGE_ID[which(SV_table_controls$identifier == x & (SV_table_controls$MANTA_GT == "1/1" | SV_table_controls$CANVAS_GT == "1/1"))])))

cases_controls_count <- data.frame(gene_list, unlist(cases_count),unlist(controls_count),unlist(cases_count_hets),unlist(controls_count_hets),unlist(cases_count_homs),unlist(controls_count_homs))
colnames(cases_controls_count) <- c("gene", "cases_count", "controls_count","cases_count_hets", "controls_count_hets","cases_count_homs", "controls_count_homs")

##Hypothesis testing and compile into table
p_val_hets_homs <- unlist(lapply(lapply(cases_controls_count$gene, function(x)fisher.test(matrix(c(cases_controls_count$cases_count[x],nrow(case_list) - cases_controls_count$cases_count[x],cases_controls_count$controls_count[x],nrow(BRIDGE_list) - cases_controls_count$controls_count[x]),nrow = 2))),"[[", 1))
q_val_hets_homs <- p.adjust(p_val_hets_homs, method = "fdr", n = nrow(gene_list))
cases_prop_hets_homs <- as.numeric(cases_count) / nrow(case_list)
controls_prop_hets_homs <- as.numeric(controls_count) / nrow(BRIDGE_list)

p_val_hets <- unlist(lapply(lapply(cases_controls_count$gene, function(x)fisher.test(matrix(c(cases_controls_count$cases_count_hets[x],nrow(case_list) - cases_controls_count$cases_count_hets[x],cases_controls_count$controls_count_hets[x],nrow(BRIDGE_list) - cases_controls_count$controls_count_hets[x]),nrow = 2))),"[[", 1))
q_val_hets <- p.adjust(p_val_hets, method = "fdr", n = nrow(gene_list))
cases_prop_hets <- as.numeric(cases_count_hets) / nrow(case_list)
controls_prop_hets <- as.numeric(controls_count_hets) / nrow(BRIDGE_list)

p_val_homs <- unlist(lapply(lapply(cases_controls_count$gene, function(x)fisher.test(matrix(c(cases_controls_count$cases_count_homs[x],nrow(case_list) - cases_controls_count$cases_count_homs[x],cases_controls_count$controls_count_homs[x],nrow(BRIDGE_list) - cases_controls_count$controls_count_homs[x]),nrow = 2))),"[[", 1))
q_val_homs <- p.adjust(p_val_homs, method = "fdr", n = nrow(gene_list))
cases_prop_homs <- as.numeric(cases_count_homs) / nrow(case_list)
controls_prop_homs <- as.numeric(controls_count_homs) / nrow(BRIDGE_list)

cases_controls_count <- cbind(cases_controls_count,cases_prop_hets,controls_prop_hets,cases_prop_homs,controls_prop_homs,cases_prop_hets_homs,controls_prop_hets_homs,p_val_hets, q_val_hets,p_val_homs,q_val_homs,p_val_hets_homs,q_val_hets_homs)

compiled_SV_table <- cases_controls_count

compiled_SV_table_sig <- compiled_SV_table[which((compiled_SV_table$q_val_hets < 0.05 & compiled_SV_table$cases_prop_hets > compiled_SV_table$controls_prop_hets) | (compiled_SV_table$q_val_homs < 0.05 & compiled_SV_table$cases_prop_homs > compiled_SV_table$controls_prop_homs) | (compiled_SV_table$q_val_hets_homs < 0.05 & compiled_SV_table$cases_prop_hets_homs > compiled_SV_table$controls_prop_hets_homs)),]

##Output results
if(args[1] == "LOW")write.csv(compiled_SV_table, paste("/home/jww39/candidate_gene_searches/TELgenes_ALL_XXX_XXX_XXX_XXX_XXX_XXX_XXX_XXX_XXX_XXX/SV_search/","low_resid_SV_table_all.csv", sep = "")) 
if(args[1] == "LOW")write.csv(compiled_SV_table_sig, paste("/home/jww39/candidate_gene_searches/TELgenes_ALL_XXX_XXX_XXX_XXX_XXX_XXX_XXX_XXX_XXX_XXX/SV_search/", "low_resid_SV_table_sig.csv", sep = "")) 


if(args[1] == "HIGH")write.csv(compiled_SV_table, paste("/home/jww39/candidate_gene_searches/TELgenes_ALL_XXX_XXX_XXX_XXX_XXX_XXX_XXX_XXX_XXX_XXX/SV_search/","high_resid_SV_table_all.csv", sep = "")) 
if(args[1] == "HIGH")write.csv(compiled_SV_table_sig, paste("/home/jww39/candidate_gene_searches/TELgenes_ALL_XXX_XXX_XXX_XXX_XXX_XXX_XXX_XXX_XXX_XXX/SV_search/", "high_resid_SV_table_sig.csv", sep = "")) 

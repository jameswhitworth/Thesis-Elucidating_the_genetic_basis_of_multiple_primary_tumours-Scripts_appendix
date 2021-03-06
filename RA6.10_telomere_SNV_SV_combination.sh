################################################################################################################################################################
##R script to count combination of single nucleotide variant, indel and structural variant calls affecting genes of interest and compare frequency vs controls##
################################################################################################################################################################

##Utilises tables of variants and counts generated by scripts for SNV/indel analysis and SV analysis

############################################################################################################
##Following commands used to launch script according to phenotypic subgroup:
#Rscript TELgenes_SNV_SV_combination.R LOW
#Rscript TELgenes_SNV_SV_combination.R HIGH
############################################################################################################
args = commandArgs(trailingOnly=TRUE)

##Read in tables of variants and combine
if(args[1] == "LOW")SNV <- read.csv("TELgenes_low_variants_per_gene_table_all.csv")
if(args[1] == "HIGH")SNV <- read.csv("TELgenes_high_variants_per_gene_table_all.csv")
SNV <- data.frame(SNV$gene, SNV$indv_per_gene_count_cases_hets, SNV$indv_per_gene_count_BRIDGES_hets, SNV$indv_per_gene_count_cases_homs, SNV$indv_per_gene_count_BRIDGES_homs)
colnames(SNV) <- c("gene", "cases_count_hets", "controls_count_hets", "cases_count_homs", "controls_count_homs")

if(args[1] == "LOW")SV <- read.csv("low_resid_SV_table_all.csv")
if(args[1] == "HIGH")SV <- read.csv("high_resid_SV_table_all.csv")
SV <- data.frame(SV$gene, SV$cases_count_hets, SV$controls_count_hets, SV$cases_count_homs, SV$controls_count_homs)
colnames(SV) <- c("gene", "cases_count_hets", "controls_count_hets", "cases_count_homs", "controls_count_homs")

library(data.table)
SNV_SV_table <- rbindlist(list(SV, SNV), fill=TRUE)

##Read in gene list
gene_list <- read.delim("TELgenes_gene_list.txt", header = F)

##Produce sums of SV and SNV counts for each gene
cases_count_hets <- unlist(lapply(gene_list$V1, function(x) sum(SNV_SV_table$cases_count_hets[which(SNV_SV_table$gene == x)])))
controls_count_hets <- unlist(lapply(gene_list$V1, function(x) sum(SNV_SV_table$controls_count_hets[which(SNV_SV_table$gene == x)])))

cases_count_homs <- unlist(lapply(gene_list$V1, function(x) sum(SNV_SV_table$cases_count_homs[which(SNV_SV_table$gene == x)])))
controls_count_homs <- unlist(lapply(gene_list$V1, function(x) sum(SNV_SV_table$controls_count_homs[which(SNV_SV_table$gene == x)])))

cases_count <- cases_count_hets + cases_count_homs
controls_count <- controls_count_hets + controls_count_homs

cases_controls_count <- data.frame(gene_list, cases_count, controls_count, cases_count_hets, controls_count_hets, cases_count_homs, controls_count_homs)
colnames(cases_controls_count)[1] <- "gene"

##Obtain relevant samples for lists of cases and controls
MPT_20170104 <- read.delim("MPT_probands_euro_20170104.txt", header = F)
BRIDGE_list <- read.delim("BRIDGE_euro_controls_20170104.txt", header = F)

if(args[1] == "LOW")case_list <- read.delim("low_resid_MPMT_samples.txt", header = F)
if(args[1] == "HIGH")case_list <- read.delim("high_resid_MPMT_samples.txt", header = F)
case_list <- data.frame(case_list[which(case_list$V1 %in% MPT_20170104$V1),])  
colnames(case_list) <- "indv"

##Hypothesis testing and compilation into table
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

cases_controls_count_sig <- cases_controls_count[which((cases_controls_count$q_val_hets < 0.05 & cases_controls_count$cases_prop_hets > cases_controls_count$controls_prop_hets) | (cases_controls_count$q_val_homs < 0.05 & cases_controls_count$cases_prop_homs > cases_controls_count$controls_prop_homs) | (cases_controls_count$q_val_hets_homs < 0.05 & cases_controls_count$cases_prop_hets_homs > cases_controls_count$controls_prop_hets_homs)),]

##Output results
write.csv(cases_controls_count, paste("/home/jww39/candidate_gene_searches/TELgenes_ALL_XXX_XXX_XXX_XXX_XXX_XXX_XXX_XXX_XXX_XXX/SNV_SV/", 
                                      paste(args[1],sep = ""),
                                      "SNV_SV_table_all.csv", sep = "")) 

write.csv(cases_controls_count_sig, paste("/home/jww39/candidate_gene_searches/TELgenes_ALL_XXX_XXX_XXX_XXX_XXX_XXX_XXX_XXX_XXX_XXX/SNV_SV/", 
                                          paste(args[1],sep = "")
                                          ,"SNV_SV_table_sig.csv", sep = "")) 


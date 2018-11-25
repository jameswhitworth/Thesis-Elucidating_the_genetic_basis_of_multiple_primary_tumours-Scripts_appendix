#############################################################################################################
##Filtering Illumina TruSight Cancer panel data for possible mosaic variants in cancer predisposition genes##
#############################################################################################################

##Read in necessary files
panel_variants <- read.table("/home/jww39/mosaicism/panel_scores_only/combined_spaces_removed.txt", header = T) ##Table of annotated variants called fro mpanel data
refseq_transcripts <- readLines("vp_tscp_refseq_no_recessive_transcripts.txt") ##List of refseq equivalents of canonical transcripts for genes of interest
trunc_consequence <- readLines("annovar_trunc_consequence.txt") ##List of truncating variant consequences as described in annotation
non_trunc_consequence <- readLines("annovar_non_trunc_consequence.txt") #List of non-truncating protein altering variant consequences as described in annotation
probands <- readLines("panel_probands.txt") ##List of included samples

##Filter variants according to criteria indicating possible mosaic variant of interest
mosaic_filtered_variants <- subset(panel_variants, panel_variants$XAcorrect == 0 &
                   panel_variants$VAF > 0.05 &
                   panel_variants$VAF < 0.3 &
                   panel_variants$Depth > 199 &
                   panel_variants$TranscriptID %in% refseq_transcripts &
                   panel_variants$X1000g2015aug_all <= 0.01 & ##1000 genomes allele frequency
                     panel_variants$sample %in% probands)

consequence_filtered_variants <- subset(mosaic_filtered_variants, 
                 (mosaic_filtered_variants$ExonicFunc.refGene %in% trunc_consequence | 
                 (mosaic_filtered_variants$ExonicFunc.refGene %in% non_trunc_consequence & (mosaic_filtered_variants$CLINSIG != "-9" | mosaic_filtered_variants$Score >= 0.75))) & ##In silico score and ClinVar annotation incorporated
                   mosaic_filtered_variants$ExonicFunc.refGene != "synonymous_SNV")

##Output filtered variants
write.csv(consequence_filtered_variants, file = "filtered_variants.csv")

###################################
##Calculating covergae statistics##
###################################

##Calculate sequencing depth per target base
for i in `cat bams.txt`; do
	samtools depth -b mosaic_exons.bed /data/Alignments/CancerPanel_MPT_hg38.bwa/${i} > ${i}.cov
done

##Collate results
ls *.cov > cov_files.txt
for i in `cat cov_files.txt`; do
  tr '\t' ',' < ${i} > ${i}.csv
done

touch mosaic_exons_cov.csv

for i in `cat cov_files.txt`; do
  cat ${i}.csv >> mosaic_exons_cov.csv
done

##R Script to calculate statistics (run in R studio)

#Read in samtools depth output combined across samples
mosaic_exons_cov <- read.csv("mosaic_exons_cov.csv", header = F)
colnames(mosaic_exons_cov) <- c("chr", "base", "depth")

#Calculate mean coverage
mean_mosaic_exons_cov <- mean(mosaic_exons_cov$depth, na.rm = T)
sd_mosaic_exons_cov <- sd(mosaic_exons_cov$depth, na.rm = T)

##Calculate percentage bases covered at 200X or above
mosaic_exons_cov_percent_at_200 <- (100/nrow(mosaic_exons_cov)) * (length(mosaic_exons_cov$depth[which(mosaic_exons_cov$depth > 199)]))

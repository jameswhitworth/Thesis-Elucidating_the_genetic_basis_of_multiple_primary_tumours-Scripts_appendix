
#########################################
##Single nucleotide variants and indels##
#########################################

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
bcftools view -R virtual_panel_exons_sorted.bed -e 'FILTER!="PASS"' -o ${i}_filtered.vcf -O v /scratch/WGS10K/data/release/20170614-A/merged-vcf/no_hgmd/${i}.vcf.gz
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

bcftools concat -Oz `cat files_to_merge.txt` -o vp_merged.vcf.gz

tabix vp_merged.vcf.gz
vcftools --gzvcf vp_merged.vcf.gz --out vp_merged_samples --recode --recode-INFO-all --keep samples_without_suffix.txt ##File containing filenames of all individuals included in analysis

##Filter merged VCF based on quality parameters (insert missing genotypes where criteria not fulfilled)
bgzip vp_merged_samples.recode.vcf
tabix vp_merged_samples.recode.vcf.gz

bcftools filter -e 'FMT/DP<10 || FMT/GQ<30 || (FMT/AD[1]/(FMT/AD[0] + FMT/AD[1]) <0.33) || (FMT/AD[1]/(FMT/AD[0] + FMT/AD[1] + FMT/AD[2]) <0.33) || (FMT/AD[2]/(FMT/AD[0] + FMT/AD[1] + FMT/AD[2]) <0.33)' -S . -o vp_merged_samples_filtered.vcf.gz -O z vp_merged_samples.recode.vcf.gz

#####################################################
##Transfer merged filtered VCF file to local server##
#####################################################

##Remove INFO fields from BioResource annotated files
bcftools annotate -x INFO vp_merged_samples_filtered.vcf -o vp_SNV_indel.vcf

##Annotate with Variant Effect Predictor
/home/jww39/ensembl-vep/vep -i vp_SNV_indel.vcf -offline --assembly GRCh37 -o VEP_out_vp_SNV_indel.vcf --pick --pick_order canonical --everything --vcf --plugin LoF --plugin CADD,/home/jww39/.vep/Plugins/whole_genome_SNVs.tsv.gz,/home/jww39/.vep/Plugins/InDels.tsv.gz --plugin ExAC,/home/jww39/.vep/Plugins/ExAC.r0.3.sites.vep.vcf.gz

##Filter with Variant Effect Predictor filter script based on annotations
/home/jww39/ensembl-vep/filter_vep -i VEP_out_vp_SNV_indel.vcf -o VEP_out_filtered_vp_SNV_indel.vcf -filter "Consequence in consequence.txt" --only_matched
/home/jww39/ensembl-vep/filter_vep -i VEP_out_filtered_vp_SNV_indel.vcf -o VEP_out_double_filtered_vp_SNV_indel.vcf -filter "(ExAC_AF < 0.01 or not ExAC_AF) and (AF < 0.01 or not AF)" --only_matched

#Prepare resulting VCF for reading into R
sed -i 's/#CHROM/CHROM/g' VEP_out_double_filtered_vp_SNV_indel.vcf
sed '/^#/ d' VEP_out_double_filtered_vp_SNV_indel.vcf > VEP_out_double_filtered_vp_SNV_indel_without_header.vcf

#######################################################################
##R script to provide further annotation to assist variant assessment##
#######################################################################

##Read in variants
variants <- read.table("VEP_out_double_filtered_vp_SNV_indel_without_header.vcf", header =T)
variants_gt_only <- variants[10:ncol(variants)]

##Read in other reference files
GOF_genes <- readLines("GOF_genes.txt") #File containing gene names of proto-oncogene cancer predisposition genes
recessive_genes <- readLines("recessive_genes.txt") #File containing gene names recessive cancer predisposition genes

##Identify samples containing each variant (heterozygotes)
samples_raw_hets <- data.frame()
for (i in 1:nrow(variants_gt_only)) {
  a <- unlist(variants_gt_only[i,])
  b <- grep("^0/1", a, perl = TRUE)
  c <- as.vector(names(a[b]))
  d <- paste(c, collapse = ';')
  samples_raw_hets <- append(samples_raw_hets, d)
} 
samples_hets <- unlist(samples_raw_hets)

##Identify samples containing each variant (homozygotes)
samples_raw_homs <- data.frame()
for (i in 1:nrow(variants_gt_only)) {
  
  a <- unlist(variants_gt_only[i,])
  b <- grep("^1/1", a, perl = TRUE)
  c <- as.vector(names(a[b]))
  d <- paste(c, collapse = ';')
  samples_raw_homs <- append(samples_raw_homs, d)

} 
samples_homs <- unlist(samples_raw_homs)

##Compile table of variants with variant and sample information 
annotated_variants <- data.frame(variants[,1:9], samples_hets, samples_homs, variants_gt_only)

##Subset table of variants to only include variants present in at least one sample
annotated_variants <- annotated_variants[which(annotated_variants$samples_hets != "" | annotated_variants$samples_homs != "" ),]

##Produce a new column to tag variants occuring in proto-oncogene cancer predisposition genes
GOF_indices <- data.frame()
for (i in 1:length(GOF_genes)) {
  a <- grep(GOF_genes[i], annotated_variants$INFO)
  GOF_indices <- append(GOF_indices, a)
}
GOF_col <- rep(0,nrow(annotated_variants))
GOF_col[unlist(GOF_indices)] <- "YES"

##Produce a new column to tag variants occuring in recessive cancer predisposition genes
recessive_indices <- data.frame(
for (i in 1:length(recessive_genes)) {
  a <- grep(recessive_genes[i], annotated_variants$INFO)
  recessive_indices <- append(recessive_indices, a)
}
recessive_col <- rep(0,nrow(annotated_variants))
recessive_col[unlist(recessive_indices)] <- "YES"

##Produce a new column to tag variants with predicted truncating consequence
trunc_indices <- grep("\\|HIGH\\|", annotated_variants$INFO)
trunc_col <- rep(0,nrow(annotated_variants))
trunc_col[trunc_indices] <- "YES"

##Produce a new column to tag variants with predicted truncating consequence an occuring in last 5% of transcript (as annotated by LOFTEE plugin)
end_trunc_indices <- grep("END_TRUNC", annotated_variants$INFO)
end_trunc_col <- rep(0,nrow(annotated_variants))
end_trunc_col[end_trunc_indices] <- "YES"

##Produce a new column to indicate variants with predicted truncating consequence that would be excluded due to being in last 5% of transcript
trunc_excl_indices <- which(trunc_col == "YES" & (GOF_col == "YES" | end_trunc_col == "YES"))
trunc_excl_col <- rep(0,nrow(annotated_variants))
trunc_excl_col[trunc_excl_indices] <- "YES"

##Produce a new column to tag variants with at least one pathogenic assertion in ClinVar database
clinvar_path_indices <- grep("pathogenic", annotated_variants$INFO)
clinvar_path_col <- rep(0,nrow(annotated_variants))
clinvar_path_col[clinvar_path_indices] <- "YES"

##Compile variant table with newly created columns to assist with assessment
annotated_variants_with_gt <- data.frame(annotated_variants[,1:11], 
                                 GOF_col, 
                                 recessive_col, 
                                 trunc_col, 
                                 end_trunc_col, 
                                 clinvar_path_col, 
                                 trunc_excl_col,
                                 annotated_variants[,12:ncol(annotated_variants)])

##Output variant table for further manual assessment
write.csv(annotated_variants_with_gt, "vp_SNV_indel_annotated_variants.csv")


#######################
##Structural variants##
#######################

##Filtered variant tables provided by NIHR BioResource Rare Diseases project exported to local server for analysis. Files used as follows:
#canvascalls_allsamples_gain.ann.filt.0.01.exons.txt
#canvascalls_allsamples_loss.ann.filt.0.01.exons.txt
#mantacalls_allsamples_bnd.ann.filt.0.01.txt
#mantacalls_allsamples_del.ann.filt.0.01.exons.txt
#mantacalls_allsamples_dup.ann.filt.0.01.exons.txt
#mantacalls_allsamples_ins.ann.filt.0.01.txt
#mantacalls_allsamples_inv.ann.filt.0.01.txt

##Subset variant table to only include samples of interest. Shown here for canvascalls_allsamples_gain.ann.filt.0.01.exons.txt. Repeated for all other files with only filename changing
fgrep -f samples_without_suffix.txt canvascalls_allsamples_gain.ann.filt.0.01.exons.txt > canvascalls_allsamples_gain.ann.filt.0.01.exons_subset.txt

##Subdivide variant table to reduce processing time
fgrep -f samples_without_suffix.txt canvascalls_allsamples_gain.ann.filt.0.01.exons_subset.txt | wc -l > number_variants.txt 
echo $(( (`cat number_variants.txt`/30) +1 )) > number_subset_files.txt #+1 to round up as bash doesn't use floating point
seq 1 `cat number_subset_files.txt` > header_file_numbers.txt
for i in `cat header_file_numbers.txt`; do
  head -1 canvascalls_allsamples_gain.ann.filt.0.01.exons.txt > variant_subset_${i}.txt
done
seq 2 30 `cat number_variants.txt` > col_1.txt
seq 31 30 `cat number_variants.txt` > col_2.txt
echo `cat number_variants.txt` >> col_2.txt
paste -d"," col_1.txt col_2.txt > variant_line_number_ranges.txt
paste variant_line_number_ranges.txt header_file_numbers.txt > variant_line_number_ranges_and_header_file_numbers.txt
IFS=$'\n'
for i in `cat variant_line_number_ranges_and_header_file_numbers.txt`;do
  COL1=`echo ${i} | awk '{ print $1 }'`
  COL2=`echo ${i} | awk '{ print $2 }'`
  sed -n ${COL1}p canvascalls_allsamples_gain.ann.filt.0.01.exons_subset.txt >> variant_subset_${COL2}.txt
done

##Resulting files named variant_subset_1.txt, variant subset_2.txt etc. retained in separate folder for each variant modality  

##R scripts for each variant modality (see below) used to identify structural variant calls fulfilling criteria. Initially used to produce a header line.
Rscript vp_SV_detection_canvas_gain.R 1
head -n 1 SV_canvas_gain_variant_table_1.txt > variants_combined.txt
rm SV_canvas_gain_variant_table_1.txt

##Subsequent execution and output to build variant table with annotation indicating if criteria fulfilled
for i in `cat header_file_numbers.txt`; do
  Rscript vp_SV_detection_canvas_gain.R ${i}
  sed '1d' SV_canvas_gain_variant_table_${i}.txt >> variants_combined.txt
done

##R scripts used 
#vp_SV_detection_canvas_gain.R
#vp_SV_detection_canvas_loss.R
#vp_SV_detection_manta_bnd.R
#vp_SV_detection_manta_del.R
#vp_SV_detection_manta_ins.R
#vp_SV_detection_manta_inv.R
#vp_SV_detection_manta_dup.R


#################################
##vp_SV_detection_canvas_gain.R##
#################################

args = commandArgs(trailingOnly=TRUE)
variant_table <- paste("variant_subset_",args[1],".txt",sep = "")
variants <- read.delim(variant_table,header = T)
elements <- read.csv("vp_genes.csv", header = T) ##File containing coordinates of cancer prediposition genes

##Adapt variant and element tables
elements$start <- as.numeric(elements$start)
elements$end <- as.numeric(elements$end)
elements$chrom <- as.character(elements$chrom)
variants$CANVAS_CHROM <- as.character(variants$CANVAS_CHROM)
variants$CANVAS_START <- as.numeric(variants$CANVAS_START)
variants$CANVAS_END <- as.numeric(variants$CANVAS_END)

##Detection of samples with predicted copy number gain involving entire gene
full_gain <- data.frame()
for (t in 1:nrow(variants)) {
  for (i in 1:nrow(elements)) {
    a <- isTRUE(
      (variants$CANVAS_CHROM[t] == elements$chrom[i]) & 
        ((variants$CANVAS_START[t] < elements$start[i]) & (variants$CANVAS_END[t] > elements$end[i]))
    )
    full_gain <- append(full_gain, a)
  }
}
full_gain <- unlist(full_gain)
full_gain <- split(full_gain, ceiling(seq_along(full_gain)/nrow(elements)))
full_gain <-t(data.frame(full_gain)) 
full_gain_hits <- data.frame()
for (i in 1:nrow(full_gain)) {
  b <- match(TRUE, full_gain[i,])
  if(is.na(b)){b <- 9999} ##9999 used as a row number that doesn't exist in the elements table
  full_gain_hits <- append(full_gain_hits, b) 
}
full_gain_hits <- unlist(full_gain_hits)
full_gain_hits <- data.frame(full_gain_hits)
full_gain_hits_IDs <- data.frame()
for (i in 1:nrow(full_gain_hits)) {
  c <- elements$identifier[full_gain_hits[i,1]] 
  d <- as.character(c)
  full_gain_hits_IDs <- append(full_gain_hits_IDs, d)
}
full_gain_hits_IDs <- cbind(full_gain_hits_IDs)
full_gain_hits_IDs <- unlist(full_gain_hits_IDs)

##Detection of samples with predicted copy number gain involving start of gene but not end
start_gain <- data.frame()
for (t in 1:nrow(variants)) {
  for (i in 1:nrow(elements)) {
    a <- isTRUE(
      (variants$CANVAS_CHROM[t] == elements$chrom[i]) & 
        ((variants$CANVAS_START[t] < elements$start[i]) & ((variants$CANVAS_END[t] > elements$start[i]) & (variants$CANVAS_END[t] < elements$end[i])))
    )
    start_gain <- append(start_gain, a)
  }
}
start_gain <- unlist(start_gain)
start_gain <- split(start_gain, ceiling(seq_along(start_gain)/nrow(elements)))
start_gain <-t(data.frame(start_gain)) 
start_gain_hits <- data.frame()
for (i in 1:nrow(start_gain)) {
  b <- match(TRUE, start_gain[i,])
  if(is.na(b)){b <- 9999} ##9999 used as a row number that doesn't exist in the elements table
  start_gain_hits <- append(start_gain_hits, b) 
}
start_gain_hits <- unlist(start_gain_hits)
start_gain_hits <- data.frame(start_gain_hits)
start_gain_hits_IDs <- data.frame()
for (i in 1:nrow(start_gain_hits)) {
  c <- elements$identifier[start_gain_hits[i,1]]
  c <- as.character(c)
  start_gain_hits_IDs <- append(start_gain_hits_IDs, c)
}
start_gain_hits_IDs <- cbind(start_gain_hits_IDs)
start_gain_hits_IDs <- unlist(start_gain_hits_IDs)

##Detection of samples with predicted copy number gain involving end of gene but not start
end_gain <- data.frame()
for (t in 1:nrow(variants)) {
  for (i in 1:nrow(elements)) {
    a <- isTRUE( 
      (variants$CANVAS_CHROM[t] == elements$chrom[i]) &     
        ((variants$CANVAS_END[t] > elements$end[i]) & ((variants$CANVAS_START[t] > elements$start[i]) & (variants$CANVAS_END[t] < elements$end[i])))
    )
    end_gain <- append(end_gain, a)
  }
}
end_gain <- unlist(end_gain)
end_gain <- split(end_gain, ceiling(seq_along(end_gain)/nrow(elements)))
end_gain <-t(data.frame(end_gain)) 
end_gain_hits <- data.frame()
for (i in 1:nrow(end_gain)) {
  b <- match(TRUE, end_gain[i,])
  if(is.na(b)){b <- 9999} ##9999 used as a row number that doesn't exist in the elements table
  end_gain_hits <- append(end_gain_hits, b)
}
end_gain_hits <- unlist(end_gain_hits)
end_gain_hits <- data.frame(end_gain_hits)
end_gain_hits_IDs <- data.frame()
for (i in 1:nrow(end_gain_hits)) {
  c <- elements$identifier[end_gain_hits[i,1]]
  c <- as.character(c)
  end_gain_hits_IDs <- append(end_gain_hits_IDs, c) 
}
end_gain_hits_IDs <- cbind(end_gain_hits_IDs)
end_gain_hits_IDs <- unlist(end_gain_hits_IDs)

##Produce final table
final_table <- data.frame(variants,full_gain_hits_IDs,start_gain_hits_IDs,end_gain_hits_IDs)

##Remove variants not fulfilling quality scores or filtering criteria
final_table <- final_table[which(final_table$CANVAS_QUAL > 29),]
final_table <- final_table[which(!is.na(final_table$full_gain_hits_IDs) | 
                                   !is.na(final_table$start_gain_hits_IDs) | 
                                   !is.na(final_table$end_gain_hits_IDs)),]


##Output annotated variant table
final_table_file <- paste("SV_canvas_gain_variant_table_",args[1],".txt",sep = "")
write.table(final_table, file = final_table_file, sep = "\t")


#################################
##vp_SV_detection_canvas_loss.R##
#################################

args = commandArgs(trailingOnly=TRUE)
variant_table <- paste("variant_subset_",args[1],".txt",sep = "")
variants <- read.delim(variant_table,header = T)
elements <- read.csv("vp_genes.csv", header = T)  ##File containing coordinates of cancer prediposition genes

##Adapt variant and element tables
elements$start <- as.numeric(elements$start)
elements$end <- as.numeric(elements$end)
elements$chrom <- as.character(elements$chrom)
variants$CANVAS_CHROM <- as.character(variants$CANVAS_CHROM)
variants$CANVAS_START <- as.numeric(variants$CANVAS_START)
variants$CANVAS_END <- as.numeric(variants$CANVAS_END)

##Detection of samples with predicted copy number loss involving entire gene
full_del <- data.frame()
for (t in 1:nrow(variants)) {
for (i in 1:nrow(elements)) {
  a <- isTRUE(
    (variants$CANVAS_CHROM[t] == elements$chrom[i]) & 
      ((variants$CANVAS_START[t] < elements$start[i]) & (variants$CANVAS_END[t] > elements$end[i]))
  )
  full_del <- append(full_del, a)
}
}
full_del <- unlist(full_del)
full_del <- split(full_del, ceiling(seq_along(full_del)/nrow(elements)))
full_del <-t(data.frame(full_del)) 
full_del_hits <- data.frame()
for (i in 1:nrow(full_del)) {
  b <- match(TRUE, full_del[i,])
  if(is.na(b)){b <- 9999} ##9999 used as a row number that doesn't exist in the elements table
  full_del_hits <- append(full_del_hits, b)
}
full_del_hits <- unlist(full_del_hits)
full_del_hits <- data.frame(full_del_hits)
full_del_hits_IDs <- data.frame()
for (i in 1:nrow(full_del_hits)) {
  c <- elements$identifier[full_del_hits[i,1]] 
  d <- as.character(c)
  full_del_hits_IDs <- append(full_del_hits_IDs, d)
}
full_del_hits_IDs <- cbind(full_del_hits_IDs)
full_del_hits_IDs <- unlist(full_del_hits_IDs)

##Detection of samples predicted copy number loss involving start of gene but not end
start_del <- data.frame()
for (t in 1:nrow(variants)) {
  for (i in 1:nrow(elements)) {
    a <- isTRUE(
      (variants$CANVAS_CHROM[t] == elements$chrom[i]) &  
           ((variants$CANVAS_START[t] < elements$start[i]) & ((variants$CANVAS_END[t] > elements$start[i]) & (variants$CANVAS_END[t] < elements$end[i])))
    )   
    start_del <- append(start_del, a)
  }
}
start_del <- unlist(start_del)
start_del <- split(start_del, ceiling(seq_along(start_del)/nrow(elements)))
start_del <-t(data.frame(start_del)) 
start_del_hits <- data.frame()
for (i in 1:nrow(start_del)) {
  b <- match(TRUE, start_del[i,])
  if(is.na(b)){b <- 9999} ##9999 used as a row number that doesn't exist in the elements table
  start_del_hits <- append(start_del_hits, b) 
}
start_del_hits <- unlist(start_del_hits)
start_del_hits <- data.frame(start_del_hits)
start_del_hits_IDs <- data.frame()
for (i in 1:nrow(start_del_hits)) { 
  c <- elements$identifier[start_del_hits[i,1]]
  c <- as.character(c)
  start_del_hits_IDs <- append(start_del_hits_IDs, c) 
}
start_del_hits_IDs <- cbind(start_del_hits_IDs)
start_del_hits_IDs <- unlist(start_del_hits_IDs)


##Detection of samples predicted copy number loss involving end of gene but not start
end_del <- data.frame()
for (t in 1:nrow(variants)) {
  for (i in 1:nrow(elements)) {
    a <- isTRUE( 
      (variants$CANVAS_CHROM[t] == elements$chrom[i]) &       
           ((variants$CANVAS_END[t] > elements$end[i]) & ((variants$CANVAS_START[t] > elements$start[i]) & (variants$CANVAS_END[t] < elements$end[i])))
    )
    end_del <- append(end_del, a)
  }
}
end_del <- unlist(end_del)
end_del <- split(end_del, ceiling(seq_along(end_del)/nrow(elements)))
end_del <-t(data.frame(end_del)) 
end_del_hits <- data.frame()
for (i in 1:nrow(end_del)) {
  b <- match(TRUE, end_del[i,])
  if(is.na(b)){b <- 9999} ##9999 used as a row number that doesn't exist in the elements table 
  end_del_hits <- append(end_del_hits, b) 
}
end_del_hits <- unlist(end_del_hits)
end_del_hits <- data.frame(end_del_hits)
end_del_hits_IDs <- data.frame()
for (i in 1:nrow(end_del_hits)) {  
  c <- elements$identifier[end_del_hits[i,1]]
  c <- as.character(c)
  end_del_hits_IDs <- append(end_del_hits_IDs, c)  
}
end_del_hits_IDs <- cbind(end_del_hits_IDs)
end_del_hits_IDs <- unlist(end_del_hits_IDs)

##Produce final table
final_table <- data.frame(variants,full_del_hits_IDs,start_del_hits_IDs,end_del_hits_IDs)

##Remove variants not fulfilling quality scores or filtering criteria
final_table <- final_table[which(final_table$CANVAS_QUAL > 29),]
final_table <- final_table[which(!is.na(final_table$full_del_hits_IDs) | 
                                   !is.na(final_table$start_del_hits_IDs) | 
                                   !is.na(final_table$end_del_hits_IDs)),]

##Output annotated variant table
final_table_file <- paste("SV_canvas_loss_variant_table_",args[1],".txt",sep = "")
write.table(final_table, file = final_table_file, sep = "\t")

###############################
##vp_SV_detection_manta_bnd.R##
###############################
args = commandArgs(trailingOnly=TRUE)
variant_table <- paste("variant_subset_",args[1],".txt",sep = "")
variants <- read.delim(variant_table,header = T)
elements <- read.csv("vp_genes.csv", header = T)

##Adapt variant and element tables
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


##Detection of samples with predicted translocations with breakpoint within gene
transloc <- data.frame()
for (t in 1:nrow(variants)) {
  for (i in 1:nrow(elements)) {
    a <- isTRUE(
      (variants$MANTA_CHROM[t] == elements$chrom[i]) &
        (((variants$MANTA_MIN_START[t] > elements$start[i]) & (variants$MANTA_MAX_START[t] < elements$end[i])) || 
           ((variants$MANTA_MIN_END[t] > elements$start[i]) & (variants$MANTA_MAX_END[t] < elements$end[i])))
    )
    transloc <- append(transloc, a)  
  }
}
transloc <- unlist(transloc)
transloc <- split(transloc, ceiling(seq_along(transloc)/nrow(elements)))
transloc <-t(data.frame(transloc)) 
transloc_hits <- data.frame()
for (i in 1:nrow(transloc)) {
  b <- match(TRUE, transloc[i,])
  if(is.na(b)){b <- 9999} ##9999 used as a row number that doesn't exist in the elements table
  transloc_hits <- append(transloc_hits, b)
}
transloc_hits <- unlist(transloc_hits)
transloc_hits <- data.frame(transloc_hits)
transloc_hits_IDs <- data.frame()
for (i in 1:nrow(transloc_hits)) {  
  c <- elements$identifier[transloc_hits[i,1]]
  c <- as.character(c)
  transloc_hits_IDs <- append(transloc_hits_IDs, c) 
}
transloc_hits_IDs <- cbind(transloc_hits_IDs)
transloc_hits_IDs <- unlist(transloc_hits_IDs)

##Produce final table
final_table <- data.frame(variants,transloc_hits_IDs)

##Remove variants not fulfilling quality scores or filtering criteria
final_table <- final_table[which(final_table$MANTA_GQ > 29),]
final_table <- final_table[which(!is.na(final_table$transloc_hits_IDs)),]

##Output annotated variant table
final_table_file <- paste("SV_manta_bnd_not_only_exons_variant_table_",args[1],".txt",sep = "")
write.table(final_table, file = final_table_file, sep = "\t")

###############################
##vp_SV_detection_manta_del.R##
###############################

args = commandArgs(trailingOnly=TRUE)
variant_table <- paste("variant_subset_",args[1],".txt",sep = "")
variants <- read.delim(variant_table,header = T)
elements <- read.csv("vp_genes.csv", header = T)

##Adapt variant and element tables
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


##Detection of samples with predicted deletion involving entire gene
full_del <- data.frame()
for (t in 1:nrow(variants)) {
  for (i in 1:nrow(elements)) {
    a <- isTRUE(  
      (variants$MANTA_CHROM[t] == elements$chrom[i]) &     
        ((variants$MANTA_MAX_START[t] < elements$start[t]) & (variants$MANTA_MIN_END[t] > elements$end[t]))  
    )   
    full_del <- append(full_del, a)    
  }
}
full_del <- unlist(full_del)
full_del <- split(full_del, ceiling(seq_along(full_del)/nrow(elements)))
full_del <-t(data.frame(full_del)) 
full_del_hits <- data.frame()
for (i in 1:nrow(full_del)) {
  b <- match(TRUE, full_del[i,])
  if(is.na(b)){b <- 9999} ##9999 used as a row number that doesn't exist in the elements table
  full_del_hits <- append(full_del_hits, b)
  
}
full_del_hits <- unlist(full_del_hits)
full_del_hits <- data.frame(full_del_hits)
full_del_hits_IDs <- data.frame()
for (i in 1:nrow(full_del_hits)) {  
  c <- elements$identifier[full_del_hits[i,1]] 
  d <- as.character(c)
  full_del_hits_IDs <- append(full_del_hits_IDs, d)  
}
full_del_hits_IDs <- cbind(full_del_hits_IDs)
full_del_hits_IDs <- unlist(full_del_hits_IDs)

##Detection of samples with predicted deletion involving start of gene but not end
start_del <- data.frame()
for (t in 1:nrow(variants)) {
  for (i in 1:nrow(elements)) {    
    a <- isTRUE(      
      (variants$MANTA_CHROM[t] == elements$chrom[i]) &        
        ((variants$MANTA_MAX_START[t] < elements$start[i]) & ((variants$MANTA_MIN_END[t] > elements$start[i]) & (variants$MANTA_MAX_END[t] < elements$end[i]))) 
    ) 
    start_del <- append(start_del, a)   
  }
}
start_del <- unlist(start_del)
start_del <- split(start_del, ceiling(seq_along(start_del)/nrow(elements)))
start_del <-t(data.frame(start_del)) 
start_del_hits <- data.frame()
for (i in 1:nrow(start_del)) {
  b <- match(TRUE, start_del[i,])
  if(is.na(b)){b <- 9999} ##9999 used as a row number that doesn't exist in the elements table
  start_del_hits <- append(start_del_hits, b)  
}
start_del_hits <- unlist(start_del_hits)
start_del_hits <- data.frame(start_del_hits)
start_del_hits_IDs <- data.frame()
for (i in 1:nrow(start_del_hits)) {  
  c <- elements$identifier[start_del_hits[i,1]]
  c <- as.character(c)
  start_del_hits_IDs <- append(start_del_hits_IDs, c)  
}

start_del_hits_IDs <- cbind(start_del_hits_IDs)
start_del_hits_IDs <- unlist(start_del_hits_IDs)

##Detection of samples with predicted deletion involving end of gene but not start
end_del <- data.frame()
for (t in 1:nrow(variants)) {
  for (i in 1:nrow(elements)) {
    a <- isTRUE(   
      (variants$MANTA_CHROM[t] == elements$chrom[i]) &        
        ((variants$MANTA_MIN_END[t] > elements$end[i]) & ((variants$MANTA_MIN_START[t] > elements$start[i]) & (variants$MANTA_MAX_END[t] < elements$end[i])))
    )    
    end_del <- append(end_del, a)    
  }
}
end_del <- unlist(end_del)
end_del <- split(end_del, ceiling(seq_along(end_del)/nrow(elements)))
end_del <-t(data.frame(end_del)) 
end_del_hits <- data.frame()
for (i in 1:nrow(end_del)) {
  b <- match(TRUE, end_del[i,])
  if(is.na(b)){b <- 9999} ##9999 used as a row number that doesn't exist in the elements table
  end_del_hits <- append(end_del_hits, b)
}
end_del_hits <- unlist(end_del_hits)
end_del_hits <- data.frame(end_del_hits)
end_del_hits_IDs <- data.frame()
for (i in 1:nrow(end_del_hits)) {  
  c <- elements$identifier[end_del_hits[i,1]]
  c <- as.character(c)
  end_del_hits_IDs <- append(end_del_hits_IDs, c)  
}
end_del_hits_IDs <- cbind(end_del_hits_IDs)
end_del_hits_IDs <- unlist(end_del_hits_IDs)


##Detection of samples with predicted deletion within gene
intra_del <- data.frame()
for (t in 1:nrow(variants)) {
  for (i in 1:nrow(elements)) {    
    a <- isTRUE(      
      (variants$MANTA_CHROM[t] == elements$chrom[i]) &         
        ((variants$MANTA_MIN_START[t] > elements$start[i]) & (variants$MANTA_MAX_END[t] < elements$end[i]))  
    )
    intra_del <- append(intra_del, a)
  }
}
intra_del <- unlist(intra_del)
intra_del <- split(intra_del, ceiling(seq_along(intra_del)/nrow(elements)))
intra_del <-t(data.frame(intra_del)) 
intra_del_hits <- data.frame()
for (i in 1:nrow(intra_del)) {
  b <- match(TRUE, intra_del[i,])
  if(is.na(b)){b <- 9999} ##9999 used as a row number that doesn't exist in the elements table
  intra_del_hits <- append(intra_del_hits, b)
}
intra_del_hits <- unlist(intra_del_hits)
intra_del_hits <- data.frame(intra_del_hits)
intra_del_hits_IDs <- data.frame()
for (i in 1:nrow(intra_del_hits)) {  
  c <- elements$identifier[intra_del_hits[i,1]]
  c <- as.character(c)
  intra_del_hits_IDs <- append(intra_del_hits_IDs, c)  
}
intra_del_hits_IDs <- cbind(intra_del_hits_IDs)
intra_del_hits_IDs <- unlist(intra_del_hits_IDs)

##Produce final table
final_table <- data.frame(variants,full_del_hits_IDs,start_del_hits_IDs,end_del_hits_IDs,intra_del_hits_IDs)

##Remove variants not fulfilling quality scores or filtering criteria
final_table <- final_table[which(final_table$MANTA_GQ > 29),]
final_table <- final_table[which(!is.na(final_table$full_del_hits_IDs) | 
                                   !is.na(final_table$start_del_hits_IDs) | 
                                   !is.na(final_table$end_del_hits_IDs) | 
                                   !is.na(final_table$intra_del_hits_IDs)),]
##Output annotated variant table
final_table_file <- paste("SV_manta_del_variant_table_",args[1],".txt",sep = "")
write.table(final_table, file = final_table_file, sep = "\t")

###############################
##vp_SV_detection_manta_ins.R##
###############################

args = commandArgs(trailingOnly=TRUE)
variant_table <- paste("variant_subset_",args[1],".txt",sep = "")
variants <- read.delim(variant_table,header = T)
elements <- read.csv("vp_genes.csv", header = T)

##Adapt variant and element tables
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

##Detection of samples with predicted insertion with breakpoints containing entire gene
full_ins <- data.frame()
for (t in 1:nrow(variants)) {
  for (i in 1:nrow(elements)) {   
    a <- isTRUE(      
      (variants$MANTA_CHROM[t] == elements$chrom[i]) &         
        ((variants$MANTA_MAX_START[t] < elements$start[t]) & (variants$MANTA_MIN_END[t] > elements$end[t]))      
    )    
    full_ins <- append(full_ins, a)    
  }
}
full_ins <- unlist(full_ins)
full_ins <- split(full_ins, ceiling(seq_along(full_ins)/nrow(elements)))
full_ins <-t(data.frame(full_ins)) 
full_ins_hits <- data.frame()
for (i in 1:nrow(full_ins)) {  
  b <- match(TRUE, full_ins[i,])
  if(is.na(b)){b <- 9999} ##9999 used as a row number that doesn't exist in the elements table
  full_ins_hits <- append(full_ins_hits, b)  
}
full_ins_hits <- unlist(full_ins_hits)
full_ins_hits <- data.frame(full_ins_hits)
full_ins_hits_IDs <- data.frame()
for (i in 1:nrow(full_ins_hits)) {  
  c <- elements$identifier[full_ins_hits[i,1]] 
  d <- as.character(c)
  full_ins_hits_IDs <- append(full_ins_hits_IDs, d)  
}
full_ins_hits_IDs <- cbind(full_ins_hits_IDs)
full_ins_hits_IDs <- unlist(full_ins_hits_IDs)

##Detection of samples with predicted insertion with breakpoints containing start but not end of gene
start_ins <- data.frame()
for (t in 1:nrow(variants)) {
  for (i in 1:nrow(elements)) {    
    a <- isTRUE(      
      (variants$MANTA_CHROM[t] == elements$chrom[i]) &         
        ((variants$MANTA_MAX_START[t] < elements$start[i]) & ((variants$MANTA_MIN_END[t] > elements$start[i]) & (variants$MANTA_MAX_END[t] < elements$end[i]))) 
    )   
    start_ins <- append(start_ins, a)   
  }
}
start_ins <- unlist(start_ins)
start_ins <- split(start_ins, ceiling(seq_along(start_ins)/nrow(elements)))
start_ins <-t(data.frame(start_ins)) 
start_ins_hits <- data.frame()
for (i in 1:nrow(start_ins)) {  
  b <- match(TRUE, start_ins[i,])
  if(is.na(b)){b <- 9999} ##9999 used as a row number that doesn't exist in the elements table
  start_ins_hits <- append(start_ins_hits, b)
}
start_ins_hits <- unlist(start_ins_hits)
start_ins_hits <- data.frame(start_ins_hits)
start_ins_hits_IDs <- data.frame()
for (i in 1:nrow(start_ins_hits)) { 
  c <- elements$identifier[start_ins_hits[i,1]]
  c <- as.character(c)
  start_ins_hits_IDs <- append(start_ins_hits_IDs, c) 
}
start_ins_hits_IDs <- cbind(start_ins_hits_IDs)
start_ins_hits_IDs <- unlist(start_ins_hits_IDs)


##Detection of samples with predicted insertion with breakpoints containing end but not start of gene
end_ins <- data.frame()
for (t in 1:nrow(variants)) {
  for (i in 1:nrow(elements)) {    
    a <- isTRUE(      
      (variants$MANTA_CHROM[t] == elements$chrom[i]) &        
        ((variants$MANTA_MIN_END[t] > elements$end[i]) & ((variants$MANTA_MIN_START[t] > elements$start[i]) & (variants$MANTA_MAX_END[t] < elements$end[i])))
    )    
    end_ins <- append(end_ins, a)    
  }
}
end_ins <- unlist(end_ins)
end_ins <- split(end_ins, ceiling(seq_along(end_ins)/nrow(elements)))
end_ins <-t(data.frame(end_ins)) 
end_ins_hits <- data.frame()
for (i in 1:nrow(end_ins)) {  
  b <- match(TRUE, end_ins[i,])
  if(is.na(b)){b <- 9999} ##9999 used as a row number that doesn't exist in the elements table
  end_ins_hits <- append(end_ins_hits, b)  
}
end_ins_hits <- unlist(end_ins_hits)
end_ins_hits <- data.frame(end_ins_hits)
end_ins_hits_IDs <- data.frame()
for (i in 1:nrow(end_ins_hits)) {  
  c <- elements$identifier[end_ins_hits[i,1]]
  c <- as.character(c)
  end_ins_hits_IDs <- append(end_ins_hits_IDs, c)  
}
end_ins_hits_IDs <- cbind(end_ins_hits_IDs)
end_ins_hits_IDs <- unlist(end_ins_hits_IDs)


##Detection of samples with predicted insertion with breakpoints within gene
intra_ins <- data.frame()
for (t in 1:nrow(variants)) {
  for (i in 1:nrow(elements)) {  
    a <- isTRUE(    
      (variants$MANTA_CHROM[t] == elements$chrom[i]) &        
        ((variants$MANTA_MIN_START[t] > elements$start[i]) & (variants$MANTA_MAX_END[t] < elements$end[i])) 
    )  
    intra_ins <- append(intra_ins, a)  
  }
}
intra_ins <- unlist(intra_ins)
intra_ins <- split(intra_ins, ceiling(seq_along(intra_ins)/nrow(elements)))
intra_ins <-t(data.frame(intra_ins)) 
intra_ins_hits <- data.frame()
for (i in 1:nrow(intra_ins)) { 
  b <- match(TRUE, intra_ins[i,])
  if(is.na(b)){b <- 9999} ##9999 used as a row number that doesn't exist in the elements table
  intra_ins_hits <- append(intra_ins_hits, b) 
}
intra_ins_hits <- unlist(intra_ins_hits)
intra_ins_hits <- data.frame(intra_ins_hits)
intra_ins_hits_IDs <- data.frame()
for (i in 1:nrow(intra_ins_hits)) {  
  c <- elements$identifier[intra_ins_hits[i,1]]
  c <- as.character(c)
  intra_ins_hits_IDs <- append(intra_ins_hits_IDs, c) 
}
intra_ins_hits_IDs <- cbind(intra_ins_hits_IDs)
intra_ins_hits_IDs <- unlist(intra_ins_hits_IDs)

##Produce final table
final_table <- data.frame(variants,full_ins_hits_IDs,start_ins_hits_IDs,end_ins_hits_IDs,intra_ins_hits_IDs)

##Remove variants not fulfilling quality scores or filtering criteria
final_table <- final_table[which(final_table$MANTA_GQ > 29),]
final_table <- final_table[which(!is.na(final_table$full_ins_hits_IDs) | 
                                   !is.na(final_table$start_ins_hits_IDs) | 
                                   !is.na(final_table$end_ins_hits_IDs) | 
                                   !is.na(final_table$intra_ins_hits_IDs)),]
                                   
##Output annotated variant table
final_table_file <- paste("SV_manta_ins_variant_table_",args[1],".txt",sep = "")
write.table(final_table, file = final_table_file, sep = "\t")

###############################
##vp_SV_detection_manta_inv.R##
###############################


args = commandArgs(trailingOnly=TRUE)
variant_table <- paste("variant_subset_",args[1],".txt",sep = "")
variants <- read.delim(variant_table,header = T)
elements <- read.csv("vp_genes.csv", header = T)

##Adapt variant and element tables
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


##Detection of samples with predicted inversion with breakpoint within gene
inv <- data.frame()
for (t in 1:nrow(variants)) {
  for (i in 1:nrow(elements)) {    
    a <- isTRUE(      
      (variants$MANTA_CHROM[t] == elements$chrom[i]) &      
        (((variants$MANTA_MIN_START[t] > elements$start[i]) & (variants$MANTA_MAX_START[t] < elements$end[i])) || 
           ((variants$MANTA_MIN_END[t] > elements$start[i]) & (variants$MANTA_MAX_END[t] < elements$end[i])))
    )    
    inv <- append(inv, a)    
  }
}
inv <- unlist(inv)
inv <- split(inv, ceiling(seq_along(inv)/nrow(elements)))
inv <-t(data.frame(inv)) 
inv_hits <- data.frame()
for (i in 1:nrow(inv)) {  
  b <- match(TRUE, inv[i,])
  if(is.na(b)){b <- 9999} ##9999 used as a row number that doesn't exist in the elements table
  inv_hits <- append(inv_hits, b) 
}
inv_hits <- unlist(inv_hits)
inv_hits <- data.frame(inv_hits)
inv_hits_IDs <- data.frame()
for (i in 1:nrow(inv_hits)) {  
  c <- elements$identifier[inv_hits[i,1]]
  c <- as.character(c)
  inv_hits_IDs <- append(inv_hits_IDs, c)  
}
inv_hits_IDs <- cbind(inv_hits_IDs)
inv_hits_IDs <- unlist(inv_hits_IDs)

##Produce final table
final_table <- data.frame(variants,inv_hits_IDs)

##Remove variants not fulfilling quality scores or filtering criteria
final_table <- final_table[which(final_table$MANTA_GQ > 29),]
final_table <- final_table[which(!is.na(final_table$inv_hits_IDs)),]

##Output annotated variant table
final_table_file <- paste("SV_manta_inv_not_only_exons_variant_table_",args[1],".txt",sep = "")
write.table(final_table, file = final_table_file, sep = "\t")

###############################
##vp_SV_detection_manta_dup.R##
###############################

args = commandArgs(trailingOnly=TRUE)
variant_table <- paste("variant_subset_",args[1],".txt",sep = "")
variants <- read.delim(variant_table,header = T)
elements <- read.csv("vp_genes.csv", header = T)

##Adapt variant and element tables
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

##Detection of samples with predicted duplication with breakpoints containing entire gene
full_dup <- data.frame()
for (t in 1:nrow(variants)) {
  for (i in 1:nrow(elements)) {    
    a <- isTRUE(      
      (variants$MANTA_CHROM[t] == elements$chrom[i]) &         
        ((variants$MANTA_MAX_START[t] < elements$start[t]) & (variants$MANTA_MIN_END[t] > elements$end[t]))      
    )    
    full_dup <- append(full_dup, a)    
  }
}
full_dup <- unlist(full_dup)
full_dup <- split(full_dup, ceiling(seq_along(full_dup)/nrow(elements)))
full_dup <-t(data.frame(full_dup)) 
full_dup_hits <- data.frame()
for (i in 1:nrow(full_dup)) {
  b <- match(TRUE, full_dup[i,])
  if(is.na(b)){b <- 9999} ##9999 used as a row number that doesn't exist in the elements table
  full_dup_hits <- append(full_dup_hits, b)
}
full_dup_hits <- unlist(full_dup_hits)
full_dup_hits <- data.frame(full_dup_hits)
full_dup_hits_IDs <- data.frame()
for (i in 1:nrow(full_dup_hits)) {  
  c <- elements$identifier[full_dup_hits[i,1]] 
  d <- as.character(c)
  full_dup_hits_IDs <- append(full_dup_hits_IDs, d)  
}
full_dup_hits_IDs <- cbind(full_dup_hits_IDs)
full_dup_hits_IDs <- unlist(full_dup_hits_IDs)


##Detection of samples with predicted duplication with breakpoints containing start but not end of gene
start_dup <- data.frame()
for (t in 1:nrow(variants)) {
  for (i in 1:nrow(elements)) {    
    a <- isTRUE(      
      (variants$MANTA_CHROM[t] == elements$chrom[i]) &         
        ((variants$MANTA_MAX_START[t] < elements$start[i]) & ((variants$MANTA_MIN_END[t] > elements$start[i]) & (variants$MANTA_MAX_END[t] < elements$end[i])))  
    )    
    start_dup <- append(start_dup, a)    
  }
}
start_dup <- unlist(start_dup)
start_dup <- split(start_dup, ceiling(seq_along(start_dup)/nrow(elements)))
start_dup <-t(data.frame(start_dup)) 
start_dup_hits <- data.frame()
for (i in 1:nrow(start_dup)) {  
  b <- match(TRUE, start_dup[i,])
  if(is.na(b)){b <- 9999} ##9999 used as a row number that doesn't exist in the elements table
  start_dup_hits <- append(start_dup_hits, b)  
}
start_dup_hits <- unlist(start_dup_hits)
start_dup_hits <- data.frame(start_dup_hits)
start_dup_hits_IDs <- data.frame()
for (i in 1:nrow(start_dup_hits)) {  
  c <- elements$identifier[start_dup_hits[i,1]]
  c <- as.character(c)
  start_dup_hits_IDs <- append(start_dup_hits_IDs, c)  
}
start_dup_hits_IDs <- cbind(start_dup_hits_IDs)
start_dup_hits_IDs <- unlist(start_dup_hits_IDs)

##Detection of samples with predicted duplication with breakpoints containing end but not start of gene
end_dup <- data.frame()
for (t in 1:nrow(variants)) {
  for (i in 1:nrow(elements)) {   
    a <- isTRUE(     
      (variants$MANTA_CHROM[t] == elements$chrom[i]) &        
        ((variants$MANTA_MIN_END[t] > elements$end[i]) & ((variants$MANTA_MIN_START[t] > elements$start[i]) & (variants$MANTA_MAX_END[t] < elements$end[i]))) 
    ) 
    end_dup <- append(end_dup, a)    
  }
}
end_dup <- unlist(end_dup)
end_dup <- split(end_dup, ceiling(seq_along(end_dup)/nrow(elements)))
end_dup <-t(data.frame(end_dup)) 
end_dup_hits <- data.frame()
for (i in 1:nrow(end_dup)) {  
  b <- match(TRUE, end_dup[i,])
  if(is.na(b)){b <- 9999} ##9999 used as a row number that doesn't exist in the elements table
  end_dup_hits <- append(end_dup_hits, b)  
}
end_dup_hits <- unlist(end_dup_hits)
end_dup_hits <- data.frame(end_dup_hits)
end_dup_hits_IDs <- data.frame()
for (i in 1:nrow(end_dup_hits)) {  
  c <- elements$identifier[end_dup_hits[i,1]]
  c <- as.character(c)
  end_dup_hits_IDs <- append(end_dup_hits_IDs, c)  
}
end_dup_hits_IDs <- cbind(end_dup_hits_IDs)
end_dup_hits_IDs <- unlist(end_dup_hits_IDs)

##Detection of samples with predicted duplication with breakpoints within gene
intra_dup <- data.frame()
for (t in 1:nrow(variants)) {
  for (i in 1:nrow(elements)) {    
    a <- isTRUE(      
      (variants$MANTA_CHROM[t] == elements$chrom[i]) &      
        ((variants$MANTA_MIN_START[t] > elements$start[i]) & (variants$MANTA_MAX_END[t] < elements$end[i])) 
    )
    intra_dup <- append(intra_dup, a)
  }
}
intra_dup <- unlist(intra_dup)
intra_dup <- split(intra_dup, ceiling(seq_along(intra_dup)/nrow(elements)))
intra_dup <-t(data.frame(intra_dup)) 
intra_dup_hits <- data.frame()
for (i in 1:nrow(intra_dup)) { 
  b <- match(TRUE, intra_dup[i,])
  if(is.na(b)){b <- 9999} ##9999 used as a row number that doesn't exist in the elements table
  intra_dup_hits <- append(intra_dup_hits, b)  
}
intra_dup_hits <- unlist(intra_dup_hits)
intra_dup_hits <- data.frame(intra_dup_hits)
intra_dup_hits_IDs <- data.frame()
for (i in 1:nrow(intra_dup_hits)) {  
  c <- elements$identifier[intra_dup_hits[i,1]]
  c <- as.character(c)
  intra_dup_hits_IDs <- append(intra_dup_hits_IDs, c)  
}
intra_dup_hits_IDs <- cbind(intra_dup_hits_IDs)
intra_dup_hits_IDs <- unlist(intra_dup_hits_IDs)

##Produce final table
final_table <- data.frame(variants,full_dup_hits_IDs,start_dup_hits_IDs,end_dup_hits_IDs,intra_dup_hits_IDs)

##Remove variants not fulfilling quality scores or filtering criteria
final_table <- final_table[which(final_table$MANTA_GQ > 29),]
final_table <- final_table[which(!is.na(final_table$full_dup_hits_IDs) | 
                                   !is.na(final_table$start_dup_hits_IDs) | 
                                   !is.na(final_table$end_dup_hits_IDs) | 
                                   !is.na(final_tableintra_dup_hits_IDs),]

##Output annotated variant table
final_table_file <- paste("SV_manta_dup_variant_table_",args[1],".txt",sep = "")
write.table(final_table, file = final_table_file, sep = "\t")



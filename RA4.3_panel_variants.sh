############################################################################################
##Filtering of VCF files containing variant calls from Illumina TruSight Cancer panel data##
############################################################################################

##Produce list of VCF files for filtering containing variant calls from Illumina TruSight Cancer panel data
ls *.vcf > samples.txt

for i in `cat samples.txt`; do
	bgzip ${i}
	tabix ${i}.gz
done

sed -r 's/.{4}$//' samples.txt > 1_samples.txt
rm samples.txt
mv 1_samples.txt samples.txt

##Filter  VCF files using BED file corresponding to coordinates of interest
for i in `cat samples.txt`; do
  bcftools view -R virtual_panel_exons_sorted.bed -o ${i}_filtered.vcf.gz -O z ${i}.vcf.gz ##BED file containing coordi
done

##Merge filtered  VCF files
for i in `cat samples.txt`; do
  gunzip ${i}_filtered.vcf.gz
done

ls *_filtered.vcf > filtered_samples.txt

for i in `cat filtered_samples.txt`; do
	vcf-sort -c ${i} > sorted_${i}
	bgzip sorted_${i}
	tabix sorted_${i}.gz
done

ls sorted_*_filtered.vcf.gz | tr "\n" " " > files_to_merge.txt
vcf-merge `cat files_to_merge.txt` | bgzip -c > vp_panels_merged.vcf.gz

##Filter merged VCF based on quality parameters (insert missing genotypes where criteria not fulfilled)
bcftools filter -e 'FMT/DP<10 || FMT/GQ<30 || (FMT/AD[1]/(FMT/AD[0] + FMT/AD[1]) <0.33)' -S . -o vp_panels_merged_filtered.vcf.gz -O z vp_panels_merged.vcf.gz

##Remove INFO fields from BioResource annotated files
bcftools annotate -x INFO vp_panels_merged_filtered.vcf.gz -o vp_SNV_indel.vcf

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

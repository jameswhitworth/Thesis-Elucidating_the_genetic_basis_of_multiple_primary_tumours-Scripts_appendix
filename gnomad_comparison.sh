###############################################
##Filtering and annotation of gnomad variants##
###############################################

##Filter gnomad VCF files for truncating variants in genes of interest (based on canonical transcripts)
/home/jww39/ensembl-vep/filter_vep -i gnomad.genomes.r2.0.2.sites.coding_only.chr1-22.vcf -o vp_gnomad.vcf -filter "(Feature in vp_transcripts.txt) and (IMPACT is HIGH)" --only_matched

##Remove INFO fields
bgzip vp_gnomad.vcf
tabix vp_gnomad.vcf.gz
bcftools annotate -x INFO vp_gnomad.vcf.gz -o vp_SNV_indel.vcf

##Re-annotate with variant effect predictor and filter
/home/jww39/ensembl-vep/vep -i vp_gnomad.vcf -offline --assembly GRCh37 -o VEP_out_vp_gnomad.vcf --pick --pick_order canonical --everything --vcf --plugin LoF --plugin CADD,/home/jww39/.vep/Plugins/whole_genome_SNVs.tsv.gz,/home/jww39/.vep/Plugins/InDels.tsv.gz --plugin ExAC,/home/jww39/.vep/Plugins/ExAC.r0.3.sites.vep.vcf.gz

/home/jww39/ensembl-vep/filter_vep -i VEP_out_vp_gnomad.vcf -o VEP_out_filtered_vp_gnomad.vcf -filter "(Feature in vp_transcripts.txt) and (IMPACT is HIGH)" --only_matched

/home/jww39/ensembl-vep/filter_vep -i VEP_out_filtered_vp_gnomad.vcf -o VEP_out_double_filtered_vp_gnomad.vcf -filter "(ExAC_AF < 0.01 or not ExAC_AF) and (AF < 0.01 or not AF)" --only_matched

##Adapt VCF file for reading into R
sed -i 's/#CHROM/CHROM/g' VEP_out_double_filtered_vp_gnomad.vcf
sed '/^#/ d' VEP_out_double_filtered_vp_gnomad.vcf > VEP_out_double_filtered_vp_gnomad_without_header.vcf

################################################################
##R script to annotate variant table to assist with exclusions##
################################################################

##Read in necessary files
variants <- read.table("VEP_out_double_filtered_vp_gnomad_without_header.vcf", header =T) ##Filtered variant table
GOF_genes <- readLines("GOF_genes.txt") ##List of proto-oncogne cancer predisposition genes
recessive_genes <- readLines("recessive_genes.txt") ##List of recessive cancer predisposition genes


##Produce a new column to tag variants occuring in proto-oncogene cancer predisposition genes
GOF_indices <- data.frame()

for (i in 1:length(GOF_genes)) {
  a <- grep(GOF_genes[i], variants$INFO)
  GOF_indices <- append(GOF_indices, a)
}

GOF_col <- rep(0,nrow(variants))
GOF_col[unlist(GOF_indices)] <- "YES"

##Produce a new column to tag variants occuring in recessive cancer predisposition genes
recessive_indices <- data.frame()

for (i in 1:length(recessive_genes)) {
  a <- grep(recessive_genes[i], variants$INFO)
  recessive_indices <- append(recessive_indices, a)
}

recessive_col <- rep(0,nrow(variants))
recessive_col[unlist(recessive_indices)] <- "YES"

##Produce a new column to tag variants with predicted truncating consequence
trunc_indices <- grep("\\|HIGH\\|", variants$INFO)
trunc_col <- rep(0,nrow(variants))
trunc_col[trunc_indices] <- "YES"

##Produce a new column to tag variants with predicted truncating consequence an occuring in last 5% of transcript (as annotated by LOFTEE plugin)
end_trunc_indices <- grep("END_TRUNC", variants$INFO)
end_trunc_col <- rep(0,nrow(variants))
end_trunc_col[end_trunc_indices] <- "YES"

##Produce a new column to indicate variants with predicted truncating consequence that would be excluded due to being in last 5% of transcript
trunc_excl_indices <- which(trunc_col == "YES" & (GOF_col == "YES" | end_trunc_col == "YES"))
trunc_excl_col <- rep(0,nrow(variants))
trunc_excl_col[trunc_excl_indices] <- "YES"

##Compile variant table with newly created columns to assist with assessment
variants_with_gt <- data.frame(variants,
                               trunc_col,
                                 GOF_col, 
                                 recessive_col, 
                                 end_trunc_col, 
                                 trunc_excl_col)

##Output variant table for further manual assessment
write.csv(variants_with_gt, "vp_gnomad.csv")

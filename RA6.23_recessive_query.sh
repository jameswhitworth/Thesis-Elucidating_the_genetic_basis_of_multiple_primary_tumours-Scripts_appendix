##############################################################################################################################################################
##Filtering of merged VCF files compiled by NIHR BioResource Rare Diseases project and stored on University of Cambridge high performance computeing cluster##
##############################################################################################################################################################

##Produce list of per chromosome files for filtering
ls /scratch/WGS10K/data/release/20170614-A/merged-vcf/no_hgmd/*.vcf.gz > 1_chromosome_files.txt
sed -r 's/.{7}$//' 1_chromosome_files.txt > 2_chromosome_files.txt
sed -r 's/.{59}//' 2_chromosome_files.txt > chromosome_files.txt
rm 1_chromosome_files.txt
rm 2_chromosome_files.txt

##Filter per chromosome VCF files using BED file corresponding to coordinates and samples of interest
for i in `cat chromosome_files.txt`; do
  bcftools view -R nexterarapidcapture_exome_targetedregions_v1.2_hg19.bed -s R014798,R014941,R014835,R014813 -e 'FILTER!="PASS"' -o ${i}_1981.vcf -O v /scratch/WGS10K/data/release/20170614-A/merged-vcf/no_hgmd/${i}.vcf.gz
done

##Merge per chromosome VCF files
ls *_1981.vcf > filtered_samples_for_bgzip_and_tabix.txt

for i in `cat filtered_samples_for_bgzip_and_tabix.txt`; do
	vcf-sort -c ${i} > sorted_${i}
	bgzip sorted_${i}
	tabix sorted_${i}.gz
done

rm filtered_samples_for_bgzip_and_tabix.txt

ls *.vcf.gz | tr "\n" " " > files_to_merge.txt

bcftools concat -Oz `cat files_to_merge.txt` -o merged.vcf.gz
tabix merged.vcf.gz

##Filter based on quality parameters (insert missing genotypes)
bcftools filter -e 'FMT/DP<10 || FMT/GQ<20 || (FMT/AD[1]/(FMT/AD[0] + FMT/AD[1]) <0.3) || (FMT/AD[1]/(FMT/AD[0] + FMT/AD[1] + FMT/AD[2]) <0.3) || (FMT/AD[2]/(FMT/AD[0] + FMT/AD[1] + FMT/AD[2]) <0.3)' -S . -o merged_filtered.vcf.gz -O z merged.vcf.gz

##Produce per individual individual VCF files 
vcftools --gzvcf merged_filtered.vcf.gz --out R014798_A --recode --recode-INFO-all --keep R014798_A_filename.txt
vcftools --gzvcf merged_filtered.vcf.gz --out R014941_A --recode --recode-INFO-all --keep R014941_A_filename.txt
vcftools --gzvcf merged_filtered.vcf.gz --out R014835_A --recode --recode-INFO-all --keep R014835_A_filename.txt
vcftools --gzvcf merged_filtered.vcf.gz --out R014813_A --recode --recode-INFO-all --keep R014813_A_filename.txt

##Extract homogygous variants from offspring for homozygous variant based analysis
bgzip R014798_A.recode.vcf
tabiz R014798_A.recode.vcf.gz
bcftools view -g hom -o R014798_A.recode_homs.vcf -O v R014798_A.recode.vcf.gz

bgzip R014941_A.recode.vcf
tabiz R014941_A.recode.vcf.gz
bcftools view -g hom -o R014941_A.recode_homs.vcf -O v R014941_A.recode.vcf.gz

##Extract heterozygous variants from parents for homozygous variant and compund heterozygote based analysis
bgzip R014835_A.recode.vcf
tabiz R014835_A.recode.vcf.gz
bcftools view -g het -o R014835_A.recode_hets.vcf -O v R014835_A.recode.vcf.gz

bgzip R014813_A.recode.vcf
tabiz R014813_A.recode.vcf.gz
bcftools view -g het -o R014813_A.recode_hets.vcf -O v R014813_A.recode.vcf.gz

##Extract heterozygous variants from offspring for compund heterozygote based analysis
bcftools view -g het -o R014798_A.recode_hets.vcf -O v R014798_A.recode.vcf.gz
bcftools view -g het -o R014941_A.recode_hets.vcf -O v R014941_A.recode.vcf.gz

##Transfer to local server

#####################################
##Homozygous variant based analysis##
#####################################

##Get list of positions shared between files
bcftools isec -n~1111 -c all R014798_A.recode_homs.vcf.gz R014941_A.recode_homs.vcf.gz R014813_A.recode_hets.vcf.gz R014835_A.recode_hets.vcf.gz -o hom_query_matches.txt

##Merge files
bcftools merge R014798_A.recode_homs.vcf.gz R014941_A.recode_homs.vcf.gz R014813_A.recode_hets.vcf.gz R014835_A.recode_hets.vcf.gz -m all -O z -o 1981_hom_merged.vcf.gz

##Extract relevant variants from merged vcf for homozygote hypothesis
tabix 1981_hom_merged.vcf.gz
bcftools view -R hom_query_matches.txt -o hom_query.vcf -O v 1981_hom_merged.vcf.gz

##Execute Variant Effect Predictor
/home/jww39/ensembl-vep/vep -i hom_query.vcf -offline --assembly GRCh37 -o VEP_out_hom_query.vcf --pick --pick_order canonical --everything --vcf --plugin LoF --plugin CADD,/home/jww39/.vep/Plugins/whole_genome_SNVs.tsv.gz,/home/jww39/.vep/Plugins/InDels.tsv.gz --vcf_info_field ANN

##Filter with Variant Effect Predictor filter script
/home/jww39/ensembl-vep/filter_vep -i VEP_out_hom_query.vcf -o VEP_out_filtered_hom_query.vcf -filter "Consequence in /home/jww39/recessive/consequence.txt" --only_matched
/home/jww39/ensembl-vep/filter_vep -i VEP_out_filtered_hom_query.vcf -o VEP_out_double_filtered_hom_query.vcf -filter "EUR_AF < 0.05 or not EUR_AF" --only_matched 

################################################
##Compound heterozygote variant based analysis##
################################################

##Get list of positions shared between files (compound homozygote hypothesis)

##Heterozygous variants in parent 1
bcftools isec -n~1110 -c all R014798_A.recode_hets.vcf.gz R014941_A.recode_hets.vcf.gz R014813_A.recode_hets.vcf.gz R014835_A.recode_hets.vcf.gz -o comp_hom_query_matches_parent_1.txt

##Heterozygous variants in parent 2
bcftools isec -n~1101 -c all R014798_A.recode_hets.vcf.gz R014941_A.recode_hets.vcf.gz R014813_A.recode_hets.vcf.gz R014835_A.recode_hets.vcf.gz -o comp_hom_query_matches_parent_2.txt

##Combine index files to give heterozygous positions that are in both offspring and one parent only
> comp_hom_query_matches_both_parents.txt
cat comp_hom_query_matches_parent_1.txt >> comp_hom_query_matches_both_parents.txt
cat comp_hom_query_matches_parent_2.txt >> comp_hom_query_matches_both_parents.txt

##Merge files
bcftools merge R014798_A.recode_hets.vcf.gz R014941_A.recode_hets.vcf.gz R014813_A.recode_hets.vcf.gz R014835_A.recode_hets.vcf.gz -m all -O z -o 1981_comp_hom_merged.vcf.gz

##Extract relevant variants from merged vcf for compound homozygote hypothesis
tabix 1981_comp_hom_merged.vcf.gz
bcftools view -R comp_hom_query_matches_both_parents.txt -o comp_hom_query.vcf -O v 1981_comp_hom_merged.vcf.gz

##Execute Variant Effect Predictor
/home/jww39/ensembl-vep/vep -i comp_hom_query.vcf -offline --assembly GRCh37 -o VEP_out_comp_hom_query.vcf --pick --pick_order canonical --everything --vcf --plugin LoF --plugin CADD,/home/jww39/.vep/Plugins/whole_genome_SNVs.tsv.gz,/home/jww39/.vep/Plugins/InDels.tsv.gz --vcf_info_field ANN

##Filter with Variant Effect Predictor filter script
/home/jww39/ensembl-vep/filter_vep -i VEP_out_comp_hom_query.vcf -o VEP_out_filtered_comp_hom_query.vcf -filter "Consequence in /home/jww39/recessive/consequence.txt" --only_matched
/home/jww39/ensembl-vep/filter_vep -i VEP_out_filtered_comp_hom_query.vcf -o VEP_out_double_filtered_comp_hom_query.vcf -filter "AF < 0.05 or not AF" --only_matched 







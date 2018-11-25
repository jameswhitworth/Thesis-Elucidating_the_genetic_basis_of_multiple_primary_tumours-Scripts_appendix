
############################################################
##Calculation of coverage in regions specified in BED file##
############################################################

##Coverage in WGS data in regions specified in virtual panel BED file
for i in `cat bams.txt`; do
	samtools depth -b vp_exons_for_cov.bed /scratch/WGS10K/data/release/20170614-A/bam/MPMT/${i} > ${i}.cov
done

ls *.bam.cov > cov_files.txt
for i in `cat cov_files.txt`; do
  tr '\t' ',' < ${i} > ${i}.csv
done

touch WGS_vp_exons_cov.csv
for i in `cat cov_files.txt`; do
  cat ${i}.csv >> WGS_vp_exons_cov.csv
done

##Coverage in WGS data in regions specified in virtual panel BED file that correspond to genes seuqenced by Illumina TruSight Cancer panel
for i in `cat bams.txt`; do
	samtools depth -b vp_exons_in_TCP_for_cov.bed /scratch/WGS10K/data/release/20170614-A/bam/MPMT/${i} > ${i}.cov
done

ls *.bam.cov > cov_files.txt
for i in `cat cov_files.txt`; do
  tr '\t' ',' < ${i} > ${i}.csv
done

touch WGS_vp_exons_in_TCP_cov.csv
for i in `cat cov_files.txt`; do
  cat ${i}.csv >> WGS_vp_exons_in_TCP_cov.csv
done

##Export coverage data from WGS to local server

##Coverage in Illumina Trusight Cancer panel data in regions specified in virtual panel BED file that correspond to genes seuqenced by Illumin TruSight Cancer panel
for i in `cat bams.txt`; do
	samtools depth -b vp_exons_in_TCP_for_cov.bed /home/jww39/cancer_panel_data/CancerPanel_hg19/${i} > ${i}.cov
done

ls *.cov > cov_files.txt
for i in `cat cov_files.txt`; do
  tr '\t' ',' < ${i} > ${i}.csv
done

touch TCP_vp_exons_in_TCP_cov.csv
for i in `cat cov_files.txt`; do
  cat ${i}.csv >> TCP_vp_exons_in_TCP_cov.csv
done


#######################################
##Calculation of coverage stattistics##
#######################################

##Coverage in WGS data in regions specified in virtual panel BED file

#Read in samtools depth output combined across samples
WGS_vp_exons_cov <- read.csv("WGS_vp_exons_cov.csv", header = F)
colnames(WGS_vp_exons_cov) <- c("chr", "base", "depth")

#Calculate mean coverage
mean_WGS_vp_exons_cov <- mean(WGS_vp_exons_cov$depth)
sd_WGS_vp_exons_cov <- sd(WGS_vp_exons_cov$depth)

#percent bases at 10X or above
WGS_vp_exons_cov_percent_at_10 <- (100/nrow(WGS_vp_exons_cov)) * (length(WGS_vp_exons_cov$depth[which(WGS_vp_exons_cov$depth > 9)]))


##Coverage in WGS data in regions specified in virtual panel BED file that correspond to genes seuqenced by Illumina TruSight Cancer panel

#Read in samtools depth output combined across samples
WGS_vp_exons_in_TCP_cov <- read.csv("WGS_vp_exons_in_TCP_cov.csv", header = F)
colnames(WGS_vp_exons_in_TCP_cov) <- c("chr", "base", "depth")

#Calculate mean coverage
mean_WGS_vp_exons_in_TCP_cov <- mean(WGS_vp_exons_in_TCP_cov$depth)
sd_WGS_vp_exons_in_TCP_cov <- sd(WGS_vp_exons_in_TCP_cov$depth)

#percent bases at 10X or above
WGS_vp_exons_in_TCP_cov_percent_at_10 <- (100/nrow(WGS_vp_exons_in_TCP_cov)) * (length(WGS_vp_exons_in_TCP_cov$depth[which(WGS_vp_exons_in_TCP_cov$depth > 9)]))


##Coverage in Illumina Trusight Cancer panel data in regions specified in virtual panel BED file that correspond to genes seuqenced by Illumin TruSight Cancer panel

#Read in samtools depth output combined across samples
TCP_vp_exons_in_TCP_cov <- read.csv("TCP_vp_exons_in_TCP_cov.csv", header = F)
colnames(TCP_vp_exons_in_TCP_cov) <- c("chr", "base", "depth")

#Calculate mean coverage
mean_TCP_vp_exons_in_TCP_cov <- mean(TCP_vp_exons_in_TCP_cov$depth)
sd_TCP_vp_exons_in_TCP_cov <- sd(TCP_vp_exons_in_TCP_cov$depth)

#percent bases at 10X or above
TCP_vp_exons_in_TCP_cov_percent_at_10 <- (100/nrow(TCP_vp_exons_in_TCP_cov)) * (length(TCP_vp_exons_in_TCP_cov$depth[which(TCP_vp_exons_in_TCP_cov$depth > 9)]))


##Output results
cov_stats_output <- t(data.frame(mean_WGS_vp_exons_cov,
                                 sd_WGS_vp_exons_cov,
                                 WGS_vp_exons_cov_percent_at_10,
                                 
                                 mean_WGS_vp_exons_in_TCP_cov,
                                 sd_WGS_vp_exons_in_TCP_cov,
                                 WGS_vp_exons_in_TCP_cov_percent_at_10,
                                 
                                 mean_TCP_vp_exons_in_TCP_cov,
                                 sd_TCP_vp_exons_in_TCP_cov,
                                 TCP_vp_exons_in_TCP_cov_percent_at_10))
write.csv(as.data.frame(cov_stats_output), "cov_stats_output.csv")

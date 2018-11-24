##############################################################################################################################################################################
##Estimate telomere length of selected BAM files provided by NIHR BioResource Rare Diseases project and stored on University of Cambridge high performance computing cluster##
##############################################################################################################################################################################

##Produce telbam files from BAM files
for i in `cat telomere_samples_MPT_BRIDGE_euro.txt`;do
  /rds/user/jww39/hpc-work/telomeres/telomerecat-3.2/telomerecat/bin/telomerecat bam2telbam ${i}
done

##Transfer length estimate files to local server
ls *_telbam.bam > telbams_for_length.txt

##Produce length estimates from telbam files. 10 estimates for each BAM file into seaparte folder
source /home/jww39/telomeres/virtualenvironment/env/bin/activate
for i in `cat telbams_for_length.txt`;do
  telomerecat telbam2length /home/jww39/telomeres/lengths_1/${i}
done

source /home/jww39/telomeres/virtualenvironment/env/bin/activate
for i in `cat telbams_for_length.txt`;do
  telomerecat telbam2length /home/jww39/telomeres/lengths_2/${i}
done

source /home/jww39/telomeres/virtualenvironment/env/bin/activate
for i in `cat telbams_for_length.txt`;do
  telomerecat telbam2length /home/jww39/telomeres/lengths_3/${i}
done

source /home/jww39/telomeres/virtualenvironment/env/bin/activate
for i in `cat telbams_for_length.txt`;do
  telomerecat telbam2length /home/jww39/telomeres/lengths_4/${i}
done

source /home/jww39/telomeres/virtualenvironment/env/bin/activate
for i in `cat telbams_for_length.txt`;do
  telomerecat telbam2length /home/jww39/telomeres/lengths_5/${i}
done

source /home/jww39/telomeres/virtualenvironment/env/bin/activate
for i in `cat telbams_for_length.txt`;do
  telomerecat telbam2length /home/jww39/telomeres/lengths_6/${i}
done

source /home/jww39/telomeres/virtualenvironment/env/bin/activate
for i in `cat telbams_for_length.txt`;do
  telomerecat telbam2length /home/jww39/telomeres/lengths_7/${i}
done

source /home/jww39/telomeres/virtualenvironment/env/bin/activate
for i in `cat telbams_for_length.txt`;do
  telomerecat telbam2length /home/jww39/telomeres/lengths_8/${i}
done

source /home/jww39/telomeres/virtualenvironment/env/bin/activate
for i in `cat telbams_for_length.txt`;do
  telomerecat telbam2length /home/jww39/telomeres/lengths_9/${i}
done

source /home/jww39/telomeres/virtualenvironment/env/bin/activate
for i in `cat telbams_for_length.txt`;do
  telomerecat telbam2length /home/jww39/telomeres/lengths_10/${i}
done

##########################################################
##R script to analyse telmomere length (run in R studio)##
##########################################################

##Combine telomere length estimates for all 10 runs
lengths_1 <- read.csv("/home/jww39/telomeres/lengths_1/collated_tel_lengths.csv")
lengths_2 <- read.csv("/home/jww39/telomeres/lengths_2/collated_tel_lengths.csv")
lengths_3 <- read.csv("/home/jww39/telomeres/lengths_3/collated_tel_lengths.csv")
lengths_4 <- read.csv("/home/jww39/telomeres/lengths_4/collated_tel_lengths.csv")
lengths_5 <- read.csv("/home/jww39/telomeres/lengths_5/collated_tel_lengths.csv")
lengths_6 <- read.csv("/home/jww39/telomeres/lengths_6/collated_tel_lengths.csv")
lengths_7 <- read.csv("/home/jww39/telomeres/lengths_7/collated_tel_lengths.csv")
lengths_8 <- read.csv("/home/jww39/telomeres/lengths_8/collated_tel_lengths.csv")
lengths_9 <- read.csv("/home/jww39/telomeres/lengths_9/collated_tel_lengths.csv")
lengths_10 <- read.csv("/home/jww39/telomeres/lengths_10/collated_tel_lengths.csv")

lengths_all <- rbind(lengths_1, lengths_2, lengths_3, lengths_4, lengths_5, lengths_6, lengths_7, lengths_8, lengths_9, lengths_10)

##Produce data frame with mean telomere length for all samples
samples <- data.frame(unique(lengths_all$Sample))
samples_with_mean_length <- apply(samples, 1, function(x) paste(x,mean(lengths_all[which(lengths_all$Sample == x),]$Length),sep = ":"))
samples_with_mean_length_matrix <- t(data.frame(strsplit(samples_with_mean_length, ":")))
samples_with_mean_length_df <- data.frame(cbind(samples_with_mean_length_matrix[,1],samples_with_mean_length_matrix[,2]))
rownames(samples_with_mean_length_df) <- seq(1:nrow(samples_with_mean_length_df))
colnames(samples_with_mean_length_df) <- c("sample", "mean_length")
samples_with_mean_length_df$sample <- gsub("_A.bam", "", samples_with_mean_length_df$sample)

##Read in ages file
ages <- read.csv("All_WGS_ages_WGSID_Feb2018v2_JW_edits.csv")

##Join lengths to ages
library(sqldf)
length_age <- sqldf("select * from samples_with_mean_length_df f1 left join ages f2 on (f1.sample == f2.WGSID)")

##Read in projects file
include_samples <- data.frame(read.delim("include_samples_20170614-A.txt", header = F))
colnames(include_samples) <- c("inc_ilumina_ID", "inc_WGSID", "inc_sample", "project", "inc_issues")

##Join length_age to project
length_age_project <- sqldf("select * from length_age f1 left join include_samples f2 on (f1.sample == f2.inc_sample)")
length_age_project <- data.frame(length_age_project$sample, length_age_project$mean_length, length_age_project$Age_at_sampling, length_age_project$project)
colnames(length_age_project) <- c("sample", "mean_length", "age", "project")

##Exclude samples where no age at sampling recorded
length_age_project_no_missing <- length_age_project[which(!(is.na(length_age_project$age))),]
length_age_project_no_missing$mean_length <- as.numeric(as.character(length_age_project_no_missing$mean_length))

##Replacing empty project columns with MPMT tag (as all were MPMT). Beware of addiing more samples if missing values in inc samples file
length_age_project_no_missing[which(is.na(length_age_project_no_missing$project)),][4] <- "MPMT"

##Tests of normality in length data
hist(length_age_project_no_missing$mean_length)
qqnorm(length_age_project_no_missing$mean_length)
qqline(length_age_project_no_missing$mean_length)
shapiro.test(length_age_project_no_missing$mean_length)

##Produce linear regression model for length (including MPMT samples) against age.

##Output residuals for further analysis
lm.1 <- lm(mean_length~age, data=length_age_project_no_missing)
summary(lm.1)


##Plot of residuals by age
plot_colour <- vector(mode = "character", length =nrow(length_age_project_no_missing))
plot_pch <- vector(mode = "numeric", length =nrow(length_age_project_no_missing))

plot_colour[length_age_project_no_missing$MPMT_col == "MPMT"] <- "red"
plot_pch[length_age_project_no_missing$MPMT_col == "MPMT"] <- 20

plot_colour[length_age_project_no_missing$MPMT_col != "MPMT"] <- "green4"
plot_pch[length_age_project_no_missing$MPMT_col != "MPMT"] <- 1
  
png("linear_model.png",
height = 450,
width = 650)

plot(mean_length~age, data=length_age_project_no_missing, 
     pch = plot_pch, 
     lwd = 0.5,
     cex = 0.5,
     ylab = "Mean length",
     xlab = "Age of participant",
     cex.axis = 1.5,
     cex.lab = 1.5,
     col = plot_colour)

abline(lm.1)

dev.off()

length_age_project_no_missing <- data.frame(length_age_project_no_missing, resid(lm.1))
colnames(length_age_project_no_missing)[5] <- "residual"

##Make MPMT and non-MPMT column
MPMT_col <- length_age_project_no_missing$project == "MPMT"
MPMT_col <-  gsub("TRUE", "MPMT",MPMT_col)
MPMT_col <-  gsub("FALSE", "non_MPMT",MPMT_col)

length_age_project_no_missing <- data.frame(length_age_project_no_missing, MPMT_col )

##Get mean ages
c(
mean(length_age_project_no_missing$age[which(length_age_project_no_missing$project == "BPD")]),
mean(length_age_project_no_missing$age[which(length_age_project_no_missing$project == "PID")]),
mean(length_age_project_no_missing$age[which(length_age_project_no_missing$project == "SPEED")]),
mean(length_age_project_no_missing$age[which(length_age_project_no_missing$project == "PMG")]),
mean(length_age_project_no_missing$age[which(length_age_project_no_missing$project == "MPMT")]),
mean(length_age_project_no_missing$age[which(length_age_project_no_missing$project == "ICP")]),
mean(length_age_project_no_missing$age[which(length_age_project_no_missing$project == "HCM")]),
mean(length_age_project_no_missing$age[which(length_age_project_no_missing$project == "CSVD")]),
mean(length_age_project_no_missing$age[which(length_age_project_no_missing$project == "NPD")]),
mean(length_age_project_no_missing$age[which(length_age_project_no_missing$project == "SRNS")])
)

##Get number of samples per project
c(
  nrow(length_age_project_no_missing[which(length_age_project_no_missing$project == "BPD"),]),
  nrow(length_age_project_no_missing[which(length_age_project_no_missing$project == "PID"),]),
  nrow(length_age_project_no_missing[which(length_age_project_no_missing$project == "SPEED"),]),
  nrow(length_age_project_no_missing[which(length_age_project_no_missing$project == "PMG"),]),
  nrow(length_age_project_no_missing[which(length_age_project_no_missing$project == "MPMT"),]),
  nrow(length_age_project_no_missing[which(length_age_project_no_missing$project == "ICP"),]),
  nrow(length_age_project_no_missing[which(length_age_project_no_missing$project == "HCM"),]),
  nrow(length_age_project_no_missing[which(length_age_project_no_missing$project == "CSVD"),]),
  nrow(length_age_project_no_missing[which(length_age_project_no_missing$project == "NPD"),]),
  nrow(length_age_project_no_missing[which(length_age_project_no_missing$project == "SRNS"),])
)

##Get list of projects
projects <- unique(length_age_project_no_missing$project)

##Make vectors with lengths by project
MPMT_residuals <- length_age_project_no_missing$residual[which(length_age_project_no_missing$project == "MPMT")]
BPD_residuals <- length_age_project_no_missing$residual[which(length_age_project_no_missing$project == "BPD")]
PID_residuals <- length_age_project_no_missing$residual[which(length_age_project_no_missing$project == "PID")]
SPEED_residuals <- length_age_project_no_missing$residual[which(length_age_project_no_missing$project == "SPEED")]
PMG_residuals <- length_age_project_no_missing$residual[which(length_age_project_no_missing$project == "PMG")]
ICP_residuals <- length_age_project_no_missing$residual[which(length_age_project_no_missing$project == "ICP")]
HCM_residuals <- length_age_project_no_missing$residual[which(length_age_project_no_missing$project == "HCM")]
CSVD_residuals <- length_age_project_no_missing$residual[which(length_age_project_no_missing$project == "CSVD")]
NPD_residuals <- length_age_project_no_missing$residual[which(length_age_project_no_missing$project == "NPD")]
SRNS_residuals <- length_age_project_no_missing$residual[which(length_age_project_no_missing$project == "SRNS")]
non_MPMT_residuals <- length_age_project_no_missing$residual[which(length_age_project_no_missing$project != "MPMT")]

##Boxplot of lengths by project
boxplot(
  MPMT_residuals,
  non_MPMT_residuals,
  BPD_residuals,
  PID_residuals,
  SPEED_residuals,
  PMG_residuals,
  ICP_residuals,
  HCM_residuals,
  CSVD_residuals,
  NPD_residuals,
  SRNS_residuals,
  
  names = c("MPMT","non_MPMT","BPD" ,"PID","SPEED","PMG","ICP","HCM","CSVD","NPD","SRNS"),
  range = 0
 
)

##Do any of the projects have significantly different residuals (ANOVA)?

##Could do one for the control samples and exclude a project that does

##Boxplot of residuals by project
png("all_projects_resid.png",
    height = 600,
    width = 1000)

boxplot(residual~project, data=length_age_project_no_missing, 
        xlab = "NIHR BioResource sub-project",
        ylab = "Residual from linear model",
        cex.lab = 1.5,
        cex.axis = 1.5)

dev.off()

##ANOVA to see if all residuals by project from same parent distribution
summary(aov(residual~project, data=length_age_project_no_missing))
##Answer is no. Probablity of coming from same parent distribution below p value threshold.

##ANOVA to see if all residuals by project from same parent distribution (non_MPMT only)
summary(aov(residual~project, data=length_age_project_no_missing[which(length_age_project_no_missing$project != "MPMT"),]))
##Answer is no. Probablity of coming from same parent distribution below p value threshold.

##ANOVA to see if all residuals by project from same parent distribution (non_MPMT and non_SPEED only)
summary(aov(residual~project, data=length_age_project_no_missing[which(length_age_project_no_missing$project != "SPEED" & length_age_project_no_missing$project != "MPMT"),]))
##Answer is no. Probablity of coming from same parent distribution below p value threshold.

##Test variance to decide between students (equal variance) and welchs (unequal variance) t test.
library(car)
bartlett.test(residual~MPMT_col, data=length_age_project_no_missing)
##p value low. Therefore reject null hypothesis that the variances are equal. Variances appear unequeal

##So use Welchs t-test (indicated by var.equal=F argument)

##Boxplot of MPMT residuals vs non-MPMT residuals
png("MPMT_resid_vs_non_MPMT_resid.png",
    height = 650,
    width = 400)


plot(residual~MPMT_col, data=length_age_project_no_missing,
     xlab = "NIHR BioResource sub-project",
     ylab = "Residual from linear model",
     names = c("MPMT", "non MPMT"),
     cex.lab = 1.5,
     cex.axis = 1.5)


dev.off()

t.test(residual~MPMT_col, data=length_age_project_no_missing, alternative = "two.sided", var.equal=F)

##Conclude that MPMT residuals are lower than non-MPMT. Could be a number of reasons for this such as chemotherapy.

##Within MPMT. Are the residuals different in cases with pathogenic variants?
length_age_project_MPMT_only <- length_age_project_no_missing[which(length_age_project_no_missing$project == "MPMT"),]

mut_pos_cases <- read.delim("mut_pos_cases_30_5_18.txt", header = F)
mut_col <- length_age_project_MPMT_only$sample %in% mut_pos_cases$V1

length_age_project_MPMT_only <- data.frame(length_age_project_MPMT_only, mut_col)

library(car)
bartlett.test(residual~mut_col, data=length_age_project_MPMT_only)
##p value high. Therefore don't reject null hypothesis that the varinaces are equal. Variances appear equeal

##So use Students t-test (indicated by var.equal=T argument)
t.test(residual~mut_col, data=length_age_project_MPMT_only, alternative = "two.sided", var.equal=T)
##Can't reject null hypothesis that no difference between the groups

####################################################################################################################
##Selection of MPMT cases with shortest and longest telomeres for assessment of variants in telomere related genes##
####################################################################################################################

##Take bottom quartile of MPMT residuals (who don't have a pathogenic variant)
non_mut_residuals <- length_age_project_MPMT_only$residual[which(length_age_project_MPMT_only$mut_col == "FALSE")]
low_resid_MPMT_samples <- length_age_project_MPMT_only$sample[which(length_age_project_MPMT_only$residual < quantile(non_mut_residuals, 0.25))]
write (as.character(low_resid_MPMT_samples), "low_resid_MPMT_samples.txt")

##Take top quartile of MPMT residuals (who don't have a pathogenic variant)
high_resid_MPMT_samples <- length_age_project_MPMT_only$sample[which(length_age_project_MPMT_only$residual > quantile(non_mut_residuals, 0.75))]
write(as.character(high_resid_MPMT_samples), "high_resid_MPMT_samples.txt")













##Read in table with record of personal and family history of tumours
MTS <- read.csv("MTS_sheet_22_5_18.csv", header = T)


################################################
##Analysis of individual independent variables##
################################################

##Prepare tables for individual variable analysis
MTS$age_ind_variable <- as.numeric(MTS$age_ind_variable)
MTS$hereditability_figure <- as.numeric(MTS$hereditability_figure)
MTS$incidence_figure <- as.numeric(MTS$incidence_figure)
MTS_with_excl <- MTS[which(MTS$exclude == "NO" & 
                             MTS$heritability_exclude == "NO" & 
                             MTS$incidence_exclude == "NO" &
                             MTS$vp_2018_proband == "YES"),]
all_IDs <- unique(MTS$WGS_ID)

##Exclude samples where an absent hertitability or incidence figure leads to the sample being excluded

##Find samples IDs where there is more than on line following exclusions and where tumour in participant is yes
two_tumour <- lapply(all_IDs, function(x) nrow(MTS_with_excl[which(MTS_with_excl$WGS_ID == x & MTS_with_excl$tumour_in_participant == "YES"),]))
two_tumour_samples <- all_IDs[which(two_tumour > 1)]

##Find sample IDs where multiple paramter = "YES"
multiple_samples <- MTS_with_excl$WGS_ID[which(MTS_with_excl$multiple == "YES" & MTS_with_excl$tumour_in_participant == "YES")]

##Combine the factors to give a list of samples that avoid exclusion
include_samples <- unique(c(as.character(two_tumour_samples), as.character(multiple_samples)))

##Produce data frame with only samples with multiple tumours (excludes those who qualify through relatives)
MTS_all <- MTS_with_excl[which(MTS_with_excl$WGS_ID %in% include_samples),]


##Subset cases from virtual panel analysis

##With family history
MTS_vp_fam <- MTS_all[which(MTS_all$cohort == "AJHG" & MTS_all$fhx == "YES"),]
vp_IDs_fam <- unique(MTS_vp_fam$WGS_ID)

##Without family history
MTS_vp_indv <- MTS_all[which(MTS_all$cohort == "AJHG" & MTS_all$tumour_in_participant == "YES"),]
vp_IDs_indv <- unique(MTS_vp_indv$WGS_ID)

##Create dataframe with variables to perform individual variable logistic regressions to assess weighting of score

##With family history
mean_age <- lapply(vp_IDs_fam, function(x) mean(MTS_vp_fam$age_ind_variable[which(MTS_vp_fam$WGS_ID == x)]))
mean_hered <- lapply(vp_IDs_fam, function(x) mean(MTS_vp_fam$hereditability_figure[which(MTS_vp_fam$WGS_ID == x)]))
mean_incid <- lapply(vp_IDs_fam, function(x) mean(MTS_vp_fam$incidence_figure[which(MTS_vp_fam$WGS_ID == x)]))

mono_var_df_vp_fam <- as.data.frame(cbind(as.character(vp_IDs_fam),mean_age, mean_hered, mean_incid))
colnames(mono_var_df_vp_fam)[1] <- "WGS_ID"

##Without family history
mean_age <- lapply(vp_IDs_indv, function(x) mean(MTS_vp_indv$age_ind_variable[which(MTS_vp_indv$WGS_ID == x)]))
mean_hered <- lapply(vp_IDs_indv, function(x) mean(MTS_vp_indv$hereditability_figure[which(MTS_vp_indv$WGS_ID == x)]))
mean_incid <- lapply(vp_IDs_indv, function(x) mean(MTS_vp_indv$incidence_figure[which(MTS_vp_indv$WGS_ID == x)]))

mono_var_df_vp_indv <- as.data.frame(cbind(as.character(vp_IDs_indv),mean_age, mean_hered, mean_incid))
colnames(mono_var_df_vp_indv)[1] <- "WGS_ID"

##Assign mutation status to data frames

##With family history
mut_status_fam <- unlist(lapply(mono_var_df_vp_fam$WGS_ID, function(x) MTS$mutation[which.max(MTS$WGS_ID == x)]))
mono_var_df_vp_fam <- data.frame(mono_var_df_vp_fam,mut_status_fam)

##Without family history
mut_status_indv <- unlist(lapply(mono_var_df_vp_indv$WGS_ID, function(x) MTS$mutation[which.max(MTS$WGS_ID == x)]))
mono_var_df_vp_indv <- data.frame(mono_var_df_vp_indv, mut_status_indv)

length(which(mono_var_df_vp_indv$mut_status_indv == "YES"))

##Perform logistic regression to assess individual variables
model_outputs_text <- vector()

##With family history - Age
mono_var_fam_age_mod <- glm(formula = (as.numeric(mut_status_fam) -1)~as.numeric(mean_age), data = mono_var_df_vp_fam, family = binomial(link = "logit"))
library(pROC)
mono_var_fam_age_roc <- roc(as.numeric(mut_status_fam) -1~as.numeric(mean_age), data = mono_var_df_vp_fam)
model_outputs_text <- append(model_outputs_text, capture.output(summary(mono_var_fam_age_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(mono_var_fam_age_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(mono_var_fam_age_roc))

##With family history - Heritability
mono_var_fam_hered_mod <- glm(formula = (as.numeric(mut_status_fam) -1)~as.numeric(mean_hered), data = mono_var_df_vp_fam, family = binomial(link = "logit"))
library(pROC)
mono_var_fam_hered_roc <- roc(as.numeric(mut_status_fam) -1~as.numeric(mean_age), data = mono_var_df_vp_fam)
model_outputs_text <- append(model_outputs_text, capture.output(summary(mono_var_fam_hered_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(mono_var_fam_hered_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(mono_var_fam_hered_roc))

##With family history - Incidence
mono_var_fam_incid_mod <- glm(formula = (as.numeric(mut_status_fam) -1)~as.numeric(mean_incid), data = mono_var_df_vp_fam, family = binomial(link = "logit"))
library(pROC)
mono_var_fam_incid_roc <- roc(as.numeric(mut_status_fam) -1~as.numeric(mean_age), data = mono_var_df_vp_fam)
model_outputs_text <- append(model_outputs_text, capture.output(summary(mono_var_fam_incid_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(mono_var_fam_incid_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(mono_var_fam_incid_roc))

##Without family history - Age
mono_var_indv_age_mod <- glm(formula = (as.numeric(mut_status_indv) -1)~as.numeric(mean_age), data = mono_var_df_vp_indv, family = binomial(link = "logit"))
library(pROC)
mono_var_indv_age_roc <- roc(as.numeric(mut_status_indv) -1~as.numeric(mean_age), data = mono_var_df_vp_indv)
model_outputs_text <- append(model_outputs_text, capture.output(summary(mono_var_indv_age_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(mono_var_indv_age_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(mono_var_indv_age_roc))

##Without family history - Heritability
mono_var_indv_hered_mod <- glm(formula = (as.numeric(mut_status_indv) -1)~as.numeric(mean_hered), data = mono_var_df_vp_indv, family = binomial(link = "logit"))
library(pROC)
mono_var_indv_hered_roc <- roc(as.numeric(mut_status_indv) -1~as.numeric(mean_hered), data = mono_var_df_vp_indv)
model_outputs_text <- append(model_outputs_text, capture.output(summary(mono_var_indv_hered_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(mono_var_indv_hered_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(mono_var_indv_hered_roc))

##Without family history - Incidence
mono_var_indv_incid_mod <- glm(formula = (as.numeric(mut_status_indv) -1)~as.numeric(mean_incid), data = mono_var_df_vp_indv, family = binomial(link = "logit"))
library(pROC)
mono_var_indv_incid_roc <- roc(as.numeric(mut_status_indv) -1~as.numeric(mean_incid), data = mono_var_df_vp_indv)
model_outputs_text <- append(model_outputs_text, capture.output(summary(mono_var_indv_incid_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(mono_var_indv_incid_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(mono_var_indv_incid_roc))


#################################################################
##Analysis of proposed scoring systems as independent variables##
#################################################################
 
##Reset tables so that cases excluded due to lack of heritability estimate or incidence for their tumours
MTS_all <- MTS[which(MTS$exclude == "NO" & MTS$vp_2018_proband == "YES"),]

##Assign scoring scores to tumours according to different proposed scoring systems

##Age 
option_1_age <- rep("empty", nrow(MTS_all))
option_1_age[which(MTS_all$age_band == "<30")] <- 4
option_1_age[which(MTS_all$age_band == "30-44")] <- 3
option_1_age[which(MTS_all$age_band == "45-59")] <- 2
option_1_age[which(MTS_all$age_band == ">59")] <- 1

option_2_age <- rep("empty", nrow(MTS_all))
option_2_age[which(MTS_all$age_band == "<30")] <- 8
option_2_age[which(MTS_all$age_band == "30-44")] <- 4
option_2_age[which(MTS_all$age_band == "45-59")] <- 2
option_2_age[which(MTS_all$age_band == ">59")] <- 1

option_3_age <- rep("empty", nrow(MTS_all))
option_3_age[which(MTS_all$age_band == "<30")] <- 27
option_3_age[which(MTS_all$age_band == "30-44")] <- 9
option_3_age[which(MTS_all$age_band == "45-59")] <- 3
option_3_age[which(MTS_all$age_band == ">59")] <- 1

option_4_age <- rep("empty", nrow(MTS_all))
option_4_age[which(MTS_all$age_band == "<30")] <- 2
option_4_age[which(MTS_all$age_band == "30-44")] <- 4
option_4_age[which(MTS_all$age_band == "45-59")] <- 8
option_4_age[which(MTS_all$age_band == ">59")] <- 16

option_5_age <- rep("empty", nrow(MTS_all))
option_5_age[which(MTS_all$age_band == "<30")] <- 2
option_5_age[which(MTS_all$age_band == "30-44")] <- 6
option_5_age[which(MTS_all$age_band == "45-59")] <- 18
option_5_age[which(MTS_all$age_band == ">59")] <- 54

option_6_age <- rep("empty", nrow(MTS_all))
option_6_age[which(MTS_all$age_band == "<30")] <- 1
option_6_age[which(MTS_all$age_band == "30-44")] <- 10
option_6_age[which(MTS_all$age_band == "45-59")] <- 20
option_6_age[which(MTS_all$age_band == ">59")] <- 30

##Heritability 
option_1_heritability <- rep("empty", nrow(MTS_all))
option_1_heritability[which(MTS_all$heritability_band == "0-25")] <- 1
option_1_heritability[which(MTS_all$heritability_band == "26-50")] <- 2
option_1_heritability[which(MTS_all$heritability_band == "51-75")] <- 3
option_1_heritability[which(MTS_all$heritability_band == "76-100")] <- 4

option_2_heritability <- rep("empty", nrow(MTS_all))
option_2_heritability[which(MTS_all$heritability_band == "0-25")] <- 1
option_2_heritability[which(MTS_all$heritability_band == "26-50")] <- 2
option_2_heritability[which(MTS_all$heritability_band == "51-75")] <- 4
option_2_heritability[which(MTS_all$heritability_band == "76-100")] <- 8

option_3_heritability <- rep("empty", nrow(MTS_all))
option_3_heritability[which(MTS_all$heritability_band == "0-25")] <- 1
option_3_heritability[which(MTS_all$heritability_band == "26-50")] <- 3
option_3_heritability[which(MTS_all$heritability_band == "51-75")] <- 9
option_3_heritability[which(MTS_all$heritability_band == "76-100")] <- 27

option_4_heritability <- rep("empty", nrow(MTS_all))
option_4_heritability[which(MTS_all$heritability_band == "0-25")] <- 2
option_4_heritability[which(MTS_all$heritability_band == "26-50")] <- 4
option_4_heritability[which(MTS_all$heritability_band == "51-75")] <- 8
option_4_heritability[which(MTS_all$heritability_band == "76-100")] <- 16

option_5_heritability <- rep("empty", nrow(MTS_all))
option_5_heritability[which(MTS_all$heritability_band == "0-25")] <- 2
option_5_heritability[which(MTS_all$heritability_band == "26-50")] <- 6
option_5_heritability[which(MTS_all$heritability_band == "51-75")] <- 18
option_5_heritability[which(MTS_all$heritability_band == "76-100")] <- 54

option_6_heritability <- rep("empty", nrow(MTS_all))
option_6_heritability[which(MTS_all$heritability_band == "0-25")] <- 1
option_6_heritability[which(MTS_all$heritability_band == "26-50")] <- 10
option_6_heritability[which(MTS_all$heritability_band == "51-75")] <- 20
option_6_heritability[which(MTS_all$heritability_band == "76-100")] <- 30

##Incidence
option_1_incidence <- rep("empty", nrow(MTS_all))
option_1_incidence[which(MTS_all$incidence_band == "0-6")] <- 4
option_1_incidence[which(MTS_all$incidence_band == "6.1-28")] <- 3
option_1_incidence[which(MTS_all$incidence_band == "28.1>50")] <- 2
option_1_incidence[which(MTS_all$incidence_band == ">50")] <- 1

option_2_incidence <- rep("empty", nrow(MTS_all))
option_2_incidence[which(MTS_all$incidence_band == "0-6")] <- 8
option_2_incidence[which(MTS_all$incidence_band == "6.1-28")] <- 4
option_2_incidence[which(MTS_all$incidence_band == "28.1>50")] <- 2
option_2_incidence[which(MTS_all$incidence_band == ">50")] <- 1

option_3_incidence <- rep("empty", nrow(MTS_all))
option_3_incidence[which(MTS_all$incidence_band == "0-6")] <- 27
option_3_incidence[which(MTS_all$incidence_band == "6.1-28")] <- 9
option_3_incidence[which(MTS_all$incidence_band == "28.1>50")] <- 3
option_3_incidence[which(MTS_all$incidence_band == ">50")] <- 1

option_4_incidence <- rep("empty", nrow(MTS_all))
option_4_incidence[which(MTS_all$incidence_band == "0-6")] <- 16
option_4_incidence[which(MTS_all$incidence_band == "6.1-28")] <- 8
option_4_incidence[which(MTS_all$incidence_band == "28.1>50")] <- 4
option_4_incidence[which(MTS_all$incidence_band == ">50")] <- 2

option_5_incidence <- rep("empty", nrow(MTS_all))
option_5_incidence[which(MTS_all$incidence_band == "0-6")] <- 54
option_5_incidence[which(MTS_all$incidence_band == "6.1-28")] <- 18
option_5_incidence[which(MTS_all$incidence_band == "28.1>50")] <- 6
option_5_incidence[which(MTS_all$incidence_band == ">50")] <- 2

option_6_incidence <- rep("empty", nrow(MTS_all))
option_6_incidence[which(MTS_all$incidence_band == "0-6")] <- 30
option_6_incidence[which(MTS_all$incidence_band == "6.1-28")] <- 20
option_6_incidence[which(MTS_all$incidence_band == "28.1>50")] <- 10
option_6_incidence[which(MTS_all$incidence_band == ">50")] <- 1



##Assign scoring scores to tumours according to original multiple tumour scoring system

original <- rep("empty", nrow(MTS_all))
original[which(MTS_all$common_tumour == "YES" & MTS_all$age_band_for_original == "<30")] <- 5
original[which(MTS_all$common_tumour == "YES" & MTS_all$age_band_for_original == "30-39")] <- 4
original[which(MTS_all$common_tumour == "YES" & MTS_all$age_band_for_original == "40-49")] <- 3
original[which(MTS_all$common_tumour == "YES" & MTS_all$age_band_for_original == "50-59")] <- 2
original[which(MTS_all$common_tumour == "YES" & MTS_all$age_band_for_original == ">59")] <- 1
original[which(MTS_all$common_tumour == "NO" & MTS_all$age_band_for_original == "<50")] <- 5
original[which(MTS_all$common_tumour == "NO" & MTS_all$age_band_for_original == "50-59")] <- 3
original[which(MTS_all$common_tumour == "NO" & MTS_all$age_band_for_original == ">59")] <- 1

##Add assigned scores to main table
MTS_all <- cbind(MTS_all,
      #
      as.numeric(option_1_age),
      as.numeric(option_2_age),
      as.numeric(option_3_age),
      as.numeric(option_4_age),
      as.numeric(option_5_age),
      as.numeric(option_6_age),
      #
      as.numeric(option_1_heritability),
      as.numeric(option_2_heritability),
      as.numeric(option_3_heritability),
      as.numeric(option_4_heritability),
      as.numeric(option_5_heritability),
      as.numeric(option_6_heritability),
      #
      as.numeric(option_1_incidence),
      as.numeric(option_2_incidence),
      as.numeric(option_3_incidence),
      as.numeric(option_4_incidence),
      as.numeric(option_5_incidence),
      as.numeric(option_6_incidence),
      
      as.numeric(original)
      
      )

colnames(MTS_all)[match("age_band_for_original",names(MTS_all)) +1:(ncol(MTS_all)-match("age_band_for_original",names(MTS_all)))] <-
  c("option_1_age",
    "option_2_age",
    "option_3_age",
    "option_4_age",
    "option_5_age",
    "option_6_age",
    
    "option_1_heritability",
    "option_2_heritability",
    "option_3_heritability",
    "option_4_heritability",
    "option_5_heritability",
    "option_6_heritability",
    
    "option_1_incidence",
    "option_2_incidence",
    "option_3_incidence",
    "option_4_incidence",
    "option_5_incidence",
    "option_6_incidence",
    
    "original")

##############################################################################
##Logistic regression based on scoring systems and assessment of performance##
##############################################################################

##Divide into training and test sets by randomly selectin individuals with and without pathognic variants separately

##With family history.
MTS_vp_fam <- MTS_all[which(MTS_all$cohort == "AJHG" & MTS_all$fhx == "YES"),]
vp_IDs_fam <- unique(MTS_vp_fam$WGS_ID)

train_samples_fam_mut <- sample(unique(MTS_vp_fam$WGS_ID[which(MTS_vp_fam$mutation == "YES")]), length(unique(MTS_vp_fam$WGS_ID[which(MTS_vp_fam$mutation == "YES")]))/2, replace=FALSE)
train_samples_fam_no_mut <- sample(unique(MTS_vp_fam$WGS_ID[which(MTS_vp_fam$mutation == "NO")]), length(unique(MTS_vp_fam$WGS_ID[which(MTS_vp_fam$mutation == "NO")]))/2, replace=FALSE)
MTS_vp_fam_train <- MTS_vp_fam[which(MTS_vp_fam$WGS_ID %in% c(as.character(train_samples_fam_mut),as.character(train_samples_fam_no_mut))),]
vp_IDs_fam_train <- unique(MTS_vp_fam_train$WGS_ID)

MTS_vp_fam_test_samples <- unique(MTS_vp_fam$WGS_ID[which(!(MTS_vp_fam$WGS_ID %in% vp_IDs_fam_train))])
MTS_vp_fam_test <- MTS_vp_fam[which(MTS_vp_fam$WGS_ID %in% MTS_vp_fam_test_samples),]
vp_IDs_fam_test <- unique(MTS_vp_fam_test$WGS_ID)

##Without family history
MTS_vp_indv <- MTS_all[which(MTS_all$cohort == "AJHG" & MTS_all$tumour_in_participant == "YES"),]
vp_IDs_indv <- unique(MTS_vp_indv$WGS_ID)

train_samples_indv_mut <- sample(unique(MTS_vp_indv$WGS_ID[which(MTS_vp_indv$mutation == "YES")]), length(unique(MTS_vp_indv$WGS_ID[which(MTS_vp_indv$mutation == "YES")]))/2, replace=FALSE)
train_samples_indv_no_mut <- sample(unique(MTS_vp_indv$WGS_ID[which(MTS_vp_indv$mutation == "NO")]), length(unique(MTS_vp_indv$WGS_ID[which(MTS_vp_indv$mutation == "NO")]))/2, replace=FALSE)

MTS_vp_indv_train <- MTS_vp_indv[which(MTS_vp_indv$WGS_ID %in% c(as.character(train_samples_indv_mut),as.character(train_samples_indv_no_mut))),]
vp_IDs_indv_train <- unique(MTS_vp_indv_train$WGS_ID)

MTS_vp_indv_test_samples <- unique(MTS_vp_indv$WGS_ID[which(!(MTS_vp_indv$WGS_ID %in% vp_IDs_indv_train))])
MTS_vp_indv_test <- MTS_vp_indv[which(MTS_vp_indv$WGS_ID %in% MTS_vp_indv_test_samples),]
vp_IDs_indv_test <- unique(MTS_vp_indv_test$WGS_ID)


##Calculate scores per proband

##All variables with family history TRAIN
option_1_score_fam_train <- lapply(vp_IDs_fam_train, function(x) sum(MTS_vp_fam_train$option_1_age[which(MTS_vp_fam_train$WGS_ID == x)]) + 
         sum(MTS_vp_fam_train$option_1_heritability[which(MTS_vp_fam_train$WGS_ID == x)]) +
         sum(MTS_vp_fam_train$option_1_incidence[which(MTS_vp_fam_train$WGS_ID == x)]))

option_2_score_fam_train <- lapply(vp_IDs_fam_train, function(x) sum(MTS_vp_fam_train$option_2_age[which(MTS_vp_fam_train$WGS_ID == x)]) + 
                                 sum(MTS_vp_fam_train$option_2_heritability[which(MTS_vp_fam_train$WGS_ID == x)]) +
                                 sum(MTS_vp_fam_train$option_2_incidence[which(MTS_vp_fam_train$WGS_ID == x)]))

option_3_score_fam_train <- lapply(vp_IDs_fam_train, function(x) sum(MTS_vp_fam_train$option_3_age[which(MTS_vp_fam_train$WGS_ID == x)]) + 
                                 sum(MTS_vp_fam_train$option_3_heritability[which(MTS_vp_fam_train$WGS_ID == x)]) +
                                 sum(MTS_vp_fam_train$option_3_incidence[which(MTS_vp_fam_train$WGS_ID == x)]))

option_4_score_fam_train <- lapply(vp_IDs_fam_train, function(x) sum(MTS_vp_fam_train$option_4_age[which(MTS_vp_fam_train$WGS_ID == x)]) + 
                                 sum(MTS_vp_fam_train$option_4_heritability[which(MTS_vp_fam_train$WGS_ID == x)]) +
                                 sum(MTS_vp_fam_train$option_4_incidence[which(MTS_vp_fam_train$WGS_ID == x)]))

option_5_score_fam_train <- lapply(vp_IDs_fam_train, function(x) sum(MTS_vp_fam_train$option_5_age[which(MTS_vp_fam_train$WGS_ID == x)]) + 
                                 sum(MTS_vp_fam_train$option_5_heritability[which(MTS_vp_fam_train$WGS_ID == x)]) +
                                 sum(MTS_vp_fam_train$option_5_incidence[which(MTS_vp_fam_train$WGS_ID == x)]))

option_6_score_fam_train <- lapply(vp_IDs_fam_train, function(x) sum(MTS_vp_fam_train$option_6_age[which(MTS_vp_fam_train$WGS_ID == x)]) + 
                                 sum(MTS_vp_fam_train$option_6_heritability[which(MTS_vp_fam_train$WGS_ID == x)]) +
                                 sum(MTS_vp_fam_train$option_6_incidence[which(MTS_vp_fam_train$WGS_ID == x)]))

original_score_fam_train <- lapply(vp_IDs_fam_train, function(x) sum(MTS_vp_fam_train$original[which(MTS_vp_fam_train$WGS_ID == x)]))
                                   
##Without incidence with family history TRAIN
option_1_score_fam_train_no_incid <- lapply(vp_IDs_fam_train, function(x) sum(MTS_vp_fam_train$option_1_age[which(MTS_vp_fam_train$WGS_ID == x)]) + 
         sum(MTS_vp_fam_train$option_1_heritability[which(MTS_vp_fam_train$WGS_ID == x)]))

option_2_score_fam_train_no_incid <- lapply(vp_IDs_fam_train, function(x) sum(MTS_vp_fam_train$option_2_age[which(MTS_vp_fam_train$WGS_ID == x)]) + 
                                              sum(MTS_vp_fam_train$option_2_heritability[which(MTS_vp_fam_train$WGS_ID == x)]))

option_3_score_fam_train_no_incid <- lapply(vp_IDs_fam_train, function(x) sum(MTS_vp_fam_train$option_3_age[which(MTS_vp_fam_train$WGS_ID == x)]) + 
                                              sum(MTS_vp_fam_train$option_3_heritability[which(MTS_vp_fam_train$WGS_ID == x)]))

option_4_score_fam_train_no_incid <- lapply(vp_IDs_fam_train, function(x) sum(MTS_vp_fam_train$option_4_age[which(MTS_vp_fam_train$WGS_ID == x)]) + 
                                              sum(MTS_vp_fam_train$option_4_heritability[which(MTS_vp_fam_train$WGS_ID == x)]))

option_5_score_fam_train_no_incid <- lapply(vp_IDs_fam_train, function(x) sum(MTS_vp_fam_train$option_5_age[which(MTS_vp_fam_train$WGS_ID == x)]) + 
                                              sum(MTS_vp_fam_train$option_5_heritability[which(MTS_vp_fam_train$WGS_ID == x)]))

option_6_score_fam_train_no_incid <- lapply(vp_IDs_fam_train, function(x) sum(MTS_vp_fam_train$option_6_age[which(MTS_vp_fam_train$WGS_ID == x)]) + 
                                              sum(MTS_vp_fam_train$option_6_heritability[which(MTS_vp_fam_train$WGS_ID == x)]))

##All variables without family history TRAIN
option_1_score_indv_train <- lapply(vp_IDs_indv_train, function(x) sum(MTS_vp_indv_train$option_1_age[which(MTS_vp_indv_train$WGS_ID == x)]) + 
                                     sum(MTS_vp_indv_train$option_1_heritability[which(MTS_vp_indv_train$WGS_ID == x)]) +
                                     sum(MTS_vp_indv_train$option_1_incidence[which(MTS_vp_indv_train$WGS_ID == x)]))

option_2_score_indv_train <- lapply(vp_IDs_indv_train, function(x) sum(MTS_vp_indv_train$option_2_age[which(MTS_vp_indv_train$WGS_ID == x)]) + 
                                     sum(MTS_vp_indv_train$option_2_heritability[which(MTS_vp_indv_train$WGS_ID == x)]) +
                                     sum(MTS_vp_indv_train$option_2_incidence[which(MTS_vp_indv_train$WGS_ID == x)]))

option_3_score_indv_train <- lapply(vp_IDs_indv_train, function(x) sum(MTS_vp_indv_train$option_3_age[which(MTS_vp_indv_train$WGS_ID == x)]) + 
                                     sum(MTS_vp_indv_train$option_3_heritability[which(MTS_vp_indv_train$WGS_ID == x)]) +
                                     sum(MTS_vp_indv_train$option_3_incidence[which(MTS_vp_indv_train$WGS_ID == x)]))

option_4_score_indv_train <- lapply(vp_IDs_indv_train, function(x) sum(MTS_vp_indv_train$option_4_age[which(MTS_vp_indv_train$WGS_ID == x)]) + 
                                     sum(MTS_vp_indv_train$option_4_heritability[which(MTS_vp_indv_train$WGS_ID == x)]) +
                                     sum(MTS_vp_indv_train$option_4_incidence[which(MTS_vp_indv_train$WGS_ID == x)]))

option_5_score_indv_train <- lapply(vp_IDs_indv_train, function(x) sum(MTS_vp_indv_train$option_5_age[which(MTS_vp_indv_train$WGS_ID == x)]) + 
                                     sum(MTS_vp_indv_train$option_5_heritability[which(MTS_vp_indv_train$WGS_ID == x)]) +
                                     sum(MTS_vp_indv_train$option_5_incidence[which(MTS_vp_indv_train$WGS_ID == x)]))

option_6_score_indv_train <- lapply(vp_IDs_indv_train, function(x) sum(MTS_vp_indv_train$option_6_age[which(MTS_vp_indv_train$WGS_ID == x)]) + 
                                     sum(MTS_vp_indv_train$option_6_heritability[which(MTS_vp_indv_train$WGS_ID == x)]) +
                                     sum(MTS_vp_indv_train$option_6_incidence[which(MTS_vp_indv_train$WGS_ID == x)]))


original_score_indv_train <- lapply(vp_IDs_indv_train, function(x) sum(MTS_vp_indv_train$original[which(MTS_vp_indv_train$WGS_ID == x)]))

##Without incidence without family history TRAIN
option_1_score_indv_train_no_incid <- lapply(vp_IDs_indv_train, function(x) sum(MTS_vp_indv_train$option_1_age[which(MTS_vp_indv_train$WGS_ID == x)]) + 
                                              sum(MTS_vp_indv_train$option_1_heritability[which(MTS_vp_indv_train$WGS_ID == x)]))

option_2_score_indv_train_no_incid <- lapply(vp_IDs_indv_train, function(x) sum(MTS_vp_indv_train$option_2_age[which(MTS_vp_indv_train$WGS_ID == x)]) + 
                                              sum(MTS_vp_indv_train$option_2_heritability[which(MTS_vp_indv_train$WGS_ID == x)]))

option_3_score_indv_train_no_incid <- lapply(vp_IDs_indv_train, function(x) sum(MTS_vp_indv_train$option_3_age[which(MTS_vp_indv_train$WGS_ID == x)]) + 
                                              sum(MTS_vp_indv_train$option_3_heritability[which(MTS_vp_indv_train$WGS_ID == x)]))

option_4_score_indv_train_no_incid <- lapply(vp_IDs_indv_train, function(x) sum(MTS_vp_indv_train$option_4_age[which(MTS_vp_indv_train$WGS_ID == x)]) + 
                                              sum(MTS_vp_indv_train$option_4_heritability[which(MTS_vp_indv_train$WGS_ID == x)]))

option_5_score_indv_train_no_incid <- lapply(vp_IDs_indv_train, function(x) sum(MTS_vp_indv_train$option_5_age[which(MTS_vp_indv_train$WGS_ID == x)]) + 
                                              sum(MTS_vp_indv_train$option_5_heritability[which(MTS_vp_indv_train$WGS_ID == x)]))

option_6_score_indv_train_no_incid <- lapply(vp_IDs_indv_train, function(x) sum(MTS_vp_indv_train$option_6_age[which(MTS_vp_indv_train$WGS_ID == x)]) + 
                                              sum(MTS_vp_indv_train$option_6_heritability[which(MTS_vp_indv_train$WGS_ID == x)]))

##Make table with each proband's sample number and assigned total score under different options

##With family history
vp_fam_train <- as.data.frame(cbind(as.character(vp_IDs_fam_train), 
      #
      option_1_score_fam_train,
      option_2_score_fam_train, 
      option_3_score_fam_train, 
      option_4_score_fam_train, 
      option_5_score_fam_train, 
      option_6_score_fam_train, 
      #
      option_1_score_fam_train_no_incid,
      option_2_score_fam_train_no_incid,
      option_3_score_fam_train_no_incid,
      option_4_score_fam_train_no_incid,
      option_5_score_fam_train_no_incid,
      option_6_score_fam_train_no_incid,
      
      original_score_fam_train
      
      ))

colnames(vp_fam_train)[1] <- "WGS_ID"

##Add pathogenic variant status to table
mut_status_vp_fam_train <- unlist(lapply(vp_fam_train$WGS_ID, function(x) MTS_all$mutation[which.max(MTS_all$WGS_ID == x)]))
vp_fam_train <- data.frame(vp_fam_train, mut_status_vp_fam_train)

##Without family history
vp_indv_train <-as.data.frame(cbind(as.character(vp_IDs_indv_train), 
      #
      option_1_score_indv_train,
      option_2_score_indv_train, 
      option_3_score_indv_train, 
      option_4_score_indv_train, 
      option_5_score_indv_train, 
      option_6_score_indv_train, 
      #
      option_1_score_indv_train_no_incid,
      option_2_score_indv_train_no_incid,
      option_3_score_indv_train_no_incid,
      option_4_score_indv_train_no_incid,
      option_5_score_indv_train_no_incid,
      option_6_score_indv_train_no_incid,
      
      original_score_indv_train
      
      ))

colnames(vp_indv_train)[1] <- "WGS_ID"

##Add pathogenic variant status to table
mut_status_vp_indv_train <- unlist(lapply(vp_indv_train$WGS_ID, function(x) MTS_all$mutation[which.max(MTS_all$WGS_ID == x)]))
vp_indv_train <- data.frame(vp_indv_train, mut_status_vp_indv_train)

##Perform logistic regressions on training sets

##With family history with incidence
train_option_fam_1_mod <- glm(formula = as.numeric(mut_status_vp_fam_train) -1~unlist(option_1_score_fam_train), data = vp_fam_train, family = binomial(link = "logit"))
library(pROC)
train_option_fam_1_roc <- roc(as.numeric(mut_status_vp_fam_train) -1 ~ unlist(option_1_score_fam_train), data = vp_fam_train)
model_outputs_text <- append(model_outputs_text, capture.output(summary(train_option_fam_1_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(train_option_fam_1_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(train_option_fam_1_roc))

train_option_fam_2_mod <- glm(formula = as.numeric(mut_status_vp_fam_train) -1~unlist(option_2_score_fam_train), data = vp_fam_train, family = binomial(link = "logit"))
library(pROC)
train_option_fam_2_roc <- roc(as.numeric(mut_status_vp_fam_train) -1 ~ unlist(option_2_score_fam_train), data = vp_fam_train)
model_outputs_text <- append(model_outputs_text, capture.output(summary(train_option_fam_2_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(train_option_fam_2_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(train_option_fam_2_roc))

train_option_fam_3_mod <- glm(formula = as.numeric(mut_status_vp_fam_train) -1~unlist(option_3_score_fam_train), data = vp_fam_train, family = binomial(link = "logit"))
library(pROC)
train_option_fam_3_roc <- roc(as.numeric(mut_status_vp_fam_train) -1 ~ unlist(option_3_score_fam_train), data = vp_fam_train)
model_outputs_text <- append(model_outputs_text, capture.output(summary(train_option_fam_3_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(train_option_fam_3_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(train_option_fam_3_roc))

train_option_fam_4_mod <- glm(formula = as.numeric(mut_status_vp_fam_train) -1~unlist(option_4_score_fam_train), data = vp_fam_train, family = binomial(link = "logit"))
library(pROC)
train_option_fam_4_roc <- roc(as.numeric(mut_status_vp_fam_train) -1 ~ unlist(option_4_score_fam_train), data = vp_fam_train)
model_outputs_text <- append(model_outputs_text, capture.output(summary(train_option_fam_4_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(train_option_fam_4_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(train_option_fam_4_roc))

train_option_fam_5_mod <- glm(formula = as.numeric(mut_status_vp_fam_train) -1~unlist(option_5_score_fam_train), data = vp_fam_train, family = binomial(link = "logit"))
library(pROC)
train_option_fam_5_roc <- roc(as.numeric(mut_status_vp_fam_train) -1 ~ unlist(option_5_score_fam_train), data = vp_fam_train)
model_outputs_text <- append(model_outputs_text, capture.output(summary(train_option_fam_5_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(train_option_fam_5_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(train_option_fam_5_roc))

train_option_fam_6_mod <- glm(formula = as.numeric(mut_status_vp_fam_train) -1~unlist(option_6_score_fam_train), data = vp_fam_train, family = binomial(link = "logit"))
library(pROC)
train_option_fam_6_roc <- roc(as.numeric(mut_status_vp_fam_train) -1 ~ unlist(option_6_score_fam_train), data = vp_fam_train)
model_outputs_text <- append(model_outputs_text, capture.output(summary(train_option_fam_6_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(train_option_fam_6_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(train_option_fam_6_roc))

train_fam_original_mod <- glm(formula = as.numeric(mut_status_vp_fam_train) -1~unlist(original_score_fam_train), data = vp_fam_train, family = binomial(link = "logit"))
library(pROC)
train_fam_original_roc <- roc(as.numeric(mut_status_vp_fam_train) -1 ~ unlist(original_score_fam_train), data = vp_fam_train)
model_outputs_text <- append(model_outputs_text, capture.output(summary(train_fam_original_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(train_fam_original_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(train_fam_original_roc))

##With family history without incidence
train_option_fam_1_no_incid_mod <- glm(formula = as.numeric(mut_status_vp_fam_train) -1~unlist(option_1_score_fam_train_no_incid), data = vp_fam_train, family = binomial(link = "logit"))
library(pROC)
train_option_fam_1_no_incid_roc <- roc(as.numeric(mut_status_vp_fam_train) -1 ~ unlist(option_1_score_fam_train_no_incid), data = vp_fam_train)
model_outputs_text <- append(model_outputs_text, capture.output(summary(train_option_fam_1_no_incid_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(train_option_fam_1_no_incid_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(train_option_fam_1_no_incid_roc))

train_option_fam_2_no_incid_mod <- glm(formula = as.numeric(mut_status_vp_fam_train) -1~unlist(option_2_score_fam_train_no_incid), data = vp_fam_train, family = binomial(link = "logit"))
library(pROC)
train_option_fam_2_no_incid_roc <- roc(as.numeric(mut_status_vp_fam_train) -1 ~ unlist(option_2_score_fam_train_no_incid), data = vp_fam_train)
model_outputs_text <- append(model_outputs_text, capture.output(summary(train_option_fam_2_no_incid_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(train_option_fam_2_no_incid_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(train_option_fam_2_no_incid_roc))

train_option_fam_3_no_incid_mod <- glm(formula = as.numeric(mut_status_vp_fam_train) -1~unlist(option_3_score_fam_train_no_incid), data = vp_fam_train, family = binomial(link = "logit"))
library(pROC)
train_option_fam_3_no_incid_roc <- roc(as.numeric(mut_status_vp_fam_train) -1 ~ unlist(option_3_score_fam_train_no_incid), data = vp_fam_train)
model_outputs_text <- append(model_outputs_text, capture.output(summary(train_option_fam_3_no_incid_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(train_option_fam_3_no_incid_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(train_option_fam_3_no_incid_roc))

train_option_fam_4_no_incid_mod <- glm(formula = as.numeric(mut_status_vp_fam_train) -1~unlist(option_4_score_fam_train_no_incid), data = vp_fam_train, family = binomial(link = "logit"))
library(pROC)
train_option_fam_4_no_incid_roc <- roc(as.numeric(mut_status_vp_fam_train) -1 ~ unlist(option_4_score_fam_train_no_incid), data = vp_fam_train)
model_outputs_text <- append(model_outputs_text, capture.output(summary(train_option_fam_4_no_incid_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(train_option_fam_4_no_incid_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(train_option_fam_4_no_incid_roc))

train_option_fam_5_no_incid_mod <- glm(formula = as.numeric(mut_status_vp_fam_train) -1~unlist(option_5_score_fam_train_no_incid), data = vp_fam_train, family = binomial(link = "logit"))
library(pROC)
train_option_fam_5_no_incid_roc <- roc(as.numeric(mut_status_vp_fam_train) -1 ~ unlist(option_5_score_fam_train_no_incid), data = vp_fam_train)
model_outputs_text <- append(model_outputs_text, capture.output(summary(train_option_fam_5_no_incid_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(train_option_fam_5_no_incid_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(train_option_fam_5_no_incid_roc))

train_option_fam_6_no_incid_mod <- glm(formula = as.numeric(mut_status_vp_fam_train) -1~unlist(option_6_score_fam_train_no_incid), data = vp_fam_train, family = binomial(link = "logit"))
library(pROC)
train_option_fam_6_no_incid_roc <- roc(as.numeric(mut_status_vp_fam_train) -1 ~ unlist(option_6_score_fam_train_no_incid), data = vp_fam_train)
model_outputs_text <- append(model_outputs_text, capture.output(summary(train_option_fam_6_no_incid_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(train_option_fam_6_no_incid_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(train_option_fam_6_no_incid_roc))

##Without family history with incidence
train_option_indv_1_mod <- glm(formula = as.numeric(mut_status_vp_indv_train) -1~unlist(option_1_score_indv_train), data = vp_indv_train, family = binomial(link = "logit"))
library(pROC)
train_option_indv_1_roc <- roc(as.numeric(mut_status_vp_indv_train) -1 ~ unlist(option_1_score_indv_train), data = vp_indv_train)
model_outputs_text <- append(model_outputs_text, capture.output(summary(train_option_indv_1_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(train_option_indv_1_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(train_option_indv_1_roc))

train_option_indv_2_mod <- glm(formula = as.numeric(mut_status_vp_indv_train) -1~unlist(option_2_score_indv_train), data = vp_indv_train, family = binomial(link = "logit"))
library(pROC)
train_option_indv_2_roc <- roc(as.numeric(mut_status_vp_indv_train) -1 ~ unlist(option_2_score_indv_train), data = vp_indv_train)
model_outputs_text <- append(model_outputs_text, capture.output(summary(train_option_indv_2_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(train_option_indv_2_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(train_option_indv_2_roc))

train_option_indv_3_mod <- glm(formula = as.numeric(mut_status_vp_indv_train) -1~unlist(option_3_score_indv_train), data = vp_indv_train, family = binomial(link = "logit"))
library(pROC)
train_option_indv_3_roc <- roc(as.numeric(mut_status_vp_indv_train) -1 ~ unlist(option_3_score_indv_train), data = vp_indv_train)
model_outputs_text <- append(model_outputs_text, capture.output(summary(train_option_indv_3_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(train_option_indv_3_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(train_option_indv_3_roc))

train_option_indv_4_mod <- glm(formula = as.numeric(mut_status_vp_indv_train) -1~unlist(option_4_score_indv_train), data = vp_indv_train, family = binomial(link = "logit"))
library(pROC)
train_option_indv_4_roc <- roc(as.numeric(mut_status_vp_indv_train) -1 ~ unlist(option_4_score_indv_train), data = vp_indv_train)
model_outputs_text <- append(model_outputs_text, capture.output(summary(train_option_indv_4_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(train_option_indv_4_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(train_option_indv_4_roc))

train_option_indv_5_mod <- glm(formula = as.numeric(mut_status_vp_indv_train) -1~unlist(option_5_score_indv_train), data = vp_indv_train, family = binomial(link = "logit"))
library(pROC)
train_option_indv_5_roc <- roc(as.numeric(mut_status_vp_indv_train) -1 ~ unlist(option_5_score_indv_train), data = vp_indv_train)
model_outputs_text <- append(model_outputs_text, capture.output(summary(train_option_indv_5_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(train_option_indv_5_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(train_option_indv_5_roc))

train_option_indv_6_mod <- glm(formula = as.numeric(mut_status_vp_indv_train) -1~unlist(option_6_score_indv_train), data = vp_indv_train, family = binomial(link = "logit"))
library(pROC)
train_option_indv_6_roc <- roc(as.numeric(mut_status_vp_indv_train) -1 ~ unlist(option_6_score_indv_train), data = vp_indv_train)
model_outputs_text <- append(model_outputs_text, capture.output(summary(train_option_indv_6_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(train_option_indv_6_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(train_option_indv_6_roc))

train_indv_original_mod <- glm(formula = as.numeric(mut_status_vp_indv_train) -1~unlist(original_score_indv_train), data = vp_indv_train, family = binomial(link = "logit"))
library(pROC)
train_indv_original_roc <- roc(as.numeric(mut_status_vp_indv_train) -1 ~ unlist(original_score_indv_train), data = vp_indv_train)
model_outputs_text <- append(model_outputs_text, capture.output(summary(train_indv_original_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(train_indv_original_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(train_indv_original_roc))

##Without family history without incidence
train_option_indv_1_no_incid_mod <- glm(formula = as.numeric(mut_status_vp_indv_train) -1~unlist(option_1_score_indv_train_no_incid), data = vp_indv_train, family = binomial(link = "logit"))
library(pROC)
train_option_indv_1_no_incid_roc <- roc(as.numeric(mut_status_vp_indv_train) -1 ~ unlist(option_1_score_indv_train_no_incid), data = vp_indv_train)
model_outputs_text <- append(model_outputs_text, capture.output(summary(train_option_indv_1_no_incid_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(train_option_indv_1_no_incid_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(train_option_indv_1_no_incid_roc))

train_option_indv_2_no_incid_mod <- glm(formula = as.numeric(mut_status_vp_indv_train) -1~unlist(option_2_score_indv_train_no_incid), data = vp_indv_train, family = binomial(link = "logit"))
library(pROC)
train_option_indv_2_no_incid_roc <- roc(as.numeric(mut_status_vp_indv_train) -1 ~ unlist(option_2_score_indv_train_no_incid), data = vp_indv_train)
model_outputs_text <- append(model_outputs_text, capture.output(summary(train_option_indv_2_no_incid_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(train_option_indv_2_no_incid_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(train_option_indv_2_no_incid_roc))

train_option_indv_3_no_incid_mod <- glm(formula = as.numeric(mut_status_vp_indv_train) -1~unlist(option_3_score_indv_train_no_incid), data = vp_indv_train, family = binomial(link = "logit"))
library(pROC)
train_option_indv_3_no_incid_roc <- roc(as.numeric(mut_status_vp_indv_train) -1 ~ unlist(option_3_score_indv_train_no_incid), data = vp_indv_train)
model_outputs_text <- append(model_outputs_text, capture.output(summary(train_option_indv_3_no_incid_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(train_option_indv_3_no_incid_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(train_option_indv_3_no_incid_roc))

train_option_indv_4_no_incid_mod <- glm(formula = as.numeric(mut_status_vp_indv_train) -1~unlist(option_4_score_indv_train_no_incid), data = vp_indv_train, family = binomial(link = "logit"))
library(pROC)
train_option_indv_4_no_incid_roc <- roc(as.numeric(mut_status_vp_indv_train) -1 ~ unlist(option_4_score_indv_train_no_incid), data = vp_indv_train)
model_outputs_text <- append(model_outputs_text, capture.output(summary(train_option_indv_4_no_incid_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(train_option_indv_4_no_incid_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(train_option_indv_4_no_incid_roc))

train_option_indv_5_no_incid_mod <- glm(formula = as.numeric(mut_status_vp_indv_train) -1~unlist(option_5_score_indv_train_no_incid), data = vp_indv_train, family = binomial(link = "logit"))
library(pROC)
train_option_indv_5_no_incid_roc <- roc(as.numeric(mut_status_vp_indv_train) -1 ~ unlist(option_5_score_indv_train_no_incid), data = vp_indv_train)
model_outputs_text <- append(model_outputs_text, capture.output(summary(train_option_indv_5_no_incid_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(train_option_indv_5_no_incid_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(train_option_indv_5_no_incid_roc))

train_option_indv_6_no_incid_mod <- glm(formula = as.numeric(mut_status_vp_indv_train) -1~unlist(option_6_score_indv_train_no_incid), data = vp_indv_train, family = binomial(link = "logit"))
library(pROC)
train_option_indv_6_no_incid_roc <- roc(as.numeric(mut_status_vp_indv_train) -1 ~ unlist(option_6_score_indv_train_no_incid), data = vp_indv_train)
model_outputs_text <- append(model_outputs_text, capture.output(summary(train_option_indv_6_no_incid_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(train_option_indv_6_no_incid_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(train_option_indv_6_no_incid_roc))

##Output performance statistics for models
write(model_outputs_text, "model_outputs_text")

##############################################################
##Assessment of best performing model (Option 3) on test set##
##############################################################

##Without family history

##Calculate scores for test set individuals and assign pathogenic variant status
option_3_score_indv_test <- lapply(vp_IDs_indv_test, function(x) sum(MTS_vp_indv_test$option_3_age[which(MTS_vp_indv_test$WGS_ID == x)]) + 
                                      sum(MTS_vp_indv_test$option_3_heritability[which(MTS_vp_indv_test$WGS_ID == x)]) +
                                      sum(MTS_vp_indv_test$option_3_incidence[which(MTS_vp_indv_test$WGS_ID == x)]))
vp_indv_test <- as.data.frame(cbind(as.character(vp_IDs_indv_test), option_3_score_indv_test))
colnames(vp_indv_test)[1] <- "WGS_ID"
mut_status_vp_indv_test <- unlist(lapply(vp_indv_test$WGS_ID, function(x) MTS_all$mutation[which.max(MTS_all$WGS_ID == x)]))
vp_indv_test <- data.frame(vp_indv_test, mut_status_vp_indv_test)

##Perform logisitic regression and assess performance
test_indv_option_3_mod <- glm(formula = as.numeric(mut_status_vp_indv_test) -1~unlist(option_3_score_indv_test), data = vp_indv_test, family = binomial(link = "logit"))
library(pROC)
test_indv_option_3_roc <- roc(as.numeric(mut_status_vp_indv_test) -1 ~ unlist(option_3_score_indv_test), data = vp_indv_test)
model_outputs_text <- append(model_outputs_text, capture.output(summary(test_indv_option_3_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(test_indv_option_3_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(test_indv_option_3_roc))


##With family history witout incidence

##Calculate scores for test set individuals and assign pathogenic variant status
option_3_score_fam_test_no_incid <- lapply(vp_IDs_fam_test, function(x) sum(MTS_vp_fam_test$option_3_age[which(MTS_vp_fam_test$WGS_ID == x)]) + 
                                              sum(MTS_vp_fam_test$option_3_heritability[which(MTS_vp_fam_test$WGS_ID == x)]))


vp_fam_test <- as.data.frame(cbind(as.character(vp_IDs_fam_test), option_3_score_fam_test_no_incid))
colnames(vp_fam_test)[1] <- "WGS_ID"
mut_status_vp_fam_test <- unlist(lapply(vp_fam_test$WGS_ID, function(x) MTS_all$mutation[which.max(MTS_all$WGS_ID == x)]))
vp_fam_test <- data.frame(vp_fam_test, mut_status_vp_fam_test)

##Perform logisitic regression and assess performance
test_fam_option_3_no_incid_mod <- glm(formula = as.numeric(mut_status_vp_fam_test) -1~unlist(option_3_score_fam_test_no_incid), data = vp_fam_test, family = binomial(link = "logit"))
library(pROC)
test_fam_option_3_no_incid_roc <- roc(as.numeric(mut_status_vp_fam_test) -1 ~ unlist(option_3_score_fam_test_no_incid), data = vp_fam_test)
model_outputs_text <- append(model_outputs_text, capture.output(summary(test_fam_option_3_no_incid_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(test_fam_option_3_no_incid_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(test_fam_option_3_no_incid_roc))

##########################################################################################################
##Also test best performing model on previous MPT series(individual information only. No family history)##
##########################################################################################################
##Best performing score (not incorporating family history was option 3

##Redefine table to include cases from previous dataset
MTS_EJHG_indv <- MTS[which(MTS$cohort == "EJHG" & 
                             MTS$tumour_in_participant == "YES" & 
                             MTS$exclude == "NO" & 
                             MTS$heritability_exclude == "NO" & 
                             MTS$incidence_exclude == "NO"),]

EJHG_IDs_indv <- unique(MTS_EJHG_indv$WGS_ID)

##Assign appropriate scoring option
option_3_age <- rep("empty", nrow(MTS_EJHG_indv))
option_3_age[which(MTS_EJHG_indv$age_band == "<30")] <- 27
option_3_age[which(MTS_EJHG_indv$age_band == "30-44")] <- 9
option_3_age[which(MTS_EJHG_indv$age_band == "45-59")] <- 3
option_3_age[which(MTS_EJHG_indv$age_band == ">59")] <- 1

option_3_heritability <- rep("empty", nrow(MTS_EJHG_indv))
option_3_heritability[which(MTS_EJHG_indv$heritability_band == "0-25")] <- 1
option_3_heritability[which(MTS_EJHG_indv$heritability_band == "26-50")] <- 3
option_3_heritability[which(MTS_EJHG_indv$heritability_band == "51-75")] <- 9
option_3_heritability[which(MTS_EJHG_indv$heritability_band == "76-100")] <- 27

option_3_incidence <- rep("empty", nrow(MTS_EJHG_indv))
option_3_incidence[which(MTS_EJHG_indv$incidence_band == "0-6")] <- 27
option_3_incidence[which(MTS_EJHG_indv$incidence_band == "6.1-28")] <- 9
option_3_incidence[which(MTS_EJHG_indv$incidence_band == "28.1>50")] <- 3
option_3_incidence[which(MTS_EJHG_indv$incidence_band == ">50")] <- 1


MTS_EJHG_indv <- cbind(MTS_EJHG_indv, as.numeric(option_3_age), as.numeric(option_3_heritability), as.numeric(option_3_incidence))
colnames(MTS_EJHG_indv)[37:39] <- c("option_3_age", "option_3_heritability", "option_3_incidence")

##Calculate scores for individuals and add pathogenic variant status
option_3_score_EJHG_test <- lapply(EJHG_IDs_indv, function(x) sum(MTS_EJHG_indv$option_3_age[which(MTS_EJHG_indv$WGS_ID == x)]) + 
                                     sum(MTS_EJHG_indv$option_3_heritability[which(MTS_EJHG_indv$WGS_ID == x)]) +
                                     sum(MTS_EJHG_indv$option_3_incidence[which(MTS_EJHG_indv$WGS_ID == x)]))
vp_EJHG_test <- as.data.frame(cbind(as.character(EJHG_IDs_indv), option_3_score_EJHG_test))
colnames(vp_EJHG_test)[1] <- "WGS_ID"
mut_status_EJHG <- unlist(lapply(vp_EJHG_test$WGS_ID, function(x) MTS_EJHG_indv$mutation[which.max(MTS_EJHG_indv$WGS_ID == x)]))
vp_EJHG_test <- data.frame(vp_EJHG_test, mut_status_EJHG)

##Perform logisitic regression and assess performance
test_EJHG_option_3_mod <- glm(formula = as.numeric(mut_status_EJHG) -1~unlist(option_3_score_EJHG_test), data = vp_EJHG_test, family = binomial(link = "logit"))
library(pROC)
test_EJHG_option_3_roc <- roc(as.numeric(mut_status_EJHG) -1 ~ unlist(option_3_score_EJHG_test), data = vp_EJHG_test)
model_outputs_text <- append(model_outputs_text, capture.output(summary(test_EJHG_option_3_mod)))
model_outputs_text <- append(model_outputs_text, capture.output(anova(test_EJHG_option_3_mod, test = "Chisq")))
model_outputs_text <- append(model_outputs_text, capture.output(test_EJHG_option_3_roc))

##Output performance statistics for models
write(model_outputs_text, "model_outputs_text")
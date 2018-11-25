##############################################################################################################################
##R Script to analyse tumour occurrences in referral and registry based multiple primary tumour series. Script run in R Studio##
##############################################################################################################################

##Read in reference file to to make tumour names consistent between series
tumour_name_conversions <- read.csv("tumour_name_conversions.csv", header = T)

##################
##MPT SERIES######
##################

##Read in data table containing tumour information
vp_tumours <- read.csv("vp_tumours.csv")

##Extract combinations from table
comb_all <- data.frame()

comb_1 <- paste(vp_tumours$Single_word_1,vp_tumours$Single_word_2, sep = "-")
comb_all <- append(comb_1, comb_all)
comb_2 <- paste(vp_tumours$Single_word_1,vp_tumours$Single_word_3, sep = "-")
comb_all <- append(comb_2, comb_all)
comb_3 <- paste(vp_tumours$Single_word_1,vp_tumours$Single_word_4, sep = "-")
comb_all <- append(comb_3, comb_all)
comb_4 <- paste(vp_tumours$Single_word_1,vp_tumours$Single_word_5, sep = "-")
comb_all <- append(comb_4, comb_all)
comb_5 <- paste(vp_tumours$Single_word_1,vp_tumours$Single_word_6, sep = "-")
comb_all <- append(comb_5, comb_all)

comb_6 <- paste(vp_tumours$Single_word_2,vp_tumours$Single_word_3, sep = "-")
comb_all <- append(comb_6, comb_all)
comb_7 <- paste(vp_tumours$Single_word_2,vp_tumours$Single_word_4, sep = "-")
comb_all <- append(comb_7, comb_all)
comb_8 <- paste(vp_tumours$Single_word_2,vp_tumours$Single_word_5, sep = "-")
comb_all <- append(comb_8, comb_all)
comb_9 <- paste(vp_tumours$Single_word_2,vp_tumours$Single_word_6, sep = "-")
comb_all <- append(comb_9, comb_all)

comb_10 <- paste(vp_tumours$Single_word_3,vp_tumours$Single_word_4, sep = "-")
comb_all <- append(comb_10, comb_all)
comb_11 <- paste(vp_tumours$Single_word_3,vp_tumours$Single_word_5, sep = "-")
comb_all <- append(comb_11, comb_all)
comb_12 <- paste(vp_tumours$Single_word_3,vp_tumours$Single_word_6, sep = "-")
comb_all <- append(comb_12, comb_all)

comb_13 <- paste(vp_tumours$Single_word_4,vp_tumours$Single_word_5, sep = "-")
comb_all <- append(comb_13, comb_all)
comb_14 <- paste(vp_tumours$Single_word_4,vp_tumours$Single_word_6, sep = "-")
comb_all <- append(comb_14, comb_all)

comb_15 <- paste(vp_tumours$Single_word_5,vp_tumours$Single_word_6, sep = "-")
comb_all <- append(comb_15, comb_all)

comb_all <- unlist(comb_all)
comb_all_split <- strsplit(comb_all, "-")

comb_all_logical <- data.frame()

for (i in 1:length(comb_all_split)) {
  
  a <- isTRUE(length(comb_all_split[[i]]) > 1)
  comb_all_logical <- append(comb_all_logical,a)
  
}

comb_all_clean <- comb_all_split[which(comb_all_logical == TRUE)]

comb_all_table <- as.data.frame(t(data.frame(comb_all_clean, header = F)))
rownames(comb_all_table) <- seq(1:nrow(comb_all_table))
comb_all_table$V1 <- as.character(comb_all_table$V1)
comb_all_table$V2 <- as.character(comb_all_table$V2)

##Convert tumour names
tumour_name_conversions$original <- as.character(tumour_name_conversions$original)
tumour_name_conversions$new <- as.character(tumour_name_conversions$new)

for (i in 1:nrow(tumour_name_conversions)) {
  
  comb_all_table$V1  <- gsub(tumour_name_conversions$original[i], tumour_name_conversions$new[i], comb_all_table$V1)
  
}

for (i in 1:nrow(tumour_name_conversions)) {
  
  comb_all_table$V2  <- gsub(tumour_name_conversions$original[i], tumour_name_conversions$new[i], comb_all_table$V2)
  
}

##Create table of combinations
comb_all_table_no_dup <- comb_all_table[which(comb_all_table$V1 != comb_all_table$V2),]

comb_all_table_no_dup_sorted <- data.frame(t(apply(comb_all_table_no_dup, 1, sort)))
comb_all_table_no_dup_sorted$X1 <- as.character(comb_all_table_no_dup_sorted$X1)
comb_all_table_no_dup_sorted$X2 <- as.character(comb_all_table_no_dup_sorted$X2)

count <- rep(1, nrow(comb_all_table_no_dup_sorted))
combination <- paste(comb_all_table_no_dup_sorted$X1, comb_all_table_no_dup_sorted$X2, sep = "-")

comb_all_table_no_dup_sorted <- data.frame(cbind(comb_all_table_no_dup_sorted$X1, comb_all_table_no_dup_sorted$X2, combination, count))
comb_all_table_no_dup_sorted$count <- as.integer(comb_all_table_no_dup_sorted$count)

##Create combination count table
combination_count <- aggregate(count~combination, data = comb_all_table_no_dup_sorted, sum)
proportion <- (100/sum(combination_count$count)) * (combination_count$count)
combination_count <- data.frame(combination_count, proportion)
combination_count <- combination_count[order(combination_count$count, decreasing = T),]

##Output table
write.csv(combination_count, "combination_count_MPT.csv")

##Check which combinations make up more than 1 percent total
combination_count_greater_one_percent_total <- combination_count[which(combination_count$proportion > 1),]

##Create chord diagram of combinations
tumors <- strsplit(as.character(combination_count$combination), "-")
tumor_1 <- sapply(tumors,function(x) x[1])
tumor_2 <- sapply(tumors,function(x) x[2])

combination_count_for_circos <- data.frame(tumor_1, tumor_2, combination_count$count, combination_count$proportion)
colnames(combination_count_for_circos) <- c("Tumor_1", "Tumor_2", "Count", "Proportion")

##Calculate characteristics of combinations before plotting chord diagram
top_5 <- unique(tumour_name_conversions$new[which(tumour_name_conversions$top_5 == "YES")])

tumour_1_or_2_not_top_5 <- !(combination_count_for_circos$Tumor_1 %in% top_5) | !(combination_count_for_circos$Tumor_2 %in% top_5)
tumour_1_and_2_not_top_5 <- !(combination_count_for_circos$Tumor_1 %in% top_5) & !(combination_count_for_circos$Tumor_2 %in% top_5)

tumour_1_or_2_not_top_5_count <- sum(combination_count_for_circos$Count[tumour_1_or_2_not_top_5])
tumour_1_or_2_not_top_5_prop <- sum(combination_count_for_circos$Count[tumour_1_or_2_not_top_5])/sum(combination_count_for_circos$Count) * 100

tumour_1_and_2_not_top_5_count <- sum(combination_count_for_circos$Count[tumour_1_and_2_not_top_5])
tumour_1_and_2_not_top_5_prop <- sum(combination_count_for_circos$Count[tumour_1_and_2_not_top_5])/sum(combination_count_for_circos$Count) * 100

high_herit <- unique(tumour_name_conversions$new[which(tumour_name_conversions$heritibility > 20)]) ## dummy number at present

tumour_1_or_2_high_herit <- !(combination_count_for_circos$Tumor_1 %in% high_herit) | !(combination_count_for_circos$Tumor_2 %in% high_herit)
tumour_1_and_2_high_herit <- !(combination_count_for_circos$Tumor_1 %in% high_herit) & !(combination_count_for_circos$Tumor_2 %in% high_herit)

tumour_1_or_2_high_herit_count <- sum(combination_count_for_circos$Count[tumour_1_or_2_high_herit])
tumour_1_or_2_high_herit_prop <- sum(combination_count_for_circos$Count[tumour_1_or_2_high_herit])/sum(combination_count_for_circos$Count) * 100

tumour_1_and_2_high_herit_count <- sum(combination_count_for_circos$Count[tumour_1_and_2_high_herit])
tumour_1_and_2_high_herit_prop <- sum(combination_count_for_circos$Count[tumour_1_and_2_high_herit])/sum(combination_count_for_circos$Count) * 100

tumour_1_or_2_not_top_5_and_high_herit_count <- sum(combination_count_for_circos$Count[tumour_1_or_2_not_top_5 == T & tumour_1_or_2_high_herit == T])
tumour_1_or_2_not_top_5_and_high_herit_prop <- sum(combination_count_for_circos$Count[tumour_1_or_2_not_top_5 == T & tumour_1_or_2_high_herit == T])/sum(combination_count_for_circos$Count) * 100

tumour_1_and_2_not_top_5_and_high_herit_count <- sum(combination_count_for_circos$Count[tumour_1_and_2_not_top_5 == T & tumour_1_and_2_high_herit == T])
tumour_1_and_2_not_top_5_and_high_herit_prop <- sum(combination_count_for_circos$Count[tumour_1_and_2_not_top_5 == T & tumour_1_and_2_high_herit == T])/sum(combination_count_for_circos$Count) * 100

##Output measures
number_discordant_combinations <- sum(combination_count_for_circos$Count)

combination_table_MPT <- as.data.frame(rbind(
number_discordant_combinations,
tumour_1_or_2_not_top_5_count,
tumour_1_or_2_not_top_5_prop,
tumour_1_and_2_not_top_5_count,
tumour_1_and_2_not_top_5_prop,
tumour_1_or_2_high_herit_count,
tumour_1_or_2_high_herit_prop,
tumour_1_and_2_high_herit_count,
tumour_1_and_2_high_herit_prop,
tumour_1_or_2_not_top_5_and_high_herit_count,
tumour_1_or_2_not_top_5_and_high_herit_prop,
tumour_1_and_2_not_top_5_and_high_herit_count,
tumour_1_and_2_not_top_5_and_high_herit_prop
))

colnames(combination_table_MPT) <- "MPT"

##Plot chord diagram
combination_count_for_circos <- combination_count_for_circos[which(combination_count_for_circos$Proportion > 0.25),]

png("circos_MPT.png",
    height = 1500,
    width = 1500)

library(circlize)
circos.par(gap.degree = 3)
chordDiagram(combination_count_for_circos, 
             transparency = 0.6, 
             annotationTrack = "grid", 
             preAllocateTracks = 1)

circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
circos.text(mean(xlim), ylim[1] + .1, 
              sector.name, facing = "clockwise", 
              niceFacing = TRUE, 
              cex = 2,
              adj = c(0, 0.5))
  
circos.axis(h = "top", 
              labels.cex = 0.6, 
              major.tick.percentage = 0.2, 
              sector.index = sector.name, 
              track.index = 2)
}, bg.border = NA)

dev.off()

circos.clear()

##Analysis of individual tumours
indv_tumours_all <- c(
as.character(vp_tumours$Single_word_1),
as.character(vp_tumours$Single_word_2),
as.character(vp_tumours$Single_word_3),
as.character(vp_tumours$Single_word_4),
as.character(vp_tumours$Single_word_5),
as.character(vp_tumours$Single_word_6)
)

indv_tumours_all <- indv_tumours_all[which(indv_tumours_all != "")]



for (i in 1:nrow(tumour_name_conversions)) {
  
  indv_tumours_all  <- gsub(tumour_name_conversions$original[i], tumour_name_conversions$new[i], indv_tumours_all)
  
}

indv_tumour_count <-rep(1, length(indv_tumours_all))
indv_tumour_table <- data.frame(indv_tumours_all, indv_tumour_count)

indv_tumour_count_table <- aggregate(indv_tumour_count~indv_tumours_all, data = indv_tumour_table,sum)

indv_tumour_count_table_ordered <- indv_tumour_count_table[order(indv_tumour_count_table$indv_tumour_count, decreasing = T),]

write.csv(indv_tumour_count_table_ordered, "indv_tumour_count_table_ordered_MPT.csv")

#############################################################
##MPT series only including tumours occurring before age 60##
#############################################################

vp_tumours <- read.csv("vp_tumours_under_60_only.csv")

##Extract combinations from table
comb_all <- data.frame()

comb_1 <- paste(vp_tumours$Single_word_1,vp_tumours$Single_word_2, sep = "-")
comb_all <- append(comb_1, comb_all)
comb_2 <- paste(vp_tumours$Single_word_1,vp_tumours$Single_word_3, sep = "-")
comb_all <- append(comb_2, comb_all)
comb_3 <- paste(vp_tumours$Single_word_1,vp_tumours$Single_word_4, sep = "-")
comb_all <- append(comb_3, comb_all)
comb_4 <- paste(vp_tumours$Single_word_1,vp_tumours$Single_word_5, sep = "-")
comb_all <- append(comb_4, comb_all)
comb_5 <- paste(vp_tumours$Single_word_1,vp_tumours$Single_word_6, sep = "-")
comb_all <- append(comb_5, comb_all)

comb_6 <- paste(vp_tumours$Single_word_2,vp_tumours$Single_word_3, sep = "-")
comb_all <- append(comb_6, comb_all)
comb_7 <- paste(vp_tumours$Single_word_2,vp_tumours$Single_word_4, sep = "-")
comb_all <- append(comb_7, comb_all)
comb_8 <- paste(vp_tumours$Single_word_2,vp_tumours$Single_word_5, sep = "-")
comb_all <- append(comb_8, comb_all)
comb_9 <- paste(vp_tumours$Single_word_2,vp_tumours$Single_word_6, sep = "-")
comb_all <- append(comb_9, comb_all)

comb_10 <- paste(vp_tumours$Single_word_3,vp_tumours$Single_word_4, sep = "-")
comb_all <- append(comb_10, comb_all)
comb_11 <- paste(vp_tumours$Single_word_3,vp_tumours$Single_word_5, sep = "-")
comb_all <- append(comb_11, comb_all)
comb_12 <- paste(vp_tumours$Single_word_3,vp_tumours$Single_word_6, sep = "-")
comb_all <- append(comb_12, comb_all)

comb_13 <- paste(vp_tumours$Single_word_4,vp_tumours$Single_word_5, sep = "-")
comb_all <- append(comb_13, comb_all)
comb_14 <- paste(vp_tumours$Single_word_4,vp_tumours$Single_word_6, sep = "-")
comb_all <- append(comb_14, comb_all)

comb_15 <- paste(vp_tumours$Single_word_5,vp_tumours$Single_word_6, sep = "-")
comb_all <- append(comb_15, comb_all)

comb_all <- unlist(comb_all)
comb_all_split <- strsplit(comb_all, "-")

comb_all_logical <- data.frame()

for (i in 1:length(comb_all_split)) {
  
  a <- isTRUE(length(comb_all_split[[i]]) > 1)
  comb_all_logical <- append(comb_all_logical,a)
  
}

comb_all_clean <- comb_all_split[which(comb_all_logical == TRUE)]

comb_all_table <- as.data.frame(t(data.frame(comb_all_clean, header = F)))
rownames(comb_all_table) <- seq(1:nrow(comb_all_table))
comb_all_table$V1 <- as.character(comb_all_table$V1)
comb_all_table$V2 <- as.character(comb_all_table$V2)

##Convert tumour names
tumour_name_conversions$original <- as.character(tumour_name_conversions$original)
tumour_name_conversions$new <- as.character(tumour_name_conversions$new)

for (i in 1:nrow(tumour_name_conversions)) {
  
  comb_all_table$V1  <- gsub(tumour_name_conversions$original[i], tumour_name_conversions$new[i], comb_all_table$V1)
  
}

for (i in 1:nrow(tumour_name_conversions)) {
  
  comb_all_table$V2  <- gsub(tumour_name_conversions$original[i], tumour_name_conversions$new[i], comb_all_table$V2)
  
}

##Create table of combinations
comb_all_table_no_dup <- comb_all_table[which(comb_all_table$V1 != comb_all_table$V2),]

comb_all_table_no_dup_sorted <- data.frame(t(apply(comb_all_table_no_dup, 1, sort)))
comb_all_table_no_dup_sorted$X1 <- as.character(comb_all_table_no_dup_sorted$X1)
comb_all_table_no_dup_sorted$X2 <- as.character(comb_all_table_no_dup_sorted$X2)

count <- rep(1, nrow(comb_all_table_no_dup_sorted))
combination <- paste(comb_all_table_no_dup_sorted$X1, comb_all_table_no_dup_sorted$X2, sep = "-")

comb_all_table_no_dup_sorted <- data.frame(cbind(comb_all_table_no_dup_sorted$X1, comb_all_table_no_dup_sorted$X2, combination, count))
comb_all_table_no_dup_sorted$count <- as.integer(comb_all_table_no_dup_sorted$count)

##Create combination count table
combination_count <- aggregate(count~combination, data = comb_all_table_no_dup_sorted, sum)
proportion <- (100/sum(combination_count$count)) * (combination_count$count)
combination_count <- data.frame(combination_count, proportion)

##Output table
combination_count <- combination_count[order(combination_count$count, decreasing = T),]
write.csv(combination_count, "combination_count_MPT_u60.csv")

##Check which combinations make up more than 1 percent total
combination_count_greater_one_percent_total <- combination_count[which(combination_count$proportion > 1),]

##Create chord diagram of combinations
tumors <- strsplit(as.character(combination_count$combination), "-")
tumor_1 <- sapply(tumors,function(x) x[1])
tumor_2 <- sapply(tumors,function(x) x[2])

combination_count_for_circos <- data.frame(tumor_1, tumor_2, combination_count$count, combination_count$proportion)
colnames(combination_count_for_circos) <- c("Tumor_1", "Tumor_2", "Count", "Proportion")


##Calculate characteristics of combinations before plotting chord diagram
top_5 <- unique(tumour_name_conversions$new[which(tumour_name_conversions$top_5 == "YES")])

tumour_1_or_2_not_top_5 <- !(combination_count_for_circos$Tumor_1 %in% top_5) | !(combination_count_for_circos$Tumor_2 %in% top_5)
tumour_1_and_2_not_top_5 <- !(combination_count_for_circos$Tumor_1 %in% top_5) & !(combination_count_for_circos$Tumor_2 %in% top_5)


tumour_1_or_2_not_top_5_count <- sum(combination_count_for_circos$Count[tumour_1_or_2_not_top_5])
tumour_1_or_2_not_top_5_prop <- sum(combination_count_for_circos$Count[tumour_1_or_2_not_top_5])/sum(combination_count_for_circos$Count) * 100

tumour_1_and_2_not_top_5_count <- sum(combination_count_for_circos$Count[tumour_1_and_2_not_top_5])
tumour_1_and_2_not_top_5_prop <- sum(combination_count_for_circos$Count[tumour_1_and_2_not_top_5])/sum(combination_count_for_circos$Count) * 100

high_herit <- unique(tumour_name_conversions$new[which(tumour_name_conversions$heritibility > 20)]) ## dummy number at present

tumour_1_or_2_high_herit <- !(combination_count_for_circos$Tumor_1 %in% high_herit) | !(combination_count_for_circos$Tumor_2 %in% high_herit)
tumour_1_and_2_high_herit <- !(combination_count_for_circos$Tumor_1 %in% high_herit) & !(combination_count_for_circos$Tumor_2 %in% high_herit)


tumour_1_or_2_high_herit_count <- sum(combination_count_for_circos$Count[tumour_1_or_2_high_herit])
tumour_1_or_2_high_herit_prop <- sum(combination_count_for_circos$Count[tumour_1_or_2_high_herit])/sum(combination_count_for_circos$Count) * 100

tumour_1_and_2_high_herit_count <- sum(combination_count_for_circos$Count[tumour_1_and_2_high_herit])
tumour_1_and_2_high_herit_prop <- sum(combination_count_for_circos$Count[tumour_1_and_2_high_herit])/sum(combination_count_for_circos$Count) * 100

tumour_1_or_2_not_top_5_and_high_herit_count <- sum(combination_count_for_circos$Count[tumour_1_or_2_not_top_5 == T & tumour_1_or_2_high_herit == T])
tumour_1_or_2_not_top_5_and_high_herit_prop <- sum(combination_count_for_circos$Count[tumour_1_or_2_not_top_5 == T & tumour_1_or_2_high_herit == T])/sum(combination_count_for_circos$Count) * 100

tumour_1_and_2_not_top_5_and_high_herit_count <- sum(combination_count_for_circos$Count[tumour_1_and_2_not_top_5 == T & tumour_1_and_2_high_herit == T])
tumour_1_and_2_not_top_5_and_high_herit_prop <- sum(combination_count_for_circos$Count[tumour_1_and_2_not_top_5 == T & tumour_1_and_2_high_herit == T])/sum(combination_count_for_circos$Count) * 100

##Output measures
number_discordant_combinations <- sum(combination_count_for_circos$Count)

combination_table_MPT_u60 <- as.data.frame(rbind(
  number_discordant_combinations,
  tumour_1_or_2_not_top_5_count,
  tumour_1_or_2_not_top_5_prop,
  tumour_1_and_2_not_top_5_count,
  tumour_1_and_2_not_top_5_prop,
  tumour_1_or_2_high_herit_count,
  tumour_1_or_2_high_herit_prop,
  tumour_1_and_2_high_herit_count,
  tumour_1_and_2_high_herit_prop,
  tumour_1_or_2_not_top_5_and_high_herit_count,
  tumour_1_or_2_not_top_5_and_high_herit_prop,
  tumour_1_and_2_not_top_5_and_high_herit_count,
  tumour_1_and_2_not_top_5_and_high_herit_prop
))

colnames(combination_table_MPT_u60) <- "MPT_u60"

##Plot chord diagram
combination_count_for_circos <- combination_count_for_circos[which(combination_count_for_circos$Proportion > 0.25),]

png("circos_MPT_u60.png",
    height = 1500,
    width = 1500)

library(circlize)
circos.par(gap.degree = 3)
chordDiagram(combination_count_for_circos, 
             transparency = 0.6, 
             annotationTrack = "grid", 
             preAllocateTracks = 1)

circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  circos.text(mean(xlim), ylim[1] + .1, 
              sector.name, facing = "clockwise", 
              niceFacing = TRUE, 
              cex = 2,
              adj = c(0, 0.5))
  
  circos.axis(h = "top", 
              labels.cex = 0.6, 
              major.tick.percentage = 0.2, 
              sector.index = sector.name, 
              track.index = 2)
}, bg.border = NA)

dev.off()

circos.clear()

##Analysis of individual tumours
indv_tumours_all <- c(
  as.character(vp_tumours$Single_word_1),
  as.character(vp_tumours$Single_word_2),
  as.character(vp_tumours$Single_word_3),
  as.character(vp_tumours$Single_word_4),
  as.character(vp_tumours$Single_word_5),
  as.character(vp_tumours$Single_word_6)
)

indv_tumours_all <- indv_tumours_all[which(indv_tumours_all != "")]

for (i in 1:nrow(tumour_name_conversions)) {
  
  indv_tumours_all  <- gsub(tumour_name_conversions$original[i], tumour_name_conversions$new[i], indv_tumours_all)
}

indv_tumour_count <-rep(1, length(indv_tumours_all))
indv_tumour_table <- data.frame(indv_tumours_all, indv_tumour_count)

indv_tumour_count_table <- aggregate(indv_tumour_count~indv_tumours_all, data = indv_tumour_table,sum)

indv_tumour_count_table_ordered <- indv_tumour_count_table[order(indv_tumour_count_table$indv_tumour_count, decreasing = T),]

write.csv(indv_tumour_count_table_ordered, "indv_tumour_count_table_ordered_MPT_u60.csv")

##############
##AVL SERIES##
##############

##Read in data table containing tumour information
vp_tumours <- read.csv("AVL_tumours.csv")

##Extract combinations from table
comb_all <- data.frame()

comb_1 <- paste(vp_tumours$Single_word_1,vp_tumours$Single_word_2, sep = "-")
comb_all <- append(comb_1, comb_all)
comb_2 <- paste(vp_tumours$Single_word_1,vp_tumours$Single_word_3, sep = "-")
comb_all <- append(comb_2, comb_all)
comb_3 <- paste(vp_tumours$Single_word_1,vp_tumours$Single_word_4, sep = "-")
comb_all <- append(comb_3, comb_all)
comb_4 <- paste(vp_tumours$Single_word_1,vp_tumours$Single_word_5, sep = "-")
comb_all <- append(comb_4, comb_all)
comb_5 <- paste(vp_tumours$Single_word_1,vp_tumours$Single_word_6, sep = "-")
comb_all <- append(comb_5, comb_all)

comb_6 <- paste(vp_tumours$Single_word_2,vp_tumours$Single_word_3, sep = "-")
comb_all <- append(comb_6, comb_all)
comb_7 <- paste(vp_tumours$Single_word_2,vp_tumours$Single_word_4, sep = "-")
comb_all <- append(comb_7, comb_all)
comb_8 <- paste(vp_tumours$Single_word_2,vp_tumours$Single_word_5, sep = "-")
comb_all <- append(comb_8, comb_all)
comb_9 <- paste(vp_tumours$Single_word_2,vp_tumours$Single_word_6, sep = "-")
comb_all <- append(comb_9, comb_all)

comb_10 <- paste(vp_tumours$Single_word_3,vp_tumours$Single_word_4, sep = "-")
comb_all <- append(comb_10, comb_all)
comb_11 <- paste(vp_tumours$Single_word_3,vp_tumours$Single_word_5, sep = "-")
comb_all <- append(comb_11, comb_all)
comb_12 <- paste(vp_tumours$Single_word_3,vp_tumours$Single_word_6, sep = "-")
comb_all <- append(comb_12, comb_all)

comb_13 <- paste(vp_tumours$Single_word_4,vp_tumours$Single_word_5, sep = "-")
comb_all <- append(comb_13, comb_all)
comb_14 <- paste(vp_tumours$Single_word_4,vp_tumours$Single_word_6, sep = "-")
comb_all <- append(comb_14, comb_all)

comb_15 <- paste(vp_tumours$Single_word_5,vp_tumours$Single_word_6, sep = "-")
comb_all <- append(comb_15, comb_all)

comb_all <- unlist(comb_all)
comb_all_split <- strsplit(comb_all, "-")

comb_all_logical <- data.frame()

for (i in 1:length(comb_all_split)) {
  
  a <- isTRUE(length(comb_all_split[[i]]) > 1)
  comb_all_logical <- append(comb_all_logical,a)
  
}

comb_all_clean <- comb_all_split[which(comb_all_logical == TRUE)]

comb_all_table <- as.data.frame(t(data.frame(comb_all_clean, header = F)))
rownames(comb_all_table) <- seq(1:nrow(comb_all_table))
comb_all_table$V1 <- as.character(comb_all_table$V1)
comb_all_table$V2 <- as.character(comb_all_table$V2)

##Convert tumour names
tumour_name_conversions$original <- as.character(tumour_name_conversions$original)
tumour_name_conversions$new <- as.character(tumour_name_conversions$new)

for (i in 1:nrow(tumour_name_conversions)) {
  
  comb_all_table$V1  <- gsub(tumour_name_conversions$original[i], tumour_name_conversions$new[i], comb_all_table$V1)
  
}

for (i in 1:nrow(tumour_name_conversions)) {
  
  comb_all_table$V2  <- gsub(tumour_name_conversions$original[i], tumour_name_conversions$new[i], comb_all_table$V2)
  
}

##Create table of combinations
comb_all_table_no_dup <- comb_all_table[which(comb_all_table$V1 != comb_all_table$V2),]

comb_all_table_no_dup_sorted <- data.frame(t(apply(comb_all_table_no_dup, 1, sort)))
comb_all_table_no_dup_sorted$X1 <- as.character(comb_all_table_no_dup_sorted$X1)
comb_all_table_no_dup_sorted$X2 <- as.character(comb_all_table_no_dup_sorted$X2)

count <- rep(1, nrow(comb_all_table_no_dup_sorted))
combination <- paste(comb_all_table_no_dup_sorted$X1, comb_all_table_no_dup_sorted$X2, sep = "-")

comb_all_table_no_dup_sorted <- data.frame(cbind(comb_all_table_no_dup_sorted$X1, comb_all_table_no_dup_sorted$X2, combination, count))
comb_all_table_no_dup_sorted$count <- as.integer(comb_all_table_no_dup_sorted$count)

##Create combination count table
combination_count <- aggregate(count~combination, data = comb_all_table_no_dup_sorted, sum)
proportion <- (100/sum(combination_count$count)) * (combination_count$count)
combination_count <- data.frame(combination_count, proportion)

combination_count <- combination_count[order(combination_count$count, decreasing = T),]

##Output table
write.csv(combination_count, "combination_count_AVL.csv")

##Check which combinations make up more than 1 percent total
combination_count_greater_one_percent_total <- combination_count[which(combination_count$proportion > 1),]

##Create chord diagram of combinations
tumors <- strsplit(as.character(combination_count$combination), "-")
tumor_1 <- sapply(tumors,function(x) x[1])
tumor_2 <- sapply(tumors,function(x) x[2])

combination_count_for_circos <- data.frame(tumor_1, tumor_2, combination_count$count, combination_count$proportion)
colnames(combination_count_for_circos) <- c("Tumor_1", "Tumor_2", "Count", "Proportion")


##Calculate characteristics of combinations before plotting chord diagram
top_5 <- unique(tumour_name_conversions$new[which(tumour_name_conversions$top_5 == "YES")])

tumour_1_or_2_not_top_5 <- !(combination_count_for_circos$Tumor_1 %in% top_5) | !(combination_count_for_circos$Tumor_2 %in% top_5)
tumour_1_and_2_not_top_5 <- !(combination_count_for_circos$Tumor_1 %in% top_5) & !(combination_count_for_circos$Tumor_2 %in% top_5)

tumour_1_or_2_not_top_5_count <- sum(combination_count_for_circos$Count[tumour_1_or_2_not_top_5])
tumour_1_or_2_not_top_5_prop <- sum(combination_count_for_circos$Count[tumour_1_or_2_not_top_5])/sum(combination_count_for_circos$Count) * 100

tumour_1_and_2_not_top_5_count <- sum(combination_count_for_circos$Count[tumour_1_and_2_not_top_5])
tumour_1_and_2_not_top_5_prop <- sum(combination_count_for_circos$Count[tumour_1_and_2_not_top_5])/sum(combination_count_for_circos$Count) * 100

high_herit <- unique(tumour_name_conversions$new[which(tumour_name_conversions$heritibility > 20)]) 

tumour_1_or_2_high_herit <- !(combination_count_for_circos$Tumor_1 %in% high_herit) | !(combination_count_for_circos$Tumor_2 %in% high_herit)
tumour_1_and_2_high_herit <- !(combination_count_for_circos$Tumor_1 %in% high_herit) & !(combination_count_for_circos$Tumor_2 %in% high_herit)

tumour_1_or_2_high_herit_count <- sum(combination_count_for_circos$Count[tumour_1_or_2_high_herit])
tumour_1_or_2_high_herit_prop <- sum(combination_count_for_circos$Count[tumour_1_or_2_high_herit])/sum(combination_count_for_circos$Count) * 100

tumour_1_and_2_high_herit_count <- sum(combination_count_for_circos$Count[tumour_1_and_2_high_herit])
tumour_1_and_2_high_herit_prop <- sum(combination_count_for_circos$Count[tumour_1_and_2_high_herit])/sum(combination_count_for_circos$Count) * 100

tumour_1_or_2_not_top_5_and_high_herit_count <- sum(combination_count_for_circos$Count[tumour_1_or_2_not_top_5 == T & tumour_1_or_2_high_herit == T])
tumour_1_or_2_not_top_5_and_high_herit_prop <- sum(combination_count_for_circos$Count[tumour_1_or_2_not_top_5 == T & tumour_1_or_2_high_herit == T])/sum(combination_count_for_circos$Count) * 100

tumour_1_and_2_not_top_5_and_high_herit_count <- sum(combination_count_for_circos$Count[tumour_1_and_2_not_top_5 == T & tumour_1_and_2_high_herit == T])
tumour_1_and_2_not_top_5_and_high_herit_prop <- sum(combination_count_for_circos$Count[tumour_1_and_2_not_top_5 == T & tumour_1_and_2_high_herit == T])/sum(combination_count_for_circos$Count) * 100

##Output measures
number_discordant_combinations <- sum(combination_count_for_circos$Count)

combination_table_AVL <- as.data.frame(rbind(
  number_discordant_combinations,
  tumour_1_or_2_not_top_5_count,
  tumour_1_or_2_not_top_5_prop,
  tumour_1_and_2_not_top_5_count,
  tumour_1_and_2_not_top_5_prop,
  tumour_1_or_2_high_herit_count,
  tumour_1_or_2_high_herit_prop,
  tumour_1_and_2_high_herit_count,
  tumour_1_and_2_high_herit_prop,
  tumour_1_or_2_not_top_5_and_high_herit_count,
  tumour_1_or_2_not_top_5_and_high_herit_prop,
  tumour_1_and_2_not_top_5_and_high_herit_count,
  tumour_1_and_2_not_top_5_and_high_herit_prop
))

colnames(combination_table_AVL) <- "AVL"

##Calculate characteristics of combinations before chord diagram
combination_count_for_circos <- combination_count_for_circos[which(combination_count_for_circos$Proportion > 0.25),]

png("circos_AVL.png",
    height = 1600,
    width = 1600)

library(circlize)
circos.par(gap.degree = 3)
chordDiagram(combination_count_for_circos, 
             transparency = 0.6, 
             annotationTrack = "grid", 
             preAllocateTracks = 1)

circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  circos.text(mean(xlim), ylim[1] + .1, 
              sector.name, facing = "clockwise", 
              niceFacing = TRUE, 
              cex = 2,
              adj = c(0, 0.5))
  
  circos.axis(h = "top", 
              labels.cex = 0.6, 
              major.tick.percentage = 0.2, 
              sector.index = sector.name, 
              track.index = 2)
}, bg.border = NA)

dev.off()

circos.clear()

##Analysis of individual tumours
indv_tumours_all <- c(
  as.character(vp_tumours$Single_word_1),
  as.character(vp_tumours$Single_word_2),
  as.character(vp_tumours$Single_word_3),
  as.character(vp_tumours$Single_word_4),
  as.character(vp_tumours$Single_word_5),
  as.character(vp_tumours$Single_word_6)
)

indv_tumours_all <- indv_tumours_all[which(indv_tumours_all != "")]

for (i in 1:nrow(tumour_name_conversions)) {
  
  indv_tumours_all  <- gsub(tumour_name_conversions$original[i], tumour_name_conversions$new[i], indv_tumours_all)
  
}

indv_tumour_count <-rep(1, length(indv_tumours_all))
indv_tumour_table <- data.frame(indv_tumours_all, indv_tumour_count)

indv_tumour_count_table <- aggregate(indv_tumour_count~indv_tumours_all, data = indv_tumour_table,sum)

indv_tumour_count_table_ordered <- indv_tumour_count_table[order(indv_tumour_count_table$indv_tumour_count, decreasing = T),]

write.csv(indv_tumour_count_table_ordered, "indv_tumour_count_table_ordered_AVL.csv")

#########################
##Dutch registry series##
#########################

##Read in data table containing tumour information
vp_tumours <- read.csv("dutch_reg_tumours.csv")

##Extract combinations from table
comb_all <- data.frame()

comb_1 <- paste(vp_tumours$Single_word_1,vp_tumours$Single_word_2, sep = "-")
comb_all <- append(comb_1, comb_all)
comb_2 <- paste(vp_tumours$Single_word_1,vp_tumours$Single_word_3, sep = "-")
comb_all <- append(comb_2, comb_all)
comb_3 <- paste(vp_tumours$Single_word_1,vp_tumours$Single_word_4, sep = "-")
comb_all <- append(comb_3, comb_all)
comb_4 <- paste(vp_tumours$Single_word_1,vp_tumours$Single_word_5, sep = "-")
comb_all <- append(comb_4, comb_all)
comb_5 <- paste(vp_tumours$Single_word_1,vp_tumours$Single_word_6, sep = "-")
comb_all <- append(comb_5, comb_all)

comb_6 <- paste(vp_tumours$Single_word_2,vp_tumours$Single_word_3, sep = "-")
comb_all <- append(comb_6, comb_all)
comb_7 <- paste(vp_tumours$Single_word_2,vp_tumours$Single_word_4, sep = "-")
comb_all <- append(comb_7, comb_all)
comb_8 <- paste(vp_tumours$Single_word_2,vp_tumours$Single_word_5, sep = "-")
comb_all <- append(comb_8, comb_all)
comb_9 <- paste(vp_tumours$Single_word_2,vp_tumours$Single_word_6, sep = "-")
comb_all <- append(comb_9, comb_all)

comb_10 <- paste(vp_tumours$Single_word_3,vp_tumours$Single_word_4, sep = "-")
comb_all <- append(comb_10, comb_all)
comb_11 <- paste(vp_tumours$Single_word_3,vp_tumours$Single_word_5, sep = "-")
comb_all <- append(comb_11, comb_all)
comb_12 <- paste(vp_tumours$Single_word_3,vp_tumours$Single_word_6, sep = "-")
comb_all <- append(comb_12, comb_all)

comb_13 <- paste(vp_tumours$Single_word_4,vp_tumours$Single_word_5, sep = "-")
comb_all <- append(comb_13, comb_all)
comb_14 <- paste(vp_tumours$Single_word_4,vp_tumours$Single_word_6, sep = "-")
comb_all <- append(comb_14, comb_all)

comb_15 <- paste(vp_tumours$Single_word_5,vp_tumours$Single_word_6, sep = "-")
comb_all <- append(comb_15, comb_all)

comb_all <- unlist(comb_all)
comb_all_split <- strsplit(comb_all, "-")

comb_all_logical <- data.frame()

for (i in 1:length(comb_all_split)) {
  
  a <- isTRUE(length(comb_all_split[[i]]) > 1)
  comb_all_logical <- append(comb_all_logical,a)
  
}

comb_all_clean <- comb_all_split[which(comb_all_logical == TRUE)]

comb_all_table <- as.data.frame(t(data.frame(comb_all_clean, header = F)))
rownames(comb_all_table) <- seq(1:nrow(comb_all_table))
comb_all_table$V1 <- as.character(comb_all_table$V1)
comb_all_table$V2 <- as.character(comb_all_table$V2)

##Convert tumour names
tumour_name_conversions$original <- as.character(tumour_name_conversions$original)
tumour_name_conversions$new <- as.character(tumour_name_conversions$new)

for (i in 1:nrow(tumour_name_conversions)) {
  
  comb_all_table$V1  <- gsub(tumour_name_conversions$original[i], tumour_name_conversions$new[i], comb_all_table$V1)
  
}

for (i in 1:nrow(tumour_name_conversions)) {
  
  comb_all_table$V2  <- gsub(tumour_name_conversions$original[i], tumour_name_conversions$new[i], comb_all_table$V2)
  
}

##Create table of combinations
comb_all_table_no_dup <- comb_all_table[which(comb_all_table$V1 != comb_all_table$V2),]

comb_all_table_no_dup_sorted <- data.frame(t(apply(comb_all_table_no_dup, 1, sort)))
comb_all_table_no_dup_sorted$X1 <- as.character(comb_all_table_no_dup_sorted$X1)
comb_all_table_no_dup_sorted$X2 <- as.character(comb_all_table_no_dup_sorted$X2)

count <- rep(1, nrow(comb_all_table_no_dup_sorted))
combination <- paste(comb_all_table_no_dup_sorted$X1, comb_all_table_no_dup_sorted$X2, sep = "-")

comb_all_table_no_dup_sorted <- data.frame(cbind(comb_all_table_no_dup_sorted$X1, comb_all_table_no_dup_sorted$X2, combination, count))
comb_all_table_no_dup_sorted$count <- as.integer(comb_all_table_no_dup_sorted$count)

##Create combination count table
combination_count <- aggregate(count~combination, data = comb_all_table_no_dup_sorted, sum)
proportion <- (100/sum(combination_count$count)) * (combination_count$count)
combination_count <- data.frame(combination_count, proportion)

combination_count <- combination_count[order(combination_count$count, decreasing = T),]

##Output table
write.csv(combination_count, "combination_count_dutch_reg.csv")

##Check which combinations make up more than 1 percent total
combination_count_greater_one_percent_total <- combination_count[which(combination_count$proportion > 1),]

##Create chord diagram of combinations
tumors <- strsplit(as.character(combination_count$combination), "-")
tumor_1 <- sapply(tumors,function(x) x[1])
tumor_2 <- sapply(tumors,function(x) x[2])

combination_count_for_circos <- data.frame(tumor_1, tumor_2, combination_count$count, combination_count$proportion)
colnames(combination_count_for_circos) <- c("Tumor_1", "Tumor_2", "Count", "Proportion")

##Calculate characteristics of combinations before chord diagram
top_5 <- unique(tumour_name_conversions$new[which(tumour_name_conversions$top_5 == "YES")])

tumour_1_or_2_not_top_5 <- !(combination_count_for_circos$Tumor_1 %in% top_5) | !(combination_count_for_circos$Tumor_2 %in% top_5)
tumour_1_and_2_not_top_5 <- !(combination_count_for_circos$Tumor_1 %in% top_5) & !(combination_count_for_circos$Tumor_2 %in% top_5)

tumour_1_or_2_not_top_5_count <- sum(combination_count_for_circos$Count[tumour_1_or_2_not_top_5])
tumour_1_or_2_not_top_5_prop <- sum(combination_count_for_circos$Count[tumour_1_or_2_not_top_5])/sum(combination_count_for_circos$Count) * 100

tumour_1_and_2_not_top_5_count <- sum(combination_count_for_circos$Count[tumour_1_and_2_not_top_5])
tumour_1_and_2_not_top_5_prop <- sum(combination_count_for_circos$Count[tumour_1_and_2_not_top_5])/sum(combination_count_for_circos$Count) * 100

high_herit <- unique(tumour_name_conversions$new[which(tumour_name_conversions$heritibility > 20)]) ## dummy number at present

tumour_1_or_2_high_herit <- !(combination_count_for_circos$Tumor_1 %in% high_herit) | !(combination_count_for_circos$Tumor_2 %in% high_herit)
tumour_1_and_2_high_herit <- !(combination_count_for_circos$Tumor_1 %in% high_herit) & !(combination_count_for_circos$Tumor_2 %in% high_herit)

tumour_1_or_2_high_herit_count <- sum(combination_count_for_circos$Count[tumour_1_or_2_high_herit])
tumour_1_or_2_high_herit_prop <- sum(combination_count_for_circos$Count[tumour_1_or_2_high_herit])/sum(combination_count_for_circos$Count) * 100

tumour_1_and_2_high_herit_count <- sum(combination_count_for_circos$Count[tumour_1_and_2_high_herit])
tumour_1_and_2_high_herit_prop <- sum(combination_count_for_circos$Count[tumour_1_and_2_high_herit])/sum(combination_count_for_circos$Count) * 100

tumour_1_or_2_not_top_5_and_high_herit_count <- sum(combination_count_for_circos$Count[tumour_1_or_2_not_top_5 == T & tumour_1_or_2_high_herit == T])
tumour_1_or_2_not_top_5_and_high_herit_prop <- sum(combination_count_for_circos$Count[tumour_1_or_2_not_top_5 == T & tumour_1_or_2_high_herit == T])/sum(combination_count_for_circos$Count) * 100

tumour_1_and_2_not_top_5_and_high_herit_count <- sum(combination_count_for_circos$Count[tumour_1_and_2_not_top_5 == T & tumour_1_and_2_high_herit == T])
tumour_1_and_2_not_top_5_and_high_herit_prop <- sum(combination_count_for_circos$Count[tumour_1_and_2_not_top_5 == T & tumour_1_and_2_high_herit == T])/sum(combination_count_for_circos$Count) * 100

##Output measures
number_discordant_combinations <- sum(combination_count_for_circos$Count)

combination_table_dutch_reg <- as.data.frame(rbind(
  number_discordant_combinations,
  tumour_1_or_2_not_top_5_count,
  tumour_1_or_2_not_top_5_prop,
  tumour_1_and_2_not_top_5_count,
  tumour_1_and_2_not_top_5_prop,
  tumour_1_or_2_high_herit_count,
  tumour_1_or_2_high_herit_prop,
  tumour_1_and_2_high_herit_count,
  tumour_1_and_2_high_herit_prop,
  tumour_1_or_2_not_top_5_and_high_herit_count,
  tumour_1_or_2_not_top_5_and_high_herit_prop,
  tumour_1_and_2_not_top_5_and_high_herit_count,
  tumour_1_and_2_not_top_5_and_high_herit_prop
))

colnames(combination_table_dutch_reg) <- "dutch_reg"

##Plot chord diagram
combination_count_for_circos <- combination_count_for_circos[which(combination_count_for_circos$Proportion > 0.25),]

png("circos_dutch_reg.png",
    height = 1500,
    width = 1500)
    
library(circlize)
circos.par(gap.degree = 3)
chordDiagram(combination_count_for_circos, 
             transparency = 0.6, 
             annotationTrack = "grid", 
             preAllocateTracks = 1)

circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  circos.text(mean(xlim), ylim[1] + .1, 
              sector.name, facing = "clockwise", 
              niceFacing = TRUE, 
              cex = 2,
              adj = c(0, 0.5))
  
  circos.axis(h = "top", 
              labels.cex = 0.6, 
              major.tick.percentage = 0.2, 
              sector.index = sector.name, 
              track.index = 2)
}, bg.border = NA)

dev.off()

circos.clear()

##Analysis of individual tumours
indv_tumours_all <- c(
  as.character(vp_tumours$Single_word_1),
  as.character(vp_tumours$Single_word_2),
  as.character(vp_tumours$Single_word_3),
  as.character(vp_tumours$Single_word_4),
  as.character(vp_tumours$Single_word_5),
  as.character(vp_tumours$Single_word_6)
)

indv_tumours_all <- indv_tumours_all[which(indv_tumours_all != "")]

for (i in 1:nrow(tumour_name_conversions)) {
  indv_tumours_all  <- gsub(tumour_name_conversions$original[i], tumour_name_conversions$new[i], indv_tumours_all)
}

indv_tumour_count <-rep(1, length(indv_tumours_all))
indv_tumour_table <- data.frame(indv_tumours_all, indv_tumour_count)

indv_tumour_count_table <- aggregate(indv_tumour_count~indv_tumours_all, data = indv_tumour_table,sum)

indv_tumour_count_table_ordered <- indv_tumour_count_table[order(indv_tumour_count_table$indv_tumour_count, decreasing = T),]

write.csv(indv_tumour_count_table_ordered, "indv_tumour_count_table_ordered_dutch_reg.csv")

###############################
##East Anglia Registry series##
###############################

##Read in data table containing tumour information
vp_tumours <- read.csv("EA_tumours.csv")

##Extract combinations from table
comb_all <- data.frame()

comb_1 <- paste(vp_tumours$Single_word_1,vp_tumours$Single_word_2, sep = "-")
comb_all <- append(comb_1, comb_all)
comb_2 <- paste(vp_tumours$Single_word_1,vp_tumours$Single_word_3, sep = "-")
comb_all <- append(comb_2, comb_all)
comb_3 <- paste(vp_tumours$Single_word_1,vp_tumours$Single_word_4, sep = "-")
comb_all <- append(comb_3, comb_all)
comb_4 <- paste(vp_tumours$Single_word_1,vp_tumours$Single_word_5, sep = "-")
comb_all <- append(comb_4, comb_all)
comb_5 <- paste(vp_tumours$Single_word_1,vp_tumours$Single_word_6, sep = "-")
comb_all <- append(comb_5, comb_all)

comb_6 <- paste(vp_tumours$Single_word_2,vp_tumours$Single_word_3, sep = "-")
comb_all <- append(comb_6, comb_all)
comb_7 <- paste(vp_tumours$Single_word_2,vp_tumours$Single_word_4, sep = "-")
comb_all <- append(comb_7, comb_all)
comb_8 <- paste(vp_tumours$Single_word_2,vp_tumours$Single_word_5, sep = "-")
comb_all <- append(comb_8, comb_all)
comb_9 <- paste(vp_tumours$Single_word_2,vp_tumours$Single_word_6, sep = "-")
comb_all <- append(comb_9, comb_all)

comb_10 <- paste(vp_tumours$Single_word_3,vp_tumours$Single_word_4, sep = "-")
comb_all <- append(comb_10, comb_all)
comb_11 <- paste(vp_tumours$Single_word_3,vp_tumours$Single_word_5, sep = "-")
comb_all <- append(comb_11, comb_all)
comb_12 <- paste(vp_tumours$Single_word_3,vp_tumours$Single_word_6, sep = "-")
comb_all <- append(comb_12, comb_all)

comb_13 <- paste(vp_tumours$Single_word_4,vp_tumours$Single_word_5, sep = "-")
comb_all <- append(comb_13, comb_all)
comb_14 <- paste(vp_tumours$Single_word_4,vp_tumours$Single_word_6, sep = "-")
comb_all <- append(comb_14, comb_all)

comb_15 <- paste(vp_tumours$Single_word_5,vp_tumours$Single_word_6, sep = "-")
comb_all <- append(comb_15, comb_all)

comb_all <- unlist(comb_all)
comb_all_split <- strsplit(comb_all, "-")

comb_all_logical <- data.frame()

for (i in 1:length(comb_all_split)) {
  
  a <- isTRUE(length(comb_all_split[[i]]) > 1)
  comb_all_logical <- append(comb_all_logical,a)
  
}

comb_all_clean <- comb_all_split[which(comb_all_logical == TRUE)]

comb_all_table <- as.data.frame(t(data.frame(comb_all_clean, header = F)))
rownames(comb_all_table) <- seq(1:nrow(comb_all_table))
comb_all_table$V1 <- as.character(comb_all_table$V1)
comb_all_table$V2 <- as.character(comb_all_table$V2)

##Convert tumour names
tumour_name_conversions$original <- as.character(tumour_name_conversions$original)
tumour_name_conversions$new <- as.character(tumour_name_conversions$new)

for (i in 1:nrow(tumour_name_conversions)) {
  comb_all_table$V1  <- gsub(tumour_name_conversions$original[i], tumour_name_conversions$new[i], comb_all_table$V1)
}

for (i in 1:nrow(tumour_name_conversions)) {
  comb_all_table$V2  <- gsub(tumour_name_conversions$original[i], tumour_name_conversions$new[i], comb_all_table$V2)
}

##Create table of combinations
comb_all_table_no_dup <- comb_all_table[which(comb_all_table$V1 != comb_all_table$V2),]

comb_all_table_no_dup_sorted <- data.frame(t(apply(comb_all_table_no_dup, 1, sort)))
comb_all_table_no_dup_sorted$X1 <- as.character(comb_all_table_no_dup_sorted$X1)
comb_all_table_no_dup_sorted$X2 <- as.character(comb_all_table_no_dup_sorted$X2)

count <- rep(1, nrow(comb_all_table_no_dup_sorted))
combination <- paste(comb_all_table_no_dup_sorted$X1, comb_all_table_no_dup_sorted$X2, sep = "-")

comb_all_table_no_dup_sorted <- data.frame(cbind(comb_all_table_no_dup_sorted$X1, comb_all_table_no_dup_sorted$X2, combination, count))
comb_all_table_no_dup_sorted$count <- as.integer(comb_all_table_no_dup_sorted$count)

##Create combination count table
combination_count <- aggregate(count~combination, data = comb_all_table_no_dup_sorted, sum)
proportion <- (100/sum(combination_count$count)) * (combination_count$count)
combination_count <- data.frame(combination_count, proportion)

combination_count <- combination_count[order(combination_count$count, decreasing = T),]

##Output table
write.csv(combination_count, "combination_count_EA.csv")

##Check which combinations make up more than 1 percent total
combination_count_greater_one_percent_total <- combination_count[which(combination_count$proportion > 1),]

##Create chord diagram of combinations
tumors <- strsplit(as.character(combination_count$combination), "-")
tumor_1 <- sapply(tumors,function(x) x[1])
tumor_2 <- sapply(tumors,function(x) x[2])

combination_count_for_circos <- data.frame(tumor_1, tumor_2, combination_count$count, combination_count$proportion)
colnames(combination_count_for_circos) <- c("Tumor_1", "Tumor_2", "Count", "Proportion")

##Calculate characteristics of combinations before chord diagram
top_5 <- unique(tumour_name_conversions$new[which(tumour_name_conversions$top_5 == "YES")])

tumour_1_or_2_not_top_5 <- !(combination_count_for_circos$Tumor_1 %in% top_5) | !(combination_count_for_circos$Tumor_2 %in% top_5)
tumour_1_and_2_not_top_5 <- !(combination_count_for_circos$Tumor_1 %in% top_5) & !(combination_count_for_circos$Tumor_2 %in% top_5)

tumour_1_or_2_not_top_5_count <- sum(combination_count_for_circos$Count[tumour_1_or_2_not_top_5])
tumour_1_or_2_not_top_5_prop <- sum(combination_count_for_circos$Count[tumour_1_or_2_not_top_5])/sum(combination_count_for_circos$Count) * 100

tumour_1_and_2_not_top_5_count <- sum(combination_count_for_circos$Count[tumour_1_and_2_not_top_5])
tumour_1_and_2_not_top_5_prop <- sum(combination_count_for_circos$Count[tumour_1_and_2_not_top_5])/sum(combination_count_for_circos$Count) * 100

high_herit <- unique(tumour_name_conversions$new[which(tumour_name_conversions$heritibility > 20)]) ## dummy number at present

tumour_1_or_2_high_herit <- !(combination_count_for_circos$Tumor_1 %in% high_herit) | !(combination_count_for_circos$Tumor_2 %in% high_herit)
tumour_1_and_2_high_herit <- !(combination_count_for_circos$Tumor_1 %in% high_herit) & !(combination_count_for_circos$Tumor_2 %in% high_herit)


tumour_1_or_2_high_herit_count <- sum(combination_count_for_circos$Count[tumour_1_or_2_high_herit])
tumour_1_or_2_high_herit_prop <- sum(combination_count_for_circos$Count[tumour_1_or_2_high_herit])/sum(combination_count_for_circos$Count) * 100

tumour_1_and_2_high_herit_count <- sum(combination_count_for_circos$Count[tumour_1_and_2_high_herit])
tumour_1_and_2_high_herit_prop <- sum(combination_count_for_circos$Count[tumour_1_and_2_high_herit])/sum(combination_count_for_circos$Count) * 100

tumour_1_or_2_not_top_5_and_high_herit_count <- sum(combination_count_for_circos$Count[tumour_1_or_2_not_top_5 == T & tumour_1_or_2_high_herit == T])
tumour_1_or_2_not_top_5_and_high_herit_prop <- sum(combination_count_for_circos$Count[tumour_1_or_2_not_top_5 == T & tumour_1_or_2_high_herit == T])/sum(combination_count_for_circos$Count) * 100

tumour_1_and_2_not_top_5_and_high_herit_count <- sum(combination_count_for_circos$Count[tumour_1_and_2_not_top_5 == T & tumour_1_and_2_high_herit == T])
tumour_1_and_2_not_top_5_and_high_herit_prop <- sum(combination_count_for_circos$Count[tumour_1_and_2_not_top_5 == T & tumour_1_and_2_high_herit == T])/sum(combination_count_for_circos$Count) * 100

##Output measures
number_discordant_combinations <- sum(combination_count_for_circos$Count)

combination_table_EA <- as.data.frame(rbind(
  number_discordant_combinations,
  tumour_1_or_2_not_top_5_count,
  tumour_1_or_2_not_top_5_prop,
  tumour_1_and_2_not_top_5_count,
  tumour_1_and_2_not_top_5_prop,
  tumour_1_or_2_high_herit_count,
  tumour_1_or_2_high_herit_prop,
  tumour_1_and_2_high_herit_count,
  tumour_1_and_2_high_herit_prop,
  tumour_1_or_2_not_top_5_and_high_herit_count,
  tumour_1_or_2_not_top_5_and_high_herit_prop,
  tumour_1_and_2_not_top_5_and_high_herit_count,
  tumour_1_and_2_not_top_5_and_high_herit_prop
))

colnames(combination_table_EA) <- "EA"

##Plot chord diagram
combination_count_for_circos <- combination_count_for_circos[which(combination_count_for_circos$Proportion > 0.25),]

png("circos_EA.png",
    height = 1500,
    width = 1500)

library(circlize)
circos.par(gap.degree = 3)
chordDiagram(combination_count_for_circos, 
             transparency = 0.6, 
             annotationTrack = "grid", 
             preAllocateTracks = 1)

circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  circos.text(mean(xlim), ylim[1] + .1, 
              sector.name, facing = "clockwise", 
              niceFacing = TRUE, 
              cex = 2,
              adj = c(0, 0.5))
  
  circos.axis(h = "top", 
              labels.cex = 0.6, 
              major.tick.percentage = 0.2, 
              sector.index = sector.name, 
              track.index = 2)
}, bg.border = NA)


dev.off()

circos.clear()

##Analysis of individual tumours
indv_tumours_all <- c(
  as.character(vp_tumours$Single_word_1),
  as.character(vp_tumours$Single_word_2),
  as.character(vp_tumours$Single_word_3),
  as.character(vp_tumours$Single_word_4),
  as.character(vp_tumours$Single_word_5),
  as.character(vp_tumours$Single_word_6)
)

indv_tumours_all <- indv_tumours_all[which(indv_tumours_all != "")]

for (i in 1:nrow(tumour_name_conversions)) {
  indv_tumours_all  <- gsub(tumour_name_conversions$original[i], tumour_name_conversions$new[i], indv_tumours_all)
}

indv_tumour_count <-rep(1, length(indv_tumours_all))
indv_tumour_table <- data.frame(indv_tumours_all, indv_tumour_count)

indv_tumour_count_table <- aggregate(indv_tumour_count~indv_tumours_all, data = indv_tumour_table,sum)

indv_tumour_count_table_ordered <- indv_tumour_count_table[order(indv_tumour_count_table$indv_tumour_count, decreasing = T),]

write.csv(indv_tumour_count_table_ordered, "indv_tumour_count_table_ordered_EA.csv")

##Create table of combination counts and charasteristics in different series
combination_table_combined <- cbind(
  combination_table_MPT,
  combination_table_MPT_u60,
  combination_table_AVL,
  combination_table_dutch_reg,
  combination_table_EA
)

##Output table
write.csv(combination_table_combined, "combination_table_combined.csv")

##Read in table (produced in excel) of combination counts in MPT series (tumours under 60 years only) and East Anglia Registry Series. Includes, for statitical tests, a count of combinations that were NOT of the type described in each line of the table 
MPT_u60_EA_comparison <- read.csv("MPT_u60_EA_comparison.csv")

##Comparison of combination frequency between MPT series (tumours under 60 years only) and East Anglia Registry Series
pvalue <- as.vector(c(), mode = "any")

for(i in 1:nrow(MPT_u60_EA_comparison)){
  
  fishtable = matrix(c(MPT_u60_EA_comparison$MPT_u60_count[i],
                       MPT_u60_EA_comparison$MPT_u60_combos_not[i],
                       MPT_u60_EA_comparison$EA_count[i],
                       MPT_u60_EA_comparison$EA_combos_not[i]), 
                     nrow = 2)
  
  fisherout <- fisher.test(fishtable)
  
  pvalue <- append(pvalue, fisherout[[1]])
  
}

pvalue <- round(pvalue, digits=5)

##Add pvalues to table and output results
MPT_u60_EA_comparison <- cbind(MPT_u60_EA_comparison, pvalue)
write.csv(MPT_u60_EA_comparison, "MPT_u60_EA_comparison.csv")


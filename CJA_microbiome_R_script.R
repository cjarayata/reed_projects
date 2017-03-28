# CJA script - last update 1/24/2017

pacman::p_load(data.table, broom, heplots)


### CAUTION: Negative control analysis is currently underway

setwd("//monell/storage/reed/Labshare/CJ Arayata/Cohen/Garth Microbiome/2017 - January/")

setwd("C:/Users/carayata/Dropbox/2017_T2R_nasal_microbiome/2000. Microbiome data/CJA")


# This reads in several files to make these: ----

# patients        = patient info (incomplete/with extra things)
# sites           = info on anatomical sampling (some amibguous codes)
# samples         = info on all samples run --> CONTAINS CRUCIAL HS NUMBERS
# tax             = taxonomic assignment for all OTU in all sinonasal samples
# counts          = count of each OTU in each sample --> columns are HS numbers
# counts.patient  = counts of each OTU in each patient, all sites summed together

################
# Patient info #
################

# This contains CRS/non-CRS info for some patients, as provided by Josh E.
patients <- read.table("final_patients_v3.txt",     header = TRUE, sep = "\t", stringsAsFactors = FALSE)

head(patients) # I might've incorrectly remembered the correct letter prefix as NSD

#############
# Site info #
#############

# This describes the sinonasal sites sampled
sites    <- read.table("final_sites.txt",        header = TRUE, sep = "\t", stringsAsFactors = FALSE)

head(sites) 
# Sensible == FALSE means unclear on anatomical site sampled

###############
# Sample info #
###############

# This is all the samples that were run (technical replicates have been removed from this)
samples  <- read.table("final_samples_trim.txt",      header = TRUE, sep = "\t", stringsAsFactors = FALSE)
samples  <- samples[ , c(3:6, 8)] # relevant columns only

head(samples) 
# HS number is the PRIMARY ACCESSION
# subsite is extra annotation, 
# is_followup == "y" means it is 1yr follow-up, not original

################################
# Taxonomy of all OTU detected #
################################

#This is the taxonomic report for every OTU in the sample
tax      <- read.table("final_OTU_taxonomy.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

tax[1:9, 1:9] 
# OTU for sorting, 
# count_ee1 is total count used in clustering, 
# count_all is using all primer-matched reads

####################################################
# Abundance of each OTU in each sample (HS number) # ----
#################################################### 

# This is the count table broken down for all patient X site combinations
counts <- read.table(file = "final_OTU_counts_HS.txt", sep = "\t", header = TRUE, row.names = 1)
# counts <- as.matrix(counts) # probably more useful as a matrix

counts[1:10, 1:10]

## Cross-reference HS number with "samples" to obtain patient, site, etc info.

# CJ: transpose "counts", then merge with "samples" in order to get site-by-site file?

counts.t <- as.data.frame(t(counts))
counts.t <- setDT(counts.t, keep.rownames = T)

site.counts <- merge(samples, counts.t, by.x = "hsid", by.y = "rn")

# finish off the site-by-site file for dani... but not before adding patient info

patient.counts <- merge(site.counts, patients, all = T)
raw.dani.file <- merge(patient.counts, sites, by.x = site, by.y = "Code")

raw.dani.file <- merge(patient.counts, sites, by.x = "site", by.y = "Code")

colnames(raw.dani.file)

raw.dani.file.test <- raw.dani.file[c(2, 489:491, 1, 3:5, 492:495, 6:488)]

write.csv(raw.dani.file.test, file = "raw.counts.for.dani.v2.csv")

############################################
# OTU MEANS - summarizing for each patient # ----
############################################

# calculate means for each patient using aggregate
table.1 <- aggregate(site.counts[c(6:488)], by = site.counts['patient_id'],
          function(x) mean=mean(x))

# OK, transpose it back so i can remove OTUs with low counts
means <- setNames(data.frame(t(table.1[,-1])), table.1[,1]) # this ensures that first column will swing out to column names

# sum it up
total_means <- data.frame(total=colSums(means))

# make new relative counts
rel_means  <- sweep(means, 2, total_means$total,"/")

# cut down to ones greater than 1%
otu_means_filt <- means[which(apply(rel_means, 1, max) > .01),]

# [josh e. line] change all zeros into some small number so our heatmaps don't freak out, since we want to take the log values (log of zero is bad)
otu_means_filt[otu_means_filt == 0] = .0001

# transpose means and patients
otu.means <- as.data.frame(t(otu_means_filt))

# keep rownames of new transposed dataframe
otu.means <- setDT(otu.means, keep.rownames = T)
setnames(otu.means, "rn", "subject.id")

# trimming white space to make sure no funny business with merge
otu.means$subject.id <- trimws(otu.means$subject.id)
patients$patient_id <- trimws(patients$patient_id)

# now, merge in patients
dem.means <- merge(patients, otu.means, by.x = "patient_id", by.y = "subject.id", all = T)

# now, merge in genotypes

genotypes <- read.csv("genotypes.csv", header = T, na.strings = NA)
genotypes$SubjectID <- trimws(genotypes$SubjectID)

############
# "full" dataset (that may be spotty; incomplete patients) #
############

mean.data <- merge(genotypes, dem.means, by.x = "SubjectID", by.y = "patient_id", all = T)

colnames(mean.data)

# rearrange
mean.data <- mean.data[c(1, 52:54, 2:51, 55:ncol(mean.data))]

# write this out just to have a copy
write.csv(mean.data, file="merged_data_means.csv")

# seems to be a lot missing; sum and check
mean.data$missing <- rowSums(is.na(mean.data))

# # cut down to those missing 4 cells or less (basically trimming off people who don't have genotype/microbiome data)
# # n = 46 who have genotypes + microbiome
# data.cut <- data[data$missing < 5, ]

# n = 46 who have clean data
means.complete <- mean.data[mean.data$missing < 3, ]

write.csv(means.complete, file="complete_data_means.csv")

########################################
# looped ANOVA - means of each patient # ----
########################################

means.complete <- read.csv("complete_data_means.csv", header = T)

#convert integers to factors for TAS2R38
# "1" = PAV+, "2" = PAV-
means.complete$TAS2R38 <- factor(means.complete$TAS2R38)
levels(means.complete$TAS2R38)

# check that patient status is a factor; OK
levels(means.complete$status)

# check colnames for use in loop
colnames <- as.data.frame(colnames(means.complete))

# create table first - need to make sure 'broom' package loaded
ANOVA_full <- tidy(anova(lm(means.complete[,56] ~ status*means.complete[,6], data=means.complete)))
ANOVA_full["OTU"]<-colnames(means.complete)[56]
ANOVA_full["SNP"]<-colnames(means.complete)[6]

#nested loop - going through each combination of SNP and OTU
for (a in 6:55){
  for (b in 56:110) {
    temp_full<-tidy(anova(lm(means.complete[,b] ~ status*means.complete[,a], data=means.complete))) 
    temp_full["OTU"]<-colnames(means.complete)[b]
    temp_full["SNP"]<-colnames(means.complete)[a]
    ANOVA_full<-rbind(ANOVA_full, temp_full)
  }
}

# now merge in taxonomy information with the results of the loop
tax <- setDT(tax, keep.rownames = T) # makes rownames first column

anova_tax <- merge(ANOVA_full, tax, by.x = "OTU", by.y = "rn") # ; 11004

# # merge in gene information for making the super-involved graph
# gene.info <- read.csv("gene.info.v3.csv", header = T)
# 
# # bonus round - match up genotypes (transpose first) and the gene.info from last time
# geno.t <- setNames(data.frame(t(genotypes[,-1])), genotypes[,1]) # this ensures that first column will swing out to column names
# geno.t <- setDT(geno.t, keep.rownames = T)
# 
# snp.list <- merge(geno.t, gene.info, by.x = "rn", by.y = "SNP", all = T)
# 
# colnames(snp.list)
# 
# snps <- subset(snp.list, select = c("rn", "Gene", "Dani_Sort"))
# 
# names(snps)[names(snps) == 'rn'] <- 'SNP'
# 
# write.csv(snps, file = "gene.info.revamp.csv")
# 
# # manually added a few missing SNPs and their gene / "sorting" info!!

# final gene info, with Chr and Mb info (guesstimated rs28581524 - TAS2R45)
snps <- read.csv("gene.info.v3.csv", header = T)


# replace the terms with more meaningful stuff
anova_tax$term[anova_tax$term == 'means.complete[, a]'] <- 'snp'
anova_tax$term[anova_tax$term == 'status:means.complete[, a]'] <- 'status*snp'

# merge with new gene info list
anova_tax_gene <- merge(anova_tax, snps, all = T)

# rearrange
anova_tax_gene <- anova_tax_gene[c(1, 31:34, 2:30)]

# results of ANOVA / SNP / OTU loop for MEAN abundance per patient
write.csv(anova_tax_gene, file = "prelim.anova.meanabundance.csv")

# then calculate -log and sort in excel for GRAPHPAD


########################
# 1/24/17 QQ stuff #####
########################


results <- read.csv("prelim.anova.meanabundance.csv", header = T)

# cut to the most abundant OTUs
## 3, 1, 2, 13, 4, 6, 9, 7, 12, 102, 21
abundant <- c(3, 1, 2, 13, 4, 6, 9, 7, 12, 102, 21)

results.11 <- results[results$OTU.1 %in% abundant, ]

#cut so only interaction term
results.interact <- results.11[which(results.11$term=="status*snp"), ]

# cut the genotypes with dani_sort = 0
qqdata <- results.interact[which(results.interact$Dani_Sort > 0), ]

#define a nice qq plot function - "coward's way" as dani so kindly put it :p
ggd.qqplot = function(pvector, main=NULL, ...) {
  o = -log10(sort(pvector,decreasing=F))
  e = -log10( 1:length(o)/length(o) )
  plot(e,o,pch=19,cex=1, main=main, ...,
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       xlim=c(0,max(e)), ylim=c(0,max(o)))
  lines(e,e,col="red")
}

# now sort and plot using fancy function
ggd.qqplot(qqdata$p.value)


# how about dani's QQ code
observed <- sort(qqdata$p.value)
lobs <- -(log10(observed))

expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))

plot(c(0,12), c(0,12), col="red", lwd=3, type="l",
     xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,5), ylim=c(0,12), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black")



########################################################################
# 1/24/17 - for manuscript, seeing which bitter markers have a main effect #####
########################################################################

# take the top 11 OTU results, cut to main effect of snp
results.bitter <- results.11[which(results.11$term=="snp"), ]

# sort on p
results.bitter <- results.bitter[with(results.bitter, order(p.value)), ]

# cut down to significant ones, p < 0.05
bitter.sig <- results.bitter[which(results.bitter$p.value < 0.05), ]

# make a list of the unique markers that are significant - to list in manuscript
sig.markers <- as.data.frame(unique(bitter.sig$SNP))

#### or am i doing this the wrong way? do i need to cut down to TAS2R38 (PAV-/PAV+), and then see which OTUs are sig? ####
tas2r38.results <- results.bitter[which(results.bitter$SNP=="TAS2R38"), ] # nothing is significant!

# how about just checking out allllll OTUs
tas2r38.all.results <- results[which(results$SNP=="TAS2R38"), ]
tas2r38.all.results <- tas2r38.all.results[with(tas2r38.all.results, order(p.value)), ]

# OTU11 is only significant one!

# need to boxplot it
library(ggplot2)
ggplot(means.complete, aes(TAS2R38, OTU_11)) +
  geom_boxplot()


#########################################
# here is a straight c/p of dani's QQ code #####
#########################################


#Dani and microbiome
library(heplots)
#from Excel file raw.counts.for.dani
#summed all raw counts per otu
#recordered columns so must abundance otus are first
#most abundant otus by this criterion are 3, 1,2, 13, 4,6,9,7, 12,102,21
setwd("C:/Users/Reed/Dropbox/2017_T2R_nasal_microbiome/2000. Microbiome data/DRR")

raw.counts.for.daniv2 <- read.csv("C:/Users/Reed/Dropbox/2017_T2R_nasal_microbiome/2000. Microbiome data/DRR/raw.counts.for.daniv2.csv")

df <- raw.counts.for.daniv2
summary(df$site)
str(df)
table(df$Type)
df1 <- df[df$Type %in% c("Biopsy", "Swab"),]
df <- df1
Y <- cbind(df$OTU_3, df$OTU_1, df$OTU_2, df$OTU_13, df$OTU_4, df$OTU6, df$OTU_9, df$OTU_7, df$OTU_12, df$OTU_102, df$OTU_21)
fit1 <- manova(Y ~ df$Type)
summary(fit1, test="Pillai")

mean <- aggregate(df[, 10:20], list(df$Type), mean)
mean

df_nobio <- df[df$Type %in% c("Swab"),]
Y <- cbind(df_nobio$OTU_3, df_nobio$OTU_1, df_nobio$OTU_2, df_nobio$OTU_13, df_nobio$OTU_4, df_nobio$OTU6, df_nobio$OTU_9, df_nobio$OTU_7, df_nobio$OTU_12, df_nobio$OTU_102, df_nobio$OTU_21)

#need topick up case/control status -- grrrr....
final_patients_v3 <- read.delim("C:/Users/Reed/Dropbox/2017_T2R_nasal_microbiome/2000. Microbiome data/CJA/final_patients_v3.txt")
df2 <- merge(final_patients_v3, df_nobio)

#reduce to essential sites A C E G I K
df3 <- df_lim_sites <- df2[df2$site %in% c("A", "C", "E", "G", "I", "R"),]
Y1 <- cbind(df_lim_sites$OTU_3, df_lim_sites$OTU_1, df_lim_sites$OTU_2, df_lim_sites$OTU_13, df_lim_sites$OTU_4, df_lim_sites$OTU6, df_lim_sites$OTU_9, df_lim_sites$OTU_7, df_lim_sites$OTU_12, df_lim_sites$OTU_102, df_lim_sites$OTU_21)
fit2 <- manova(Y1 ~ df3$site + df3$status)
summary(fit2, test="Pillai")
etasq(fit2)

#fit <- lm(OTU_3 ~ site, data=df)
#summary(fit) # show results

###############
# intro stats - CJA checked 1/23/2017 #####
###############

df <- raw.dani.file.test

table(df$Type)

# cuts out one Unknown
df <- df[df$Type %in% c("Biopsy", "Swab"),]

# let's cut down to just OTUs

# sum them up
colnames(df)

otus.only <- df[c(13:ncol(df))]
sums <- as.data.frame(colSums(otus.only))

# yep, 3, 1, 2, 13, 4, 6, 9, 7, 12 102, 21
dvs.raw <- cbind(df$OTU_3, df$OTU_1, df$OTU_2, df$OTU_13, df$OTU_4, df$OTU6, df$OTU_9, df$OTU_7, df$OTU_12, df$OTU_102, df$OTU_21)

manova.1 <- manova(dvs.raw ~ df$Type)

# for first swab v. biospy 'intro' stats
summary(manova.1, test="Pillai")

# now, focus on swabs only
df_nobio <- df[df$Type %in% c("Swab"),]

# cut down to certain sites
df_nobio <- df_nobio[df_nobio$site %in% c("A", "C", "E", "G", "I", "R"),]

# create list of DVs
dvs.nobio <- cbind(df_nobio$OTU_3, df_nobio$OTU_1, df_nobio$OTU_2, df_nobio$OTU_13, df_nobio$OTU_4, df_nobio$OTU6, df_nobio$OTU_9, df_nobio$OTU_7, df_nobio$OTU_12, df_nobio$OTU_102, df_nobio$OTU_21)

manova.2 <- manova(dvs.nobio ~ df_nobio$site + df_nobio$status)

summary(manova.2, test="Pillai")

library(heplots)

etasq(manova.2)


###############
# make a two-way table for pseudomones result, using TAS2R38 and CRS status ####
###############

table <- table(means.complete$TAS2R38, means.complete$status)

prop.table(table)

summary(table)

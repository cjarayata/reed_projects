# Adip20 R script - last update Aug 30th 2016 #

# Developed by Charles J Arayata, Danielle R Reed and Cailu Lin
# Monell Chemical Senses Center


## Prepping the session: Packages and Data Files ## ----

#use 'pacman' package to install and load packages needed
if (!require("pacman")) install.packages("pacman", "CRAN")
library(pacman)
pacman::p_load(psych, lubridate, plyr, dplyr, broom, reshape, data.table, lsr)

#Read in data
Data_near_final <- read.csv("Data_near_final_current_8_29_2016.csv", header=TRUE)

#Read in SNP info
snp <- read.csv("snp.csv", header=TRUE)
snp$SNPs <- tolower(snp$SNPs)



## Table 1 ## ----
## Characteristics of N=2,137 male mice used in Adip20 mapping studies ##
## Drop Inbrefd HF - n=19, final N=2,118 ##


#start with n's for each strain
n1 <- summary(Data_near_final$Table1_Mapping_population)
n1 <- data.frame(n1)
setDT(n1, keep.rownames=TRUE)
names(n1) <- c("Table1_Mapping_population", "n")

#get rest of descriptives
tbl1.min <- ddply(Data_near_final, "Table1_Mapping_population", summarize, min=min(Days))
tbl1.max <- ddply(Data_near_final, "Table1_Mapping_population", summarize, max=max(Days))
tbl1.mean <- ddply(Data_near_final, "Table1_Mapping_population", summarize, mean=mean(Days))
tbl1.sd <- ddply(Data_near_final, "Table1_Mapping_population", summarize, sd=sd(Days))

#merge together
table.1 <- Reduce(function(...) merge(..., by="Table1_Mapping_population", all=TRUE), list(n1, tbl1.min, tbl1.max, tbl1.mean, tbl1.sd))

#need to get START and END dates
#first subset data to cut out genotypes; then make all variable names lowercase for ease of variable reference
cj.data1 <- Data_near_final[ , 1:25]
colnames(cj.data1) <- tolower(names(cj.data1))

#parse out birthdate into new variable - tacks onto end of data frame (which is why we just created a new data frame)
cj.data1$birthdate.fix <- parse_date_time(cj.data1$birthdate, orders="mdy")

#get min and max dates for each strain
t1.bday.min <- ddply(cj.data1, "table1_mapping_population", summarize, start=min(birthdate.fix))
t1.bday.max <- ddply(cj.data1, "table1_mapping_population", summarize, end=max(birthdate.fix))
t1.bday <- merge(t1.bday.min, t1.bday.max)

#merge descriptives and start/end dates together - now complete
table.1 <- merge(table.1, t1.bday, by.x="Table1_Mapping_population", by.y="table1_mapping_population", all=TRUE)

## Table 1 END ##



## Table 2 ## ----
## Characteristics of N=1,293 congenic mice by N=22 strains ##

#for n's
n2 <- summary(Data_near_final$Table_2_Substrains)
n2 <- data.frame(n2)
setDT(n2, keep.rownames=TRUE)
names(n2) <- c("Table_2_Substrains", "n")

#for 'Day' descriptives
tbl2.min <- ddply(Data_near_final, "Table_2_Substrains", summarize, min=min(Days))
tbl2.max <- ddply(Data_near_final, "Table_2_Substrains", summarize, max=max(Days))
tbl2.mean <- ddply(Data_near_final, "Table_2_Substrains", summarize, mean=mean(Days))
tbl2.sd <- ddply(Data_near_final, "Table_2_Substrains", summarize, sd=sd(Days))

#merge all together
table.2 <- Reduce(function(...) merge(..., by="Table_2_Substrains", all=TRUE), list(n2, tbl2.min, tbl2.max, tbl2.mean, tbl2.sd))

## Table 2 END ##



## Table 3 ## ----
## Congenic gonadal depot weight comparing host and full-length donor regions by strain ##

#n's - count the number of host vs donor mice (full length donor)
tbl3_n <- as.data.frame(xtabs(~ Table_2_Substrains + Full_length_donor, data=Data_near_final))

#subset based on strains to actually use
table.3 <- Data_near_final[Data_near_final$Table_2_Substrains %in% c(1, 1.1, 3, 3.1, "3.1.1", "3.1.1.1", "3.1.1.2", "3.1.1.4", "3.1.4.1", 4, 4.1, 4.4), ]

#there are unused levels lurking , 'refactoring' the variable will drop all empty strain groups
levels(table.3$Table_2_Substrains)
table.3$Table_2_Substrains <- factor(table.3$Table_2_Substrains)
levels(table.3$Table_2_Substrains)

# ANOVAs #
#set full donor status as factor for analysis
full.donor <- factor(table.3$Full_length_donor, levels = c("yes","host"))

#build table first - aov with partial eta, put in strain info, and rearrange table
aov.out3 <- aov(Gonadal ~ full.donor + BW, data=table.3, subset=Table_2_Substrains=="1")
aov.out3 <- tidy(etaSquared(aov.out3, type=1, anova=TRUE))
aov.out3["strain"] <- "1"
aov.out3 <- aov.out3[, c(9, 1, 4:8, 2, 3)]

#loop through all strains
for (strain in levels(table.3$Table_2_Substrains)){
  temp <- aov(Gonadal ~ full.donor + BW, data=table.3, subset=Table_2_Substrains==strain)
  temp <- tidy(etaSquared(temp, type=1, anova=TRUE))
  temp["strain"] <- strain
  temp <- temp[, c(9, 1, 4:8, 2, 3)]
  aov.out3 <- rbind(aov.out3, temp)
}

#cut first 3 rows - duplicates from initial table set up
aov.out3 <- aov.out3[4:nrow(aov.out3), ]

## Table 3 - END ##



## Supplemental Table 5 ## ----
## Congenic gonadal depot weight comparing host and donor regions (partial or full) by strain ##
# adds additional substrains 1.1.1, 1.2, 3.1.1.3, 3.1.2, 3.1.3, 4.1a, 4.2, 4.3, 4.5 

#count n's
tbl5_n <- as.data.frame(xtabs(~ Table_2_Substrains + Any_length_donor, data=Data_near_final))

#ANOVAs
#set partial donor status as factor for analysis
any.donor <- factor(Data_near_final$Any_length_donor, levels = c("donor","host"))

#build table first - aov with partial eta, put in strain info, and rearrange table
aov.out5 <- aov(Gonadal ~ any.donor + BW, data=Data_near_final, subset=Table_2_Substrains=="1")
aov.out5 <- tidy(etaSquared(aov.out5, type=1, anova=TRUE))
aov.out5["strain"] <- "1"
aov.out5 <- aov.out5[, c(9, 1, 4:8, 2, 3)]

#this loop breaks at strain 3.1.3, because 0/2 mice in groups
for (strain in levels(Data_near_final$Table_2_Substrains)){
temp <- aov(Gonadal ~ any.donor + BW, data=Data_near_final, subset=Table_2_Substrains==strain)
temp <- tidy(etaSquared(temp, type=1, anova=TRUE))
temp["strain"] <- strain
temp <- temp[, c(9, 1, 4:8, 2, 3)]
aov.out5 <- rbind(aov.out5, temp)
}

#continue loop where it left off
for (strain in levels(Data_near_final$Table_2_Substrains)[14:23]){
  temp <- aov(Gonadal ~ any.donor + BW, data=Data_near_final, subset=Table_2_Substrains==strain)
  temp <- tidy(etaSquared(temp, type=1, anova=TRUE))
  temp["strain"] <- strain
  temp <- temp[, c(9, 1, 4:8, 2, 3)]
  aov.out5 <- rbind(aov.out5, temp)
}

#cut first 3 rows - duplicates from initial table set up
aov.out5 <- aov.out5[4:nrow(aov.out5), ]

## Supplemental Table 5 - END ##



## Supplemental Table 4 ## ----
## Weight of gonadal depot and body weight by congenic strain and genotype group ##

#n's - same code as above for Table 3 (not necessary to run below)
tbl3_n <- as.data.frame(xtabs(~ Table_2_Substrains + Full_length_donor, data=Data_near_final))

#get means and SDs
s4_means <- aggregate(Data_near_final[, 23:25], list(Data_near_final$Fulldonor_vs_Host), mean)
s4_SD <- aggregate(Data_near_final[, 23:25], list(Data_near_final$Fulldonor_vs_Host), sd)

#column bind both tables, then cut/rearrange/rename variables
sup.4 <- cbind(s4_means, s4_SD)
sup.4 <- sup.4[, c(1, 4, 8, 2, 6)]
names(sup.4) <- c("Strain", "Gonadal.Mean", "Gonadal.SD", "Body.Weight.Mean", "Body.Weight.SD")
View(sup.4)

## Supplemental Table 4 - END ##



## In-text analyses - ANCOVA ## ----

#Dani method
#Please note that default is type 1 Sum of Squares and this is different than default in STATISTICA
##Clarification - "glm" call outputs Type VI (unique), but using 'anova' to print output of F-test will give Type I (sequential)
B6.129_Chr9 <- Data_near_final[Data_near_final$Groups %in% c("B6.129-Chr9", "C57BL/6ByJ"),]
Model1 <- glm(Gonadal ~ Groups + BW, data=B6.129_Chr9)
anova(Model1, test = "F")
Consomics <- Data_near_final[Data_near_final$Groups %in% c("129.B6-Chr9", "129P3/J"),]
Model2 <- glm(Gonadal ~ Groups + BW, data=Consomics)
anova(Model2, test = "F")

#CJA method
#note use of 'aov' - default is Type I (Sequential) Sums of Squares
b6_chr9 <- Data_near_final[which(Data_near_final$Groups=='B6.129-Chr9' | Data_near_final$Groups=='C57BL/6ByJ'), ]
cj.model1 <- aov(Gonadal ~ Groups + BW, data=b6_chr9)
summary(cj.model1)

cj.consomics <- Data_near_final[which(Data_near_final$Groups=='129.B6-Chr9' | Data_near_final$Groups=='129P3/J'), ]
cj.model2 <- aov(Gonadal ~ Groups + BW, data=cj.consomics)
summary(cj.model2)

## In-text analyses - ANCOVA END ##



## Figure 2 Behemoth ## ----
## -log values across gemone, plus gonadal data overlayed ##

# step 1: -log values #

#load in data
cj.figure2 <- Data_near_final
colnames(cj.figure2) <- tolower(names(cj.figure2))
names(cj.figure2)

## below is ORIGINAL code used to produce p-values

#create table first
aov.out2 <- tidy(aov(lm(gonadal ~ cj.figure2[,26] + bw, data=cj.figure2, subset=table1_mapping_population=='F2')))
aov.out2["SNP"] <- colnames(cj.figure2)[26]
aov.out2["Strain"] <- "F2"

#write one big loop - loop through SNPs for each strain subsetted, then bind
for (i in 27:length(names(cj.figure2))){
  temp.f2<-tidy(aov(lm(gonadal~cj.figure2[,i] + bw, data=cj.figure2, subset=table1_mapping_population=='F2')))
  temp.f2["SNP"] <- colnames(cj.figure2)[i]
  temp.f2["Strain"] <- "F2"
  aov.out2 <- rbind(aov.out2, temp.f2)
  
  temp.f2.hf<-tidy(aov(lm(gonadal~cj.figure2[,i] + bw, data=cj.figure2, subset=table1_mapping_population=='F2-HF')))
  temp.f2.hf["SNP"] <- colnames(cj.figure2)[i]
  temp.f2.hf["Strain"] <- "F2-HF"
  aov.out2 <- rbind(aov.out2, temp.f2.hf)
  
  temp.f1.129<-tidy(aov(lm(gonadal~cj.figure2[,i] + bw, data=cj.figure2, subset=table1_mapping_population=='N2 (F1 x 129)')))
  temp.f1.129["SNP"] <- colnames(cj.figure2)[i]
  temp.f1.129["Strain"] <- "N2 (F1 x 129)"
  aov.out2 <- rbind(aov.out2, temp.f1.129)
  
  temp.f1.b6<-tidy(aov(lm(gonadal~cj.figure2[,i] + bw, data=cj.figure2, subset=table1_mapping_population=='N2 (F1 x B6)')))
  temp.f1.b6["SNP"] <- colnames(cj.figure2)[i]
  temp.f1.b6["Strain"] <- "N2 (F1 x B6)"
  aov.out2 <- rbind(aov.out2, temp.f1.b6)
  
  temp.n3<-tidy(aov(lm(gonadal~cj.figure2[,i] + bw, data=cj.figure2, subset=table1_mapping_population=='N3')))
  temp.n3["SNP"] <- colnames(cj.figure2)[i]
  temp.n3["Strain"] <- "N3"
  aov.out2 <- rbind(aov.out2, temp.n3)
  
  temp.n4<-tidy(aov(lm(gonadal~cj.figure2[,i] + bw, data=cj.figure2, subset=table1_mapping_population=='N4')))
  temp.n4["SNP"] <- colnames(cj.figure2)[i]
  temp.n4["Strain"] <- "N4"
  aov.out2 <- rbind(aov.out2, temp.n4)
  
  temp.n6<-tidy(aov(lm(gonadal~cj.figure2[,i] + bw, data=cj.figure2, subset=table1_mapping_population=='N6')))
  temp.n6["SNP"] <- colnames(cj.figure2)[i]
  temp.n6["Strain"] <- "N6"
  aov.out2 <- rbind(aov.out2, temp.n6)
  
}

#merge in SNP info
fig2.log <- merge(aov.out2, snp, by.x="SNP", by.y="SNPs")

#rearrange results file
fig2.log <- fig2.log[, c(2:8, 1, 9, 10)]

#calculate -log (base 10)
fig2.log$neg.logp <- (-log10(fig2.log$p.value))


## edited 9/30/16 to get etaSquared values from ANOVA test ## ----
## HOWEVER, results of this ANOVA test yields '0' for extremely small p-values, therefore does not calcuate -log properly!


#create table first
aov.out2 <- aov(gonadal ~ cj.figure2[,27] + bw, data=cj.figure2, subset=table1_mapping_population=='F2')
aov.out2 <- tidy(etaSquared(aov.out2, type=1, anova=T))
aov.out2["SNP"] <- colnames(cj.figure2)[27]
aov.out2["Strain"] <- "F2"

#write one big loop - loop through SNPs for each strain subsetted, then bind
for (i in 27:length(names(cj.figure2))){
  temp.f2 <- aov(gonadal~cj.figure2[,i] + bw, data=cj.figure2, subset=table1_mapping_population=='F2')
  temp.f2 <- tidy(etaSquared(temp.f2, type=1, anova=T))
  temp.f2["SNP"] <- colnames(cj.figure2)[i]
  temp.f2["Strain"] <- "F2"
  aov.out2 <- rbind(aov.out2, temp.f2)
  
  temp.f2.hf <- aov(gonadal~cj.figure2[,i] + bw, data=cj.figure2, subset=table1_mapping_population=='F2-HF')
  temp.f2.hf <- tidy(etaSquared(temp.f2.hf, type=1, anova=T))
  temp.f2.hf["SNP"] <- colnames(cj.figure2)[i]
  temp.f2.hf["Strain"] <- "F2-HF"
  aov.out2 <- rbind(aov.out2, temp.f2.hf)
  
  temp.f1.129 <- aov(gonadal~cj.figure2[,i] + bw, data=cj.figure2, subset=table1_mapping_population=='N2 (F1 x 129)')
  temp.f1.129 <- tidy(etaSquared(temp.f1.129, type=1, anova=T))
  temp.f1.129["SNP"] <- colnames(cj.figure2)[i]
  temp.f1.129["Strain"] <- "N2 (F1 x 129)"
  aov.out2 <- rbind(aov.out2, temp.f1.129)
  
  temp.f1.b6 <- aov(gonadal~cj.figure2[,i] + bw, data=cj.figure2, subset=table1_mapping_population=='N2 (F1 x B6)')
  temp.f1.b6 <- tidy(etaSquared(temp.f1.b6, type=1, anova=T))
  temp.f1.b6["SNP"] <- colnames(cj.figure2)[i]
  temp.f1.b6["Strain"] <- "N2 (F1 x B6)"
  aov.out2 <- rbind(aov.out2, temp.f1.b6)
  
  temp.n3 <- aov(gonadal~cj.figure2[,i] + bw, data=cj.figure2, subset=table1_mapping_population=='N3')
  temp.n3 <- tidy(etaSquared(temp.n3, type=1, anova=T))
  temp.n3["SNP"] <- colnames(cj.figure2)[i]
  temp.n3["Strain"] <- "N3"
  aov.out2 <- rbind(aov.out2, temp.n3)
  
  temp.n4 <- aov(gonadal~cj.figure2[,i] + bw, data=cj.figure2, subset=table1_mapping_population=='N4')
  temp.n4 <- tidy(etaSquared(temp.n4, type=1, anova=T))
  temp.n4["SNP"] <- colnames(cj.figure2)[i]
  temp.n4["Strain"] <- "N4"
  aov.out2 <- rbind(aov.out2, temp.n4)
  
  temp.n6 <- aov(gonadal~cj.figure2[,i] + bw, data=cj.figure2, subset=table1_mapping_population=='N6')
  temp.n6 <- tidy(etaSquared(temp.n6, type=1, anova=T))
  temp.n6["SNP"] <- colnames(cj.figure2)[i]
  temp.n6["Strain"] <- "N6"
  aov.out2 <- rbind(aov.out2, temp.n6)

}

#merge in SNP info
fig2.log <- merge(aov.out2, snp, by.x="SNP", by.y="SNPs")

#rearrange results file
fig2.log <- fig2.log[, c(2, 6, 5, 7, 8, 9, 3:4, 10, 11, 12, 1)]


#calculate -log (base 10)
fig2.log$neg.logp <- (-log10(fig2.log$p))

# step 1 END #

# step 2: gonadal weights by strain and genotype

#build table first
## need to use "i" even though just represents one SNP (column 27) so that colnames will match during loop

for (i in 27:27) {
fig2.table <- ddply(cj.figure2, .(cj.figure2[,2], cj.figure2[,i]), summarize, gonadal.mean=mean(gonadal, na.rm=T))
fig2.table["SNP"] <- colnames(cj.figure2)[i]
}

#loop and rbind results together

for (i in 28:ncol(cj.figure2)) {
  fig2.temp <- ddply(cj.figure2, .(cj.figure2[,2], cj.figure2[,i]), summarize, gonadal.mean=mean(gonadal, na.rm=T))
  fig2.temp["SNP"] <- colnames(cj.figure2)[i]
  fig2.table <- rbind(fig2.temp, fig2.table)
}

#table now built, but need to cast for GRAPHPAD
#define columns; cast data
colnames(fig2.table) <- c("strain", "genotype", "gonadal.mean", "snp")
fig2.cast <- reshape(fig2.table,
                      timevar="genotype",
                      idvar=c("strain", "snp"),
                      direction="wide")

#merge casted table with snp info
fig2.cast.snp <- merge(fig2.cast, snp, by.x="snp", by.y="SNPs")

#rearrange and export - ready for graphpad once sorted/exported
fig2.cast.snp <- fig2.cast.snp[, c(2, 1, 7, 8, 3:6)]

## Figure 2 - END ##



## Figure 3 ## ----
## Graph of in-text analyses for consomics - gonadal weight ##

#start with n's for each strain
n3 <- summary(Data_near_final$Groups)
n3 <- data.frame(n3)
setDT(n3, keep.rownames=TRUE)
names(n3) <- c("Groups", "n")

#get rest of descriptives
fig3.min <- ddply(Data_near_final, "Groups", summarize, min=min(Gonadal))
fig3.max <- ddply(Data_near_final, "Groups", summarize, max=max(Gonadal))
fig3.mean <- ddply(Data_near_final, "Groups", summarize, mean=mean(Gonadal))
fig3.sd <- ddply(Data_near_final, "Groups", summarize, sd=sd(Gonadal))

#merge all together - ready for graphpad
figure.3 <- Reduce(function(...) merge(..., by="Groups", all=TRUE), list(n3, fig3.min, fig3.max, fig3.mean, fig3.sd))

## Figure 3 END ##



## August 16, 2016 - CJA Figure 2 tweaking ##
## combine N3, N4, N6 into one panel - n=159 ## ----

#subset
n_comb <- Data_near_final[Data_near_final$Table1_Mapping_population %in% c("N3", "N4", "N6"), ]
colnames(n_comb) <- tolower(names(n_comb))

## below is ORIGINAL code

# step 1 - -log values
#create table first
aov.ncomb <- tidy(aov(lm(gonadal ~ n_comb[,26] + bw, data=n_comb)))
aov.ncomb["SNP"] <- colnames(n_comb)[26]

#write one big loop - loop through SNPs, then bind
for (i in 27:length(names(n_comb))){
  temp.ncomb<-tidy(aov(lm(gonadal~n_comb[,i] + bw, data=n_comb)))
  temp.ncomb["SNP"] <- colnames(n_comb)[i]
  aov.ncomb <- rbind(aov.ncomb, temp.ncomb)
}

#merge in SNP info
ncomb.results <- merge(aov.ncomb, snp, by.x="SNP", by.y="SNPs")

#rearrange results file
ncomb.results <- ncomb.results[, c(2:9, 1)]

#calculate -log (base 10)
ncomb.results$neg.logp <- (-log10(ncomb.results$p.value))


## edited 9/30/16 to get etaSquared from ANOVA tests - but p-value issue ##


# step 1 - -log values
#create table first
aov.ncomb <- aov(gonadal ~ n_comb[,27] + bw, data=n_comb)
aov.ncomb <- tidy(etaSquared(aov.ncomb, type=1, anova=T))
aov.ncomb["SNP"] <- colnames(n_comb)[27]

#write one big loop - loop through SNPs, then bind
for (i in 28:length(names(n_comb))){
  temp.ncomb <- aov(gonadal~n_comb[,i] + bw, data=n_comb)
  temp.ncomb <- tidy(etaSquared(temp.ncomb, type=1, anova=T))
  temp.ncomb["SNP"] <- colnames(n_comb)[i]
  aov.ncomb <- rbind(aov.ncomb, temp.ncomb)
}

#merge in SNP info
ncomb.results <- merge(aov.ncomb, snp, by.x="SNP", by.y="SNPs")

#rearrange results file
ncomb.results <- ncomb.results[, c(2, 6, 5, 7, 8, 9, 3:4, 10, 11, 1)]

#calculate -log (base 10)
ncomb.results$neg.logp <- (-log10(ncomb.results$p))


# step 2: gonadal weights by genotype only
#build table first with "i" so column names match
for (i in 27:27) {
  fig2.table.ncomb <- ddply(n_comb, .(n_comb[,i]), summarize, gonadal.mean=mean(gonadal, na.rm=T))
  fig2.table.ncomb["SNP"] <- colnames(n_comb)[i]
}

#loop and rbind results together
for (i in 28:ncol(n_comb)) {
  fig2.temp.ncomb <- ddply(n_comb, .(n_comb[,i]), summarize, gonadal.mean=mean(gonadal, na.rm=T))
  fig2.temp.ncomb["SNP"] <- colnames(n_comb)[i]
  fig2.table.ncomb <- rbind(fig2.temp.ncomb, fig2.table.ncomb)
}

#table built; redefine columns and cast for GRAPHPAD
colnames(fig2.table.ncomb) <- c("genotype", "gonadal.mean", "snp")

fig2.ncomb.cast <- reshape(fig2.table.ncomb,
                     timevar="genotype",
                     idvar="snp",
                     direction="wide")

#merge casted table with snp info
fig2.ncomb.cast.snp <- merge(fig2.ncomb.cast, snp, by.x="snp", by.y="SNPs")

#rearrange - ready for graphpad once sorted
fig2.ncomb.cast.snp <- fig2.ncomb.cast.snp[, c(1, 5, 6, 2:4)]

## N combined END ##


## August 16, 2016 - CJA Figure 2 tweaking ##
## run same pooled analyses for all congenic mice; n=1293 ## ----

#subset
congenic <- Data_near_final[Data_near_final$Table1_Mapping_population=="Congenic", ]
colnames(congenic) <- tolower(names(congenic))

## below is ORIGINAL code

#step 1 - -log values
#create table first
aov.congenic <- tidy(aov(lm(gonadal ~ congenic[,26] + bw, data=congenic)))
aov.congenic["SNP"] <- colnames(congenic)[26]

#loop through SNPs
for (i in 27:length(names(congenic))){
  temp.congenic<-tidy(aov(lm(gonadal~congenic[,i] + bw, data=congenic)))
  temp.congenic["SNP"] <- colnames(congenic)[i]
  aov.congenic <- rbind(aov.congenic, temp.congenic)
}

#merge in SNP info; rearrange, -log
congenic.results <- merge(aov.congenic, snp, by.x="SNP", by.y="SNPs")
congenic.results <- congenic.results[, c(2:9, 1)]
congenic.results$neg.logp <- (-log10(congenic.results$p.value))


## edited 9/30/16 to get etaSquared from ANOVA tests - but p-value issue ##

#step 1 - -log values
#create table first
aov.congenic <- aov(gonadal ~ congenic[,27] + bw, data=congenic)
aov.congenic <- tidy(etaSquared(aov.congenic, type=1, anova=T))
aov.congenic["SNP"] <- colnames(congenic)[27]

#loop through SNPs
for (i in 28:length(names(congenic))){
  temp.congenic <- aov(gonadal~congenic[,i] + bw, data=congenic)
  temp.congenic <- tidy(etaSquared(temp.congenic, type=1, anova=T))
  temp.congenic["SNP"] <- colnames(congenic)[i]
  aov.congenic <- rbind(aov.congenic, temp.congenic)
}

#merge in SNP info; rearrange, -log
congenic.results <- merge(aov.congenic, snp, by.x="SNP", by.y="SNPs")
congenic.results <- congenic.results[, c(2, 6, 5, 7, 8, 9, 3:4, 10, 11, 1)]

congenic.results$neg.logp <- (-log10(congenic.results$p))

# step 2: gonadal weights by strain and genotype

#same method as above; build, then loop and bind
for (i in 27:27) {
  fig2.table.congenic <- ddply(congenic, .(congenic[,i]), summarize, gonadal.mean=mean(gonadal, na.rm=T))
  fig2.table.congenic["SNP"] <- colnames(congenic)[i]
}

for (i in 28:ncol(congenic)) {
  fig2.temp.congenic <- ddply(congenic, .(congenic[,i]), summarize, gonadal.mean=mean(gonadal, na.rm=T))
  fig2.temp.congenic["SNP"] <- colnames(congenic)[i]
  fig2.table.congenic <- rbind(fig2.temp.congenic, fig2.table.congenic)
}

#table built; redefine columns and cast for GRAPHPAD
colnames(fig2.table.congenic) <- c("genotype", "gonadal.mean", "snp")
fig2.congenic.cast <- reshape(fig2.table.congenic,
                           timevar="genotype",
                           idvar="snp",
                           direction="wide")

#merge casted table with snp info
fig2.congenic.cast.snp <- merge(fig2.congenic.cast, snp, by.x="snp", by.y="SNPs")

#rearrange and export - ready for graphpad once sorted
fig2.congenic.cast.snp <- fig2.congenic.cast.snp[, c(1, 6, 7, 2:5)]

# 9/30/16 - why are there NAs there?!?!

## consomics END
## August 16, 2016 - CJA Figure 2 tweaking END##


## Write all files into new 'checking' working directory ## ----
setwd("S:/bachmanov/Labshare/Publications&Presentations/Monell-InPrep/Adip20/7.  Near submission/CJ Effects Sizes")

write.csv(table.1, file="Table_1.csv")
write.csv(table.2, file="Table_2.csv")
write.csv(tbl3_n, file="Table_3_n.csv")
write.csv(aov.out3, file="Table_3_ANOVAs.csv")
write.csv(tbl5_n, file="Supplemental_Table_5_n.csv")
write.csv(aov.out5, file="Supplemental_Table_5_ANOVAS.csv")
write.csv(sup.4, file="Supplemental_Table_4_means.csv")
write.csv(figure.3, file="Figure_3_means.csv")

#Figure 2: original
write.csv(fig2.log, "Figure_2_log 9-30-16.csv")
write.csv(fig2.cast.snp, file="Figure_2_means 9-30-16.csv")

#Figure 2: combined N and consomics
write.csv(ncomb.results, "Figure_2_log_N_combined - 9-30-16.csv")
write.csv(fig2.ncomb.cast.snp, file="Figure_2_Ncombined_means 9-30-16.csv")

write.csv(congenic.results, "Figure_2_log_congenics - 9-30-16.csv")
write.csv(fig2.congenic.cast.snp, file="Figure_2_congenic_means 9-30-16.csv")


## Codebook preparation - Aug 30 2016 ## ----

var.class <- as.data.frame(sapply(Data_near_final, class))
var.levels <- as.data.frame(sapply(Data_near_final, nlevels))
vars <- cbind(var.class, var.levels)
colnames(vars) <- c("var.class", "num.levels")

#write file; need to fill variable labels and actual variable levels MANUALLY
write.csv(vars, file="codebook_start.csv")

#read in completed codebook
codebook <- read.csv("Adip20_Codebook - Aug 30 2016.csv", header=T, na.string="")


# Adip20 Script END #
# Burly1 R script - 2018-02-05 #

# Developed by Charles J Arayata, Danielle R Reed and Cailu Lin
# Monell Chemical Senses Center

## Prepping the session: Packages and Data Files ####

# Use 'pacman' package to install and load packages needed
if (!require("pacman")) install.packages("pacman", "CRAN")

pacman::p_load(psych, lubridate, plyr, dplyr, broom, reshape, data.table, lsr, scatterplot3d, multcompView, lsmeans, ggplot2, cowplot, ggpubr, agricolae, DescTools)

# Set to your working directory
setwd("~/GitHub/reed_projects/Burly1")

# Read in data
Data_burly1 <- read.csv("Data/Data_burly1.csv", header=TRUE)
# Data_burly1 <- read.csv("Data/Data_burly1.csv", header=TRUE, stringsAsFactors = F)

# Read in SNP info
# Duplicated as S2 Table: Polymorphic mChr2 markers
snp <- read.csv("Data/snp.csv", header=TRUE)
snp$Marker <- trimws(tolower(snp$Marker))

# Standard Error function
sem <- function(x){
  sd(x)/sqrt(length(x))
}

# Negate "not in" function
"%not%" <- Negate("%in%")

## Prep END ##

## Table 1: Characteristics of N=2,053 male mice used in Burly1 mapping studies ####

# Massage dates via lubridate for start/end dates
Data_burly1$birthdate.fix <- parse_date_time(Data_burly1$Birth.Date, orders="mdy")
Data_burly1$endpoint.fix <- parse_date_time(Data_burly1$Endpoint.Date, orders="mdy")

# Create table
table.1 <- ddply(Data_burly1, ~Table1_Mapping_population, summarise,
                 n = length(Mouse.ID),
                 min.age = min(Age),
                 max.age = max(Age),
                 mean.MRage = mean(Age),
                 sd.MRage = sd(Age),
                 start = min(birthdate.fix),
                 end = max(birthdate.fix))
                 # end = max(endpoint.fix)) # is endpoint correct?

write.table(table.1, file="Results/Table1.csv", sep=",", row.names=F)
## Table 1 END ##


## Table 2: Characteristics of N=1,293 congenic mice by N=22 strains ####

## Exclude mice
# n=1 mouse from substrain 2, endpoint = 100
# n=2 mice from 2.2, endpoint 120
bad.ids <-
  c(as.character(Data_burly1$Mouse.ID[which(Data_burly1$Table_2_Substrains == "2" &
                                              Data_burly1$Age == "100")]),
    as.character(Data_burly1$Mouse.ID[which(Data_burly1$Table_2_Substrains == "2.2" &
                                              Data_burly1$Age == "120")]))

# Subset data, NOT bad.ids = good IDs
table2.data <- Data_burly1[which(Data_burly1$Mouse.ID %not% bad.ids), ]

# Create table
table.2 <- ddply(table2.data, ~Table_2_Substrains, summarise,
                 n = length(Mouse.ID),
                 min.age = min(Age),
                 max.age = max(Age),
                 mean.DEXA = round(mean(Age), digits = 0),
                 sd.DEXA = round(sd(Age), digits = 0))


write.table(table.2, file="Results/Table2.csv", sep=",", row.names=F)
## Table 2 END ##

## Figure 2: 3D scatter plot of body weight and lean body mass ####

data.3D <- read.csv("Data/3dcorrelation.csv")

scatterplot3d(data.3D, xlim=c(15,45), ylim=c(15,45), zlim=c(15,45), xlab="Lean by DEXA, g", zlab="Body weight, g", ylab="Lean by MR, g", cex.lab = 2, cex.axis = 1.2)
# savePlot(filename ="Figure_2_3d.tif",type ="tiff", device = dev.cur())
## Figure 2 END ##

## Figure 3: The genomic location of mouse QTL Burly1 identified in multiple mapping populations ####
# Copy data
figure3.dat <- Data_burly1[,c(1:142)] # cut the birth and endpoint dates
colnames(figure3.dat) <- tolower(names(figure3.dat))
names(figure3.dat)

## F2: Bodyweight plus Age covariate ####
## Part 1: -logp values across Chromosome 2

aov.out <-
  tidy(aov(lm(
      bodyweight ~ figure3.dat[, 23] + age,
      data = figure3.dat,
      subset = figure3.dat$table1_mapping_population == 'F2'
    )
  ))
aov.out["SNP"] <- colnames(figure3.dat)[23]
aov.out["Strain"] <- "F2"

for(i in 24:length(names(figure3.dat))) {
  temp.f2 <-
    tidy(aov(lm(
        bodyweight ~ figure3.dat[, i] + age,
        data = figure3.dat,
        subset = table1_mapping_population == 'F2'
      )
    ))
  temp.f2["SNP"] <- colnames(figure3.dat)[i]
  temp.f2["Strain"] <- "F2"
  aov.out <- rbind(aov.out, temp.f2)
}

# Only take main effect of genotype, which is listed in output as the looped term
f2.results <- aov.out[grep("figure3.dat", aov.out$term), ]

f2.results %<>% 
  .[, c(8, 7, 2:6)] %>% # reorder
  merge(., snp, by.x = "SNP", by.y = "Marker", all = T) %>% # merge with SNP info
  arrange(., Mark.Position.Mb)

# Calculate negative logp
f2.results$neg.logp <- (-log10(f2.results$p.value))

# WRITE OUT for part one of graphpad figure 3
# Mark.Position and neg.logp are the relevant columns
write.table(f2.results, file="Results/Figure3_F2_association.csv", sep=",", row.names=F)

## Part 2: Bar Chart
# which SNP is the top hit?
as.character(f2.results[which.max(f2.results$neg.logp), ][1]) # rs3023694

# get the n, body weight, and standard deviation for top snp
raw.table <-
  figure3.dat %>% 
  filter(., table1_mapping_population2 == "F2") %>%
  ddply(., ~rs3023694, summarise,
        n = length(mouse.id),
        BW.mean = mean(bodyweight),
        BW.sd = sd(bodyweight),
        BW.se = sem(bodyweight),
        Bw.max = max(bodyweight),
        Bw.min = min(bodyweight))


## Least Square Means
# We do this to obtain "corrected" (covariate-adjusted) means and standard deviations for our DV, so that we can make meaningful post-hoc comparisons

# Run manually on the top hit
manual.fit <- aov(lm(bodyweight ~ rs3023694 + age,
                     data = figure3.dat,
                     subset = figure3.dat$table1_mapping_population == 'F2'))

# Partial eta squared for the model
etaSquared(manual.fit, type=1)

# Run lsmeans to generate new fit object
leastsquare.fit <- lsmeans(manual.fit,
                        pairwise ~ rs3023694,
                        data = figure3.dat,
                        subset = figure3.dat$table1_mapping_population == 'F2')

# Compact Letter Display of pairwise comparisons
ls.table <- cld(leastsquare.fit,
                alpha=.05, 
                Letters=letters)

# Combine with raw
full.table <- merge(raw.table, ls.table)

# Calculate least squares SD from SE
full.table$ls.sd <- full.table$SE * sqrt(full.table$n)

# Cut down to relevant; This is graphpad part 2
cohen <- full.table[, c(1, 2, 8, 14)]
write.table(cohen, file="Results/Figure3_F2_bar.csv", sep=",", row.names=F)

# Exclude genotype B (129/129 in this case)
cohen <- cohen[which(full.table$rs3023694 != "B"), ]

# Obtain effects size; this is in column D of graphpad (part 3)
effectsize <- numeric(0)

effectsize[1] <- 
((cohen$lsmean[1] - cohen$lsmean[2])
/sqrt(
  (cohen$ls.sd[1]^2 + cohen$ls.sd[2]^2)
  /2
))

effectsize[1] # 0.2664629

## F2_Second: Select set of genotypes; Lean + Total ####
## Part 1: -logp values across Chromosome 2

aov.out <- tidy(aov(lm(lean ~ figure3.dat[,23]+total, data=figure3.dat, subset=table1_mapping_population=='F2_second')))
aov.out["SNP"] <- colnames(figure3.dat)[23]
aov.out["Strain"] <- "F2_second"


for (i in c(34:52, 72:81, 115:133)){  
  temp.f2.hf<-tidy(aov(lm(lean~figure3.dat[,i] + total, data=figure3.dat, subset=table1_mapping_population=='F2_second')))
  temp.f2.hf["SNP"] <- colnames(figure3.dat)[i]
  temp.f2.hf["Strain"] <- "F2_second"
  aov.out <- rbind(aov.out, temp.f2.hf)
}


# Main Effect
f2.hf.results <- aov.out[grep("figure3.dat", aov.out$term), ]

f2.hf.results %<>% 
  .[, c(8, 7, 2:6)] %>% # reorder
  merge(., snp, by.x = "SNP", by.y = "Marker", all = T) %>% # merge with SNP info
  arrange(., Mark.Position.Mb)

# Negative logp for graphpad
f2.hf.results$neg.logp <- (-log10(f2.hf.results$p.value))

# WRITE OUT for part one of graphpad figure 3
# Mark.Position and neg.logp are the relevant columns
write.table(f2.hf.results, file="Results/Figure3_F2_second_association.csv", sep=",", row.names=F)

# Part 2: Bar Chart
# which SNP is the top hit?
as.character(f2.hf.results[which.max(f2.hf.results$neg.logp), ][1]) # d2mit285

# get the n, body weight, and standard deviation for top snp
raw.table <-
  figure3.dat %>% 
  filter(., table1_mapping_population == "F2_second") %>%
  ddply(., ~d2mit285, summarise,
        n = length(mouse.id),
        lean.mean = mean(lean),
        lean.sd = sd(lean),
        lean.se = sem(lean),
        lean.max = max(lean),
        lean.min = min(lean))


## Least Square Means
# run manually on the top hit
manual.fit <- aov(lm(
      lean ~ d2mit285 + total,
      data = figure3.dat,
      subset = figure3.dat$table1_mapping_population == 'F2_second'))

# Partial eta squared for the model
etaSquared(manual.fit, type=1)

# lsmeans to generate new fit object
leastsquare.fit <- lsmeans(manual.fit,
                           pairwise ~ d2mit285,
                           # adjust=,
                           data = figure3.dat,
                           subset = figure3.dat$table1_mapping_population == 'F2_second')

# Pairwise comparisons
ls.table <- cld(leastsquare.fit,
                alpha=.05, 
                Letters=letters)

# Combine w. raw
full.table <- merge(raw.table, ls.table)

# Calculate least squares SD from SE
full.table$ls.sd <- full.table$SE * sqrt(full.table$n)

# Graphpad part 2
cohen <- full.table[, c(1, 2, 8, 14)]
write.table(cohen, file="Results/Figure3_F2_second_bar.csv", sep=",", row.names=F)

# Exclude genotype B (129/129 in this case)
cohen <- cohen[which(full.table$d2mit285 != "B"), ]

# Obtain effects size; this is in column D of graphpad (part 3)
effectsize[2] <- 
  ((cohen$lsmean[1] - cohen$lsmean[2])
   /sqrt(
     (cohen$ls.sd[1]^2 + cohen$ls.sd[2]^2)
     /2
   ))

effectsize[2] # 0.5176211


## Backcross_129: Lean w. Total covariate ####
## Part 1: -logp values across Chromosome 2

aov.out <- tidy(aov(lm(lean ~ figure3.dat[,22]+total, data=figure3.dat, subset=table1_mapping_population2=='Backcross_129')))
aov.out["SNP"] <- colnames(figure3.dat)[22]
aov.out["Strain"] <- "Backcross_129"

for (i in 23:(length(names(figure3.dat)))-1){ 
  temp.Backcross_129<-tidy(aov(lm(lean~figure3.dat[,i] + total, data=figure3.dat, subset=table1_mapping_population2=='Backcross_129')))
  temp.Backcross_129["SNP"] <- colnames(figure3.dat)[i]
  temp.Backcross_129["Strain"] <- "Backcross_129"
  aov.out <- rbind(aov.out, temp.Backcross_129)
}

# Main effect of genotype
backcross.129.results <- aov.out[grep("figure3.dat", aov.out$term), ]

backcross.129.results %<>% 
  .[, c(8, 7, 2:6)] %>% # reorder
  merge(., snp, by.x = "SNP", by.y = "Marker", all = T) %>% # merge with SNP info
  arrange(., Mark.Position.Mb)

# Negative logp
backcross.129.results$neg.logp <- (-log10(backcross.129.results$p.value))

# WRITE OUT for part one of graphpad figure 3
# Mark.Position and neg.logp are the relevant columns
write.table(backcross.129.results, file="Results/Figure3_Backcross_129_association.csv", sep=",", row.names=F)

# Part 2: Bar Chart
# which SNP is the top hit?
as.character(backcross.129.results[which.max(backcross.129.results$neg.logp), ][1]) # rs3687512

# get the n, body weight, and standard deviation for top snp
raw.table <-
  figure3.dat %>% 
  filter(., table1_mapping_population2 == "Backcross_129") %>%
  ddply(., ~rs3687512, summarise,
        n = length(mouse.id),
        BW.mean = mean(lean),
        BW.sd = sd(lean),
        BW.se = sem(lean),
        Bw.max = max(lean),
        Bw.min = min(lean))


## LS Means
# run manually on the top hit
manual.fit <- aov(lm(lean ~ rs3687512 + total,
                     data = figure3.dat,
                     subset = figure3.dat$table1_mapping_population2 == 'Backcross_129'))

# Partial eta squared for the model
etaSquared(manual.fit, type=1)

# lsmeans to generate new fit object
leastsquare.fit <- lsmeans(manual.fit,
                           pairwise ~ rs3687512,
                           data = figure3.dat,
                           subset = figure3.dat$table1_mapping_population == 'Backcross_129')

# Pairwise comparisons
ls.table <- cld(leastsquare.fit,
                alpha=.05, 
                Letters=letters)

# Combine with raw
full.table <- merge(raw.table, ls.table)

# Calculate least squares SD from SE
full.table$ls.sd <- full.table$SE * sqrt(full.table$n)

# Graphpad part 2
cohen <- full.table[, c(1, 2, 8, 14)]
write.table(cohen, file="Results/Figure3_Backcross_129_bar.csv", sep=",", row.names=F)

# Obtain effects size; this is in column D of graphpad (part 3)

effectsize[3] <- 
  ((cohen$lsmean[1] - cohen$lsmean[2])
   /sqrt(
     (cohen$ls.sd[1]^2 + cohen$ls.sd[2]^2)
     /2
   ))

effectsize[3] # -0.3202797


## Backcross_B6: Lean w. Total ####
## Part 1: -logp values across Chromosome 2

aov.out <- tidy(aov(lm(lean ~ figure3.dat[,22]+total, data=figure3.dat, subset=table1_mapping_population2=='Backcross_B6')))
aov.out["SNP"] <- colnames(figure3.dat)[22]
aov.out["Strain"] <- "Backcross_B6"

for (i in 23:(length(names(figure3.dat)))-1){ 
  temp.Backcross_B6<-tidy(aov(lm(lean~figure3.dat[,i] + total, data=figure3.dat, subset=table1_mapping_population2=='Backcross_B6')))
  temp.Backcross_B6["SNP"] <- colnames(figure3.dat)[i]
  temp.Backcross_B6["Strain"] <- "Backcross_B6"
  aov.out <- rbind(aov.out, temp.Backcross_B6)
}

# Main effect of genotype
backcross.b6.results <- aov.out[grep("figure3.dat", aov.out$term), ]

backcross.b6.results %<>% 
  .[, c(8, 7, 2:6)] %>% # reorder
  merge(., snp, by.x = "SNP", by.y = "Marker", all = T) %>% # merge with SNP info
  arrange(., Mark.Position.Mb)

# Negative logp
backcross.b6.results$neg.logp <- (-log10(backcross.b6.results$p.value))

# WRITE OUT for part one of graphpad figure 3
# Mark.Position and neg.logp are the relevant columns
write.table(backcross.b6.results, file="Results/Figure3_Backcross_B6_association.csv", sep=",", row.names=F)

# Part 2: Bar Chart
# which SNP is the top hit?
as.character(backcross.b6.results[which.max(backcross.b6.results$neg.logp), ][1]) # rs27350529, but several are same

# get the n, body weight, and standard deviation for top snp
raw.table <-
  figure3.dat %>% 
  filter(., table1_mapping_population2 == "Backcross_B6") %>%
  ddply(., ~rs27350529, summarise,
        n = length(mouse.id),
        BW.mean = mean(lean),
        BW.sd = sd(lean),
        BW.se = sem(lean),
        Bw.max = max(lean),
        Bw.min = min(lean))


## LS Means

# run manually on the top hit
manual.fit <- aov(lm(lean ~ rs27350529 + total,
                     data = figure3.dat,
                     subset = figure3.dat$table1_mapping_population2 == 'Backcross_B6'))

# Partial eta squared for the model
etaSquared(manual.fit, type=1)

# lsmeans to generate new fit object
leastsquare.fit <- lsmeans(manual.fit,
                           pairwise ~ rs27350529,
                           data = figure3.dat,
                           subset = figure3.dat$table1_mapping_population2 == 'Backcross_B6')

# Pairwise comparisons
ls.table <- cld(leastsquare.fit,
                alpha=.05, 
                Letters=letters)

# combine
full.table <- merge(raw.table, ls.table)

# Calculate least squares SD from SE
full.table$ls.sd <- full.table$SE * sqrt(full.table$n)

# Graphpad part 2
cohen <- full.table[, c(1, 2, 8, 14)]
write.table(cohen, file="Results/Figure3_Backcross_B6_bar.csv", sep=",", row.names=F)

# Obtain effects size; this is in column D of graphpad (part 3)

effectsize[4] <- 
  ((cohen$lsmean[1] - cohen$lsmean[2])
   /sqrt(
     (cohen$ls.sd[1]^2 + cohen$ls.sd[2]^2)
     /2
   ))

effectsize[4] # -0.3202797

## Congenics: Log transformed Lean + Total, subset SNPs ####
## Part 1: -logp values across Chromosome 2
aov.out <- tidy(aov(lm(log(lean) ~ figure3.dat[,60]+total, data=figure3.dat, subset=table1_mapping_population2=='Congenic')))
aov.out["SNP"] <- colnames(figure3.dat)[60]
aov.out["Strain"] <- "Congenic"

for (i in 61:length(names(figure3.dat))){ 
  temp.Congenic<-tidy(aov(lm(log(lean)~figure3.dat[,i] + total, data=figure3.dat, subset=table1_mapping_population2=='Congenic')))
  temp.Congenic["SNP"] <- colnames(figure3.dat)[i]
  temp.Congenic["Strain"] <- "Congenic"
  aov.out <- rbind(aov.out, temp.Congenic)
}

# Main effect of genotype
congenic.results <- aov.out[grep("figure3.dat", aov.out$term), ]

congenic.results %<>% 
  .[, c(8, 7, 2:6)] %>% # reorder
  merge(., snp, by.x = "SNP", by.y = "Marker", all = T) %>% # merge with SNP info
  arrange(., Mark.Position.Mb)

# Negative logp
congenic.results$neg.logp <- (-log10(congenic.results$p.value))

# WRITE OUT for part one of graphpad figure 3
# Mark.Position and neg.logp are the relevant columns
write.table(congenic.results, file="Results/Figure3_Congenic_association.csv", sep=",", row.names=F)

# Part 2: Bar Chart
# which SNP is the top hit?
as.character(congenic.results[which.max(congenic.results$neg.logp), ][1]) # rs27346400, but several are same

# get the n, body weight, and standard deviation for top snp
raw.table <-
  figure3.dat %>% 
  filter(., table1_mapping_population2 == "Congenic") %>%
  ddply(., ~rs27346400, summarise,
        n = length(mouse.id),
        BW.mean = mean(lean),
        BW.sd = sd(lean),
        BW.se = sem(lean),
        Bw.max = max(lean),
        Bw.min = min(lean))

## LS Means
# run manually on the top hit
manual.fit <- aov(lm(lean ~ rs27346400 + total,
                     data = figure3.dat,
                     subset = figure3.dat$table1_mapping_population2 == 'Congenic'))

# Partial eta squared for the model
etaSquared(manual.fit, type=1)

# lsmeans to generate new fit object
leastsquare.fit <- lsmeans(manual.fit,
                           pairwise ~ rs27346400,
                           data = figure3.dat,
                           subset = figure3.dat$table1_mapping_population2 == 'Congenic')

# Pairwise comparisons
ls.table <- cld(leastsquare.fit,
                alpha=.05, 
                Letters=letters)

# Combine
full.table <- merge(raw.table, ls.table)

# Calculate least squares SD from SE
full.table$ls.sd <- full.table$SE * sqrt(full.table$n)

# Graphpad part 2
cohen <- full.table[, c(1, 2, 8, 14)]
write.table(cohen, file="Results/Figure3_Congenic_bar.csv", sep=",", row.names=F)

# obtain effects size; this is in column D of graphpad (part 3)

effectsize[5] <- 
  ((cohen$lsmean[1] - cohen$lsmean[2])
   /sqrt(
     (cohen$ls.sd[1]^2 + cohen$ls.sd[2]^2)
     /2
   ))

effectsize[5] # 0.6855266

## Part 3: Save effect sizes of all analyses
strains <- c("F2_First", "F2_Second", "Backcross_129", "Backcross_B6", "Congenic")
effect.sizes <- data.frame("Strains" = strains,
                           "Effect Size " = effectsize)

write.table(effect.sizes, file="Results/Figure3_Effect_Sizes.csv", sep=",", row.names=F)

## Figure 3 END ##

## Figure 4: Reciprocal Consomics ####
## Part 1: -logp for 129.B6_chr2 across Chr2
aov.out <- tidy(aov(lm(lean ~ figure3.dat[,22]+sex+total, data=figure3.dat, subset=table1_mapping_population=='129.B6-Chr2')))
aov.out["SNP"] <- colnames(figure3.dat)[22]
aov.out["Strain"] <- "129.B6_Chr2"
  
for (i in 23:length(names(figure3.dat))){ 
  temp.129.B6_Chr2<-tidy(aov(lm(lean~figure3.dat[,i]+sex+total, data=figure3.dat, subset=table1_mapping_population=='129.B6-Chr2')))
  temp.129.B6_Chr2["SNP"] <- colnames(figure3.dat)[i]
  temp.129.B6_Chr2["Strain"] <- "129.B6-Chr2"
  aov.out <- rbind(aov.out, temp.129.B6_Chr2)
}

# Main Effect
results.129.B6_Chr2 <- aov.out[grep("figure3.dat", aov.out$term), ]

results.129.B6_Chr2 %<>% 
  .[, c(8, 7, 2:6)] %>% # reorder
  merge(., snp, by.x = "SNP", by.y = "Marker", all = T) %>% # merge with SNP info
  arrange(., Mark.Position.Mb)

# Negative logp for graphpad
results.129.B6_Chr2$neg.logp <- (-log10(results.129.B6_Chr2$p.value))

# WRITE OUT for Graphpad figure 4
write.table(results.129.B6_Chr2, file="Results/Figure4_129.B6-Chr2_association.csv", sep=",", row.names=F)

# For Bar Charts: Which SNP is the top hit?
as.character(results.129.B6_Chr2[which.max(results.129.B6_Chr2$neg.logp), ][1]) # "rs3681694"

## Figure 4 END##

## Figure 5 & Figure S3: Lean ####
## Part 1: -logp values across Chromosome 2

aov.out <- tidy(aov(lm(lean ~ figure3.dat[,60]+total, data=figure3.dat, subset=table1_mapping_population2=='Congenic')))
aov.out["SNP"] <- colnames(figure3.dat)[60]
aov.out["Strain"] <- "Congenic"

for (i in 61:length(names(figure3.dat))){ 
  temp.Congenic<-tidy(aov(lm(lean ~ figure3.dat[,i] + total, data=figure3.dat, subset=table1_mapping_population2=='Congenic')))
  temp.Congenic["SNP"] <- colnames(figure3.dat)[i]
  temp.Congenic["Strain"] <- "Congenic"
  aov.out <- rbind(aov.out, temp.Congenic)
}

congenic.lean.results <- aov.out[grep("figure3.dat", aov.out$term), ]

congenic.lean.results %<>% 
  .[, c(8, 7, 2:6)] %>% # reorder
  merge(., snp, by.x = "SNP", by.y = "Marker", all = T) %>% # merge with SNP info
  arrange(., Mark.Position.Mb)

# calculate negative logp for graphpad
congenic.lean.results$neg.logp <- (-log10(congenic.lean.results$p.value))

write.table(congenic.lean.results, file="Results/Figure5_Congenic_lean_association.csv", sep=",", row.names=F)

## Figure 5/S3 END ##

## Figure S1: Correlations among three measures (BW/Lean by DEXA and MR) ####
Large <- read.csv("Data/corr_burly1.csv")
Large<- Large[-which(is.na(Large[,13])), ]
p1 <- ggplot(Large, aes(Lean_by_DEXA,  Lean_by_MR)) + geom_point(shape=1)+ geom_smooth(method=lm) + facet_grid(. ~ Strain) + labs (x="Lean by DEXA, g", y="Lean by MR, g")

Large1 <- read.csv("Data/corr_burly1_2.csv")
Large1<- Large1[-which(is.na(Large1[,11])), ]
p2 <- ggplot(Large1, aes(Lean_body_mass,Body_weight, color=method)) + geom_point(shape=1)+geom_smooth(method=lm) + facet_grid(. ~ Strain) + labs(x="Lean body mass, g", y="Body weight, g")

ggarrange(ggarrange(p1,  ncol = 2, labels = c(""), widths=c(1, 0.9)), p2, nrow = 2, labels = c("a","b") )

## Figure S4: No effect on Fat ####
# Part 1: -logp across Chr2

aov.out <- tidy(aov(lm(fat ~ figure3.dat[,60]+total, data=figure3.dat, subset=table1_mapping_population2=='Congenic')))
aov.out["SNP"] <- colnames(figure3.dat)[60]
aov.out["Strain"] <- "Congenic"

for (i in 61:length(names(figure3.dat))){
  temp.Congenic<-tidy(aov(lm(fat~figure3.dat[,i] + total, data=figure3.dat, subset=table1_mapping_population2=='Congenic')))
  temp.Congenic["SNP"] <- colnames(figure3.dat)[i]
  temp.Congenic["Strain"] <- "Congenic"
  aov.out <- rbind(aov.out, temp.Congenic)
}

congenic.fat.results <- aov.out[grep("figure3.dat", aov.out$term), ]

congenic.fat.results %<>%
  .[, c(8, 7, 2:6)] %>% # reorder
  merge(., snp, by.x = "SNP", by.y = "Marker", all = T) %>% # merge with SNP info
  arrange(., Mark.Position.Mb)

# Negative logp for graphpad
congenic.fat.results$neg.logp <- (-log10(congenic.fat.results$p.value))

write.table(congenic.fat.results, file="Results/FigureS4_Congenic_fat_association.csv", sep=",", row.names=F)




## S4 Table correlation within each mapping population (CJ NOT checked yet) ####
# 
# DEXA=subset(Large, method=="DEXA")
# # DEXA <- Large
# a=subset(DEXA, Strain=="129.B6-Chr2")
# cor.test(a$Lean_body_mass,a$Body_weight, method="pearson" )
# # cor.test(a$Lean_by_DEXA,a$Body.weight..g, method="pearson" )
# 
# b=subset(DEXA, Strain=="B6.129-Burly1")
# cor.test(b$Lean_body_mass,b$Body_weight, method="pearson" )
# 
# 
# c=subset(DEXA, Strain=="B6.129-Chr2")
# cor.test(c$Lean_body_mass,c$Body_weight, method="pearson" )
# 
# d=subset(DEXA, Strain=="C57BL/6ByJ")
# cor.test(d$Lean_body_mass,d$Body_weight, method="pearson" )
# 
# e=subset(DEXA, Strain=="F1 x 129")
# cor.test(e$Lean_body_mass,e$Body_weight, method="pearson" )
# 
# f=subset(DEXA, Strain=="F1 x B6")
# cor.test(f$Lean_body_mass,f$Body_weight, method="pearson" )
# 
# g=subset(DEXA, Strain=="F2_Second")
# cor.test(g$Lean_body_mass,g$Body_weight, method="pearson" )
# 
# 

## Supplementary Figure 2: Monthly lean body mass in male mice from two populations: mice with B6background (a; n=319) and mice with 129 background
# where is this?

## Supplementary Table 6: Congenics strains included and excluded from analyses? ####
# NEED TO MAKE SURE TO EXCLUDE ALL MICE ALREADY EXCLUDED ELSEWHERE
# then, substitute filtered dataset into below call
table(Data_burly1$Table_2_Substrains, Data_burly1$HasFragment)
# table(s7.data$Table_2_Substrains, s7.data$HasFragment)

## Supplementary Table 7: Post Hoc comparision between mice w. and without donor region within congenic strains ####
# Recreated to match Statistica done by Cailu

# Substrains to exclude, plus NA
# From Statistica output:
# Exclude condition: v3="ss.1.1" OR v3="ss.1.11" OR v3="ss.1.14" OR v3="ss1.6" OR v3="ss.2.2" OR v3="ss.1.3" OR v3="ss.1.9"
exclude <- c("1.1", "1.11", "1.14", "1.6", "2.2", "1.3", "1.9", NA)

# Subset data
s7.data <- 
  Data_burly1 %>% 
  .[which(.$Table_2_Substrains %not% exclude), ] %>% 
  .[which(.$Table1_Mapping_population2 != "Congenic_not"), ]

# Make factor for Type I ANOVA
s7.data$HasFragment <- as.factor(s7.data$HasFragment)

table(s7.data$Table_2_Substrains, s7.data$HasFragment)

# ANOVA on Lean
anov.obj <- aov(s7.data$Lean ~ s7.data$Table_2_Substrains +
                  s7.data$HasFragment +
                  s7.data$Table_2_Substrains:s7.data$HasFragment)

summary(anov.obj) # matches Statistica
#                                                   Df Sum Sq Mean Sq F value   Pr(>F)    
# s7.data$Table_2_Substrains                        17    648   38.09   9.104  < 2e-16 ***
#   s7.data$HasFragment                              1    303  302.61  72.326  < 2e-16 ***
#   s7.data$Table_2_Substrains:s7.data$HasFragment  17    216   12.71   3.038 3.33e-05 ***
#   Residuals                                      920   3849    4.18   

## Post-Hoc Pairwise Table: LSD
# Take ANOVA object, calculate LSD pairwise comparisons, and turn into Data Table
posthocs <- 
  anov.obj %>% 
  PostHocTest(., method = "lsd") %>% 
  .$`s7.data$Table_2_Substrains:s7.data$HasFragment` %>% 
  data.frame(.) %>% 
  setDT(., keep.rownames = T)

write.table(posthocs, file="Results/TableS7_Congenic_LSD_posthoc.csv", sep=",", row.names=F)

## Table S7 END ##

## Labmaster: Data Processing ####
## This script was drafted by Cailu Lin in 2016_08 for processing the labmaster data 

## create statistical table
Table_For_analyses<-file("Results/table_processed.csv", "w")
cat("Animal.No.", "Average.RER", "Hr.H.2", "Hr.XT.YT.Z", "Dy.Feed", "Dy.Drink", "Hr.VO2.2", "Hr.VCO2.2", file=Table_For_analyses, "\n" ,sep=", " )


id <- as.character(read.csv(file = "Data/labmaster.id.csv", header = TRUE)[,1])

for (i in 1:length(id))
{
  #### process labmaster data aninaml by animal
  requested_animal_no <- id[i] ##### type animal ID #####
  
  # Load records of animal measurements.
  records <- read.csv(file="Data/labmaster.records.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
  
  ##########################################################
  
  # Extract rows for the animal "requested_animal_no".
  animal_records <- records[records$Animal.No == requested_animal_no,]
  
  # Sort records for the animal on the "Date" and "Time" columns.
  sorted_animal_records <- animal_records[order(animal_records$Date, animal_records$Time),]
  
  # Slice to only contain the "Date" and "Time" columns.
  sorted_animal_date_times <- sorted_animal_records[,c("Date", "Time")]
  
  # Print the sorted (Date, Time) pairs for the animal.
  sorted_animal_date_times
  
  ##########################################################
  # Use procedural statements to find the minimum and maximum
  # date_time string in rows where "Animal.No." is
  # "requested_animal_no".
  
  min_date_time_string <- NA
  max_date_time_string <- NA
  for (row in 1 : dim(animal_records)[1])
  {
    
    # Skip this row if either "Animal.No" is not available (NA) or
    # the row is not for the requested animal.
    if (!is.na(animal_records[row, "Animal.No."]) &&
        animal_records[row, "Animal.No."] == requested_animal_no)
    {
      # Extract the "Date" and "Time" value from this row and
      # store as a single string.
      date_time <- animal_records[row, c("Date", "Time")]
      date_time_string <- paste(date_time[1], date_time[2])
      
      # Update the minimum and maximum date-time string seen so far.
      if (is.na(min_date_time_string) || date_time_string < min_date_time_string)
      {
        min_date_time_string = date_time_string
      }
      else if (is.na(max_date_time_string) || date_time_string > max_date_time_string)
      {
        max_date_time_string = date_time_string
      }
    }
  }
  
  # Print the minimum and maximum date time strings that were found for the
  # requested animal.
  # print "minimum date time string:"
  min_date_time_string
  #print "maximum date time string:"
  max_date_time_string
  
  #### function to extract <hour> from a string in the format of "hour:minute...".
  
  extract_hour <- function( hour_time_string )
  {
    as.numeric( strsplit( hour_time_string, ":" )[[1]][[1]] )
  }
  
  ##### function to determine if a given row from records is empty.
  # The row is considered empty when the Date value is missing or is empty.
  
  is_empty_row <- function( row_data )
  {
    is.na( row_data[ "Date" ] ) || row_data[ "Date" ] == ""
  }
  
  # Load records of measurements for a given animal.
  
  requested_animal_no <- id[i]
  
  all_records <- read.csv( file="Data/labmaster.records.csv", header=TRUE, sep=",", stringsAsFactors=FALSE )
  records <- all_records[ all_records$"Animal.No." == requested_animal_no, ]
  
  ##### Update rows in "records" with statistics.
  
  first_data_row <- 2
  last_data_row <- dim( records )[ 1 ]
  first_current_hour_row <- NA
  rows_with_hourly_statistics <- NULL
  
  for( row in first_data_row : ( last_data_row - 1 ) )
  {
    # Skip empty rows.
    
    if( is_empty_row( records[ row, ] ) )
    {
      next
    }
    
    if( is_empty_row( records[ row + 1, ] ) )
    {
      next
    }
    
    ##### Compute statistics over the next and the current row.
    
    current_drink_value <- as.numeric( records[ row, "Drink1" ] )
    next_drink_value <- as.numeric( records[ row + 1, "Drink1" ] )
    records[ row, "Diff.Drink" ] <- next_drink_value - current_drink_value
    
    current_feed_value <- as.numeric( records[ row, "Feed1" ] )
    next_feed_value <- as.numeric( records[ row + 1, "Feed1" ] )
    records[ row, "Diff.Feed" ] <- next_feed_value - current_feed_value
    
    # Compute statistics over the hour for the current row.
    #
    # This computation is performed over a number of rows that have the same hour.
    # The computation is triggered when the next row has a different hour.
    
    if( is.na( first_current_hour_row ) )
    {
      first_current_hour_row = row;
    }
    
    current_hour <- extract_hour( records[ row, "Time" ] )
    next_hour <- extract_hour( records[ row + 1, "Time" ] )
    
    if( current_hour != next_hour )
    {
      last_current_hour_row = row;
      
      # Compute statistics over the first_current_hour_row : last_current_hour_row rows
      # (the rows for the current hour) and record with the first_current_hour_row row.
      
      records[ first_current_hour_row, "Average.RER" ] <-
        mean( as.numeric( records[ first_current_hour_row : last_current_hour_row, "RER" ] ) )
      records[ first_current_hour_row, "Hr.H.2." ] <-
        mean( as.numeric( records[ first_current_hour_row : last_current_hour_row, "H.2." ] ) )
      records[ first_current_hour_row, "Hr.XT.YT.Z" ] <-
        2*(mean( as.numeric( records[ first_current_hour_row : last_current_hour_row, "XT.YT.Z" ] ) ))
      records[ first_current_hour_row, "Hr.Drink" ] <-
        2*(mean( as.numeric( records[ first_current_hour_row : last_current_hour_row, "Diff.Drink" ] )))
      records[ first_current_hour_row, "Hr.Feed" ] <- 
        2*(mean( as.numeric( records[ first_current_hour_row : last_current_hour_row, "Diff.Feed" ] )))
      records[ first_current_hour_row, "Hr.VO2.2" ] <-
        mean( as.numeric( records[ first_current_hour_row : last_current_hour_row, "VO2.2." ] ) )
      records[ first_current_hour_row, "Hr.VCO2.2" ] <-
        mean( as.numeric( records[ first_current_hour_row : last_current_hour_row, "VCO2.2." ] ) )
      records[ first_current_hour_row, "Hour" ] <- current_hour
      
      rows_with_hourly_statistics <- c( rows_with_hourly_statistics, first_current_hour_row )
      
      # Start new hour block.
      
      first_current_hour_row = row + 1;
    }
  }
  
  # Plot hourly statistics.
  
  plot( records[ rows_with_hourly_statistics, "Hour" ], records[ rows_with_hourly_statistics, "Average.RER" ] )
  plot( records[ rows_with_hourly_statistics, "Hour" ], records[ rows_with_hourly_statistics, "Hr.H.2." ] )
  plot( records[ rows_with_hourly_statistics, "Hour" ], records[ rows_with_hourly_statistics, "Hr.XT.YT.Z" ] )
  plot( records[ rows_with_hourly_statistics, "Hour" ], records[ rows_with_hourly_statistics, "Hr.Feed" ] )
  plot( records[ rows_with_hourly_statistics, "Hour" ], records[ rows_with_hourly_statistics, "Hr.Drink" ] )
  plot( records[ rows_with_hourly_statistics, "Hour" ], records[ rows_with_hourly_statistics, "Hr.VO2.2" ] )
  plot( records[ rows_with_hourly_statistics, "Hour" ], records[ rows_with_hourly_statistics, "Hr.VCO2.2" ] )
  
  a = mean(records[ rows_with_hourly_statistics, "Average.RER" ])
  b = mean(records[ rows_with_hourly_statistics, "Hr.H.2." ])
  c = mean(records[ rows_with_hourly_statistics, "Hr.XT.YT.Z" ] )
  d = 24*(mean(records[ rows_with_hourly_statistics, "Hr.Feed" ]))
  e = 24*(mean(records[ rows_with_hourly_statistics, "Hr.Drink"]))
  f = mean(records[ rows_with_hourly_statistics, "Hr.VO2.2" ]) 
  g = mean(records[ rows_with_hourly_statistics, "Hr.VCO2.2" ])
  
  cat(requested_animal_no, a, b, c, d, e, f, g, file=Table_For_analyses, "\n" ,sep=", " )
}


## Labmaster: T-tests on processed data (CJ not tweaked yet) ####

geno <- read.csv("Data/labmaster.id.csv", header = TRUE)
metab <- read.csv ("Results/table_processed.csv", header = TRUE)
table <- merge(geno, metab, by="Animal.No.")

##t.test for Average.RER
t.test(table$Average.RER ~table$geno)

##t.test for Hr.XT.YT.Z
t.test(table$Hr.XT.YT.Z ~table$geno)

##t.test for Hr.VO2.2
t.test(table$Hr.VO2.2 ~table$geno)

##t.test for Hr.VCO2.2
t.test(table$Hr.VCO2.2 ~table$geno)

##t.test for Hr.H.2
t.test(table$Hr.H.2 ~table$geno)


##t.test for Dy.Feed
t.test(((table$Dy.Feed*1000)/table$lean) ~table$geno)

##t.test for Dy.Drink
t.test(((table$Dy.Drink*1000)/table$lean) ~table$geno)



## Burly1 script END ##
setwd("S:/reed/Labshare/CJ Arayata/Mennella/Genotyping/16s Saliva Bacteria/Check out FUT2 - May 2016")

#get melted data, get genotype data, merge together
bact_melt <- read.delim("S:/reed/Labshare/CJ Arayata/Mennella/Genotyping/16s Saliva Bacteria/Check out FUT2 - May 2016/2016_0112 - 16s Statistica Stacked.txt", header=TRUE)
View(bact_melt)

fut2_genotypes <- read.delim("S:/reed/Labshare/CJ Arayata/Mennella/Genotyping/16s Saliva Bacteria/Check out FUT2 - May 2016/2016_0517 - FUT2 genotypes.txt", header=TRUE)
View(fut2_genotypes)

bact_fut <- merge(bact_melt, fut2_genotypes, by.x="SubjectID", by.y="SubjectID", all.x=TRUE)

#subset based on study; also need to drop unused levels!!
mbf_bact_fut <- bact_fut[which(bact_fut$Study=="MBF"), ]
mbf_bact_fut <- droplevels(mbf_bact_fut)
gro_bact_fut <- bact_fut[which(bact_fut$Study=="GRO"), ]
gro_bact_fut <- droplevels(gro_bact_fut)


#cross-tabs
library(broom)
cross_tabs_mbf <- tidy(xtabs(~ rs601338 + FeedingType, mbf_bact_fut))
cross_tabs_gro <- tidy(xtabs(~ rs601338 + FeedingType, gro_bact_fut))


#summary statistics i guess
library(psych)
describe(bact_fut)
describe.by(bact_fut, group=bact_fut$Study)


#run a mixed ANOVA - study (MBF/GRO) x feeding type (???) x genotype x round (3 within)
anova_mbf <- aov(Baseline.Change ~ (FeedingType*rs601338*RoundNumber) + Error((SubjectID/RoundNumber)+(FeedingType*rs601338)), data=mbf_bact_fut)
summary(anova_mbf)
anova_gro <- aov(Baseline.Change ~ (FeedingType*rs601338*RoundNumber) + Error((SubjectID/RoundNumber)+(FeedingType*rs601338)), data=gro_bact_fut)
summary(anova_gro)

#plots - baseline change
interaction.plot(mbf_bact_fut$RoundNumber, mbf_bact_fut$rs601338, mbf_bact_fut$Baseline.Change, fun=mean, legend=TRUE)

interaction.plot(gro_bact_fut$RoundNumber, gro_bact_fut$rs601338, gro_bact_fut$Baseline.Change, fun=mean, legend=TRUE)

#plots - mean fold change for MBF - FUT2 genotype for breastfeeders
interaction.plot(mbf_bact_fut$RoundNumber, mbf_bact_fut$rs601338, mbf_bact_fut$Mean.Fold.Change, fun=mean, legend=TRUE)

#plots - mean fold change for GRO - feed group for formula feeders
interaction.plot(gro_bact_fut$RoundNumber, gro_bact_fut$FeedingType, gro_bact_fut$Mean.Fold.Change, fun=mean, legend=TRUE)

#what about nicer ggplots?
library(ggplot2)




#ideas for future poking - 16s, FUT2, and anthro data - would need to tie it all together and perhaps there is a story to be told?!?!
#bacteria is good? transfer via breastmilk and FUT2 (MBF); do they have better growth outcomes
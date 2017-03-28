setwd("//monell/storage/reed/Labshare/CJ Arayata/Cohen/Dec 2016 - Jen Gene Expression/send to dani")

# from raw access file, i added two columns for row/col of plate (eg. 'A1' turns to 'A' '1')
# turned 'Undetermined' cycle counts into NAs (ct_na) AND 40 (ct_40) to calculate both ways

###############
# STEP 2: START
###############

data <- read.csv("jen_exp.csv", header = T, na.strings = c("", NA))

# to lower
colnames(data) <- tolower(colnames(data))

# make dates actual dates
data$date <- as.Date(data$date, format = "%m/%d/%Y")

table(data$subjectid) #something is up, need to strip whitespace

data$subjectid <- as.character(data$subjectid)
data$subjectid <- trimws(data$subjectid) # fancy, its in base R now!
data$target.name <- trimws(data$target.name)
table(data$target.name) # check target names
data$target.name <- gsub("rs102\\b", "rs10246939", data$target.name) # replace exact string of "rs102"
  
# OKKKKK..... sum up data as summary of
## subjectid, plate, sample type, targetname, reporter, MEAN and MEDIAN
## need to run on 'ct_na' column as well as 'ct_40' column

library(plyr)

# for spot checking
fun.check <- data[data$subjectid %in% c('FUN-1'), ]

# here it is - within each plate
summary.perplate <- ddply(data, .(subjectid, plate, date, sample.type, target.name, reporter), summarize,
                 replicates = length(ct_na),
                 num.undetermined = sum(is.na(ct_na)),
                 mean.remove.na = mean(ct_na, na.rm=T),
                 mean.setas40 = mean(ct_40),
                 median.remove.na = median(ct_na, na.rm=T),
                 median.setas40 = median(ct_40))

# # collapse plates / dates
# summary.acrossplate <- ddply(data, .(subjectid, sample.type, target.name, reporter), summarize,
#                           replicates = length(ct_na),
#                           num.undetermined = sum(is.na(ct_na)),
#                           mean.remove.na = mean(ct_na, na.rm=T),
#                           mean.setas40 = mean(ct_40),
#                           median.remove.na = median(ct_na, na.rm=T),
#                           median.setas40 = median(ct_40))



# write.csv(summary.acrossplate, file = "mean_medians_acrossplate.csv")

# ROUTE ONE: for per-plate, need to cut down to most recent plate per person

## paste together subjectid - sample type - marker - reporter
summary.perplate$full <- paste(summary.perplate$subjectid, "-", summary.perplate$sample.type, "-", summary.perplate$target.name, "-", summary.perplate$reporter)

# move full to first position
summary.perplate <- summary.perplate[, c(13, 1:12)]

# sort by date - newest first
summary.perplate <- summary.perplate[order(summary.perplate$full, summary.perplate$date, summary.perplate$plate, decreasing = T), ]

# take first instance
best.avail <- summary.perplate[!duplicated(summary.perplate$full), ]

# then i will have one 'run' per person that should represent 'best' data

means.best <- ddply(best.avail, .(sample.type, target.name, reporter), summarise,
               samples = length(mean.setas40),
               num.undeter = sum(is.na(mean.remove.na)),
               mean.all.remove = mean(mean.remove.na, na.rm=T),
               mean.all.as40 = mean(mean.setas40),
               median.all.remove = median(median.remove.na, na.rm=T),
               median.all.as40 = median(median.setas40))

# now merge to have one row per person/target/reporter/etc
best.merge <- merge(best.avail, means.best)

# reorganize
best.merge <- best.merge[, c(5, 1:3, 6:ncol(best.merge))]

# export
write.csv(best.merge, file = "best_replicates.csv")

#########
# STEP 3: calculate delta Cts in EXCEL
#########

# IGNORE THIS
# # ROUTE TWO: collapse and summarise across plates
# ## one run per person that is mean of all replicates
# 
# # THEN:
# ## calculate mean of each target/reporter across ALL SAMPLES (but within each sample type?)
# ## can i do this with ddply into another table, and then MERGE so that appropriate means are duplicated?
# 
# # looking at this is helpful to figure out the unique combinations of sample/target/reporter
# means <- ddply(summary.acrossplate, .(sample.type, target.name, reporter), summarise,
#                samples = length(mean.setas40),
#                num.undeter = sum(is.na(mean.remove.na)),
#                mean.all.remove = mean(mean.remove.na, na.rm=T),
#                mean.all.as40 = mean(mean.setas40),
#                median.all.remove = median(median.remove.na, na.rm=T),
#                median.all.as40 = median(median.setas40))
# 
# # this works
# test.merge <- merge(summary.acrossplate, means)
# 
# # going to export to excel before i try to do anything
# write.csv(test.merge, file = "all_replicates.csv")
# # reorganize
# colnames(test.merge)
# test.merge <- test.merge[, c(4, 1:3, 5, 11, 6, 12, 7:10, 13:16)]

# THEN:
## formula to calculate delta-delta of each TAS2R38 marker against:
## GAPDH
## GNAT3 for taste, CFTR for nasal



## delta delta Ct:
## 2^[(meansampleCt housekeeping - meansampleCt target) - (meanallsamplesCt housekeeping - meanallsamplesCt target)]


#######
# STEP 4: trying R to give me TAS2R38 taster/nontaster stuff
#######

# read in sheet 3 from excel
stuff <- read.csv("tas2r38.csv", header = T, na.strings = c("", NA))

# checking - need to edit (should have done this earlier)
table(stuff$target.name)
stuff$target.name <- gsub("rs102\\b", "rs10246939", stuff$target.name)

# drop unused rows (ie. select TAS2R38 markers only)
tas2r38 <- stuff[stuff$target.name %in% c("rs10246939", "rs1726866", "rs713598"), ]

# using 'count' to find subjects who have 12 values, ie. they have data for all three markers
library(plyr)
bleh <- count(tas2r38$subjectid)

# merge the count data
tas2r38 <- merge(tas2r38, bleh, by.x = "subjectid", by.y = "x", all = T)

# select only those with 12
tas2r38 <- tas2r38[tas2r38$freq == 12, ]

# sort it
tas2r38.sort <- tas2r38[order(tas2r38$subjectid, tas2r38$sample.type, tas2r38$target.name), ]

# create composite variable target and reporter
tas2r38.sort$target.reporter <- paste(tas2r38.sort$target.name, "-", tas2r38.sort$reporter)

# create taster and nontaster means, medians, sums
# 713598 FAM, 1726866 FAM, 10246939 VIC - taster
# 713598 VIC, 1726866 VIC, 10246939 FAM - nontaster

table(tas2r38.sort$target.reporter)

tas2r38.sort$taster[tas2r38.sort$target.reporter == "rs713598 - FAM" |
                      tas2r38.sort$target.reporter == "rs1726866 - FAM" |
                      tas2r38.sort$target.reporter == "rs10246939 - VIC"] <- "taster"

tas2r38.sort$taster[tas2r38.sort$target.reporter == "rs713598 - VIC" |
                      tas2r38.sort$target.reporter == "rs1726866 - VIC" |
                      tas2r38.sort$target.reporter == "rs10246939 - FAM"] <- "nontaster"

tas2r38.table <- ddply(tas2r38.sort, .(subjectid, sample.type, taster), summarise,
                       mean.na = mean(mean.remove.na, na.rm = T),
                       mean.40 = mean(mean.setas40),
                       median.na = median(median.remove.na, na.rm = T),
                       median.40 = median(median.setas40),
                       sum.mean.na = sum(mean.remove.na, na.rm = T),
                       sum.mean.40 = sum(mean.setas40),
                       sum.median.na = sum(median.remove.na, na.rm = T),
                       sum.median.40 = sum(median.setas40))

tas2r38.table.means <- ddply(tas2r38.sort, .(sample.type, taster), summarise,
                    mean.all.remove = mean(mean.remove.na, na.rm=T),
                    mean.all.as40 = mean(mean.setas40),
                    median.all.remove = median(median.remove.na, na.rm=T),
                    median.all.as40 = median(median.setas40))

# merge and reorganize
tas2r38.all <- merge(tas2r38.table, tas2r38.table.means)
tas2r38.all <- tas2r38.all[, c(3, 1:2, 4:7, 12:15, 8:11)]

# export
write.csv(tas2r38.all, file = "taster.means.csv")




##################################################
# START HERE FOR PREPARING DATASET ###############
## really need to make a better file to start with
##################################################


#####
# STEP 5: excel to create deltas for these new summaries
deltas <- read.csv("deltas.csv", header = T, na.strings = c("", NA, "NA"))


#####
# 1/3/2017 - edits to access query; jen gave me DOBs so added into access first, re-exported, then re-churned here


#####
# STEP 6: access query!
## note: i tweaked a bit manually in excel (such as biopsy date and taste test date) to make further down steps easier
taste.dems <- read.csv("query.export3.csv", header = T, na.strings = c("", NA, "NA")) 

library(reshape2)
taste.cast <- dcast(taste.dems, ... ~ Marker, value.var = "Value") # ... is shorthand for all other variables

#######
# STEP 7: prepare analysis-ready dataset
prelim.dataset <- merge(deltas, taste.cast, by.x = "subjectid", by.y = "SubjectID", all = T)

# sort
prelim.dataset <- prelim.dataset[order(prelim.dataset$subjectid, prelim.dataset$sample.type, prelim.dataset$reporter), ]

## qualitative variable - expressed or not expressed (NA / 40 cycles)
prelim.dataset$expressed[prelim.dataset$mean.remove.na >= 0] <- "expressed"
prelim.dataset$expressed[is.na(prelim.dataset$mean.remove.na)] <- "not expressed"

# reorder
prelim.dataset <- prelim.dataset[, c(1, 31:34, 2:4, 72, 5, 7:30, 35:67, 70, 69, 68, 71)]

# # rename 'date' to 'date.taste' because that's what it actually represents
# names(prelim.dataset)[names(prelim.dataset)=='Date'] <- "Date.taste"

names(prelim.dataset)[names(prelim.dataset)=='Biposy.Time'] <- "Biopsy.Time"

# convert biopsy and taste testing dates to actual dates
prelim.dataset$Biopsy.Date <- as.Date(prelim.dataset$Biopsy.Date, "%m/%d/%Y")
prelim.dataset$Date.taste <- as.Date(prelim.dataset$Date.taste, "%m/%d/%Y")

# TIME: step 1 - convert to character
prelim.dataset$Biopsy.Time <- as.character(prelim.dataset$Biopsy.Time)

# TIME: step 1.5 - strptime function check - works
## adds today's date in arbitrarily but OK for our purposes
## %I is for 12-hr time, %p used in conjunction to denote AM/PM - result is 24hr time
time.test <- as.data.frame(strptime(prelim.dataset$Biopsy.Time, format = "%I:%M %p"))

# TIME: step 2 - convert in place
prelim.dataset$Biopsy.Time <- strptime(prelim.dataset$Biopsy.Time, format = "%I:%M %p")

# plot just to see whats up - time looks good!
library(ggplot2)
ggplot(prelim.dataset, aes(x = Biopsy.Time, y = mean.remove.na)) +
  geom_point()

# NOW i need to figure out how to calculate age; where is that super useful function i had??

# change class to date so i can calculate age
class(prelim.dataset$DOB)
prelim.dataset$DOB <- as.Date(prelim.dataset$DOB, format = "%m/%d/%Y")

## need to subtract subject DOB from swab date?

# define an age function for the control subjects
age <- function(from, to) {
  from_lt = as.POSIXlt(from)
  to_lt = as.POSIXlt(to)
  
  age = to_lt$year - from_lt$year
  
  ifelse(to_lt$mon < from_lt$mon |
           (to_lt$mon == from_lt$mon & to_lt$mday < from_lt$mday),
         age - 1, age)
}

# now calculate testing age
prelim.dataset["testing.age"] <- age(from = prelim.dataset$DOB, to = prelim.dataset$Biopsy.Date)


# write just in case

write.csv(prelim.dataset, file = "prelim.dataset.csv")

# 1/26/2017 - picking it back up

# cut down to FUN only
fun <- prelim.dataset[grep("FUN", prelim.dataset$subjectid), ]

# cut down missing gene expression data
fun.exp <- fun[which(!is.na(fun$mean.setas40)), ]

fun.exp.taste <- fun.exp[which(!is.na(fun.exp$MilliQWaterIntensity1)), ]

# see who is missing what; can only take complete data

# what is the final n??
fun.exp.taste$subjectid <- factor(fun.exp.taste$subjectid)
unique(fun.exp.taste$subjectid) # n = 26
sum(is.na(fun.exp.taste$rs713598)) # zero missing genotypes, that's good!

# i'm going to write out this cleaned file juuuuust in case!

write.csv(fun.exp.taste, file="prelim.clean.csv")




#######################
### you can start here....
# but need to re-do DOB, Biopsy Date, Biopsy Time, Date Taste....
# missing testing age:
## FUN-12 (missing DOB), FUN-17 (missing biopsy date), FUN-31 (missing DOB)
###########################

fun.exp.taste <- read.csv("prelim.clean.csv", header = T)

# lets check out taster gene expression correlated with taste test data
gene.exp.taste <- fun.exp.taste[grep("tas2r38", fun.exp.taste$target.name), ]

# calculate average PTC
gene.exp.taste$avg.ptc <- (gene.exp.taste$PTCIntensity4 + gene.exp.taste$PTCIntensity12)/2

################################################################
# CORRELATION TESTING - subsets FP only but not very useful ####
# cut to FP
# fp.taste <- gene.exp.taste[which(gene.exp.taste$sample.type == 'FP'), ]
# fp.taste <- fp.taste[which(fp.taste$target.name == 'tas2r38.taste'), ]
# 
# # ***note i am using *mean.setas40* just to get some code running. can also use median, etc etc.
# 
# # then, let's correlate the tongue vs. nose expression
# with(fp.taste, cor.test(mean.setas40, avg.ptc))
# 
# fp.taste$TAS2R38_Diplotype <- factor(fp.taste$TAS2R38_Diplotype)
# library(plyr)
# 
# # cut down to three columns?? there's gotta be a better way
# 
# # why do i need to convert biopsy time? something for ddply to work
# fp.taste$Biopsy.Time <- as.POSIXct(fp.taste$Biopsy.Time)
# 
# correlation.table <- ddply(fp.taste, .(TAS2R38_Diplotype, target.name), summarise,
#       Cor = cor(mean.setas40, avg.ptc))
# 
# # 2/6/17 - messing around trying to get additional stuff in ddply table - this works-ish!
# blahblah <- ddply(fp.taste, .(TAS2R38_Diplotype, target.name), summarise,
#                   n = length(unique(subjectid)),
#                   cor = cor(mean.setas40, avg.ptc),
#                   pvalue = cor.test(mean.setas40, avg.ptc)$estimate)
# 
# blah <- tidy(cor.test(fp.taste$mean.setas40, fp.taste$avg.ptc))
# cor(fp.taste$mean.setas40, fp.taste$avg.ptc)
####################################################################


# ###################
# # CORRELATIONS 2/8/2017 ####
# ###################

# why do i need to convert biopsy time? something for ddply to work
gene.exp.taste$Biopsy.Time <- as.POSIXct(gene.exp.taste$Biopsy.Time)

library(plyr)

# # original correlation for posterity
# correlation.taste <- ddply(gene.exp.taste, .(TAS2R38_Diplotype, target.name, sample.type), summarise,
#                            n = length(unique(subjectid)),
#                            cor = cor(mean.setas40, avg.ptc))
# 
# write.csv(correlation.taste, file="cor.table.exp.ptc.csv")

# if you use TRY it gets around the errors but fills resulting data.frame with crap
correlation.revamp <- ddply(gene.exp.taste, .(TAS2R38_Diplotype, target.name, sample.type), summarise,
                          n = length(unique(subjectid)),
                          cor = try(cor.test(mean.setas40, avg.ptc)$estimate),
                          pvalue = try(cor.test(mean.setas40, avg.ptc)$p.value))

# replace crap with NA - VICTORY IS MINE
correlation.revamp[correlation.revamp == "Error in cor.test.default(mean.setas40, avg.ptc) : \n  not enough finite observations\n"] <- NA

write.csv(correlation.blah, file="cor.table.exp.ptc.REVAMP.csv")

# just take a peek at genotype breakdown
table <- table(gene.exp.taste$TAS2R38_Diplotype)

# i should also figure out how to cast/melt the data so that i have
# subjectid        target       FP      Nasal
# FUN-1           taster        x         y
# FUN-2           taster        x         y

library(reshape2)

data.2 <- dcast(gene.exp.taste, subjectid + target.name + TAS2R38_Diplotype + avg.ptc ~ sample.type, value.var = "mean.setas40")

# cut it to FP and Biopsy
FP.biopsy <- data.2[which(!is.na(data.2$`Nasal Biopsy`)), ]

names(FP.biopsy)[names(FP.biopsy) == 'Nasal Biopsy'] <- 'nasal.biopsy'

# correlation between FP and nasal biopsies
correlation.biopsy <- ddply(FP.biopsy, .(TAS2R38_Diplotype, target.name), summarise,
                            n = length(unique(subjectid)),
                            cor = try(cor.test(FP, nasal.biopsy)$estimate),
                            pvalue = try(cor.test(FP, nasal.biopsy)$p.value))

correlation.biopsy[correlation.biopsy == "Error in cor.test.default(FP, nasal.biopsy) : \n  not enough finite observations\n"] <- NA

write.csv(correlation.biopsy, file="cor.table.FP.biopsy.REVAMP.csv")


# cut it to FP and Brush
FP.brush <- data.2[which(!is.na(data.2$`Nasal Brush`)), ]

names(FP.brush)[names(FP.brush) == 'Nasal Brush'] <- 'nasal.brush'

# correlation between FP and nasal brushes
correlation.brush <- ddply(FP.brush, .(TAS2R38_Diplotype, target.name), summarise,
                           n = length(unique(subjectid)),
                           cor = try(cor.test(FP, nasal.brush)$estimate),
                           pvalue = try(cor.test(FP, nasal.brush)$p.value))

correlation.brush[correlation.brush == "Error in cor.test.default(FP, nasal.brush) : \n  not enough finite observations\n"] <- NA

write.csv(correlation.brush, file="cor.table.FP.brush.REVAMP.csv")

# useless below
# # and can then correlate (perhaps facet it by TAS2R38?? 3 plots?)
# 
# library(ggplot2)
# 
# # scatter of mean.as40 values for different sample types.
# ## plot FP vs. Biopsy
# ## is this helpful? actually not, the jitter gives a false scatterplot!! need to rethink what needs to be graphed here
# ggplot(data.2, aes(x = 'FP', y = 'Nasal Biopsy')) +
#   geom_point() +
#   # geom_smooth(method=lm) + 
#   facet_wrap(~TAS2R38_Diplotype + target.name)

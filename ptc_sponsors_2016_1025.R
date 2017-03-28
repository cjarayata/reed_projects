# made queries in access; hopefully i got all the right data!!

setwd("S:/reed/Labshare/CJ Arayata/Twinsburg/Fall 2016 Sponsors Meeting")

library(pacman)
pacman::p_load(reshape2, tm, qdap, plyr, lsr, broom, ggplot2, gridExtra, psych)

# laborious master file creation: run as chunk  ----

# read in TWA - all data with PTC and genotypes
twa <- read.csv("twa.csv", header=T, na.string="")

# read in CNTL - cohen controls in 2015
control <- read.csv("control.csv", header=T, na.string="")

# start formatting; deleting columns; preparing for merge 

# to lower
colnames(twa) <- tolower(colnames(twa))
colnames(control) <- tolower(colnames(control))

# TWA - delete unnecessary variables
twa <- twa[, c(1:3, 10, 4:5, 7, 8, 9, 11, 14, 22, 19, 27, 28)]

# convert dates to proper format
twa$dob <- as.Date(twa$dob, format = "%m/%d/%Y")
twa$testdate <- as.Date(twa$testdate, format = "%m/%d/%Y")

# CONTROL

# copy over testing date so i can add testing day variable
control["testing.day"] <- control$testingdate

# use gsub to replace strings with other strings
control$testing.day <- gsub('8/8/2015', 'Saturday', control$testing.day)
control$testing.day <- gsub('8/9/2015', 'Sunday', control$testing.day)

# change class to date so i can calculate age
class(control$dob)
control$dob <- as.Date(control$dob, format = "%m/%d/%Y")

class(control$testingdate)
control$testingdate <- as.Date(control$testingdate, format = "%m/%d/%Y")


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
control["testing.age"] <- age(from = control$dob, to = control$testingdate)

# normalize the control taste data to be out of 7.6 instead of 12
control["ptc.norm"] <- round(control$ptcmeanintensity * (7.6/12), digits = 1)
control["water.norm"] <- round(control$h2omeanintensity * (7.6/12), digits = 1)

# reorder datasets
twa.clean <- twa[, c(1:9, 13, 11, 12, 14, 15)]
control.clean <- control[, c(1:3, 15, 4:6, 14, 7, 9, 16, 17, 12, 13)]

# check names
var.names <- cbind(as.data.frame(colnames(twa.clean)), as.data.frame(colnames(control.clean)))

# make the column names match each other
colnames(twa.clean) <- c("subjectid", "uniqueid", "dob", "testing.age", "sex", "race", "testing.date", "testing.day", "investigator", "sp.code", "ptc.intensity", "water.intensity", "snp", "genotype")
colnames(control.clean) <- colnames(twa.clean)

# now row bind both TWA and CNTL subjects together

ptc.full <- rbind(twa.clean, control.clean)

str(ptc.full)

# need to cast genotypes into separate variables, plus diplotype tomfoolery

# cast out
ptc.cast <- dcast(ptc.full, ... ~ snp, value.var = "genotype")

str(ptc.cast)

# dip and hap are two different variables; combine them and clean
ptc.cast$tas2r38 <- paste(ptc.cast$TAS2R38_Diplotype, ptc.cast$TAS2R38_Haplotype) # this paste introduces NA strings
ptc.cast$tas2r38 <- gsub("NA", "", ptc.cast$tas2r38) # remove NA strings
ptc.cast$tas2r38 <- stripWhitespace(ptc.cast$tas2r38) # take out whitespace just in case
ptc.cast$tas2r38 <- strip(ptc.cast$tas2r38, char.keep = "/") # takes out the () in P(A)V etc.
ptc.cast$tas2r38 <- toupper(ptc.cast$tas2r38) # changes pav to PAV
ptc.cast$tas2r38 <- replace(ptc.cast$tas2r38, ptc.cast$tas2r38=="", NA) # changes blanks to NA
ptc.cast$tas2r38 <- as.factor(ptc.cast$tas2r38) # factorize

# make sure the blanks are proper NAs
table(ptc.cast$tas2r38)

# checking out the file to be sure what's going on
str(ptc.cast)

# cut extra tas2r38 columns
ptc.cast<- ptc.cast[, c(1:15, 18)]

# turn genotype x's into NAs - dumb code - should have done this pre-casting

ptc.cast$rs10246939 <- gsub("X", NA, ptc.cast$rs10246939)
ptc.cast$rs10246939 <- as.factor(ptc.cast$rs10246939)
levels(ptc.cast$rs10246939)

ptc.cast$rs1726866 <- gsub("X", NA, ptc.cast$rs1726866)
ptc.cast$rs1726866 <- as.factor(ptc.cast$rs1726866)
levels(ptc.cast$rs1726866)

ptc.cast$rs713598 <- gsub("X", NA, ptc.cast$rs713598)
ptc.cast$rs713598 <- as.factor(ptc.cast$rs713598)
levels(ptc.cast$rs713598)

str(ptc.cast)

# start to change race to factor
class(ptc.cast$race)

ptc.cast$race <- replace(ptc.cast$race, ptc.cast$race==1, 'European American')
ptc.cast$race <- replace(ptc.cast$race, ptc.cast$race==4, 'African American')

# check investigator names - kind of a mess - need to fix
table(ptc.cast$investigator)
ptc.cast$investigator <- tolower(ptc.cast$investigator)
ptc.cast$investigator <- replace(ptc.cast$investigator, ptc.cast$investigator=='bard', 'brad')
ptc.cast$investigator <- replace(ptc.cast$investigator, ptc.cast$investigator=='danielle', 'dani')
ptc.cast$investigator <- replace(ptc.cast$investigator, ptc.cast$investigator=='laura + michelle', 'laura')
ptc.cast$investigator <- replace(ptc.cast$investigator, ptc.cast$investigator=='mathew', 'matthew')
ptc.cast$investigator <- droplevels(as.factor(ptc.cast$investigator))
ptc.cast$investigator <- factor(ptc.cast$investigator)

table(ptc.cast$investigator)

# this is convoluted but eventually anonymizes tester names into letters
levels(ptc.cast$investigator)
ptc.cast$investigator.anon <- ptc.cast$investigator
levels(ptc.cast$investigator.anon) <- c(LETTERS[1:nlevels(ptc.cast$investigator)])

# checking counts.... OK
table(ptc.cast$investigator)
table(ptc.cast$investigator.anon)

# sex clean
table(ptc.cast$sex)
levels(ptc.cast$sex)

# OK well i guess finally ready for first dataset!
str(ptc.cast)
visit.prep <- ptc.cast

# FIRST VISIT prep: run as chunk - file is "first.bw.dip" ----
head(visit.prep)
str(visit.prep)

# sort by date so that i can then take first instance of unique - MUST be done first!
visit.prep <- visit.prep[with(visit.prep, order(testing.date)), ]

# now get first instance of unique
first <- visit.prep[!duplicated(visit.prep$uniqueid), ]

# cut race to only white and black - then factor
first.bw <- first[first$race %in% c('European American', 'African American'), ]
class(first.bw$race)
levels(first.bw$race)
first.bw$race <- factor(first.bw$race)

# SP code prep
table(first.bw$sp.code)

# make 0 and 4 NA, then when plotting, drop it
first.bw$sp.code <- replace(first.bw$sp.code, first.bw$sp.code==0, NA)
first.bw$sp.code <- replace(first.bw$sp.code, first.bw$sp.code==4, NA)
first.bw$sp.code <- factor(first.bw$sp.code)
levels(first.bw$sp.code)
class(first.bw$sp.code)

# cut to three main diplotypes - for whatever reason, %in% works but [which] does not!!
first.bw.dip <- first.bw[first.bw$tas2r38 %in% c('PAV/PAV', 'AVI/AVI', 'AVI/PAV'), ]

# checking out dip counts just to be sure
table(first.bw$tas2r38) # 'full' dataset with all diplotype groups
table(first.bw.dip$tas2r38)

# drop unused levels (as.factor doesn't work but 'factor' does?)
first.bw.dip$tas2r38 <- factor(first.bw.dip$tas2r38)

# check out n's
xtabs(~ tas2r38 + race, data=first.bw.dip)

# median split on AGE - then using "young" vs. "old" factor
median(na.omit(first.bw.dip$testing.age)) # median is 33
count(first.bw.dip$testing.age < 33) # 868 / 869
count(first.bw.dip$testing.age >= 33) # 869 / 868; OK
first.bw.dip$age.grp <- first.bw.dip$testing.age # copy

first.bw.dip$age.grp[first.bw.dip$age.grp < 33] <- 'young'
first.bw.dip$age.grp[first.bw.dip$age.grp >= 33 & first.bw.dip$age.grp != 'young'] <- 'old'

first.bw.dip$age.grp <- factor(first.bw.dip$age.grp)

# make parallel dataset for LAST VISIT ----

# count visits
visits <- as.data.frame(table(visit.prep$uniqueid))
colnames(visits) <- c("uniqueid", "num.visits")

# tie number of visits onto full dataset as "visits" variable

visits.2 <- merge(visit.prep, visits)
visits.2$experienced <- visits.2$num.visits

# sort by date so that i can then take LAST (most recent) instance of unique - MUST be done first!
visits.2 <- visits.2[with(visits.2, rev(order(testing.date))), ]

# now get LAST instance of unique
last <- visits.2[!duplicated(visits.2$uniqueid), ]

# recode number of visits
last$experienced[last$experienced == '1'] <- 'naive'
last$experienced[last$experienced == '2' | last$experienced == '3'] <- 'experienced'
last$experienced[last$experienced != 'naive' & last$experienced != 'experienced'] <- 'very experienced'
last$experienced <- factor(last$experienced)

# subset on race
last.bw <- last[last$race %in% c('European American', 'African American'), ]
class(last.bw$race)
last.bw$race <- factor(last.bw$race)
levels(last.bw$race)

# SP code prep
table(last.bw$sp.code)

# make 0 and 4 NA, then when plotting, drop it
last.bw$sp.code <- replace(last.bw$sp.code, last.bw$sp.code==0, NA)
last.bw$sp.code <- replace(last.bw$sp.code, last.bw$sp.code==4, NA)
last.bw$sp.code <- factor(last.bw$sp.code)
levels(last.bw$sp.code)
class(last.bw$sp.code)

# cut to three main diplotypes - for whatever reason, %in% works but [which] does not!!
last.bw.dip <- last.bw[last.bw$tas2r38 %in% c('PAV/PAV', 'AVI/AVI', 'AVI/PAV'), ]

# checking out dip counts just to be sure
table(last.bw$tas2r38) # 'full' dataset with all diplotype groups
table(last.bw.dip$tas2r38)

# drop unused levels (as.factor doesn't work but 'factor' does?)
last.bw.dip$tas2r38 <- factor(last.bw.dip$tas2r38)

# need to order experience levels
last.bw.dip$experienced <- ordered(last.bw.dip$experienced, levels = c("naive", "experienced", "very experienced"))

# check out n's
xtabs(~ tas2r38 + race, data=last.bw.dip)

# FOR DANI 9-22-16: look within groups age, sex, race, experience ----

# first, EVERYONE anova that shows the relation between genotype - phenotype = 'base level' eta squared
{
# anova TAS2R38 (3 main) - 69% variance
aov <- aov(ptc.intensity ~ tas2r38, data = first.bw.dip)
aov <- tidy(etaSquared(aov, type=1, anova=TRUE))

# simple plot - could be prettier
ggplot(data=first.bw.dip, aes(x = tas2r38, y = ptc.intensity)) +
  geom_boxplot()
}
# then, AGE GROUP (median split at 33) - not significant
{
aov.age <- aov(ptc.intensity ~ age.grp, data = first.bw.dip)
aov.age <- tidy(etaSquared(aov.age, type=1, anova=T))

# order age groups
first.bw.dip$age.grp <- ordered(first.bw.dip$age.grp, levels = c("young", "old"))

# plot RAW AGE - continuous variable
ggplot(first.bw.dip, aes(x = testing.age, y = ptc.intensity)) +
  geom_point() +
  facet_wrap(~tas2r38)

# plot AGE GROUP (dichotomized) variable
ggplot(first.bw.dip[!is.na(first.bw.dip$age.grp), ], aes(x = age.grp, y = ptc.intensity)) +
  geom_boxplot() +
  facet_wrap(~tas2r38)

xtabs(~ tas2r38 + age.grp, data=first.bw.dip)
xtabs(~ tas2r38 + race, data=last.bw.dip)

## difference in anova when subsetting old/young
aov.age.o <- aov(ptc.intensity ~ tas2r38, data = first.bw.dip[first.bw.dip$age.grp=='old', ])
aov.age.o <- tidy(etaSquared(aov.age.o, type=1, anova=T))
# old - 67%

aov.age.y <- aov(ptc.intensity ~ tas2r38, data = first.bw.dip[first.bw.dip$age.grp=='young', ])
aov.age.y <- tidy(etaSquared(aov.age.y, type=1, anova=T))
# young - 72%
}
# then, SEX - significant but small eta.sq
{
aov.sex <- aov(ptc.intensity ~ sex, data = first.bw.dip)
aov.sex <- tidy(etaSquared(aov.sex, type=1, anova=T))

ggplot(first.bw.dip[!is.na(first.bw.dip$sex), ], aes(x = sex, y = ptc.intensity)) +
  geom_boxplot() +
  facet_wrap(~tas2r38)

## difference in anova when subsetting male/female
aov.sex.m <- aov(ptc.intensity ~ tas2r38, data = first.bw.dip[first.bw.dip$sex=='M', ])
aov.sex.m <- tidy(etaSquared(aov.sex.m, type=1, anova=T))
# men - 66%

aov.sex.f <- aov(ptc.intensity ~ tas2r38, data = first.bw.dip[first.bw.dip$sex=='F', ])
aov.sex.f <- tidy(etaSquared(aov.sex.f, type=1, anova=T))
# women - 70% - slightly better than men?
}
# then, RACE - significant but small eta.sq
{
aov.race <- aov(ptc.intensity ~ race, data = first.bw.dip)
aov.race <- tidy(etaSquared(aov.race, type=1, anova=T))

ggplot(first.bw.dip[!is.na(first.bw.dip$race), ], aes(x = race, y = ptc.intensity)) +
  geom_boxplot() +
  facet_wrap(~tas2r38)

## difference in anova when subsetting white/black
aov.race.w <- aov(ptc.intensity ~ tas2r38, data = first.bw.dip[first.bw.dip$race=='european american', ])
aov.race.w <- tidy(etaSquared(aov.race.w, type=1, anova=T))
# white - 71%

# plot white
w <- ggplot(first.bw.dip[first.bw.dip$race=='european american', ], aes(x = tas2r38, y = ptc.intensity)) +
  geom_boxplot() +
  ggtitle("TAS2R38 -> PTC Intensity for Euro Americans - 71%")

aov.race.b <- aov(ptc.intensity ~ tas2r38, data = first.bw.dip[first.bw.dip$race=='african american', ])
aov.race.b <- tidy(etaSquared(aov.race.b, type=1, anova=T))
# black - 27% - staggering difference!

# plot black
b <- ggplot(first.bw.dip[first.bw.dip$race=='african american', ], aes(x = tas2r38, y = ptc.intensity)) +
  geom_boxplot() +
  ggtitle("TAS2R38 -> PTC Intensity for African Am. - 26%")

# package 'gridExtra' can arrange multiple objects onto one plot. kind of like faceting but can be much more complex
grid.arrange(w, b, ncol=2)
}
# then, TESTERS - not significant
{
# make a data frame that is the number of tests each person has administered
n.tested <- as.data.frame(table(first.bw.dip$investigator.anon))

# merge with full dataset; then subset
testmerge <- merge(first.bw.dip, n.tested, by.x = 'investigator.anon', by.y = 'Var1')
testers <- testmerge[testmerge$Freq > 50, ]

# then, global test for effect of 'trained' testers (those >50)
aov.testers <- aov(ptc.intensity ~ investigator.anon, data = testers)
aov.testers <- tidy(etaSquared(aov.testers, type=1, anova=T))

# i guess i didn't do this before, but let's factorize (and drop unused levels)
levels(testers$investigator.anon)
testers$investigator.anon <- factor(testers$investigator.anon)
levels(testers$investigator.anon)

# plot with filtered testers
ggplot(data=testers, aes(x = investigator.anon, y = ptc.intensity)) +
  geom_boxplot() +
  facet_wrap(~ tas2r38)



# use 'psych' package to get a nice summary table, note the 'mat=TRUE' call to output to matrix/data frame
tester.summary <- describeBy(testers$ptc.intensity, group = testers$investigator.anon, mat=T)

## difference in anova when testing best/worst tester? (sort tester.summary by SKEW, a hunch)

# lowest skew - tester V
aov.tester.v <- aov(ptc.intensity ~ tas2r38, data = testers[testers$investigator.anon=='V', ])
aov.tester.v <- tidy(etaSquared(aov.tester.v, type=1, anova=T))
# tester V - 70% variance explained

# highest skew - tester W
aov.tester.w <- aov(ptc.intensity ~ tas2r38, data = testers[testers$investigator.anon=='W', ])
aov.tester.w <- tidy(etaSquared(aov.tester.w, type=1, anova=T))
# tester W - 19% variance explained?? wowzers that's bad

# dani, you will like this. just because i am curious let's get a tester loop going!!

#build table first - aov with partial eta, put in tester info, and rearrange table
aov.tester.loop <- aov(ptc.intensity ~ tas2r38, data=testers[testers$investigator.anon=='A', ])
aov.tester.loop <- tidy(etaSquared(aov.tester.loop, type=1, anova=TRUE))
aov.tester.loop["tester"] <- "A"
aov.tester.loop <- aov.tester.loop[, c(9, 1, 4:8, 2, 3)]

#loop through all testers (A will be in there twice, sorry)
for (i in levels(testers$investigator.anon)){
  temp <- aov(ptc.intensity ~ tas2r38, data=testers[testers$investigator.anon==i, ])
  temp <- tidy(etaSquared(temp, type=1, anova=TRUE))
  temp["tester"] <- i
  temp <- temp[, c(9, 1, 4:8, 2, 3)]
  aov.tester.loop <- rbind(aov.tester.loop, temp)
}

# now just by viewing/sorting 'aov.tester.loop', we can check out the % variance explained for each tester!
}
# then, EXPERIENCE - on the cusp
{
## 1 - naive, 2 or 3 tests - experienced, 4+ tests are very experienced (probably could have subsetted better based on 'years' participated but oh well)

# global anova
aov.exp <- aov(ptc.intensity ~ experienced, data = last.bw.dip)
aov.exp <- tidy(etaSquared(aov.exp, type = 1, anova = T))

# plot it
ggplot(data=last.bw.dip, aes(x = tas2r38, y = ptc.intensity)) +
  geom_boxplot() +
  facet_wrap(~ experienced)

ggplot(data=last.bw.dip, aes(x = experienced, y = ptc.intensity)) +
  geom_boxplot() +
  facet_wrap(~ tas2r38)

# build table w. naive subjects - aov with partial eta, put in experience info, and rearrange table
aov.exp.loop <- aov(ptc.intensity ~ tas2r38, data = last.bw.dip[last.bw.dip$experienced=='naive', ])
aov.exp.loop <- tidy(etaSquared(aov.exp.loop, type=1, anova=TRUE))
aov.exp.loop["experience"] <- "naive"
aov.exp.loop <- aov.exp.loop[, c(9, 1, 4:8, 2, 3)]

#loop through 3 experience levels (naive will be in there twice, sorry)
for (i in levels(last.bw.dip$experienced)){
  temp <- aov(ptc.intensity ~ tas2r38, data=last.bw.dip[last.bw.dip$experienced==i, ])
  temp <- tidy(etaSquared(temp, type=1, anova=TRUE))
  temp["experience"] <- i
  temp <- temp[, c(9, 1, 4:8, 2, 3)]
  aov.exp.loop <- rbind(aov.exp.loop, temp)
}
}
# then, TEST DAY - not significant
{
table(first.bw.dip$testing.day)

aov.day <- aov(ptc.intensity ~ testing.day, data = first.bw.dip)
aov.day <- tidy(etaSquared(aov.day, type=1, anova=TRUE))

ggplot(data=first.bw.dip, aes(x = testing.day, y = ptc.intensity)) +
  geom_boxplot() +
  facet_wrap(~ tas2r38)

## subset saturday vs. sunday testing
# need to do

}
# then, SP CODE - not significant (small numbers for "3")
{
aov.sp <- aov(ptc.intensity ~ sp.code, data = first.bw.dip)
aov.sp <- tidy(etaSquared(aov.sp, type=1, anova=TRUE))

# need to drop NA!
ggplot(data=first.bw.dip[!is.na(first.bw.dip$tas2r38) & !is.na(first.bw.dip$sp.code), ], aes(x = sp.code, y = ptc.intensity)) +
  geom_boxplot() +
  facet_wrap(~ tas2r38)
}

# # (old) main analyses: TAS2R38, age, sex, race ----
# to uncomment a whole section, highlight and CONTROL + SHIFT + C (thanks yusuke)
# str(first.bw.dip)
# 
# # anova TAS2R38 (3 main) - 69% variance
# aov <- aov(ptc.intensity ~ tas2r38, data = first.bw.dip)
# aov <- tidy(etaSquared(aov, type=1, anova=TRUE))
# 
# # equivalent calls but just checking output to be sure
# lm <- lm(ptc.intensity ~ tas2r38, data=first.bw.dip)
# summary(lm)
# anova(lm)
# 
# glm <- glm(ptc.intensity ~ tas2r38, data=first.bw.dip)
# summary(glm)
# anova(glm)
# 
# # plot - could be prettier
# ggplot(data=first.bw.dip, aes(x = tas2r38, y = ptc.intensity)) +
#   geom_boxplot()
# 
# # anova TAS2R38 (type I) plus AGE (continuous)
# aov.age <- aov(ptc.intensity ~ (tas2r38 * testing.age), data = first.bw.dip)
# aov.age <- tidy(etaSquared(aov.age, type=1, anova=TRUE))

# editing the facet labels with a helper
haplotypes <- list(
  'AVI/AVI'="AVI",
  'AVI/PAV'="AVI/PAV",
  'PAV/PAV'="PAV"
)
hap_label <- function(variable,value){
  return(haplotypes[value])
}

# plot - again, could be prettier
ggplot(data=first.bw.dip, aes(x=testing.age, y=ptc.intensity)) +
  geom_jitter() +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 7.6),
                     limits = c(0, 7.6)) +
  labs(x = "Testing Age",
       y = "PTC Intensity\n") +
  geom_smooth(method = 'lm', span = 0.5) +
  facet_wrap(~ tas2r38, labeller = hap_label) +
  theme(text = element_text(size = 26))
  
# 
# # anova plus RACE
# aov.agerace <- aov(ptc.intensity ~ (tas2r38 * testing.age * race), data = first.bw.dip)
# aov.agerace <- tidy(etaSquared(aov.agerace, type=1, anova=TRUE))
# 
# # simple plot of RACE
# ggplot(data=first.bw.dip, aes(x = tas2r38, y = ptc.intensity)) +
#   geom_boxplot() +
#   facet_wrap(~ race)
# 
# # plot AGE, RACE, and TAS2R38
# ggplot(data=first.bw.dip, aes(x=testing.age, y=ptc.intensity)) +
#   geom_jitter() +
#   geom_smooth(span = 0.5) +
#   geom_hline(yintercept = 7.6) + # line drawn at 7.6 to show max, but not sure if that's useful
#   facet_wrap(~ race + tas2r38)

# 10-18-16: dotplot things for sponsors slides----

# all data graph 
all <- ggplot(first.bw.dip, aes(x = 1, y = ptc.intensity)) +
  geom_dotplot(binaxis = "y",
               binwidth = 0.1,
               stackdir="centerwhole",
               dotsize = 0.3) +
  stat_summary(fun.y = "median",
               fun.ymin = "median",
               fun.ymax = "median",
               geom = "crossbar",
               width = 0.4,
               color = "red") +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 7.6),
                     limits = c(0, 7.6)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 26)) +
  labs(x = "\nAll Participants",
       y = "PTC Intensity\n")

all

all.2 <- ggplot(first.bw.dip, aes(x = 1, y = ptc.intensity)) +
  geom_dotplot(binaxis = "y",
               binwidth = 0.1,
               stackdir="centerwhole",
               dotsize = 0.3) +
  stat_summary(fun.y = "median",
               fun.ymin = "median",
               fun.ymax = "median",
               geom = "crossbar",
               width = 0.4,
               color = "red") +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 7.6),
                     limits = c(0, 7.6)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 26)) +
  labs(x = "\nAll Participants",
       y = "PTC Intensity Ratings, visual analog scale\n")

all.2

# count n's just in case
length(unique(first.bw.dip$uniqueid)) # n = 1755

# sex graph
{
sex <- ggplot(data = subset(first.bw.dip, !is.na(sex)),
              aes(x = sex, y = ptc.intensity, color = tas2r38)) +
  geom_dotplot(binaxis = "y",
               binwidth = 0.1,
               binpositions="all",
               stackdir="centerwhole",
               dotsize = 0.5) +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 7.6),
                     limits = c(0, 7.6)) +
  stat_summary(fun.y = "median",
               fun.ymin = "median",
               fun.ymax = "median",
               geom = "crossbar",
               width = 0.4,
               color = "red") +
               #position = position_dodge(width = 1)) +
  theme(legend.position = "none",
        text = element_text(size = 26)) +
  labs(x = "Sex",
       y = "PTC Intensity\n") +
  facet_wrap(~ tas2r38, labeller=hap_label)

sex
}

# i cant figure out how to get the median crossbar length to be representative of the actual count of the median value
count(subset(first.bw.dip, sex=='M')$ptc.intensity) == 4.2))
length(which(subset(first.bw.dip, sex=='M')$ptc.intensity == median(subset(first.bw.dip, sex=='M')$ptc.intensity, na.rm=T)))


# race
{
race <- ggplot(data = subset(first.bw.dip, !is.na(race)),
              aes(x = race, y = ptc.intensity, color = tas2r38)) +
  geom_dotplot(binaxis = "y",
               binwidth = 0.1,
               binpositions="all",
               stackdir="centerwhole",
               dotsize = 0.5) +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 7.6),
                     limits = c(0, 7.6)) +
  stat_summary(fun.y = "median",
               fun.ymin = "median",
               fun.ymax = "median",
               geom = "crossbar",
               width = 0.4,
               color = "red") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12), 
        text = element_text(size = 26)) +
  labs(x = "\nRace",
       y = "PTC Intensity\n") +
  facet_wrap(~ tas2r38, labeller = hap_label)

race
}

genotype <- ggplot(data = subset(first.bw.dip, !is.na(tas2r38)),
                   aes(x = tas2r38, y = ptc.intensity, color = tas2r38)) +
  geom_dotplot(binaxis = "y",
               binwidth = 0.1,
               binpositions="all",
               stackdir="centerwhole",
               dotsize = 0.5) +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 7.6),
                     limits = c(0, 7.6)) +
  stat_summary(fun.y = "median",
               fun.ymin = "median",
               fun.ymax = "median",
               geom = "crossbar",
               width = 0.5,
               color = "red") +
  theme(legend.position = "none",
        text = element_text(size = 26)) +
  labs(x = "TAS2R38 Genotype",
       y = "PTC Intensity\n")

genotype

# tester dotplot - looks bad
tester <- ggplot(data = subset(testers, !is.na(race)),
               aes(x = investigator.anon, y = ptc.intensity, color = investigator.anon)) +
  geom_dotplot(binaxis = "y",
               binwidth = 0.1,
               binpositions="all",
               stackdir="centerwhole",
               dotsize = 0.8) +
  coord_cartesian(ylim=c(0, 7.6)) +
  stat_summary(fun.y = "median",
               fun.ymin = "median",
               fun.ymax = "median",
               geom = "crossbar",
               width = 0.4,
               color = "red") +
  theme(legend.position = "none") +
  labs(x = "Tester",
       y = "PTC Intensity") +
  facet_wrap(~ tas2r38)

tester

# another try - still bad - very sparse
test.2 <- ggplot(data = testers[testers$investigator.anon %in% c('J', 'W'), ],
                 aes(x = investigator.anon, y = ptc.intensity, color = investigator.anon)) +
  geom_dotplot(binaxis = "y",
               binwidth = 0.1,
               binpositions="all",
               stackdir="centerwhole",
               dotsize = 0.3) +
  coord_cartesian(ylim=c(0, 7.6)) +
  stat_summary(fun.y = "median",
               fun.ymin = "median",
               fun.ymax = "median",
               geom = "crossbar",
               width = 0.4,
               color = "red") +
  theme(legend.position = "none") +
  labs(x = "Tester",
       y = "PTC Intensity") +
  facet_wrap(~ tas2r38)

# here's the boxplot for tester J (best tester) vs. tester W (worst tester)
ggplot(data=testers[testers$investigator.anon %in% c('J', 'W'), ],
       aes(x = investigator.anon, y = ptc.intensity)) +
  geom_boxplot() +
  labs(x = "Tester",
       y = "PTC Intensity\n") +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 7.6),
                     limits = c(0, 7.6)) +
  theme(text = element_text(size = 26)) +
  facet_wrap(~ tas2r38, labeller = hap_label)

# arrange in a layout such that 'all data' takes up 1/3 of page, and sex takes up 2/3s
layout <- grid.arrange(all, sex, ncol = 2,
                       layout_matrix = rbind(c(1, 2, 2)))

layout <- grid.arrange(all, sex, ncol = 2,
                       widths = c(1, 2))

# arrange in a layout such that 'all data' takes up 1/3 of page, and race takes up 2/3s
layout2 <- grid.arrange(all, race, ncol = 2,
                       layout_matrix = rbind(c(1, 2, 2)))

# now all and genotype
layout3 <- grid.arrange(all, genotype, ncol = 2,
                        layout_matrix = rbind(c(1, 2, 2)))

# this lines up x-axes but don't know how to control panel split
grid::grid.draw(gridExtra:::cbind_gtable(ggplotGrob(all),ggplotGrob(sex)))



# 10-25-16: comparing twinsburg subjects and cohen controls
## slight misnomer - control data collected during twinsburg 2015 (but using clinic method)

# create a 'study' variable so i can compare both methods
first.bw.dip$study <- first.bw.dip$subjectid # copy over subject id
first.bw.dip$study <- gsub('\\d', '', first.bw.dip$study) # strip id numbers
first.bw.dip$study <- strip(first.bw.dip$study) # get rid of '-' dash
first.bw.dip$study <- gsub('cntl', 2, first.bw.dip$study) # rename for plotting labels
first.bw.dip$study <- gsub('twa', 1, first.bw.dip$study)

# factorize
class(first.bw.dip$study)
first.bw.dip$study <- factor(first.bw.dip$study)
levels(first.bw.dip$study)

table(first.bw.dip$study)

# study ANOVA - significant but small eta
  aov.study <- aov(ptc.intensity ~ study, data = first.bw.dip)
  aov.study <- tidy(etaSquared(aov.study, type=1, anova=T))

# boxplot
  ggplot(first.bw.dip[!is.na(first.bw.dip$study), ], aes(x = study, y = ptc.intensity)) +
    geom_boxplot() +
    facet_wrap(~tas2r38)
  
  ## difference in anova when subsetting
  aov.study.twa <- aov(ptc.intensity ~ tas2r38, data = first.bw.dip[first.bw.dip$study=='Twinsburg', ])
  aov.study.twa <- tidy(etaSquared(aov.study.twa, type=1, anova=T))
  # variance is 70%
  
  aov.study.cntl <- aov(ptc.intensity ~ tas2r38, data = first.bw.dip[first.bw.dip$study=='Clinic', ])
  aov.study.cntl <- tidy(etaSquared(aov.study.cntl, type=1, anova=T))
  # variance is 66%
  
# dotplot like all of the rest
  study.graph <- ggplot(data = first.bw.dip,
                   aes(x = study, y = ptc.intensity, color = tas2r38)) +
    geom_dotplot(binaxis = "y",
                 binwidth = 0.1,
                 binpositions="all",
                 stackdir="centerwhole",
                 dotsize = 0.5) +
    scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 7.6),
                       limits = c(0, 7.6)) +
    stat_summary(fun.y = "median",
                 fun.ymin = "median",
                 fun.ymax = "median",
                 geom = "crossbar",
                 width = 0.4,
                 color = "red") +
    theme(legend.position = "none",
          text = element_text(size = 26)) +
    labs(x = "Testing Method",
         y = "PTC Intensity\n") +
    facet_wrap(~ tas2r38)
  
  study.graph
  
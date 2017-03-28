# twinsburg milk stuff for sponsors - 10-26-16
# edit 11-16-2016 for twin 1 vs. twin 2

setwd("//monell/storage/reed/Labshare/CJ Arayata/Twinsburg/2016 Milk Fall Sponsors - Extracurricular Time of Day")

# load packages ----
library(pacman)
pacman::p_load(lubridate, ggplot2, psych, gridExtra, plyr, dplyr, broom)

# read in raw - from my own access query
twins <- read.csv("twin_milk_10-26-16.csv", header=T, na.string="")

# newer version that has twin ID - for splitting zygosity
twins <- read.csv("twin_milk_11-16-16.csv", header=T, na.string = "")

# formatting colnames
colnames(twins) <- tolower(colnames(twins))
names(twins)[16] <- "total.correct"
names(twins)[17] <- "correct.fat.first"
names(twins)[18] <- "correct.fat.second"

# real quick for cailu
ggplot(twins, aes(x = howoftendrinkmilk, y = total.correct)) +
  geom_bar(stat="identity")

# mess with dates
twins$testdate <- as_date(mdy(twins$testdate))

# checking for weird misspellings or NA's, etc.... everything looks good
table(twins$tester)
table(twins$sex)
table(twins$race)

# cut down to first visit - first order by test date, then by subject id ----
twin.prep <- twins[with(twins, order(testdate, subjectid)), ]

# take first instance of unique ID
first <- twin.prep[!duplicated(twin.prep$uniqueid), ]

# check length
length(unique(first$uniqueid)) # n = 463

# define a summary table function - this is modified from a ddply call, since i have no grouping variables to subset by
cj.summary <- function(subjectid, age, sex, race){
  c(n.test.total = length(unique(subjectid)), # total tested
    age.mean = mean(age, na.rm=T), # mean age
    age.sd = sd(age, na.rm=T), # age sd
    n.white = as.numeric(table(race)[names(table(race)) == 1]), # num white participants
    perc.white = as.numeric(table(race)[names(table(race)) == 1]) / length(unique(subjectid)), # percent white
    n.female = as.numeric(table(sex)[names(table(sex)) == 'F']), # number female
    perc.female = as.numeric(table(sex)[names(table(sex)) == 'F']) / length(unique(subjectid)) # percent female
  )
}

# run summary function with required subjectid, age, sex, race arguments, and transpose result
table <- t(cj.summary(
  first$subjectid,
  first$testing.age,
  first$sex,
  first$race))

# not necessary but would irk me otherwise
row.names(table) <- "2016 unique summary"

# use psych package to get another table - semi useful but calculates for ALL variables
blah <- describe(first)


# test-retest - means need to make a dataset of sat-sunday repeats
# KEEP duplicates of unique id - both SAT and SUN - n = 254 data points
duplicates <- twin.prep[twin.prep$uniqueid %in% twin.prep$uniqueid[duplicated(twin.prep$uniqueid)], ]

length(unique(duplicates$uniqueid)) # n = 127 repeat customers

# 463 one-visits + 127 repeats = 590 total subjects = checks out

# correlations between sat and sunday testing

# create subsets
sat <- duplicates[which(duplicates$testing.day == 'Saturday'), ]
sun <- duplicates[which(duplicates$testing.day == 'Sunday'), ]

# sort both by unique id - crucial so correlation is paired properly!!
sat <- sat[order(sat$uniqueid), ]
sun <- sun[order(sun$uniqueid), ]

# plot histogram (unique) and scatter (reliability) ----
# histogram
ggplot(data = first, aes(total.correct)) +
  geom_histogram(binwidth=0.5, fill = "black") +
  labs(x = "Total Correct",
       y = "Frequency") +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  geom_vline(aes(xintercept = mean(total.correct, na.rm=T)),
             color = "red", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = 5),
             color = "red", linetype = "solid", size = 1) +
  stat_function(fun = function(x, mean, sd, n){
    n * dnorm(x = x, mean = mean, sd = sd)
  }, 
  args = with(first, c(mean = mean(total.correct), sd = sd(total.correct), n
                       = length(total.correct)))) +
  theme_bw() +
  theme(text = element_text(size = 26))
  

# one-sample t-test against 5 (better than chance?)
t <- t.test(first$total.correct, mu = 5)

t
# t(462) = 7.9135, p = 1.869e-14
# mean of x = 5.812095


# run correlation ----
cor.test(sat$total.correct, sun$total.correct, method = 'pearson')
# Pearson's product-moment correlation
# 
# data:  sat$total.correct and sun$total.correct
# t = 7.4077, df = 125, p-value = 1.676e-11
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.4183747 0.6627715
# sample estimates:
#       cor 
# 0.5523316 

# jittered scatterplot
ggplot(data = data.frame(x = sat$total.correct, y = sun$total.correct), aes(x, y)) +
  geom_jitter() +
  geom_smooth(method = lm, se = T) +
  labs(x = "Saturday Correct",
       y = "Sunday Correct") +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10),
                     limits = c(0, 10)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10),
                     limits = c(0, 10)) +
  theme(text = element_text(size = 26))

# boxplot
ggplot(data = duplicates, aes(x = testing.day, y = total.correct)) +
  geom_boxplot()



# working with unique dataset, do SEX, RACE, AGE effects ----

# SEX
sex.plot <- ggplot(first, aes(x = sex, y = total.correct, fill = sex)) + 
  geom_violin(aes(fill=sex)) +
  geom_dotplot(binaxis='y',
               binwidth = 0.05,
               stackdir='center',
               dotsize = 1.2) + 
  scale_y_continuous(breaks=1:10) + 
  labs(x = "Sex",
       y = "Total Correct") +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 26))
  
sex.plot

t.test(first$total.correct ~ first$sex)
# Welch Two Sample t-test
# 
# data:  first$total.correct by first$sex
# t = -0.36408, df = 270.45, p-value = 0.7161
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.4984927  0.3428970
# sample estimates:
#   mean in group F mean in group M 
# 5.790419        5.868217 

sex.table <- ddply(first, ~ sex, summarise,
                    n.test.total = length(unique(subjectid)), # number tested
                    total.mean = mean(total.correct, na.rm=T), # mean age
                    total.sd = sd(total.correct, na.rm=T)) # sd age


# RACE

# recode so that 1 = euro, 4 = african, anything else = other
first$race <- replace(first$race, first$race==1, 'European American')
first$race <- replace(first$race, first$race==4, 'African American')
first$race <- replace(first$race, first$race!= 'European American' & first$race!= 'African American', 'Other')
table(first$race)
class(first$race)


# make race a factor
first$race <- factor(first$race)
class(first$race)

# ordering is for plotting purposes i swear!!
first$race <- ordered(first$race, levels = c("European American", "African American", "Other"))


race.plot <- ggplot(data = subset(first, !is.na(race)),
                    aes(x = race, y = total.correct, fill = race)) + 
  geom_violin(aes(fill=race)) +
  geom_dotplot(binaxis='y',
               binwidth = 0.05,
               stackdir='centerwhole',
               dotsize = 1.2) + 
  scale_y_continuous(breaks=1:10) + 
  labs(x = "Race",
       y = "Total Correct") +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 24))

race.plot

# subset for a race t-test
first.bw <- first[first$race %in% c('European American', 'African American'), ]
t.test(first.bw$total.correct ~ first.bw$race)

# Welch Two Sample t-test
# 
# data:  first.bw$total.correct by first.bw$race
# t = -0.038858, df = 69.011, p-value = 0.9691
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.6035226  0.5804605
# sample estimates:
#   mean in group European American  mean in group African American 
# 5.831606                        5.843137 

# wait, i could just run an anova on white/black/other
race.anova <- aov(total.correct ~ race, data = first)
summary(race.anova)
#               Df Sum Sq Mean Sq F value Pr(>F)
# race          2    4.3   2.134   0.439  0.645
# Residuals   457 2222.3   4.863               
# 3 observations deleted due to missingness


race.table <- ddply(first, ~ race, summarise,
                 n.test.total = length(unique(subjectid)), # number tested
                 total.mean = mean(total.correct, na.rm=T), # mean age
                 total.sd = sd(total.correct, na.rm=T)) # sd age

# SP CODE

# prep
table(first$sp.code)
first$sp.code <- replace(first$sp.code, first$sp.code==1, 'Excellent')
first$sp.code <- replace(first$sp.code, first$sp.code==2, 'Good')
first$sp.code <- replace(first$sp.code, first$sp.code==3, 'Poor')
table(first$sp.code)

first$sp.code <- factor(first$sp.code)

# plot
sp.plot <- ggplot(data = subset(first, !is.na(sp.code)),
                    aes(x = sp.code, y = total.correct, fill = sp.code)) + 
  geom_violin() +
  geom_dotplot(binaxis='y',
               binwidth = 0.05,
               stackdir='centerwhole',
               dotsize = 1.2) + 
  scale_y_continuous(breaks=1:10) + 
  labs(x = "Attention",
       y = "Total Correct") +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 26))

sp.plot

sp.anova <- aov(total.correct ~ sp.code, data = first)
summary(sp.anova)
#               Df Sum Sq Mean Sq F value Pr(>F)
# sp.code       2    3.9   1.936   0.396  0.673
# Residuals   460 2248.8   4.889

sp.table <- ddply(first, ~ sp.code, summarise,
                    n.test.total = length(unique(subjectid)), # number tested
                    total.mean = mean(total.correct, na.rm=T), # mean age
                    total.sd = sd(total.correct, na.rm=T)) # sd age


# AGE

ggplot(data = first, aes(x = testing.age, y = total.correct)) +
  geom_jitter() +
  geom_smooth(method = lm, se = T) +
  labs(x = "Age",
       y = "Total Correct") +
  scale_x_continuous(breaks = c(20, 30, 40, 50, 60, 70, 80)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10),
                     limits = c(0, 10)) +
  theme(text = element_text(size = 26))

cor.test(first$testing.age, first$total.correct)
# Pearson's product-moment correlation
# 
# data:  first$testing.age and first$total.correct
# t = -2.4162, df = 461, p-value = 0.01607
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# -0.20091083 -0.02091089
# sample estimates:
# cor 
# -0.111828 

# TEST DAY - this is using ALL DATA (ie. repeats included)

day.plot <- ggplot(twins, aes(x = testing.day, y = total.correct, fill = testing.day)) + 
  geom_violin() +
  geom_dotplot(binaxis='y',
               binwidth = 0.05,
               stackdir='center',
               dotsize = 1.2) + 
  scale_y_continuous(breaks=1:10) + 
  labs(x = "Testing Day",
       y = "Total Correct") +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 26))

day.plot

t.test(twins$total.correct ~ twins$testing.day)
# Welch Two Sample t-test
# 
# data:  twins$total.correct by twins$testing.day
# t = 0.8732, df = 385.22, p-value = 0.3831
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.2164747  0.5623801
# sample estimates:
#   mean in group Saturday   mean in group Sunday 
# 5.809769               5.636816 

day.table <- ddply(twins, ~ testing.day, summarise,
                   n.test.total = length(unique(subjectid)), # number tested
                   total.mean = mean(total.correct, na.rm=T), # mean age
                   total.sd = sd(total.correct, na.rm=T)) # sd age


# TEST DAY - this is using UNIQUE DATA (ie. only first-visits) - same code, different dataset
day.plot <- ggplot(first, aes(x = testing.day, y = total.correct, fill = testing.day)) + 
  geom_violin() +
  geom_dotplot(binaxis='y',
               binwidth = 0.05,
               stackdir='center',
               dotsize = 1.2) + 
  scale_y_continuous(breaks=1:10) + 
  labs(x = "Testing Day",
       y = "Total Correct") +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 26))

day.plot

t.test(first$total.correct ~ first$testing.day)
# Welch Two Sample t-test
# 
# data:  first$total.correct by first$testing.day
# t = -0.05056, df = 100.37, p-value = 0.9598
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.5856910  0.5565796
# sample estimates:
#   mean in group Saturday   mean in group Sunday 
# 5.809769               5.824324 

day.table <- ddply(first, ~ testing.day, summarise,
                   n.test.total = length(unique(subjectid)), # number tested
                   total.mean = mean(total.correct, na.rm=T), # mean age
                   total.sd = sd(total.correct, na.rm=T)) # sd age

######################
# 11-16-2016 - Twin 1 v. Twin 2 correlations ----
######################

# start with unique dataset - "first"

# format and check
table(first$status_selfreported)

# clean up - make "don't know"s into NA
first$status_selfreported <- gsub("Don't know", NA, first$status_selfreported) # remove NA strings
first$status_selfreported <- gsub("Don't Know", NA, first$status_selfreported) # remove NA strings
first$status_selfreported <- gsub("don't know maybe identical", NA, first$status_selfreported) # remove NA strings
table(first$status_selfreported) # 390 MZ, 67 DZ

# giving something a shot from stack overflow
## http://stackoverflow.com/questions/8214303/conditional-replacement-of-values-in-a-data-frame

# copy variable over
first$best.twin.status <- first$status_genotyped

# then ifelse statement: if genotyped is NA, then use self-report. else, use genotyped.
twins.test <- transform(first, best.twin.status = ifelse(!is.na(status_genotyped), status_genotyped, status_selfreported))

# check my work - why is it turning them into numbers?
twins.zygo.check <- twins.test[, c("status_selfreported", "status_genotyped", "best.twin.status")]

# replace number manually - now done
twins.test$best.twin.status <- replace(twins.test$best.twin.status, twins.test$best.twin.status==1, "DZ")
twins.zygo.check <- twins.test[, c("status_selfreported", "status_genotyped", "best.twin.status")]

# OK, now sort by twin ID
## cast out and have duplicate values based on twin ID - need to match by co-twin ID
twins.sort <- twins.test[order(twins.test$twin.id, twins.test$subjectid), ]

# now what if i cut into two different dfs, then merge together?
twin.one <- twins.sort[which(twins.sort$twin.id == 1), ]
twin.two <- twins.sort[which(twins.sort$twin.id == 2), ]

# give it a shot - this merge cuts out a triplet set
twin.pairs <- merge(twin.one, twin.two, by.x = "co.twin.id", by.y = "subjectid")

# cut it down to twin IDs, number correct, and zygosity
colnames(twin.pairs)
twin.pared <- twin.pairs[, c("subjectid", "dob.x", "co.twin.id", "dob.y", "best.twin.status.x", "total.correct.x", "total.correct.y")]


# correlation by zygosity
twin.pared.mz <- twin.pared[which(twin.pared$best.twin.status.x == "MZ"), ]
twin.pared.dz <- twin.pared[which(twin.pared$best.twin.status.x == "DZ"), ]

mz.cor <- tidy(cor.test(twin.pared.mz$total.correct.x, twin.pared.mz$total.correct.y))
row.names(mz.cor) <- "MZ"
dz.cor <- tidy(cor.test(twin.pared.dz$total.correct.x, twin.pared.dz$total.correct.y))
row.names(dz.cor) <- "DZ"

cor.table <- rbind(mz.cor, dz.cor)

# why can't i get r.sq to work?
colnames(cor.table) <- c("correlation", "t", "p.value", "df", "conf.low", "conf.high", "method")
cor.table$r.sq <- cor.table$correlation*cor.table$correlation


ggplot(twin.pared.mz, aes(x = total.correct.x, y = total.correct.y)) +
  geom_jitter() +
  geom_smooth(method = lm, se = F) +
  labs(x = "Twin 1 Correct",
       y = "Twin 2 Correct") +
  ggtitle("MZ twins") +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10),
                     limits = c(0, 10)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10),
                     limits = c(0, 10))

ggplot(twin.pared.dz, aes(x = total.correct.x, y = total.correct.y)) +
  geom_jitter() +
  geom_smooth(method = lm, se = F) +
  labs(x = "Twin 1 Correct",
       y = "Twin 2 Correct") +
  ggtitle("DZ twins") +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10),
                     limits = c(0, 10)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10),
                     limits = c(0, 10))


ggplot(data = data.frame(x = sat$total.correct, y = sun$total.correct), aes(x, y)) +
  geom_jitter() +
  geom_smooth(method = lm, se = T) +
  labs(x = "Saturday Correct",
       y = "Sunday Correct") +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10),
                     limits = c(0, 10)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10),
                     limits = c(0, 10)) +
  theme(text = element_text(size = 26))


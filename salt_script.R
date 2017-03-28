# twinsburg - salt consortium data - 2009 to 2016
# last edit 10-17-16

##############################
# Prep session; read in and format data ----
##############################


setwd("//monell/storage/reed/Labshare/CJ Arayata/Twinsburg/2016 Salt Consortium Progress Report")

# the salt.csv file is just the .csv version of the Access export
salt <- read.csv("salt.csv", header=T, na.strings="")

# i just took this from the PTC sponsors stuff; probably won't be using all packages
library(pacman)
pacman::p_load(reshape2, tm, qdap, plyr, lsr, broom, ggplot2, gridExtra, psych, data.table)

# lower varnames for ease of reference
colnames(salt) <- tolower(colnames(salt))
colnames(salt)

# change dob and test date to actual date format
salt$dob <- as.Date(salt$dob, format = "%m/%d/%Y")
salt$testdate <- as.Date(salt$testdate, format = "%m/%d/%Y")

# extract testing year
salt$testyear <- as.numeric(format(salt$testdate, '%Y'))

# recode races
salt$race <- replace(salt$race, salt$race==1, 'european american')
salt$race <- replace(salt$race, salt$race==4, 'african american')

# check tester names
table(salt$tester)
salt$tester <- tolower(salt$tester)

# check out var names
varnames <- as.data.frame(colnames(salt))


##############################
# 10-17-16: New directives----
## objective: make summary table of n's per year, repeats, sex age race, etc for 2009 - 2016
##############################


# i didn't realize i could calculate a bunch of things in one ddply call
table.1 <- ddply(salt, ~ testyear, summarise,
                 n.test.total = length(unique(subjectid)), # number tested
                 n.white = as.numeric(table(race)[names(table(race)) == 'european american']), # number white
                 perc.white = as.numeric(table(race)[names(table(race)) == 'european american']) / length(unique(subjectid)), # percent white
                 n.female = as.numeric(table(sex)[names(table(sex)) == 'F']), # number female
                 perc.female = as.numeric(table(sex)[names(table(sex)) == 'F']) / length(unique(subjectid)), # percent female
                 age.mean = mean(testing.age, na.rm=T), # mean age
                 age.sd = sd(testing.age, na.rm=T)) # sd age

summary <- function(subjectid, age, sex, race){
  c(n.test.total = length(unique(subjectid)),
    age.mean = mean(age, na.rm=T), # mean age
    age.sd = sd(age, na.rm=T),
    n.white = as.numeric(table(race)[names(table(race)) == 1]),
    perc.white = as.numeric(table(race)[names(table(race)) == 1]) / length(unique(subjectid)),
    n.female = as.numeric(table(sex)[names(table(sex)) == 'F']), # number female
    perc.female = as.numeric(table(sex)[names(table(sex)) == 'F'])
    )
}

# how to calcuate those who did both saturday and sunday - 'duplicated' ??
# first, paste the unique and test year together (e.g. "TW-1 - 2009)
salt$dup.prep <- paste(salt$uniqueid, "-", salt$testyear)

# then, call up (KEEP) duplicates of those instances only (i could not figure out how to use two variables as criteria for duplicates, hence creating one in code above)
duplicates <- salt[salt$dup.prep %in% salt$dup.prep[duplicated(salt$dup.prep)], ]

# export so i can check n's manually
write.csv(duplicates, file="duplicate_checking.csv")

# calculate table of how many participants were repeat customers *within one year*
dup.table <- ddply(duplicates, ~ testyear, summarise,
                   n.test.bothdays = length(unique(uniqueid)))

# merge and rearrange - now complete
full.table <- merge(table.1, dup.table)
full.table <- full.table[, c(1, 2, 9, 3:8)]

write.csv(full.table, file="salt_table.csv")

##############################
# ANALYSES FOR FINAL REPORT????? thursday, 10/13/16----
##############################


# get a quick summary of n's and means - caution: contains repeat customers
salt.summary <- describe(salt)
setDT(salt.summary, keep.rownames=T)

# FIRST VISIT prep:

# sort by date so that i can then take first instance of unique - MUST be done first!
salt.prep <- salt[with(salt, order(testdate)), ]

# now get first instance of unique
first <- salt.prep[!duplicated(salt.prep$uniqueid), ]

# cut race to only white and black - then factor
first.bw <- first[first$race %in% c('european american', 'african american'), ]
class(first.bw$race)
first.bw$race <- factor(first.bw$race)
levels(first.bw$race)

# first.bw is now streamlined first-visit

# NOW do summary means - has all solutions and age for first-timers
first.summary <- describe(first.bw)
setDT(first.summary, keep.rownames=T)

# check n's for sex and race
table(first.bw$sex)
table(first.bw$race)

## solutions of interest:
### water and PTC - controls
### NaCl and KCl - are they related?
### CaCl and CsCl - bonus

# need to check if low-concentration NaCl and KCl show similar distributions to 'full' concentration
## 'normalizing' and combining data?

## ratings of interest:
### liking, bitter, intensity, salitness


# HISTOGRAMS - loop through all possible solution/taste quality graphs


# first, messing with indexing to get ylab to work below - returns 'n' from salt.summary table based on matching of solution variable
as.numeric(salt.summary[which(salt.summary$rn == 'nacl_like3'), ]$n) # returns 3924
as.numeric(salt.summary[which(salt.summary$rn == colnames(first.bw[37])), ]$n) # returns 3924


## i wish i knew a way to store each graph...

for (phenotype in 15:ncol(first.bw)){
  graph <- ggplot(data=first.bw,
                  aes(x=first.bw[,phenotype])) +
    geom_histogram(binwidth=0.2) +
    xlab(colnames(first.bw)[phenotype]) +
    ylab(paste("Frequency (total n:", as.numeric(salt.summary[which(salt.summary$rn == colnames(first.bw[phenotype])), ]$n), ")")) +
  ggtitle(paste("Histogram of", as.character(colnames(first.bw)[phenotype])))
  print(graph)
}


# what if i faceted each graph by sex (and/or race)?
for (phenotype in 15:15){
  graph <- ggplot(data=first.bw,
                  aes(x=first.bw[,phenotype])) +
    geom_histogram(binwidth=0.2) +
    xlab(colnames(first.bw)[phenotype]) +
    ylab(paste("Frequency (total n:", as.numeric(salt.summary[which(salt.summary$rn == colnames(first.bw[phenotype])), ]$n), ")")) +
    ggtitle(paste("Histogram of", as.character(colnames(first.bw)[phenotype]))) +
    facet_wrap(~ sex)
  print(graph)
}


# how about some scatterplots?
ggplot(first.bw, aes(x = nacl_salty3, y = x0.75_nacl_salty20)) +
  geom_point() +
  geom_smooth(method=lm, se=F)

x <- first.bw$nacl_salty3
y <- first.bw$x0.75_nacl_salty20
tidy(cor.test(x, y)) # r = 0.32

summary(aov(x0.75_nacl_salty20 ~ nacl_salty3, data=first.bw)) # gives same p-value as correlation

ggplot(first.bw, aes(x = kcl_salty10, y = x0.75_kcl_salty21)) +
  geom_point() +
  geom_smooth(method=lm, se=F)


ggplot(first.bw, aes(x = nacl_salty3, y = kcl_salty10)) +
  geom_point() +
  geom_smooth(method=lm, se=F)

x <- first.bw$nacl_salty3
y <- first.bw$kcl_salty10
cor.test(x, y)

ggplot(first.bw, aes(x = nacl_like3, y = kcl_like10)) +
  geom_point() +
  geom_smooth(method=lm, se=F)

x <- first.bw$nacl_like3
y <- first.bw$kcl_like10
cor.test(x, y)


# heatmaps - what exactly?

# reliability - test-retest
## can i take a leaf out of my mennella analyses, and graph longitudinal trends for each DV, controlling for repeated measures?
## or do it 'old school' way where i split dataset into saturday-sunday, and correlate


# zygosity
# need to split original dataset into twin 1 and twin 2 - then correlate
## need to make sure that twins LINE UP in both datasets
## then, can i run correlations between, say, PTC bitter (twin 1) and PTC bitter (twin 2)?


# BATCH UNCOMMENT by highlighting section, then CTRL + SHIFT + C

# junk stuff to fix later
# collapse twin status variables - this doesn't work!
# copy genotyped status into new 'twin_status' variable
# salt$twin_status <- salt$status_genotyped
# 
# # now, if twin_status is NA, set value to self-reported status
# salt$twin_status[is.na(salt$twin_status)] <- salt$status_selfreported

# reorganize
# salt <- salt[, c(1:7, 51, 8:50)]

# cbind.fill <- function(...){
#   nm <- list(...)
#   nm <- lapply(nm, as.matrix)
#   n <- max(sapply(nm, nrow))
#   do.call(cbind, lapply(nm, function (x)
#     rbind(x, matrix(, n-nrow(x), ncol(x)))))
# }

# # various subsets - might not be needed
# sodium <- subset(salt, !is.na(nacl_bitter3)) # only choosing bitter because i dont know fuller criteria
# sodium.low <- subset(salt, !is.na(x0.75_nacl_bitter20))
# potassium <- subset(salt, !is.na(kcl_bitter10))
# potassium.low <- subset(salt, !is.na(x0.75_kcl_bitter21))
# calcium <- subset(salt, !is.na(cacl2_bitter11))
# cesium <- subset(salt, !is.na(cscl_bitter8))
# lpar <- subset(salt, !is.na(lpar1agonistbitter9)) # only in 2015
# twin.2016 <- subset(salt, testdate == '2016-08-06' | testdate == '2016-08-07') # just NaCl




#make sex blanks NAs
  
salt$sex[salt$sex=="BLANK"] <- NA
salt$sex <- factor(salt$sex)
table.sex <- as.data.frame(table(salt$testyear, salt$sex))


# 11-16-16 - checking if zygosity collapse works
# check
table(salt$status_selfreported)

# clean up - make "don't know"s into NA - there's gotta be a better way to do this
salt$status_selfreported <- strip(salt$status_selfreported) # takes out the () in P(A)V etc.
salt$status_selfreported <- replace(salt$status_selfreported, salt$status_selfreported=="", NA)
salt$status_selfreported <- gsub("dont know", NA, salt$status_selfreported) # remove NA strings
salt$status_selfreported <- gsub("dont know maybe identical", NA, salt$status_selfreported) # remove NA strings
salt$status_selfreported <- gsub("likely mz not sure", "mz", salt$status_selfreported) # remove NA strings
salt$status_selfreported <- gsub("mz mirror", "mz", salt$status_selfreported) # remove NA strings
salt$status_selfreported <- gsub("mzmirror", "mz", salt$status_selfreported) # remove NA strings
salt$status_selfreported <- gsub("both", NA, salt$status_selfreported) # kick triplets, sorry
salt$status_selfreported <- gsub("not sure most likely mz", "mz", salt$status_selfreported) # remove NA strings

# giving something a shot from stack overflow - kinda works
## http://stackoverflow.com/questions/8214303/conditional-replacement-of-values-in-a-data-frame

# copy variable over
salt$best.twin.status <- salt$status_genotyped

# then ifelse statement: if genotyped is *not* NA, then use genotype. else, use self-report.
salt.test <- transform(salt, best.twin.status = ifelse(!is.na(status_genotyped), status_genotyped, status_selfreported))

# check my work - why does it change it into numbers?
twins.zygo.check <- salt.test[, c(6, 7, 52)]

# ugh change the numbers into appropriate things, ugh
salt.test$best.twin.status <- replace(salt.test$best.twin.status, salt.test$best.twin.status==1, NA)
salt.test$best.twin.status <- replace(salt.test$best.twin.status, salt.test$best.twin.status==2, "dz")
salt.test$best.twin.status <- replace(salt.test$best.twin.status, salt.test$best.twin.status==3, "mz")



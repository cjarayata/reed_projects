# Burly1 R script - 2017-1-18 #

# Adapted from the Adip20 developed by Charles J Arayata, Danielle R Reed and Cailu Lin
# Monell Chemical Senses Center

## Prepping the session: Packages and Data Files ## ----

#use 'pacman' package to install and load packages needed
if (!require("pacman")) install.packages("pacman", "CRAN")
library(pacman)
pacman::p_load(psych, lubridate, plyr, dplyr, broom, reshape, data.table, lsr)

#Read in data
Data_burly1 <- read.csv("Data_burly1.csv", header=TRUE)

#Read in SNP info
snp <- read.csv("snp.csv", header=TRUE)
snp$SNPs <- tolower(snp$SNPs)


## Table 1 ## ----
## Characteristics of N=2,074 male mice used in burly1 mapping studies ##


#start with n's for each strain
n1 <- summary(Data_burly1$Table1_Mapping_population)
n1 <- data.frame(n1)
setDT(n1, keep.rownames=TRUE)
names(n1) <- c("Table1_Mapping_population", "n")

#get rest of descriptives
tbl1.min <- ddply(Data_burly1, "Table1_Mapping_population", summarize, min=min(Days))
tbl1.max <- ddply(Data_burly1, "Table1_Mapping_population", summarize, max=max(Days))
tbl1.mean <- ddply(Data_burly1, "Table1_Mapping_population", summarize, mean=mean(Days))
tbl1.sd <- ddply(Data_burly1, "Table1_Mapping_population", summarize, sd=sd(Days))

#merge together
table.1 <- Reduce(function(...) merge(..., by="Table1_Mapping_population", all=TRUE), list(n1, tbl1.min, tbl1.max, tbl1.mean, tbl1.sd))

#need to get START and END dates
#first subset data to cut out genotypes; then make all variable names lowercase for ease of variable reference
cj.data1 <- Data_burly1[ , 1:15]
colnames(cj.data1) <- tolower(names(cj.data1))

#parse out birthdate into new variable - tacks onto end of data frame (which is why we just created a new data frame)
cj.data1$birthdate.fix <- parse_date_time(cj.data1$birth.date, orders="mdy")

#get min and max dates for each strain
t1.bday.min <- ddply(cj.data1, "table1_mapping_population", summarize, start=min(birthdate.fix))
t1.bday.max <- ddply(cj.data1, "table1_mapping_population", summarize, end=max(birthdate.fix))
t1.bday <- merge(t1.bday.min, t1.bday.max)

#merge descriptives and start/end dates together - now complete
table.1 <- merge(table.1, t1.bday, by.x="Table1_Mapping_population", by.y="table1_mapping_population", all=TRUE)
write.table(table.1, file="Table1.csv", sep=",", row.names=F)
## Table 1 END ##



## Table 2 ## ----
## Characteristics of N=1,293 congenic mice by N=22 strains ##

#for n's
n2 <- summary(Data_burly1$Table_2_Substrains)
n2 <- data.frame(n2)
setDT(n2, keep.rownames=TRUE)
names(n2) <- c("Table_2_Substrains", "n")

#for 'Day' descriptives
tbl2.min <- ddply(Data_burly1, "Table_2_Substrains", summarize, min=min(Days))
tbl2.max <- ddply(Data_burly1, "Table_2_Substrains", summarize, max=max(Days))
tbl2.mean <- ddply(Data_burly1, "Table_2_Substrains", summarize, mean=mean(Days))
tbl2.sd <- ddply(Data_burly1, "Table_2_Substrains", summarize, sd=sd(Days))

#merge all together
table.2 <- Reduce(function(...) merge(..., by="Table_2_Substrains", all=TRUE), list(n2, tbl2.min, tbl2.max, tbl2.mean, tbl2.sd))

write.table(table.2, file="Table2.csv", sep=",", row.names=F)
## Table 2 END ##


## Figure 3, 5, S3 Figure- Behemoth ## ----
#load in data
cj.figure2 <- Data_burly1
colnames(cj.figure2) <- tolower(names(cj.figure2))
names(cj.figure2)

## below is ORIGINAL code used to produce p-values
#create table first
rm(aov.out2)
##F2_First
aov.out2 <- tidy(aov(lm(bodyweight ~ cj.figure2[,22]+age, data=cj.figure2, subset=table1_mapping_population=='F2')))
aov.out2["SNP"] <- colnames(cj.figure2)[22]
aov.out2["Strain"] <- "F2"

for(i in 23:length(names(cj.figure2))){
  temp.f2<-tidy(aov(lm(bodyweight~ cj.figure2[,i]+age, data=cj.figure2, subset=table1_mapping_population=='F2')))
  temp.f2["SNP"] <- colnames(cj.figure2)[i]
  temp.f2["Strain"] <- "F2"
  aov.out2 <- rbind(aov.out2, temp.f2)
}
write.table(aov.out2, file="F2__association.csv", sep=",", row.names=F)

#####
rm(aov.out2)
##F2_Second

aov.out2 <- tidy(aov(lm(lean ~ cj.figure2[,22]+total, data=cj.figure2, subset=table1_mapping_population=='F2_second')))
aov.out2["SNP"] <- colnames(cj.figure2)[22]
aov.out2["Strain"] <- "F2_HF"

for (i in c(33:51, 71:80, 114:133)){  
  temp.f2.hf<-tidy(aov(lm(lean~cj.figure2[,i] + total, data=cj.figure2, subset=table1_mapping_population=='F2_second')))
  temp.f2.hf["SNP"] <- colnames(cj.figure2)[i]
  temp.f2.hf["Strain"] <- "F2_HF"
  aov.out2 <- rbind(aov.out2, temp.f2.second)
}
write.table(aov.out2, file="F2_second_association.csv", sep=",", row.names=F)

########
  rm(aov.out2)
###Backcross_129
aov.out2 <- tidy(aov(lm(lean ~ cj.figure2[,22]+total, data=cj.figure2, subset=table1_mapping_population2=='Backcross_129')))
aov.out2["SNP"] <- colnames(cj.figure2)[22]
aov.out2["Strain"] <- "Backcross_129"
for (i in 23:(length(names(cj.figure2))-1)){ 
  temp.Backcross_129<-tidy(aov(lm(lean~cj.figure2[,i] + total, data=cj.figure2, subset=table1_mapping_population2=='Backcross_129')))
  temp.Backcross_129["SNP"] <- colnames(cj.figure2)[i]
  temp.Backcross_129["Strain"] <- "Backcross_129"
  aov.out2 <- rbind(aov.out2, temp.Backcross_129)
}
write.table(aov.out2, file="Backcross_129_withoutsonsomic_association1.csv", sep=",", row.names=F)


rm(aov.out2)
###Backcross_B6
aov.out2 <- tidy(aov(lm(lean ~ cj.figure2[,22]+total, data=cj.figure2, subset=table1_mapping_population2=='Backcross_B6')))
aov.out2["SNP"] <- colnames(cj.figure2)[22]
aov.out2["Strain"] <- "Backcross_B6"
for (i in 23:(length(names(cj.figure2))-1)){ 
  temp.Backcross_B6<-tidy(aov(lm(lean~cj.figure2[,i] + total, data=cj.figure2, subset=table1_mapping_population2=='Backcross_B6')))
  temp.Backcross_B6["SNP"] <- colnames(cj.figure2)[i]
  temp.Backcross_B6["Strain"] <- "Backcross_B6"
  aov.out2 <- rbind(aov.out2, temp.Backcross_B6)
}
write.table(aov.out2, file="Backcross_B6_withoutN7_association1.csv", sep=",", row.names=F)

rm(aov.out2)

#######
  
rm(aov.out2)


##For 129.B6_chr2
aov.out2 <- tidy(aov(lm(lean ~ cj.figure2[,21]+sex+total, data=cj.figure2, subset=table1_mapping_population=='129.B6-Chr2')))
aov.out2["SNP"] <- colnames(cj.figure2)[21]
aov.out2["Strain"] <- "129.B6_Chr2"
  
for (i in 22:length(names(cj.figure2))){ 
  temp.129.B6_Chr2<-tidy(aov(lm(lean~cj.figure2[,i]+sex+total, data=cj.figure2, subset=table1_mapping_population=='129.B6-Chr2')))
  temp.129.B6_Chr2["SNP"] <- colnames(cj.figure2)[i]
  temp.129.B6_Chr2["Strain"] <- "129.B6-Chr2"
  aov.out2 <- rbind(aov.out2, temp.129.B6_Chr2)
}
write.table(aov.out2, file="129.B6-Chr2_association1.csv", sep=",", row.names=F) 

 
#######
rm(aov.out2)
###Congenics
aov.out2 <- tidy(aov(lm(log(lean) ~ cj.figure2[,60]+total, data=cj.figure2, subset=table1_mapping_population2=='Congenic')))
aov.out2["SNP"] <- colnames(cj.figure2)[60]
aov.out2["Strain"] <- "Congenic"
for (i in 61:length(names(cj.figure2))){ 
  temp.Congenic<-tidy(aov(lm(log(lean)~cj.figure2[,i] + total, data=cj.figure2, subset=table1_mapping_population2=='Congenic')))
  temp.Congenic["SNP"] <- colnames(cj.figure2)[i]
  temp.Congenic["Strain"] <- "Congenic"
  aov.out2 <- rbind(aov.out2, temp.Congenic)
}
write.table(aov.out2, file="Congenic_v2_association1.csv", sep=",", row.names=F)


###Figure 4####
#CJA method
#note use of 'aov' - default is Type I (Sequential) Sums of Squares
b6_chr9 <- Data_near_final[which(Data_near_final$Groups=='B6.129-Chr2' | Data_near_final$Groups=='C57BL/6ByJ'), ]
cj.model1 <- aov(lean ~ Groups + BW, data=b6_chr2)
summary(cj.model1)

##For 129.B6_chr2
aov.out2 <- tidy(aov(lm(lean ~ cj.figure2[,21]+sex+total, data=cj.figure2, subset=table1_mapping_population=='129.B6-Chr2')))
aov.out2["SNP"] <- colnames(cj.figure2)[21]
aov.out2["Strain"] <- "129.B6_Chr2"
  
for (i in 22:length(names(cj.figure2))){ 
  temp.129.B6_Chr2<-tidy(aov(lm(lean~cj.figure2[,i]+sex+total, data=cj.figure2, subset=table1_mapping_population=='129.B6-Chr2')))
  temp.129.B6_Chr2["SNP"] <- colnames(cj.figure2)[i]
  temp.129.B6_Chr2["Strain"] <- "129.B6-Chr2"
  aov.out2 <- rbind(aov.out2, temp.129.B6_Chr2)
}
write.table(aov.out2, file="129.B6-Chr2_association1.csv", sep=",", row.names=F) 

## analyses - consomics END ##

##S4 Figure for Fat
load in data
cj.figure2 <- Data_burly1
colnames(cj.figure2) <- tolower(names(cj.figure2))
names(cj.figure2)


#######
rm(aov.out2)
###Congenics
aov.out2 <- tidy(aov(lm(fat ~ cj.figure2[,60]+total, data=cj.figure2, subset=table1_mapping_population2=='Congenic')))
aov.out2["SNP"] <- colnames(cj.figure2)[60]
aov.out2["Strain"] <- "Congenic"
for (i in 61:length(names(cj.figure2))){ 
  temp.Congenic<-tidy(aov(lm(fat~cj.figure2[,i] + total, data=cj.figure2, subset=table1_mapping_population2=='Congenic')))
  temp.Congenic["SNP"] <- colnames(cj.figure2)[i]
  temp.Congenic["Strain"] <- "Congenic"
  aov.out2 <- rbind(aov.out2, temp.Congenic)
}
write.table(aov.out2, file="Congenic_fat_association1.csv", sep=",", row.names=F)


### S1 Figure ggplot2 for correlation within each population
library("ggplot2")

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


Large <- read.csv("corr_burly1.csv")
Large<- Large[-which(is.na(Large[,13])), ]
p1 <- ggplot(Large, aes(Lean_by_DEXA,  Lean_by_MR)) + geom_point(shape=1)+ geom_smooth(method=lm) + facet_grid(. ~ Strain) + labs (x="Lean by DEXA, g", y="Lean by MR, g")
#p1 <- ggplot(Large, aes(Lean_by_DEXA,  Lean_by_MR)) + geom_point(shape=1, aes(fill=as.factor(variable))+ geom_smooth(method=lm) + facet_grid(. ~ Strain) + labs (x="Lean by DEXA, g", y="Lean by MR, g")+scale_fill_manual(values=c("red","blue"),guide=guide_legend(title="Variable"))

Large1 <- read.csv("Large1corr_burly1_2.csv")
Large1<- Large1[-which(is.na(Large1[,11])), ]

p2 <- ggplot(Large1, aes(Lean_body_mass,Body_weight, color=method)) + geom_point(shape=1)+geom_smooth(method=lm) + facet_grid(. ~ Strain) + labs(x="Lean body mass, g", y="Body weight, g")

##p2 <- ggplot(Large1, aes(Lean_body_mass,Body_weight, color=method)) + geom_point(shape=1)+geom_smooth(method=lm) + facet_grid(. ~ Strain) + theme(legend.position=c(1, 5),legend.justification=c(1,5))+labs(x="Lean body mass, g", y="Body weight, g")

ggdraw() +
  draw_plot(p1, 0, .5, .8, .5) +
  draw_plot(p2, 0, 0, 1, .5)+ 
   draw_plot_label(c("A", "B"), c(0, 0), c(1, .5))
multiplot(p1,p2)

##another option to arrange the plots
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
library(ggpubr)
ggarrange(ggarrange(p1,  ncol = 2, labels = c(""), widths=c(1, 0.9)), p2, nrow = 2, labels = c("A","B") ) 


###S4 Table correlation within each mapping population
DEXA=subset(Large, method=="DEXA")
a=subset(DEXA, Strain=="129.B6-Chr2")
cor.test(a$Lean_body_mass,a$Body_weight, method="pearson" )


b=subset(DEXA, Strain=="B6.129-Burly1")
cor.test(b$Lean_body_mass,b$Body_weight, method="pearson" )


c=subset(DEXA, Strain=="B6.129-Chr2")
cor.test(c$Lean_body_mass,c$Body_weight, method="pearson" )

d=subset(DEXA, Strain=="C57BL/6ByJ")
cor.test(d$Lean_body_mass,d$Body_weight, method="pearson" )

e=subset(DEXA, Strain=="F1 x 129")
cor.test(e$Lean_body_mass,e$Body_weight, method="pearson" )

f=subset(DEXA, Strain=="F1 x B6")
cor.test(f$Lean_body_mass,f$Body_weight, method="pearson" )

g=subset(DEXA, Strain=="F2_Second")
cor.test(g$Lean_body_mass,g$Body_weight, method="pearson" )


###Figure 2 3-d-correlation
library("scatterplot3d")
setwd("D:/Data/B6.129/B6.129-congenics/B6.129-Burly1/ggplot2")
dat=read.csv("3dcorrelation.csv")
scatterplot3d(dat, xlim=c(15,45), ylim=c(15,45), zlim=c(15,45), xlab="Lean by DEXA, g", zlab="Body weight, g", ylab="Lean by MR, g", cex.lab = 2, cex.axis = 1.2)
savePlot(filename ="Figure_2_3d.tif",type ="tiff", device = dev.cur())

##Burly1 end###
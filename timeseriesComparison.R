#mAyo Time series Dpm vs dpnm analysis. Or cor/Inc analysis. 
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
## Norms the data within specified groups in a data frame; it normalizes each
## subject (identified by idvar) so that they have the same mean, within each group
## specified by betweenvars.
##   data: a data frame.
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   na.rm: a boolean that indicates whether to ignore NA's
normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
  library(plyr)
  
  # Measure var on left, idvar + between vars on right of formula.
  data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
                         .fun = function(xx, col, na.rm) {
                           c(subjMean = mean(xx[,col], na.rm=na.rm))
                         },
                         measurevar,
                         na.rm
  )
  
  # Put the subject means with original data
  data <- merge(data, data.subjMean)
  
  # Get the normalized data in a new column
  measureNormedVar <- paste(measurevar, "_norm", sep="")
  data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
    mean(data[,measurevar], na.rm=na.rm)
  
  # Remove this subject mean column
  data$subjMean <- NULL
  
  return(data)
}
## Summarizes data, handling within-subjects variables by removing inter-subject variability.
## It will still work if there are no within-S variables.
## Gives count, un-normed mean, normed mean (with same between-group mean),
##   standard deviation, standard error of the mean, and confidence interval.
## If there are within-subject variables, calculate adjusted values using method from Morey (2008).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   withinvars: a vector containing names of columns that are within-subjects variables
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  
  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
                       FUN=is.factor, FUN.VALUE=logical(1))
  
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  
  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL
  
  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)
  
  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")
  
  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                                  FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
  
  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor
  
  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}

library("broom", lib.loc="~/R/R-3.3.2/library")
library("tidyr", lib.loc="~/R/R-3.3.2/library")
library("ggplot2", lib.loc="~/R/R-3.3.1/library")
library("gridExtra", lib.loc="~/R/R-3.3.1/library")
library("plyr")
library("stringr")
library("dplyr", lib.loc="~/R/R-3.3.1/library")
library("data.table", lib.loc="~/R/R-3.3.2/library")

stud <- "" #this one for test phase
#stud <- "Study_" #this one if doing study phase sanity check

#vector for different comparisons (DPM and left/right brain)
#Change which is commented out for DPM or cor/inc analysis
#vects <- c("DPM_l", 'DPM_r',"DPNM_l", 'DPNM_r')
vects <- c("correct_l", 'correct_r',"incorrect_l", 'incorrect_r')


# ROIs <- c("acc","ant_hipp","AT","caud","erc","fusi","MD","pcc","phc","post_hipp","prc","precun") # "all_hipp_"
 ROIs <- c("ant_hipp","AT","erc","fusi","MD","phc","post_hipp","prc")
 # ROIs <- c("prc") # "all_hipp_"
# ROIs <- c("post_hipp") # "all_hipp_"
# ROIs <- c("phc") # "all_hipp_"
# ROIs <- c("MD") # "all_hipp_"
#ROIs <- c("fusi") # "all_hipp_"
# ROIs <- c("erc") # "all_hipp_"
# ROIs <- c("AT") # "all_hipp_"
# ROIs <- c("ant_hipp") # "all_hipp_"



#Subjects applying code to. Right now skipping subjects 3,6,11, 18 bc. don't have ET data or have error)
y <- c(seq(from = 1, to = 2, by =1),seq(from = 4, to = 5, by =1),seq(from = 7, to = 10, by =1),seq(from = 12, to = 17, by =1),seq(from = 19, to = 28, by =1)) #Subjects applying code to. Right now skipping subjects 3,6,11, 18 bc. don't have ET data or have error)

#making empty dataframe to fill with loop
All <- data.frame(matrix(ncol = 5, nrow = 0))

#Double for loop (subj# and condition) to load in files and add to 1 dataframe
mylistnew <- list()
for (i in y){
  for (t in vects){
  for (q in ROIs){
    b = i + 1000
  b
  data_dir = 'C:/Users/kgeier/Desktop/SuperTeamROIoutput' 
  filename1 = paste0(stud,"timeseries_",t, "_", q,"_")
  filename2 =  b
  fileextension = ".csv"
  filename = paste0(filename1,filename2, fileextension)
  filepath = paste(data_dir, '/', filename, sep = ''); filepath
  
  df.21 <- fread(filepath)
  
  #Part that adds SID and Condition info 
  SID <- rep(b,21)
  Condition <- rep(t,21)
  Region <- rep(q,21)
  SIDCon <- cbind(SID,Condition, Region)
  SIDCon <- as.data.frame(SIDCon)
  
  withSIDCon <- cbind(SIDCon, df.21)
  
  All <- rbind(All,withSIDCon, fill = F)
  }
  }
}

#Renaming columns, removing NA fills, and adding TP list for each of the 21 seconds
df <- All[,c("SID", "Condition", "Region","V1")]
df1 <- subset(df, subset = is.na(df$SID)==F)
Name_of_TP <- rep(seq(1:21),96) #this part is done semi manually so double check each subject/condition has exactly 21 time points etc. 
df2 <- cbind(df1,Name_of_TP)

df2
df2$SID <- factor(df2$SID)
detach("package:dplyr", unload=TRUE)
df2_stat <- summarySEwithin(df2, measurevar = "V1", withinvars = c("Condition","Region","Name_of_TP"), idvar = c("SID"),na.rm = FALSE, conf.interval = .95)

head(df2_stat)


library("dplyr", lib.loc="~/R/R-3.3.1/library")
#Using ddply to summarize mean, sd, se of the data for graphing. 
TPdata <- ddply(df2, c("Condition", "Region", "Name_of_TP"), summarise,
               N    = length(Condition),
               mean = mean(V1),
               sd   = sd(V1),
               se   = sd / sqrt(N)
)

head(TPdata)

#Disecting the condition column (for within subject error bar version) so can colour by DPM_DPNM and facet by side of brain
df2_stat$brainhemi <- str_sub(df2_stat$Condition,start = -1)
df2_stat$Match.NonMatch <- str_sub(df2_stat$Condition,end = -3)

#Disecting the condition column so can colour by DPM_DPNM and facet by side of brain
TPdata$brainhemi <- str_sub(TPdata$Condition,start = -1)
TPdata$Match.NonMatch <- str_sub(TPdata$Condition,end = -3)

# #Making x axis numeric so jitter will work
# TPdata$Name_of_TP <- as.numeric(as.character(TPdata$Name_of_TP))
head(df2_stat)
#making graph 
ggplot(df2_stat, aes(x=Name_of_TP, y=V1, group=Condition, 
                   colour=Match.NonMatch)) +
  geom_line(size=1.5) +
  geom_errorbar(width=0, aes(ymin=V1-se, ymax=V1+se), colour="red", size=1.5) +
  geom_errorbar(width=0, aes(ymin=V1-se, ymax=V1+se), data=df2_stat, size=1.5) +
  geom_point(size=2) + 
  scale_color_brewer(palette='Set2')+ 
  scale_x_discrete(limits=seq(1:21)) +
  facet_wrap(~Region~brainhemi)+
  xlab("Time (sec)") +
  ylab("Percent Signal Change?")

#between subject not within subject error bars
ggplot(TPdata, aes(x=Name_of_TP, y=mean, group=Condition, 
                     colour=Match.NonMatch)) +
  geom_line(size=1.5) +
  geom_errorbar(width=0, aes(ymin=mean-se, ymax=mean+se), colour="red", size=1.5) +
  geom_errorbar(width=0, aes(ymin=mean-se, ymax=mean+se), data=TPdata, size=1.5) +
  geom_point(size=2) + 
  scale_color_brewer(palette='Set2')+ 
  scale_x_discrete(limits=seq(1:21)) +
  facet_wrap(~Region~brainhemi)+
  xlab("Time (sec)") +
  ylab("Percent Signal Change?")

check <- ( (TPdata$se-df2_stat$se )/ TPdata$se )*100
range(check)
mean(check)

# #renaming columns so 1-21 seconds are now timepoint 1-21 so can use dplyr on them (inshallah)
#df5 <- spread(df2, Name_of_TP, V1)
# namevec <- "Tp1"
# for (q in 2:21){
#   w<- paste0("Tp",q)
#   namevec <- c(namevec,w)
# }
# colnames(df5) <- c("SID", "Condition", namevec)
# 
# tapply(df5$, df$group, summary)

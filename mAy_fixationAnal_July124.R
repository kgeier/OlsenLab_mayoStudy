#mAy fixation Analysis version 2
 
library("broom", lib.loc="~/R/R-3.3.2/library")
library("tidyr", lib.loc="~/R/R-3.3.2/library")
library("ggplot2", lib.loc="~/R/R-3.3.1/library")
library("gridExtra", lib.loc="~/R/R-3.3.1/library")
library("plyr")
library("dplyr", lib.loc="~/R/R-3.3.1/library")
library("data.table", lib.loc="~/R/R-3.3.2/library")
library(lme4)

#Find How_to document either in the Dropbox\OlsenLabDocs\Useful Documents\How-Tos folder called README_RcodeMRIonsets
#Will also be copied at the bottom of this code as a comment. 


#GOALS OF THIS R CODE
#0. Loading data, renaming, and removing null trials
#1. Clean data so only trials with eyes open/functioning data remain
#2. Hanulla analysis: which trials get looked at more. IF only look at those for MRI does it show up?
#More specifically:
#a) bin the viewing by 500 ms (0-3000)
#b) Anova to check if any bins have significantly more viewing of correctly chosen faces than incorrectly chosen faces
#c) Divide into 2 categories of trial. 1)-DPMatch at least 10% more viewing of matched (matched means target face) face than next highest face. Compare DPM to those with 10% higher viewing to WRONG face. Equal viewing trials are ignored (so partly comparing to "inappropriate memory" not to "no memory")
#d) Potentially add that a new column to onset files and run more GLMs. Test if hippo activity at cue is proportional to DPM vs DPNM
#e) hannula got 62% cor, 25% inc, 12% "dont know". 
#f) Functional connectivity shows correlation with explicit correct answers not eye-mvmts. So test this as well.
#g) Hannula said a cool 'future analysis' is compare hippo activity w/ DPM correct vs DPNM incorrect. 
#3. Repeat of IA analysis. Do they look more at x, does fixation predict accuracy? 


#0 loading Steph's scripts for SE from Objfam
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
  datac$ci<- datac$se * ciMult
  
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


#0)############################################################################################################################################################################################################
#import csv for eye-tracking. Don't forget to change the file name and path if you're using a different set of data!
df_study_IAreport <- fread("C:/Users/kgeier/Google Drive/OlsenLab_cloud/All Experiments/MAYO/R studio analysis/1to28_studyreport_fixationAnalysis.csv")
#uncomment if want full test IP rather than just IP with 3 faces on screen# df_test_IAreport <- fread("C:/Users/kgeier/Google Drive/OlsenLab_cloud/All Experiments/MAYO/R studio analysis/1to25_testreport_fixationAnalysis.csv")
df_test_IAreport <- fread("C:/Users/kgeier/Google Drive/OlsenLab_cloud/All Experiments/MAYO/R studio analysis/1to28_test3FacesReport_fixationAnalysis.csv")

#P1 -Correct IPs ####
#following method will take these two dataframes, add on a column with ordered numbers, 
# subset them so that study phase IP just has study phase trials and that the df with test IP just has test trials,
# then add them together and resort with the ordered numbers from original. This will get the correct fixation data 
# for each trial but combined appropriately.  See practice below for a simpler idea of what is going on. 

#appending column with seq 1:length of rows
df_study_withcount <- cbind(df_study_IAreport, seq(1:nrow(df_study_IAreport)))
df_test_withcount <-cbind(df_test_IAreport, seq(1:nrow(df_test_IAreport)))

#making sure the new column is named 
names(df_study_withcount)[names(df_study_withcount) == 'V2'] <- 'orderer'
names(df_test_withcount)[names(df_test_withcount) == 'V2'] <- 'orderer'

#and numeric
df_study_withcount$orderer <- as.numeric(as.character(df_study_withcount$orderer))
df_test_withcount$orderer <- as.numeric(as.character(df_test_withcount$orderer))

#subsetting 
subset_study <-  df_study_withcount[grepl("1", df_study_withcount$task),]
subset_test <-  df_test_withcount[grepl("2", df_test_withcount$task),]

#rbind 
df_rebound <- rbind(subset_study, subset_test, fill = TRUE) 

#re-ordering
df_final <- df_rebound[order(df_rebound$orderer),]

#Cleaning trials and blocks manually that have been noted as bad. ####
#getting rid of blocks that didn't trigger, etc.
df_final <- df_final[df_final$RECORDING_SESSION_LABEL != "may14kg6"]
df_final <- df_final[df_final$RECORDING_SESSION_LABEL != "may05kg4"]

#changing TRIAL_LABEL from Trial: 19 to 19. 
df_final$TRIAL_LABEL <- gsub("Trial: ", "", df_final$TRIAL_LABEL) 

###- REMOVING 7 trials (both study and test) that had innapropriately had a mix of older/yonger/male/female in the lures. Ex. OF11, YF15, OF23. 
#trials were determined to be 
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may16kg1" | df_final$`TRIAL_LABEL` !="15",)
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may16kg4" | df_final$`TRIAL_LABEL` !="42",)
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may18kg3" | df_final$`TRIAL_LABEL` !="25",)
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may04kg1" | df_final$`TRIAL_LABEL` !="19",)
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may06kg3" | df_final$`TRIAL_LABEL` !="35",)
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may11kg4" | df_final$`TRIAL_LABEL` !="24",)
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may23kg4" | df_final$`TRIAL_LABEL` !="36",)
#removing an 8th and 9th trial b.c. lure face matched target face
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may22kg6" | df_final$`TRIAL_LABEL` !="9",)
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may13kg6" | df_final$`TRIAL_LABEL` !="19",)

df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may16kg1" | df_final$`TRIAL_LABEL` !="60",)
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may16kg4" | df_final$`TRIAL_LABEL` !="54",)
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may18kg3" | df_final$`TRIAL_LABEL` !="55",)
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may04kg1" | df_final$`TRIAL_LABEL` !="64",)
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may06kg3" | df_final$`TRIAL_LABEL` !="60",)
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may11kg4" | df_final$`TRIAL_LABEL` !="50",)
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may23kg4" | df_final$`TRIAL_LABEL` !="53",)
#removing an 8th and 9thtrial b.c. lure face matched target face
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may22kg6" | df_final$`TRIAL_LABEL` !="62",)
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may13kg6" | df_final$`TRIAL_LABEL` !="53",)

#adding Subject ID to both study and test. 
df_final$SLabNoBlock <- df_final$RECORDING_SESSION_LABEL #makes a copy of session label to get SID from
df_final$SLabNoBlock <- substring(df_final$SLabNoBlock,1, 7) #selects first 7 charachters to remove block number
df_final$SID <- gsub("[^0-9]","",df_final$SLabNoBlock) #removes all non-numbers
df_final$SID <- gsub("146","14",df_final$SID) #subject 14 has issue with not being same #char. So this fixes it.
df_final$SID <- paste0("10",df_final$SID)
#still not working... use lines below   df_final$SID <- ifelse(grepl("may05bk", df_final$RECORDING_SESSION_LABEL), gsub("1005 ", "2001", df_final$SID), df_final$SID <- df_final$SID)
#replace next 6 lines with regular expressions when you have more time. 
df_final[df_final$SID == "1005" & df_final$RECORDING_SESSION_LABEL == "may05kg1", "SID"] <- "2001" #subject5a should be 2001 (artist).
df_final[df_final$SID == "1005" & df_final$RECORDING_SESSION_LABEL == "may05kg2", "SID"] <- "2001" #subject5a should be 2001 (artist).
df_final[df_final$SID == "1005" & df_final$RECORDING_SESSION_LABEL == "may05kg3", "SID"] <- "2001" #subject5a should be 2001 (artist).
df_final[df_final$SID == "1005" & df_final$RECORDING_SESSION_LABEL == "may05kg4", "SID"] <- "2001" #subject5a should be 2001 (artist).
df_final[df_final$SID == "1005" & df_final$RECORDING_SESSION_LABEL == "may05kg5", "SID"] <- "2001" #subject5a should be 2001 (artist).
df_final[df_final$SID == "1005" & df_final$RECORDING_SESSION_LABEL == "may05kg6", "SID"] <- "2001" #subject5a should be 2001 (artist).


#Getting rid of Subjects we're excluding.
df_final <- subset(df_final, subset = df_final$SID != "1003")
df_final <- subset(df_final, subset = df_final$SID != "1006")
df_final <- subset(df_final, subset = df_final$SID != "1011")
df_final <- subset(df_final, subset = df_final$SID != "1018")
df_final <- subset(df_final, subset = df_final$SID != "2001")

unique(df_final$SID) #tests to see if this error prone process really worked (depends on if file name is properly done)
df_temmpcheck2 <- df_final %>% group_by(SID, block) %>% summarise(length(task))



############################################################################################################################################################################################################
#Checking accuracy###

#convert response registered
df_final$ConvertedTarget <- ifelse(grepl("0", df_final$FCRTargetIndex), "2",
                                       ifelse(grepl("1", df_final$FCRTargetIndex), "1",
                                              ifelse(grepl("2", df_final$FCRTargetIndex), "3","error")))
#Calculate Accuracy
df_final$Accuracy <- ifelse((df_final$ConvertedTarget==df_final$FCRResponse), "correct", 
                                ifelse((df_final$FCRResponse==-1), "no_response", 
                                       ifelse((df_final$FCRResponse==4), "NoGuess", "incorrect")))


# #making subset with just test rows.#This is also where trial 65s with "." as values get deleted 
# df_test_full <- df_test_full[df_test_full$task==2,] 
# 
# #Getting percentage correct for each subject. Still not ideal because done partially manually so can't be faceted by block later :(
# grouped_df_bySID <- group_by(df_test_full, SID)
# temp.table.ga <- (table(grouped_df_bySID$SID, grouped_df_bySID$Accuracy))
# temp.table.ga <- as.data.frame(temp.table.ga)
# temp.wide.ga <- spread(temp.table.ga, Var2, Freq) 
# temp.wide.ga$freqsum <- temp.wide.ga$correct + temp.wide.ga$incorrect +temp.wide.ga$NoGuess +temp.wide.ga$no_response
# temp.wide.ga$Percent.Correct <- (temp.wide.ga$correct / temp.wide.ga$freqsum) * 100
# temp.wide.ga$Percent.incorrect <- (temp.wide.ga$incorrect / temp.wide.ga$freqsum) * 100
# temp.wide.ga$Percent.guess <- (temp.wide.ga$NoGuess / temp.wide.ga$freqsum) * 100
# temp.wide.ga$Percent.no_response <- (temp.wide.ga$no_response / temp.wide.ga$freqsum) * 100
# temp.wide.ga$Perc.cor.ifguessingall <- (temp.wide.ga$Percent.Correct + (temp.wide.ga$Percent.guess/3 + temp.wide.ga$Percent.no_response/3)) 
# wide.df.percentAcc <- gather(temp.wide.ga, Percent.Correct, Percent.incorrect, Percent.guess, Percent.no_response, key = "Accuracy.type", value = "Percent")

#temp.wide.ga
#install.packages("xlsx")
#library(xlsx) #load excel saving package
#write.xlsx(x = temp.wide.ga, "C:/Users/kgeier/Google Drive/OlsenLab_cloud/All Experiments at the Olsen Lab/MAYO/R studio analysis/Accuracyv4.xlsx",
#           sheetName = "summary", row.names = FALSE)

# #graphs faceted by SID with total count.
# pAccuracy <- ggplot(df_test_full, aes(factor(df_test_full$Accuracy))) +
#   geom_bar() +
#   facet_grid(.~ SID) +
#   geom_hline(aes(yintercept=137.48, linetype="dashed", size=.02)) +
#   theme_classic() +
#   theme(axis.line.x= element_line(colour = "black"), axis.line.y = element_line(color= "black")) +
#   theme(plot.title = element_text(size=36), axis.text = element_text(size=24), 
#         axis.title = element_text(size=32),
#         strip.text.x = element_text(size=16)) +
#   theme(legend.position="none") +
#   ggtitle("placeholder") +
#   xlab("Type of Accuracy") +
#   ylab("Count Correct")
# #pAccuracy




#1) Cleaning data#########################################################################################################################################################################################################

###Code to label all trials without eye-tracking as 'bad'###

#Make a new column with SID, block, and trial so there is a unique identifier for all Interest Areas to group by. 
df_final$Subj_trial <- paste0(df_final$SID, "_", df_final$block, "_", df_final$TRIAL_LABEL)

#This code makes a new df that has piped into it all the trials (all interest areas) that do NOT have 0 fixations on full screne. This means can have 0 fixations for part of the image, but if 0 for the whole screen means eyes were closed or worse.
#Note: Current set-up tags 0 fixations. May also want to tag with more specificity (eg. fewer than 5, more than 30??). (code credit to akrun on stackoverflow) 

#df_final2 made to just have trials we want to tag.   
df_final2 <-  df_final %>%
  group_by(Subj_trial) %>%
  filter(any(IA_FIXATION_COUNT==0 & IA_LABEL == "FullSceneInterestArea")) 
#df_final3 made to just have trials we want to keep as eye-tracking working.
df_final3 <- df_final %>%
  group_by(Subj_trial) %>%
  filter(!any(IA_FIXATION_COUNT==0 & IA_LABEL == "FullSceneInterestArea")) 

#giving them different values for new column. 1 for bad(ignore these trials in GLMs) 0 for okay. 
df_final2$noET <- 1 # 1 means bad - no eye-t
df_final3$noET <- 0

df_finalbound.withETtags <- rbind(df_final2, df_final3)
df_final <- df_finalbound.withETtags[order(df_finalbound.withETtags$orderer),]

#c) Divide into DPM DPNM#####################################
#Divide into 2 categories of trial. 1)-DPMatch at least 10% more viewing of matched (matched means target face) face than next highest face. Compare DPM to those with 10% higher viewing to WRONG face. Equal viewing trials are ignored (so partly comparing to "inappropriate memory" not to "no memory")

#if statement to make new column with whether or not that row was CORRECT interest area
df_final$CorrectIA <- ifelse((df_final$IA_LABEL=="FullSceneInterestArea"), "Full",
                                   ifelse((df_final$ConvertedTarget == 1 & df_final$IA_LABEL=="LeftFaceInterestArea"), "corIA",
                                          ifelse((df_final$ConvertedTarget == 2 & df_final$IA_LABEL=="TopMiddleFaceInterestArea"), "corIA",
                                                 ifelse((df_final$ConvertedTarget == 3 & df_final$IA_LABEL=="RightFaceInterestArea"), "corIA", "incorIA"))))


#if statement to make new column with whether or not that row was GUESSED interets area
df_final$ChosenIA <- ifelse((df_final$IA_LABEL=="FullSceneInterestArea"), "Full",
                                  ifelse((df_final$FCRResponse == 4), "noguess",
                                         ifelse((df_final$FCRResponse == -1), "noresponse",
                                                ifelse((df_final$FCRResponse == 1 & df_final$IA_LABEL=="LeftFaceInterestArea"), "chosenIA",
                                                       ifelse((df_final$FCRResponse == 2 & df_final$IA_LABEL=="TopMiddleFaceInterestArea"), "chosenIA",
                                                              ifelse((df_final$FCRResponse == 3 & df_final$IA_LABEL=="RightFaceInterestArea"), "chosenIA", "NotChosenIA"))))))

#converting dwell time to numeric
df_final$`IA_DWELL_TIME_%` <- as.numeric(as.character(df_final$`IA_DWELL_TIME_%`)) 
class(df_final$`IA_DWELL_TIME_%`)


##############'''''''''''''''
#copy for subqacc #
df_subqacc <- df_final


#making a new subset without 'full' screne IA so can compare the remaining IA's easier. 
subset_df_noFull <-  df_final[!grepl("Full", df_final$CorrectIA),] 
#making a new subset with just test phase as that's only phase that had the relevant eye-movements.  
subset_df_noFull_t2 <-  subset_df_noFull[grepl(2, subset_df_noFull$task),]




#Nest step is to label each trial as disproportiante viewing (dp) vs not dp. After that can combine with match vs non-match
#To do this must order within each trial most to least viewing. Then if first most is 10% greater than second most label it as DP

#Thank you stack overflow for this easy bit of code. Makes a new df with ordered by Subj_block_trial THEN ordered by dwell time %
df.final2 <-subset_df_noFull_t2[order(subset_df_noFull_t2$Subj_trial,subset_df_noFull_t2$`IA_DWELL_TIME_%`),]

#multiplying by 100 for %
df.final2$IA_DwellTime_Percent <- df.final2$`IA_DWELL_TIME_%` *100

df.final2$rankedposition <- rep(1:3, (length(df.final2$IA_DwellTime_Percent)/3))

#overly complicated for loop in order to compare IAs with most viewing to trials with 2nd most viewing 
# and abel them as disproportionate or not depending on if view by more tahn 10 percentage points. 
df.final2$dp_ndp <- 'placeholder'
for (i in 2:length(df.final2$dp_ndp)){
  if(df.final2[i,39] == 1) {
    1 -> df.final2[i, 40]
}
  else if(df.final2[i,39]==2) {
    2 -> df.final2[i, 40]
}
    else if(is.na(df.final2[i,38])) {
    "TrialNotViewed" -> df.final2[i, 40]
    }
  else if(df.final2[i,35]==1) {
    "TrialNotViewed" -> df.final2[i, 40]
  }
  else if(df.final2[i,39]== 3 & df.final2[i, 38] > (10+(df.final2[i-1, 38]))) {
    'dp' -> df.final2[i, 40]
  }
    else if(df.final2[i,39]== 3 & df.final2[i, 38] <= (10+(df.final2[i-1, 38]))) {
      df.final2[i, 40] <- "ndp" 
  } else
      df.final2[i, 40] <- "ERROR" 
    }
1 -> df.final2[1, 40]   

#Please double check this, esp. if you have edited the code as it is set up to deal with column number 
#not column name which is very sensitive to change. 


df.final2$dpm_dpnm <- ifelse((df.final2$dp_ndp<3), 'NotRelevant', 
                          ifelse((df.final2$dp_ndp == 'TrialNotViewed'), 'TrialNotViewed',
                            ifelse((df.final2$dp_ndp == 'ndp'), 'NDPV', #ndpv is Not DisProportional Viewing trial
                              ifelse((df.final2$dp_ndp == 'dp' & df.final2$CorrectIA =="corIA"), "DPM", #DPM is DisProportional viewing to Match (i.e. correct face)
                                ifelse((df.final2$dp_ndp == 'dp' & df.final2$CorrectIA =="incorIA"), "DPNM", 'Error'))))) #DPNM is DisProportional viewing to Non-Match (i.e. incorrect face). Note this isn't absence of eye-mvmt memory but 'inappropriate' memory to wrong face. 
                                   
                                  
#sanity check #Count of DPM and DPNM and NDPV
a<-nrow(df.final2[grepl("DPNM", df.final2$dpm_dpnm),])
b<-nrow(df.final2[grepl("DPM", df.final2$dpm_dpnm),])
c<-nrow(df.final2[grepl("NDPV", df.final2$dpm_dpnm),])
d<-nrow(df.final2[grepl("TrialNotViewed", df.final2$dpm_dpnm),])
sum(a,b,c,d)

#tmrw Kirk, please make a table with number of each trial type for each subject. 

df.final.dpBins <- df.final2[!grepl("NotRelevant", df.final2$dpm_dpnm),]
df.final.dpBins <- df.final.dpBins[,c('Subj_trial','dpm_dpnm')]

write.csv(x = df.final.dpBins, "C:/Users/kgeier/ViewingType28subj.csv", row.names = FALSE)

#ARG. BC google drive not syncing, won't save automatically to drive. So saving to files then moving manually to drive. 
#Add Google Drive/OlsenLab_cloud/All Experiments at the Olsen Lab/MAYO/R studio analysis/ to path to get it to gdrive. 
)


#Data Cleaning#### Kind of confusing bc did study and test phase separately (dif. box plots) but now fuse with dif. dataframes
#test#Weird names bc copied exactly from the study cleaning
df_subq <- subset(df_subqacc, subset = df_subqacc$CorrectIA == "Full")

df_stud <- subset(df_subq, subset = df_subq$task == 2)

#cleaning just df_stud based on study fixations (tukey boxplot method)
#Find and remove all trials with 0 fixations. 
df.no0fix <- df_stud[df_stud$IA_FIXATION_COUNT !=0 ,] 

#checking trial counts for each subject with 0 fixation trials removed
table(df.no0fix$SID)

#Calculating quartiles and median #REMOVE "IP_LABEL" if want to group by suject only not subject and study/test
df.sum <- ddply(df.no0fix, .(SID, IP_LABEL), function(df.no0fix) quantile(df.no0fix$IA_FIXATION_COUNT, na.rm =T))

#calculating pos(ible) and prob(able) ranges for outliers
df.sum$BotWhisk_pos <- (df.sum$`25%`-(1.5* (df.sum$`75%`-df.sum$`25%`))) 
df.sum$TopWhisk_pos <- (df.sum$`75%`+ (1.5* (df.sum$`75%`-df.sum$`25%`)))
df.sum$BotWhisk_prob <- (df.sum$`25%`-(3 * (df.sum$`75%`-df.sum$`25%`)))
df.sum$TopWhisk_prob <- (df.sum$`75%`+(3 * (df.sum$`75%`-df.sum$`25%`)))

#Making new columns and merging boxplot extremes with core data
df.no0fix$Outlier_pos <- "placeholder"
df.no0fix$Outlier_prob <- "ph"
df.no0fix <- merge(df.no0fix, df.sum, by = c("SID")) 

#marking all trials that have fixations outside of the whiskers as yes or no. 
df.no0fix$Outlier_pos <- ifelse((df.no0fix$IA_FIXATION_COUNT < df.no0fix$BotWhisk_pos | df.no0fix$IA_FIXATION_COUNT > df.no0fix$TopWhisk_pos), "yes", "no")
df.no0fix$Outlier_prob <- ifelse((df.no0fix$IA_FIXATION_COUNT < df.no0fix$BotWhisk_prob | df.no0fix$IA_FIXATION_COUNT > df.no0fix$TopWhisk_prob), "yes", "no")

#counting number of trials as eliminated or kept. 
table(df.no0fix$Outlier_pos)
table(df.no0fix$Outlier_prob)

df.no0fix4m <- df.no0fix[, c("Subj_trial", "Outlier_prob", "Outlier_pos")]
df.mergeOutliers <- merge(df.final2, df.no0fix4m, by = "Subj_trial", all.y = F)
#in this case removing all probable outliers to not bias the data. 
#Feel free to adjust. 
df.clean <- df.mergeOutliers[df.mergeOutliers$Outlier_prob == "no",]  

df.final2 <-df.clean
#############################################################################################################

#study (for subq acc)

#subsetting then merging study and test phase to match accuracy
df_subq <- subset(df_subqacc, subset = df_subqacc$CorrectIA == "Full")

df_stud <- subset(df_subq, subset = df_subq$task == 1)

#cleaning just df_stud based on study fixations (tukey boxplot method)
#Find and remove all trials with 0 fixations. 
df.no0fix <- df_stud[df_stud$IA_FIXATION_COUNT !=0 ,] 

#checking trial counts for each subject with 0 fixation trials removed
table(df.no0fix$SID)

#Calculating quartiles and median #REMOVE "IP_LABEL" if want to group by suject only not subject and study/test
df.sum <- ddply(df.no0fix, .(SID, IP_LABEL), function(df.no0fix) quantile(df.no0fix$IA_FIXATION_COUNT, na.rm =T))

#calculating pos(ible) and prob(able) ranges for outliers
df.sum$BotWhisk_pos <- (df.sum$`25%`-(1.5* (df.sum$`75%`-df.sum$`25%`))) 
df.sum$TopWhisk_pos <- (df.sum$`75%`+ (1.5* (df.sum$`75%`-df.sum$`25%`)))
df.sum$BotWhisk_prob <- (df.sum$`25%`-(3 * (df.sum$`75%`-df.sum$`25%`)))
df.sum$TopWhisk_prob <- (df.sum$`75%`+(3 * (df.sum$`75%`-df.sum$`25%`)))

#Making new columns and merging boxplot extremes with core data
df.no0fix$Outlier_pos <- "placeholder"
df.no0fix$Outlier_prob <- "ph"
df.no0fix <- merge(df.no0fix, df.sum, by = c("SID")) 

#marking all trials that have fixations outside of the whiskers as yes or no. 
df.no0fix$Outlier_pos <- ifelse((df.no0fix$IA_FIXATION_COUNT < df.no0fix$BotWhisk_pos | df.no0fix$IA_FIXATION_COUNT > df.no0fix$TopWhisk_pos), "yes", "no")
df.no0fix$Outlier_prob <- ifelse((df.no0fix$IA_FIXATION_COUNT < df.no0fix$BotWhisk_prob | df.no0fix$IA_FIXATION_COUNT > df.no0fix$TopWhisk_prob), "yes", "no")

#counting number of trials as eliminated or kept. 
table(df.no0fix$Outlier_pos)
table(df.no0fix$Outlier_prob)

#in this case removing all probable outliers to not bias the data. 
#Feel free to adjust. 
df.clean <- df.no0fix[df.no0fix$Outlier_prob == "no",]  

df_stud <- df.clean

df_test <- subset(df_subq, subset = df_subq$task == 2)

df_subqaccM <- merge(df_stud, df_test, by =c('SID', "block", "facesimg"), all.x = F, all.y = F)

#plotting fixation count by accuracy.y

# df_subqaccM %>% group_by(Accuracy.y) %>% summarise(IA_FIX = mean(IA_FIXATION_COUNT.x))
# 
# plot_fixbySubA <- df_subqaccM %>% group_by(Accuracy.y) %>% summarise(IA_FIX = mean(IA_FIXATION_COUNT.x))  %>% ggplot(aes(x=Accuracy.y, y=IA_FIX)) +
#   #geom_errorbar(width=0, aes(ymin=FIXATION_COUNT-ci, ymax=FIXATION_COUNT+ci), colour="red",position=pd, size=1.5) +
#  # geom_errorbar(width=0, aes(ymin=FIXATION_COUNT-ci, ymax=FIXATION_COUNT+ci), data=d_sum, position=pd, size=1.5) +
#   geom_point(size=3) 
#  # geom_smooth(method='lmer',formula=y~x + (1|Subject.Number) + (1|studyphase), data = summary.df_eachviewing)
#   #geom_text(x = 0.4, y = 8, label = lm_eqn(dfp), parse = TRUE)
#   plot_fixbySubA


  pd <- position_dodge(0.1) # move them .05 to the left and right
  
#finding SE
  
  d_sum = ddply(df_subqaccM, ~Accuracy.y + SID,
                summarise,AVERAGE_FIX_COUNT=mean(IA_FIXATION_COUNT.x)) 
  
  detach("package:dplyr", unload=TRUE)
  
  d_sum <- d_sum[!is.na(d_sum$AVERAGE_FIX_COUNT),]
  
  d_sum = summarySEwithin(d_sum, 
                          measurevar="AVERAGE_FIX_COUNT", 
                          withinvars="Accuracy.y",
                          idvar="SID") 
  d_sum
  
  plot2_fixbySubA <- ggplot(d_sum, aes(x=Accuracy.y, y=AVERAGE_FIX_COUNT)) +
   # geom_line(position=pd, size=1.5) +
    geom_errorbar(width=0, aes(ymin=AVERAGE_FIX_COUNT-ci, ymax=AVERAGE_FIX_COUNT+ci), colour="red",position=pd, size=1.5) +
    geom_errorbar(width=0, aes(ymin=AVERAGE_FIX_COUNT-ci, ymax=AVERAGE_FIX_COUNT+ci), data=d_sum, position=pd, size=1.5) +
    geom_point(size=5, position=pd) + 
    scale_color_brewer(palette='Set1')
  plot2_fixbySubA
  
  d_sum2 = ddply(df_subqaccM, ~ Accuracy.y  + SID,
                 summarise,AVERAGE_FIX_COUNT=mean(IA_FIXATION_COUNT.x))
  #trying Anova instead of linear model bc. categorical on X axis. 
  res_SA <- aov(AVERAGE_FIX_COUNT ~ Accuracy.y, d_sum2)
  summary(res_SA)
  #nothing significant so no need for post hoc testing. 
  
  # #THIS IS RANDOME INTERCEPT FIXED EFFECTS. SO ASSUMING EACH SUBJECT MAY have different total number of fixations but a similar change due to subq acc.
  # #May need to try other models and compare for best fit?
  # res = lmer(AVERAGE_FIX_COUNT ~ Accuracy.y  + (1|SID),
  #             data=d_sum2)
  # summary(res) #no res1


#Graphs of dwell time to correct vs incorrect faces####

  #df.final2$ChosCor4Opt <- paste0(df.final2$CorrectIA, "_", df.final2$ChosenIA)
  
  #finding SE
  
  d_sum3 = ddply(df.final2, ~CorrectIA*ChosenIA + SID,
                summarise,Average_Dwell_Time=mean(IA_DwellTime_Percent)) 
  
  d_sum3 <- d_sum3[!is.na(d_sum3$Average_Dwell_Time),]
  
  detach("package:dplyr", unload=TRUE)
  d_sum3 = summarySEwithin(d_sum3, 
                          measurevar="Average_Dwell_Time", 
                          withinvars=c("CorrectIA","ChosenIA"),
                          idvar="SID") 
  d_sum3
  
  #Remove this next line if itnerested in "no response and no guess"
  # d_sum3 <- subset(d_sum3, subset = ChosenIA %in%  c("chosenIA","NotChosenIA"))
  # d_sum3
  
plot3_dwelltime <- ggplot(d_sum3, aes(x=CorrectIA, y=Average_Dwell_Time, group = ChosenIA, fill = ChosenIA, color = ChosenIA)) +
    geom_line(position=pd, size=1.5) +
    geom_errorbar(width=0, aes(ymin=Average_Dwell_Time-ci, ymax=Average_Dwell_Time+ci), colour="red",position=pd, size=1.5) +
    geom_errorbar(width=0, aes(ymin=Average_Dwell_Time-ci, ymax=Average_Dwell_Time+ci), data=d_sum3, position=pd, size=1.5) +
  #  geom_point() + 
    scale_color_brewer(palette='Set1')
plot3_dwelltime
  
  d_sum4 = ddply(df.final2, ~ CorrectIA * ChosenIA  + SID,
                 summarise,Average_Dwell_Time=mean(IA_DwellTime_Percent))
  
  #Remove following lines if want NoGuess and noresponse. Then use the post-hoc testing below
  # d_sum4 <- subset(d_sum4, subset = ChosenIA %in%  c("chosenIA","NotChosenIA"))
  # d_sum4
 #again trying anova (2 way) bc x is categorical (but this time 2 way bc it's correct vs incor AND Chosen vs not chosen)
res_DT <- aov(Average_Dwell_Time~CorrectIA*ChosenIA,d_sum4)  
summary(res_DT)  

with(d_sum4, pairwise.t.test(Average_Dwell_Time, ChosenIA, 
                                    p.adj="holm", paired=F))


#THIS IS RANDOM INTERCEPT FIXED EFFECTS. SO ASSUMING EACH SUBJECT MAY have different total number of fixations but a similar change due to subq acc. 
  #May need to try other models and compare for best fit?
  res2 = lmer(Average_Dwell_Time ~ CorrectIA * ChosenIA + (1|SID),
              data=d_sum4)
  summary(res2) #no res1 

#Find out a very large difference due to chosen vs. not chosen. Subtle difference between correct vs incorrect. 
#Maybe should filter out no response and no guess   


  
  
  
  
#OLD CODE to be scrapped for pieces####
# 
# 
# #subset only "full" ia because then only 1 row per trial.
# subset_df_onlyFull <-  df_testcompact[grepl("FullSceneInterestArea", df_testcompact$IA_LABEL),]
# 
# )}
# 
# subset_df_onlyFullclean1 <-  subset_df_onlyFull[!grepl("no_response", subset_df_onlyFull$Accuracy),]
# subset_df_onlyFullclean2 <-  subset_df_onlyFullclean1[!grepl("guess", subset_df_onlyFullclean1$Accuracy),]
# subset_df_onlyFullClean <-  subset_df_onlyFullclean2
# 
# 
# #Goal: group by subject and count Correct vs. incorrect. Uses Full, not cleaned so will show how many 
# #no responses there are.
# df_count<- subset_df_onlyFull %>% 
#   group_by(Accuracy) %>%
#   summarise(countAcc = length(Accuracy))
# 
# summary_by_subj <-subset_df_onlyFull %>% group_by(cb, Accuracy)
# 
# #playing around trying to make accuracy a factor, count numeric, and the whole thing a dataframe not matrix. 
# #Not sure how much of this is neccessary
# df_count$countAcc <- as.numeric(as.character(df_count$countAcc))
# df_count$Accuracy <- as.factor(df_count$Accuracy)
# df_count <- as.data.frame(df_count)
# df_count
# 
# #graphs of accuracy 
# #horizontal line is just 34.37. I.e. number would get correct if never guessed 4 but randomly went 1,2,3. 
# p9 <- ggplot(subset_df_onlyFull, aes(factor(subset_df_onlyFull$Accuracy))) +
#   geom_bar() +
#   facet_wrap(~cb) +
#   geom_hline(aes(yintercept=34.37, linetype="dashed", size=.02)) +
#   theme_classic() +
#   theme(axis.line.x= element_line(colour = "black"), axis.line.y = element_line(color= "black")) +
#   theme(plot.title = element_text(size=36), axis.text = element_text(size=24), 
#         axis.title = element_text(size=32),
#         strip.text.x = element_text(size=16)) +
#   theme(legend.position="none") +
#   ggtitle("placeholder") +
#   xlab("Type of Accuracy") +
#   ylab("Count Correct")
# p9
# 
# 
# #################################P3- More fixations for correct Interest Area? 
# 
# #if statement to make new column with whether or not that row was CORRECT interest area
# df_testcompact$CorrectIA <- ifelse((df_testcompact$IA_LABEL=="FullSceneInterestArea"), "Full",
#                                    ifelse((df_testcompact$ConvertedTarget == 1 & df_testcompact$IA_LABEL=="LeftFaceInterestArea"), "corIA",
#                                           ifelse((df_testcompact$ConvertedTarget == 2 & df_testcompact$IA_LABEL=="TopMiddleFaceInterestArea"), "corIA",
#                                                  ifelse((df_testcompact$ConvertedTarget == 3 & df_testcompact$IA_LABEL=="RightFaceInterestArea"), "corIA", "incorIA"))))
# 
# 
# #if statement to make new column with whether or not that row was GUESSED interets area
# df_testcompact$ChosenIA <- ifelse((df_testcompact$IA_LABEL=="FullSceneInterestArea"), "Full",
#                                   ifelse((df_testcompact$FCRResponse == 4), "noguess",
#                                          ifelse((df_testcompact$FCRResponse == -1), "noresponse",
#                                                 ifelse((df_testcompact$FCRResponse == 1 & df_testcompact$IA_LABEL=="LeftFaceInterestArea"), "chosenIA",
#                                                        ifelse((df_testcompact$FCRResponse == 2 & df_testcompact$IA_LABEL=="TopMiddleFaceInterestArea"), "chosenIA",
#                                                               ifelse((df_testcompact$FCRResponse == 3 & df_testcompact$IA_LABEL=="RightFaceInterestArea"), "chosenIA", "NotChosenIA"))))))
# 
# #converting dwell time to numeric
# df_testcompact$`IA_DWELL_TIME_%` <- as.numeric(as.character(df_testcompact$`IA_DWELL_TIME_%`)) 
# class(df_testcompact$`IA_DWELL_TIME_%`)
# 
# ##############'''''''''''''''
# 
# #making a new subset without 'full' screne IA so can compare the remaining IA's easier. 
# subset_df_noFull <-  df_workingETonly[!grepl("Full", df_workingETonly$CorrectIA),]
# 
# #boxplot dwell time subset without full
# p01 <- ggplot(subset_df_noFull, aes(factor(subset_df_noFull$CorrectIA), subset_df_noFull$`IA_DWELL_TIME_%`)) +
#   geom_boxplot() +
#   theme_classic() +
#   theme(axis.line.x= element_line(colour = "black"), axis.line.y = element_line(color= "black")) +
#   theme(plot.title = element_text(size=36), axis.text = element_text(size=24), 
#         axis.title = element_text(size=32),
#         strip.text.x = element_text(size=16)) +
#   theme(legend.position="none") +
#   ggtitle("Implicit Memory for Correct Face") +
#   xlab("Correct or Not Correct IA") +
#   ylab("Dwell Time %")
# p01
# 
# p02 <- ggplot(subset_df_noFull, aes(factor(subset_df_noFull$ChosenIA), subset_df_noFull$`IA_DWELL_TIME_%`)) +
#   geom_boxplot() +
#   theme_classic() +
#   theme(axis.line.x= element_line(colour = "black"), axis.line.y = element_line(color= "black")) +
#   theme(plot.title = element_text(size=36), axis.text = element_text(size=24), 
#         axis.title = element_text(size=32),
#         strip.text.x = element_text(size=16)) +
#   theme(legend.position="none") +
#   ggtitle("Implicit Memory for Correct Face") +
#   xlab("Correct or Not Correct IA") +
#   ylab("Dwell Time %")
# p02
# 
# grid.arrange(p01,p02)
# 
# 














# #Violin plot for dwell time
# p6 <- ggplot(df_testcompact, aes(factor(df_testcompact$CorrectIA), df_testcompact$`IA_DWELL_TIME_%`)) +
#   geom_violin() +
#   geom_point(aes(), size= 2, stroke=1.3, position =position_jitter(width=.3, height=0)) +
#   theme_classic() +
#   theme(axis.line.x= element_line(colour = "black"), axis.line.y = element_line(color= "black")) +
#   theme(plot.title = element_text(size=36), axis.text = element_text(size=24), 
#         axis.title = element_text(size=32),
#         strip.text.x = element_text(size=16)) +
#   theme(legend.position="none") +
#   ggtitle("Implicit Memory for Correct Face") +
#   xlab("Correct or Not Correct IA") +
#   ylab("Dwell Time %")
# p6
# 
# p7 <- ggplot(df_testcompact, aes(factor(df_testcompact$ChosenIA), df_testcompact$`IA_DWELL_TIME_%`)) +
#   geom_violin() +
#   geom_point(aes(), size= 2, stroke=1.3, position =position_jitter(width=.3, height=0)) +
#   theme_classic() +
#   theme(axis.line.x= element_line(colour = "black"), axis.line.y = element_line(color= "black")) +
#   theme(plot.title = element_text(size=36), axis.text = element_text(size=24), 
#         axis.title = element_text(size=32),
#         strip.text.x = element_text(size=16)) +
#   theme(legend.position="none") +
#   ggtitle("Implicit Memory for Correct Face") +
#   xlab("Correct or Not Correct IA") +
#   ylab("Dwell Time %")
# p7
# 
# grid.arrange(p6,p7)
  
  
  
  
  
  
  
  
  
  
  #By time period####
  
  #import csv
  df_tp1 <- fread("C:/Users/kgeier/Google Drive/OlsenLab_cloud/All Experiments/MAYO/R studio analysis/TP1_Dwell.csv")
  df_tp2 <- fread("C:/Users/kgeier/Google Drive/OlsenLab_cloud/All Experiments/MAYO/R studio analysis/TP2_Dwell.csv")
  df_tp3 <- fread("C:/Users/kgeier/Google Drive/OlsenLab_cloud/All Experiments/MAYO/R studio analysis/TP3_Dwell.csv")
  df_tp4 <- fread("C:/Users/kgeier/Google Drive/OlsenLab_cloud/All Experiments/MAYO/R studio analysis/TP4_Dwell.csv")
  df_tp5 <- fread("C:/Users/kgeier/Google Drive/OlsenLab_cloud/All Experiments/MAYO/R studio analysis/TP5_Dwell.csv")
  df_tp6 <- fread("C:/Users/kgeier/Google Drive/OlsenLab_cloud/All Experiments/MAYO/R studio analysis/TP6_Dwell.csv")
  
df_alltp <- rbind(df_tp1,df_tp2,df_tp3,df_tp4,df_tp5,df_tp6)
  
#adding orderer column
df_test_withcount <-cbind(df_alltp, seq(1:nrow(df_alltp)))
names(df_test_withcount)[names(df_test_withcount) == 'V2'] <- 'orderer'


#re-ordering
df_final <- df_test_withcount

#Cleaning trials and blocks manually that have been noted as bad. ####
#getting rid of blocks that didn't trigger, etc.
df_final <- df_final[df_final$RECORDING_SESSION_LABEL != "may14kg6"]
df_final <- df_final[df_final$RECORDING_SESSION_LABEL != "may05kg4"]

#changing TRIAL_LABEL from Trial: 19 to 19. 
df_final$TRIAL_LABEL <- gsub("Trial: ", "", df_final$TRIAL_LABEL) 

###- REMOVING 7 trials (both study and test) that had innapropriately had a mix of older/yonger/male/female in the lures. Ex. OF11, YF15, OF23. 
#trials were determined to be 
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may16kg1" | df_final$`TRIAL_LABEL` !="15",)
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may16kg4" | df_final$`TRIAL_LABEL` !="42",)
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may18kg3" | df_final$`TRIAL_LABEL` !="25",)
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may04kg1" | df_final$`TRIAL_LABEL` !="19",)
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may06kg3" | df_final$`TRIAL_LABEL` !="35",)
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may11kg4" | df_final$`TRIAL_LABEL` !="24",)
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may23kg4" | df_final$`TRIAL_LABEL` !="36",)
#removing an 8th and 9th trial b.c. lure face matched target face
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may22kg6" | df_final$`TRIAL_LABEL` !="9",)
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may13kg6" | df_final$`TRIAL_LABEL` !="19",)

df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may16kg1" | df_final$`TRIAL_LABEL` !="60",)
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may16kg4" | df_final$`TRIAL_LABEL` !="54",)
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may18kg3" | df_final$`TRIAL_LABEL` !="55",)
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may04kg1" | df_final$`TRIAL_LABEL` !="64",)
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may06kg3" | df_final$`TRIAL_LABEL` !="60",)
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may11kg4" | df_final$`TRIAL_LABEL` !="50",)
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may23kg4" | df_final$`TRIAL_LABEL` !="53",)
#removing an 8th and 9thtrial b.c. lure face matched target face
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may22kg6" | df_final$`TRIAL_LABEL` !="62",)
df_final <- subset(df_final, df_final$RECORDING_SESSION_LABEL != "may13kg6" | df_final$`TRIAL_LABEL` !="53",)

#adding Subject ID to both study and test. 
df_final$SLabNoBlock <- df_final$RECORDING_SESSION_LABEL #makes a copy of session label to get SID from
df_final$SLabNoBlock <- substring(df_final$SLabNoBlock,1, 7) #selects first 7 charachters to remove block number
df_final$SID <- gsub("[^0-9]","",df_final$SLabNoBlock) #removes all non-numbers
df_final$SID <- gsub("146","14",df_final$SID) #subject 14 has issue with not being same #char. So this fixes it.
df_final$SID <- paste0("10",df_final$SID)
#still not working... use lines below   df_final$SID <- ifelse(grepl("may05bk", df_final$RECORDING_SESSION_LABEL), gsub("1005 ", "2001", df_final$SID), df_final$SID <- df_final$SID)
#replace next 6 lines with regular expressions when you have more time. 
df_final[df_final$SID == "1005" & df_final$RECORDING_SESSION_LABEL == "may05kg1", "SID"] <- "2001" #subject5a should be 2001 (artist).
df_final[df_final$SID == "1005" & df_final$RECORDING_SESSION_LABEL == "may05kg2", "SID"] <- "2001" #subject5a should be 2001 (artist).
df_final[df_final$SID == "1005" & df_final$RECORDING_SESSION_LABEL == "may05kg3", "SID"] <- "2001" #subject5a should be 2001 (artist).
df_final[df_final$SID == "1005" & df_final$RECORDING_SESSION_LABEL == "may05kg4", "SID"] <- "2001" #subject5a should be 2001 (artist).
df_final[df_final$SID == "1005" & df_final$RECORDING_SESSION_LABEL == "may05kg5", "SID"] <- "2001" #subject5a should be 2001 (artist).
df_final[df_final$SID == "1005" & df_final$RECORDING_SESSION_LABEL == "may05kg6", "SID"] <- "2001" #subject5a should be 2001 (artist).


#Getting rid of Subjects we're excluding.
df_final <- subset(df_final, subset = df_final$SID != "1003")
df_final <- subset(df_final, subset = df_final$SID != "1006")
df_final <- subset(df_final, subset = df_final$SID != "1011")
df_final <- subset(df_final, subset = df_final$SID != "1018")
df_final <- subset(df_final, subset = df_final$SID != "2001")

unique(df_final$SID) #tests to see if this error prone process really worked (depends on if file name is properly done)
df_temmpcheck2 <- df_final %>% group_by(SID, block) %>% summarise(length(task))



############################################################################################################################################################################################################
#Checking accuracy###

#convert response registered
df_final$ConvertedTarget <- ifelse(grepl("0", df_final$FCRTargetIndex), "2",
                                   ifelse(grepl("1", df_final$FCRTargetIndex), "1",
                                          ifelse(grepl("2", df_final$FCRTargetIndex), "3","error")))
#Calculate Accuracy
df_final$Accuracy <- ifelse((df_final$ConvertedTarget==df_final$FCRResponse), "correct", 
                            ifelse((df_final$FCRResponse==-1), "no_response", 
                                   ifelse((df_final$FCRResponse==4), "NoGuess", "incorrect")))



#1) Cleaning data#########################################################################################################################################################################################################

###Code to label all trials without eye-tracking as 'bad'###

#Make a new column with SID, block, and trial so there is a unique identifier for all Interest Areas to group by. 
df_final$Subj_trial <- paste0(df_final$SID, "_", df_final$block, "_", df_final$TRIAL_LABEL)

#This code makes a new df that has piped into it all the trials (all interest areas) that do NOT have 0 fixations on full screne. This means can have 0 fixations for part of the image, but if 0 for the whole screen means eyes were closed or worse.
#Note: Current set-up tags 0 fixations. May also want to tag with more specificity (eg. fewer than 5, more than 30??). (code credit to akrun on stackoverflow) 

#df_final2 made to just have trials we want to tag.   
df_final2 <-  df_final %>%
  group_by(Subj_trial) %>%
  filter(any(IA_FIXATION_COUNT==0 & IA_LABEL == "FullSceneInterestArea")) 
#df_final3 made to just have trials we want to keep as eye-tracking working.
df_final3 <- df_final %>%
  group_by(Subj_trial) %>%
  filter(!any(IA_FIXATION_COUNT==0 & IA_LABEL == "FullSceneInterestArea")) 

#giving them different values for new column. 1 for bad(ignore these trials in GLMs) 0 for okay. 
df_final2$noET <- 1 # 1 means bad - no eye-t
df_final3$noET <- 0

df_finalbound.withETtags <- rbind(df_final2, df_final3)
df_final <- df_finalbound.withETtags[order(df_finalbound.withETtags$orderer),]

#c) Divide into DPM DPNM#####################################
#Divide into 2 categories of trial. 1)-DPMatch at least 10% more viewing of matched (matched means target face) face than next highest face. Compare DPM to those with 10% higher viewing to WRONG face. Equal viewing trials are ignored (so partly comparing to "inappropriate memory" not to "no memory")

#if statement to make new column with whether or not that row was CORRECT interest area
df_final$CorrectIA <- ifelse((df_final$IA_LABEL=="FullSceneInterestArea"), "Full",
                             ifelse((df_final$ConvertedTarget == 1 & df_final$IA_LABEL=="LeftFaceInterestArea"), "corIA",
                                    ifelse((df_final$ConvertedTarget == 2 & df_final$IA_LABEL=="TopMiddleFaceInterestArea"), "corIA",
                                           ifelse((df_final$ConvertedTarget == 3 & df_final$IA_LABEL=="RightFaceInterestArea"), "corIA", "incorIA"))))


#if statement to make new column with whether or not that row was GUESSED interets area
df_final$ChosenIA <- ifelse((df_final$IA_LABEL=="FullSceneInterestArea"), "Full",
                            ifelse((df_final$FCRResponse == 4), "noguess",
                                   ifelse((df_final$FCRResponse == -1), "noresponse",
                                          ifelse((df_final$FCRResponse == 1 & df_final$IA_LABEL=="LeftFaceInterestArea"), "chosenIA",
                                                 ifelse((df_final$FCRResponse == 2 & df_final$IA_LABEL=="TopMiddleFaceInterestArea"), "chosenIA",
                                                        ifelse((df_final$FCRResponse == 3 & df_final$IA_LABEL=="RightFaceInterestArea"), "chosenIA", "NotChosenIA"))))))

#converting dwell time to numeric
df_final$`IA_DWELL_TIME_%` <- as.numeric(as.character(df_final$`IA_DWELL_TIME_%`)) 
class(df_final$`IA_DWELL_TIME_%`)

#copy for subqacc #
df_subqacc <- df_final

#making a new subset without 'full' screne IA so can compare the remaining IA's easier. 
subset_df_noFull <-  df_final[!grepl("Full", df_final$CorrectIA),] 
#making a new subset with just test phase as that's only phase that had the relevant eye-movements.  
subset_df_noFull_t2 <-  subset_df_noFull[grepl(2, subset_df_noFull$task),]


#Nest step is to label each trial as disproportiante viewing (dp) vs not dp. After that can combine with match vs non-match
#To do this must order within each trial most to least viewing. Then if first most is 10% greater than second most label it as DP

#Thank you stack overflow for this easy bit of code. Makes a new df with ordered by Subj_block_trial THEN ordered by dwell time %
df.final2 <-subset_df_noFull_t2[order(subset_df_noFull_t2$IP_LABEL,subset_df_noFull_t2$Subj_trial, subset_df_noFull_t2$`IA_DWELL_TIME_%`),]

#multiplying by 100 for %
df.final2$IA_DwellTime_Percent <- df.final2$`IA_DWELL_TIME_%` *100

df.final2$rankedposition <- rep(1:3, (length(df.final2$IA_DwellTime_Percent)/3))

#overly complicated for loop in order to compare IAs with most viewing to trials with 2nd most viewing 
# and abel them as disproportionate or not depending on if view by more tahn 10 percentage points. 
df.final2$dp_ndp <- 'placeholder'
for (i in 2:length(df.final2$dp_ndp)){
  if(df.final2[i,39] == 1) {
    1 -> df.final2[i, 40]
  }
  else if(df.final2[i,39]==2) {
    2 -> df.final2[i, 40]
  }
  else if(is.na(df.final2[i,38])) {
    "TrialNotViewed" -> df.final2[i, 40]
  }
  else if(df.final2[i,35]==1) {
    "TrialNotViewed" -> df.final2[i, 40]
  }
  else if(df.final2[i,39]== 3 & df.final2[i, 38] > (10+(df.final2[i-1, 38]))) {
    'dp' -> df.final2[i, 40]
  }
  else if(df.final2[i,39]== 3 & df.final2[i, 38] <= (10+(df.final2[i-1, 38]))) {
    df.final2[i, 40] <- "ndp" 
  } else
    df.final2[i, 40] <- "ERROR" 
}
1 -> df.final2[1, 40]   

#Please double check this, esp. if you have edited the code as it is set up to deal with column number 
#not column name which is very sensitive to change. 


df.final2$dpm_dpnm <- ifelse((df.final2$dp_ndp<3), 'NotRelevant', 
                             ifelse((df.final2$dp_ndp == 'TrialNotViewed'), 'TrialNotViewed',
                                    ifelse((df.final2$dp_ndp == 'ndp'), 'NDPV', #ndpv is Not DisProportional Viewing trial
                                           ifelse((df.final2$dp_ndp == 'dp' & df.final2$CorrectIA =="corIA"), "DPM", #DPM is DisProportional viewing to Match (i.e. correct face)
                                                  ifelse((df.final2$dp_ndp == 'dp' & df.final2$CorrectIA =="incorIA"), "DPNM", 'Error'))))) #DPNM is DisProportional viewing to Non-Match (i.e. incorrect face). Note this isn't absence of eye-mvmt memory but 'inappropriate' memory to wrong face. 


#sanity check #Count of DPM and DPNM and NDPV
a<-nrow(df.final2[grepl("DPNM", df.final2$dpm_dpnm),])
b<-nrow(df.final2[grepl("DPM", df.final2$dpm_dpnm),])
c<-nrow(df.final2[grepl("NDPV", df.final2$dpm_dpnm),])
d<-nrow(df.final2[grepl("TrialNotViewed", df.final2$dpm_dpnm),])
sum(a,b,c,d)

#tmrw Kirk, please make a table with number of each trial type for each subject. 

df.final.dpBins <- df.final2[!grepl("NotRelevant", df.final2$dpm_dpnm),]
df.final.dpBins <- df.final.dpBins[,c('Subj_trial','dpm_dpnm')]

write.csv(x = df.final.dpBins, "C:/Users/kgeier/28subj_timeperiods_viewingtype.csv", row.names = FALSE)


#############################################################################################################


pd <- position_dodge(0.1) # move them .05 to the left and right


#Graphs of dwell time to correct vs incorrect faces####

#df.final2$ChosCor4Opt <- paste0(df.final2$CorrectIA, "_", df.final2$ChosenIA)

#finding SE

df.final2chosen <- subset(df.final2, subset = df.final2$ChosenIA == 'chosenIA') 
#Removing NAs to avoid having them make issues with average
df.final2chosen <- df.final2chosen[!is.na(df.final2chosen$IA_DwellTime_Percent),]

d_sum3 = ddply(df.final2chosen, ~CorrectIA*IP_LABEL + SID,
               summarise, Average_Dwell_Time=mean(IA_DwellTime_Percent)) 
d_sum3all = ddply(df.final2chosen, ~CorrectIA + SID,
               summarise, Average_Dwell_Time=mean(IA_DwellTime_Percent)) 

d_sum3
#mmmm
d_sum3tp1 <- subset(d_sum3, subset = d_sum3$IP_LABEL  == "tp1")
d_sum3tp2<- subset(d_sum3, subset = d_sum3$IP_LABEL  == "tp2")
d_sum3tp3<- subset(d_sum3, subset = d_sum3$IP_LABEL  == "tp3")
d_sum3tp4<- subset(d_sum3, subset = d_sum3$IP_LABEL  == "tp4")
d_sum3tp5<- subset(d_sum3, subset = d_sum3$IP_LABEL  == "tp5")
d_sum3tp6<- subset(d_sum3, subset = d_sum3$IP_LABEL  == "tp6")
d_sum3tp1

detach("package:dplyr", unload=TRUE)
d_sum3tp1 = summarySEwithin(d_sum3tp1, 
                         measurevar=c("Average_Dwell_Time"), 
                         withinvars=c("CorrectIA"),
                         idvar=c("SID"))
d_sum3tp2 = summarySEwithin(d_sum3tp2, 
                         measurevar=c("Average_Dwell_Time"), 
                         withinvars=c("CorrectIA"),
                         idvar=c("SID"))
d_sum3tp3 = summarySEwithin(d_sum3tp3, 
                         measurevar=c("Average_Dwell_Time"), 
                         withinvars=c("CorrectIA"),
                         idvar=c("SID"))
d_sum3tp4 = summarySEwithin(d_sum3tp4, 
                         measurevar=c("Average_Dwell_Time"), 
                         withinvars=c("CorrectIA"),
                         idvar=c("SID"))
d_sum3tp5 = summarySEwithin(d_sum3tp5, 
                         measurevar=c("Average_Dwell_Time"), 
                         withinvars=c("CorrectIA"),
                         idvar=c("SID"))
d_sum3tp6 = summarySEwithin(d_sum3tp6, 
                         measurevar=c("Average_Dwell_Time"), 
                         withinvars=c("CorrectIA"),
                         idvar=c("SID"))
d_sum3all = summarySEwithin(d_sum3all, 
                            measurevar=c("Average_Dwell_Time"), 
                            withinvars=c("CorrectIA"),
                            idvar=c("SID"))

d_sum3tp1
d_sum3 <- rbind(d_sum3tp1, d_sum3tp2,d_sum3tp3,d_sum3tp4,d_sum3tp5,d_sum3tp6, d_sum3all)
TP_label <- c("tp1","tp1","tp2","tp2","tp3","tp3","tp4","tp4","tp5","tp5","tp6","tp6","all","all")
avgVStc <- c("tc","tc","tc","tc","tc","tc","tc","tc","tc","tc","tc","tc","all","all")
d_sum3 <-cbind(d_sum3, TP_label, avgVStc)
d_sum3

plot3_dwelltime <- ggplot(d_sum3, aes(x=TP_label, y=Average_Dwell_Time, group = interaction(CorrectIA, avgVStc), color = CorrectIA)) +
  geom_line(position=pd, size=1.5) +
  geom_errorbar(width=0, aes(ymin=Average_Dwell_Time-ci, ymax=Average_Dwell_Time+ci), colour="red",position=pd, size=1.5) +
 geom_errorbar(width=0, aes(ymin=Average_Dwell_Time-ci, ymax=Average_Dwell_Time+ci), data=d_sum3, position=pd, size=1.5) +
  #  geom_point() + 
  geom_hline(aes(yintercept=33.333), linetype= "dashed", size = 0.2) +
  scale_color_brewer(palette='Set1') +
  annotate("text", x = "tp2", y = 49, label = "*", size = 8)+
  annotate("text", x = "all", y = 41, label = "*", size = 8)+
  annotate("text", x = "tp3", y = 45, label = "+", size = 8)
plot3_dwelltime

df1 <- data.frame(a = c(1, "tp2", 3))


d_sum4 = ddply(df.final2chosen, ~ CorrectIA   + SID*IP_LABEL,
               summarise,Average_Dwell_Time=mean(IA_DwellTime_Percent))
d_sum4
#Remove following lines if want NoGuess and noresponse. Then use the post-hoc testing below
# d_sum4 <- subset(d_sum4, subset = ChosenIA %in%  c("chosenIA","NotChosenIA"))
# d_sum4
#again trying anova (2 way) bc x is categorical (but this time 2 way bc it's correct vs incor AND Chosen vs not chosen)
res_DT <- aov(Average_Dwell_Time~CorrectIA*IP_LABEL,d_sum4)  
summary(res_DT)  

#is this working? 
with(d_sum4, pairwise.t.test(Average_Dwell_Time, ChosenIA, 
                             p.adj="holm", paired=F))

#time for t teeeeessts! Mult comp correction : avg. will be uncorrected. tp1-tp6 will be Holm-Bonferroni corrected

incorDwell <- subset(d_sum4, subset = CorrectIA == "incorIA")
incorDwell
corDwell <- subset(d_sum4, subset = CorrectIA == "corIA")
corDwell
avgt <- t.test(incorDwell$Average_Dwell_Time,corDwell$Average_Dwell_Time,paired=T)
avgt


#do with s apply instead of repating 6 times...
incorDwell1 <- subset(d_sum4, subset = CorrectIA == "incorIA"& IP_LABEL == "tp1")
corDwell1 <- subset(d_sum4, subset = CorrectIA == "corIA" & IP_LABEL == "tp1" )
tp1t <- t.test(incorDwell1$Average_Dwell_Time,corDwell1$Average_Dwell_Time,paired=TRUE)

 
incorDwell2 <- subset(d_sum4, subset = CorrectIA == "incorIA"& IP_LABEL == "tp2")
corDwell2 <- subset(d_sum4, subset = CorrectIA == "corIA" & IP_LABEL == "tp2" )
tp2t <- t.test(incorDwell2$Average_Dwell_Time,corDwell2$Average_Dwell_Time,paired=TRUE)
 

incorDwell3 <- subset(d_sum4, subset = CorrectIA == "incorIA"& IP_LABEL == "tp3")
corDwell3 <- subset(d_sum4, subset = CorrectIA == "corIA" & IP_LABEL == "tp3" )
tp3t <- t.test(incorDwell3$Average_Dwell_Time,corDwell3$Average_Dwell_Time,paired=TRUE)
 

incorDwell4 <- subset(d_sum4, subset = CorrectIA == "incorIA"& IP_LABEL == "tp4")
corDwell4 <- subset(d_sum4, subset = CorrectIA == "corIA" & IP_LABEL == "tp4" )
tp4t <- t.test(incorDwell4$Average_Dwell_Time,corDwell4$Average_Dwell_Time,paired=TRUE)
 

incorDwell5 <- subset(d_sum4, subset = CorrectIA == "incorIA"& IP_LABEL == "tp5")
corDwell5 <- subset(d_sum4, subset = CorrectIA == "corIA" & IP_LABEL == "tp5" )
tp5t <- t.test(incorDwell5$Average_Dwell_Time,corDwell5$Average_Dwell_Time,paired=TRUE)


incorDwell6 <- subset(d_sum4, subset = CorrectIA == "incorIA"& IP_LABEL == "tp6")
corDwell6 <- subset(d_sum4, subset = CorrectIA == "corIA" & IP_LABEL == "tp6" )
tp6t <- t.test(incorDwell6$Average_Dwell_Time,corDwell6$Average_Dwell_Time,paired=TRUE)

tp1t
tp2t
tp3t
tp4t
tp5t 
tp6t 

#Holm-bonferroni correction (not done here just did it manually but can add soon) resulted in significance for tp2 only. 


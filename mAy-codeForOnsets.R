#making onsets for mAy! 
library("broom", lib.loc="~/R/R-3.3.2/library")
library("tidyr", lib.loc="~/R/R-3.3.2/library")
library("ggplot2", lib.loc="~/R/R-3.3.1/library")
library("gridExtra", lib.loc="~/R/R-3.3.1/library")
library("plyr")
library("dplyr", lib.loc="~/R/R-3.3.1/library")
library("data.table", lib.loc="~/R/R-3.3.2/library")



#Find How_to document either in the Dropbox\OlsenLabDocs\Useful Documents\How-Tos folder called README_RcodeMRIonsets
#Will also be copied at the bottom of this code as a comment. 

#import csv for eye-tracking. Don't forget to change the file name and path if you're using a different set of data!
df_study_full <- fread("C:/Users/kgeier/Google Drive/OlsenLab_cloud/All Experiments/MAYO/R studio analysis/1to28_studyreport_forOnsetfiles.csv")
df_test_full <- fread("C:/Users/kgeier/Google Drive/OlsenLab_cloud/All Experiments/MAYO/R studio analysis/1to28_testreport_forOnsetfiles.csv")

#Also need another dataframe with sync times from edf converter, 
synctime_df <- fread("C:/Users/kgeier/Google Drive/OlsenLab_cloud/All Experiments/MAYO/R studio analysis/Sync_time_subj1to28.csv")
#and EPI and etc scan numbers that need to be checked manually from Tech notes. 
ET_SB_conversion_df <- fread("C:/Users/kgeier/Google Drive/OlsenLab_cloud/All Experiments/MAYO/R studio analysis/ET_to_scanblock_conversion1to28.csv")
#and trial type as DPM, DPNM,NDPV, or TrialNotViewed. This gets generated by the mAy-FixationAnalysis_inR_Version2.R script  
viewingType_df <- fread("C:/Users/kgeier/Google Drive/OlsenLab_cloud/All Experiments/MAYO/R studio analysis/ViewingType28subj.csv")

#getting rid of blocks that didn't trigger, etc.
df_study_full <- df_study_full[df_study_full$RECORDING_SESSION_LABEL != "may14kg6"]
df_test_full <- df_test_full[df_test_full$RECORDING_SESSION_LABEL != "may14kg6"]
df_study_full <- df_study_full[df_study_full$RECORDING_SESSION_LABEL != "may05kg4"]
df_test_full <- df_test_full[df_test_full$RECORDING_SESSION_LABEL != "may05kg4"]

###- REMOVING 7 trials (both study and test) that had innapropriately had a mix of older/yonger/male/female in the lures. Ex. OF11, YF15, OF23. 
#trials were determined to be 
df_study_full2 <- subset(df_study_full, df_study_full$RECORDING_SESSION_LABEL != "may16kg1" | df_study_full$block != 1 | df_study_full$TRIAL_LABEL !="Trial: 15",)
df_study_full2 <- subset(df_study_full2, df_study_full2$RECORDING_SESSION_LABEL != "may16kg4" | df_study_full2$block != 4 | df_study_full2$TRIAL_LABEL !="Trial: 42",)
df_study_full2 <- subset(df_study_full2, df_study_full2$RECORDING_SESSION_LABEL != "may18kg3" | df_study_full2$block != 3 | df_study_full2$TRIAL_LABEL !="Trial: 25",)
df_study_full2 <- subset(df_study_full2, df_study_full2$RECORDING_SESSION_LABEL != "may04kg1" | df_study_full2$block != 1 | df_study_full2$TRIAL_LABEL !="Trial: 19",)
df_study_full2 <- subset(df_study_full2, df_study_full2$RECORDING_SESSION_LABEL != "may06kg3" | df_study_full2$block != 3 | df_study_full2$TRIAL_LABEL !="Trial: 35",)
df_study_full2 <- subset(df_study_full2, df_study_full2$RECORDING_SESSION_LABEL != "may11kg4" | df_study_full2$block != 4 | df_study_full2$TRIAL_LABEL !="Trial: 24",)
df_study_full2 <- subset(df_study_full2, df_study_full2$RECORDING_SESSION_LABEL != "may23kg4" | df_study_full2$block != 4 | df_study_full2$TRIAL_LABEL !="Trial: 36",)
#removing an 8th and 9th trial b.c. lure face matched target face
df_study_full2 <- subset(df_study_full2, df_study_full2$RECORDING_SESSION_LABEL != "may22kg6" | df_study_full2$block != 6 | df_study_full2$TRIAL_LABEL !="Trial: 9",)
df_study_full2 <- subset(df_study_full2, df_study_full2$RECORDING_SESSION_LABEL != "may13kg6" | df_study_full2$block != 6 | df_study_full2$TRIAL_LABEL !="Trial: 19",)
df_study_full <- df_study_full2

df_test_full2 <- subset(df_test_full, df_test_full$RECORDING_SESSION_LABEL != "may16kg1" | df_test_full$block != 1 | df_test_full$TRIAL_LABEL !="Trial: 60",)
df_test_full2 <- subset(df_test_full2, df_test_full2$RECORDING_SESSION_LABEL != "may16kg4" | df_test_full2$block != 4 | df_test_full2$TRIAL_LABEL !="Trial: 54",)
df_test_full2 <- subset(df_test_full2, df_test_full2$RECORDING_SESSION_LABEL != "may18kg3" | df_test_full2$block != 3 | df_test_full2$TRIAL_LABEL !="Trial: 55",)
df_test_full2 <- subset(df_test_full2, df_test_full2$RECORDING_SESSION_LABEL != "may04kg1" | df_test_full2$block != 1 | df_test_full2$TRIAL_LABEL !="Trial: 64",)
df_test_full2 <- subset(df_test_full2, df_test_full2$RECORDING_SESSION_LABEL != "may06kg3" | df_test_full2$block != 3 | df_test_full2$TRIAL_LABEL !="Trial: 60",)
df_test_full2 <- subset(df_test_full2, df_test_full2$RECORDING_SESSION_LABEL != "may11kg4" | df_test_full2$block != 4 | df_test_full2$TRIAL_LABEL !="Trial: 50",)
df_test_full2 <- subset(df_test_full2, df_test_full2$RECORDING_SESSION_LABEL != "may23kg4" | df_test_full2$block != 4 | df_test_full2$TRIAL_LABEL !="Trial: 53",)
#removing an 8th and 9thtrial b.c. lure face matched target face
df_test_full2 <- subset(df_test_full2, df_test_full2$RECORDING_SESSION_LABEL != "may22kg6" | df_test_full2$block != 6 | df_test_full2$TRIAL_LABEL !="Trial: 62",)
df_test_full2 <- subset(df_test_full2, df_test_full2$RECORDING_SESSION_LABEL != "may13kg6" | df_test_full2$block != 6 | df_test_full2$TRIAL_LABEL !="Trial: 53",)
df_test_full <- df_test_full2


#adding Subject ID to both study and test. 
df_test_full$SLabNoBlock <- df_test_full$RECORDING_SESSION_LABEL #makes a copy of session label to get SID from
df_test_full$SLabNoBlock <- substring(df_test_full$SLabNoBlock,1, 7) #selects first 7 charachters to remove block number
df_test_full$SID <- gsub("[^0-9]","",df_test_full$SLabNoBlock) #removes all non-numbers
df_test_full$SID <- gsub("146","14",df_test_full$SID) #subject 14 has issue with not being same #char. So this fixes it.
df_test_full$SID <- paste0("10",df_test_full$SID)
#still not working... use lines below   df_test_full$SID <- ifelse(grepl("may05bk", df_test_full$RECORDING_SESSION_LABEL), gsub("1005 ", "2001", df_test_full$SID), df_test_full$SID <- df_test_full$SID)
#replace next 6 lines with regular expressions when you have more time. 
df_test_full[df_test_full$SID == "1005" & df_test_full$RECORDING_SESSION_LABEL == "may05kg1", "SID"] <- "2001" #subject5a should be 2001 (artist).
df_test_full[df_test_full$SID == "1005" & df_test_full$RECORDING_SESSION_LABEL == "may05kg2", "SID"] <- "2001" #subject5a should be 2001 (artist).
df_test_full[df_test_full$SID == "1005" & df_test_full$RECORDING_SESSION_LABEL == "may05kg3", "SID"] <- "2001" #subject5a should be 2001 (artist).
df_test_full[df_test_full$SID == "1005" & df_test_full$RECORDING_SESSION_LABEL == "may05kg4", "SID"] <- "2001" #subject5a should be 2001 (artist).
df_test_full[df_test_full$SID == "1005" & df_test_full$RECORDING_SESSION_LABEL == "may05kg5", "SID"] <- "2001" #subject5a should be 2001 (artist).
df_test_full[df_test_full$SID == "1005" & df_test_full$RECORDING_SESSION_LABEL == "may05kg6", "SID"] <- "2001" #subject5a should be 2001 (artist).

unique(df_test_full$SID) #tests to see if this error prone process really worked (depends on if file name is properly done)
df_temmpcheck2 <- df_test_full %>% group_by(SID, block) %>% summarise(length(task))

#study
df_study_full$SLabNoBlock <- df_study_full$RECORDING_SESSION_LABEL #makes a copy of session label to get SID from
df_study_full$SLabNoBlock <- substring(df_study_full$SLabNoBlock,1, 7) #selects first 7 charachters to remove block number
df_study_full$SID <- gsub("[^0-9]","",df_study_full$SLabNoBlock) #removes all non-numbers
df_study_full$SID <- gsub("146","14",df_study_full$SID) #subject 14 has issue with not being same #char. So this fixes it.
df_study_full$SID <- paste0("10",df_study_full$SID)
df_study_full[df_study_full$SID == "1005" & df_study_full$RECORDING_SESSION_LABEL == "may05kg1", "SID"] <- "2001" #subject5b should be 2001 (artist).
df_study_full[df_study_full$SID == "1005" & df_study_full$RECORDING_SESSION_LABEL == "may05kg2", "SID"] <- "2001" #subject5b should be 2001 (artist).
df_study_full[df_study_full$SID == "1005" & df_study_full$RECORDING_SESSION_LABEL == "may05kg3", "SID"] <- "2001" #subject5b should be 2001 (artist).
df_study_full[df_study_full$SID == "1005" & df_study_full$RECORDING_SESSION_LABEL == "may05kg4", "SID"] <- "2001" #subject5b should be 2001 (artist).
df_study_full[df_study_full$SID == "1005" & df_study_full$RECORDING_SESSION_LABEL == "may05kg5", "SID"] <- "2001" #subject5b should be 2001 (artist).
df_study_full[df_study_full$SID == "1005" & df_study_full$RECORDING_SESSION_LABEL == "may05kg6", "SID"] <- "2001" #subject5b should be 2001 (artist).

unique(df_study_full$SID) #tests to see if this error prone process really worked (depends on if file name is properly done)
df_temmpcheck <- df_study_full %>% group_by(SID, block) %>% summarise(length(task))

###First checking accuracy ###
#convert response registered  (use this grepl syntax with reg exp for step above changing to 2001)
df_test_full$ConvertedTarget <- ifelse(grepl("0", df_test_full$FCRTargetIndex), "2",
                                     ifelse(grepl("1", df_test_full$FCRTargetIndex), "1",
                                            ifelse(grepl("2", df_test_full$FCRTargetIndex), "3","error")))
#Calculate Accuracy
df_test_full$Accuracy <- ifelse((df_test_full$ConvertedTarget==df_test_full$FCRResponse), "correct", 
                              ifelse((df_test_full$FCRResponse==-1), "no_response", 
                                     ifelse((df_test_full$FCRResponse==4), "NoGuess", "incorrect")))

#making a copy for later recomobination with study phase
df_test_full2 <- df_test_full

#making subset with just test rows.#This is also where trial 65s with "." as values get deleted 
df_test_full <- df_test_full[df_test_full$task==2,] 

#Getting percentage correct for each subject. Still not ideal because done partially manually so can't be faceted by block later :(
grouped_df_bySID <- group_by(df_test_full, SID)
temp.table.ga <- (table(grouped_df_bySID$SID, grouped_df_bySID$Accuracy))
temp.table.ga <- as.data.frame(temp.table.ga)
temp.wide.ga <- spread(temp.table.ga, Var2, Freq) 
temp.wide.ga$freqsum <- temp.wide.ga$correct + temp.wide.ga$incorrect +temp.wide.ga$NoGuess +temp.wide.ga$no_response
temp.wide.ga$Percent.Correct <- (temp.wide.ga$correct / temp.wide.ga$freqsum) * 100
temp.wide.ga$Percent.incorrect <- (temp.wide.ga$incorrect / temp.wide.ga$freqsum) * 100
temp.wide.ga$Percent.guess <- (temp.wide.ga$NoGuess / temp.wide.ga$freqsum) * 100
temp.wide.ga$Percent.no_response <- (temp.wide.ga$no_response / temp.wide.ga$freqsum) * 100
temp.wide.ga$Perc.cor.ifguessingall <- (temp.wide.ga$Percent.Correct + (temp.wide.ga$Percent.guess/3 + temp.wide.ga$Percent.no_response/3)) 
wide.df.percentAcc <- gather(temp.wide.ga, Percent.Correct, Percent.incorrect, Percent.guess, Percent.no_response, key = "Accuracy.type", value = "Percent")

#temp.wide.ga
#install.packages("xlsx")
library(xlsx) #load excel saving package
write.xlsx(x = temp.wide.ga, "C:/Users/kgeier/Google Drive/OlsenLab_cloud/All Experiments/MAYO/R studio analysis/Accuracyv6.xlsx",
           sheetName = "summary", row.names = FALSE)



####GraphingAccuracy####
# #graphs faceted by SID with total count.
pAccuracy <- ggplot(df_test_full, aes(factor(df_test_full$Accuracy))) +
  geom_bar() +
  facet_grid(.~ SID) +
  geom_hline(aes(yintercept=34.37, linetype="dashed", size=.02)) +
  theme_classic() +
  theme(axis.line.x= element_line(colour = "black"), axis.line.y = element_line(color= "black")) +
  theme(plot.title = element_text(size=36), axis.text = element_text(size=24),
        axis.title = element_text(size=32),
        strip.text.x = element_text(size=16)) +
  theme(legend.position="none") +
  ggtitle("placeholder") +
  xlab("Type of Accuracy") +
  ylab("Count Correct")
pAccuracy
 
#Getting rid of 2001
df_test_onlyYoungers <- df_test_full[df_test_full$SID != '2001',]
df_percent_acc <- temp.wide.ga[temp.wide.ga$Var1 != '2001',]

#changing names for graph to be prettier (surely theres a better way)
names(df_percent_acc)[names(df_percent_acc) == 'Percent.Correct'] <- 'Correct'
names(df_percent_acc)[names(df_percent_acc) == 'Percent.incorrect'] <- 'Incorrect'
names(df_percent_acc)[names(df_percent_acc) == 'Percent.guess'] <- "NoGuess"
names(df_percent_acc)[names(df_percent_acc) == 'Percent.no_response'] <- 'NoResponse'

#Converting to long format
df_percent_acc <- gather(df_percent_acc, Correct, Incorrect, NoGuess, NoResponse, key = "Accuracy.type", value = "Percent")

#plotting accuracy
  df_percent_acc %>%
  #gather(c(`R_MD`:`L_hipp`), key = "Name_of_ROI", value = "Tvalue") %>%
  #group_by(Var1) %>%
  ggplot(aes(x=factor(Accuracy.type), y= Percent, fill = (factor(Accuracy.type)))) + 
    geom_boxplot() +
   # geom_point() +
  #geom_point() +
  #geom_smooth(method = 'lm')+ 
  scale_x_discrete(limits=c("Correct","Incorrect","NoGuess", "NoResponse")) +
  #facet_wrap(~Var1, scales = 'free')
  geom_hline(aes(size = 12, yintercept=33.33, linetype="dashed")) +
    theme_classic() +
    #theme(axis.line.x= element_line(colour = "black"), axis.line.y = element_line(color= "black")) +
    #theme(plot.title = element_text(size=30), axis.text = element_text(size=20),
    #      axis.title = element_text(size=20)) +
    theme(legend.position ="none") +
    ggtitle("Response Distribution") +
    xlab("Response Type") +
    ylab("Percent of Responses") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size=56))




#this is framework code for later saving graphs, 1 for each subject faceted by block.
# jpeg(file = “C://R//SAVEHERE//myplot.jpeg”)
# #to build tomorrow into a for loop to go through all SID, make a faceted graph for each block for each subject to see how well they do.
# for(i in 15:22){
#   x = rnorm(500, 105, 10)
#   y = beta0 + beta1[i]*x + 15*rnorm(500)
#   
#   mypath <- file.path("C:","R","SAVEHERE",paste("myplot_", names[i], ".jpg", sep = ""))
#   
#   jpg(file=mypath)
#   mytitle = paste("my title is", names[i])
#   plot(x,y, main = mytitle)
#   dev.off()
# }


#### ONSET PART !!!!!!!! WOOOT#######################################################
#4- matching study and test IPs of double study and triple test dfs

#following method will take these two dataframes, add on a column with ordered numbers, 
# subset them so that study phase IP just has study phase trials and that the df with test IP just has test trials,
# then add them together and resort with the ordered numbers from original. This will get the correct fixation data 
# for each trial but combined appropriately.  See practice below for a simpler idea of what is going on. 

#filling in empty accuracy and target columns for study in case that's needed.
df_study_full <- cbind(df_study_full, 'NA')
df_study_full <- cbind(df_study_full, 'NA')

#appending column with seq 1:length of rows
df_study_withcount<- cbind(df_study_full, seq(1:nrow(df_study_full)))
df_test_withcount <-cbind(df_test_full2, seq(1:nrow(df_test_full2)))

#making sure the new column is named 
colnames(df_study_withcount)[23] <- "orderer"
colnames(df_study_withcount)[21] <- "ConvertedTarget"
colnames(df_study_withcount)[22] <- "Accuracy"
colnames(df_test_withcount)[23] <- "orderer"
#and numeric
df_study_withcount$orderer <- as.numeric(as.character(df_study_withcount$orderer))
df_test_withcount$orderer <- as.numeric(as.character(df_test_withcount$orderer))

#subsetting 
subset_study <-  df_study_withcount[grepl("1", df_study_withcount$task),]
subset_test <-  df_test_withcount[grepl("2", df_test_withcount$task),]

# df with 2 rows for each study phase and 3 for each test phase. This is for each stage of onset time. (ex. scene showing and face showing)
df_study_double <- subset_study[rep(1:nrow(subset_study),each=2),] 
df_test_triple <- subset_test[rep(1:nrow(subset_test),each=3),]

#rbind 
df_rebound <- rbind(df_study_double, df_test_triple) 

#re-ordering#############start here for regenerating df.final if screw it up
df.final <- df_rebound[order(df_rebound$orderer),]

##############
###5 Making column 3 trial
#Need to a) make name trial instead of Trial_Label and b) replace "Trial: #" with "#". 
#a)
df.final$TRIAL_LABEL <- gsub("Trial: ", "", df.final$TRIAL_LABEL) 
#b)
names(df.final)[names(df.final) == 'TRIAL_LABEL'] <- 'trial'

###6 Turning task into task phase (1 to study, 2 to test)
df.final$task <- gsub("1", "study", df.final$task) # 1 to study
df.final$task <- gsub("2", "test", df.final$task) #2 to test
names(df.final)[names(df.final) == 'task'] <- 'task_phase' #changing name of column to task_phase


###7 Making a new column as column 4. a) name it "cond" b) Concatenate together task"_scene_"face/prev/delay
#a)
df.final <- as.data.frame(append(df.final, list(cond = NA), after = 3))

#b)  
#Making new column with prev/face/delay
#Need an if statement to do ---- If task is 1, and value of trial row above is different, prev, if same, "face", elseif task is 2, and value of trial above is different, "prev", one above is same but two above is different, "delay", 2 above is same but 3 above is different, "face" 

#Really need a way to reference row by number and column by name... OTherwise very susceptible to changing columns. 
#elaborate for loop because I can't figure out how to reference cells 1 and 2 above otherwise
df.final$prev_face_delay <- 'placeholder'
for (i in 2:length(df.final$prev_face_delay)){
  if(df.final[i, 5]== 'study') {
    if(df.final[i, 3] == df.final[(i-1), 3]) { #line to check if trial is same as one above
      "face" -> df.final[i, 25]
    } else 
        "prev" -> df.final[i,25]
  } else if (df.final[i,5]== 'test') {
      if(df.final[i,3]== df.final[(i-2), 3]) {
        "face" -> df.final[i,25]
      } else if(df.final[i,3] == df.final[(i-1),3]) {
          "delay" -> df.final[i,25]
      } else
          "prev" -> df.final[i,25]
}}
df.final[1,25] <- "prev"
  
#Appending the 3 parts to make cond column like study_scene_face
df.final$cond <- paste0(df.final$task_phase,"_scene_", df.final$prev_face_delay)

###8 rename sceneimg to scene
names(df.final)[names(df.final) == 'scenesimg'] <- 'scene'

###9 faceimg to target_face
names(df.final)[names(df.final) == 'facesimg'] <- 'target_face' 

###10 rename lure faces to lure_face#
names(df.final)[names(df.final) == 'lure1'] <- "lure_face1"
names(df.final)[names(df.final) == 'lure2'] <- "lure_face2"

###11 rename response to plaus_rating and _REACTION_time_ to plaus rt(ms)
names(df.final)[names(df.final) == 'Response'] <- "plaus_rating"
names(df.final)[names(df.final) == 'X__Reaction_Time__1'] <- "plaus rt (ms)" 

###12 rename gender to sex
names(df.final)[names(df.final) == 'gender'] <- "sex" 

###13 Yay! Fun part. Need to match sync-time as a new column from the synctime_df by subject and block. 

df.final$Sync_time <- "dontknowyet" #clearing it to reduce confusion (sync-time there to start with )

#changing names for synctime_df so they match those of df.final
names(synctime_df)<- c("SID", "block","Sync_time") 
#merging the 2 dataframes by a vector of both column names in common, so that the new column, sync_time is matched to where subject and block are in common!
df.merged <- merge(df.final, synctime_df, by = c('SID', "block"))
df.final <- df.merged #################Another good spot for starting again with unscrewed up df.final####

###14 Calculating initial onset time of preview by subtracting sync time from IP start time

df.final$onset.time <- df.final$IP_START_TIME - df.final$Sync_time.y #note it's sync_time.y because there are 2 synctimes after merge  

#check for weird values by selecting all first and last trials 
df.subsetcheck <- df.final[(df.final$trial== "1"),c(1:4, 18, 27)] 
df.subsetcheck64 <- df.final[df.final$trial=="64",c(1:4, 18, 27)]
df.check1norep <- unique(df.subsetcheck)
df.check64norep <- unique(df.subsetcheck64)

range(df.check1norep$onset.time) #sholuld be all very close to 14000 ms
range(df.check64norep$onset.time) #sholuld be all very close to 589 000 ms
#Okay. Subject 16 Block 5 would not get converted from edf to asc. So rather than a true sync time, it has the average synctime of all other blocks. i.e. 14100ms. So may be up to 160 ms off. Waiting to get the original edf from Annette to see if that will work with ASC converter.  

###15 Getting duration in ms of each trial by type (ex. study_prev -> 1000ms and test_delay -> 7000ms)
df.final$Duration <- ifelse(grepl("study_scene_prev", df.final$cond), 1000,
                            ifelse(grepl("study_scene_face", df.final$cond), 2000,
                                   ifelse(grepl("test_scene_prev", df.final$cond), 1000,
                                          ifelse(grepl("test_scene_delay", df.final$cond), 7000,
                                                 ifelse(grepl("test_scene_face", df.final$cond), 3000, "error")))))

###16 Calculating onset time of each part of trial. Ex. study_prev at 14075 ms and study_face at 15075ms.  
df.final$OnsetByTrialPart <- ifelse(grepl("study_scene_prev", df.final$cond), (df.final$onset.time),
                            ifelse(grepl("study_scene_face", df.final$cond), (df.final$onset.time +1000),
                                   ifelse(grepl("test_scene_prev", df.final$cond), (df.final$onset.time),
                                          ifelse(grepl("test_scene_delay", df.final$cond), (df.final$onset.time + 1000),
                                                 ifelse(grepl("test_scene_face", df.final$cond), (df.final$onset.time + 8000), "error"))))) #8000ms because it's from display time of prev. so add 1000ms of preview display and 7000ms of delay = 8000s

#checking for errors
grep("error", df.final$Duration)
grep("error", df.final$OnsetByTrialPart)                    

###17 Calculating Offset time
#Converting new columns to numeric
df.final$OnsetByTrialPart <- as.numeric(as.character(df.final$OnsetByTrialPart))
df.final$Duration <- as.numeric(as.character(df.final$Duration))

df.final$offset <- (df.final$OnsetByTrialPart + df.final$Duration)

###18 face memory response and reaction time. a) change names b) NA to replace all studyphase cells and c) no_response to replace -1 
#a) Changing names of columns to match onset file format
names(df.final)[names(df.final)== 'FCRResponse'] <- 'faceMem_resp'
names(df.final)[names(df.final)== 'FCR_RT'] <- 'faceMem rt (ms)'

#b) NA to replace response and rt values for study #AGHHH. Should work but putting in values 1 higher 4 some reason. 

#Converting to numeric. -"why?" -"because for some damn reason I can't figure out if left as a factor the value has 1 added to it in the following code"
df.final$faceMem_resp <- as.numeric(as.character(df.final$faceMem_resp))
#Converting all non-"test" to "NA"
df.final$faceMem_resp <- ifelse(grepl('test', df.final$task_phase), (df.final$faceMem_resp), 'NA')

#Now for reaction time: Again convert to numeric for some reason..
df.final$`faceMem rt (ms)` <- as.numeric(as.character(df.final$`faceMem rt (ms)`))
#Converting all non-'test' to NA
df.final$`faceMem rt (ms)` <- ifelse(grepl('test', df.final$task_phase), (df.final$`faceMem rt (ms)`), 'NA')

#c) no_response to replace -1
#Response
df.final$faceMem_resp <- gsub("-1","no_response",df.final$faceMem_resp) 

#Reaction Time
df.final$`faceMem rt (ms)` <- gsub("-1","no_response",df.final$`faceMem rt (ms)`) 


###19 MATCH ACCURACY OF TEST PHASE TO STUDY PHASE. This is so you can tell which study phase trials they eventually got right and compare signals.
#2/3 of study trials will still be "NA" or "not tested"

#making new copy dataframe so don't screw up original. 
df.final2 <- df.final
#adding an "orderer" columns so can re-sort after dividing and appending. 
df.final2 <- cbind(df.final2, seq(1:nrow(df.final2)))
names(df.final2)[names(df.final2) == 'seq(1:nrow(df.final2))'] <- 'Orderer2'

#Making a df with just study
df.finalstudy <- df.final2[(df.final2$task_phase== "study"),]

#Making a df with just test and only the rows relevent to match with the study phase as well as Accuracy which is what we're trying to find out
df.finaltest <- df.final2[(df.final2$task_phase== "test"),]
df.finaltest2 <- df.finaltest[,c('SID', 'block', 'scene', 'Accuracy')]

#merge function. This is what we use to match. Takes the study phase and matches it to the test phase by Subj ID, block, and the scene shown. This will make a new accuracy column and add the value for accuracy within the test phase df to the same scene for the study phase df. 
#All.x=true is crucial to prevent deletion of study phase trials that weren't tested on. However it creates duplicates of all where they match. 
df.merged4acc <- merge(df.finalstudy, df.finaltest2, by = c("SID", 'block', 'scene'),  all.x = TRUE)

#Making the names match so can rbind study with subsequent accuracy to test. 
df.finaltest$Accuracy.y <- df.finaltest$Accuracy
names(df.finaltest)[names(df.finaltest) == 'Accuracy'] <- 'Accuracy.x'
names(df.finaltest)[names(df.finaltest) == 'Accuracy.y'] <- 'SubqAcc'
names(df.merged4acc)[names(df.merged4acc) == 'Accuracy.y'] <- 'SubqAcc'

#rbind to combine study and test
df.merged2 <- rbind(df.merged4acc, df.finaltest)
#Reorder by Orderer2 column to sort by SID, block, trial instead of all study then all test. 
df.final2 <- df.merged2[order(df.merged2$Orderer2),]
#remove duplicate rows created during merge function based on duplicate values in Orderer2 column
#note: still uses simpler version of [,31] which is more susceptible to issues with edits in other parts of code so maybe should change. 
df.duplicatesremoved <- subset(df.final2, !duplicated(df.final2[,31])) 

#Making it df.final again
df.final <-df.duplicatesremoved

df.studysumcheck <-df.final[(df.final$task_phase== "study"),]
df.testsumcheck <- df.final[(df.final$task_phase== "test"),]
df.studysumcheck %>% group_by(SubqAcc) %>% summarise(mean = length(SubqAcc))
df.testsumcheck %>% group_by(SubqAcc) %>% summarise(mean = length(SubqAcc))

#Replace all NAs in study phase with "not_tested".  
df.final$SubqAcc <- as.character(df.final$SubqAcc) #Need to change to string for some reason..
df.final$SubqAcc[is.na(df.final$SubqAcc)] <- "not_tested"

###20.Making a "Bad" column with 0 for okay trials and 1 for NoResponse trials or if something else is wrong with that trial
df.final$bad <- 0
df.final$bad[df.final$SubqAcc == "no_response"] <- "1"

###Adding EPI, ET and Scan block
#Making copy so don't screw up original
df.final2 <- df.final

#Eye tracking block
df.final2$`E-T_Block` <- df.final$block

#renaming convesion sheet to SID from subject so can merge
names(ET_SB_conversion_df)[names(ET_SB_conversion_df) == 'subject'] <- "SID"

#merging conversion sheet to df.final to match epi and scan blocks 
df.merged3 <- merge(df.final2, ET_SB_conversion_df, by = c('SID', "E-T_Block"))

#some checks to make sure nothing is duplicated nor lost by the merge
n_occur <- data.frame(table(df.merged3$SID, df.merged3$`E-T_Block`))
n_occur2 <- data.frame(table(df.final2$SID, df.final2$`E-T_Block`))
check <- cbind(n_occur, n_occur2) 

#Return to df.final
df.final <- df.merged3

#DPM####
#matching trial viewing type to each trial (DPM vs DPNM)
df.final2 <- df.final
df.final2$orderer3 <- 1:length(df.final2$SID)
df.final2$Subj_trial <- paste0(df.final2$SID, "_", df.final2$`E-T_Block`, "_", df.final2$trial)
#Removing sub18block6 from ViewingType_df. This isn't bc there's anything wrong with this block, but simply bc it's missing from the onset df and it has to be consistent or it will add new rows. Also. Subject 18 performed so poorly they are being excluded so it's not currently worth my time to save that block by re-generating a trial report. 
viewingType_df <- viewingType_df[!grepl("1018_6", viewingType_df$Subj_trial),]
df.final3 <- merge(df.final2, viewingType_df, by = 'Subj_trial', all = TRUE)
df.final3RO <- df.final3[order(df.final3$orderer3),]
df.final <- df.final3RO

###Final- Reorder to have exact spelling and order of columns. Specifically... 
#selecting only the right variables in the right order. But may have slightly different names
df.final2 <- df.final[,c("RECORDING_SESSION_LABEL", "SID", "E-T_Block", "Scan_Block", "EPI", "trial", "cond","task_phase","scene","target_face","lure_face1","lure_face2","plaus_rating","plaus rt (ms)","sex","age","OnsetByTrialPart","Duration","offset","faceMem_resp","faceMem rt (ms)","SubqAcc","bad", 'dpm_dpnm')]

#correcting the names that aren't quite right

names(df.final2)[names(df.final2) == 'SID'] <- "Subject"
names(df.final2)[names(df.final2) == 'OnsetByTrialPart'] <- 'onset'
names(df.final2)[names(df.final2) == 'offset'] <- "Offset Time (ms)"

#making an "other" column which is a SubqAcc column with 'other' instead of "no guess" and "no response
df.final2$AccWithOther <- ifelse(grepl("no_response", df.final2$SubqAcc), "other", 
                            ifelse(grepl("NoGuess", df.final2$SubqAcc), "other", df.final2$SubqAcc))

#Making a different other column with 1s for other and 0 for not other
df.final2$other <- ifelse(grepl("no_response", df.final2$SubqAcc), 1, 
                                 ifelse(grepl("NoGuess", df.final2$SubqAcc), 1, 0))

df.final<-df.final2 

#Replacing NA with 'study' in DPM column
df.final$dpm_dpnm[is.na(df.final$dpm_dpnm)] <- 'study'

#DID IT!
###JK still need to save each subject into it's own .csv file 
#write out portion ####
df.list = list(0)
#Making a new dataframe for each subject
for (i in 1001:1028) {
  assign(paste0("df.", i), df.final[df.final$Subject == i,])
}

#and 2001 #Next time use a repeat loop with a break function that can get to 1-24 and 2001 in one loop. But for now this works to save time
assign(paste0("df.", 2001), df.final[df.final$Subject == 2001,])

#Making a loop to save into excel files. 

#tediously adding each df to a list. Ideally add this to a loop as well but not enough time right now. 
df.list <- list(df.1001,df.1002, df.1003, df.1004, df.1005, df.1006, df.1007, df.1008, df.1009, df.1010, df.1011, df.1012, df.1013, df.1014, df.1015,df.1016,df.1017,df.1018,df.1019,df.1020,df.1021,df.1022,df.1023,df.1024,df.1025,df.1026, df.1027,df.1028,df.2001)

j = 1 
for (i in 1001:1028){
  mainDir <- "C:/Users/kgeier/Google Drive/OlsenLab_cloud/All Experiments/MAYO/R studio analysis/design.files"
  subDir <- i
  path <- file.path(mainDir,subDir)
  dir.create(path, showWarnings = FALSE)
  write.csv(x=df.list[[j]], file= paste0(path,"/","design.csv"), row.names = FALSE)
  j= j+1
}
#subj2001
mainDir <- "C:/Users/kgeier/Google Drive/OlsenLab_cloud/All Experiments/MAYO/R studio analysis/design.files"
subDir <- 2001
path <- file.path(mainDir,subDir)
dir.create(path, showWarnings = FALSE)
write.csv(x=df.list[[25]], file= paste0(path,"/","design.csv"), row.names = FALSE)

#writing file for df.final (design file not separated by subj)
mainDir <- "C:/Users/kgeier/Google Drive/OlsenLab_cloud/All Experiments/MAYO/R studio analysis/design.files_aug23"
subDir <- "all28"
path <- file.path(mainDir,subDir)
dir.create(path, showWarnings = FALSE)
write.csv(x=df.final, file= paste0(path,"/","designAll.csv"), row.names = FALSE)


#done writing out designs!####


#Making a graph to check onset time is as we expect ###### Everything works but the loop doesn't make the graphs... 
j = 1 
for (i in 1001:1028){
  df.temp4onsetplot <- df.list[[j]]
  df.temp4onsetplot$trial <- as.numeric(as.character(df.temp4onsetplot$trial)) 
  OnsetPlot <- ggplot(df.temp4onsetplot, aes(x = trial, y = onset)) +
    facet_wrap(~`E-T_Block`) +
    geom_point() +
    theme_classic() +
    theme(axis.line.x= element_line(colour = "black"), axis.line.y = element_line(color= "black")) +
    theme(plot.title = element_text(size=36), axis.text = element_text(size=12), 
          axis.title = element_text(size=32),
          strip.text.x = element_text(size=16)) +
    theme(legend.position="none") +
    ggtitle(j) 
  OnsetPlot
  j= j+1
}






















#README FILE COPY WITHOUT SCREENSHOTS
# Hey! 
#   So you're about to process the eye-tracking and behavioural data to generate onset times. 
# This is a big project. Some things may need to be done by hand which will be directed below. 
# 1.	Dataviewer:
#   a)	Open dataviewer and use file>import data> multiple data files to import the files of all blocks and subjects of your study. 
# b)	We're going to generate 2 reports that will be nearly identical but with different interest periods.  studyFull: Study_display to DISPLAY_fixation_cross
# testFull: test1_display_scene to DISPLAY_fixation_cross[2]
# Include the following variables in the following order.  
# RECORDING_SESSION_LABEL DATA_FILE IA_LABEL IA_DWELL_TIME IA_DWELL_TIME_% IA_FIRST_FIXATION_DURATION IA_FIXATION_% IA_FIXATION_COUNT IA_LAST_FIXATION_DURATION IP_LABEL IP_START_TIME FCRFaceList FCRResponse FCRTargetIndex FCR_RT Response __Reaction_Time__1 __Trial_Count__1 age block cb facesimg gender lure1 lure2 old_new sceneTime targetface task
# 
# 
# Save the generated reports as .csv into C:/Users/kgeier/Google Drive/OlsenLab_cloud/All Experiments/MAYO/R studio analysis/       1to24_studyreport_forOnsetfiles.csv  OR  1to24_testreport_forOnsetfiles.csv
# #If you have another path then change the path read in by the fread() function. 
# 
# 2.	You will also need to generate a file with sync times for each block. This will be used to subtract from IP start time to get the exact time the experiment was triggered (needed to align fMRI time to experiment). For mAy this needs to be done manually. If you have an Experiment Builder file that does this already you will need to EDIT this R script. Change the merge function under part #13.  
# A)	To generate this file follow the template img below with column names (subj, block, SyncTime). Then have every subject and block number in the first 2 columns. 
# B)	To find SyncTime use the SR research Utility called 'Visual EDF2ASC'. Open it and tediously load in each individual edf file (ex. May01kg1) and "convert". Once converted return to the original folder that contained the edf and it will now also contain and ASC file.
# C)	Find the row with "0 KEYBOARD_STUDY".  Record the associated time in the SyncTime column.  NOTE: This is very error prone process so a) be careful b)thoroughly double check later that the onset times match what you expect (for mAy study return values of 14000 ms from trigger until IP study trial 1. This was confirmed to be what the participant sees with a slow motion video). 
# 
# D)	Save this file as C:/Users/kgeier/Google Drive/OlsenLab_cloud/All Experiments/MAYO/R studio analysis/Sync_time_subj1to24.csv  (or change the path to your file appropriately). 
# 3.	Step 3 is to make a csv file that can link E-T block to scan block and EPI (example below). Get this information from the note the technologist gives you and double check the DICOM files on the server (once you download them and unzip them from RRiNId).  Save this file to C:/Users/kgeier/Google Drive/OlsenLab_cloud/All Experiments/MAYO/R studio analysis/ET_to_scanblock_conversion1to24.csv
# 
# 
# 4.	Edit R code! This is if you are doing new subjects or there is something different about your experiment. One tip is to change the subject numbers where relevant. And to decide which trials or blocks may need to be excluded and altering the first part after loading in the dataframes.  
# 5.	Now run the R code! 
#   a)	This will produce: a file called AccuracyV4.xlsx, and 25 folders from 1001:1024 and 2001 with design files inside, and a designAll.csv file with the design files of all subjects. 
# b)	You will also need to do double checks on the data. 
# For instance:
#   -how many blocks are missing
# -how many trials have no-responses
# -are any lures duplicated? Do lures match face? Are there inappropriate lures (ex. Old man with 2 young men)
# -do test targets match study scene face pairs?
# -Do faces or scenes repeat across blocks
# -are onset times what we expect. To do this, use the code at the very end of the R script. The loop isn't currently working to return the graphs to the plot area so either alter or go through one at a time for all of your subjects. Onset times should follow pattern below. 
# 
# 
# Complete!
  

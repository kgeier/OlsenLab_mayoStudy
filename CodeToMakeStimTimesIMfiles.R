#To transform onset files to stim_times files...

library("broom", lib.loc="~/R/R-3.3.2/library")
library("tidyr", lib.loc="~/R/R-3.3.2/library")
library("ggplot2", lib.loc="~/R/R-3.3.1/library")
library("gridExtra", lib.loc="~/R/R-3.3.1/library")
library("plyr")
library("stringr")
library("dplyr", lib.loc="~/R/R-3.3.1/library")
library("data.table", lib.loc="~/R/R-3.3.2/library")


TypeOfAccList <- c("SubqAcc", 'dpm_dpnm')
CondList <- c("study_scene_prev","test_scene_prev",  "test_scene_delay", "test_scene_face")
AccListResp <- c("not_tested", "correct", "incorrect", "NoGuess")
AccListEM <- c("study", "DPM", "DPNM", "NDPV")

#fake dataframe to fill empty df to be deleted later...
dfFake <- c(1:4)

#how many subjects do you have?
y <- 28

for (i in 1:y){
  data_dir = 'C:/Users/kgeier/Google Drive/OlsenLab_cloud/All Experiments/MAYO/R studio analysis/design.files.all.July11' 
  sub_dir = i + 1000 #our naming convention starts at 1001
  filename1 = 'design'
  fileextension = ".csv"
  filename = paste0(filename1, fileextension)
  filepath = paste(data_dir, '/', sub_dir, "/", filename, sep = ''); filepath
  
  df.1 <- fread(filepath)

  #Loops through whether we are classifying accuracy based on response or eye-movement
  for (j in TypeOfAccList){
    ifelse(j == "SubqAcc", m <- AccListResp, m <- AccListEM) #Sets m based on accuracy type to change which columns are used later

    #Getting rid of unneeded columns
    df.short <- df.1[,c('onset','cond','Scan_Block')]
    ifelse(j == "SubqAcc",  df.short <- cbind(df.short, df.1[,22]), df.short <- cbind(df.short, df.1[,24]))    
           
    #Subsetting just relevant condition
  for (k in CondList){  
   
    for (l in m){
   
    df.cond <- subset(df.short, subset = cond == k & get(j) == l)
    df.cond

    df.timeblock <- df.cond[,c('onset','Scan_Block')]
    #Milliseconds to Seconds
    df.timeblock$onset <- (df.timeblock$onset /1000)
    

    
    #Transposing to have each time set horizontal with new row for each block
    if(nrow(df.timeblock) != 0){
      df.timeblock%>% 
        group_by(Scan_Block) %>%
        mutate(time=seq_along(Scan_Block)) %>% spread(time,onset) -> df.sp

    #Getting rid of Scan_Block column
    df.final <- subset(df.sp, select = -(Scan_Block)) 
    
    
      
    #Saving as text file
    pathway4save <- 'C:/Users/kgeier/Desktop/StimTimeFiles/' # 'C:/Users/kgeier/Google Drive/OlsenLab_cloud/All Experiments/MAYO/R studio analysis/StimTimeFiles/'
    w = i + 1000
    filename4save <- paste0(w,"_",k, "_", l, ".1D")
    outputPath <- paste0(pathway4save, filename4save)
   # if(nrow(df.timeblock) != 0){
    write.table(df.final, outputPath, col.names = F, quote = F, row.names = F, sep="\t", na ="") 
} else 
    NULL
    }
    }
  }
  }
  

  



#Practice on one subject to see if it works before making loop

# #Reading in DF
# data_dir = 'C:/Users/kgeier/Google Drive/OlsenLab_cloud/All Experiments/MAYO/R studio analysis/design.files.all.July11' 
# sub_dir = '1001'
# filename1 = 'design'
# fileextension = ".csv"
# filename = paste0(filename1, fileextension)
# filepath = paste(data_dir, '/', sub_dir, "/", filename, sep = ''); filepath
# 
# df.1 <- fread(filepath)
# 
# #Getting rid of unneeded columns
# df.short <- df.1[,c('onset','cond','Scan_Block','SubqAcc','dpm_dpnm')]
# 
# #Making subset based on condition (ADD HERE ALL THE CONDITIONS)
# df.testcueCor <- subset(df.short, subset = cond == 'test_scene_prev' & SubqAcc == 'correct')
# df.timeblock <- df.testcueCor[,c('onset','Scan_Block')]
# 
# #Milliseconds to Seconds
# df.timeblock$onset <- (df.timeblock$onset /1000)
# 
# 
# #Transposing to have each time set horizontal with new row for each block
# df.timeblock%>% 
#   group_by(Scan_Block) %>%
#   mutate(time=seq_along(Scan_Block)) %>% spread(time,onset) -> df.sp
# 
# #Getting rid of Scan_Block column
# df.final <- subset(df.sp, select = -(Scan_Block)) 
# 
# 
# #Saving as text file
# write.table(df.final, 'C:/Users/kgeier/Google Drive/OlsenLab_cloud/All Experiments/MAYO/R studio analysis/StimTimeFiles/try.txt', col.names = F, quote = F, row.names = F,  ,sep=" ")

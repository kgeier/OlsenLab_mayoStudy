mayo_accuracy <- read.csv("C:/Users/kgeier/Dropbox/OlsenLabDocs/All Experiments at the Olsen Lab/MAYO/ay_accuracy_2subjectsRemoved.csv")
mayo_accuracy$Row.Labels <- factor(mayo_accuracy$Row.Labels)
library(tidyr)
library(ggplot2)
df_long <- gather(mayo_accuracy, condition, measurement, correct:NoResponse,factor_key = TRUE )
#df_long <- subset(df_long, Row.Labels!=1005)

df_long <- subset(df_long, condition!='GrandTotal')

ggplot(df_long, aes(x=condition, y=measurement)) + 
  geom_violin(aes(mapping=df_long$condition)) + 
  geom_point(aes(colour=df_long$condition), size= 2,
             stroke=1.7,
             position=position_jitter(width=.2,height=0))  +
  theme_classic() +
  theme(axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black")) +
  theme(plot.title = element_text(size=36),
        axis.text = element_text(size=24), 
        axis.title = element_text(size=32),
        strip.text.x = element_text(size=16)) +
  ylab("Number of trials") +
  xlab("Face-Scene accuracy") +
  #scale_y_discrete(breaks=seq(0, 100, 10)) +
  theme(legend.position="none")


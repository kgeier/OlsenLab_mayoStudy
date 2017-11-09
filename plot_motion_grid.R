# script to read in and plot motion
library(ggplot2)
library(dplyr)
library(tidyr)

df_scan1 <- read.table("/Volumes/172.24.4.65/1004/epi/rscan01-motpar.1D", quote="\"", comment.char="")
df_scan1$x <- seq(1,length(df_scan1$V1))
df_scan1_long <- gather(data=df_scan1, key=MotPar, value=shift, V1:V6)

df_scan2 <- read.table("/Volumes/172.24.4.65/1004/epi/rscan02-motpar.1D", quote="\"", comment.char="")
df_scan2$x <- seq(1,length(df_scan2$V1))
df_scan2_long <- gather(data=df_scan2, key=MotPar, value=shift, V1:V6)

df_scan3 <- read.table("/Volumes/172.24.4.65/1004/epi/rscan03-motpar.1D", quote="\"", comment.char="")
df_scan3$x <- seq(1,length(df_scan3$V1))
df_scan3_long <- gather(data=df_scan3, key=MotPar, value=shift, V1:V6)

df_scan4 <- read.table("/Volumes/172.24.4.65/1004/epi/rscan04-motpar.1D", quote="\"", comment.char="")
df_scan4$x <- seq(1,length(df_scan4$V1))
df_scan4_long <- gather(data=df_scan4, key=MotPar, value=shift, V1:V6)

df_scan5 <- read.table("/Volumes/172.24.4.65/1004/epi/rscan05-motpar.1D", quote="\"", comment.char="")
df_scan5$x <- seq(1,length(df_scan5$V1))
df_scan5_long <- gather(data=df_scan5, key=MotPar, value=shift, V1:V6)

df_scan6 <- read.table("/Volumes/172.24.4.65/1004/epi/rscan06-motpar.1D", quote="\"", comment.char="")
df_scan6$x <- seq(1,length(df_scan6$V1))
df_scan6_long <- gather(data=df_scan6, key=MotPar, value=shift, V1:V6)

df_scan1_long$MotPar <- factor(df_scan1_long)

p1 <- ggplot(df_scan1_long, aes(x=x, y=shift, color=MotPar)) +
  geom_point() +
  facet_wrap(~MotPar) +
  scale_fill_discrete() +
  geom_hline(aes(yintercept=3)) +
  geom_hline(aes(yintercept=-3))

p2 <- ggplot(df_scan2_long, aes(x=x, y=shift, color=MotPar)) +
  geom_point() +
  facet_wrap(~MotPar) +
  scale_fill_discrete() +
  geom_hline(aes(yintercept=3)) +
  geom_hline(aes(yintercept=-3))

p3 <- ggplot(df_scan3_long, aes(x=x, y=shift, color=MotPar)) +
  geom_point() +
  facet_wrap(~MotPar) +
  scale_fill_discrete() +
  geom_hline(aes(yintercept=3)) +
  geom_hline(aes(yintercept=-3))

p4 <- ggplot(df_scan4_long, aes(x=x, y=shift, color=MotPar)) +
  geom_point() +
  facet_wrap(~MotPar) +
  scale_fill_discrete() +
  geom_hline(aes(yintercept=3)) +
  geom_hline(aes(yintercept=-3))

p5 <- ggplot(df_scan5_long, aes(x=x, y=shift, color=MotPar)) +
  geom_point() +
  facet_wrap(~MotPar) +
  scale_fill_discrete() +
  geom_hline(aes(yintercept=3)) +
  geom_hline(aes(yintercept=-3))

p6 <- ggplot(df_scan6_long, aes(x=x, y=shift, color=MotPar)) +
  geom_point() +
  facet_wrap(~MotPar) +
  scale_fill_discrete() +
  geom_hline(aes(yintercept=3)) +
  geom_hline(aes(yintercept=-3))



p1

p2

p3

p4

p5

p6

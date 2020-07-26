install.packages("ctsem")
library(ctsem)

setwd("~/Benslab/CENT/veilahti")

df <- read.csv('./nfb/data/CENT_DB_2013-08-27/Sessions/auto_block.csv', sep=";", header=1)
str(df)

gudf <- df[df$reject_filter==0,]

boxplot(score ~ patient, data=df)
boxplot(score ~ patient, data=gudf)
boxplot(adj_score ~ patient, data=gudf)

boxplot(score ~ gametype, data=df)
boxplot(score ~ gametype, data=gudf)
boxplot(adj_score ~ gametype, data=gudf)

boxplot(score ~ trialtype, data=df)
boxplot(score ~ trialtype, data=gudf)
boxplot(adj_score ~ trialtype, data=gudf)


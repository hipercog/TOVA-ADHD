library(tidyverse)
library(viridis)
library(RColorBrewer)
source('R/znbnz_visuals.R')
source('R/LC_analysis_functions.R')
source('R/LC_visuals_functions.R')

odir <- "/home/bcowley/Dropbox/project_CENT/Publications/WIP - NFB Learning Curves/Figures/"

dat <- read.csv(file.path('data', 'scoreXregimeXlearner.csv'), sep = '\t')
df <- dat %>% filter(session_num > 2)

aggv <- 'norm_score'
byv <- 'session_num'
colv <- brewer.pal(n = 12, name = "Paired")
colv <- paste0(colv, "70")

lrnTB <- getTrialsMCI(filter(df, Regime == "TB", Learner == "Learner"), byv, aggv)
lrnSMR <- getTrialsMCI(filter(df, Regime == "SMR", Learner == "Learner"), byv, aggv)
nlrnTB <- getTrialsMCI(filter(df, Regime == "TB", Learner == "NonLearner"), byv, aggv)
nlrnSMR <- getTrialsMCI(filter(df, Regime == "SMR", Learner == "NonLearner"), byv, aggv)

# plot_4_meanCI(x1s, ml1, ci1, 
#               x2s, ml2, ci2, 
#               x3s, ml3, ci3, 
#               x4s, ml4, ci4, 
#               xlabel, ylabel, lgdstr, colgrp)
ylab <- "adjusted score"
xs <- seq(3,40)
plot_4_meanCI(xs, lrnTB$M[[aggv]], lrnTB$CI, 
              xs, nlrnTB$M[[aggv]], nlrnTB$CI, 
              xs, lrnSMR$M[[aggv]], lrnSMR$CI, 
              xs, nlrnSMR$M[[aggv]], nlrnSMR$CI, 
              "Session", ylab, 
              c("Learner:TB", "nonLearner:TB", "Learner:SMR", "nonLearner:SMR"), 
              colv[c(2, 1, 4, 3, 10, 9, 12, 11)],
              SMTH = TRUE, smthk = 3)
ggsave(paste0(odir, "learnerXregimeXperf.svg"))

library(tidyverse)
library(ggridges)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(RColorBrewer)
source('R/znbnVisuals.R')

odir <- "/home/bcowley/Dropbox/project_CENT/Publications/WIP - NFB Learning Curves/Figures"

colv <- brewer.pal(n = 12, name = "Paired")

#---- FIGURE FOR DES x LC SLOPE
des <- read.csv(file.path('data', 'ctsem_paper', 'DESxLC.csv'), sep = '\t')
dat <- read.csv(file.path('data', 'ctsem_paper', 'selfrepXLC.csv'))
srp <- read.csv(file.path('data', 'ctsem_paper', 'selfrep_TrainingPlaceboASRS.csv'))
srp <- select(srp, -c(1, 13:24))
srp$learner <- des$learner
srp$regime <- des$regime
srp$slope.score <- des$slope.score
srp$placebo_a1 <- srp$X40_Placebo_A1 - srp$PTre_Placebo_A1
srp$placebo_a2 <- srp$X40_Placebo_A2 - srp$PTre_Placebo_A2
srp$placebo_a3 <- srp$X40_Placebo_A3 - srp$PTre_Placebo_A3
srp$placebo_a4 <- srp$X40_Placebo_A3 - srp$PTre_Placebo_A4
srp$placebo_b1 <- srp$X40_Placebo_B1 - srp$PTre_Placebo_B1
srp$placebo_b2 <- srp$X40_Placebo_B2 - srp$PTre_Placebo_B2

t.test(Communication ~ learner, srp)
t.test(Technical ~ learner, srp)
t.test(Debriefing ~ learner, srp)
t.test(Stimulation_to_learning ~ learner, srp)
t.test(Patient_Communication ~ learner, srp)
t.test(Patient_Attitude.Behaviour ~ learner, srp)

t.test(Communication ~ regime, srp)
t.test(Technical ~ regime, srp)
t.test(Debriefing ~ regime, srp)
t.test(Stimulation_to_learning ~ regime, srp)
t.test(Patient_Communication ~ regime, srp)
t.test(Patient_Attitude.Behaviour ~ regime, srp)

corrMatrix("Communication", srp$Communication,
          "Technical", srp$Technical,
          "Debriefing", srp$Debriefing,
          "Stimulation_to_learning", srp$Stimulation_to_learning,
          "Patient_Communication", srp$Patient_Communication,
          "Patient_Attitude.Behaviour", srp$Patient_Attitude.Behaviour,
          "slope.score", srp$slope.score,
          title = "Patient-Trainer self-rep")

t.test(placebo_a1 ~ learner, srp)
t.test(placebo_a2 ~ learner, srp)
t.test(placebo_a3 ~ learner, srp)
t.test(placebo_a4 ~ learner, srp)
t.test(placebo_b1 ~ learner, srp)
t.test(placebo_b2 ~ learner, srp)

t.test(placebo_a1 ~ regime, srp)
t.test(placebo_a2 ~ regime, srp)
t.test(placebo_a3 ~ regime, srp)
t.test(placebo_a4 ~ regime, srp)
t.test(placebo_b1 ~ regime, srp)
t.test(placebo_b2 ~ regime, srp)

corrMatrix("Placebo A1", srp$placebo_a1,
           "Placebo A2", srp$placebo_a2,
           "Placebo A3", srp$placebo_a3,
           "Placebo A4", srp$placebo_a4,
           "Placebo B1", srp$placebo_b1,
           "Placebo B2", srp$placebo_b2,
           "slope.score", srp$slope.score,
           title = "Placebo 40-pre self-rep")

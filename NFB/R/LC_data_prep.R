#LIBRARIES for data management
library(plyr)
library(dplyr)
library(tidyr)
library(data.table)

setwd('~/Dropbox/PROJECT_CENT/Publications/WIP - Learning Curves/R')
source('./LC_analysis_functions.R')
source('./LC_visuals_functions.R')

pth <- '~/Dropbox/PROJECT_CENT/CENT_patients/CENT_DB_2013-08-27/Sessions/'
wkspace <- "LCanal.RData"
if(file.exists(wkspace)){
  load(wkspace)
}else{
  # Load the session, meta, daily-session, and outcome data
  ssn <- read.csv(paste0(pth, 'auto_session.csv'), sep=",")
  ssn <- ssn[order(ssn[,1]),]
  mta <- read.csv(paste0(pth, 'meta_data.csv'), sep=",")
  mta <- mta[order(mta[,1]),]
  dsl <- read.csv(paste0(pth, 'daily_session_wLatencyCalculation.csv'), sep=",")
  dsl <- dsl[order(dsl[,1]),]
  nfb <- read.csv('test_outcomes.csv', sep="\t")
  nfb <- nfb[order(nfb[,1]),]
  
  # Derive some indexing values from the datasets
  n.NOR <- 38
  n.INV <- 16
  n.TRA <- 8
  SBJS <- unique(nfb$Part_number)
  n.SBJS <- NCOL(SBJS) * NROW(SBJS)
  TB <- unique(subset(SBJS, nfb$TB.1.SMR1 == -1))
  SMR <- unique(subset(SBJS, nfb$TB.1.SMR1 == 1))
  MEN <- unique(subset(SBJS, nfb$Female0.Male2 == 2))
  WMN <- unique(subset(SBJS, nfb$Female0.Male2 == 0))
  
  # Get the most important data: trials/blocks, and make needed variables
  tr.raw <- read.csv(paste0(pth, 'auto_block.csv'), sep=";")
  tr.blk <- droplevels(subset(tr.raw, session > 1 & session < 41 & score != 0 & trainer_says_no != 1 & patient %in% SBJS))
  tr.blk <- tr.blk[order(tr.blk[,1]),]
  tr.blk$secs <- as.numeric(strptime(tr.blk$date_time, format='%Y:%m:%d %H:%M:%S'))
  
  # Subset the block data for analysis. Patients may have different start sessions for inverse, transfer.
  source('./LC_analysis_functions.R')
  trials <- getAllTrials(tr.blk)
  tr.nor <- getNorTrials(tr.blk, min_trial = 1) # Normal trials only, cut sessions with 1 trial
  tr.not <- getNotInvTrials(tr.blk, min_trial = 1) # not-inverse trials only, cut sessions with 1 trial
  tr.inv <- getInvTrials(tr.blk, min_trial = 1) # Inverse trials only, cut sessions with 1 trial
  tr.tra <- getTraTrials(tr.blk, min_trial = 1) # transfer trials only, cut sessions with 1 trial
  tr.nor0 <- getNorTrials(tr.blk) # Normal trials only
  tr.not0 <- getNotInvTrials(tr.blk) # not-inverse trials only
  tr.inv0 <- getInvTrials(tr.blk) # Inverse trials only
  tr.tra0 <- getTraTrials(tr.blk) # transfer trials only
  # tr.nor1 <- getNorTrials(tr.blk, n.NOR, 1) # Normal trials only, pad/trim, cut 1-trial sessions
  # tr.not1 <- getNotInvTrials(tr.blk, n.NOR, 1) # Normal trials only, pad/trim, cut 1-trial sessions
  # tr.inv1 <- getInvTrials(tr.blk, n.INV, 1) # Inverse trials only, pad/trim, cut 1-trial sessions
  # tr.tra1 <- getTraTrials(tr.blk, n.TRA, 1) # transfer trials only, pad/trim, cut 1-trial sessions
  
  # Create some subsets to work with individually
  ssn.tb.men <- droplevels(subset(ssn, patient%in%TB & patient%in%MEN))
  ssn.smr.men <- droplevels(subset(ssn, patient%in%SMR & patient%in%MEN))
  ssn.tb.wmn <- droplevels(subset(ssn, patient%in%TB & patient%in%WMN))
  ssn.smr.wmn <- droplevels(subset(ssn, patient%in%SMR & patient%in%WMN))
  
  dsl.tb <- droplevels(subset(dsl, patient%in%TB))
  dsl.smr <- droplevels(subset(dsl, patient%in%SMR))
  dsl.tb.men <- droplevels(subset(dsl, patient%in%TB & patient%in%MEN))
  dsl.smr.men <- droplevels(subset(dsl, patient%in%SMR & patient%in%MEN))
  dsl.tb.wmn <- droplevels(subset(dsl, patient%in%TB & patient%in%WMN))
  dsl.smr.wmn <- droplevels(subset(dsl, patient%in%SMR & patient%in%WMN))
  
  vigi_in <- read.table("~/Dropbox/PROJECT_CENT/Analysis/Vigilance/VIGALL_Intake.csv", header=TRUE, sep="\t", row.names=1)
  vigiPT <- subset(vigi_in, Group=="PT" & as.numeric(row.names(vigi_in)) %in% SBJS)
  vigiPWL <- subset(vigi_in, Group=="PWL")
  vigiC <- subset(vigi_in, Group=="C")

  DV <- 'adj_score'
  corvar <- 'secs'
  
  # centered to 0, scaled to [-1..1], geometric mean scores per session
  trial.scor <- getTrialCSGmean(trials, DV, cs = 'none')
  tr.nor.scr <- getTrialCSGmean(tr.nor, DV, cs = 'none')
  tr.not.scr <- getTrialCSGmean(tr.not, DV, cs = 'none')
  tr.inv.scr <- getTrialCSGmean(tr.inv, DV, cs = 'none')
  tr.tra.scr <- getTrialCSGmean(tr.tra, DV, cs = 'none')
  # testi <- gather(as.data.frame(tr.nor.scr), session, adj_score, 2:41, factor_key=TRUE)
  
  # Kendall correlations of per-trial adjusted score with order - Kendall is used for small N
  trial.cors <- getTrialCors(trials, DV, corvar)
  tr.nor.cor <- getTrialCors(tr.nor, DV, corvar)
  tr.not.cor <- getTrialCors(tr.not, DV, corvar)
  tr.inv.cor <- getTrialCors(tr.inv, DV, corvar)
  tr.tra.cor <- getTrialCors(tr.tra, DV, corvar)
}


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ---
## ## ---- TEST DATA CONSISTENCY ----
for (tst.num in SBJS){
  print(paste0("**** **** **** **** **** **** **** **** **** SUBJECT: ", tst.num))
  # tst.dsl <- subset(dsl, dsl$patient == tst.num)
  # tst.blk <- subset(tr.blk, tr.blk$patient == tst.num)
  # test <- as.data.frame(table(tst.blk$session))$Freq - tst.dsl[tst.dsl$trials > 0,]$trials
  
  # test <- length(unique(tr.tra[tr.tra$patient == tst.num,]$session))
  
  test <- tr.nor[tr.nor$patient == tst.num, c("patient", "session", "score", "nom_session")]
  print("#### #### NORMALS")
  print(test)
  test <- tr.inv[tr.inv$patient == tst.num, c("patient", "session", "score", "nom_session")]
  print("#### #### INVERSE")
  print(test)
  test <- tr.tra[tr.tra$patient == tst.num, c("patient", "session", "score", "nom_session")]
  print("#### #### TRANSFER")
  print(test)
  test <- tr.not[tr.not$patient == tst.num, c("patient", "session", "score", "nom_session")]
  print("#### #### NOT INVERSE (NORMAL+INVERSE)")
  print(test)
  # test <- tr.tra[tr.tra$patient == tst.num, c("patient", "session", "score", "nom_session")]
  # print(test)
}

testOnes(tr.nor)
testOnes(tr.nor0)
testOnes(tr.inv)
testOnes(tr.inv0)
testOnes(tr.tra)
testOnes(tr.tra0)
testOnes(tr.not)
testOnes(tr.not0)

ct.nor0 <- getTrialCount(tr.nor0)
ct.inv0 <- getTrialCount(tr.inv0)
ct.tra0 <- getTrialCount(tr.tra0)
ct.not0 <- getTrialCount(tr.not0)
ct.nor <- getTrialCount(tr.nor)
ct.inv <- getTrialCount(tr.inv)
ct.tra <- getTrialCount(tr.tra)
ct.not <- getTrialCount(tr.not)

apply(!(ct.nor0 == 1 | is.na(ct.nor0)), 1, sum)
apply(!(ct.inv0 == 1 | is.na(ct.inv0)), 1, sum)
apply(!(ct.tra0 == 1 | is.na(ct.tra0)), 1, sum)
apply(!(ct.nor == 1 | is.na(ct.nor)), 1, sum)
apply(!(ct.inv == 1 | is.na(ct.inv)), 1, sum)
apply(!(ct.tra == 1 | is.na(ct.tra)), 1, sum)


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ---
## ## ---- REPORT SIMPLE STATS OF LCs -------------------
# For daily session
summary(dsl.tb$adj_score)
summary(dsl.tb$score)
summary(dsl.smr$adj_score)
summary(dsl.smr$score)

# Print summary of ALL trials in ALL categories
reportTrials(trials, tr.nor, tr.inv, tr.tra, "adj_score", TB, SMR)



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ---
## ## ---- CALCULATE NON-PARAMETRIC LCs --------------------

## ## ---- MONOTONICITY of LEARNING ----
# Aggregation of Kendall correlations of per-trial adjusted score with order
trials.MLC <- getTrialAggCor(trials, DV, corvar)
tr.nor.MLC <- getTrialAggCor(tr.nor, DV, corvar)
tr.inv.MLC <- getTrialAggCor(tr.inv, DV, corvar)
tr.tra.MLC <- getTrialAggCor(tr.tra, DV, corvar)
tr.not.MLC <- getTrialAggCor(tr.not, DV, corvar)

# BY SESSIONS: Kendall correlations of daily session adjusted score with session number
ssn.MLC <- setDT(subset(dsl, patient %in% SBJS))[
  , .(cors = cor(get(DV), session, use="pair", method="k")), by=patient]$cors

ssn.trl.MLC <- getSsnMLC(trials, DV)
ssn.nor.MLC <- getSsnMLC(tr.nor0, DV)
ssn.inv.MLC <- getSsnMLC(tr.inv0, DV)
ssn.tra.MLC <- getSsnMLC(tr.tra0, DV)
ssn.not.MLC <- getSsnMLC(tr.not0, DV)

## MONOTONIC LEARNING
MonoLC <- aggMLC(trials.MLC, ssn.trl.MLC)
MLC.nor <- aggMLC(tr.nor.MLC, ssn.nor.MLC)
MLC.inv <- aggMLC(tr.inv.MLC, ssn.inv.MLC)
MLC.tra <- aggMLC(tr.tra.MLC, ssn.tra.MLC)
MLC.tra.fill <- aggMLC(tr.tra.MLC, ssn.tra.MLC, TRUE)
MLC.not <- aggMLC(tr.not.MLC, ssn.not.MLC)


## ## ---- CONSISTENCY of LEARNING ----

# BY TRIALS: cosine similarity of trial cors to an arbitrary hypothetical LC (AHLC)
IdealLC <- getTrialCosSim(trial.cors[,-(1)], trial.cors[,1], AHLC = "Fitts")
ILC.nor <- getTrialCosSim(tr.nor.cor[,-(1)], tr.nor.cor[,1], AHLC = "Fitts")
ILC.inv <- getTrialCosSim(tr.inv.cor[,-(1)], tr.inv.cor[,1], AHLC = "Fitts")
ILC.tra <- getTrialCosSim(tr.tra.cor[,-(1)], tr.tra.cor[,1], AHLC = "Fitts")
ILC.not <- getTrialCosSim(tr.not.cor[,-(1)], tr.not.cor[,1], AHLC = "Fitts")

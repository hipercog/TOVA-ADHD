library(tidyr)
library(dplyr)
library(data.table)

setwd('~/Dropbox/PROJECT_CENT/Publications/WIP - Learning Curves/R')
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
  tr.blk <- droplevels(subset(tr.raw, session_num > 1 & score != 0 & trainer_says_no != 1 & patient %in% SBJS))
  tr.blk <- tr.blk[order(tr.blk[,1]),]
  tr.blk$secs <- as.numeric(strptime(tr.blk$date_time, format='%Y:%m:%d %H:%M:%S'))
  
  # Subset the block data for analysis. Patients may have different start sessions for inverse, transfer.
  source('./LC_analysis_functions.R')
  trials <- getAllTrials(tr.blk)
  tr.nor <- getNorTrials(tr.blk, min_trial = 1) # Normal trials only, cut sessions with 1 trial
  tr.inv <- getInvTrials(tr.blk, min_trial = 1) # Inverse trials only, cut sessions with 1 trial
  tr.tra <- getTraTrials(tr.blk, min_trial = 1) # transfer trials only, cut sessions with 1 trial
  tr.nor0 <- getNorTrials(tr.blk) # Normal trials only
  tr.inv0 <- getInvTrials(tr.blk) # Inverse trials only
  tr.tra0 <- getTraTrials(tr.blk) # transfer trials only
  # tr.nor1 <- getNorTrials(tr.blk, n.NOR, 1) # Normal trials only, pad/trim, cut 1-trial sessions
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
}


## TEST DATA CONSISTENCY
for (tst.num in SBJS){
  print(paste0("**** **** **** **** **** **** **** **** **** SUBJECT: ", tst.num))
  # tst.dsl <- subset(dsl, dsl$patient == tst.num)
  # tst.blk <- subset(tr.blk, tr.blk$patient == tst.num)
  # test <- as.data.frame(table(tst.blk$session_num))$Freq - tst.dsl[tst.dsl$trials > 0,]$trials
  
  # test <- length(unique(tr.tra[tr.tra$patient == tst.num,]$session_num))
  
  test <- tr.nor[tr.nor$patient == tst.num, c("patient", "session_num", "score", "nom_session")]
  print("#### #### NORMALS")
  print(test)
  test <- tr.inv[tr.inv$patient == tst.num, c("patient", "session_num", "score", "nom_session")]
  print("#### #### INVERSE")
  print(test)
  test <- tr.tra[tr.tra$patient == tst.num, c("patient", "session_num", "score", "nom_session")]
  print("#### #### TRANSFER")
  print(test)
  # test <- tr.tra[tr.tra$patient == tst.num, c("patient", "session_num", "score", "nom_session")]
  # print(test)
}

testOnes(tr.nor)
testOnes(tr.nor0)
testOnes(tr.inv)
testOnes(tr.inv0)
testOnes(tr.tra)
testOnes(tr.tra0)

ct.nor0 <- getTrialCount(tr.nor0)
ct.inv0 <- getTrialCount(tr.inv0)
ct.tra0 <- getTrialCount(tr.tra0)
ct.nor <- getTrialCount(tr.nor)
ct.inv <- getTrialCount(tr.inv)
ct.tra <- getTrialCount(tr.tra)

apply(!(ct.nor0 == 1 | is.na(ct.nor0)), 1, sum)
apply(!(ct.inv0 == 1 | is.na(ct.inv0)), 1, sum)
apply(!(ct.tra0 == 1 | is.na(ct.tra0)), 1, sum)
apply(!(ct.nor == 1 | is.na(ct.nor)), 1, sum)
apply(!(ct.inv == 1 | is.na(ct.inv)), 1, sum)
apply(!(ct.tra == 1 | is.na(ct.tra)), 1, sum)


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ---
## ## ## ## ## --------------------- Simple stats of LCs --------------------
# For daily session
summary(dsl.tb$adj_score)
summary(dsl.tb$score)
summary(dsl.smr$adj_score)
summary(dsl.smr$score)

# Print summary of ALL trials in ALL categories
reportTrials(trials, tr.nor, tr.inv, tr.tra, "adj_score", TB, SMR)


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ---
## ## ## ## ## --------------------- non-para correlations of individual LCs --------------------

DV <- 'adj_score'
corvar <- 'secs'


### MONOTONICITY of LEARNING
# Aggregation of Kendall correlations of per-trial adjusted score with order
trials.MLC <- getTrialAggCor(trials, DV, corvar)
tr.nor.MLC <- getTrialAggCor(tr.nor, DV, corvar)
tr.inv.MLC <- getTrialAggCor(tr.inv, DV, corvar)
tr.tra.MLC <- getTrialAggCor(tr.tra, DV, corvar)

# BY SESSIONS: Kendall correlations of daily session adjusted score with session number
ssn.MLC <- setDT(subset(dsl, patient %in% SBJS))[
  , .(cors = cor(get(DV), session_num, use="pair", method="k")), by=patient]$cors

ssn.trl.MLC <- getSsnMLC(trials, DV)
ssn.nor.MLC <- getSsnMLC(tr.nor0, DV)
ssn.inv.MLC <- getSsnMLC(tr.inv0, DV)
ssn.tra.MLC <- getSsnMLC(tr.tra0, DV)

## MONOTONIC LEARNING
MonoLC <- aggMLC(trials.MLC, ssn.trl.MLC)
MLC.nor <- aggMLC(tr.nor.MLC, ssn.nor.MLC)
MLC.inv <- aggMLC(tr.inv.MLC, ssn.inv.MLC)
MLC.tra <- aggMLC(tr.tra.MLC, ssn.tra.MLC)
MLC.tra.fill <- aggMLC(tr.tra.MLC, ssn.tra.MLC, TRUE)


### CONSISTENCY of LEARNING
# Kendall correlations of per-trial adjusted score with order - Kendall is used for small N
trial.cors <- getTrialCors(trials, DV, corvar)
tr.nor.cor <- getTrialCors(tr.nor, DV, corvar)
tr.inv.cor <- getTrialCors(tr.inv, DV, corvar)
tr.tra.cor <- getTrialCors(tr.tra, DV, corvar)

# BY TRIALS: cosine similarity of trial cors to an arbitrary hypothetical LC (AHLC)
IdealLC <- getTrialCosSim(trial.cors[,-(1)], trial.cors[,1], AHLC = "Fitts")
ILC.nor <- getTrialCosSim(tr.nor.cor[,-(1)], tr.nor.cor[,1], AHLC = "Fitts")
ILC.inv <- getTrialCosSim(tr.inv.cor[,-(1)], tr.inv.cor[,1], AHLC = "Fitts")
ILC.tra <- getTrialCosSim(tr.tra.cor[,-(1)], tr.tra.cor[,1], AHLC = "Fitts")


# COMPARISON of MLC & ILC
cor(trials.MLC, IdealLC[,2])
cor(tr.nor.MLC, ILC.nor[,2])
cor(tr.inv.MLC, ILC.inv[,2])
cor(tr.tra.MLC[tr.tra.cor[,1] %in% ILC.tra[,1]], ILC.tra[,2])



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ---
## ## ## ## ## --------------------- TESTING ----

# MONOTONIC LEARNING
aov.MLC <- aov(MonoLC[,2] ~ nfb$TB.1.SMR1 * nfb$Female.2.Male2)
summary(aov.MLC)
lm.MLC <- lm(MonoLC[,2] ~ nfb$TB.1.SMR1 * nfb$Female.2.Male2)
summary(lm.MLC)
lm.MLC.nor <- lm(MLC.nor[,2] ~ nfb$TB.1.SMR1 * nfb$Female.2.Male2)
summary(lm.MLC.nor)
lm.MLC.inv <- lm(MLC.inv[,2] ~ nfb$TB.1.SMR1 * nfb$Female.2.Male2)
summary(lm.MLC.inv)
idx <- nfb$Part_number %in% MLC.tra[,1]
lm.MLC.tra <- lm(MLC.tra[,2] ~ nfb[idx,]$TB.1.SMR1 * nfb[idx,]$Female.2.Male2)
summary(lm.MLC.tra)

lm.ssn.MLC <- lm(ssn.trl.MLC[,2] ~ nfb$TB.1.SMR1 * nfb$Female.2.Male2)
summary(lm.ssn.MLC)
lm.trials.MLC <- lm(trials.MLC[,2] ~ nfb$TB.1.SMR1 * nfb$Female.2.Male2)
summary(lm.trials.MLC)


# IDEAL LEARNING
aov.ILC <- aov(IdealLC[,2] ~ nfb$TB.1.SMR1 * nfb$Female.2.Male2)
summary(aov.ILC)
lm.ILC <- lm(IdealLC[,2] ~ nfb$TB.1.SMR1 * nfb$Female.2.Male2)
summary(lm.ILC)
lm.ILC.nor <- lm(ILC.nor[,2] ~ nfb$TB.1.SMR1 * nfb$Female.2.Male2)
summary(lm.ILC.nor)
lm.ILC.inv <- lm(ILC.inv[,2] ~ nfb$TB.1.SMR1 * nfb$Female.2.Male2)
summary(lm.ILC.inv)
idx <- nfb$Part_number %in% ILC.tra[,1]
lm.ILC.tra <- lm(ILC.tra[,2] ~ nfb[idx,]$TB.1.SMR1 * nfb[idx,]$Female.2.Male2)
summary(lm.ILC.tra)



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ---
## ## ## ## ## --------------------- Extract mean and CIs for each data subset ----

# Organising variables
group_by_var <- 'nom_session'
agregate_var <- 'adj_score'
agg_var_str <- "adjusted score"


# GET MEAN AND CIs FOR ALL DATA...
allT <- getTrialsMCI(getAllTrials(tr.blk), group_by_var, agregate_var) # All trials
norT <- getTrialsMCI(getNorTrials(tr.blk), group_by_var, agregate_var) # Normal trials
invT <- getTrialsMCI(getInvTrials(tr.blk), group_by_var, agregate_var) # Inverse trials
traT <- getTrialsMCI(getTraTrials(tr.blk), group_by_var, agregate_var) # Transfer trials

# GET MEAN AND CIs FOR TB PROTOCOL...
df.tb <- droplevels(subset(tr.blk, patient%in%TB))
tb.allT <- getTrialsMCI(getAllTrials(df.tb), group_by_var, agregate_var) # All trials
tb.norT <- getTrialsMCI(getNorTrials(df.tb), group_by_var, agregate_var) # Normal trials
tb.invT <- getTrialsMCI(getInvTrials(df.tb), group_by_var, agregate_var) # Inverse trials
tb.traT <- getTrialsMCI(getTraTrials(df.tb), group_by_var, agregate_var) # Transfer trials

# ...AND FOR SMR!
df.smr <- droplevels(subset(tr.blk, patient%in%SMR))
smr.allT <- getTrialsMCI(getAllTrials(df.smr), group_by_var, agregate_var) # All trials
smr.norT <- getTrialsMCI(getNorTrials(df.smr), group_by_var, agregate_var) # Normal trials
smr.invT <- getTrialsMCI(getInvTrials(df.smr), group_by_var, agregate_var) # Inverse trials
smr.traT <- getTrialsMCI(getTraTrials(df.smr), group_by_var, agregate_var) # Transfer trials


# DERIVE MEAN AND CIs FOR THE MEN...
df.men <- droplevels(subset(tr.blk, patient%in%MEN))
men.allT <- getTrialsMCI(getAllTrials(df.men), group_by_var, agregate_var) # All trials
men.norT <- getTrialsMCI(getNorTrials(df.men), group_by_var, agregate_var) # Normal trials
men.invT <- getTrialsMCI(getInvTrials(df.men), group_by_var, agregate_var) # Inverse trials
men.traT <- getTrialsMCI(getTraTrials(df.men), group_by_var, agregate_var) # Transfer trials

# ...AND THE WOMEN!
df.wmn <- droplevels(subset(tr.blk, patient%in%TB & patient%in%WMN))
wmn.allT <- getTrialsMCI(getAllTrials(df.wmn), group_by_var, agregate_var) # All trials
wmn.norT <- getTrialsMCI(getNorTrials(df.wmn), group_by_var, agregate_var) # Normal trials
wmn.invT <- getTrialsMCI(getInvTrials(df.wmn), group_by_var, agregate_var) # Inverse trials
wmn.traT <- getTrialsMCI(getTraTrials(df.wmn), group_by_var, agregate_var) # Transfer trials


# DERIVE MEAN AND CIs BY PROTOCOL, FOR THE MEN...
df.tb.men <- droplevels(subset(tr.blk, patient%in%TB & patient%in%MEN))
men.tb.allT <- getTrialsMCI(getAllTrials(df.tb.men), group_by_var, agregate_var) # All trials
men.tb.norT <- getTrialsMCI(getNorTrials(df.tb.men), group_by_var, agregate_var) # Normal trials
men.tb.invT <- getTrialsMCI(getInvTrials(df.tb.men), group_by_var, agregate_var) # Inverse trials
men.tb.traT <- getTrialsMCI(getTraTrials(df.tb.men), group_by_var, agregate_var) # Transfer trials

df.smr.men <- droplevels(subset(tr.blk, patient%in%SMR & patient%in%MEN))
men.smr.allT <- getTrialsMCI(getAllTrials(df.smr.men), group_by_var, agregate_var) # All trials
men.smr.norT <- getTrialsMCI(getNorTrials(df.smr.men), group_by_var, agregate_var) # Normal trials
men.smr.invT <- getTrialsMCI(getInvTrials(df.smr.men), group_by_var, agregate_var) # Inverse trials
men.smr.traT <- getTrialsMCI(getTraTrials(df.smr.men), group_by_var, agregate_var) # Transfer trials

# ...AND THE WOMEN!
df.tb.wmn <- droplevels(subset(tr.blk, patient%in%TB & patient%in%WMN))
wmn.tb.allT <- getTrialsMCI(getAllTrials(df.tb.wmn), group_by_var, agregate_var) # All trials
wmn.tb.norT <- getTrialsMCI(getNorTrials(df.tb.wmn), group_by_var, agregate_var) # Normal trials
wmn.tb.invT <- getTrialsMCI(getInvTrials(df.tb.wmn), group_by_var, agregate_var) # Inverse trials
wmn.tb.traT <- getTrialsMCI(getTraTrials(df.tb.wmn), group_by_var, agregate_var) # Transfer trials

df.smr.wmn <- droplevels(subset(tr.blk, patient%in%SMR & patient%in%WMN))
wmn.smr.allT <- getTrialsMCI(getAllTrials(df.smr.wmn), group_by_var, agregate_var) # All trials
wmn.smr.norT <- getTrialsMCI(getNorTrials(df.smr.wmn), group_by_var, agregate_var) # Normal trials
wmn.smr.invT <- getTrialsMCI(getInvTrials(df.smr.wmn), group_by_var, agregate_var) # Inverse trials
wmn.smr.traT <- getTrialsMCI(getTraTrials(df.smr.wmn), group_by_var, agregate_var) # Transfer trials





## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ---
## ## ## ## ## --------------------- quadratic fit for men vs women, TB vs SMR protocols ----

quads <- matrix(NaN, nrow=23, ncol=39)
rownames(quads) <- pnums
sex = c("female", "male")
protocol = c("TB","null", "SMR")

# subplots for each participant, from daily session scores
par(mfrow=c(5,5))

df <- droplevels(subset(dsl.tb.men, exclude==0))
pnums <- unique(df$patient)
for (i in 1:length(pnums)){
  p <- pnums[i]
  tmp <- df #outgoners(df, df$adj_score, thr=3)
  tmp <- droplevels(subset(tmp, patient==p))
  quad <- with(tmp, fitted(lm(adj_score ~ session_num + I(session_num^2))))
  plot(tmp, main=p, col=2)
  lines(quad, col=3, lwd=2)
}
plot.new()

# subplots for each participant, from trial scores
par(mfrow=c(5,5))

df <- droplevels(subset(df.smr.wmn, reject_filter==0 & session_num>5 & norm_notInv==1 & trialtype==0))
pnums <- unique(df$patient)
# quads <- matrix(NaN, nrow=length(pnums), ncol=39)
for (i in 1:length(pnums)){
  p <- pnums[i]
  tmp <- df #outgoners(df, df$adj_score, thr=3)
  tmp <- ddply(subset(tmp, patient==p), .(session_num), summarize, adj_score=mean(adj_score))
  quad <- with(tmp, fitted(lm(adj_score ~ session_num + I(session_num^2))))
  # ttl = paste("Participant", p, ",", sex[nfb[SBJS==p,]$Gender..1.F..2.M.], ",", protocol[nfb[SBJS==p,]$TB1.SMR3])
  plot(tmp, main=p, xlab="", ylab="", col=2)
  lines(quad, col=3, lwd=2)
}
plot.new()



library(dplyr)
library(ggplot2)
library(plotrix)
library(plotly)

setwd('~/Dropbox/PROJECT_CENT/Publications/WIP - Learning Curves/R')
pth <- '~/Dropbox/PROJECT_CENT/CENT_patients/CENT_DB_2013-08-27/Sessions/'
wkspace <- "LCanal.Rdata"
if(file.exists(wkspace)){
  load(wkspace)
}else{
  blk <- read.csv(paste0(pth, 'auto_block.csv'), sep=";")
  blk <- droplevels(subset(blk, score!=0)) # & session_num<41))
  blk$secs <- as.numeric(strptime(blk$date_time, format='%Y:%m:%d %H:%M:%S'))
  blk <- blk[order(blk[,1]),]
  ssn <- read.csv(paste0(pth, 'auto_session.csv'), sep=",")
  ssn <- ssn[order(ssn[,1]),]
  mta <- read.csv(paste0(pth, 'meta_data.csv'), sep=",")
  mta <- mta[order(mta[,1]),]
  dsl <- read.csv(paste0(pth, 'daily_session_wLatencyCalculation.csv'), sep=",")
  dsl <- dsl[order(dsl[,1]),]
  nfb <- read.csv('test_outcomes.csv', sep="\t")
  nfb <- nfb[order(nfb[,1]),]
  sbjs <- unique(nfb$Part_number)
  tb <- unique(subset(nfb$Part_number, nfb$TB.1.SMR1 == -1))
  smr <- unique(subset(nfb$Part_number, nfb$TB.1.SMR1 == 1))
  men <- unique(subset(nfb$Part_number, nfb$Female0.Male2 == 2))
  wmn <- unique(subset(nfb$Part_number, nfb$Female0.Male2 == 0))
  
  ssn.tb.men <- droplevels(subset(ssn, patient%in%tb & patient%in%men))
  ssn.smr.men <- droplevels(subset(ssn, patient%in%smr & patient%in%men))
  ssn.tb.wmn <- droplevels(subset(ssn, patient%in%tb & patient%in%wmn))
  ssn.smr.wmn <- droplevels(subset(ssn, patient%in%smr & patient%in%wmn))
  
  dsl.tb <- droplevels(subset(dsl, patient%in%tb))
  dsl.smr <- droplevels(subset(dsl, patient%in%smr))
  dsl.tb.men <- droplevels(subset(dsl, patient%in%tb & patient%in%men))
  dsl.smr.men <- droplevels(subset(dsl, patient%in%smr & patient%in%men))
  dsl.tb.wmn <- droplevels(subset(dsl, patient%in%tb & patient%in%wmn))
  dsl.smr.wmn <- droplevels(subset(dsl, patient%in%smr & patient%in%wmn))
}

## TEST DATA CONSISTENCY
for (tst.num in sbjs){
  tst.dsl <- subset(dsl, dsl$patient == tst.num)
  tst.blk <- subset(blk, blk$patient == tst.num)
  test <- as.data.frame(table(tst.blk$session_num))$Freq - tst.dsl[tst.dsl$trials > 0,]$trials
  print(tst.num)
  print(test)
}



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ----
## ## ## ## ## --------------------- simple stats of daily session LCs --------------------
summary(dsl.tb$adj_score)
summary(dsl.tb$score)
summary(dsl.smr$adj_score)
summary(dsl.smr$score)



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ----
## ## ## ## ## --------------------- non-para correlations of individual LCs --------------------
# df <- droplevels(subset(dsl, exclude==0))
# Spearman correlations of daily session adjusted score with session number
ssn.cors <- vector("numeric", 23)
for (idx in 1:23){
  ssn.cors[idx] <- with(dsl[dsl$patient==sbjs[idx],], cor(adj_score, session_num, use="pair", method="spear"))
  # print(paste(ssn.cors[idx], sbjs[idx], nfb[nfb$Part_number==sbjs[idx],]$Gender..1.F..2.M.))
}
# ddply(dsl, dsl$patient, cor, x="adj_score", y="session_num", use="pairwise.complete.obs", method="spearman")
# Spearman correlations of per-trial adjusted score with order
use_days <- 37
trial.cors <- matrix(NA, 23, use_days)
print(paste("idx", "subj", "sessions", sep=" . . . "))

for (idx in 1:23){
  subj <- nfb$Part_number[idx]
  pins <- mta[mta$patient==subj,]$num_sessions
  print(paste(idx, subj, pins, sep=" . . . "))

  for (step in 1:use_days){
    day <- seq(pins-(use_days-1), pins)[step]
    # if (with(dsl[dsl$patient==subj & dsl$session_num==day,], trials - rejected_trials > 0)){
    # trials <- count(blk[blk$patient==subj & blk$session_num==day,])
    # print(trials)
    if (count(blk[blk$patient==subj & blk$session_num==day,]) > 0){
      trial.cors[idx, step] <- with(blk[blk$patient==subj & blk$session_num==day,], cor(adj_score, secs, use="pair", method="spear"))
    }
  }
}

trial.cosines <- apply(trial.cors, 1, cos.na, y=rep(1, use_days))[2,]
trial.angles <- apply(trial.cors, 1, angle.na, y=rep(1, use_days))

lrng <- (trial.cosines + ssn.cors) / 2


## Exploration Plots
par(mfrow=c(3,2))
LC <- ssn.cors
ttl <- "Session Spearman"
# DEMOGRAPHICS
plot(LC ~ nfb$Age, main = paste(ttl, "x age"))
boxplot(LC ~ nfb$Gender..1.F..2.M., main = paste(ttl, "x sex"), names = c("F", "M"))
boxplot(LC ~ nfb$TB.1.SMR1, main = paste(ttl, "x protocol"), names = c("TB", "SMR"))
plot(LC ~ apply(rbind(nfb$BAS.Drive, nfb$BAS.Fun.Seeking, nfb$BAS.Reward.Responsiveness), 2, mean), main = paste(ttl, "x BAS"))
plot(LC ~ nfb$BIS, main = paste(ttl, "x BIS"))
boxplot(LC ~ nfb$scales.COMORBID, main = paste(ttl, "x comorbidities"))#, names = c("0", "1", "2"))

# ADHD SELF-REPORTS
boxplot(LC ~ nfb$diagnosis..1.ADHD..2.ADD., main = paste(ttl, "x diagnosis"), names = c("ADHD", "ADD"))
boxplot(LC ~ nfb$diagnostic.COMORBID, main = paste(ttl, "x diagnostic comorbidity"))#, names = c("0", "1", "2"))
plot(LC ~ nfb$ASRS_total, main = paste(ttl, "x ASRS"))
plot(LC ~ nfb$ASRS_inattention_score, main = paste(ttl, "x ASRS-I"))
plot(LC ~ nfb$ASRS_hyperactivity.impulsivity_score, main = paste(ttl, "x ASRS-H"))
plot(LC ~ nfb$BADDS, main = paste(ttl, "x BADDS"))

# IQ
plot(LC ~ nfb$FSIQ, main = paste(ttl, "x IQ"))
plot(LC ~ nfb$VIQ, main = paste(ttl, "x verbal IQ"))
plot(LC ~ nfb$PIQ, main = paste(ttl, "x performance IQ"))
plot(LC ~ nfb$Span_total_standardized, main = paste(ttl, "x digit span"))
plot(LC ~ nfb$Span_fwd, main = paste(ttl, "x forward WM span"))
plot(LC ~ nfb$Span_bwd, main = paste(ttl, "x backward WM span"))

# VIGILANCE
# idx <- as.numeric(row.names(vigiPT)) %in% sbjs
# vigiPT <- subset(vigiPT, idx)
idx <- sbjs %in% as.numeric(row.names(vigiPT))
plot(LC[idx] ~ vigiPT$Index_W, main = paste(ttl, "x Vigi-lability"))
plot(LC[idx] ~ vigiPT$prW, main = paste(ttl, "x Vigi Waking-stage %"))
plot(LC[idx] ~ vigiPT$prA, main = paste(ttl, "x Vigi A-stage %"))
plot(LC[idx] ~ vigiPT$prB, main = paste(ttl, "x Vigi B-stage %"))
plot(LC[idx] ~ vigiPT$prC, main = paste(ttl, "x Vigi C-stage %"))
plot(LC[idx] ~ vigiPT$Art, main = paste(ttl, "x Vigi atefacts"))



# POLAR PLOT FOR VIEWING VECTOR ANGLES
mag <- rep(1, 23)
directionLabels <- c("")
colors <- c("black", "red", "green", "blue", "orange", "purple", "pink")
par(cex.axis=0.7)
par(cex.lab=0.5)

radial.plot(c(0, mag), 
            c(0, trial.angles), 
            lwd=2, line.col=(LC+1)*128,
            labels=directionLabels,
            radial.lim = c(0,1), #range of grid circle
            main="circular diagram!",
            show.grid.label=1, #put the concentric circle labels going down
            show.radial.grid=TRUE
)


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ----
## ## ## ## ## --------------------- BWRC COLOURS --------------------

# colors: 1 orange - 2 blue - 3 dark blue - 4 darker gray - 5 lighter gray - 6 green - 7 red
bwrc <- palette(c('#F0AA00','#00A0BE','#003C78','#828282','#F0F0F0','#B4BE00','#FA4100'))
bwrcol <- colorRampPalette(c('#828282','#F0F0F0','#00A0BE','#003C78'))
bwrcol <- colorRampPalette(c('#F0AA00','#00A0BE','#003C78'))
bwrcol <- colorRampPalette(c('#00A0BE','#F0AA00'))


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ----
## ## ## ## ## --------------------- quadratic fit for men vs women, TB vs SMR protocols ----

quads <- matrix(NaN, nrow=23, ncol=39)
rownames(quads) <- pnums
sex = c("female", "male")
protocol = c("TB","null", "SMR")

# subplots for protocol X sex
par(mfrow=c(2,2))

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
  # ttl = paste("Participant", p, ",", sex[nfb[nfb$Part_number==p,]$Gender..1.F..2.M.], ",", protocol[nfb[nfb$Part_number==p,]$TB1.SMR3])
  plot(tmp, main=p, xlab="", ylab="", col=2)
  lines(quad, col=3, lwd=2)
}
plot.new()

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ----
## ## ## ## ## --------------------- trend of men vs women, TB vs SMR protocol LCs with CIs ----

df.tb.men <- droplevels(subset(blk, patient%in%tb & patient%in%men))
# All trials
df <- droplevels(subset(df.tb.men, reject_filter==0 & session_num>1))
men.tb.mnsc <- ddply(df, .(session_num), summarize, sesc=mean(adj_score))
men.tb.scCI <- ddply(df, .(session_num), bootin, varname="adj_score")
# Normal trials
df <- droplevels(subset(df.tb.men, reject_filter==0 & session_num>1 & norm_notInv==1 & trialtype==0))
men.nor.tb.mnsc <- ddply(df, .(session_num), summarize, sesc=mean(adj_score))
men.nor.tb.scCI <- ddply(df, .(session_num), bootin, varname="adj_score")
# Inverse trials
df <- droplevels(subset(df.tb.men, reject_filter==0 & session_num>22 & norm_notInv==0 & trialtype==0))
men.inv.tb.mnsc <- ddply(df, .(session_num), summarize, sesc=mean(adj_score))
men.inv.tb.scCI <- ddply(df, .(session_num), bootin, varname="adj_score")
# Transfer trials
df <- droplevels(subset(df.tb.men, reject_filter==0 & session_num>30 & norm_notInv==1 & trialtype>0))
men.tfr.tb.mnsc <- ddply(df, .(session_num), summarize, sesc=mean(adj_score))
men.tfr.tb.scCI <- ddply(df, .(session_num), bootin, varname="adj_score")


df.smr.men <- droplevels(subset(blk, patient%in%smr & patient%in%men))
# All trials
df <- droplevels(subset(df.smr.men, reject_filter==0 & session_num>1))
men.smr.mnsc <- ddply(df, .(session_num), summarize, sesc=mean(adj_score))
men.smr.scCI <- ddply(df, .(session_num), bootin, varname="adj_score")
# Normal trials
df <- droplevels(subset(df.smr.men, reject_filter==0 & session_num>1 & norm_notInv==1 & trialtype==0))
men.nor.smr.mnsc <- ddply(df, .(session_num), summarize, sesc=mean(adj_score))
men.nor.smr.scCI <- ddply(df, .(session_num), bootin, varname="adj_score")
# Inverse trials
df <- droplevels(subset(df.smr.men, reject_filter==0 & session_num>22 & norm_notInv==0 & trialtype==0))
men.inv.smr.mnsc <- ddply(df, .(session_num), summarize, sesc=mean(adj_score))
men.inv.smr.scCI <- ddply(df, .(session_num), bootin, varname="adj_score")
# Transfer trials
df <- droplevels(subset(df.smr.men, reject_filter==0 & session_num>30 & norm_notInv==1 & trialtype>0))
men.tfr.smr.mnsc <- ddply(df, .(session_num), summarize, sesc=mean(adj_score))
men.tfr.smr.scCI <- ddply(df, .(session_num), bootin, varname="adj_score")


# ALTOGETHER NOW - master plot
xtk = seq(2,40)
plot(range(xtk),
     range(c(men.nor.smr.scCI$V4, men.nor.smr.scCI$V5, men.inv.smr.scCI$V4, men.inv.smr.scCI$V5, 
             men.nor.tb.scCI$V4, men.nor.tb.scCI$V5, men.inv.tb.scCI$V4, men.inv.tb.scCI$V5)), 
     type="n", xaxt="n", xlab="all sessions, males", ylab="adjusted score")
axis(1, at=xtk)

polygon(c(xtk, rev(xtk)), c(men.nor.tb.scCI$V5, rev(men.nor.tb.scCI$V4)), col=2, border=NA)
polygon(c(xtk, rev(xtk)), c(men.nor.smr.scCI$V5, rev(men.nor.smr.scCI$V4)), col=2, border=NA)

lines(men.nor.tb.mnsc$session_num, men.nor.tb.mnsc$sesc, type="l", lwd=2, lty=2, col=3)
lines(men.nor.smr.mnsc$session_num, men.nor.smr.mnsc$sesc, type="l", lwd=2, col=3)

lines(men.inv.tb.mnsc$session_num, men.inv.tb.mnsc$sesc, type="l", lwd=2, lty=2, col='black')
lines(men.inv.smr.mnsc$session_num, men.inv.smr.mnsc$sesc, type="l", lwd=2, col='black')

with(men.inv.tb.scCI, lines(session_num, V4, type="l", lwd=1, lty=2, col='black'))
with(men.inv.tb.scCI, lines(session_num, V5, type="l", lwd=1, lty=2, col='black'))
with(men.inv.smr.scCI, lines(session_num, V4, type="l", lwd=1, col='black'))
with(men.inv.smr.scCI, lines(session_num, V5, type="l", lwd=1, col='black'))

legend("topleft", inset=c(0.05,0.0), y.intersp=1, cex=1, bty="n",
       legend=c("normal TB", "normal SMR", "inverse TB", "inverse SMR"), 
       lty=c(2,1,2,1), lwd=2, col=c(3,3,'black','black'))


 # ...AND THE WOMEN!
df.tb.wmn <- droplevels(subset(blk, patient%in%tb & patient%in%wmn))
# All trials
df <- droplevels(subset(df.tb.wmn, reject_filter==0 & session_num>1))
wmn.tb.mnsc <- ddply(df, .(session_num), summarize, sesc=mean(adj_score))
wmn.tb.scCI <- ddply(df, .(session_num), bootin, varname="adj_score")
# Normal trials
df <- droplevels(subset(df.tb.wmn, reject_filter==0 & session_num>1 & norm_notInv==1 & trialtype==0))
wmn.nor.tb.mnsc <- ddply(df, .(session_num), summarize, sesc=mean(adj_score))
wmn.nor.tb.scCI <- ddply(df, .(session_num), bootin, varname="adj_score")
# Inverse trials
df <- droplevels(subset(df.tb.wmn, reject_filter==0 & session_num>22 & norm_notInv==0 & trialtype==0))
wmn.inv.tb.mnsc <- ddply(df, .(session_num), summarize, sesc=mean(adj_score))
wmn.inv.tb.scCI <- ddply(df, .(session_num), bootin, varname="adj_score")
# Transfer trials
df <- droplevels(subset(df.tb.wmn, reject_filter==0 & session_num>30 & norm_notInv==1 & trialtype>0))
wmn.tfr.tb.mnsc <- ddply(df, .(session_num), summarize, sesc=mean(adj_score))
wmn.tfr.tb.scCI <- ddply(df, .(session_num), bootin, varname="adj_score")


df.smr.wmn <- droplevels(subset(blk, patient%in%smr & patient%in%wmn))
# All trials
df <- droplevels(subset(df.smr.wmn, reject_filter==0 & session_num>1))
wmn.smr.mnsc <- ddply(df, .(session_num), summarize, sesc=mean(adj_score))
wmn.smr.scCI <- ddply(df, .(session_num), bootin, varname="adj_score")
# Normal trials
df <- droplevels(subset(df.smr.wmn, reject_filter==0 & session_num>1 & norm_notInv==1 & trialtype==0))
wmn.nor.smr.mnsc <- ddply(df, .(session_num), summarize, sesc=mean(adj_score))
wmn.nor.smr.scCI <- ddply(df, .(session_num), bootin, varname="adj_score")
# Inverse trials
df <- droplevels(subset(df.smr.wmn, reject_filter==0 & session_num>22 & norm_notInv==0 & trialtype==0))
wmn.inv.smr.mnsc <- ddply(df, .(session_num), summarize, sesc=mean(adj_score))
wmn.inv.smr.scCI <- ddply(df, .(session_num), bootin, varname="adj_score")
# Transfer trials
df <- droplevels(subset(df.smr.wmn, reject_filter==0 & session_num>30 & norm_notInv==1 & trialtype>0))
wmn.tfr.smr.mnsc <- ddply(df, .(session_num), summarize, sesc=mean(adj_score))
wmn.tfr.smr.scCI <- ddply(df, .(session_num), bootin, varname="adj_score")


# ALTOGETHER NOW - master plot
xtk = seq(2,40)
plot(range(xtk),
     range(c(wmn.nor.smr.scCI$V4, wmn.nor.smr.scCI$V5, wmn.inv.smr.scCI$V4, wmn.inv.smr.scCI$V5, 
             wmn.nor.tb.scCI$V4, wmn.nor.tb.scCI$V5, wmn.inv.tb.scCI$V4, wmn.inv.tb.scCI$V5)), 
     type="n", xaxt="n", xlab="all sessions, females", ylab="adjusted score")
axis(1, at=xtk)

polygon(c(xtk, rev(xtk)), c(wmn.nor.tb.scCI$V5, rev(wmn.nor.tb.scCI$V4)), col=2, border=NA)
polygon(c(xtk, rev(xtk)), c(wmn.nor.smr.scCI$V5, rev(wmn.nor.smr.scCI$V4)), col=2, border=NA)

lines(wmn.nor.tb.mnsc$session_num, wmn.nor.tb.mnsc$sesc, type="l", lwd=2, lty=2, col=3)
lines(wmn.nor.smr.mnsc$session_num, wmn.nor.smr.mnsc$sesc, type="l", lwd=2, col=3)

lines(wmn.inv.tb.mnsc$session_num, wmn.inv.tb.mnsc$sesc, type="l", lwd=2, lty=2, col='black')
lines(wmn.inv.smr.mnsc$session_num, wmn.inv.smr.mnsc$sesc, type="l", lwd=2, col='black')

with(wmn.inv.tb.scCI, lines(session_num, V4, type="l", lwd=1, lty=2, col='black'))
with(wmn.inv.tb.scCI, lines(session_num, V5, type="l", lwd=1, lty=2, col='black'))
with(wmn.inv.smr.scCI, lines(session_num, V4, type="l", lwd=1, col='black'))
with(wmn.inv.smr.scCI, lines(session_num, V5, type="l", lwd=1, col='black'))

legend("topleft", inset=c(0.05,0.0), y.intersp=1, cex=1, bty="n",
       legend=c("normal TB", "normal SMR", "inverse TB", "inverse SMR"), 
       lty=c(2,1,2,1), lwd=2, col=c(3,3,'black','black'))



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ----
## ## ## ## ## --------------------- trend of TB vs SMR protocol LCs with CIs ----

df.tb <- droplevels(subset(blk, patient%in%tb))
# All trials
df <- droplevels(subset(df.tb, reject_filter==0 & session_num>1))
tb.mnsc <- ddply(df, .(session_num), summarize, sesc=mean(adj_score))
tb.scCI <- ddply(df, .(session_num), bootin, varname="adj_score")
# Normal trials
df <- droplevels(subset(df.tb, reject_filter==0 & session_num>1 & norm_notInv==1 & trialtype==0))
nor.tb.mnsc <- ddply(df, .(session_num), summarize, sesc=mean(adj_score))
nor.tb.scCI <- ddply(df, .(session_num), bootin, varname="adj_score")
# Inverse trials
df <- droplevels(subset(df.tb, reject_filter==0 & session_num>22 & norm_notInv==0 & trialtype==0))
inv.tb.mnsc <- ddply(df, .(session_num), summarize, sesc=mean(adj_score))
inv.tb.scCI <- ddply(df, .(session_num), bootin, varname="adj_score")
# Transfer trials
df <- droplevels(subset(df.tb, reject_filter==0 & session_num>30 & norm_notInv==1 & trialtype>0))
tfr.tb.mnsc <- ddply(df, .(session_num), summarize, sesc=mean(adj_score))
tfr.tb.scCI <- ddply(df, .(session_num), bootin, varname="adj_score")

# master plot
xtk = seq(2,40)
plot(range(xtk),
     range(c(nor.tb.scCI$V4, nor.tb.scCI$V5, inv.tb.scCI$V4, inv.tb.scCI$V5)), 
     type="n", xaxt="n", xlab="all TB sessions", ylab="adjusted score")
axis(1, at=xtk)
polygon(c(xtk, rev(xtk)), c(nor.tb.scCI$V5, rev(nor.tb.scCI$V4)), col=5, border=NA)
lines(nor.tb.mnsc$session_num, nor.tb.mnsc$sesc, type="l", lwd=2, col=3)
lines(inv.tb.mnsc$session_num, inv.tb.mnsc$sesc, type="l", lwd=2, col='black')

with(nor.tb.scCI, lines(session_num, V4, type="l", lwd=1, col=2))
with(nor.tb.scCI, lines(session_num, V5, type="l", lwd=1, col=2))
with(inv.tb.scCI, lines(session_num, V4, type="l", lwd=1, col=4))
with(inv.tb.scCI, lines(session_num, V5, type="l", lwd=1, col=4))

legend("topleft", inset=c(0.05,0.05), y.intersp=1.5, cex=1, bty="n",
       legend=c("normal", "inverse"), lty=1, lwd=2, col=c(3,'black'))

# lines(tfr.tb.mnsc$session_num, tfr.tb.mnsc$sesc, type="l", lwd=2, col=1)
# with(tfr.tb.scCI, lines(session_num, V4, type="l", lwd=1, col=1))
# with(tfr.tb.scCI, lines(session_num, V5, type="l", lwd=1, col=1))


df.smr <- droplevels(subset(blk, patient%in%smr))
# All trials
df <- droplevels(subset(df.smr, reject_filter==0 & session_num>1))
smr.mnsc <- ddply(df, .(session_num), summarize, sesc=mean(adj_score))
smr.scCI <- ddply(df, .(session_num), bootin, varname="adj_score")
# Normal trials
df <- droplevels(subset(df.smr, reject_filter==0 & session_num>1 & norm_notInv==1 & trialtype==0))
nor.smr.mnsc <- ddply(df, .(session_num), summarize, sesc=mean(adj_score))
nor.smr.scCI <- ddply(df, .(session_num), bootin, varname="adj_score")
# Inverse trials
df <- droplevels(subset(df.smr, reject_filter==0 & session_num>22 & norm_notInv==0 & trialtype==0))
inv.smr.mnsc <- ddply(df, .(session_num), summarize, sesc=mean(adj_score))
inv.smr.scCI <- ddply(df, .(session_num), bootin, varname="adj_score")
# Transfer trials
df <- droplevels(subset(df.smr, reject_filter==0 & session_num>30 & norm_notInv==1 & trialtype>0))
tfr.smr.mnsc <- ddply(df, .(session_num), summarize, sesc=mean(adj_score))
tfr.smr.scCI <- ddply(df, .(session_num), bootin, varname="adj_score")

# master plot
xtk = seq(2,40)
plot(range(xtk),
     range(c(nor.smr.scCI$V4, nor.smr.scCI$V5, inv.smr.scCI$V4, inv.smr.scCI$V5)), 
     type="n", xaxt="n", xlab="all SMR sessions", ylab="adjusted score")
axis(1, at=xtk)
polygon(c(xtk, rev(xtk)), c(nor.smr.scCI$V5, rev(nor.smr.scCI$V4)), col=5, border=NA)
lines(nor.smr.mnsc$session_num, nor.smr.mnsc$sesc, type="l", lwd=2, col=3)
lines(inv.smr.mnsc$session_num, inv.smr.mnsc$sesc, type="l", lwd=2, col='black')

with(nor.smr.scCI, lines(session_num, V4, type="l", lwd=1, col=2))
with(nor.smr.scCI, lines(session_num, V5, type="l", lwd=1, col=2))
with(inv.smr.scCI, lines(session_num, V4, type="l", lwd=1, col=4))
with(inv.smr.scCI, lines(session_num, V5, type="l", lwd=1, col=4))

legend("topleft", inset=c(0.05,0.05), y.intersp=1.5, cex=1, bty="n",
       legend=c("normal", "inverse"), lty=1, lwd=2, col=c(3,'black'))

# lines(tfr.smr.mnsc$session_num, tfr.smr.mnsc$sesc, type="l", lwd=2, col=1)
# with(tfr.smr.scCI, lines(session_num, V4, type="l", lwd=1, col=1))
# with(tfr.smr.scCI, lines(session_num, V5, type="l", lwd=1, col=1))

# ALTOGETHER NOW - master plot
xtk = seq(2,40)
plot(range(xtk),
     range(c(nor.smr.scCI$V4, nor.smr.scCI$V5, inv.smr.scCI$V4, inv.smr.scCI$V5, 
             nor.tb.scCI$V4, nor.tb.scCI$V5, inv.tb.scCI$V4, inv.tb.scCI$V5)), 
     type="n", xaxt="n", xlab="all sessions", ylab="adjusted score")
axis(1, at=xtk)

polygon(c(xtk, rev(xtk)), c(nor.tb.scCI$V5, rev(nor.tb.scCI$V4)), col=2, border=NA)
polygon(c(xtk, rev(xtk)), c(nor.smr.scCI$V5, rev(nor.smr.scCI$V4)), col=2, border=NA)

lines(nor.tb.mnsc$session_num, nor.tb.mnsc$sesc, type="l", lwd=2, lty=2, col=3)
lines(nor.smr.mnsc$session_num, nor.smr.mnsc$sesc, type="l", lwd=2, col=3)

lines(inv.tb.mnsc$session_num, inv.tb.mnsc$sesc, type="l", lwd=2, lty=2, col='black')
lines(inv.smr.mnsc$session_num, inv.smr.mnsc$sesc, type="l", lwd=2, col='black')

with(inv.tb.scCI, lines(session_num, V4, type="l", lwd=1, lty=2, col='black'))
with(inv.tb.scCI, lines(session_num, V5, type="l", lwd=1, lty=2, col='black'))
with(inv.smr.scCI, lines(session_num, V4, type="l", lwd=1, col='black'))
with(inv.smr.scCI, lines(session_num, V5, type="l", lwd=1, col='black'))

legend("topleft", inset=c(0.05,0.0), y.intersp=1, cex=1, bty="n",
       legend=c("normal TB", "normal SMR", "inverse TB", "inverse SMR"), 
       lty=c(2,1,2,1), lwd=2, col=c(3,3,'black','black'))


# with(nor.smr.scCI, lines(session_num, V4, type="l", lwd=1, col=2))
# with(nor.smr.scCI, lines(session_num, V5, type="l", lwd=1, col=2))



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ----
## ## ## ## ## --------------------- trend of all LCs with CIs ----

# All trials
df <- droplevels(subset(blk, reject_filter==0 & session_num>1))
mnsc <- ddply(df, .(session_num), summarize, sesc=mean(adj_score))
scCI <- ddply(df, .(session_num), bootin, varname="adj_score")
# Normal trials
df <- droplevels(subset(blk, reject_filter==0 & session_num>1 & norm_notInv==1 & trialtype==0))
nor.mnsc <- ddply(df, .(session_num), summarize, sesc=mean(adj_score))
nor.scCI <- ddply(df, .(session_num), bootin, varname="adj_score")
# Inverse trials
df <- droplevels(subset(blk, reject_filter==0 & session_num>19 & norm_notInv==0 & trialtype==0))
inv.mnsc <- ddply(df, .(session_num), summarize, sesc=mean(adj_score))
inv.scCI <- ddply(df, .(session_num), bootin, varname="adj_score")
# Transfer trials
df <- droplevels(subset(blk, reject_filter==0 & session_num>30 & norm_notInv==1 & trialtype>0))
tfr.mnsc <- ddply(df, .(session_num), summarize, sesc=mean(adj_score))
tfr.scCI <- ddply(df, .(session_num), bootin, varname="adj_score")

# master plot
xtk = seq(2,40)
plot(range(xtk),
      range(c(nor.scCI$V4, nor.scCI$V5, inv.scCI$V4, inv.scCI$V5)), 
      type="n", xaxt="n", xlab="all sessions", ylab="adjusted score")
axis(1, at=xtk)
polygon(c(xtk, rev(xtk)), c(nor.scCI$V5, rev(nor.scCI$V4)), col=5, border=NA)
lines(nor.mnsc$session_num, nor.mnsc$sesc, type="l", lwd=2, col=3)
lines(inv.mnsc$session_num, inv.mnsc$sesc, type="l", lwd=2, col='black')

with(nor.scCI, lines(session_num, V4, type="l", lwd=1, col=2))
with(nor.scCI, lines(session_num, V5, type="l", lwd=1, col=2))
with(inv.scCI, lines(session_num, V4, type="l", lwd=1, col=4))
with(inv.scCI, lines(session_num, V5, type="l", lwd=1, col=4))

legend("topleft", inset=c(0.05,0.05), y.intersp=1.5, cex=1, bty="n",
       legend=c("normal", "inverse"), lty=1, lwd=2, col=c(3,'black'))

# lines(tfr.mnsc$session_num, tfr.mnsc$sesc, type="l", lwd=2, col=1)
# with(tfr.scCI, lines(session_num, V4, type="l", lwd=1, col=1))
# with(tfr.scCI, lines(session_num, V5, type="l", lwd=1, col=1))

# All trials
with(df, plot(range(session_num), range(c(scCI$V4, scCI$V5)), type="n", xaxt="n", xlab="all trials", ylab="adj_score"))
xtk = sort(unique(mnsc$session_num))
axis(1, at=xtk)
polygon(c(xtk, rev(xtk)), c(scCI$V5, rev(scCI$V4)), col=5, border=NA)
with(df, lines(xtk, mnsc$sesc, type="l", lwd=2, col=3))
with(scCI, lines(xtk, V4, type="l", lwd=1, col=2))
with(scCI, lines(xtk, V5, type="l", lwd=1, col=2))

# Normal trials
with(df, plot(range(session_num), range(c(nor.scCI$V4, nor.scCI$V5)), type="n", xaxt="n", xlab="normal trials", ylab="adj_score"))
xtk = sort(unique(nor.mnsc$session_num))
axis(1, at=xtk)
polygon(c(xtk, rev(xtk)), c(nor.scCI$V5, rev(nor.scCI$V4)), col=5, border=NA)
with(df, lines(xtk, nor.mnsc$sesc, type="l", lwd=2, col=3))
with(nor.scCI, lines(xtk, V4, type="l", lwd=1, col=2))
with(nor.scCI, lines(xtk, V5, type="l", lwd=1, col=2))

# Inverse trials
with(df, plot(range(session_num), range(c(inv.scCI$V4, inv.scCI$V5)), type="n", xaxt="n", xlab="inverse trials", ylab="adj_score"))
xtk = sort(unique(inv.mnsc$session_num))
axis(1, at=xtk)
polygon(c(xtk, rev(xtk)), c(inv.scCI$V5, rev(inv.scCI$V4)), col=5, border=NA)
with(df, lines(xtk, inv.mnsc$sesc, type="l", lwd=2, col=3))
with(inv.scCI, lines(xtk, V4, type="l", lwd=1, col=2))
with(inv.scCI, lines(xtk, V5, type="l", lwd=1, col=2))

# Transfer trials
with(df, plot(range(session_num), range(c(tfr.scCI$V4, tfr.scCI$V5)), type="n", xaxt="n", xlab="transfer trials", ylab="adj_score"))
xtk = sort(unique(tfr.mnsc$session_num))
axis(1, at=xtk)
polygon(c(xtk, rev(xtk)), c(tfr.scCI$V5, rev(tfr.scCI$V4)), col=5, border=NA)
with(df, lines(xtk, tfr.mnsc$sesc, type="l", lwd=2, col=3))
with(tfr.scCI, lines(xtk, V4, type="l", lwd=1, col=2))
with(tfr.scCI, lines(xtk, V5, type="l", lwd=1, col=2))
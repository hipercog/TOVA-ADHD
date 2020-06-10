# library(lattice)
# library(car)
# library(visreg)
# library(foreign)
# library(MASS)
# library(lme4)
# # library(nlme)
# library(lmerTest)
# library(gvlma)
# library(gee)
library(geepack)
library(QICpack)
library(moments)
library(colorspace)
library(ROCR)
# library(pROC)
# library(caTools)
library(boot)
library(lattice)
library(gplots)

chkldpkg(c("lattice","car","visreg","foreign","MASS","lme4","lmerTest","gvlma",
           "geepack","moments","plyr"))

## useful indices
ids <- substr(unique(blocs$ID), 8, 9)
rules <- unique(blocs$currentrule)
globs <- unique(blocs$featg)
locas <- unique(blocs$featl)


## ------------ generic subsetting per hypothesis ---------------------------
# get a temp data structure to hold variables of interest
tmpbk <- merge(blocs, cortimes, sort=FALSE)
tmpbk <- merge(tmpbk, search, sort=FALSE)

# from this we extract a subset for a hypothesis, used in all model fits

# H1a: performance during integration task blocks (global processing) will be
# unaffected by target stimulus condition, comparing faces, letters, colours.
{
hydat<-subset(hydat, currentrule=="G1" & featl=="shape")  # 1.1
hydat<-subset(hydat, currentrule=="G1" & featl=="letter") # 1.2
hydat<-subset(hydat, currentrule=="G1" & featl=="patch")  # 1.3a
hydat<-subset(hydat, currentrule=="G1" & featl!="patch")  # 1.3b
hydat<-subset(hydat, currentrule=="G1" )                  # 1.4
hydat<-subset(hydat, currentrule=="G2" & featl=="shape")  # 1.1.2
hydat<-subset(hydat, currentrule=="G2" & featl=="letter") # 1.2.2
hydat<-subset(hydat, currentrule=="G2" & featl=="patch")  # 1.3.2a
hydat<-subset(hydat, currentrule=="G2" & featl!="patch")  # 1.3.2b
hydat<-subset(hydat, currentrule=="G2" )                  # 1.4.2
hydat<-subset(hydat, focus=="G" & featl!="patch")         # 1.5
hydat<-subset(hydat, focus=="G")                          # 1.6
}
# simple diff test - with no conditioning on subject, use normalised times
{
hydat <- outgoners(hydat, hydat$discrim)
with(hydat, tapply(discrim, featg, shapiro.test))
wilcox.test(discrim ~ featg, data=hydat)
boxplot(discrim ~ featg, data=hydat, main="discrimination score")


hydat <- outgoners(hydat, hydat$mn)
with(hydat, tapply(mn, featg, shapiro.test))
wilcox.test(mn ~ featg, data=hydat)
boxplot(mn ~ featg, data=hydat, main="normalised mean msrt")


hydat <- outgoners(hydat, hydat$varn)
with(hydat, tapply(varn, featg, shapiro.test))
wilcox.test(varn ~ featg, data=hydat)
boxplot(varn ~ featg, data=hydat, main="normalised msrtv")
}

# H1b: performance during focal attention blocks (local processing) will be 
# unaffected by target stim condition, comparing shapes, letters, orientation
{
hydat<-subset(hydat, currentrule=="L1" & featg=="face")   # 2.1
hydat<-subset(hydat, currentrule=="L1" & featg=="letter") # 2.2
hydat<-subset(hydat, currentrule=="L1" & featg=="noise")  # 2.3a
hydat<-subset(hydat, currentrule=="L1" & featg!="noise")  # 2.3b
hydat<-subset(hydat, currentrule=="L1" )                  # 2.4
hydat<-subset(hydat, currentrule=="L2" & featg=="face")   # 2.1.2
hydat<-subset(hydat, currentrule=="L2" & featg=="letter") # 2.2.2
hydat<-subset(hydat, currentrule=="L2" & featg=="noise")  # 2.3.2a
hydat<-subset(hydat, currentrule=="L2" & featg!="noise")  # 2.3.2b
hydat<-subset(hydat, currentrule=="L2" )                  # 2.4.2
hydat<-subset(hydat, focus=="L" & featg!="noise")         # 2.5
hydat<-subset(hydat, focus=="L")                          # 2.6
}
# simple diff test - with no conditioning on subject, use normalised times
{
hydat <- outgoners(hydat, hydat$discrim)
wilcox.test(discrim ~ featl, data=hydat)
boxplot(discrim ~ featl, data=hydat, main="Rule=L2", ylab="discrim score")


hydat <- outgoners(hydat, hydat$mn)
wilcox.test(mn ~ featl, data=hydat)


hydat <- outgoners(hydat, hydat$varn)
wilcox.test(varn ~ featl, data=hydat)
}

hydat<-droplevels(subset(tmpbk, featg!="noise" & featl!="patch"))     # 3

## ---- PREP FOR EXPLORATION AND FITTING LINEAR MODELS ----

# trim 5% of outliers
hydat <- outgoners(hydat, hydat$mt, 1.96)
hydat <- outgoners(hydat, hydat$vart, 1.96)
hydat <- outgoners(hydat, hydat$discrim, 1.96)
# ...for search response times
hydat <- outgoners(hydat, hydat$srt, 1.96)
hydat <- outgoners(hydat, hydat$srtv, 1.96)

# when blocks will be modelled, define a linearly incrementing block count ----
hydat$bloco <- NA
# tapply works for a single level of partitioning
temp <- with(hydat, tapply(block, subj, rank))
hydat$bloco <- unlist(temp)
rm(temp)
# ddply for more complex partitioning
hydat <- ddply(hydat, .(set, subj), blocord)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ----
## ## ## ## ## --------------------- BWRC COLOURS --------------------

# colors: 1 orange - 2 blue - 3 dark blue - 4 darker gray - 5 lighter gray - 6 green - 7 red
bwrc <- palette(c('#F0AA00','#00A0BE','#003C78','#828282','#F0F0F0','#B4BE00','#FA4100'))
bwrcRGB <- hex2RGB(c('#F0AA00','#00A0BE','#003C78','#828282','#F0F0F0','#B4BE00','#FA4100'))
bwrcLAB <- as(bwrcRGB, "LAB")
bwrcHSV <- as(bwrcRGB, "HSV")
bwrcDst <- matrix(NaN, nrow=length(bwrc), ncol=length(bwrc))
# rownames(bwrcDst) <- ...
# colnames(bwrcDst) <- ...
for (i in 1:6){
  for (j in (i+1):7){
    L1 <- bwrcLAB@coords[i,"L"]
    L2 <- bwrcLAB@coords[j,"L"]
    a1 <- bwrcLAB@coords[i,"A"]
    a2 <- bwrcLAB@coords[j,"A"]
    b1 <- bwrcLAB@coords[i,"B"]
    b2 <- bwrcLAB@coords[j,"B"]
    bwrcDst[i,j] <- (L1-L2)^2 + (a1-a2)^2 + (b1-b2)^2
  }
}
bwrcol <- colorRampPalette(c('#828282','#F0F0F0','#00A0BE','#003C78'))
bwrcol <- colorRampPalette(c('#F0AA00','#00A0BE','#003C78'))
bwrcol <- colorRampPalette(c('#00A0BE','#F0AA00'))

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ----
## ## ## ## ## --------------------- EXPLORATION --------------------
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
hist(hydat$mt)
hist(sqrt(hydat$mt))
hist(hydat$vart)
hist(sqrt(hydat$vart))
hist(hydat$discrim)
hist(sqrt(hydat$discrim))
hist(hydat$discrim^2)

spacer <- c(1,2,3,4,6,7,8,9,11,12,13,14,16,17,18,19)

tmp <- outgoners(hydat, hydat$mt, 1.96)
with(tmp, boxplot(mt ~ featg*featl*currentrule, col=c(1,2,3,5), xaxt="n",
                  main="4 stimulus conditions x Rule", ylab="RT (sec)", ylim=c(0,2),
                  notch=TRUE, pch="-",
                  at=spacer))
axis(1, at=spacer, labels = FALSE)

mars <- par("mar")
mars[1] <- 5.1
par(mar=mars)

tmp <- outgoners(hydat, hydat$vart, 1.96)
with(tmp, boxplot(vart ~ featg*featl*currentrule, col=c(1,2,3,5), xaxt="n", 
                  pch="-", ylab="RTV (sec)", ylim=c(0,1.5),
                  notch=TRUE,
                  at=spacer))

axis(1, at=spacer, labels = FALSE)

lbls <- sort(unique(tmp$currentrule))
text(c(3,8,13,18), par("usr")[3] - 0.4, adj = 1, labels = lbls, xpd = TRUE, cex=1.1)#seq(3,15,4)

## LABELS WITH COMPLETE NAMES
lbls <- c("face", "letter", "face", "letter", "face", "letter", "face", "letter", 
          "face", "letter", "face", "letter", "face", "letter", "face", "letter")
text(spacer, par("usr")[3] - 0.1, srt = 45, adj = 0.5, labels = lbls, xpd = TRUE, cex=0.8)
lbls <- c("shape", "letter", "shape", "letter", "shape", "letter", "shape", "letter")
text(c(2,4,7,9,12,14,17,19), par("usr")[3] - 0.25, adj = 0.8, labels = lbls, xpd = TRUE)#seq(2,16,2)

## LABELS WITH CAPITALS ONLY
lbls <- c("F", "L", "F", "L", "F", "L", "F", "L", "F", "L", "F", "L", "F", "L", "F", "L")
text(1:16, par("usr")[3] - 0.1, adj = 0.5, labels = lbls, xpd = TRUE, cex=0.8)
lbls <- c("s", "lr", "s", "lr", "s", "lr", "s", "lr")
text(seq(2,16,2), par("usr")[3] - 0.2, adj = 2, labels = lbls, xpd = TRUE)


# tmp <- outgoners(hydat, hydat$discrim, 1.96)
# with(hydat, boxplot(discrim ~ featg*featl*currentrule, col=c(1,2,3,5),
#                     xaxt="n", pch="-",  main="4 stimulus conditions x Rule", 
#                     ylab="d (discrimination score)", ylim=c(0,1)))

rm(tmp)

# ---- # ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ----
## -- SUBJECT META-DATA (e.g. DIFFICULTY) VERSUS BLOCK-LEVEL DATA ----
df <- merge(hydat, meta, by="ID")


# ---- # ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ----
## --------------- RESPONSE TIME SWITCH COST ----
# (label block where previous block had different focus)
hydat$switsh <- logical(nrow(hydat))
for(i in 2:nrow(hydat)){
  if (hydat[i-1,]$focus!=hydat[i,]$focus){
    hydat$switsh[i] <- TRUE
  }
}
with(hydat, boxplot(mt ~ switsh, ylim=c(0,2.5), ylab="RT (sec)"))
with(hydat, boxplot(vart ~ switsh, ylim=c(0,1), ylab="RTV (sec)"))
with(hydat, boxplot(discrim ~ switsh, ylim=c(0,1), ylab="discrimination"))
# ------- General Estimating Equations
gee_switch_i <- geeglm(vart ~ 1 + block + switsh, 
                       data=hydat, id=ID, family=Gamma, corstr="exchangeable" )
summary(model <- gee_switch_i)
qic(model)
# ------- ROC curve
# ROCR approach
pred <- prediction(hydat$discrim, hydat$switsh)
perf <- performance(pred, "tpr", "fpr")
plot(perf, lwd=2, main="switch cost block  :  change level v. no change", 
     xlab=paste("fpr /  no change discrimination"),
     ylab=paste("tpr / change level discrimination"))

auc <- performance(pred, "auc")
auc <- auc@y.values[[1]]

A <- round(auc,2)
P <- round(aucp(auc, 1204, 1196), 2)
E <- round(aucES(auc, 1204, 1196), 2)
legend("bottomright", inset=c(0.05,0.05), y.intersp=2, cex=1.3, bty="n",
       legend=bquote(atop(atop("AUC = "~.(A),
                               italic(p)~" = "~.(P)),
                          atop("ES = "~.(E), ""))))


# ---- # ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ----
# -------- SEARCH ERRORS GLOBAL VERSUS LOCAL ----

# matrix plot of errors?
Lrt <- hydat[hydat$focus=="L",]
Lrt <- Lrt[with(Lrt, order(omerr)), ]
Grt <- hydat[hydat$focus=="G",]
Grt <- Grt[with(Grt, order(omerr)), ]
# plot(Lrt$omerr, Grt$omerr)
# abline(0,1)
symx <- c(table(Lrt$omerr))#[1:15]
symy <- c(table(Grt$omerr))#[1:15]
sym <- c(table(Grt$omerr, Lrt$omerr))
sysq <- matrix(sym, nrow = 16, byrow = TRUE)
sysq=sysq[1:7,1:7]
sysq[sysq==0] <- -Inf
# errors per block count
levelplot(sysq, xlab="local", ylab="global", col.regions=bwrcol(300), cex.axis=3)

# ---- # ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ----
## -------- SEARCH TIME PERFORMANCE VERSUS RULE FOUND PERFORMANCE ----

# ...for rule-found response times
hydat <- outgoners(hydat, hydat$rfrt, 1.96)
hydat <- outgoners(hydat, hydat$rfrtv, 1.96)
hydat <- hydat[hydat$omerr!=0 & hydat$omerr<8,]

# exploration plots
with(hydat, plot(srt, mt))
with(hydat[!is.na(hydat$srt),], cor(mt, srt))
with(hydat, plot(srt, rfrt))
with(hydat[!is.na(hydat$srt),], cor(rfrt, srt))

with(hydat, plot(srtv, vart, col=currentrule))
with(hydat[!is.na(hydat$srtv),], cor(vart, srtv))
with(hydat, plot(srtv, rfrtv))
with(hydat[!is.na(hydat$srt),], cor(rfrtv, srt))

with(hydat, boxplot(rfrt ~ omerr, add=TRUE))
with(hydat[!is.na(hydat$rfrt),], cor(rfrt, omerr))
summary(with(hydat[!is.na(hydat$rfrt),], lm(rfrt ~ omerr)))

with(hydat[!is.na(hydat$omerr),], cor(comer, omerr))

# mean RT per omerr level, with bootstrapped 95% CIs
ormnrt <- with(hydat, tapply(rfrt, omerr, mean))
bCI <- ddply(hydat, .(omerr), bootin, varname="rfrt")
xidx <- sort(unique(hydat$omerr))

with(hydat, plot(range(omerr), range(c(bCI$V4, bCI$V5)), type="n", xaxt="n", 
                 xlab="",
                 ylab="'rule found' mean RT (sec)"))
axis(1, at=c(1:9))
polygon(c(xidx, rev(xidx)), c(bCI$V5, rev(bCI$V4)), col=5, border=NA)
with(hydat, lines(xidx, ormnrt, type="l", lwd=2, col=3))
with(bCI, lines(xidx, V4, type="l", lwd=1, col=2))
with(bCI, lines(xidx, V5, type="l", lwd=1, col=2))


# colors: 1 orange - 2 blue - 3 dark blue - 4 darker gray - 5 lighter gray - 6 green - 7 red
# # mean RT per attention focus (G=global L=local), bootstrap 95% CIs
Gormnrt <- with(hydat[hydat$focus=="G",], tapply(rfrt, omerr, mean))
Lormnrt <- with(hydat[hydat$focus=="L",], tapply(rfrt, omerr, mean))
bCIg <- ddply(hydat[hydat$focus=="G",], .(omerr), bootin, varname="rfrt")
bCIl <- ddply(hydat[hydat$focus=="L",], .(omerr), bootin, varname="rfrt")

# with(hydat, interaction.plot(omerr, focus, rfrt, col=3, lwd=2, ylim=c(0.7,1.8),
#                              xlab="'search phase' error count", xaxt="n",
#                              ylab="'rule found' mean RT (sec)"))

with(hydat, plot(range(omerr), range(c(bCIl$V5, bCIg$V4)), type="n", xaxt="n", 
                 xlab="'search phase' error count",
                 ylab="'rule found' mean RT (sec)"))
axis(1, at=c(1:9))
polygon(c(xidx, rev(xidx)), c(bCIg$V5, rev(bCIg$V4)), col=5, border=NA)

with(hydat, lines(xidx, Gormnrt, type="l", lwd=2, col=1))
with(hydat, lines(xidx, Lormnrt, type="l", lwd=2, col='black', lty=2))

with(bCIg, lines(xidx, V4, type="l", lwd=1, col=1))
with(bCIg, lines(xidx, V5, type="l", lwd=1, col=1))
with(bCIl, lines(xidx, V4, type="l", lwd=1, col='black', lty=2))
with(bCIl, lines(xidx, V5, type="l", lwd=1, col='black', lty=2))

legend("topleft", inset=c(0.05,0.05), y.intersp=1.5, cex=1, bty="n",
       legend=c("global", "local"), lty=c(1,2), lwd=3, col=c(1,'black'))


# ------- General Estimating Equations
gee_search_i <- geeglm(rfrt ~ 1 + block + srt + omerr, 
                       data=hydat[!is.nan(hydat$srt),], 
                       id=ID, family=Gamma, corstr="exchangeable" )
summary(model <- gee_search_i)
qic(model)


# ---- # ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ----
## -------- BLOCK PERFORMANCE VERSUS ALL-TEST PERFORMANCE ----
df <- outgoners(cortrials, cortrials$t_reaction)
df <- droplevels(subset(df, featg!="noise" & featl!="patch"))


bkxrt <- with(df, tapply(t_reaction, block, mean))
stxbkrt <- ddply(df, .(set, block), summarize, rt=mean(t_reaction))
stxbkrt$bpset <- seq(1, 12)

# HEATMAP OF BLOCK VERSUS SET MEAN TIMES
sxbmat <- stxbkrt[,c(1,3,4)]
sxbmat1 <- reshape(sxbmat, idvar="set", timevar="bpset", direction="wide")
sxbmat2 <- data.matrix(sxbmat1[,2:13])
rownames(sxbmat2) <- sxbmat1$set

# following code limits the lowest and highest color to 5%, and 95% of range
# quantile.range <- quantile(sxbmat2, probs = seq(0, 1, 0.01))
# palette.breaks <- seq(quantile.range["1%"], quantile.range["99%"], 0.1)
heatmap.2(
  sxbmat2,
  Rowv=NA, Colv=NA,
  dendrogram = "none",
  scale = "none",
  trace = "none",
  col = bwrcol,
  #   key        = FALSE,
  keysize = 3,
  key.title = "",
  density.info=c("density"),
  denscol="black",
  key.xlab = "RT (sec)",
  key.ylab = "",
  key.ytick = "",
  #   labRow     = NA,
  labCol = 1:12,
  srtCol = 0,
  xlab = "block / set",
  ylab = "full-stimulus set"
)
# UNUSED SHITTY INTERACTION PLOTS
# ckCI <- ddply(df, .(block), bootin, varname="t_reaction")
# with(df, plot(range(block), range(c(0,2)), type="n", xaxt="n", 
#                  xlab="block", ylab="mean RT (sec)"))
# axis(1, at=c(1:120))
# with(df, lines(1:120, bkxrt, type="l", lwd=2, col=3))
# with(ckCI, lines(1:120, V4, type="l", lwd=1, col=2))
# with(ckCI, lines(1:120, V5, type="l", lwd=1, col=2))
# with(stxbkrt, interaction.plot(bpset, set, rt, col=3, lwd=2,
#               xlab="block", ylab="mean RT (sec)"))

# WITHIN BLOCK TREND, ALL TRIALS ---
df <- droplevels(subset(bedat$logs, featg!="noise" & featl!="patch"))
df <- outgoners(df, df$trial, 1.96)
df <- outgoners(df, df$t_reaction, 1.96)
df <- ddply(df, .(ID, block), trlord)
hist(df$trial)
trl <- ddply(df, .(trial), summarize, rt=mean(t_reaction))
trCI <- ddply(df, .(trial), bootin, varname="t_reaction")
with(df, plot(range(trial), range(c(trCI$V4, trCI$V5)), type="n", xaxt="n", 
              xlab="trial", ylab="mean RT (sec)"))
xtk = c(1:length(unique(df$trial)))
axis(1, at=xtk)
polygon(c(xtk, rev(xtk)), c(trCI$V5, rev(trCI$V4)), col=5, border=NA)
with(df, lines(xtk, trl$rt, type="l", lwd=2, col=3))
with(trCI, lines(xtk, V4, type="l", lwd=1, col=2))
with(trCI, lines(xtk, V5, type="l", lwd=1, col=2))

# WITHIN BLOCK TREND, CORRECT TRIALS ---
df <- droplevels(subset(cortrials, featg!="noise" & featl!="patch"))
df <- outgoners(df, df$trial, 1.96)
df <- outgoners(df, df$t_reaction, 1.96)
df <- ddply(df, .(ID, block), trlord)
hist(df$trial)
trl <- ddply(df, .(trial), summarize, rt=mean(t_reaction))
trCI <- ddply(df, .(trial), bootin, varname="t_reaction")
with(df, plot(range(trial), range(c(trCI$V4, trCI$V5)), type="n", xaxt="n", 
              xlab="correct trials", ylab="mean RT (sec)"))
xtk = c(1:length(unique(df$trial)))
axis(1, at=xtk)
polygon(c(xtk, rev(xtk)), c(trCI$V5, rev(trCI$V4)), col=5, border=NA)
with(df, lines(xtk, trl$rt, type="l", lwd=2, col=3))
with(trCI, lines(xtk, V4, type="l", lwd=1, col=2))
with(trCI, lines(xtk, V5, type="l", lwd=1, col=2))



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ----
## ## ## ## ## --------------------- REGRESSION ANALYSIS --------------------
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# A four-factor contrast matrix
contr.mat <- matrix(c(c(-1,1,0,0),c(0,-1,1,0),c(0,0,-1,1)),4)
# A four-factor sliding-differences (repeated contrasts) matrix
contr.mat <- matrix(c(c(-0.75,0.25,0.25,0.25),
                      c(-0.5,-0.5,0.5,0.5),
                      c(-0.25,-0.25,-0.25,0.75)),4)
# Contrast names
colnames(contr.mat) <- c("G1","G2","L2")
contrasts(cortimes$currentrule) <- contr.mat

# GEE significance: computed from output of summary(gee(...)) which prints
# either a "z" or a "t" depending in the "family" option. z follows, under
# the null hypothesis, a normal distribution N(0, 1), you have the
# corresponding P-value with (for a two-tailed test):
#   
#   2 * (1 - pnorm(abs(z)))
# 
# t follows a 'Student' distribution with df degrees of freedom given by N- k
# - 1, where N is the number of observations, and k is the number of
# estimated paramaters. I think, but am not definitely sure, that N is
# counted among all clusters, and k is the number of parameters in the GLM
# eventually included the estimated scale (correlation parameters are not
# counted). As above, you have the P-value with:
#   
#   2 * (1 - pnorm(abs(t), df))

# ---- # ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ----
# ---- TEST DISCRIMINATION ----
# shapiro test normality
shapiro.test(hydat$discrim)
shapiro.test(sqrt(hydat$discrim))

# --- ANOVA ----
aov_d <- aov(discrim ~ featg*bloco + Error(ID), hydat)

# Linear Mixed model fitting
lmer_H1G_d <- lmer(discrim ~ 0 + featg*featl + bloco + (1|ID), hydat, FALSE)
lmer_d_null <- lmer(discrim ~ 0 + (1|ID), hydat, FALSE)
model <- lmer_H1G_d
nullmdl <- lmer_d_null

# ------- General Estimating Equations -------
## model that works for testing stimulus condition effects per rule
gee_d_feat_i <- geeglm(discrim ~ 1 + bloco + featg*featl, data=hydat, id=subj,
                       subset=currentrule=="L2", family=Gamma, corstr="exchangeable" )
summary(model <- gee_d_feat_i)
qic(model)
## model that works for testing stimulus condition effects per rule
gee_d_rule_i <- geeglm(discrim ~ 1 + bloco + currentrule, data=hydat, id=subj,
                       family=Gamma, corstr="exchangeable" )
summary(model <- gee_d_rule_i)
qic(model)
## model to test difference between focal levels
gee_d_focus_i <- geeglm(discrim ~ 1 + bloco + focus, data=hydat, id=subj,
                        family=Gamma, corstr="exchangeable" )
summary(model <- gee_d_focus_i)
qic(model)
## models that works for testing stimulus condition effects per rule
tmp <- droplevels(subset(hydat, focus=="G"))
gee_d_ruleG_i <- geeglm(discrim ~ 1 + bloco + currentrule, data=tmp, id=subj,
                        family=Gamma, corstr="exchangeable" )
summary(model <- gee_d_ruleG_i)
qic(model)
tmp <- droplevels(subset(hydat, focus=="L"))
gee_d_ruleL_i <- geeglm(discrim ~ 1 + bloco + currentrule, data=tmp, id=subj,
                        family=Gamma, corstr="exchangeable" )
summary(model <- gee_d_ruleL_i)
qic(model)
rm(tmp)

# ---- # ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ----
# ---- TEST RESPONSE TIME ----
# shapiro test normality
shapiro.test(hydat$mt)
shapiro.test(sqrt(hydat$mt))

# --- ANOVA ----
aov.mt.G <- aov(mt ~ featg + Error(ID), hydat)
aov.mt.L <- aov(mt ~ featl + Error(ID), hydat)
aov.mt.GxL <- aov(mt ~ featg*featl + Error(ID), hydat)
aov.mt.GxLxB <- aov(mt ~ featg*featl + bloco + Error(ID), hydat)
summary(aov.mt.G)
summary(aov.mt.L)
summary(aov.mt.GxL)
summary(aov.mt.GxLxB)

with(hydat, interaction.plot(bloco, featg, mt, ylab="rt ms", main="full-stim only"))
with(hydat[hydat$featg=="letter",], abline(lm(mt ~ bloco), lty=2))
with(hydat[hydat$featg=="face",], abline(lm(mt ~ bloco), lty=3))

# ----- Linear Mixed model fitting ------
lmer_H1G_d <- lmer(discrim ~ 0 + featg*featl + bloco + (1|ID), hydat, FALSE)
lmer_d_null <- lmer(discrim ~ 0 + (1|ID), hydat, FALSE)

# ------- General Estimating Equations -------
## model that works for testing stimulus condition effects per rule
gee_mt_feat_i <- geeglm(mt ~ 1 + bloco + featg*featl, data=hydat, id=subj,
                        subset=currentrule=="G1", family=Gamma, corstr="exchangeable" )
summary(model <- gee_mt_feat_i)

## model that works for testing stimulus condition effects per rule
gee_mt_rule_i <- geeglm(mt ~ 1 + bloco + currentrule, data=hydat, id=subj,
                        family=Gamma, corstr="exchangeable" )
summary(model <- gee_mt_rule_i)

## model to test difference between focal levels
gee_mt_focus_i <- geeglm(mt ~ 1 + bloco + focus, data=hydat, id=subj,
                         family=Gamma, corstr="exchangeable" )
summary(model <- gee_mt_focus_i)

## models that works for testing stimulus condition effects per rule
tmp <- droplevels(subset(hydat, focus=="G"))
gee_mt_ruleG_i <- geeglm(mt ~ 1 + bloco + currentrule, data=tmp, id=subj,
                         family=Gamma, corstr="exchangeable" )
summary(model <- gee_mt_ruleG_i)
qic(model)
tmp <- droplevels(subset(hydat, focus=="L"))
gee_mt_ruleL_i <- geeglm(mt ~ 1 + bloco + currentrule, data=tmp, id=subj,
                         family=Gamma, corstr="exchangeable" )
summary(model <- gee_mt_ruleL_i)
qic(model)
rm(tmp)


## null model - needed?
gee_mt_zi_null <- geeglm(mt ~ 1, data=hydat, id=subj, 
                         waves=bloco, family=Gamma, corstr="exchangeable" )
nullmdl <- gee_mt_zi_null

# ---- # ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ----
# ---- TEST RESPONSE TIME VARIABILITY ----
# shapiro test normality
shapiro.test(hydat$vart)
shapiro.test(sqrt(hydat$vart))

# ----- Linear Mixed model fitting -----
lme_rul_rtv <- lmer(vart ~ 0 + currentrule + (1|ID), hydat, FALSE,
                    contrasts=contrasts(cortimes$currentrule))
model <- lme_rul_rtv
# all conditions predict response time variability
lme_cnd_rtv <- lmer(vart ~ 0 + currentrule + featg + featl + (1|ID), hydat,
                    FALSE, contrasts=contrasts(cortimes$currentrule) )
model <- lme_cnd_rtv

lmer_rtv <- lmer(varn ~ 0 + targetstim + (1|bloco) + (1|ID), hydat, FALSE)
lmer_rtvFB <- lmer(varn ~ 0 + targetstim + bloco + (1|ID), hydat, FALSE)
lmer_rtvNB <- lmer(varn ~ 0 + targetstim + (1|ID), hydat, FALSE)
lmer_rtv_null <- lmer(varn ~ 0 + (1|bloco) + (1|ID), hydat, FALSE)
lmer_rtv_nullB <- lmer(varn ~ 0 + (1|ID), hydat, FALSE)
# lme_H1G1_rtv <- lme(vart ~ 0 + featg + block, hydat, random=~I(block-12)|ID)
# model <- lme_H1G1_rtv

# ------- General Estimating Equations -------
## model that works for testing stimulus condition effects per rule
gee_rtv_feat_i <- geeglm(vart ~ 1 + bloco + featg*featl, data=hydat, id=subj,
                         subset=currentrule=="G1", family=Gamma, corstr="exchangeable" )
summary(model <- gee_rtv_feat_i)
qic(model)
## model that works for testing stimulus condition effects per rule
gee_rtv_rule_i <- geeglm(vart ~ 1 + bloco + currentrule, data=hydat, id=subj,
                         family=Gamma, corstr="exchangeable" )
summary(model <- gee_rtv_rule_i)
qic(model)
## model to test difference between focal levels
gee_rtv_focus_i <- geeglm(vart ~ 1 + bloco + focus, data=hydat, id=subj,
                          family=Gamma, corstr="exchangeable" )
summary(model <- gee_rtv_focus_i)
qic(model)
## models that works for testing stimulus condition effects per rule
tmp <- droplevels(subset(hydat, focus=="G"))
gee_rtv_ruleG_i <- geeglm(vart ~ 1 + bloco + currentrule, data=tmp, id=subj,
                          family=Gamma, corstr="exchangeable" )
summary(model <- gee_rtv_ruleG_i)
qic(model)
tmp <- droplevels(subset(hydat, focus=="L"))
gee_rtv_ruleL_i <- geeglm(vart ~ 1 + bloco + currentrule, data=tmp, id=subj,
                          family=Gamma, corstr="exchangeable" )
summary(model <- gee_rtv_ruleL_i)
qic(model)
rm(tmp)


# ---- # ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ----
# ---- TEST RIGHT VS LEFT PRESENTATION OF TARGET-MATCH STIMULI ----
df <- bedat$logs[bedat$logs$correct==1,]
df <- droplevels(df[df$target_pos=="right" | df$target_pos=="left",])

# ---- visuals
with(df, boxplot(t_reaction ~ target_pos*focus, 
                 main="LvR presented target", ylab="RT (sec)"))
with(df, boxplot(t_reaction ~ target_pos*currentrule, 
                 main="LvR presented target", ylab="RT (sec)"))
with(df, boxplot(t_reaction ~ target_pos*featg*featl, 
                 main="LvR presented target", ylab="RT (sec)"))
with(df, boxplot(t_reaction ~ target_pos*featg*featl*focus, 
                 main="LvR presented target", ylab="RT (sec)"))
with(df, boxplot(t_reaction ~ target_pos*featg*featl*currentrule, 
                 main="LvR presented target", ylab="RT (sec)"))
# --- ANOVA
aov.mt.RvL <- aov(t_reaction ~ target_pos*featg*featl*currentrule + Error(ID), df)
summary(aov.mt.RvL)
# ------- General Estimating Equations
gee_rt_RvL_i <- geeglm(t_reaction ~ 1 + block + target_pos*featg + target_pos*featl + target_pos*currentrule, 
                       data=df, id=ID, family=Gamma, corstr="exchangeable" )
summary(model <- gee_rt_RvL_i)
qic(model)


# ---- # ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ----
# ------- Checking and summarising any model -------------
# lm.goodness(model, hydat)

# gvmodel <- gvlma(model)
# summary(gvmodel)

visreg(model)
anova(nullmdl, model)

ggCaterpillar(ranef(model, condVar = T), QQ = T, likeDotplot = F)

qqnorm(resid(model), ylab="Residuals")
qqline(resid(model))


# ---- # ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ----
## ------------------------ ROC ANALYSIS ------------------------


# colAUC(su_blocs$discrim, class, plotROC=T, alg=c("Wilcoxon","ROC"))

## -------- ENVIRONMENTS FOR ROC CURVES BY SUBSET ---------
su_pred <- new.env()
su_perf <- new.env()

## --------------- TEST ROC CURVES ALL DATA ---------------
aucs_errs <- c(1:nerr)
rox_errs <- new.env()

for (e in 1:nerr){
  name <- paste0(diff[e], '_', e+2, '_consec_correct')
  
  ## ROCR approach
  pred <- prediction(blocs$discrim, blocs$cterr[,e])
  perf <- performance(pred, "tpr", "fpr")
  auc <- performance(pred, "auc")
  aucs_errs[e] <- auc@y.values[[1]]
  plot(perf, main=paste(name, paste0('auc=',round(auc@y.values[[1]],2))))
  
  #     assign(name, pred, envir=su_pred)
  #     assign(name, perf, envir=su_perf)
  
  ## pROC approach
  rox <- roc(blocs$cterr[,e], blocs$discrim, ci=TRUE)
  assign(name, rox, envir=rox_errs)
}



## --------------- ROC CURVES BY RULExSTIMULUS CONDITION ----
# rox_rule <- new.env()
aucs_hyps <- c(1:7)
hyps <- c("H1a", "H1b", "H2a", "H2b", "H2c", "H2d", "H3c")
names(aucs_hyps) <- hyps
aucs_subj <- matrix(NaN, nrow=length(ids), ncol=length(hyps))
dimnames(aucs_subj) <- list(ids, hyps)

r <- c("G1", "L1", "G2", "L2", "G2", "L2", "L1")
f <- c("G", "L", "G", "L", "L", "G", "L")
feat <- c("letter", "letter", "face", "shape", "shape", "face")

# ROC hypothesis testing data sets ----
hy.ps <- c(1:7)
names(hy.ps) <- hyps
hy.es <- c(1:7)
names(hy.es) <- hyps
# for (i in 1:6){
#   hy.ps[i] <- aucp(aucs_hyps[i], 300, 300)
#   hy.es[i] <- aucES(aucs_hyps[i], 300, 300)
# }
# hy.ps[7] <- aucp(aucs_hyps[7], 600, 600)
# hy.es[7] <- aucES(aucs_hyps[7], 600, 600)
# 
# subj.Hy.p <- matrix(NaN, nrow=length(ids), ncol=length(hyps))
# rownames(subj.Hy.p) <- ids
# colnames(subj.Hy.p) <- hyps
# subj.Hy.ES <- matrix(NaN, nrow=length(ids), ncol=length(hyps))
# rownames(subj.Hy.ES) <- ids
# colnames(subj.Hy.ES) <- hyps
# sz <- 12
# for (j in hyps){
#   if (j==hyps[7]){ sz <- 24 }
#   for (i in ids){ 
#     subj.Hy.p[i,j] <- aucp(aucs_subj[i,j], sz, sz)
#     subj.Hy.ES[i,j] <- aucES(aucs_subj[i,j], sz, sz)
#   }
# }

h <- 6
# one rule at a time, compare stimulus conditions ----
if (h<7){
  su_blocs <- subset(blocs, currentrule==r[h] & featg!="noise" & featl!="patch")
  if (f[h]=="G"){
    class <- su_blocs$featg==feat[h]
    nf1 <- unique(su_blocs[class,]$featg)
    nf2 <- unique(su_blocs[!class,]$featg)
    name <- paste(hyps[h], "  :   ", r[h], nf1, "v", nf2)
  }else{
    class <- su_blocs$featl==feat[h]
    nf1 <- unique(su_blocs[class,]$featl)
    nf2 <- unique(su_blocs[!class,]$featl)
    name <- paste(hyps[h], "  :   ", r[h], nf1, "v", nf2)
  }
}else{
  # ---- two rules at a time,compare rules ----
  su_blocs <- subset(blocs, focus==f[h] & featg!="noise" & featl!="patch")
  class <- su_blocs$currentrule==r[h]
  nf1 <- unique(su_blocs[class,]$currentrule)
  nf2 <- unique(su_blocs[!class,]$currentrule)
  name <- paste(hyps[h], nf1, "v", nf2)
}
# Make the main ROC curve ----
if (length(unique(class))==2){    
  # ROCR approach
  predm <- prediction(su_blocs$discrim, class)
  perfm <- performance(predm, "tpr", "fpr")
  plot(perfm, main=name, lwd=4, xlab=paste("fpr / ", nf2, "discrimination"),
       ylab=paste("tpr / ", nf1,"discrimination"))
  
  auc <- performance(predm, "auc")
  aucs_hyps[hyps[h]] <- auc@y.values[[1]]
  hy.ps[h] <- aucp(aucs_hyps[h], 300, 300)
  hy.es[h] <- aucES(aucs_hyps[h], 300, 300)
  
  A <- round(aucs_hyps[hyps[h]],2)
  P <- round(hy.ps[h], 2)
  E <- round(hy.es[h], 2)
  if(hy.ps[h] >= 0.001){
    legend("bottomright", inset=c(0.05,0.05), y.intersp=2, cex=1.3, bty="n",
           legend=bquote(atop(atop("AUC = "~.(A),
                                   italic(p)~" = "~.(P)),
                              atop("ES = "~.(E), ""))))
  }else{
    legend("bottomright", inset=c(0.05,0.05), y.intersp=2, cex=1.3, bty="n",
           legend=bquote(atop(atop("AUC = "~.(A),
                                   italic(p)~" < 0.001"),
                              atop("ES = "~.(E), ""))))  
  }
  
  assign(name, predm, envir=su_pred)
  assign(name, perfm, envir=su_perf)
}
## --------------- ROC CURVES BY SUBJECT ------------------
for (i in ids){
  if (h<7){
    su_blocs <- subset(blocs, currentrule==r[h] & featg!="noise" & 
                         featl!="patch" & ID==paste0('reknow0',i))
    if (f[h]=="G"){
      class <- su_blocs$featg==feat[h]
    }else{
      class <- su_blocs$featl==feat[h]
    }
  }else{
    su_blocs <- subset(blocs, focus==f[h] & featg!="noise" & 
                         featl!="patch" & ID==paste0('reknow0',i))
    class <- su_blocs$currentrule==r[h]
  }
  if (length(unique(class))==2){
    # ROCR approach
    pred <- prediction(su_blocs$discrim, class)
    perf <- performance(pred, "tpr", "fpr")
    plot(perf, add=TRUE, lty=2, col=1, lwd=1)
    auc <- performance(pred, "auc")
    aucs_subj[i,hyps[h]] <- auc@y.values[[1]]
  }
}
# -- Clean up main curve and overlay diagonal
abline(0,1,lwd=2)
plot(perfm, add=TRUE, lwd=4)



## --------------- ROC CURVES BY RULExSTIMULUS CONDITION ----
with(rox_rule, print(roc.test(G1, G2)))
with(rox_rule, print(roc.test(G1, L1)))
with(rox_rule, print(roc.test(G1, L2)))
with(rox_rule, print(roc.test(G2, L1)))
with(rox_rule, print(roc.test(G2, L2)))
with(rox_rule, print(roc.test(L1, L2)))

# aucs_rust <- array(NaN, c(nerr, length(rules), length(globs), length(locas)))
# dimnames(aucs_rust) <- list(diff, rules, globs, locas)
# rox_rust <- new.env()
# 
# for (e in 1:nerr){
#   for (r in rules){
#     for (g in globs){
#       for (l in locas){
#         su_blocs <- subset(blocs, currentrule==r & featg==g & featl==l)
#         
#         if (length(unique(su_blocs$cterr[,e]))==2){
#           name <- paste0(diff[e], '_', r, '_', g, '_', l)
#           
#           # ROCR approach
#           pred <- prediction(su_blocs$discrim, su_blocs$cterr[,e])
#           perf <- performance(pred, "tpr", "fpr")
#           plot(perf, main=name)
#           auc <- performance(pred, "auc")
#           aucs_rust[e,r,g,l] <- auc@y.values[[1]]
#           assign(name, pred, envir=su_pred)
#           assign(name, perf, envir=su_perf)
#     
#           # pROC approach
#           rox <- roc(su_blocs$cterr[,e], su_blocs$discrim)
#           assign(name, rox, envir=rox_rust)
#         }
#       }
#     }
#   }
# }

for (d in diff){
  print(d)
  print(aucs_rust[d,,,])
}

# pROC curve tests
with(rox_rust, { # G1 rule
  # CO-LEVEL conditions
  print(roc.test(pf_G1_face_shape, pf_G1_letter_shape))
  print(roc.test(pf_G1_face_letter, pf_G1_letter_letter))
  print(roc.test(pf_G1_face_patch, pf_G1_letter_patch))
  # OTHER-LEVEL conditions
  print(roc.test(pf_G1_face_shape, pf_G1_face_letter))
  print(roc.test(pf_G1_letter_shape, pf_G1_letter_letter))
})

with(rox_rust, { # G2 rule
  # CO-LEVEL conditions
  print(roc.test(pf_G2_face_shape, pf_G2_letter_shape))
  print(roc.test(pf_G2_face_letter, pf_G2_letter_letter))
  print(roc.test(pf_G2_face_patch, pf_G2_letter_patch))
  # OTHER-LEVEL conditions
  print(roc.test(pf_G2_face_shape, pf_G2_face_letter))
  print(roc.test(pf_G2_letter_shape, pf_G2_letter_letter))
})

with(rox_rust, { # L1 rule
  # CO-LEVEL conditions
  print(roc.test(pf_L1_face_shape, pf_L1_face_letter))
  print(roc.test(pf_L1_letter_shape, pf_L1_letter_letter))
  print(roc.test(pf_L1_noise_shape, pf_L1_noise_letter))
  # OTHER-LEVEL conditions
  print(roc.test(pf_L1_face_shape, pf_L1_letter_shape))
  print(roc.test(pf_L1_face_letter, pf_L1_letter_letter))
})

with(rox_rust, { # L2 rule
  # CO-LEVEL conditions
  print(roc.test(pf_L2_face_shape, pf_L2_face_letter))
  print(roc.test(pf_L2_letter_shape, pf_L2_letter_letter))
  print(roc.test(pf_L2_noise_shape, pf_L2_noise_letter))
  # OTHER-LEVEL conditions
  print(roc.test(pf_L2_face_shape, pf_L2_face_letter))
  print(roc.test(pf_L2_letter_shape, pf_L2_letter_letter))
})


## ------- ROC CURVES BY SUBJECTxRULExSTIMULUS CONDITION ----
aucs_srst <- 
  array(NaN, c(nerr, length(ids), length(rules), length(globs), length(locas)))
dimnames(aucs_srst) <- list(diff, ids, rules, globs, locas)

for (e in 1:nerr){
  for (i in ids){
    for (r in rules){
      for (g in globs){
        for (l in locas){
          su_blocs <- subset(blocs,ID==paste0('reknow0',i) &
                               currentrule==r & featg==g & featl==l)
          
          if (length(unique(su_blocs$cterr[,e]))==2){
            name <- paste0(diff[e], '_r0', i, '_', r, '_', g, '_', l)
            
            # ROCR approach
            pred <- prediction(su_blocs$discrim, su_blocs$cterr[,e])
            perf <- performance(pred, "tpr", "fpr")
            plot(perf, add=TRUE, main=name)
            #             auc <- performance(pred, "auc")
            #             aucs_srst[e,i,r,g,l] <- auc@y.values[[1]]
            
            #             assign(name, pred, envir=su_pred)
            #             assign(name, perf, envir=su_perf)
            
            # pROC approach
            #           rox <- roc(su_blocs$cterr[,e], su_blocs$discrim)
          }
        }
      }
    }
  }
}


## ------------- CLEAN UP ---------------------------- 
rm(i, r, g, l, su_blocs, name, pred, perf, auc, rox)
rm(ids, globs, locas, rules)
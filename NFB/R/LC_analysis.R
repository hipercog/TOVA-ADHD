library(psych)
library(ggplot2)
library(nlme) #for mixed effects models
library(broom)
# FOLLOWING ARE NEEDED, SHOULD BE ALREADY LOADED FROM LC_data_prep.R
# library(data.table) #for fast data management
# library(plyr) #for data management
# library(dplyr) #for data management

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ---
## ## ## ## ## --------------------- TESTING & MODEL FITTING ----
# Gather some data in long format
scorelong <- gather(as.data.frame(tr.not.scr), session, adj_score, 2:40, convert = TRUE, factor_key=TRUE, na.rm = TRUE)
kndlcorlg <- gather(as.data.frame(tr.not.cor), session, kendall_cor, 2:40, convert = TRUE, factor_key=TRUE, na.rm = TRUE)
#converting to a data.table object
scorelong.dt <- data.table(scorelong)

# * Method, part 1 - derive different aspects of learning by quantitative measures
#    * Consistency - derived as score monotonicity, using Kendall rank-order correlation
#    * Magnitude - derived as score mean, using outlier-resistant geometric mean
## Exploration: facet plot of score and kendall features for each subject
ggplot(scorelong, aes(x=session, y=adj_score)) +
  geom_jitter(alpha=0.4) +
  geom_smooth(method="lm", formula=y~x, se=F, alpha=0.8, linetype="dashed", size = 0.5) +
  facet_wrap(~patient) +
  theme_bw() +
  labs(x = "session#", y = "adjusted score")
ggplot(kndlcorlg, aes(x=session, y=kendall_cor)) +
  geom_jitter(alpha=0.4) +
  geom_smooth(method="lm", formula=y~x, se=F, alpha=0.8, linetype="dashed", size = 0.5) +
  facet_wrap(~patient) +
  theme_bw() +
  labs(x = "session#", y = "intra-session kendall correlation")


# ---- TESTING A RANDOM SINGLE SUBJECT ----
subj_1019 <- scorelong[which(scorelong$patient == 1019), ]
kdl_s1019 <- kndlcorlg[which(kndlcorlg$patient == 1019), ]
# subj_1019 <- dsl[which(dsl$patient == 1019), c("patient", "session", "adj_score")]
#plotting intraindividual change
ggplot(data = subj_1019, aes(x = session, y = adj_score, group = patient)) +
  geom_point() + 
  geom_line() +
  geom_smooth(method=lm, se=FALSE,colour="red", size=1) +
  xlab("session") + 
  ylab("NF Score") + ylim(min(subj_1019$adj_score),max(subj_1019$adj_score))
  # scale_x_continuous(breaks=seq(1,40,by=1))

#regress adj_score on session
linear_1019 <- lm(formula = adj_score ~ 1 + session, 
                  data = subj_1019, 
                  na.action=na.exclude)
summary(linear_1019)
lm1019.reg.linear <- as.list(coef(linear_1019))
lm1019.rse.linear <- summary(linear_1019)$sigma #to obtain residual standard error

# ---- REGRESSION FOR ALL SUBJECTS ----
#collecting regression output by patient 
indiv.reg <- scorelong.dt[,c(reg.1 = as.list(coef(lm(adj_score ~ session))), 
                            reg.1.sigma = as.list(summary(lm(adj_score ~ session))$sigma)), by=patient]
describe(as.data.frame(indiv.reg)[-1])
cor(as.data.frame(indiv.reg)[-1], use="complete.obs", method="pearson")
pairs.panels(as.data.frame(indiv.reg)[-1]) #SIGMA =  residual std error =  residual std dev
#making intraindividual change plot
ggplot(data = scorelong, aes(x = session, y = adj_score, group = patient)) +
  geom_smooth(method=lm,se=FALSE,colour="red", size=1) +
  xlab("session") + 
  ylab("NF Score") + ylim(5,50)




## ## ---- RQ1 - MAGNITUDE: LINEAR GROWTH MODELS ----

# ---- UNCONDITIONAL MEANS MODEL - BASE COMPARISON MODEL ----
um.fit <- lme(fixed = adj_score ~ 1, 
              random = ~ 1|patient, 
              data = scorelong,
              na.action = na.exclude)
summary(um.fit)
VarCorr(um.fit)
RandomEffects <- as.numeric(VarCorr(um.fit)[,1])
ICC_between <- RandomEffects[1]/(RandomEffects[1]+RandomEffects[2]) # between-person variance
# within-person variance = 100 - (ICC_between * 100)
scorelong$pred.um <- predict(um.fit)
scorelong$resid.um <- residuals(um.fit)
#plotting PREDICTED intraindividual change; overlay PROTOTYPE (average individual)
#create the function for the prototype
fun.um <- function(x) {
  as.numeric(um.fit$coefficients$fixed) + 0*x
}
#add the prototype as an additional layer
ggplot(data = scorelong, aes(x = session, y = pred.um, group = patient)) +
  ggtitle("Unconditional Means Model") +
  #  geom_point() + 
  geom_line() +
  xlab("session") + 
  ylab("PREDICTED NF Score") + ylim(0,50) +
  stat_function(fun=fun.um, color="red", size = 2)
#plotting RESIDUAL intraindividual change
ggplot(data = scorelong, aes(x = session, y = resid.um, group = patient)) +
  ggtitle("Unconditional Means Model") +
  #  geom_point() + 
  geom_line() +
  xlab("session") + 
  ylab("RESIDUAL NF Score")

# ---- FIXED LINEAR RANDOM INTERCEPT GROWTH MODEL (SESSION AS TIME)
fl.ri.fit <- lme(fixed = adj_score ~ 1 + session, 
                 random = ~ 1|patient, 
                 data = scorelong,
                 na.action = na.exclude)
summary(fl.ri.fit)
#Place individual predictions and residuals into the dataframe
scorelong$pred.fl.ri <- predict(fl.ri.fit)
scorelong$resid.fl.ri <- residuals(fl.ri.fit)
#Create a function for the prototype
fun.fl.ri <- function(x) {
  as.numeric(fl.ri.fit$coefficients$fixed[1]) + as.numeric(fl.ri.fit$coefficients$fixed[2])*x
}
#plotting PREDICTED intraindividual change
ggplot(data = scorelong, aes(x = session, y = pred.fl.ri, group = patient)) +
  ggtitle("Fixed Linear, Random Intercept") +
  #  geom_point() + 
  geom_line() +
  xlab("session") + 
  ylab("PREDICTED NF Score") + 
  stat_function(fun=fun.fl.ri, color="red", size = 2)
#plotting RESIDUAL intraindividual change
ggplot(data = scorelong, aes(x = session, y = resid.fl.ri, group = patient)) +
  ggtitle("Fixed Linear, Random Intercept") +
  #  geom_point() + 
  geom_line() +
  xlab("session") + 
  ylab("RESIDUAL NF Score")

# ---- RANDOM LINEAR FIXED INTERCEPT GROWTH MODEL (SESSION AS TIME)
# rl.fi.fit <- lme(fixed = adj_score ~ 1, 
#                  random = ~ 1 + session|patient, 
#                  data = scorelong,
#                  na.action = na.exclude)
# summary(rl.fi.fit)
# #Place individual predictions and residuals into the dataframe
# scorelong$pred.rl.fi <- predict(rl.fi.fit)
# scorelong$resid.rl.fi <- residuals(rl.fi.fit)
# #Create a function for the prototype
# fun.rl.fi <- function(x) {
#   as.numeric(rl.fi.fit$coefficients$fixed[1]) + 0*x
# }
# #plotting PREDICTED intraindividual change
# ggplot(data = scorelong, aes(x = session, y = pred.rl.fi, group = patient)) +
#   ggtitle("Random Linear, Fixed Intercept") +
#   geom_line() +
#   xlab("session") + 
#   ylab("PREDICTED NF Score") + 
#   stat_function(fun=fun.rl.fi, color="red", size = 2)
# #plotting RESIDUAL intraindividual change
# ggplot(data = scorelong, aes(x = session, y = resid.rl.fi, group = patient)) +
#   ggtitle("Random Linear, Fixed Intercept") +
#   geom_line() +
#   xlab("session") + 
#   ylab("RESIDUAL NF Score")
  
  # ---- RANDOM LINEAR SLOPES AND INTERCEPTS ----
  rl.ri.fit <- lme(fixed = adj_score ~ 1 + session,
                   random = ~ 1 + session|patient,
                   data = scorelong,
                   na.action = na.exclude)
  summary(rl.ri.fit)
  intervals(rl.ri.fit)
  #Place individual predictions and residuals into the dataframe
  scorelong$pred.rl.ri <- predict(rl.ri.fit)
  scorelong$resid.rl.ri <- residuals(rl.ri.fit)
  #Create a function for the prototype
  fun.rl.ri <- function(x) {
    as.numeric(rl.ri.fit$coefficients$fixed[1]) + as.numeric(rl.ri.fit$coefficients$fixed[2])*x
  }
  #plotting PREDICTED intraindividual change
  ggplot(data = scorelong, aes(x = session, y = pred.rl.ri, group = patient)) +
    ggtitle("Random Linear, Random Intercept") +
    #  geom_point() + 
    geom_line() +
    xlab("session") + 
    ylab("PREDICTED NF Score") +
    stat_function(fun=fun.rl.ri, color="red", size = 2)
  #plotting RESIDUAL intraindividual change
  ggplot(data = scorelong, aes(x = session, y = resid.rl.ri, group = patient)) +
    ggtitle("Random Linear, Random Intercept") +
    #  geom_point() + 
    geom_line() +
    xlab("session") + 
    ylab("RESIDUAL NF Score")
  # significance of random slopes: compare models by anova() for difference in fit between two nested models
  anova(um.fit,fl.ri.fit)
  anova(fl.ri.fit,rl.ri.fit)





## ## ---- RQ2 - plateau of performance improvement reached? ----
## ---- Method - quadratic fit for each subject
nfOrd <- setorder(data.table(nfb), TB.1.SMR1, Female0.Male2)
subOrd <- nfOrd$Part_number
quads <- matrix(NaN, nrow=n.SBJS, ncol=39)
quad.lms <- list()
quad.x2s <- vector(length = n.SBJS)
rownames(quads) <- subOrd
sex = c("F", "M")
protocol = c("TB","null", "SMR")
# subplots for each participant, from trial scores
par(mfrow=c(6,4),mar=c(1,2,1,1))
for (i in seq_along(subOrd)){
  p <- subOrd[i]
  tmp <- tr.not #outgoners(tr.not, tr.not$adj_score, thr=3)
  tmp <- ddply(subset(tmp, patient==p), .(session), summarize, adj_score=mean(adj_score))
  ttl = paste0(p, ", ", sex[nfOrd[subOrd==p,]$Gender..1.F..2.M.], ",", protocol[nfOrd[subOrd==p,]$TB.1.SMR1 + 2])
  if (i < length(subOrd) - 3){
    plot(tmp, main=ttl, xlab="", ylab="", col=2, xaxt='n')
  }else{
    plot(tmp, main=ttl, xlab="", ylab="", col=2)
  }
  tmp.lm <- with(tmp, lm(adj_score ~ session + I(session^2)))
  quad.lms[[i]] <- tidy(tmp.lm)
  quad.x2s[i] <- as.numeric(tmp.lm$coefficients[3])
  quads[i,1:nrow(tmp)] <- fitted(tmp.lm)
  lines(na.omit(quads[i, 1:nrow(tmp)]), col=3, lwd=2)
}
plot.new()
par(mfrow=c(1,1))
names(quad.x2s) <- subOrd

# compare fitted model curve direction to learning, i.e. PLATEAU v LEARN
cor.test(rl.ri.fit$coefficients$random$patient[,2], quad.x2s[sort(names(quad.x2s))])

## ---- Method - quadratic growth model of all subjects
rq.ri.fit <- lme(fixed = adj_score ~ 1 + session + I(session^2),
                 random = ~ 1 + session|patient,
                 data = scorelong,
                 na.action = na.exclude)
summary(rq.ri.fit)
#Place individual predictions and residuals into the dataframe
scorelong$pred.rq.ri <- predict(rq.ri.fit)
scorelong$resid.rq.ri <- residuals(rq.ri.fit)
#Create a function for the prototype
fun.rq.ri <- function(x) {
  as.numeric(rq.ri.fit$coefficients$fixed[1]) + 
    as.numeric(rq.ri.fit$coefficients$fixed[2])*x + 
    as.numeric(rq.ri.fit$coefficients$fixed[3])*x^2
}
#plotting PREDICTED intraindividual change
ggplot(data = scorelong, aes(x = session, y = pred.rq.ri, group = patient)) +
  ggtitle("Random Quadratic, Random Intercept") +
  geom_line() +
  xlab("session") + 
  ylab("PREDICTED NF Score") +
  stat_function(fun=fun.rq.ri, color="red", size = 2)
#plotting RESIDUAL intraindividual change
ggplot(data = scorelong, aes(x = session, y = resid.rq.ri, group = patient)) +
  ggtitle("Random Quadratic, Random Intercept") +
  geom_line() +
  xlab("session") + 
  ylab("RESIDUAL NF Score")

anova(rl.ri.fit,rq.ri.fit)



## ## ---- RQ3 - what kind of learning? ----
## ---- Parametric modelling to test Power law/exponential vs skill acquisition/phases.
# Method - fit curves from separate families
# * POWER LAW: is linear in log-log space
ggplot(scorelong, aes(x=log(session), y=log(adj_score))) + #, color=hrs_since_sleep)) + 
  geom_jitter(alpha=0.4) +
  geom_smooth(method="lm", formula=y~x, se=F, alpha=0.8, linetype="dashed", size = 0.5) +
  # scale_color_gradient(low="green", high="red") +
  facet_wrap(~patient) +
  theme_bw() +
  labs(x = "log(session#)", y = "log(adjusted score)") #, color = "hours awake")
#collecting regression output by patient 
idv.loglog.reg <- scorelong.dt[,c(reg.1 = as.list(coef(lm(log(adj_score) ~ log(session)))), 
                             reg.1.sigma = as.list(summary(lm(log(adj_score) ~ log(session)))$sigma)), by=patient]
describe(as.data.frame(idv.loglog.reg)[-1])
cor(as.data.frame(idv.loglog.reg)[-1], use="complete.obs", method="pearson")
pairs.panels(as.data.frame(idv.loglog.reg)[-1]) #SIGMA =  residual std error =  residual std dev
# OR, fit a linear model to all patients (so to compare curves with ANOVA)
loglog.rl.ri.fit <- lme(fixed = log(adj_score) ~ 1 + log(session),
                    random = ~ 1 + log(session)|patient,
                    data = scorelong,
                    na.action = na.exclude)
summary(loglog.rl.ri.fit)

# * EXPONENTIAL: linear in semi-log space
ggplot(scorelong, aes(x=session, y=log(adj_score))) + #, color=hrs_since_sleep)) + 
  geom_jitter(alpha=0.4) +
  geom_smooth(method="lm", formula=y~x, se=F, alpha=0.8, linetype="dashed", size = 0.5) +
  facet_wrap(~patient) +
  theme_bw() +
  labs(x = "session#", y = "log(adjusted score)")
#collecting regression output by patient 
idv.loglin.reg <- scorelong.dt[,c(reg.1 = as.list(coef(lm(log(adj_score) ~ session))), 
                                  reg.1.sigma = as.list(summary(lm(log(adj_score) ~ session))$sigma)), by=patient]
describe(as.data.frame(idv.loglin.reg)[-1])
cor(as.data.frame(idv.loglin.reg)[-1], use="complete.obs", method="pearson")
pairs.panels(as.data.frame(idv.loglin.reg)[-1]) #SIGMA =  residual std error =  residual std dev
# OR, fit a linear model to all patients (so to compare curves with ANOVA)
loglin.rl.ri.fit <- lme(fixed = log(adj_score) ~ 1 + session,
                        random = ~ 1 + session|patient,
                        data = scorelong,
                        na.action = na.exclude)
summary(loglin.rl.ri.fit)

# * sigmoid: linear in WHAT space??
#   * piecewise powerlaw

# * Tests - measure fit of linear regression with r²
prettyPrintModelFit("Random slope+intercept linear model:", rl.ri.fit)
prettyPrintModelFit("Random slope+intercept log-log model:", loglog.rl.ri.fit)
prettyPrintModelFit("Random slope+intercept log-linear model:", loglin.rl.ri.fit)



## ## ---- RQ4 - NON-PARAMETRIC METHODS ----
# * RQ4 - what kind of learning? Non-parametric methods to test skill acquisition model.
# * Method - fit these session-wise scores to models of possible learning trajectories: 
#                    Fitts-Posner, power-law/exponential, linear(?); with following fitting methods:

#    * cosine similarity: scaled distribution of learning characteristic values
# Fit of the chosen Fitts-Posner learning vector (0, 1, 0.5) to the Kendall-correlations, by cosine similarity
CosSimFit <- function(x){
  # Sigma_(1..n) 1 - |x_n| / n 
  mean(1 - abs(na.omit(x)))
}
CosSimFit(IdealLC[,2])
CosSimFit(ILC.nor[,2])
CosSimFit(ILC.not[,2])
CosSimFit(ILC.inv[,2])
CosSimFit(ILC.tra[,2])
#    * difference of cosine value distribution mean from 0

wilcox.test(IdealLC[,2], exact = TRUE, conf.int = TRUE)
wilcox.test(ILC.nor[,2], exact = TRUE, conf.int = TRUE)
wilcox.test(ILC.not[,2], exact = TRUE, conf.int = TRUE)
wilcox.test(ILC.inv[,2], exact = TRUE, conf.int = TRUE)
wilcox.test(ILC.tra[,2], exact = TRUE, conf.int = TRUE)


# Plot of Kendall-correlations in log-log space with linear fits, showing if close to power law
ggplot(kndlcorlg, aes(x=log(session), y=log(kendall_cor))) +
  geom_jitter(alpha=0.4) +
  geom_smooth(method="lm", formula=y~x, se=F, alpha=0.8, linetype="dashed", size = 0.5) +
  facet_wrap(~patient) +
  theme_bw() +
  labs(x = "log(session#)", y = "log(intra-session kendall correlation)")

# Linear non-parametric fits: monotonic LC vectors
CosSimFit(MonoLC[,2])
CosSimFit(MLC.nor[,2])
CosSimFit(MLC.not[,2])
CosSimFit(MLC.inv[,2])
CosSimFit(MLC.tra[,2])
CosSimFit(ssn.trl.MLC[,2])
CosSimFit(ssn.nor.MLC[,2])
CosSimFit(ssn.not.MLC[,2])
CosSimFit(ssn.inv.MLC[,2])
CosSimFit(ssn.tra.MLC[,2])
CosSimFit(trials.MLC[,2])
CosSimFit(tr.nor.MLC[,2])
CosSimFit(tr.not.MLC[,2])
CosSimFit(tr.inv.MLC[,2])
CosSimFit(tr.tra.MLC[,2])

# * Tests - test significance of distribution and model fit
# * __TODO__ - finalise output of cosine approach
wilcox.test(MonoLC[,2], exact = TRUE, conf.int = TRUE)
wilcox.test(MLC.nor[,2], exact = TRUE, conf.int = TRUE)
wilcox.test(MLC.not[,2], exact = TRUE, conf.int = TRUE)
wilcox.test(MLC.inv[,2], exact = TRUE, conf.int = TRUE)
wilcox.test(MLC.tra[,2], exact = TRUE, conf.int = TRUE)


# ---- COMPARISON of MLC & ILC ----
cor(trials.MLC$ret.agg, IdealLC[,2], use = "n")
cor(tr.nor.MLC$ret.agg, ILC.nor[,2], use = "n")
cor(tr.not.MLC$ret.agg, ILC.not[,2], use = "n")
cor(tr.inv.MLC$ret.agg, ILC.inv[,2], use = "n")
cor(tr.tra.MLC[tr.tra.cor[,1] %in% ILC.tra[,1],]$ret.agg, ILC.tra[,2], use = "n")

# ---- MONOTONIC LEARNING ----
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


# ---- IDEAL LEARNING ----
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


# * __TODO__ - Theil-Sen estimator: non-parametric regression
#    * model fit to be estimated from Theil-Sen regression
library(mblm)
#collecting regression output by patient 
TS.reg.1019 <- mblm(adj_score ~ session, subj_1019, repeated = FALSE)
TS.reg.1019 <- mblm(kendall_cor ~ session, kdl_s1019, repeated = FALSE)
summary(TS.reg.1019)
idv.TS.reg <- scorelong.dt[,c(reg.1 = as.list(coef(mblm(adj_score ~ session))), 
                              reg.1.sigma = as.list(summary(mblm(adj_score ~ session))$sigma)), by=patient]
describe(as.data.frame(idv.TS.reg)[-1])
cor(as.data.frame(idv.TS.reg)[-1], use="complete.obs", method="pearson")
pairs.panels(as.data.frame(idv.TS.reg)[-1])



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ---
## ## ## ## ## --------------------- Extract mean and CIs for each data subset ----

# Organising variables
group_by_var <- 'nom_session'
agregate_var <- 'adj_score'
agg_var_str <- "adjusted score"


# GET MEAN AND CIs FOR ALL DATA...
allT <- getTrialsMCI(getAllTrials(tr.blk), group_by_var, agregate_var) # All trials
norT <- getTrialsMCI(getNorTrials(tr.blk), group_by_var, agregate_var) # Normal trials
notT <- getTrialsMCI(getNotInvTrials(tr.blk), group_by_var, agregate_var) # Not-Inverse trials
invT <- getTrialsMCI(getInvTrials(tr.blk), group_by_var, agregate_var) # Inverse trials
traT <- getTrialsMCI(getTraTrials(tr.blk), group_by_var, agregate_var) # Transfer trials

# GET MEAN AND CIs FOR TB PROTOCOL...
df.tb <- droplevels(subset(tr.blk, patient%in%TB))
tb.allT <- getTrialsMCI(getAllTrials(df.tb), group_by_var, agregate_var) # All trials
tb.norT <- getTrialsMCI(getNorTrials(df.tb), group_by_var, agregate_var) # Normal trials
tb.notT <- getTrialsMCI(getNotInvTrials(df.tb), group_by_var, agregate_var) # Not-Inverse trials
tb.invT <- getTrialsMCI(getInvTrials(df.tb), group_by_var, agregate_var) # Inverse trials
tb.traT <- getTrialsMCI(getTraTrials(df.tb), group_by_var, agregate_var) # Transfer trials

# ...AND FOR SMR!
df.smr <- droplevels(subset(tr.blk, patient%in%SMR))
smr.allT <- getTrialsMCI(getAllTrials(df.smr), group_by_var, agregate_var) # All trials
smr.norT <- getTrialsMCI(getNorTrials(df.smr), group_by_var, agregate_var) # Normal trials
smr.notT <- getTrialsMCI(getNotInvTrials(df.smr), group_by_var, agregate_var) # Not-Inverse trials
smr.invT <- getTrialsMCI(getInvTrials(df.smr), group_by_var, agregate_var) # Inverse trials
smr.traT <- getTrialsMCI(getTraTrials(df.smr), group_by_var, agregate_var) # Transfer trials


# DERIVE MEAN AND CIs FOR THE MEN...
df.men <- droplevels(subset(tr.blk, patient%in%MEN))
men.allT <- getTrialsMCI(getAllTrials(df.men), group_by_var, agregate_var) # All trials
men.norT <- getTrialsMCI(getNorTrials(df.men), group_by_var, agregate_var) # Normal trials
men.notT <- getTrialsMCI(getNotInvTrials(df.men), group_by_var, agregate_var) # Not-Inverse trials
men.invT <- getTrialsMCI(getInvTrials(df.men), group_by_var, agregate_var) # Inverse trials
men.traT <- getTrialsMCI(getTraTrials(df.men), group_by_var, agregate_var) # Transfer trials

# ...AND THE WOMEN!
df.wmn <- droplevels(subset(tr.blk, patient%in%TB & patient%in%WMN))
wmn.allT <- getTrialsMCI(getAllTrials(df.wmn), group_by_var, agregate_var) # All trials
wmn.norT <- getTrialsMCI(getNorTrials(df.wmn), group_by_var, agregate_var) # Normal trials
wmn.notT <- getTrialsMCI(getNotInvTrials(df.wmn), group_by_var, agregate_var) # Not-Inverse trials
wmn.invT <- getTrialsMCI(getInvTrials(df.wmn), group_by_var, agregate_var) # Inverse trials
wmn.traT <- getTrialsMCI(getTraTrials(df.wmn), group_by_var, agregate_var) # Transfer trials


# DERIVE MEAN AND CIs BY PROTOCOL, FOR THE MEN...
df.tb.men <- droplevels(subset(tr.blk, patient%in%TB & patient%in%MEN))
men.tb.allT <- getTrialsMCI(getAllTrials(df.tb.men), group_by_var, agregate_var) # All trials
men.tb.norT <- getTrialsMCI(getNorTrials(df.tb.men), group_by_var, agregate_var) # Normal trials
men.tb.notT <- getTrialsMCI(getNotInvTrials(df.tb.men), group_by_var, agregate_var) # Not-Inverse trials
men.tb.invT <- getTrialsMCI(getInvTrials(df.tb.men), group_by_var, agregate_var) # Inverse trials
men.tb.traT <- getTrialsMCI(getTraTrials(df.tb.men), group_by_var, agregate_var) # Transfer trials

df.smr.men <- droplevels(subset(tr.blk, patient%in%SMR & patient%in%MEN))
men.smr.allT <- getTrialsMCI(getAllTrials(df.smr.men), group_by_var, agregate_var) # All trials
men.smr.norT <- getTrialsMCI(getNorTrials(df.smr.men), group_by_var, agregate_var) # Normal trials
men.smr.notT <- getTrialsMCI(getNotInvTrials(df.smr.men), group_by_var, agregate_var) # Not-Inverse trials
men.smr.invT <- getTrialsMCI(getInvTrials(df.smr.men), group_by_var, agregate_var) # Inverse trials
men.smr.traT <- getTrialsMCI(getTraTrials(df.smr.men), group_by_var, agregate_var) # Transfer trials

# ...AND THE WOMEN!
df.tb.wmn <- droplevels(subset(tr.blk, patient%in%TB & patient%in%WMN))
wmn.tb.allT <- getTrialsMCI(getAllTrials(df.tb.wmn), group_by_var, agregate_var) # All trials
wmn.tb.norT <- getTrialsMCI(getNorTrials(df.tb.wmn), group_by_var, agregate_var) # Normal trials
wmn.tb.notT <- getTrialsMCI(getNotInvTrials(df.tb.wmn), group_by_var, agregate_var) # Not-Inverse trials
wmn.tb.invT <- getTrialsMCI(getInvTrials(df.tb.wmn), group_by_var, agregate_var) # Inverse trials
wmn.tb.traT <- getTrialsMCI(getTraTrials(df.tb.wmn), group_by_var, agregate_var) # Transfer trials

df.smr.wmn <- droplevels(subset(tr.blk, patient%in%SMR & patient%in%WMN))
wmn.smr.allT <- getTrialsMCI(getAllTrials(df.smr.wmn), group_by_var, agregate_var) # All trials
wmn.smr.norT <- getTrialsMCI(getNorTrials(df.smr.wmn), group_by_var, agregate_var) # Normal trials
wmn.smr.notT <- getTrialsMCI(getNotInvTrials(df.smr.wmn), group_by_var, agregate_var) # Not-Inverse trials
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



library(ggplot2)
library(plotrix)
library(plotly)



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ---
## ## ## ## ## --------------------- BWRC COLOURS --------------------

# colors: 1 orange - 2 blue - 3 dark blue - 4 darker gray - 5 lighter gray - 6 green - 7 red
bwrc <- palette(c('#F0AA00','#00A0BE','#003C78','#828282','#F0F0F0','#B4BE00','#FA4100'))
bwrcol <- colorRampPalette(c('#828282','#F0F0F0','#00A0BE','#003C78'))
bwrcol <- colorRampPalette(c('#F0AA00','#00A0BE','#003C78'))
bwrcol <- colorRampPalette(c('#00A0BE','#F0AA00'))


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ---
## ## ## ## ## ---- Plot LCs in log-log space with linear regression line, showing power law fit ----
ggplot(dsl, aes(x=log(session_num), y=log(adj_score), color=hrs_since_sleep)) + 
  geom_jitter(alpha=0.4) +
  geom_smooth(method="lm", formula=y~x, se=F, alpha=0.8, linetype="dashed", size = 0.5) +
  scale_color_gradient(low="green", high="red") +
  facet_wrap(~patient) +
  theme_bw() +
  labs(x = "log(session#)", y = "log(adjusted score)", color = "hours awake")

library(standardize)
dsl$TBR.n <- scale_by(TB_ratio ~ patient, data=dsl)
ggplot(dsl, aes(x=log(session_num), y=log(adj_score), color=TBR.n)) + 
  geom_jitter(alpha=0.4) +
  geom_smooth(method="lm", formula=y~x, se=F, alpha=0.8, linetype="dashed", size = 0.5) +
  scale_color_gradient(low="green", high="red") +
  facet_wrap(~patient) +
  theme_bw() +
  labs(x = "log(session#)", y = "log(adjusted score)", color = "z(TB ratio)")

dsl.smr$TBR.n <- scale_by((theta_base / beta_base) ~ patient, data=dsl.smr)
ggplot(dsl.smr, aes(x=log(session_num), y=log(adj_score), color=TBR.n)) + 
  geom_jitter(alpha=0.4) +
  geom_smooth(method="lm", formula=y~x, se=F, alpha=0.8, linetype="dashed", size = 0.5) +
  scale_color_gradient(low="red", high="green") +
  facet_wrap(~patient) +
  theme_bw() +
  labs(x = "log(session#)", y = "log(adjusted score)", color = "z(TB ratio)")

dsl.tb$TBR.n <- scale_by((theta_base / beta_base) ~ patient, data=dsl.tb)
ggplot(dsl.tb, aes(x=log(session_num), y=log(adj_score), color=TBR.n)) + 
  geom_jitter(alpha=0.4) +
  geom_smooth(method="lm", formula=y~x, se=F, alpha=0.8, linetype="dashed", size = 0.5) +
  scale_color_gradient(low="red", high="green") +
  facet_wrap(~patient) +
  theme_bw() +
  labs(x = "log(session#)", y = "log(adjusted score)", color = "z(TB ratio)")

dsl$theta.n <- scale_by(theta_base ~ patient, data=dsl)
ggplot(dsl, aes(x=log(session_num), y=log(adj_score), color=theta.n)) + 
  geom_jitter(alpha=0.4) +
  geom_smooth(method="lm", formula=y~x, se=F, alpha=0.8, linetype="dashed", size = 0.5) +
  scale_color_gradient(low="red", high="green") +
  facet_wrap(~patient) +
  theme_bw() +
  labs(x = "log(session#)", y = "log(adjusted score)", color = "z(theta BL)")

dsl$beta.n <- scale_by(beta_base ~ patient, data=dsl)
ggplot(dsl, aes(x=log(session_num), y=log(adj_score), color=beta.n)) + 
  geom_jitter(alpha=0.4) +
  geom_smooth(method="lm", formula=y~x, se=F, alpha=0.8, linetype="dashed", size = 0.5) +
  scale_color_gradient(low="green", high="red") +
  facet_wrap(~patient) +
  theme_bw() +
  labs(x = "log(session#)", y = "log(adjusted score)", color = "z(beta BL)")

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ---
## ## ## ## ## ---------------------  EXPLORATION PLOTS ----
df <- as.data.frame(cbind(MonoLC[,2], IdealLC[,2], 
                          MLC.nor[,2], ILC.nor[,2], 
                          MLC.inv[,2], ILC.inv[,2],
                          c(MLC.tra[,2], rep(NA, 6)), c(ILC.tra[,2], rep(NA, 10))))
names(df) <- c("MonoLC", "IdealLC", "MLC-nor", "ILC-nor", "MLC-inv", "ILC-inv", "MLC-tra", "ILC-tra")
ggplot(melt(df), aes(factor(variable), value)) + geom_boxplot() + 
  geom_point(aes(color = factor(variable)), position = position_dodge(width = 0.5))

require(corrgram)
corrMatrix("MonoLC", MonoLC[,2],
           "session.MLC", ssn.trl.MLC[,2], 
           "trial.MLC", trials.MLC[,2],
           "IdealLC", IdealLC[,2],
           title = "Ideal, Monotonic learning")

corrMatrix("MonoLC", MonoLC[,2],
           "IdealLC", IdealLC[,2],
           "MLC.nor", MLC.nor[,2],
           "MLC.inv", MLC.inv[,2],
           "ILC.nor", ILC.nor[,2],
           "ILC.inv", ILC.inv[,2],
           "TB-1/SMR1", nfb$TB.1.SMR1,
           title = "Ideal, Monotonic x Training Modes, Protocol")

corrMatrix("MonoLC", MonoLC[,2],
           "IdealLC", IdealLC[,2],
           "age", nfb$Age,
           "sex", nfb$Gender..1.F..2.M.,
           "BAS", apply(rbind(nfb$BAS.Drive, nfb$BAS.Fun.Seeking, nfb$BAS.Reward.Responsiveness), 2, mean), 
           "BIS", nfb$BIS, 
           "comorbid", nfb$scales.COMORBID,
           title = "Ideal, Monotonic LCs vs background")

corrMatrix("MonoLC", MonoLC[,2],
           "IdealLC", IdealLC[,2],
           "diagnosis", nfb$diagnosis..1.ADHD..2.ADD.,
           "diag. comorbid", nfb$diagnostic.COMORBID,
           "ASRS", nfb$ASRS_total,
           "ASRS-I", nfb$ASRS_inattention_score, 
           "ASRS-H", nfb$ASRS_hyperactivity.impulsivity_score, 
           "BADDS", nfb$BADDS,
           title = "Ideal, Monotonic LCs vs diagnosis")

corrMatrix("MonoLC", MonoLC[,2],
           "IdealLC", IdealLC[,2],
           "IQ", nfb$FSIQ,
           "verbal IQ", nfb$VIQ,
           "performance IQ", nfb$PIQ,
           "digit span", nfb$Span_total_standardized, 
           "forward WM span", nfb$Span_fwd, 
           "backward WM span", nfb$Span_bwd,
           title = "Ideal, Monotonic LCs vs IQ")

idx <- SBJS %in% as.numeric(row.names(vigiPT))
corrMatrix("MonoLC", MonoLC[idx,2],
           "IdealLC", IdealLC[idx,2],
           "lability", vigiPT$Index_W,
           "waking-stage %", vigiPT$prW,
           "A-stage %", vigiPT$prA,
           "B-stage %", vigiPT$prB, 
           "C-stage %", vigiPT$prC, 
           "artefacts", vigiPT$Art,
           title = "Ideal, Monotonic LCs vs Vigilance")

par(mfrow=c(1,3))
idx <- 2
LC <- rbind(MonoLC[,2], IdealLC[,2])
ttl <- c("Monotonic", "Ideal")
boxplot(LC[idx,] ~ nfb$Gender..1.F..2.M., main = paste(ttl[idx], "v sex"), names = c("F", "M"))
boxplot(LC[idx,] ~ nfb$TB.1.SMR1, main = paste(ttl[idx], "v protocol"), names = c("TB", "SMR"))
boxplot(LC[idx,] ~ nfb$Gender..1.F..2.M. * nfb$TB.1.SMR1, main = paste(ttl[idx], "v \nprotocol x sex"), names = c("F.TB", "M.TB", "F.SMR", "M.SMR"))





## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ---
## ## ## ## ## --------------------- Plot mean trendlines of all LCs with CIs ----

# Plots of means + CIs
par(mfrow=c(2,2))
plot_meanCI(allT$M[[group_by_var]], allT$M[[agregate_var]], allT$CI, "sessions, all trials", agg_var_str, c(3, 2, 5))
plot_meanCI(norT$M[[group_by_var]], norT$M[[agregate_var]], norT$CI, "sessions, normal trials", agg_var_str, c(3, 2, 5))
plot_meanCI(notT$M[[group_by_var]], notT$M[[agregate_var]], notT$CI, "sessions, normal trials", agg_var_str, c(3, 2, 5))
plot_meanCI(invT$M[[group_by_var]], invT$M[[agregate_var]], invT$CI, "sessions, inverse trials", agg_var_str, c(3, 2, 5))
plot_meanCI(traT$M[[group_by_var]], traT$M[[agregate_var]], traT$CI, "sessions, transfer trials", agg_var_str, c(3, 2, 5))

# master plot - FIXME SEQUENCE LENGTHS
plot_2_meanCI(seq(1,38), norT$M[[agregate_var]], norT$CI, seq(22,38), invT$M[[agregate_var]], invT$CI, 
              "all sessions", agg_var_str, c("normal", "inverse"), c(3, 'black', 2, 4, 5))


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ---
## ## ## ## ## --------------------- Plot mean trendlines of male, female LCs with CIs ----

# men nor-inv plot
plot_2_meanCI(seq(1,38), men.norT$M[[agregate_var]], men.norT$CI, seq(22,38), men.invT$M[[agregate_var]], men.invT$CI, 
              "all M sessions", agg_var_str, c("normal, men", "inverse, men"), c(3, 'black', 2, 4, 5))

# women nor-inv plot
plot_2_meanCI(seq(1,38), wmn.norT$M[[agregate_var]], wmn.norT$CI, seq(22,38), wmn.invT$M[[agregate_var]], wmn.invT$CI, 
              "all F sessions", agg_var_str, c("normal, women", "inverse, women"), c(3, 'black', 2, 4, 5))


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ---
## ## ## ## ## --------------------- Plot mean trendlines of TB, SMR protocol LCs with CIs ----

# TB nor-inv plot
plot_2_meanCI(seq(1,38), tb.norT$M[[agregate_var]], tb.norT$CI, seq(22,38), tb.invT$M[[agregate_var]], tb.invT$CI, 
              "all TB sessions", agg_var_str, c("normal TB", "inverse TB"), c(3, 'black', 2, 4, 5))

# SMR nor-inv plot
plot_2_meanCI(seq(1,38), smr.norT$M[[agregate_var]], smr.norT$CI, seq(22,38), smr.invT$M[[agregate_var]], smr.invT$CI, 
              "all SMR sessions", agg_var_str, c("normal SMR", "inverse SMR"), c(3, 'black', 2, 4, 5))


# ALTOGETHER NOW - master plots
par(mfrow=c(1,2))

plot_4_meanCI(seq(1,38), men.tb.norT$M[[agregate_var]], men.tb.norT$CI, 
              seq(1,38), men.smr.norT$M[[agregate_var]], men.smr.norT$CI, 
              seq(22,38), men.tb.invT$M[[agregate_var]], men.tb.invT$CI, 
              seq(22,38), men.smr.invT$M[[agregate_var]], men.smr.invT$CI, 
              "MEN - all sessions, by protocol", agg_var_str, 
              c("normal TB", "normal SMR", "inverse TB", "inverse SMR"), c(2, 8, 6, 5, 1, 4, 3, 5))

plot_4_meanCI(seq(1,38), wmn.tb.norT$M[[agregate_var]], wmn.tb.norT$CI, 
              seq(1,38), wmn.smr.norT$M[[agregate_var]], wmn.smr.norT$CI, 
              seq(22,38), wmn.tb.invT$M[[agregate_var]], wmn.tb.invT$CI, 
              seq(22,38), wmn.smr.invT$M[[agregate_var]], wmn.smr.invT$CI, 
              "WMN - all sessions, by protocol", agg_var_str, 
              c("normal TB", "normal SMR", "inverse TB", "inverse SMR"), c(2, 8, 6, 5, 1, 4, 3, 5))


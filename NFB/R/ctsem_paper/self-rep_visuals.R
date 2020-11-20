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


## BIS x LC ----
eqn <- lm_eqn(dat, lm(slope.score ~ BIS, data = dat))
eqTB <- lm_eqn(dat, lm(slope.score ~ BIS, data = filter(dat, regime == "TBR")))
eqSM <- lm_eqn(dat, lm(slope.score ~ BIS, data = filter(dat, regime == "SMR")))

ggplot(dat, aes(x = BIS, y = slope.score)) +
  geom_point(aes(pch = regime)) +
  geom_smooth(method = "lm") +
  ylab("Slope of learning curve") +
  ylim(-0.25, 0.6) +
  annotate("text",
           x = max(dat$BIS),
           y = max(dat$slope.score),
           hjust = 1,
           label = eqn,
           parse = TRUE) +
  theme_minimal()
ggsave(file.path(odir, "LCxBIS.svg"))

ggplot(dat, aes(x = BIS, y = slope.score, color = learner)) +
  geom_point(aes(pch = learner)) +
  geom_smooth(method = "lm") +
  ylab("Slope of learning curve") +
  scale_colour_manual(values = colv[c(1, 2)]) +
  theme_minimal()

ggplot(dat, aes(x = BIS, y = slope.score, color = regime)) +
  geom_point(aes(pch = regime)) +
  geom_smooth(method = "lm") +
  ylab("Slope of learning curve") +
  ylim(-0.25, 0.7) +
  scale_colour_manual(values = colv[c(9, 10)]) +
  annotate("text",
           x = max(dat$BIS),
           y = max(dat$slope.score) + 0.17,
           hjust = 1,
           label = eqSM,
           parse = TRUE, 
           color = colv[9]) +
  annotate("text",
           x = max(dat$BIS),
           y = max(dat$slope.score) + 0.1,
           hjust = 1,
           label = eqTB,
           parse = TRUE, 
           color = colv[10]) +
  theme_minimal()
ggsave(file.path(odir, "LCxBISxRegime.svg"))


## DES x LC ----
eqn <- lm_eqn(dat, lm(slope.score ~ DES, data = dat))
eqTB <- lm_eqn(dat, lm(slope.score ~ DES, data = filter(dat, regime == "TBR")))
eqSM <- lm_eqn(dat, lm(slope.score ~ DES, data = filter(dat, regime == "SMR")))

ggplot(dat, aes(x = DES, y = slope.score)) +
  geom_point(aes(pch = regime)) +
  geom_smooth(method = "lm") +
  ylab("Slope of learning curve") +
  ylim(-0.25, 0.6) +
  annotate("text",
           x = max(dat$DES),
           y = max(dat$slope.score),
           hjust = 1,
           label = eqn,
           parse = TRUE) +
  theme_minimal()
ggsave(file.path(odir, "LCxDES.svg"))

ggplot(dat, aes(x = DES, y = slope.score, color = learner)) +
  geom_point(aes(pch = learner)) +
  geom_smooth(method = "lm") +
  ylab("Slope of learning curve") +
  scale_colour_manual(values = colv[c(1, 2)]) +
  theme_minimal()

ggplot(dat, aes(x = DES, y = slope.score, color = regime)) +
  geom_point(aes(pch = regime)) +
  geom_smooth(method = "lm") +
  ylab("Slope of learning curve") +
  ylim(-0.25, 0.7) +
  scale_colour_manual(values = colv[c(9, 10)]) +
  annotate("text",
           x = max(dat$DES),
           y = max(dat$slope.score) + 0.17,
           hjust = 1,
           label = eqSM,
           parse = TRUE, 
           color = colv[9]) +
  annotate("text",
           x = max(dat$DES),
           y = max(dat$slope.score) + 0.1,
           hjust = 1,
           label = eqTB,
           parse = TRUE, 
           color = colv[10]) +
  theme_minimal()
ggsave(file.path(odir, "LCxDESxRegime.svg"))


## GAD x LC ----
eqn <- lm_eqn(dat, lm(slope.score ~ GAD, data = dat))
eqTB <- lm_eqn(dat, lm(slope.score ~ GAD, data = filter(dat, regime == "TBR")))
eqSM <- lm_eqn(dat, lm(slope.score ~ GAD, data = filter(dat, regime == "SMR")))

ggplot(dat, aes(x = GAD, y = slope.score)) +
  geom_point(aes(pch = regime)) +
  geom_smooth(method = "lm") +
  ylab("Slope of learning curve") +
  ylim(-0.25, 0.6) +
  annotate("text",
           x = max(dat$GAD),
           y = max(dat$slope.score),
           hjust = 1,
           label = eqn,
           parse = TRUE) +
  theme_minimal()
ggsave(file.path(odir, "LCxGAD.svg"))

ggplot(dat, aes(x = GAD, y = slope.score)) +
  geom_point(aes(pch = learner)) +
  geom_smooth(method = "lm") +
  ylab("Slope of learning curve") +
  scale_colour_manual(values = colv[c(1, 2)]) +
  theme_minimal()

ggplot(dat, aes(x = GAD, y = slope.score, color = regime)) +
  geom_point(aes(pch = regime)) +
  geom_smooth(method = "lm") +
  ylab("Slope of learning curve") +
  ylim(-0.25, 0.7) +
  scale_colour_manual(values = colv[c(9, 10)]) +
  annotate("text",
           x = max(dat$GAD),
           y = max(dat$slope.score) + 0.17,
           hjust = 1,
           label = eqSM,
           parse = TRUE, 
           color = colv[9]) +
  annotate("text",
           x = max(dat$GAD),
           y = max(dat$slope.score) + 0.1,
           hjust = 1,
           label = eqTB,
           parse = TRUE, 
           color = colv[10]) +
  theme_minimal()
ggsave(file.path(odir, "LCxGADxRegime.svg"))


#### ASRS DATA ----
asrs <- read.csv(file.path('data', 'ctsem_paper', 'ASRSxLC.csv'), sep = '\t')
asrs$DES <- des$DES < 25
df <- asrs %>%
  pivot_longer(cols = contains(c("inattention", "hyper")), names_to = "ASRS")
df$ASRS <- as.factor(df$ASRS)
df$Learner <- as.factor(df$Learner)
levels(df$Learner) <- c("nonLearner", "Learner")
df$DES <- as.factor(df$DES)
levels(df$DES) <- c("DES≥25", "DES<25")
df$status.ASRS <- with(df, Learner:ASRS)
df$DES.ASRS <- with(df, DES:ASRS)

#---- FIGURE FOR ASRS x LEARNER STATUS ----
ggplot(df, aes(x = status.ASRS, y = value)) +
  geom_jitter(aes(fill=Regime), width = 0.2, height = 0.1, size=2, shape=21, stroke=0) +
  scale_colour_manual(values = c(colv[9], colv[10]), aesthetics = c("color", "fill")) +
  geom_boxplot(width=0.2, outlier.shape = NA, alpha = 0.4) +
  labs(title="ASRS factors before and after 30 sessions of treatment", x="ASRS", y = "Score") +
  theme_minimal() +
  theme(
    # legend.position="none",
    axis.title.x = element_text(size=11),
    axis.title.y= element_text(size=11),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave(file.path(odir, "LCxASRS.svg"))


#---- FIGURE FOR ASRS x DES LEVEL ----
ggplot(df, aes(x = DES.ASRS, y = value)) +
  geom_jitter(aes(fill=Regime), width = 0.2, height = 0.1, size=2, shape=21, stroke=0) +
  scale_colour_manual(values = c(colv[9], colv[10]), aesthetics = c("color", "fill")) +
  geom_boxplot(width=0.2, outlier.shape = NA, alpha = 0.4) +
  labs(title="ASRS factors before and after 30 sessions of treatment", x="ASRS", y = "Score") +
  theme_minimal() +
  theme(
    # legend.position="none",
    axis.title.x = element_text(size=11),
    axis.title.y= element_text(size=11),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave(file.path(odir, "DESxASRS.svg"))

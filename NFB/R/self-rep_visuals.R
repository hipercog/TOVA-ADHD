library(tidyverse)
library(ggridges)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(RColorBrewer)
source('R/znbnz_visuals.R')

odir <- "/home/bcowley/Dropbox/project_CENT/Publications/WIP - NFB Learning Curves/Figures/"

colv <- brewer.pal(n = 12, name = "Paired")

#---- FIGURE FOR DES x LC SLOPE
dat <- read.csv(file.path('data', 'DESxLC.csv'), sep = '\t')

eqn <- lm_eqn(dat, lm(slope.score ~ DES, data = dat))
eqTB <- lm_eqn(dat, lm(slope.score ~ DES, data = filter(dat, regime == "TB")))
eqSM <- lm_eqn(dat, lm(slope.score ~ DES, data = filter(dat, regime == "SMR")))
eqLR <- lm_eqn(dat, lm(slope.score ~ DES, data = filter(dat, learner == "Learner")))
eqNL <- lm_eqn(dat, lm(slope.score ~ DES, data = filter(dat, learner == "Non-Learner")))

ggplot(dat, aes(x = DES, y = slope.score)) +
  geom_point(aes(pch = dat$regime)) +
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
ggsave(paste0(odir, "LCxDES.svg"))

ggplot(dat, aes(x = DES, y = slope.score, color = regime)) +
  geom_point(aes(pch = dat$regime)) +
  geom_smooth(method = "lm") +
  ylab("Slope of learning curve") +
  ylim(-0.25, 0.6) +
  scale_colour_manual(values = colv[c(9, 10)]) +
  annotate("text",
           x = max(dat$DES),
           y = max(dat$slope.score) + 0.07,
           hjust = 1,
           label = eqTB,
           parse = TRUE) +
  annotate("text",
           x = max(dat$DES),
           y = max(dat$slope.score),
           hjust = 1,
           label = eqSM,
           parse = TRUE) +
  theme_minimal()
ggsave(paste0(odir, "LCxDESxRegime.svg"))

ggplot(dat, aes(x = DES, y = slope.score, color = learner)) +
  geom_point(aes(pch = dat$learner)) +
  geom_smooth(method = "lm") +
  ylab("Slope of learning curve") +
  ylim(-0.25, 0.6) +
  scale_colour_manual(values = colv[c(9, 10)]) +
  annotate("text",
           x = max(dat$DES),
           y = max(dat$slope.score) + 0.07,
           hjust = 1,
           label = eqLR,
           parse = TRUE) +
  annotate("text",
           x = max(dat$DES),
           y = max(dat$slope.score),
           hjust = 1,
           label = eqNL,
           parse = TRUE) +
  theme_minimal()
ggsave(paste0(odir, "LCxDESxLearner.svg"))

#---- FIGURE FOR ASRS x LEARNER STATUS
dat <- read.csv(file.path('data', 'ASRSxLC.csv'), sep = '\t')
df <- dat %>%
  pivot_longer(cols = contains(c("inattention", "hyper")), names_to = "ASRS")
df$ASRS <- as.factor(df$ASRS)
df$Learner <- as.factor(df$Learner)
levels(df$Learner) <- c("nonLearner", "Learner")
df$status.ASRS <- with(df, Learner:ASRS)

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
ggsave(paste0(odir, "LCxASRS.svg"))

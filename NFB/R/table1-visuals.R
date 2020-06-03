library(tidyverse)
library(ggridges)
library(ggplot2)
library(viridis)
library(hrbrthemes)

df <- read.csv('data/table1.csv', sep = "\t")
df$Patient <- as.factor(df$Patient)
df$Plateau <- as.factor(df$Quadratic.sign > 0)
df$status.regime <- with(df, Status:Regime)
df$status.regime.plateau <- with(df, Status:Regime:Plateau)
str(df)
n <- nrow(df)

x <- df$z.growth.model.slope
xlbl <- "z_slope"
# y <- df$z.inverse.slope
# ylbl <- "z_inv_slope"
x <- df$Quadratic.sign
xlbl <- "plateau"
y <- df$Change.in.BPR
ylbl <- "BPR change"
LBL <- df$Regime

plot(x,y, type='p', xlab = xlbl, ylab = ylbl)

source('R/znbnz_visuals.R')
clustered_scatter(x, y, xlab = xlbl, ylab = ylbl, ptnames = LBL)

df %>% select(c(2, 4, 6, 7, 11, 12, 13)) %>% pairs()



odir <- "/home/bcowley/Dropbox/project_CENT/Publications/WIP - NFB Learning Curves/Figures/"
# Plot
ggplot(df, aes(x = z.growth.model.slope, y = status.regime, fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE,
    quantiles = c(0.025, 0.5, 0.975),
    quantile_lines = TRUE
  ) +
  scale_fill_manual(
    name = "Probability", values = c("#666666A0", "#E0E0E0A0", "#E0E0E0A0", "#999999A0"),
    labels = c("(0, 0.025]", "(0.025, 0.5]", "(0.5, 0.975]", "(0.975, 1]")
  ) +
  xlab('Growth curve slope of normal trials') +
  ylab('Status x Regime') +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    axis.title.x = element_text(size=11),
    axis.title.y= element_text(size=11),
    axis.text.y = element_text(angle = 45, hjust = 1)
  )
ggsave(paste0(odir, "learnXregime_normal.svg"))

ggplot(df, aes(x = z.inverse.slope, y = status.regime, fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE,
    quantiles = c(0.025, 0.5, 0.975),
    quantile_lines = TRUE
  ) +
  scale_fill_manual(
    name = "Probability", values = c("#666666A0", "#E0E0E0A0", "#E0E0E0A0", "#999999A0"),
    labels = c("(0, 0.025]", "(0.025, 0.5]", "(0.5, 0.975]", "(0.975, 1]")
  ) +
  xlab('Growth curve slope of inverse trials') +
  ylab('Status x Regime') +
  theme_ipsum()+
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    axis.title.x = element_text(size=11),
    axis.title.y= element_text(size=11),
    axis.text.y = element_text(angle = 45, hjust = 1)
  )
ggsave(paste0(odir, "learnXregime_inverse.svg"))

ggplot(df, aes(x = Gain, y = status.regime, fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE,
    quantiles = c(0.025, 0.5, 0.975),
    quantile_lines = TRUE
  ) +
  scale_fill_manual(
    name = "Probability", values = c("#666666A0", "#E0E0E0A0", "#E0E0E0A0", "#999999A0"),
    labels = c("(0, 0.025]", "(0.025, 0.5]", "(0.5, 0.975]", "(0.975, 1]")
  ) +
  xlab('Gain of normal trials') +
  ylab('Status x Regime') +
  theme_ipsum()+
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    axis.title.x = element_text(size=11),
    axis.title.y= element_text(size=11),
    axis.text.y = element_text(angle = 45, hjust = 1)
  )
ggsave(paste0(odir, "learnXregime_gain.svg"))

ggplot(df, aes(x = Change.in.BPR, y = status.regime, fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE,
    quantiles = c(0.025, 0.5, 0.975),
    quantile_lines = TRUE
  ) +
  scale_fill_manual(
    name = "Probability", values = c("#666666A0", "#E0E0E0A0", "#E0E0E0A0", "#999999A0"),
    labels = c("(0, 0.025]", "(0.025, 0.5]", "(0.5, 0.975]", "(0.975, 1]")
  ) +
  xlab('Change in BPR') +
  ylab('Status x Regime') +
  theme_ipsum()+
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    axis.title.x = element_text(size=11),
    axis.title.y= element_text(size=11),
    axis.text.y = element_text(angle = 45, hjust = 1)
  )
ggsave(paste0(odir, "learnXregime_BPR.svg"))

ggplot(df, aes(x = BIS, y = status.regime, fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE,
    quantiles = c(0.025, 0.5, 0.975),
    quantile_lines = TRUE
  ) +
  scale_fill_manual(
    name = "Probability", values = c("#666666A0", "#E0E0E0A0", "#E0E0E0A0", "#999999A0"),
    labels = c("(0, 0.025]", "(0.025, 0.5]", "(0.5, 0.975]", "(0.975, 1]")
  ) +
  xlab('Behavioural Inhibition Scale') +
  ylab('Status x Regime') +
  theme_ipsum()+
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    axis.title.x = element_text(size=11),
    axis.title.y= element_text(size=11),
    axis.text.y = element_text(angle = 45, hjust = 1)
  )
ggsave(paste0(odir, "learnXregime_BIS.svg"))

ggplot(df, aes(x = DES, y = status.regime, fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE,
    quantiles = c(0.025, 0.5, 0.975),
    quantile_lines = TRUE
  ) +
  scale_fill_manual(
    name = "Probability", values = c("#666666A0", "#E0E0E0A0", "#E0E0E0A0", "#999999A0"),
    labels = c("(0, 0.025]", "(0.025, 0.5]", "(0.5, 0.975]", "(0.975, 1]")
  ) +
  xlab('Dissociative Experiences Scale') +
  ylab('Status x Regime') +
  theme_ipsum()+
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    axis.title.x = element_text(size=11),
    axis.title.y= element_text(size=11),
    axis.text.y = element_text(angle = 45, hjust = 1)
  )
ggsave(paste0(odir, "learnXregime_DES.svg"))
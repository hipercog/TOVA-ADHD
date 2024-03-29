library(tidyverse)
library(ggplot2)
library(gghalves)
library(ggridges)
library(ggpubr)
library(viridis)
library(hrbrthemes)
library(patchwork)
library(here)

odir <- file.path(here(), 'figures')
SAVE <- FALSE
comps = list( c("ADHD", "control") )
pvals <- c(0.046, 0.024, 0.037)

dat <- readxl::read_xlsx(file.path(here(), 'data', 'tova_parsed__2013-08-22_12-15-08.xlsx'))

df <-
  dat %>% filter(SESNUM == 1 & SEITOTAL < 3) %>% 
  select(-SESNUM, -GroupName_prepost, -DOB, -DoB) #PatientHealthy == 2 | TDATE < "2013-01-01")
df$PatientHealthy <- as.factor(df$PatientHealthy)
levels(df$PatientHealthy) <- c("Control", "ADHD")
df$TestWLControl <- as.factor(df$TestWLControl)
df <- rename(df, Group = GroupName)
df$Group <- as.factor(df$Group)
df$GENDER <- as.factor(df$GENDER)


# Plot after simple group comparisons
violnCOM <- ggplot(data = df, aes(x = PatientHealthy, y = COMSSH1, fill = PatientHealthy)) +
  geom_violin() + 
  geom_point(pch = 21, position = position_jitter(width = 0.15, height = 0.05)) +
  xlab(NULL) +
  ylab('H1 Commission errors\n (standard scores)') +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12), 
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "cm"))

violnDPR <- ggplot(data = df, aes(x = PatientHealthy, y = DPRSSH1, fill = PatientHealthy)) +
  geom_violin() + 
  geom_point(pch = 21, position = position_jitter(width = 0.15, height = 0.05)) +
  xlab(NULL) +
  ylab('H1 D\' (standard scores)') +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12), 
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "cm"))


violnVAR <- ggplot(data = df, aes(x = PatientHealthy, y = VARSSH2, fill = PatientHealthy)) +
  geom_violin() + 
  geom_point(pch = 21, position = position_jitter(width = 0.15, height = 0.05)) +
  xlab(NULL) +
  ylab('H2 RTV (standard scores)') +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12), 
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "cm"))



# Plot after time x group comparisons - COMERRH1, COMERRH2, OMERRH1, OMERRH2, DPRIMEH1, DPRIMEH2

# First dataframe includes all DVs but they don't filter to same subset so are handled individually below
# dfl <-
#   df %>% select(1:8, contains("OMSSH") | starts_with("DPRSSH")) %>% 
#   pivot_longer(cols = ends_with(c("H1", "H2")), 
#                names_to = c(".value", "Half"), 
#                names_pattern = "([A-Z]*)H([12])") %>%
#   mutate(Half = factor(Half, labels = c("H1", "H2"))) %>%
#   mutate(group.time = PatientHealthy:Half) %>%
#   filter(OMSS > mean(OMSS) - sd(OMSS) * 3 & OMSS < mean(OMSS) + sd(OMSS) * 3 &
#            COMSS > mean(COMSS) - sd(COMSS) * 3 & COMSS < mean(COMSS) + sd(COMSS) * 3 &
#            DPRSS > mean(DPRSS) - sd(DPRSS) * 3 & DPRSS < mean(DPRSS) + sd(DPRSS) * 3)

dfl <- 
  df %>% select(1:8, starts_with("OMSSH")) %>% 
  pivot_longer(cols = ends_with(c("H1", "H2")), 
               names_to = c(".value", "Half"), 
               names_pattern = "([A-Z]*)H([12])") %>%
  mutate(Half = factor(Half, labels = c("H1", "H2"))) %>%
  mutate(group.time = PatientHealthy:Half) %>%
  filter(OMSS > mean(OMSS) - (sd(OMSS)) & OMSS < mean(OMSS) + (sd(OMSS)))

ridgeOMS <- ggplot(dfl, aes(x = OMSS, y = group.time, fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE,
    quantiles = c(0.025, mean(dfl$OMSS < mean(dfl$OMSS)), 0.975),
    quantile_lines = TRUE,
    vline_color = "red3",
    jittered_points = TRUE,
    position = position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 2.5, point_alpha = 0.8
  ) +
  scale_fill_manual(
    name = "Probability", values = c("#666666A0", "#E0E0E0A0", "#E0E0E0A0", "#999999A0")
  ) +
  xlab('Omission errors (standard scores)') +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 12),
    plot.margin = unit(c(0,0,0,0), "cm")
  )

dfl <-
  df %>% select(1:8, starts_with("COMSSH")) %>% 
  pivot_longer(cols = ends_with(c("H1", "H2")), 
               names_to = c(".value", "Half"), 
               names_pattern = "([A-Z]*)H([12])") %>%
  mutate(Half = factor(Half, labels = c("H1", "H2"))) %>%
  mutate(group.time = PatientHealthy:Half) %>%
  filter(COMSS > mean(COMSS) - sd(COMSS) * 3 & COMSS < mean(COMSS) + sd(COMSS) * 3)

ridgeCOM <- ggplot(dfl, aes(x = COMSS, y = group.time, fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE,
    quantiles = c(0.025, mean(dfl$COMSS < mean(dfl$COMSS)), 0.975),
    quantile_lines = TRUE,
    vline_color = "red3",
    jittered_points = TRUE,
    position = position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 2.5, point_alpha = 0.7
  ) +
  scale_fill_manual(
    name = "Probability", values = c("#666666A0", "#E0E0E0A0", "#E0E0E0A0", "#999999A0")
  ) +
  xlab('Commission errors (standard scores)') +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 12),
    plot.margin = unit(c(0,0,0,0), "cm")
  )


dfl <-
  df %>% select(1:8, starts_with("DPRSSH")) %>% 
  pivot_longer(cols = ends_with(c("H1", "H2")), 
               names_to = c(".value", "Half"), 
               names_pattern = "([A-Z]*)H([12])") %>%
  mutate(Half = factor(Half, labels = c("H1", "H2"))) %>%
  mutate(group.time = PatientHealthy:Half) %>%
  filter(DPRSS > mean(DPRSS) - sd(DPRSS) * 3 & DPRSS < mean(DPRSS) + sd(DPRSS) * 3) 

ridgeDPR <- ggplot(dfl, aes(x = DPRSS, y = group.time, fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE,
    quantiles = c(0.025, mean(dfl$DPRSS < mean(dfl$DPRSS)), 0.975),
    quantile_lines = TRUE,
    vline_color = "red3",
    jittered_points = TRUE,
    position = position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 2.5, point_alpha = 0.8
  ) +
  scale_fill_manual(
    name = "Probability", values = c("#666666A0", "#E0E0E0A0", "#E0E0E0A0", "#999999A0")
  ) +
  xlab('D\' (standard scores)') +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 12),
    plot.margin = unit(c(0,0,0,0), "cm")
  )


if(SAVE){
  ggsave(file.path(odir, "comerr_ss_h1.svg"), plot = violnCOM)
  ggsave(file.path(odir, "dprime_ss_h1.svg"), plot = violnDPR)
  ggsave(file.path(odir, "rtv_ss_h2.svg"), plot = violnVAR)
  ggsave(file.path(odir, "omerr_ss_ridge.svg"), plot = ridgeOMS)
  ggsave(file.path(odir, "comerr_ss_ridge.svg"), plot = ridgeCOM)
  ggsave(file.path(odir, "dprime_ss_ridge.svg"), plot = ridgeDPR)
}

(violnCOM + violnDPR + violnVAR + ridgeCOM + ridgeDPR + ridgeOMS) + 
  plot_layout(guides="collect") + plot_annotation(tag_levels = 'A')

library(tidyverse)

df <- read.csv('data/table1.csv', sep = "\t")
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

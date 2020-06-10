## ## ## ## ## --------------------- PLOTTING FUNCTIONS ----

corrMatrix <- function(..., title = "Correlation matrix")
{
  require(corrgram)
  x <- list(...)
  # assume data and names have been passed as name-value pairs
  df <- as.data.frame(x[seq(2,length(x),2)])
  colnames(df) <- x[seq(1,length(x),2)]
  
  corrgram(df, order = TRUE, 
           lower.panel = panel.ellipse,
           upper.panel = panel.conf, 
           text.panel = panel.txt,
           diag.panel = panel.minmax, 
           main = title,
           gap = 1)
}

# Plot one mean and CI LC
plot_meanCI <- function(xmarks, mline, cilines, xlabel, ylabel, colgrp)
{
  
  plot(range(xmarks), range(c(cilines$V4, cilines$V5)), type="n", xaxt="n", xlab=xlabel, ylab=ylabel)
  axis(1, at = xmarks)
  polygon(c(xmarks, rev(xmarks)), c(cilines$V5, rev(cilines$V4)), col=colgrp[3], border=NA)
  lines(xmarks, t(mline), type="l", lwd=2, col=colgrp[1])
  with(cilines, lines(xmarks, V4, type="l", lwd=1, col=colgrp[2]))
  with(cilines, lines(xmarks, V5, type="l", lwd=1, col=colgrp[2]))
}

# Plot two means and their CIs
plot_2_meanCI <- function(x1s, ml1, ci1, 
                          x2s, ml2, ci2, 
                          xlabel, ylabel, lgdstr, colgrp)
{
  xmarks <- unique(c(x1s, x2s))
  plot(range(xmarks),
       range(c(ci1$V4, ci1$V5, ci2$V4, ci2$V5)), 
       type="n", xaxt="n", xlab=xlabel, ylab=ylabel)
  axis(1, at=xmarks)
  
  polygon(c(x1s, rev(x1s)), c(ci1$V5, rev(ci1$V4)), col=colgrp[5], border=NA)
  with(ci1, lines(x1s, V4, type="l", lwd=1, col=colgrp[3]))
  with(ci1, lines(x1s, V5, type="l", lwd=1, col=colgrp[3]))
  lines(x1s, ml1, type="l", lwd=2, col=colgrp[1])
  
  with(ci2, lines(x2s, V4, type="l", lwd=1, col=colgrp[4]))
  with(ci2, lines(x2s, V5, type="l", lwd=1, col=colgrp[4]))
  lines(x2s, ml2, type="l", lwd=2, col=colgrp[2])
  
  legend("topleft", inset=c(0.05,0.05), y.intersp=1.5, cex=1, bty="n",
         legend=lgdstr, lty=1, lwd=2, col=colgrp[1:2])
}

# Plot three means and their CIs
plot_3_meanCI <- function(x1s, ml1, ci1, 
                          x2s, ml2, ci2, 
                          x3s, ml3, ci3, 
                          xlabel, ylabel, lgdstr, colgrp)
{
  xmarks <- unique(c(x1s, x2s, x3s))
  plot(range(xmarks),
       range(c(ci1$V4, ci1$V5, ci2$V4, ci2$V5, ci3$V4, ci3$V5)), 
       type="n", xaxt="n", xlab=xlabel, ylab=ylabel)
  axis(1, at=xmarks)
  
  polygon(c(x1s, rev(x1s)), c(ci1$V5, rev(ci1$V4)), col=colgrp[2], border=NA)
  lines(x1s, ml1, type="l", lwd=2, col=colgrp[1])
  
  with(ci2, lines(x2s, V4, type="l", lwd=1, col=colgrp[4]))
  with(ci2, lines(x2s, V5, type="l", lwd=1, col=colgrp[4]))
  lines(x2s, ml2, type="l", lwd=2, col=colgrp[3])
  
  with(ci3, lines(x3s, V4, type="l", lwd=1, col=colgrp[6]))
  with(ci3, lines(x3s, V5, type="l", lwd=1, col=colgrp[6]))
  lines(x3s, ml3, type="l", lwd=2, col=colgrp[5])
  
  legend("topleft", inset=c(0.05,0.05), y.intersp=1.5, cex=1, bty="n",
         legend=lgdstr, lty=1, lwd=2, col=colgrp[c(1,3,5)])
}


# Plot four means and their CIs
plot_4_meanCI <- function(x1s, ml1, ci1, 
                          x2s, ml2, ci2, 
                          x3s, ml3, ci3, 
                          x4s, ml4, ci4, 
                          xlabel, ylabel, lgdstr, colgrp)
{
  xmarks <- unique(c(x1s, x2s, x3s, x4s))
  plot(range(xmarks),
       range(c(ci1$V4, ci1$V5, ci2$V4, ci2$V5, ci3$V4, ci3$V5, ci4$V4, ci4$V5)), 
       type="n", xaxt="n", xlab=xlabel, ylab=ylabel)
  axis(1, at=xmarks)
  
  polygon(c(x1s, rev(x1s)), c(ci1$V5, rev(ci1$V4)), col=colgrp[2], border=NA)
  lines(x1s, ml1, type="l", lwd=2, col=colgrp[1])
  
  polygon(c(x2s, rev(x2s)), c(ci2$V5, rev(ci2$V4)), col=colgrp[4], border=NA)
  lines(x2s, ml2, type="l", lwd=2, col=colgrp[3])
  
  with(ci3, lines(x3s, V4, type="l", lwd=1, col=colgrp[6]))
  with(ci3, lines(x3s, V5, type="l", lwd=1, col=colgrp[6]))
  lines(x3s, ml3, type="l", lwd=2, col=colgrp[5])
  
  with(ci4, lines(x4s, V4, type="l", lwd=1, col=colgrp[8]))
  with(ci4, lines(x4s, V5, type="l", lwd=1, col=colgrp[8]))
  lines(x4s, ml4, type="l", lwd=2, col=colgrp[7])
  
  legend("topleft", inset=c(0.05,0.05), y.intersp=1.5, cex=1, bty="n",
         legend=lgdstr, lty=1, lwd=2, col=colgrp[c(1,3,5,7)])
}


# Create text that can be added to ggplot LM ablines
# function to create the text equation
lm_eqn <- function(df, lm_object) {
  eq <-
    substitute(
      italic(y) == a + b %.% italic(x) * "," ~  ~ italic(r) ^ 2 ~ "=" ~ r2,
      list(
        a = format(coef(lm_object)[1], digits = 2),
        b = format(coef(lm_object)[2], digits = 2),
        r2 = format(summary(lm_object)$r.squared, digits = 3),
        digits = 2
      )
    )
  as.character(as.expression(eq))
}


# VISUALIZING CLUSTERS IN 2D
# From https://www.r-bloggers.com/visualizing-clusters/

within_var = function(I){
  I0=which(I==0)
  I1=which(I==1)
  xbar0=mean(x[I0])
  xbar1=mean(x[I1])
  ybar0=mean(y[I0])
  ybar1=mean(y[I1])
  w=sum(I0)*sum( (x[I0]-xbar0)^2+(y[I0]-ybar0)^2 )+
    sum(I1)*sum( (x[I1]-xbar1)^2+(y[I1]-ybar1)^2 )
  return(c(I,w))
}

base2 = function(z,n=10){
  Base.b=rep(0,n)
  ndigits=(floor(logb(z, base=2))+1)
  for(i in 1:ndigits){
    Base.b[ n-i+1]=(z%%2)
    z=(z%/%2)}
  return(Base.b)
}

cluster_viz = function(x, y, indices, xlbl, ylbl, ptnm){
  library(RColorBrewer)
  
  CL2palette=rev(brewer.pal(n = 9, name = "RdYlBu"))
  CL2f=CL2palette[c(1,9)]
  
  plot(x, y, pch=19, xlab=xlbl, ylab=ylbl, cex=2, col=CL2f[1+indices])
  text(x, y, label=ptnm, cex = 0.5)
  
  CL2cl=CL2palette[c(3,7)]
  I0=which(indices==0)
  I1=which(indices==1)
  xbar0=mean(x[I0])
  xbar1=mean(x[I1])
  ybar0=mean(y[I0])
  ybar1=mean(y[I1])
  segments(x[I0],y[I0],xbar0,ybar0,col=CL2cl[1])
  segments(x[I1],y[I1],xbar1,ybar1,col=CL2cl[2])
  points(xbar0,ybar0,pch=19,cex=1.5,col=CL2cl[1])
  points(xbar1,ybar1,pch=19,cex=1.5,col=CL2cl[2])
}

clustered_scatter = function(x, y, xlab, ylab, ptnames, clustn = 2, clustrng = "kmeans", meth = ""){
  
  xn = length(x)
  yn = length(y)
  if (xn != yn){
    stop("x and y are not of equal length")
  }else{
    n = xn
  }
  
  if (clustrng == "min var")
  {
    L = function(x) within_var(base2(x, n))
    S = sapply(1:(2^10),L)
    I = S[1:n,which.min(S[n+1,])]
  }else if (clustrng == "hier")
  {
    if (meth == ""){
      meth = "complete" # also ward.D or single
    }
    H = hclust(dist(cbind(x,y)), method = meth)
    I = cutree(H, k=clustn)-1
  }else if (clustrng == "kmeans")
  {
    if (meth == ""){
      meth = "Lloyd"
    }
    km = kmeans(cbind(x,y), centers=clustn, nstart = 1, algorithm = meth)
    I = km$cluster-1
  }
  
  cluster_viz(x, y, I, xlab, ylab, ptnm = ptnames)
}


# Line plotting useful for MWEs
plotlines <- function(x, tvec, ...) {
  lines(tvec,x[,"lo0"],lty="dotted",...)
  lines(tvec,x[,"up0"],lty="dotted",...)
  lines(tvec,x[,"lo"],lty="dashed",lwd=1.5,...)
  lines(tvec,x[,"up"],lty="dashed",lwd=1.5,...)
  lines(tvec,x[,"mean0"],lty="solid",lwd=1.5,...)
}

plotlines_comp <- function(x,tvec,...) {
  lines(tvec,x[,"lo0"],lty="dashed",...)
  lines(tvec,x[,"up0"],lty="dashed",...)
  lines(tvec,x[,"mean0"],lty="solid",lwd=1.5,...)
}

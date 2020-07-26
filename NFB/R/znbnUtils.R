outgoners <- function(x, var, thr=2.5){
  # remove z-score outliers of a given variable at some multiple of sigma
  tz <- scale(var)
  outgoners <- x[tz<thr,]
}


chkldpkg <- function(x){
  for( i in x ){
    #  require returns TRUE invisibly if it was able to load package
    if( ! require( i , character.only = TRUE ) ){
      #  If package was not able to be loaded then re-install
      install.packages( i , dependencies = TRUE )
      #  Load package after installing
      require( i , character.only = TRUE )
    }
  }
}

detach_package <- function(pkg, character.only = FALSE)
{
  if(!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}

zdiff <- function(lm1, lm2){
  # test the slopes of two linear models for significant difference
  # the model IV should be after Intercept - 2nd in the coefficients list
  slmF <- summary(lm1)
  slmL <- summary(lm2)
  zdiff <- (slmF$coefficients[2,1] - slmL$coefficients[2,1]) / 
    sqrt((slmF$coefficients[2,1]*slmF$coefficients[2,2])^2 + 
           (slmL$coefficients[2,1]*slmL$coefficients[2,2])^2)
}


#bootstrapped CIs for differences
bootstrapCIs <- function( df, varname, runs ){
  require( boot )
  set.seed(42)
  boot.fun <- function( x, ind ){
    c("mean" = mean(x[ind]), "median" = median(x[ind]), "var" = var(x[ind]))    
  }
  
  if (nrow(df) > 1){
    b <- boot( df[[varname]], boot.fun, runs)
    bci <- boot.ci(b, type = "bca")
  }else{
    bci <- list(bca = list(conf = 1, 1, 0, df[[varname]], df[[varname]]))
  }
}

bootin <- function(df, varname){
  bCI <- bootstrapCIs(df, varname, 5000)
  data.frame(bCI$bca)
}

quadfit <- function(df, dependent, factor){

  df$quads <- fitted(lm(df[[dependent]]~factor+I(factor^2)))
  df
}

#define function that removes NA before calculating the cosine similarity
library(lsa)
cos.na <- function(x,y)
{
  cosine(na.omit(cbind(x,y)))
}

#...same for angle
angle.na <- function(x,y){
  tmp <- na.omit(cbind(x,y))
  x <- tmp[,1]
  y <- tmp[,2]
  dot.prod <- x%*%y 
  norm.x <- norm(x,type="2")
  norm.y <- norm(y,type="2")
  theta <- acos(dot.prod / (norm.x * norm.y))
  as.numeric(theta)
}


#function to average Pearson correlations, from Alexander (1990), Bull Psych Soc, 28 (4). Formula 5
#fails for small samples n<4, and works badly for rank correlations, due to large likelihood of r = +-1
agg.pearson.cors <- function(rs.ns)
{
  na.idx <- is.na(rs.ns[,1])
  rs <- rs.ns[!na.idx,1]
  ns <- rs.ns[!na.idx,2]
  k <- length(rs)
  ns.k <- sum(ns) - k
  r.curly <- rs + (rs * (1 - rs ^ 2)) / (2 * (ns - 3))
  
  aggcors <- sum((ns - 1) * r.curly) / ns.k
  return(aggcors)
}

#function computes aggregated rank correlation of multiple samples,
agg.rank.cors <- function(rs.ns)
{
  na.idx <- is.na(rs.ns[,1]) | is.na(rs.ns[,2]) # NA index
  rs <- rs.ns[!na.idx,1]
  ns <- rs.ns[!na.idx,2]
  
  aggcors <- agg.cors(rs, ns)
  
  return(aggcors)
}


# function to aggregate correlations by Hunter-Schmidt method,
# see Zhang & Wang (2014), Multivariate Behavioral Research, 49:2, 130-148
agg.cors <- function(rs, ns)
{
  rhat <- sum(rs * ns)
  
  return( rhat / sum(ns) )
}

# GM means due to here: https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
# geometric mean, handling NA and zeros
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# geometric mean, handling NA and optionally propagating zeros
gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}

center_scale <- function(df, na.rm = TRUE){
  
  if (length(df[!is.na(df)]) == 1){
    df[!is.na(df)] = 0
    return(df)
  }
  mx_df <- max(df, na.rm = na.rm)
  mn_df <- min(df, na.rm = na.rm)
  half_dist <- (mx_df - mn_df) / 2
  cs_df <- (df - (mn_df + half_dist)) / half_dist
  
  return(cs_df)
}

prettyPrintModelFit <- function(text, model){
  require(broom)
  tmp <- capture.output(glance(model))
  cat(text, "\n", tmp[1], "\n", tmp[2])
}



## ----
# Functions to generate Maximum Width Envelope analysis
# from K. Puolamäki
#
#

#' Priority queue
#' @param p order of column, i.e., output of order (x)
#' @return Lexical closure which contains the following O(1) methods:
#'  $first(), $second(), $last(), $last2(), and $remove(i) as documented
#'  in Korpela et al. (2014), appendix A.
#' @export
PQ <- function(p) {
  n <- length(p)
  idx <- 1+order(p)  # index at which item i can be found
  p <- c(NA,p,NA)      # head is at p[1] and tail at p[n+2]
  nxt <- c(2:(n+2),NA) # pointer to next
  prv <- c(NA,1:(n+1)) # pointer to previous
  first <- function() p[nxt[1]]
  second <- function() p[nxt[nxt[1]]]
  last <- function() p[prv[n+2]]
  last2 <- function() p[prv[prv[n+2]]]
  history <- rep(NA,n)
  pos <- maxpos <- 0
  remove <- function(i) {
    pos <<- maxpos <<- pos+1 # update position
    history[pos] <<- i       # add this to history
    j <- idx[i]
    prv[nxt[j]] <<- prv[j] # previous of the next is previous of the current 
    nxt[prv[j]] <<- nxt[j] # next of the previous is next of the current
    pos
  }
  ## remove/unremove previously removed items. This is not really needed
  ## but costs nothing to include here at this stage...
  goto <- function(newpos) { 
    if(newpos<pos) { # go backward
      for(i in pos:(newpos+1)) {
        j <- idx[history[i]]
        prv[nxt[j]] <<- j # unremove
        nxt[prv[j]] <<- j # unremove
      }
    } else if(newpos>pos) { # go forward
      for(i in pos:(newpos-1)) {
        j <- idx[history[i]]
        prv[nxt[j]] <<- prv[j] # remove
        nxt[prv[j]] <<- nxt[j] # remove
      }
    }
    pos <<- newpos
    pos
  }
  show <- function() list(idx=idx,p=p,nxt=nxt,prv=prv,pos=pos,
                          history=if(maxpos==0) c() else history[1:maxpos]) # for debugging
  list(first=first,second=second,last=last,last2=last2,remove=remove,goto=goto,show=show)
}

#' Finds MWE efficiently using separate validation data as documented in Korpela et al.
#' (2014). The algorithm has time complexity of O(m*n*log(n)), where m is the length of 
#' time series and n is the rows in the training data or validation data, whichever is
#' larger.
#' @param max_k integer, maximum value of k searched
#' @param max_alpha real number, fraction of validation data samples out of MWE
#' @param data_tr nXm matrix, training data
#' @param data_va n'Xm matrix, validation data
#' @return Returns a list that contains number of time series removed k, final alpha 
#'  (should be the largest alpha which is at most max_alpha), and lower and upper 
#'  bounds.
#' @export
find_mwe <- function(max_k,max_alpha,data_tr,data_va) {
  m <- dim(data_tr)[2]
  # use priority queues for both training and validation set. This finds it efficient
  # to find the largest and smallest non-removed curves.
  q_tr <- apply(data_tr,2,function(x) PQ(order(x))) 
  q_va <- apply(data_va,2,function(x) PQ(order(x)))
  lo0 <- up0 <- NULL
  k <- alpha <- -1
  removed <- 0
  while(k<max_k && removed/dim(data_va)[1]<=max_alpha) {
    idx <- matrix(c(rep(1:m,2),
                    sapply(q_tr,function(q) q$first()), sapply(q_tr,function(q) q$last()),
                    sapply(q_tr,function(q) q$second()),sapply(q_tr,function(q) q$last2())),2*m,3)
    lo <- data_tr[idx[  1:m ,c(2,1)]] # current lower envelope
    up <- data_tr[idx[-(1:m),c(2,1)]] # current upper envelope
    for(i in 1:m) { # column/time index i
      j <- q_va[[i]]$first() # index of the smallest item in validation set
      while(data_va[j,i]<lo[i] && removed/dim(data_va)[1]<=max_alpha) { 
        ## remove if below lower limit
        sapply(q_va,function(q) q$remove(j))
        j <- q_va[[i]]$first()
        removed <- removed+1
      }
      j <- q_va[[i]]$last() # index of the largest item in validation set
      while(up[i]<data_va[j,i] && removed/dim(data_va)[1]<=max_alpha) { 
        ## remove if above upper limit
        sapply(q_va,function(q) q$remove(j))
        j <- q_va[[i]]$last()
        removed <- removed+1
      }
    }
    if(removed/dim(data_va)[1]<=max_alpha) {
      lo0 <- lo
      up0 <- up
      k <- k+1
      alpha <- removed/dim(data_va)[1]
      
      ## greedy MWE algorithm: remove time series which decreases the MWE most:
      a <- aggregate(x=data.frame(gain=abs(data_tr[idx[,c(2,1)]]-data_tr[idx[,c(3,1)]])),
                     by=list(i=idx[,2]),
                     FUN=sum)
      j <- a[which.max(a[,"gain"]),"i"]
      sapply(q_tr,function(q) q$remove(j))
    }
  }
  list(k=k,alpha=alpha,lo=lo0,up=up0)
}

readdataset <- function(filename) {
  data <- read.csv(filename,header=FALSE,sep=";",na.strings="NA")[,-1]
  if(any(apply(data,2,function(x) any(is.na(x)) & any(!is.na(x))))) stop("readdataset: NAs detected")
  if(dim(data)[1]!=41) stop("readdataset: wrong dim")
  idx <- which(apply(data,2,function(x) all(!is.na(x))))
  list(data=t(as.matrix(data[,idx])),idx=idx)
}

findcurves <- function(data,mtr=5000,mva=5000,alpha=0.05,signflip=FALSE) {
  n <- dim(data)[1]
  m <- dim(data)[2]
  if(signflip) {
    samples_tr <- t(replicate(mtr,colMeans((sample(c(-1,1),size=n,replace=TRUE) %o% rep(1,m))*data)))
    samples_va <- t(replicate(mva,colMeans((sample(c(-1,1),size=n,replace=TRUE) %o% rep(1,m))*data)))
  } else {
    samples_tr <- t(replicate(mtr,colMeans(data[sample.int(n,replace=TRUE),]))) # training set
    samples_va <- t(replicate(mva,colMeans(data[sample.int(n,replace=TRUE),]))) # validation set
  } 
  a <- find_mwe(mtr-2,alpha,samples_tr,samples_va)
  q <- apply(samples_tr,2,function(x) quantile(x,probs=c(alpha/2,1-alpha/2)))
  data.frame(mean0=colMeans(data),lo=a$lo,up=a$up,k=a$k,alpha=a$alpha,lo0=q[1,],up0=q[2,])
}

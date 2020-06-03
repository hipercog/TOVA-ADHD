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
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
  b <- boot( df[[varname]], boot.fun, runs)
  bci <- boot.ci(b, type = "bca")
}

bootin <- function(df, varname){
  bCI <- bootstrapCIs(df, varname, 5000)
  data.frame(bCI$bca)
}

quadfit <- function(df, dependent, factor){

  df$quads <- fitted(lm(df[[dependent]]~factor+I(factor^2)))
  df
}
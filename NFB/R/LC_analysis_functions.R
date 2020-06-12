## ----
# Functions to assist LC analysis of CENT data
# B. Cowley
#
source('R/znbnUtils.R')

testOnes <- function(df, ttl = deparse(substitute(df)))
{
  tmp <- df %>% group_by(patient,session_num) %>% count()
  # hist(tmp[tmp$n == 1,]$session_num, breaks = 1:41)
  hist(tmp[tmp$n == 1,]$patient, breaks = 1001:1081, main = ttl)
}

# Print summaries for a list of data.frames of trial data
reportTrials <- function(df, ...)
{
  if (inherits(df, "list"))
  {
    L2 <- lapply(df, function(x) as.data.frame(x, stringsAsFactors = FALSE))
    lapply(seq_along(L2), function(i){ printTrials(L2[[i]], names(L2)[i]) } )
    # for (idx in 1:length(df))
    # {
    #   printTrials(df[idx], ...)
    # }
  }else if (inherits(df, "data.frame")){
    printTrials(df, ...)
  }
}

printTrials <- function(df, colonom, group1 = NA, group2 = NA)
{
  print(paste0("Data: ", deparse(substitute(df)), "$", colonom))
  df <- as.data.frame(df)
  print(summary(df[,colonom]))
  
  if (!is.na(group1) && !is.na(group2))
  {
    print(paste("... subset: ", deparse(substitute(group1))))
    print(summary(subset(df, patient%in%group1)[,colonom]))
    print(paste("... subset: ", deparse(substitute(group2))))
    print(summary(subset(df, patient%in%group2)[,colonom]))
  }
}

## ## ## ## ## --------------------- DERIVE FEATURES FUNCTIONS ----
#
getSSN <- function(df, DV, summary_method = "median")
{
  f <- match.fun(summary_method)
  ssn <- setDT(df)[, .(f(get(DV))), by = .(patient, nom_session)]
  setnames(ssn, names(ssn)[3], DV)
  
  return(ssn)
}

getSsnMLC <- function(df, DV, precomp = FALSE, ...)
{
  if (precomp){
    ssn <- df
  }else{
    ssn <- getSSN(df, DV, ...)
  }
  
  ssn.MLC <- setDT(ssn)[
    , .(cors = cor(get(DV), session, use="pair", method="k"))
    , by = patient]$cors
  
  ssn.count <- ssn %>% group_by(patient) %>% summarise(ssn.count = max(session))
  
  return(cbind(ssn.count[,1], ssn.MLC, ssn.count[,2]))
}

getTrialAggCor <- function(mtx, corx, cory, id = TRUE)
{
  require(abind)
  trial.cors <- getTrialCors(mtx, corx, cory, id)
  trial.count <- getTrialCount(mtx, FALSE)
  ssn.count <- apply(!is.na(trial.count), 1, sum)
  
  if (id)
  {
    id <- trial.cors[,1]
    trial.cors <- trial.cors[,-(1)]
  }
  agg.cor <- apply(abind(trial.cors, trial.count, rev.along = 0), 1, agg.rank.cors)
  
  return(as.data.frame(cbind(id, agg.cor, ssn.count)))
}

aggMLC <- function(tr, sn, fill = FALSE)
{
  if (nrow(tr) < nrow(sn))
  {
    ix <- sn[,1] %in% tr[,1]
    aggMLC <- apply(abind(cbind(tr[,2], sn[ix,2]), cbind(tr[,3], sn[ix,3]), rev.along = 0)
                    , 1, agg.rank.cors)
    aggMLC <- as.data.frame(cbind(tr[,1], aggMLC))
    if (fill)
    {
      names(aggMLC) <- names(sn)[1:2]
      aggMLC <- rbind( aggMLC, sn[!ix,1:2] )
      aggMLC <- aggMLC[order(aggMLC[,1]),]
    }
  }else{
    aggMLC <- apply(abind(cbind(tr[,2], sn[,2]), cbind(tr[,3], sn[,3]), rev.along = 0)
                    , 1, agg.rank.cors)
    aggMLC <- cbind(tr[,1], aggMLC)
  }
  return(as.data.frame(aggMLC))
}

getTrialCosSim <- function(trial.cors, labels, AHLC = "Fitts", phase123 = c(0, 1, 0.5))
{
  # BY TRIALS: cosine similarity of trial cors to an arbitrary hypothetical LC (AHLC)
  N <- ncol(trial.cors)
  # AHLC by linear increment from 0 to 1
  if (AHLC == "Fitts")
  {
    # phase123 <- c(0, 1, 0.5)
    one3rd <- ceiling(N / 3)
    arhyLC <- c(rep(phase123[1], one3rd), rep(phase123[2], one3rd), rep(phase123[3], one3rd))[2:(N+1)]
    min.ssn <- 3
  }else if (AHLC == "increment")
  {
    arhyLC <- seq(0, 1, 1 / (N - 1))
    min.ssn <- 2
  }else if (AHLC == "flat")
  {
    arhyLC <- rep(1, N)
    min.ssn <- 1
  }
  
  idx <- apply(!is.na(trial.cors), 1, sum) >= min.ssn
  trial.MLC <- apply(trial.cors[idx, ], 1, cos.na, y = arhyLC)
  
  return(cbind(labels[idx], trial.MLC[2,]))
}

getTrialCors <- function(df, X, Y, id = TRUE)
{
  require(data.table)
  tr.cors <- as.matrix(spread(
    setDT(df)[
      ,.(cors = cor(round(x = get(X)), round(y = get(Y)), use="pair", method="k"))
      , by = .(patient, nom_session)]
      , nom_session, `cors`))
  if (!id)
  {
    tr.cors <- tr.cors[,-(1)]
  }
  return(tr.cors)
}

getTrialCSGmean <- function(df, X, id = TRUE, cs = '-1to1')
{
  require(data.table)
  tr.gmean <- as.matrix(spread(
    setDT(df)[
      , .(geom_mean = gm_mean(x = get(X)))
      , by = .(patient, nom_session)]
      , nom_session, `geom_mean`))
  # center and scale
  tmp <- tr.gmean[,-(1)]
  if (cs == '-1to1')
  {
    tmp <- t(apply(tmp, 1, center_scale)) #center_scale is from znbnUtils.R
  }else if (cs == 'SD'){
    tmp <- apply(tmp, 1, scale)
  }
  tr.gmean[,-(1)] <- tmp
  # apply column names & remove patient ID column if requested
  # rownames(tr.gmean) <- c('patient', X)
  if (!id)
  {
    subset(tr.gmean, select = -c('patient'))
  }
  
  return(tr.gmean)
}

getTrialCount <- function(df, id = TRUE)
{
  tr.cnt <- as.matrix(spread(
    df %>% group_by(patient, nom_session) %>% count()
    , nom_session, `n`))
  if (!id)
  {
    tr.cnt <- tr.cnt[,-(1)]
  }
  return(tr.cnt)
}



## ## ## ## ## --------------------- EXTRACT DATA FUNCTIONS ----
# Trial mean & CI
getTrialsMCI <- function(x, byvar, jcol)
{
  require(data.table)
  fM <- setDT(x)[, lapply(.SD, mean), by=byvar, .SDcols=jcol]
  fM <- fM[order(fM[,1]),]
  fCI <- setDT(x)[, bootin(.SD, varname = jcol), by=byvar]
  fCI <- fCI[order(fCI[,1]),]
  list(M=fM, CI=fCI)
}

# Remove any sessions containing trials <= min_trial
cut1TrialSsn <- function(df, min_trial = 1)
{
  nr <- nrow(df)
  if (nr <= min_trial)
  {
    df <- df[0,]
  }
  return(df)
}

# IF subject sessions are less than canon_s :> pad by last-observation-carried-forward
# IF subject sessions greater than canon_s :> trim away earliest sessions 
padTrimOnce <- function(df, canon_s)
{
  num_ssn <- length(unique(df$session_num))
  if (num_ssn < canon_s)
  {
    tile <- df[df$session_num == max(df$session_num),]
    for (ridx in 1:(canon_s - num_ssn))
    {
      tile$session_num <- tile$session_num + 1
      df <- rbind(df, tile)
    }
    df$session_num <- df$session_num - (canon_s - num_ssn)
  }else if(num_ssn > canon_s)
  {
    df <- df[!(df$session_num %in% unique(df$session_num)[1:(num_ssn - canon_s)]),]
  }

  return(df)
}

tidyTrials <- function(df, canon_s = 0, min_trial = 0)
{
  # remove any sessions with trials < min_trial - because e.g. cannot calculate correlations from 1 point
  if (min_trial > 0)
  {
    df <- df %>% group_by(patient,session_num) %>% do(cut1TrialSsn(., min_trial))
  }
  if (canon_s > 0)
  {
    # Trim, or pad by repeating last session, to match required number of sessions
    df <- df %>% group_by(patient) %>% do(padTrimOnce(., canon_s))
    # create an index of session numbers that aligns by the LAST completed session
    df <- with(df, df %>% group_by(patient) %>%
                 mutate(nom_session = session_num - (max(session_num) - canon_s)))
  }else{
    # create an index of session numbers that goes from 1:number_of_sessions
    df <- with(df, df %>% group_by(patient) %>% 
                 mutate(nom_session = as.numeric(as.factor(rank(session_num))) ))
  }
  return(df)
}

# All trials
getAllTrials <- function(df, ...)
{
  df <- tidyTrials(df, ...)
  getAllTrials <- droplevels(subset(df, nom_session > 0))
}

### Normal trials ----
getNorTrials <- function(df, ...)
{
  norTrials <- droplevels(subset(df, norm_notInv==1 & trialtype==0))
  norTrials <- tidyTrials(norTrials, ...)
  return(norTrials)
}

### Inverse trials ----
getInvTrials <- function(df, ...)
{
  invTrials <- droplevels(subset(df, norm_notInv==0 & trialtype==0))
  invTrials <- tidyTrials(invTrials, ...)
  return(invTrials)
}
# Transfer trials
getTraTrials <- function(df, ...)
{
  traTrials <- droplevels(subset(df, norm_notInv==1 & trialtype>0))
  traTrials <- tidyTrials(traTrials, ...)
  return(traTrials)
}

### Not-inverse trials ----
getNotInvTrials <- function(df, ...)
{
  notiTrials <- droplevels(subset(df, norm_notInv==1))
  notiTrials <- tidyTrials(notiTrials, ...)
  return(notiTrials)
}

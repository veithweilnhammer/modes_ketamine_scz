
fitMetaDprime   <- function (data, model = "ML", nInits = 5, nRestart = 3, .parallel = FALSE, 
          n.cores = NULL) 
{
  if (!model %in% c("ML", "F")) {
    stop("model must be either 'ML' or 'F'")
  }
  if (length(unique(data$stimulus)) != 2) {
    stop("There must be exacltly two different possible values of stimulus")
  }
  if (!is.factor(data$stimulus)) {
    data$stimulus <- factor(data$stimulus)
    warning("stimulus is transformed to a factor!")
  }
  if (!is.factor(data$rating)) {
    data$rating <- factor(data$rating)
    warning("rating is transformed to a factor!")
  }
  if (!all(data$correct %in% c(0, 1))) 
    stop("correct should be 1 or 0")
  if (!any(data$correct == 0)) 
    stop("There should be at least one erroneous response")
  if (!any(data$correct == 1)) 
    stop("There should be at least one correct response")
  nRatings <- length(unique(data$rating))
  outnames <- c("model", "participant", "dprime", "c", "metaD", 
                "Ratio")
  call_fitfct <- function(X) {
    cur_model <- model[X[1]]
    cur_sbj <- X[2]
    participant <- NULL
    data_part <- subset(data, participant == cur_sbj)
    res <- int_fitMetaDprime(ratings = data_part$rating, 
                             stimulus = data_part$stimulus, correct = data_part$correct, 
                             ModelVersion = cur_model, nInits = nInits, nRestart = nRestart)
    res$model <- cur_model
    res$participant <- cur_sbj
    res[outnames[!(outnames %in% names(res))]] <- NA
    res <- res[, outnames]
    return(res)
  }
  subjects <- unique(data$participant)
  nJobs <- length(model) * length(subjects)
  jobs <- expand.grid(model = 1:length(model), sbj = subjects)
  if (.parallel) {
    listjobs <- list()
    for (i in 1:nrow(jobs)) {
      listjobs[[i]] <- c(model = jobs[["model"]][i], sbj = jobs[["sbj"]][i])
    }
    if (is.null(n.cores)) 
      n.cores <- detectCores() - 1
    cl <- makeCluster(type = "SOCK", n.cores)
    clusterExport(cl, c("data", "model", "outnames", "call_fitfct", 
                        "nInits", "nRestart"), envir = environment())
    on.exit(try(stopCluster(cl), silent = TRUE))
    res <- clusterApplyLB(cl, listjobs, fun = call_fitfct)
    stopCluster(cl)
  }
  else {
    res <- apply(X = jobs, 1, FUN = call_fitfct)
  }
  res <- do.call(rbind, res)
  res <- res[, apply(res, 2, function(X) any(!is.na(X)))]
  return(res)
}





int_fitMetaDprime   <- function(ratings, stimulus, correct,
                                ModelVersion = "ML",
                                nInits = 5, nRestart = 3){
  
  # prepare data
  
  nRatings <-  length(levels(ratings))
  nCriteria <- nRatings * 2 - 1
  abj_f <- 1 /(nRatings*2)
  abs_corrects <-  table(ratings[correct == 1], stimulus[correct == 1]) + abj_f
  abs_errors <- table(ratings[correct == 0], stimulus[correct == 0]) + abj_f
  
  abs_S1 <- c(rev(abs_errors[,2]),abs_corrects[,2])
  ratingHrs <- qnorm(1 - cumsum(abs_S1)/sum(abs_S1))
  abs_S2 <-  c(rev(abs_corrects[,1]), abs_errors[,1] )
  ratingFrs <-  qnorm(1 - cumsum(abs_S2)/sum(abs_S2))
  finits <- is.finite(ratingHrs) & is.finite(ratingFrs)
  ratingHrs <- as.vector(ratingHrs[finits])
  ratingFrs <- as.vector(ratingFrs[finits])
  
  nC_rS1 <- rev(as.vector(abs_corrects[,1]))
  nI_rS1 <- rev(as.vector(abs_errors[,2]))
  nC_rS2 <- as.vector(abs_corrects[,2])
  nI_rS2 <- as.vector(abs_errors[,1])
  
  # compute type 1 parameters based on formulae
  
  dprime <- ratingHrs[nRatings] - ratingFrs[nRatings]
  cs <- (-.5 * ( ratingHrs  + ratingFrs))
  cprime <- cs[nRatings]/dprime
  
  # use a coarse grid search to identify the most promising starting parameters
  temp <- expand.grid(d = seq(0, 5, length.out = 10),
                      tauMin =  seq(.1,2,length.out=10),  # position of the most conservative confidence criterion related to stimulus A
                      tauRange = seq(0.5,5,length.out=10))  # range of rating criteria stimulus B  #  position of the most liberal confidence criterion with respect to thet
  
  inits <- data.frame(matrix(data=NA, nrow= nrow(temp), ncol = 1 + (nRatings-1)*2))
  inits[,1] <- temp$d #  qnorm((temp$d + 10)/ 20)
  inits[,2:(nRatings-1)] <-
    log(t(mapply(function(tauRange) rep(tauRange/(nRatings-1), nRatings-1),
                 temp$tauRange)))
  inits[,nRatings:(nRatings+1)] <- rep(log(temp$tauMin),2)
  inits[,(nRatings+2):(nRatings*2-1)] <-
    log(t(mapply(function(tauRange) rep(tauRange/(nRatings-1), nRatings-1),
                 temp$tauRange)))
  
  if(ModelVersion == "ML"){
    logL <- apply(inits, MARGIN = 1,
                  function(p) try(negLoglMetaD(p, nC_rS1,nI_rS1, nC_rS2,nI_rS2,nRatings, cprime), silent = TRUE))
  }
  if (ModelVersion == "F"){
    logL <- apply(inits, MARGIN = 1,
                  function(p) try(negLoglFleming(p, nC_rS1,nI_rS1, nC_rS2,nI_rS2,nRatings, cs[nRatings]), silent = TRUE))
  }
  
  logL <- as.numeric(logL)
  inits <- inits[order(logL),]
  noFitYet <- TRUE
  
  if(ModelVersion == "ML"){
    for (j in 1:nInits){
      m <- try(optim(par = inits[j,], fn = negLoglMetaD, gr = NULL,
                     nC_rS1 = nC_rS1, nI_rS1 = nI_rS1, nC_rS2 = nC_rS2, nI_rS2 = nI_rS2,
                     nRatings = nRatings, cprime = cprime,
                     control = list(maxit = 10^6, reltol = 10^-8)), silent=T)
      for(i in 2:nRestart){
        try(m <- try(optim(par = m$par, fn = negLoglMetaD, gr = NULL,
                           nC_rS1 = nC_rS1, nI_rS1 = nI_rS1, nC_rS2 = nC_rS2, nI_rS2 = nI_rS2,
                           nRatings = nRatings, cprime = cprime,
                           control = list(maxit = 10^6, reltol = 10^-8)), silent=T),
            silent=T)
      }
      if ((!inherits(m, "try-error"))){
        if (noFitYet) {
          fit <- m
          noFitYet <- FALSE
        } else {
          if (m$value < fit$value) fit <- m
        }
      }
    }
  }
  
  if(ModelVersion == "F"){
    for (j in 1:nInits){
      m <- try(optim(par = inits[j,], fn = negLoglFleming, gr = NULL,
                     nC_rS1 = nC_rS1, nI_rS1 = nI_rS1, nC_rS2 = nC_rS2, nI_rS2 = nI_rS2,
                     nRatings = nRatings, type1_c = cs[nRatings],
                     control = list(maxit = 10^6, reltol = 10^-8)), silent=T)
      for(i in 2:nRestart){
        try(m <- try(optim(par = m$par, fn = negLoglFleming, gr = NULL,
                           nC_rS1 = nC_rS1, nI_rS1 = nI_rS1, nC_rS2 = nC_rS2, nI_rS2 = nI_rS2,
                           nRatings = nRatings, type1_c = cs[nRatings],
                           control = list(maxit = 10^6, reltol = 10^-8)), silent=T),
            silent=T)
      }
      if ((!inherits(m, "try-error"))){
        if (noFitYet) {
          fit <- m
          noFitYet <- FALSE
        } else {
          if (m$value < fit$value) fit <- m
        }
      }
    }
  }
  
  result <- data.frame(
    dprime = dprime,
    cprime = cprime,
    c = dprime * cprime,
    metaD = NA, Ratio = NA, ModelVersion = ModelVersion)
  
  if(exists("fit")){
    if(is.list(fit)){
      result$metaD = fit$par[1]
      result$Ratio =  result$metaD / dprime
      
    }
  }
  
  result
}

negLoglMetaD <- function(parameters, nC_rS1,nI_rS1, nC_rS2,nI_rS2,nRatings, cprime){
  metadprime <- parameters[1] # pnorm(parameters[1])*20 - 10
  S1mu <- -metadprime/2
  S2mu <- metadprime/2
  meta_c <- metadprime*cprime
  t2_rS1 <- c(-Inf, meta_c - rev(cumsum(exp(parameters[2:nRatings]))), meta_c)
  t2_rS2 <- c(meta_c, meta_c + cumsum(exp(parameters[(nRatings+1):length(parameters)])), Inf)
  
  prC_rS1 <- (pnorm(t2_rS1[2:(nRatings+1)],S1mu) -
                pnorm(t2_rS1[1:nRatings],S1mu)) /pnorm(meta_c,S1mu)
  prI_rS1 <- (pnorm(t2_rS1[2:(nRatings+1)], S2mu) -
                pnorm(t2_rS1[1:nRatings],S2mu) ) / pnorm(meta_c,S2mu)
  prC_rS2 <- ((1- pnorm(t2_rS2[1:nRatings],S2mu)) -
                (1- pnorm(t2_rS2[2:(nRatings+1)], S2mu))) /
    (1 - pnorm(meta_c, S2mu))
  prI_rS2 <- ((1- pnorm(t2_rS2[1:nRatings],S1mu)) -
                (1- pnorm(t2_rS2[2:(nRatings+1)],S1mu))) /
    (1 - pnorm(meta_c, S1mu))
  
  logL <-  -sum(nC_rS1*log(prC_rS1),nI_rS1*log(prI_rS1),
                nC_rS2*log(prC_rS2),nI_rS2*log(prI_rS2))
  prC_rS1[is.na( prC_rS1) | is.nan(prC_rS1) | prC_rS1 < 10^-10] <- 10^-10
  prI_rS1[is.na(prI_rS1) | is.nan(prI_rS1) | prI_rS1 < 10^-10] <- 10^-10
  prC_rS2[is.na( prC_rS2) | is.nan( prC_rS2) |  prC_rS2 < 10^-10] <- 10^-10
  prI_rS2[is.na(prI_rS2) | is.nan(prI_rS2) | prI_rS2 < 10^-10] <- 10^-10
  logL
}

negLoglFleming <- function(parameters, nC_rS1,nI_rS1, nC_rS2,nI_rS2,nRatings, type1_c){
  metadprime <- parameters[1] # pnorm(parameters[1])*20 - 10
  S1mu <- -metadprime/2
  S2mu <- metadprime/2
  t2_rS1 <- c(-Inf, type1_c - rev(cumsum(exp(parameters[2:nRatings]))), type1_c)
  t2_rS2 <- c(type1_c, type1_c + cumsum(exp(parameters[(nRatings+1):length(parameters)])), Inf)
  
  prC_rS1 <- (pnorm(t2_rS1[2:(nRatings+1)],S1mu) -
                pnorm(t2_rS1[1:nRatings],S1mu)) /pnorm(type1_c,S1mu)
  prI_rS1 <- (pnorm(t2_rS1[2:(nRatings+1)], S2mu) -
                pnorm(t2_rS1[1:nRatings],S2mu) ) / pnorm(type1_c,S2mu)
  prC_rS2 <- ((1- pnorm(t2_rS2[1:nRatings],S2mu)) -
                (1- pnorm(t2_rS2[2:(nRatings+1)], S2mu))) /
    (1 - pnorm(type1_c, S2mu))
  prI_rS2 <- ((1- pnorm(t2_rS2[1:nRatings],S1mu)) -
                (1- pnorm(t2_rS2[2:(nRatings+1)],S1mu))) /
    (1 - pnorm(type1_c, S1mu))
  
  logL <-  -sum(nC_rS1*log(prC_rS1),nI_rS1*log(prI_rS1),
                nC_rS2*log(prC_rS2),nI_rS2*log(prI_rS2))
  prC_rS1[is.na( prC_rS1) | is.nan(prC_rS1) | prC_rS1 < 10^-10] <- 10^-10
  prI_rS1[is.na(prI_rS1) | is.nan(prI_rS1) | prI_rS1 < 10^-10] <- 10^-10
  prC_rS2[is.na( prC_rS2) | is.nan( prC_rS2) |  prC_rS2 < 10^-10] <- 10^-10
  prI_rS2[is.na(prI_rS2) | is.nan(prI_rS2) | prI_rS2 < 10^-10] <- 10^-10
  logL
}
bonferroni <- function(C) {
  sign = which(C$coefficients[,ncol(C$coefficients)] * nrow(C$coefficients) < 0.05)
  C$c_sign_coefficients = C$coefficients[sign,]
  
  C$c_coefficients = C$coefficients
  C$coefficients[,ncol(C$coefficients)] = C$coefficients[,ncol(C$coefficients)] * nrow(C$coefficients)
  C$coefficients[C$coefficients[,ncol(C$coefficients)] > 1,ncol(C$coefficients)] = 1
  return(C)
}

checkBinaryTrait = function(v, naVal = "NA") {
  if (!is.numeric(v))
    stop("Only numeric vectors are accepted.")
  vSet = unique(v)
  if (!missing(naVal))
    vSet[vSet == naVal] = NA
  vSet = vSet[!is.na(vSet)]
  if (any(as.integer(vSet) != vSet))
    "con"
  else if (length(vSet) > 2)
    "cat"
  else
    "bin"
}


frequency.spectrum <- function(X.k, xlimits=c(0,length(X.k))) {
  cbind(0:(length(X.k)-1), Mod(X.k))
}

compute_AIC <- function(LL, K, n, correction) {
AIC = -2*LL + 2*K + correction * ((2*K*(K+1)/(n-K-1)))
return(AIC)
}



exclude_3SD = function(x) {
  v = x
  average = mean(v, na.rm = TRUE)
  standard_dev = sd(v, na.rm = TRUE)
  
  v[v > average + 3*standard_dev] = NA
  v[v < average - 3*standard_dev] = NA
  return(v)
}


getmode <- function(v) {
  uniqv <- unique(v, na.rm = TRUE)
  uniqv[which.max(tabulate(match(v, uniqv)))]
} 

#' Type 1 SDT for a 2AFC design
type_1_sdt <- function(df, stimulus = NULL, response = NULL,
                       counts = total, s = 1, add_constant = TRUE) {
  
  stim_col <- dplyr::enquo(stimulus)
  count_col <- dplyr::enquo(counts)
  resp_col <- dplyr::enquo(response)
  
  df <- dplyr::group_by(df, !!stim_col)
  
  if (add_constant) {
    df <- dplyr::mutate(df, proportions = ( (!!count_col) + 1 /
                                              n()) / (sum( (!!count_col)) + 1))
  } else{
    df <- dplyr::mutate(df, proportions = (!!count_col) / sum( (!!count_col)))
  }
  
  reps_only <- dplyr::filter(df, (!!resp_col) == 1)
  
  s1_HR <- reps_only$proportions[[1]]
  s1_FA <- reps_only$proportions[[2]]
  
  d_prime <- (1 / s) * qnorm(s1_HR) - qnorm(s1_FA)
  c_raw <- (-1 / (1 + s)) * (qnorm(s1_HR) + qnorm(s1_FA))
  c_prime <- c_raw / d_prime
  data.frame(d_prime, c_raw, c_prime, s1_HR, s1_FA)
}


fit_meta_d_SSE <- function(nR_S1, nR_S2, s = 1, d_min = -5, d_max = 5,
                           d_grain = .01, add_constant = TRUE) {
  
  if (add_constant) {
    nR_S1 <- nR_S1 + (1 / length(nR_S1))
    nR_S2 <- nR_S2 + (1 / length(nR_S2))
  }
  
  n_ratings <- length(nR_S1) / 2
  n_criteria <- 2 * n_ratings - 1
  
  S1_HR <- sum(nR_S2[(n_ratings + 1):(n_ratings * 2)]) / sum(nR_S2)
  S1_FA <- sum(nR_S1[(n_ratings + 1):(n_ratings * 2)]) / sum(nR_S1)
  
  d_1 <- (1 / s) * qnorm(S1_HR) - qnorm(S1_FA)
  c_1 <- (-1 / (1 + s)) * (qnorm(S1_HR) + qnorm(S1_FA))
  c_prime <- c_1 / d_1
  
  obs_HR2_rS1 <- NULL
  obs_FAR2_rS1 <- NULL
  
  for (i in 1:(n_ratings - 1)) {
    obs_HR2_rS1[i] <- sum(nR_S1[1:i]) / sum(nR_S1[1:n_ratings])
    obs_FAR2_rS1[i] <- sum(nR_S2[1:i]) / sum(nR_S2[1:n_ratings])
  }
  
  obs_HR2_rS2 <- NULL
  obs_FAR2_rS2 <- NULL
  for (i in (n_ratings + 2):(2 * n_ratings)) {
    obs_HR2_rS2[(i - (n_ratings + 2) + 1)] <-
      sum(nR_S2[i:(n_ratings * 2)]) /
      sum(nR_S2[(n_ratings + 1):(n_ratings * 2)])
    
    obs_FAR2_rS2[(i - (n_ratings + 2) + 1)] <-
      sum(nR_S1[i:(n_ratings * 2)]) /
      sum(nR_S1[(n_ratings + 1):(n_ratings * 2)])
  }
  
  d_grid <- seq(d_min, d_max, by = d_grain)
  c_grid <- c_prime * d_grid
  
  S1mu <- -d_grid / 2
  S2mu <- d_grid / 2
  S1sd <- 1
  S2sd <- 1 / s
  
  bounds <- 5 * max(S1sd, S2sd)
  SSEmin <- Inf
  
  param_space <- map2(S1mu, S2mu,
                      ~seq(.x - bounds, .y + bounds, by = .001)) # param space
  
  min_z <- map2(param_space, c_grid, ~min(abs(.x - .y)))
  c_ind <- map2(param_space, c_grid, ~which.min(abs(.x - .y)))
  
  HRs <- map2(param_space, S2mu,
              ~ 1 - (pnorm(.x, .y, S2sd)))
  FARs <- map2(param_space, S1mu,
               ~ 1 - (pnorm(.x, .y, S1sd)))
  
  # fit type 2 data for S1 responses
  est_HR2s_rS1 <- map2(FARs, c_ind, ~(1 - .x[1:.y]) / (1 - .x[.y]))
  est_FAR2s_rS1 <- map2(HRs, c_ind, ~(1 - .x[1:.y]) / (1 - .x[.y]))
  
  SSE <- NULL
  SSE_rS1 <- matrix(data = NaN, nrow = length(d_grid), ncol = n_ratings - 1)
  rS1_ind <- matrix(data = NaN, nrow = length(d_grid), ncol = n_ratings - 1)
  
  for (n in (1:(n_ratings - 1))) {
    SSE <- map2(est_HR2s_rS1, est_FAR2s_rS1,
                ~ (.x - obs_HR2_rS1[[n]]) ^ 2 + (.y - obs_FAR2_rS1[[n]]) ^ 2)
    SSE_rS1[, n] <- map_dbl(SSE, min)
    inds <- unlist(map(SSE, which.min))
    rS1_ind[1:length(inds), n] <- inds
  }
  
  # fit type 2 data for S2 responses
  est_HR2s_rS2 <- map2(HRs, c_ind, ~ .x[.y:length(.x)] / .x[.y])
  est_FAR2s_rS2 <- map2(FARs, c_ind, ~ .x[.y:length(.x)] / .x[.y])
  
  SSE_rS2 <- matrix(data = NaN, nrow = length(d_grid), ncol = n_ratings - 1)
  rS2_ind <- matrix(data = NaN, nrow = length(d_grid), ncol = n_ratings - 1)
  for (n in (1:(n_ratings - 1))) {
    SSE <- map2(est_HR2s_rS2, est_FAR2s_rS2,
                ~(.x - obs_HR2_rS2[[n]]) ^ 2 + (.y - obs_FAR2_rS2[[n]]) ^ 2)
    SSE_rS2[, n] <- map_dbl(SSE, min)
    inds <- unlist(map(SSE, which.min))
    rS2_ind[1:length(inds), n] <- inds
  }
  
  # update analysis
  SSEtot <- rowSums(SSE_rS1) + rowSums(SSE_rS2)
  SSEmin <- Inf
  
  meta_d <- vector("numeric", 1)
  meta_c <- vector("numeric", 1)
  t2c_rS1 <- NULL
  t2c_rS2 <- NULL
  
  min_SSE <- which.min(SSEtot)
  meta_d <- d_grid[[min_SSE]]
  meta_c <- c_grid[[min_SSE]]
  t2c_rS1 <- param_space[[min_SSE]][rS1_ind[[min_SSE]]]
  t2c_rS2 <- param_space[[min_SSE]][c_ind[[min_SSE]] + rS2_ind[[min_SSE]] - 1]
  est_HR2_rS1  <- est_HR2s_rS1[[min_SSE]][rS1_ind[[min_SSE]]]
  est_FAR2_rS1 <- est_FAR2s_rS1[[min_SSE]][rS1_ind[[min_SSE]]]
  est_HR2_rS2  <- est_HR2s_rS2[[min_SSE]][rS2_ind[[min_SSE]]]
  est_FAR2_rS2 <- est_FAR2s_rS2[[min_SSE]][rS2_ind[[min_SSE]]]
  
  out <- data.frame(da = d_1,
                    meta_da = meta_d,
                    M_diff = meta_d - d_1,
                    M_ratio = meta_d / d_1,
                    meta_ca = meta_c,
                    s = s,
                    t2ca_rS1 = t2c_rS1,
                    t2ca_rS2 = t2c_rS2,
                    SSE = SSEtot[[min_SSE]],
                    est_HR2_rS1 = est_HR2_rS1,
                    obs_HR2_rS1 = obs_HR2_rS1,
                    est_HR2_rS2 = est_HR2_rS2,
                    obs_HR2_rS2 = obs_HR2_rS2,
                    est_FAR2_rS1 = est_FAR2_rS1,
                    obs_FAR2_rS1 = obs_FAR2_rS1,
                    est_FAR2_rS2 = est_FAR2_rS2,
                    obs_FAR2_rS2 = obs_FAR2_rS2)
  return(out)
}

#' Internal function for fitting meta_d_plus
#' @param x0 Starting parameters.
#' @param th Theta.
#' @param hp Hit rate for positive responses
#' @param fp False positive rate for positive responses.

fit_metad_plus <- function(x0, th, hp, fp) {
  y1 <- (1 - pnorm(x0[[1]], x0[[2]], 1)) /
    (1 - pnorm(th * x0[[2]], x0[[2]], 1)) - hp
  y2 <- (1 - pnorm(x0[[1]], 0, 1)) / (1 - pnorm(th * x0[[2]], 0, 1)) - fp
  c(y1, y2)
}

#' Internal function for fitting meta_d_minus
#' @param x0 Starting parameters.
#' @param th Theta
#' @param hm Hit rate for negative responses
#' @param fm False positive rate for negative responses

fit_metad_minus <- function(x0, th, hm, fm) {
  y1 <- pnorm(x0[[1]], 0, 1) / pnorm(th * x0[[2]], 0, 1) - hm
  y2 <- pnorm(x0[[1]], x0[[2]], 1) / pnorm(th * x0[[2]], x0[[2]], 1) - fm
  c(y1, y2)
}

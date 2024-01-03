#' @title Estimating derivative of effects of incremental propensity score
#' interventions
#'
#' @description \code{ipsi_deriv} is used to estimate the derivative of effects 
#' of incremental propensity score interventions, i.e., estimates of the mean 
#' outcomes if the odds of receiving treatment were increased infitesimally
#' from those corresponding to incremental parameter delta (which are themselves
#' the odds of treatment if the natural odds of receiving treatment were 
#' multiplied by delta).
#'
#' @usage ipsi_deriv(y, a, x.trt, x.out, id, delta.seq, nsplits, 
#'  ci_level = 0.95, progress_bar = TRUE, return_ifvals = FALSE, fit,
#'  sl.lib=c("SL.earth","SL.gam","SL.glm","SL.glmnet","SL.glm.interaction", 
#'  "SL.mean","SL.ranger","rpart"))
#'
#' @param y Outcome of interest
#' @param a Binary treatment.
#' @param x.trt Covariate matrix for treatment regression.
#' @param x.out Covariate matrix for outcome regression.
#' @param id Subject identifier.
#' @param delta.seq Sequence of delta increment values for incremental
#'  propensity score intervention.
#' @param nsplits Integer number of sample splits for nuisance estimation. If
#'  \code{nsplits = 1}, sample splitting is not used, and nuisance functions are
#'  estimated n full sample (in which case validity of standard errors and
#'  confidence intervals requires empirical process conditions). Otherwise must
#'  have \code{nsplits > 1}.
#' @param ci_level A \code{numeric} value giving the level (1 - alpha) of the
#'  confidence interval to be computed around the point estimate.
#' @param progress_bar A \code{logical} value indicating whether to print a
#'  customized progress bar as various stages of computation reach completion.
#'  The default is \code{TRUE}, printing a progress bar to inform the user.
#' @param return_ifvals A \code{logical} indicating whether the estimated
#'  observation-level values of the influence function ought to be returned as
#'  part of the output object. The default is \code{FALSE} as these values are
#'  rarely of interest in standard usage.
#' @param fit How nuisance functions should be estimated. Options are "rf" for
#'  random forests via the \code{ranger} package, or "sl" for super learner.
#' @param sl.lib sl.lib algorithm library for SuperLearner.
#' Default library includes "earth", "gam", "glm", "glmnet", "glm.interaction",
#' "mean", "ranger", "rpart.
#'
#' @section Details:
#' Treatment and covariates are expected to be one timepoint. Therefore if
#' \code{n}, then \code{a}, \code{y}, and \code{id} should all be vectors of 
#' length \code{n}, and \code{x.trt} and \code{x.out} should be matrices with 
#' \code{n} rows. The subject ordering should be consistent  across function 
#' inputs, based on the ordering specified by \code{id}. See example below for
#' an illustration.
#'
#' @return A list containing the following components:
#' \item{res}{ estimates/SEs and uniform CIs for population means.}
#' \item{res.ptwise}{ estimates/SEs and pointwise CIs for population means.}
#' \item{calpha}{ multiplier bootstrap critical value.}
#' \item{ifvals}{ the influence function values is optional.}
#'
#' @importFrom stats qnorm as.formula
#' @importFrom ranger ranger
#'
#' @export
#'
#' @examples
#' n <- 500
#'
#' id <- rep(1:n)
#' x.trt <- matrix(rnorm(n * 5), nrow = n)
#' x.out <- matrix(rnorm(n * 5), nrow = n)
#' a <- rbinom(n , 1, .5)
#' y <- rnorm(mean = 1, n)
#'
#' d.seq <- seq(0.1, 5, length.out = 10)
#'
#' results <- ipsi_deriv(y, a, x.trt, x.out, id, d.seq)
#' @references McClean et al. Nonparametric estimation of conditional 
#' incremental effects.
#' \href{https://arxiv.org/abs/2212.03578}{arxiv:2212.03578}
#
ipsi_deriv <- function(y, a, x.trt, x.out, id, delta.seq, nsplits = 2, 
                       ci_level = 0.95, progress_bar = TRUE, 
                       return_ifvals = FALSE, fit = "rf",
                       sl.lib = c("SL.earth", "SL.gam", "SL.glm", 
                                  "SL.glm.interaction", "SL.mean", "SL.ranger",
                                  "SL.rpart")) {
  
  library("SuperLearner")
  library("earth")
  library("gam")
  library("ranger")
  library("rpart")
  
  # Setup storage
  n <- length(unique(id))
  dat <- data.frame(id = id, a = a, y = y,
                    ps = rep(NA, n), mu0 = rep(NA, n), mu1 = rep(NA, n))
  k <- length(delta.seq)
  ifvals <- matrix(nrow = n, ncol = k)
  est.eff <- rep(NA, k)
  x.trt <- data.frame(x.trt)
  x.out <- data.frame(x.out); x.out$a <- a
  
  if (progress_bar) {
    pb <- txtProgressBar(
      min = 0, max = nsplits + 2 * k + 1,
      style = 3
    )
  }
  
  s <- sample(rep(seq_len(nsplits), ceiling(n / nsplits))[seq_len(n)])
  
  if (progress_bar) {
    pbcount <- 0
  }
  
  # Counterfactual case for treatment for A = 1 and A = 0 
  newx1 <- x.out; newx1$a <- 1; newx0 <- x.out; newx0$a <- 0
  
  # Estimate nuisance functions using sample splitting
  for (split in seq_len(nsplits)) {
    
    if (progress_bar) {
      Sys.sleep(0.1)
      setTxtProgressBar(pb, pbcount)
      pbcount <- pbcount + 1
    }
    
    # Fit treatment model
    if (fit == "rf"){
      if (nsplits == 1) {
        trtmod <- ranger::ranger(stats::as.formula("a ~ ."),
                                 dat = cbind(x.trt, a = dat$a))
        dat$ps <- predict(trtmod, data = x.trt)$predictions        
      } else {
        trtmod <- ranger::ranger(stats::as.formula("a ~ ."),
                                 dat = cbind(x.trt, a = dat$a)[s != split, ])
        dat$ps[s == split] <- predict(trtmod, data = x.trt[s == split, ])$predictions
      }
    }
    
    if (fit == "sl"){
      if (nsplits == 1) {
        trtmod <- SuperLearner(dat$a, 
                               x.trt,
                               SL.library = sl.lib, 
                               family = binomial)
        
        dat$ps <- trtmod$SL.predict
        
      } else {
        trtmod <- SuperLearner(dat$a[s != split], x.trt[s != split,],
                               newX = x.trt[s == split, ], 
                               SL.library = sl.lib, 
                               family = binomial)
        
        dat$ps[s == split] <- trtmod$SL.predict
      }
    }
    
    # Fit outcome model
    split_length <- nrow(dat[s == split,])
    if (fit == "rf"){
      if (nsplits == 1) {
        outmod <- ranger::ranger(stats::as.formula("y ~ ."), dat = cbind(x.out, y))
        
        dat$mu1 <- predict(outmod, data = newx1)$predictions
        dat$mu0 <- predict(outmod, data = newx0)$predictions
        
      } else {
        outmod <- ranger::ranger(stats::as.formula("y ~ ."),
                                 dat = cbind(x.out, y)[s != split, ])
        
        dat$mu1[s == split] <- predict(outmod, data = newx1[s == split, ])$predictions
        dat$mu0[s == split] <- predict(outmod, data = newx0[s == split, ])$predictions
      }
    }
      
    if (fit == "sl"){
      if (nsplits == 1) {
        outmod <- SuperLearner(y, x.out,
                               SL.library = sl.lib,
                               newX = rbind(newx1, newx0))
        
        dat$mu1 <- outmod$SL.predict[1:n]
        dat$mu0 <- outmod$SL.predict[(n + 1):(2 * n)]
        
      } else {
        outmod <- SuperLearner(y[s != split],
                               x.out[s != split,],
                               SL.library = sl.lib,
                               newX = rbind(newx1[s == split, ], 
                                            newx0[s == split, ]))
        
        dat$mu1[s == split] <- outmod$SL.predict[1:split_length]
        dat$mu0[s == split] <- outmod$SL.predict[(split_length + 1):(2 * split_length)]
      }
    }
        
  }
  
  # Compute influence function values
  for (j in seq_len(k)) {
    
    if (progress_bar) {
      Sys.sleep(0.1)
      setTxtProgressBar(pb, pbcount)
      pbcount <- pbcount + 1
    }
    
    delta <- delta.seq[j]

    ifvals[, j] <- dat$ps * (1 - dat$ps) / (delta * dat$ps + 1 - dat$ps)^2 *  
      (dat$mu1 - dat$mu0 + (dat$a / dat$ps) * (dat$y - dat$mu1) -
         ((1 - dat$a) / (1 - dat$ps)) * (dat$y - dat$mu0)) +
      (dat$a - dat$ps) * (dat$mu1 - dat$mu0) * 
      (1 / (delta * dat$ps + 1 - dat$ps)^2 - 
         2 * delta * dat$ps / (delta * dat$ps + 1 - dat$ps)^3)
    
  }
  
  # Compute estimates
  for (j in seq_len(k)) {
    
    if (progress_bar) {
      Sys.sleep(0.1)
      setTxtProgressBar(pb, pbcount)
      pbcount <- pbcount + 1
    }
    
    est.eff[j] <- mean(ifvals[, j])
  }
  
  # Compute asymptotic variance
  sigma <- sqrt(apply(ifvals, 2, var))
  ci_norm_bounds <- abs(stats::qnorm(p = (1 - ci_level) / 2))
  eff.ll <- est.eff - ci_norm_bounds * sigma / sqrt(n)
  eff.ul <- est.eff + ci_norm_bounds * sigma / sqrt(n)
  
  # Multiplier bootstrap
  if (progress_bar) {
    Sys.sleep(0.1)
    setTxtProgressBar(pb, pbcount)
    pbcount <- pbcount + 1
  }
  
  eff.mat <- matrix(rep(est.eff, n), nrow = n, byrow = TRUE)
  sig.mat <- matrix(rep(sigma, n), nrow = n, byrow = TRUE)
  ifvals2 <- (ifvals - eff.mat) / sig.mat
  nbs <- 10000
  mult <- matrix(2 * rbinom(n * nbs, 1, 0.5) - 1, nrow = n, ncol = nbs)
  maxvals <- sapply(seq_len(nbs), function(col) {
    max(abs(apply(mult[, col] * ifvals2, 2, sum) / sqrt(n)))
  })
  calpha <- quantile(maxvals, ci_level)
  eff.ll2 <- est.eff - calpha * sigma / sqrt(n)
  eff.ul2 <- est.eff + calpha * sigma / sqrt(n)
  
  if (progress_bar) {
    Sys.sleep(0.1)
    setTxtProgressBar(pb, pbcount)
    close(pb)
  }
  
  res <- data.frame(
    increment = delta.seq, est = est.eff, se = sigma,
    ci.ll = eff.ll2, ci.ul = eff.ul2
  )
  res2 <- data.frame(
    increment = delta.seq, est = est.eff, se = sigma,
    ci.ll = eff.ll, ci.ul = eff.ul
  )
  
  # Output
  if (return_ifvals) {
    return(invisible(list(
      res = res, res.ptwise = res2, calpha = calpha,
      ifvals = (ifvals - est.eff)
    )))
  } else {
    return(invisible(list(res = res, res.ptwise = res2, calpha = calpha)))
  }
}

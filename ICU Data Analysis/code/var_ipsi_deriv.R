#' @title Estimating the variance of derivatives of effects of incremental 
#' propensity score interventions
#'
#' @description \code{var_ipsi_deriv} is used to estimate the variance of the 
#' derivative of effects of incremental propensity score interventions, i.e., 
#' estimates of the variance of the change in outcomes if the odds of receiving 
#' treatment were increased infitesimally from those corresponding to 
#' incremental parameter delta (which are themselves the odds of treatment if 
#' the natural odds of receiving treatment were multiplied by delta). This 
#' estimand is a summarize of treatment effect heterogeneity, like the variance
#' of the conditional average treatment effect.
#'
#' @usage var_ipsi_deriv(y, a, x.trt, x.out, cond.vars id, delta.seq, nsplits, 
#'  ci_level = 0.95, progress_bar = TRUE, return_ifvals = FALSE, fit,
#'  sl.lib=c("SL.earth","SL.gam","SL.glm","SL.glmnet","SL.glm.interaction", 
#'  "SL.mean","SL.ranger","rpart"))
#'
#' @param y Outcome of interest measured at end of study.
#' @param a Binary treatment.
#' @param x.trt Covariate matrix for treatment regression.
#' @param x.out Covariate matrix for outcome regression.
#' @param cond.vars Integer list of covariates over which to calculate variance.
#' Note --- these must be indexes in \code{x.trt}, not the full dataset.
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
#' @param sl.lib sl.lib algorithm library for SuperLearner. Default library 
#'  includes "earth", "gam", "glm", "glmnet", "glm.interaction", "mean", 
#'  "ranger", and "rpart.
#'
#' @section Details:
#' Treatment and covariates are expected to be one timepoint. Therefore if
#' \code{n}, then \code{a}, \code{y}, and \code{id} should all be vectors of 
#' length \code{n}, and \code{x.trt} and \code{x.out} should be matrices with 
#' \code{n} rows. The subject ordering should be consistent  across function 
#' inputs, based on the ordering specified by \code{id}. The covariate ordering
#' should be consistent across \code{x.trt} and \code{x.out}. \code{cond.vars} 
#' is assumed to contain a list of integers, corresponding to columns in
#' \code{x.trt} 
#'
#' @return A list containing the following components:
#' \item{res}{ estimates/SEs and pointwise CIs for population means.}
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
#' d.seq <- seq(0.5, 2, length.out = 10)
#'
#' results <- var_ipsi_deriv(y, a, x.trt, x.out, cond.vars = c(1,2), id, d.seq)
#' @references McClean et al. Nonparametric estimation of conditional 
#' incremental effects.
#' \href{https://arxiv.org/abs/2212.03578}{arxiv:2212.03578}
#
var_ipsi_deriv <- function(y, a, x.trt, x.out, cond.vars, id, delta.seq, 
                           nsplits = 2, ci_level = 0.95, progress_bar = TRUE, 
                           return_ifvals = FALSE, fit = "rf",
                           sl.lib = c("SL.earth", "SL.gam", "SL.glm", 
                                      "SL.glm.interaction", "SL.mean", 
                                      "SL.ranger", "SL.rpart")) {
  
  library("SuperLearner")
  library("earth")
  library("gam")
  library("ranger")
  library("rpart")
  
  # Setup
  n <- length(unique(id))
  dat <- data.frame(id = id, a = a, y = y,
                    ps = rep(NA, n), mu0 = rep(NA, n), mu1 = rep(NA, n),
                    cide = rep(NA, n))
  k <- length(delta.seq)
  ifvals <- matrix(nrow = n, ncol = k)
  ifvals_t1 <- ifvals_t2 <- ifvals 
  est.eff <- rep(NA, k)
  x.trt <- data.frame(x.trt)
  x.out <- data.frame(x.out); x.out$a <- a
  
  # Print names of conditioning variables; stop if index out of bounds
  if (max(cond.vars) > ncol(x.trt) || min(cond.vars) < 1) {
    stop("Conditioning variable subscript out of bounds. Please check and re-run.")
  } else {
    cat("Conditioning variables:", colnames(x.trt)[cond.vars], "\n")
  }
  
  if (progress_bar) {
    pb <- txtProgressBar(
      min = 0, max = nsplits + k + 1,
      style = 3
    )
  }
  
  s <- sample(rep(seq_len(nsplits), ceiling(n / nsplits))[seq_len(n)])
  
  if (progress_bar) {
    pbcount <- 0
  }

  # Counterfactual case for treatment for A = 1 and A = 0 
  newx1 <- x.out; newx1$a <- 1; newx0 <- x.out; newx0$a <- 0
  
  # Estimate propensity score and outcome regression using sample splitting
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
  
  # Iterate through delta values to calculate conditional incremental
  # derivative effects and then influence function values
  for (j in seq_len(k)) {
      
    if (progress_bar) {
      setTxtProgressBar(pb, pbcount)
      pbcount <- pbcount + 1
    }
    
    delta <- delta.seq[j]

    # Calculate pseudo-outcomes for estimating conditional incremental 
    # derivative effect
    omega <- dat$ps * (1 - dat$ps) / ((delta * dat$ps + 1 - dat$ps)^2)
    tau <- dat$mu1 - dat$mu0
    eif_omega <- (dat$a - dat$ps) * 
      (1 / (delta * dat$ps + 1 - dat$ps)^3 -
         (2 * delta * dat$ps) / (delta * dat$ps + 1 - dat$ps)^2)
    eif_tau <- (dat$a / dat$ps) * (dat$y - dat$mu1) + 
      ((1 - dat$a) / (1 - dat$ps)) * (dat$y - dat$mu0)
    
    pseudo <- omega * tau + eif_omega * tau + omega * eif_tau
    
    # If conditionining on all covariates or if not sample splitting, use
    # the pseudo-value
    if (length(cond.vars) == ncol(x.trt) || nsplits == 1) {
      
      dat$cide <- pseudo 
      
    }
    
    # Otherwise, regress the pseudo-value against conditioning covariates
    if (length(cond.vars) < ncol(x.trt) & nsplits > 1) {
      
      for (split in seq_len(nsplits)) {
        
        # Run model on training data --- estimate on held out split.
        if (fit == "rf") {
          cidemod <- 
            ranger::ranger(
              stats::as.formula("pseudo ~ ."),
              dat = cbind(x.trt[, cond.vars, drop = F], pseudo)[s != split, ]
            )
          dat$cide[s == split] <- 
            predict(cidemod, data = x.trt[s == split, cond.vars, drop = F])$predictions
        }
        
        if (fit == "sl") {
          cidemod <- SuperLearner(pseudo[s != split],
                                  x.trt[s != split, cond.vars, drop = F],
                                  newX = x.trt[s == split, cond.vars, drop = F],
                                  SL.library = sl.lib)
          dat$cide[s == split] <- cidemod$SL.predict
        }
      
      }
    }
      
    # Compute influence function values
    ifvals_t1[, j] <- (dat$cide)^2 + 2 * dat$cide *
      (omega * tau + omega * eif_tau + eif_omega * tau - (dat$cide)^2)
    
    ifvals_t2[, j] <- mean(omega * tau + omega * eif_tau + eif_omega * tau) *
      (omega * tau + omega * eif_tau + eif_omega * tau)
    
    ifvals[, j] <- ifvals_t1[, j] - ifvals_t2[, j]
    
  }
  
  # Compute estimates
  for (j in seq_len(k)) {
    est.eff[j] <- mean(ifvals[, j])
  }

  
  if (progress_bar) {
    setTxtProgressBar(pb, pbcount)
    pbcount <- pbcount + 1
  }
  
  # Compute asymptotic variance 
  sigma <- sqrt(apply(ifvals, 2, var))
  ci_norm_bounds <- abs(stats::qnorm(p = (1 - ci_level) / 2))
  eff.ll <- est.eff - ci_norm_bounds * sigma / sqrt(n)
  eff.ul <- est.eff + ci_norm_bounds * sigma / sqrt(n)
  
  # Compute asymptotic variance for when the estimand is on the boundary of
  # the parameter space (VCIDE = 0)
  sigma2 <- sqrt(apply(ifvals_t1, 2, var) + apply(ifvals_t2, 2, var))
  ci_norm_bounds <- abs(stats::qnorm(p = (1 - ci_level) / 2))
  eff.ll2 <- est.eff - ci_norm_bounds * sigma2 / sqrt(n)
  eff.ul2 <- est.eff + ci_norm_bounds * sigma2 / sqrt(n)
  
  if (progress_bar) {
    setTxtProgressBar(pb, pbcount)
    close(pb)
  }
  
  res <- data.frame(
    increment = delta.seq, est = est.eff, se = pmax(sigma, sigma2),
    ci.ll = pmin(eff.ll, eff.ll2), ci.ul = pmax(eff.ul, eff.ul2)
  )

  # Output
  if (return_ifvals) {
    return(invisible(list(
      res = res, ifvals = (ifvals - est.eff)
    )))
  } else {
    return(invisible(list(res = res)))
  }
}

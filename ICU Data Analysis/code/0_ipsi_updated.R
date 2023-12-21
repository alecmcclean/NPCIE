### This is the IPSI function from npcausal, but with some added commands to
### return the propensity scores and the split information
ipsi_update <- function (y, a, x.trt, x.out, time, id, delta.seq, nsplits = 2, 
                         ci_level = 0.95, progress_bar = TRUE, return_ifvals = FALSE, 
                         fit = "rf", sl.lib = c("SL.earth", "SL.gam", 
                                                "SL.glm", "SL.glm.interaction", "SL.mean", 
                                                "SL.ranger", "SL.rpart")) 
{
  require("SuperLearner")
  require("earth")
  require("gam")
  require("ranger")
  require("rpart")
  ntimes <- length(table(time))
  end <- max(time)
  n <- length(unique(id))
  ynew <- rep(NA, n * ntimes)
  ynew[time == end] <- y
  dat <- data.frame(time = time, id = id, y = ynew, a = a)
  k <- length(delta.seq)
  ifvals <- matrix(nrow = n, ncol = k)
  q_scores <- matrix(nrow = n, ncol = k)
  est.eff <- rep(NA, k)
  wt <- matrix(nrow = n * ntimes, ncol = k)
  cumwt <- matrix(nrow = n * ntimes, ncol = k)
  rt <- matrix(nrow = n * ntimes, ncol = k)
  vt <- matrix(nrow = n * ntimes, ncol = k)
  x.trt <- data.frame(x.trt)
  x.out <- data.frame(x.out)
  x.out$a <- a
  if (progress_bar) {
    pb <- txtProgressBar(min = 0, max = 2 * nsplits * length(delta.seq) + 
                           3, style = 3)
  }
  s <- sample(rep(seq_len(nsplits), ceiling(n/nsplits))[seq_len(n)])
  slong <- rep(s, rep(ntimes, n))
  if (progress_bar) {
    pbcount <- 0
  }
  for (split in seq_len(nsplits)) {
    if (progress_bar) {
      Sys.sleep(0.1)
      setTxtProgressBar(pb, pbcount)
      pbcount <- pbcount + 1
    }
    if (fit == "rf") {
      trtmod <- ranger::ranger(stats::as.formula("a ~ ."), 
                               dat = cbind(x.trt, a = dat$a)[slong != split, 
                               ])
      dat$ps <- predict(trtmod, data = x.trt)$predictions
    }
    if (fit == "sl") {
      trtmod <- SuperLearner(Y = dat$a[slong != split], 
                             X = x.trt[slong != split, ], 
                             newX = x.trt, 
                             SL.library = sl.lib, 
                             family = binomial)
      dat$ps <- trtmod$SL.predict
    }
    for (j in seq_len(k)) {
      if (progress_bar) {
        Sys.sleep(0.1)
        setTxtProgressBar(pb, pbcount)
        pbcount <- pbcount + 1
      }

      delta <- delta.seq[j]
      
      # Save shifted propensity scores
      q_scores[s == split, j] <- delta * dat$ps[s == split] / 
        (delta * dat$ps[s == split] + 1 - dat$ps[s == split])
      
      wt[, j] <- (delta * dat$a + 1 - dat$a)/(delta * dat$ps + 1 - dat$ps)
      cumwt[, j] <- as.numeric(t(aggregate(wt[, j], by = list(dat$id), 
                                           cumprod)[, -1]))
      vt[, j] <- (1 - delta) * (dat$a * (1 - dat$ps) - 
                                  (1 - dat$a) * delta * dat$ps)/delta
      outmod <- vector("list", ntimes)
      rtp1 <- dat$y[dat$time == end]
      if (progress_bar) {
        Sys.sleep(0.1)
        setTxtProgressBar(pb, pbcount)
        pbcount <- pbcount + 1
      }
      for (i in seq_len(ntimes)) {
        t <- rev(unique(dat$time))[i]
        newx1 <- x.out[dat$time == t, ]
        newx1$a <- 1
        newx0 <- x.out[dat$time == t, ]
        newx0$a <- 0
        if (fit == "rf") {
          outmod[[i]] <- ranger::ranger(stats::as.formula("rtp1 ~ ."), 
                                        dat = cbind(x.out, rtp1)[dat$time == t & 
                                                                   slong != split, ])
          m1 <- predict(outmod[[i]], data = newx1)$predictions
          m0 <- predict(outmod[[i]], data = newx0)$predictions
        }
        if (fit == "sl") {
          print(c(i, j))
          flush.console()
          outmod[[i]] <- SuperLearner(rtp1[s != split], 
                                      x.out[dat$time == t & slong != split, ], 
                                      SL.library = sl.lib, newX = rbind(newx1, newx0))
          m1 <- outmod[[i]]$SL.predict[1:dim(newx1)[1]]
          m0 <- outmod[[i]]$SL.predict[(dim(newx1)[1] + 
                                          1):(dim(newx1)[1] + dim(newx0)[1])]
        }
        pi.t <- dat$ps[dat$time == t]
        rtp1 <- (delta * pi.t * m1 + (1 - pi.t) * m0)/(delta * 
                                                         pi.t + 1 - pi.t)
        rt[dat$time == t, j] <- rtp1
      }
      ifvals[s == split, j] <- ((cumwt[, j] * dat$y)[dat$time == 
                                                       end] + aggregate(cumwt[, j] * vt[, j] * rt[, 
                                                                                                  j], by = list(dat$id), sum)[, -1])[s == split]
    }
  }
  for (j in seq_len(k)) {
    est.eff[j] <- mean(ifvals[, j])
  }
  sigma <- sqrt(apply(ifvals, 2, var))
  ci_norm_bounds <- abs(stats::qnorm(p = (1 - ci_level)/2))
  eff.ll <- est.eff - ci_norm_bounds * sigma/sqrt(n)
  eff.ul <- est.eff + ci_norm_bounds * sigma/sqrt(n)

  if (progress_bar) {
    Sys.sleep(0.1)
    setTxtProgressBar(pb, pbcount)
    pbcount <- pbcount + 1
  }
  eff.mat <- matrix(rep(est.eff, n), nrow = n, byrow = TRUE)
  sig.mat <- matrix(rep(sigma, n), nrow = n, byrow = TRUE)
  ifvals2 <- (ifvals - eff.mat)/sig.mat
  nbs <- 10000
  mult <- matrix(2 * rbinom(n * nbs, 1, 0.5) - 1, nrow = n, 
                 ncol = nbs)
  maxvals <- sapply(seq_len(nbs), function(col) {
    max(abs(apply(mult[, col] * ifvals2, 2, sum)/sqrt(n)))
  })
  calpha <- quantile(maxvals, ci_level)
  eff.ll2 <- est.eff - calpha * sigma/sqrt(n)
  eff.ul2 <- est.eff + calpha * sigma/sqrt(n)
  if (progress_bar) {
    Sys.sleep(0.1)
    setTxtProgressBar(pb, pbcount)
    close(pb)
  }
  res <- data.frame(increment = delta.seq, est = est.eff, se = sigma, 
                    ci.ll = eff.ll2, ci.ul = eff.ul2)
  res2 <- data.frame(increment = delta.seq, est = est.eff, 
                     se = sigma, ci.ll = eff.ll, ci.ul = eff.ul)
  if (return_ifvals) {
    return(invisible(list(res = res, 
                          res.ptwise = res2, 
                          calpha = calpha, 
                          ifvals = (ifvals - est.eff), 
                          splits = s,
                          q_scores = q_scores)))
  }
  else {
    return(invisible(list(res = res, res.ptwise = res2, calpha = calpha)))
  }
}




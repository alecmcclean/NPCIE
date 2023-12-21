###############################################################################
### Author: Alec McClean
### Purpose: Functions for simulation study on conditional incremental effects
###############################################################################


#################################
### Helper functions

# QTransform: Calculate transformed propensity score
# Input: delta - odds multiplier, pi - the propensity score
# Output: Transformed propensity score
QTransform <- function(delta, pi) {
  return(delta * pi / (delta * pi + 1 - pi))
}

# PhiTransform: Calculate efficient influence function for incremental 
# propensity score (labelled phi in McClean et al. 2023)
# Input: delta - odds multiplier  A - treatment values, pi - propensity score
# Output: EIF for incremental propensity score
PhiTransform <- function(delta, A, pi) {
  return(delta * (A - pi) / ((delta * pi + 1 - pi)^2))
}


# IFValsCIE: Calculate efficient influence function values for incremental effect
# Input: Y - outcome, A - treatment, pi - the propensity score,
#        mu1 - outcome regression when A = 1, 
#        mu0 - outcome regression when A = 0
#        q - transformed pi (see QTransform),
#        phi - influence function values for q (see PhiTransform)
# Output: Incremental effect efficient influence function values
IFValsCIE <- function(Y, A, pi, mu1, mu0, q, phi) {
  return(
    (Y - (A * mu1 + (1 - A) * mu0)) * # Residual
      (A * q + (1 - A) * (1 - q)) / (A * pi + (1 - A) * (1 - pi)) + # Inverse weight
      mu1 * q + mu0 * (1 - q) + # Plug-in 
      (mu1 - mu0) * phi # Uncertainty estimating q
  )
}

# expit: Logistic transformation
expit <- function(x){ exp(x) / (1 + exp(x)) }

# logit: Inverse logistic transformation
logit <- function(x){ log(x / (1 - x)) }


# CreatePiMu0: Create propensity score and mu0 based on covariate X
CreatePiMu0 <- function(dat) {
  # Simple propensity score model
  dat$pi <- expit(dat$X / 2)
  
  # Complicated outcome regression model (sim 1 from Polley & VdL)
  dat %<>% mutate(
    mu0 = (X < -3) * -2 +
      (X > -2) * 2.55 +
      (X > 0) * -2 +
      (X > 2) * 4 +
      (X > 3) * -1
  )
  
  return(dat)
}

# CreateAYVariables: Generate treatment (A) and outcome (Y) variables based on 
# propensity score (pi) and mu0, mu1
CreateAYVariables <- function(dat) {
  assert_that("pi" %in% colnames(dat))
  assert_that("mu0" %in% colnames(dat))
  assert_that("mu1" %in% colnames(dat))
  
  N <- nrow(dat)
  
  # Generate A
  dat$A <- rbinom(N, 1, prob = dat$pi)
  
  # Generate Y
  dat$Y <- dat$A * dat$mu1 + (1 - dat$A) * dat$mu0 + rnorm(N)
  
  return(dat)
}

# CreateCICE: Create CICE-specific variable tau_cice
CreateCICE <- function(dat) {
  dat$tau_cice <- 1 + 0.5 * dat$X - 0.2 * (dat$X)^2
  return(dat)
}


###################################
### Generate data for main figures

# generate_data: Generate synthetic data for the simulation study
# Input: N - the number of samples, 
#        seed - random seed,
#        delta_up and delta_down - \delta_u and \delta_l for CICE
#        pos_viol - flag for creating positivity violations
# Output: Synthetic dataset
generate_data <- function(N, seed, delta_up, delta_down) {
  set.seed(seed)
  
  if (is.numeric(N)) {
    # Create covariate
    dat <- data.frame(X = runif(n = N, min = -4, max = +4))
  } else {
    dat <- N
  }
  
  # Calculate propensity score (sim 2 from Polley * VdL 2007) and mu0
  dat <- CreatePiMu0(dat) 
  
  # Calculate q_up and q_down
  dat$q_up <- QTransform(delta_up, dat$pi)
  dat$q_down <- QTransform(delta_down, dat$pi)
  
  # Define tau_cate
  dat <- CreateCICE(dat)
  
  # Calculate tau_cate from tau_cice
  dat$tau_cate <- dat$tau_cice / (dat$q_up - dat$q_down)
  
  # Calculate mu1 from tau_CATE
  dat$mu1 <- dat$mu0 + dat$tau_cate
  
  # Generate A & Y
  dat <- CreateAYVariables(dat)
  
  return(dat)
}

#####################################
### Helper function: construct outcomes

# ConstructOutcomes: Construct outcomes for the simulation study
# Input: dat - the synthetic dataset, 
#        rate_pi - convergence rate propensity score estimator
#        rate_mu - convergence rate outcome regression estimator 
#        delta_up and delta_down - \delta_u and \delta_l for CICE
# Output: Dataset with constructed outcomes
ConstructOutcomes <- function(dat, rate_pi, rate_mu, delta_up, delta_down) {
  N <- nrow(dat)
  
  ##############################################
  ### Construct nuisance function "estimates"
  ### Essentially, add Gaussian noise
  
  dat$pihat <- expit(logit(dat$pi) + rnorm(N, 
                                           mean = 1 / (N^rate_pi), 
                                           sd = 1 / (N^rate_pi)))
  
  ### Outcome models
  span_mu1 <- max(dat$mu1) - min(dat$mu1)
  span_mu0 <- max(dat$mu0) - min(dat$mu0)
  dat$mu1hat <- dat$mu1 + rnorm(N, 
                                mean = span_mu1 / (N^rate_mu), 
                                sd = span_mu1 / (N^rate_mu))
  dat$mu0hat <- dat$mu0 + rnorm(N, 
                                mean = span_mu0 / (N^rate_mu), 
                                sd = span_mu0 / (N^rate_mu))
  
  ########################################
  ### Calculate pseudo outcomes
  
  ### CIE
  # Oracle pseudo outcomes
  dat %<>%
    mutate(
      phi_up = PhiTransform(delta_up, A, pi),
      phi_down = PhiTransform(delta_down, A, pi),
      
      # Pseudo values for regression with higher delta
      pseudo_up = IFValsCIE(Y, A, pi, mu1, mu0, q = q_up, phi = phi_up),
      
      # Pseudo values for regression with lower delta
      pseudo_down = IFValsCIE(Y, A, pi, mu1, mu0, q = q_down, phi = phi_down),
      
      # Pseudo values for contrast
      pseudo_cice_oracle = pseudo_up - pseudo_down
    )
  
  dat %<>% select(-phi_up, -phi_down, -pseudo_up, -pseudo_down)
  
  # Estimated pseudo outcomes
  dat %<>%
    mutate(
      
      qhat_up = QTransform(delta_up, pihat),
      qhat_down = QTransform(delta_down, pihat),
      
      phihat_up = PhiTransform(delta_up, A, pihat),
      phihat_down = PhiTransform(delta_down, A, pihat),
      
      # Pseudo values for regression with higher delta
      pseudo_up = IFValsCIE(Y, A, pihat, mu1 = mu1hat, mu0 = mu0hat, q = qhat_up, phi = phihat_up),
      
      # Pseudo values for regression with lower delta
      pseudo_down = IFValsCIE(Y, A, pihat, mu1 = mu1hat, mu0 = mu0hat, q = qhat_down, phi = phihat_down),
    )
  
  # Pseudo values for contrast
  dat$pseudo_cice = dat$pseudo_up - dat$pseudo_down
  dat %<>% select(-phihat_up, -phihat_down, -pseudo_up, -pseudo_down)
  
  return(dat)
}

#########################################
### EstimateContrast function

# EstimateContrast: Estimate contrast values and their mean squared errors (MSEs) for a range of parameters
# Input: N - the number of samples, 
#        seed - random seed,
#        rate_seq - convergence rates for propensity score outcome regression estimators,
#        delta_up and delta_down - \delta_u and \delta_l for CICE
#        pos_violations - flag for creating positivity violations
# Output: Data frame with MSE values for different contrast estimation methods
EstimateContrast <- function(N, seed, rate_seq, delta_up, delta_down) {
  # Generate dataset
  dat <- generate_data(N, seed, delta_up, delta_down)
  
  # Create output data frame
  mse_output <- expand.grid(rate_pi = rate_seq, rate_mu = rate_seq)
  mse_output$cice_baseline <- NA; mse_output$cice_dr <- NA; mse_output$cice_or <- NA
  mse_output$intercept_in <- NA; mse_output$beta1_in <- NA; mse_output$beta2_in <- NA
  
  # Estimate nuisance functions and return MSE of each estimator
  for (index in 1:nrow(mse_output)) {
    
    cat("\nRound ", index, " out of ", nrow(mse_output))
    
    rate_pi <- mse_output$rate_pi[index]
    rate_mu <- mse_output$rate_mu[index]
    
    # Construct outcomes
    dat <- ConstructOutcomes(dat, rate_pi, rate_mu, delta_up, delta_down)
    
    ###################################
    ### Construct estimators
    
    # CIE plug-in  
    dat %<>% mutate(baseline_cice = (qhat_up - qhat_down) * (mu1hat - mu0hat))
    
    # CIE DR
    dat$drl_cice <- predict(smooth.spline(x = dat$X, y = dat$pseudo_cice), x = dat$X)$y
    
    # CIE oracle
    dat$oracle_cice <- predict(smooth.spline(x = dat$X, y = dat$pseudo_cice_oracle), x = dat$X)$y
    
    ### Estimate a working model assuming tau_cice = beta_0 + beta_1 X + beta_2 X^2
    working_mod_cice <- lm(data = dat, formula = pseudo_cice ~ X + I(X^2))
    
    COEFFICIENT_SDs <- sqrt(diag(vcovHC(working_mod_cice, type = "HC")))
    
    # Apply Bonferroni correction for testing three coefficients
    scaling_factor <- qnorm(1 - 0.025)
    
    LOWER <- coef(working_mod_cice) - scaling_factor * COEFFICIENT_SDs
    UPPER <- coef(working_mod_cice) + scaling_factor * COEFFICIENT_SDs
    
    mse_output$intercept_in[index] <- LOWER[[1]] < 1 & UPPER[[1]] > 1
    mse_output$beta1_in[index] <- LOWER[[2]] < 0.5 & UPPER[[2]] > 0.5
    mse_output$beta2_in[index] <- LOWER[[3]] < -0.2 & UPPER[[3]] > -0.2
    
    
    ###############################
    ### Calculate MSEs
    
    mse_output$cice_baseline[index] <- mean((dat$baseline_cice - dat$tau_cice)^2)
    mse_output$cice_dr[index] <- mean((dat$drl_cice - dat$tau_cice)^2)
    mse_output$cice_or[index] <- mean((dat$oracle_cice - dat$tau_cice)^2)
  }
  
  return(mse_output)
}

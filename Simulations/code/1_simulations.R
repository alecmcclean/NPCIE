###############################################################################
### Author: Alec McClean
### Purpose: Estimate CICE on simulated data with I-DR-Learner and 
### Projection-Learner and nuisance function estimators converging 
### at different rates
###############################################################################

# Number of iterations for the simulation
NUM_ITERS <- 1000

# Initialize an empty data frame to store the results
bind_output <- data.frame()

# Run simulations at n = 1,000 and n = 10,000
for (sample_size in c(1000, 10000)) {
  
  # Loop through multiple iterations for the simulation
  for (i in 1:NUM_ITERS) {
    
    cat("\n------------------------------------",
        "\n--- Run: ", i, " out of ", NUM_ITERS,
        "\n------------------------------------")
    
    # Estimate contrast models and calculate Mean Squared Errors (MSE)
    output <- EstimateContrast(N = sample_size,
                               seed = 20210118 + sample_size + i,
                               rate_seq = seq(0.1, 0.5, 0.1),
                               delta_up = 5,
                               delta_down = 1/5)
    
    # Add an iteration indicator to the output
    output$iter <- i
    output$sample_size <- sample_size
    
    # Append the output to the overall result data frame
    bind_output %<>% bind_rows(output)
    
  }
}


###############################################################################
### Author: Alec McClean
### Purpose: Master script for conditional incremental effects analysis.  
### Loads packages and calls other scripts.
###
### Last run: 2022-04-18
###############################################################################

if (!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}

p_load(assertthat,
       GGally,
       ggplot2,
       ggthemes,
       magrittr,
       mgcv,
       sandwich,
       tidyverse)

options(stringsAsFactors = F)

source("0_simulation_functions.R")
source("1_simulations.R")
source("2_simulations_figures.R")

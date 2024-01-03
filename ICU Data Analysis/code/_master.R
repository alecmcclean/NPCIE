###############################################################################
### Author: Alec McClean
### Purpose: Master script for conditional incremental effects analysis.  
### Loads packages and calls other scripts.
###
### Last run: 2024-01-03
###############################################################################

if (!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}

p_load(assertthat,
       e1071,
       earth,
       gbm,
       GGally,
       ggplot2,
       ggthemes,
       glmnet,
       latex2exp,
       lubridate,
       magrittr,
       mgcv,
       modelr,
       nnet,
       np,
       pander,
       plotly,
       ranger,
       rpart,
       sandwich,
       stargazer,
       SuperLearner,
       tidyverse,
       xgboost)

options(stringsAsFactors = F)

source("1_icu_data_analysis.R")
source("2_cide_data_analysis.R")

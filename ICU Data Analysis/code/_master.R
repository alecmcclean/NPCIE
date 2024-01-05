###############################################################################
### Author: Alec McClean
### Purpose: Master script for conditional incremental effects analysis.  
### Loads packages and calls other scripts.
###
### Last run: 2024-01-04
###############################################################################

if (!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}

p_load(e1071,
       earth,
       gbm,
       GGally,
       ggplot2,
       ggthemes,
       glmnet,
       magrittr,
       mgcv,
       nnet,
       np,
       pander,
       plotly,
       ranger,
       rpart,
       sandwich,
       SuperLearner,
       tidyverse,
       xgboost)

options(stringsAsFactors = F)
set.seed(20240104)

source("1_icu_data_analysis.R")
source("2_cide_data_analysis.R")

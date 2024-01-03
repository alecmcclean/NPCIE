###############################################################################
### Author: Alec McClean
### Purpose: Data analysis with ICU data
###############################################################################

source("ipsi_deriv.R")
source("var_ipsi_deriv.R")

### Column names from original analysis
COLS <- c("site", "dead28", "age", "male", "sepsis_dx",
          "periarrest", "weekend", "winter", "icnarc_score", "news_score",
          "sofa_score", "icu_bed", "open_bin", "v_cc1", "v_cc2", "v_cc3", 
          "v_cc4", "v_cc5", "v_cc_r1", "v_cc_r2", "v_cc_r3", "v_cc_r4", 
          "v_cc_r5")  

#####################################
### Load and clean data

# icu <- read.csv("../data/icuData.csv") 
### ^^^ The data used in the original analysis is not publicly available.
### Below, we use a similar, publicly available dataset
icu <- read.csv("../data/icu_pseudo_data.csv") 

# Select columns from original analysis
icu %<>% select(all_of(COLS))

# Create id column
icu$id <- seq_len(nrow(icu))

# Change variables to factors
icu %<>% mutate_at(vars(site, male, sepsis_dx:winter, v_cc1:v_cc_r5), as.factor)

# Deltas
DELTAS <- seq(0.2, 0.9, 0.1)
DELTAS <- c(DELTAS, 1, rev(1 / DELTAS))


##############################################
### Calculate Average incremental 
### derivative effect IF values

results <- ipsi_deriv(y = icu$dead28,
                      a = icu$icu_bed,
                      id = icu$id,
                      x.trt = icu %>% select(-id, -dead28, -icu_bed),
                      x.out = icu %>% select(-id, -dead28, -icu_bed),
                      fit = "rf",
                      delta.seq = DELTAS,
                      nsplits = 2,
                      return_ifvals = TRUE)

ifvals <- data.frame(ifval = as.vector(results$ifvals))
ifvals$delta <- rep(DELTAS, each = nrow(icu))
ifvals$icnarc_score <- rep(icu$icnarc_score, times = length(DELTAS))

##########################################
### Second stage regressions for CIDE

for (DELTA in unique(ifvals$delta)) {
  
  cat("\nDelta: ", DELTA)
    
  # Second stage model with smoothing spline using mgcv package
  mod <- gam(ifval ~ s(icnarc_score), data = ifvals %>% filter(delta == DELTA))
  
  # Predicted value
  ifvals$pred[ifvals$delta == DELTA] <- predict(mod)
  
  # Upper bound of pointwise 95% CI
  ifvals$upr[ifvals$delta == DELTA] <- 
    ifvals$pred[ifvals$delta == DELTA] + 1.96 * predict(mod, se.fit = TRUE)$se.fit
  
  # Lower bound of pointwise 95% CI
  ifvals$lwr[ifvals$delta == DELTA] <- 
    ifvals$pred[ifvals$delta == DELTA] - 1.96 * predict(mod, se.fit = TRUE)$se.fit
    
}


####################################
### Calculate V-CIDE across delta

x.trt.colnames <- colnames(icu %>% select(-id, -dead28, -icu_bed))
results <- var_ipsi_deriv(y = icu$dead28,
                          a = icu$icu_bed,
                          id = icu$id,
                          x.trt = icu %>% select(-id, -dead28, -icu_bed),
                          x.out = icu %>% select(-id, -dead28, -icu_bed),
                          cond.vars = which(x.trt.colnames %in% "icnarc_score"),
                          fit = "rf",
                          delta.seq = DELTAS,
                          nsplits = 2,
                          return_ifvals = TRUE)


##############################
### Create plots

### CIDE
p1 <- ifvals %>% filter(delta %in% c(0.2, 0.5, 1, 2, 5)) %>%
  ggplot(aes(x = icnarc_score, y = pred)) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  facet_wrap(~ paste0("Delta: ", substr(delta, 1, 4)),
             scales = "free_y") +
  labs(x = "ICNARC", y = "Conditional Incremental Derivative Effect Estimate") +
  theme_clean()

ggsave(plot = p1, filename = "../figures/cide.png",
       width = 8, height = 6)

### VCIDE
p2 <- ggplot(data = results$res %>% filter(increment > 0.2), 
             aes(x = as.factor(round(increment, 2)), group = 1, y = est)) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = ci.ll, ymax = ci.ul), alpha = 0.2) +
  labs(x = "Delta", y = "Variance of the Conditional Incremental Derivative Effect Estimate") +
  theme_clean()

ggsave(plot = p2, filename = "../figures/vcide.png")

rm(list = ls(all = T))
gc()


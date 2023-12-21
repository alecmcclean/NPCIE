###############################################################################
### Author: Alec McClean
### Purpose: Data analysis with ICU data
###############################################################################

#####################################
### Load and clean data

# icu <- read.csv("../data/icuData.csv") 
### ^^^ This data is not publicly available.
### There is a similar dataset available (same size, sampled with replacement)
### in ivmodel::icu.data.

# Use row number as ID
icu %<>% rename(id = X)

# Outcome variable: dead28
icu %<>% select(-dead7, -dead90) 

# Change variables to factors
icu %<>% mutate_at(vars(site, male, sepsis_dx:winter, v_cc1:v_cc_r5), as.factor)

# Deltas
DELTAS <- seq(0.2, 0.9, 0.1)
DELTAS <- c(DELTAS, 1, rev(1 / DELTAS))


##############################################
### Calculate Average incremental 
### derivative effect IF values

FOLDS <- 2
icu$fold <- sample(1:FOLDS, size = nrow(icu), replace = T)
output <- data.frame()

for (FOLD in 1:FOLDS) {
  test <- icu %>% filter(fold == FOLD)
  train <- icu %>% filter(fold != FOLD)
  
  pimod <- ranger(icu_bed ~ ., dat = train %>% select(-id, -dead28, -fold))
  mumod <- ranger(dead28 ~ ., dat = train %>% select(-id, -fold))
  
  test$pihat <- predict(pimod, data = test %>% select(-id, -dead28, -fold))$predictions
  test$mu1hat <- predict(mumod, 
                         data = test %>% select(-id, -fold) %>% mutate(icu_bed = 1))$predictions
  test$mu0hat <- predict(mumod, 
                         data = test %>% select(-id, -fold) %>% mutate(icu_bed = 0))$predictions

  output %<>% bind_rows(test)
}

### Calculate IF values at each delta
ifvals <- data.frame()
for (DELTA in DELTAS) {
  temp <- output %>% mutate(delta = DELTA)

  temp %<>% mutate(
    omega = pihat * (1 - pihat) / ((DELTA * pihat + 1 - pihat)^2),
    tau = mu1hat - mu0hat,
    eif_omega = (icu_bed - pihat) * (1 / (DELTA * pihat + 1 - pihat)^3 - (2 * DELTA * pihat) / (DELTA * pihat + 1 - pihat)^2),
    eif_tau = (icu_bed / pihat) * (dead28 - mu1hat) + ((1 - icu_bed) / (1 - pihat)) * (dead28 - mu0hat),
    plugin = omega * tau,
    eif_terms = omega * eif_tau + eif_omega * tau,
    ifval_cide = eif_terms + plugin,
  )

  temp$ifval_vcide_t2 <-  mean(temp$eif_terms + temp$plugin) * (temp$eif_terms + temp$plugin)

  ifvals %<>% bind_rows(temp)
}


##########################################
### Second stage regressions for CIDE

temp <- data.frame()
for (DELTA in unique(ifvals$delta)) {
  
  for (FOLD in 1:max(ifvals$fold)) {
    cat("\nDelta: ", DELTA)
    dat <- ifvals %>% filter(delta == DELTA) %>% filter(fold == FOLD)
    
    # Second stage model with smoothing spline using mgcv package
    mod <- gam(ifval_cide ~ s(icnarc_score), data = dat)
    dat$pred <- predict(mod)
    dat$upr <- dat$pred + 1.96 * predict(mod, se.fit = TRUE)$se.fit
    dat$lwr <- dat$pred - 1.96 * predict(mod, se.fit = TRUE)$se.fit
    
    temp %<>% bind_rows(dat)
  }
    
}

ifvals <- temp
rm(temp)
gc()

####################################
### Calculate V-CIDE across delta

### Calculate IF values 
ifvals %<>% mutate(
  ifval_vcide_t1 = (pred * pred) + 2 * pred * (eif_terms + plugin - (pred * pred)),
  ifval_vcide = ifval_vcide_t1 - ifval_vcide_t2
)

vcide <- ifvals %>% group_by(delta) %>%
  summarize(pt_est = mean(ifval_vcide),
            sd_est = sd(ifval_vcide) / sqrt(n()),
            ad_hoc = sqrt((var(ifval_vcide_t1) + var(ifval_vcide_t2)) / n())) %>%
  ungroup() %>%
  mutate(lower = pt_est - 1.96 * pmax(sd_est, ad_hoc),
         upper = pt_est + 1.96 * pmax(sd_est, ad_hoc))

##############################
### Create plots

### CIDE
p1 <- ifvals %>% filter(fold == 1, delta %in% c(0.2, 0.5, 1, 2, 5)) %>%
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
p2 <- ggplot(data = vcide, 
             aes(x = as.factor(round(delta, 2)), group = 1, y = pt_est)) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  labs(x = "Delta", y = "Variance of the Conditional Incremental Derivative Effect Estimate") +
  theme_clean()

ggsave(plot = p2, filename = "../figures/vcide.png")

rm(list = ls(all = T))
gc()


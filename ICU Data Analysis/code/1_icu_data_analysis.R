###############################################################################
### Author: Alec McClean
### Purpose: Data analysis with ICU data
###############################################################################

source("0_ipsi_updated.R") # An amended version of npcausal::ipsi() that outputs
# estimated propensity scores and influence function values.

#####################################
### Load and clean data

### The data used in the original analysis is not publicly available.
### We use a similar, publicly available dataset
icu <- read.csv("../data/icu_pseudo_data.csv")

# Create id variable, remove unused outcome variables (7- and 90-day mortality)
icu %<>% mutate(id = seq_len(nrow(icu))) %>% select(-X, -dead7, -dead90)


# Change variables to factors
icu %<>% mutate_at(vars(site, male, sepsis_dx:winter, v_cc1:v_cc_r5), as.factor)

# Deltas
DELTAS <- seq(0.2, 0.9, 0.1)
DELTAS <- c(DELTAS, 1, rev(1 / DELTAS))

##########################
### Calculate CIE

### First stage regressions
results <- ipsi_update(y = icu$dead28,
                       a = icu$icu_bed,
                       id = icu$id,
                       x.trt = icu %>% select(-id, -dead28, -icu_bed),
                       x.out = icu %>% select(-id, -dead28, -icu_bed),
                       time = rep(1, nrow(icu)),
                       fit = "rf",
                       delta.seq = DELTAS,
                       nsplits = 2,
                       return_ifvals = TRUE)

# Second stage regressions
ifvals <- as.data.frame(results$ifvals)
colnames(ifvals) <- DELTAS
ifvals$id <- icu$id 
ifvals %<>% left_join(icu %>% select(id, icnarc_score))
ifvals$split <- results$splits

ifvals %<>% 
  gather(delta, ifval, `0.2`:`5`) %>% 
  mutate(delta = round(as.numeric(delta), 2))

# Center ifvals at mean
point_estimates <- results$res %>% 
  select(delta = increment, pt = est) %>%
  mutate(delta = round(delta, 2))

ifvals %<>% left_join(point_estimates) %>%
  mutate(pseudo = pt + ifval)

for (DELTA in unique(ifvals$delta)) {
  
  cat("\nDelta: ", DELTA)
  dat <- ifvals %>% filter(delta == DELTA)
    
  # Second stage model with smoothing spline using mgcv package
  mod <- gam(pseudo ~ s(icnarc_score), data = dat)
  dat$pred <- predict(mod)
  dat$upr <- dat$pred + 1.96 * predict(mod, se.fit = TRUE)$se.fit
  dat$lwr <- dat$pred - 1.96 * predict(mod, se.fit = TRUE)$se.fit
    
  ifvals %<>% filter(delta != DELTA) %>% bind_rows(dat)
    
}


##############################
### Create plots

ifvals %<>% mutate(delta = as.factor(delta),
                   delta = reorder(delta, as.numeric(delta)))

### Propensity scores
q_scores <- results$q_scores
colnames(q_scores) <- DELTAS
q_scores %<>% as.data.frame %>% gather(delta, q, `0.2`:`5`) 
q_scores %<>% mutate(delta = factor(delta),
                     delta = reorder(delta, as.numeric(delta)))
q_scores$id <- rep(icu$id, length(unique(q_scores$delta)))

p1 <- q_scores %>% filter(delta == 1) %>%
  left_join(select(icu, id, icnarc_score)) %>% 
  ggplot(aes(x = icnarc_score, y = q, group = icnarc_score)) +
  geom_boxplot() +
  theme_clean() +
  labs(x = "ICNARC Score", y = "Estimated Propensity Score")

ggsave(filename = "../figures/diagnostics.png", plot = p1,
       width = 6, height = 4)

### 3d plot of CIE
p2 <- ifvals %>%
  ggplot(aes(x = as.factor(delta), y = icnarc_score, fill = pred)) +
  geom_tile() +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000",
                       midpoint = 0.5) + 
  labs(y = "ICNARC", x = "Delta", fill = "CIE Estimate") +
  theme_clean() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(plot = p2, filename = "../figures/3d_results_icnarc.png",
       width = 5, height = 4)

p3 <- ifvals %>%
  filter(icnarc_score %in% c(0, 15, 30, 40)) %>%
  mutate(icnarc_score = reorder(as.factor(icnarc_score), desc(icnarc_score))) %>%
  ggplot(aes(x = as.factor(delta), y = pred, group = icnarc_score)) + 
  geom_point() + 
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = factor(icnarc_score)), alpha = 0.2) +
  labs(x = "Delta", y = "CIE Estimate", fill = "ICNARC Score") +
  theme_clean() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


ggsave(plot = p3, filename = "../figures/select_results_icnarc.png",
       width = 6, height = (10/3))


### Clean
rm(list = ls(all = T))
gc()

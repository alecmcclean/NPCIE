###############################################################################
### Author: Alec McClean
### Purpose: Create figures for simulations section
###############################################################################

##############################
### Data generation figures

dat <- generate_data(N = 1000, seed = 20220516, delta_up = 5, delta_down = 1/5)

p0 <- dat %>% 
  gather(var, value, mu0, mu1, tau_cate, tau_cice, pi) %>%
  mutate(
    panel = case_when(
      var == "pi" ~ "Propensity Scores",
      var == "tau_cice" ~ "CICE",
      T ~ "Outcome Regressions and CATE"),
    var = case_when(
      var == "mu0" ~ "E(Y | X, A = 0)",
      var == "mu1" ~ "E(Y | X, A = 1)",
      var == "tau_cate" ~ "CATE",
      var == "tau_cice" ~ "CICE",
      T ~ "expit(X / 2)"
    ),
    panel = factor(panel, levels = c("Propensity Scores", "Outcome Regressions and CATE", "CICE")),
    var = factor(var, levels = c("expit(X / 2)",
                                 "E(Y | X, A = 0)", "E(Y | X, A = 1)", "CATE",
                                 "CICE"))
  ) %>%
  ggplot(aes(x = X, y = value, color = var)) +
  geom_point() + 
  labs(x = "X", y = "", color = "") +
  scale_color_colorblind() +
  facet_wrap(~ panel, scales = "free_y", ncol = 1) +
  theme_classic() 

ggsave(p0, filename = "../figures/data.png", width = 8, height = 6)



#########################################################################
### Compare I-DR-Learner, Oracle I-DR-Learner and Baseline Learner 

NUM_ITERS <- 1000
data <- bind_output

# Calculate means and SDs
data %<>%
  gather(var, value, cice_baseline:cice_or) %>%
  group_by(rate_pi, rate_mu, var, sample_size) %>%
  summarize(avg = mean(value),
            error = sd(value) * 1.96 / sqrt(NUM_ITERS)) %>%
  ungroup()

p1 <- data %>% 
  mutate(rate_mu = paste0("Estimation rate mu: ", rate_mu),
         var = case_when(
           var == "cice_dr" ~ "I-DR-Learner",
           var == "cice_or" ~ "Oracle I-DR-Learner",
           T ~ "Baseline CICE"
         ),
         var = factor(var, levels = c("Baseline CICE", "I-DR-Learner", "Oracle I-DR-Learner")),
         sample_size = paste0("N = ", sample_size)) %>%
  ggplot(aes(x = rate_pi, y = avg, color = var)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = pmax(avg - error, 0.00001), ymax = avg + error)) +
  scale_y_log10() +
  scale_color_colorblind() +
  theme_classic() +
  facet_grid(sample_size ~ rate_mu) +
  labs(x = latex2exp::TeX("Estimation rate pi ($\\alpha$)"),
       y = "MSE for estimator",
       color = "Estimator")

ggsave(p1, filename = "../figures/mse.png")


#################################################################
### Show Projection-Learner coefficient estimates coverage

data <- bind_output

p2 <- data %>%
  group_by(rate_pi, rate_mu, sample_size) %>%
  summarize_at(vars(contains("_in")), mean) %>%
  ungroup() %>%
  gather(var, value, intercept_in, beta1_in, beta2_in) %>%
  mutate(var = case_when(
    var == "intercept_in" ~ "Intercept",
    var == "beta1_in" ~ "Beta 1",
    T ~ "Beta 2"
  ),
  rate_mu = paste0("Estimation rate mu: ", rate_mu),
  sample_size = paste0("N = ", sample_size)
  ) %>%
  ggplot(aes(x = rate_pi, y = value, color = var)) +
  geom_line() +
  geom_point() +
  scale_color_colorblind() +
  theme_classic() +
  facet_grid(sample_size ~ rate_mu) +
  labs(x = latex2exp::TeX("Estimation rate pi ($\\alpha$)"),
       y = "Coverage over 1,000 simulations",
       color = "Coefficient in model") +
  geom_hline(yintercept = 0.95, linetype = "dashed")

ggsave(p2, filename = "../figures/coverage.png")


####################
### Clean

rm(list = ls(all = T))
gc()

############################
########### PS 5 ###########
############################


#load packages#
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(xtable)
library(fixest)
library(broom)
library(sandwich)
library(purrr)
library(knitr)


## Question 2 ## 

# Generate 1000 df:s # 

set.seed(123)
n_sims <- 10000
n_obs  <- 100

for (i in 1:n_sims) {
  
  x_i = rbern(n_obs, prob = 0.1) 
  e_i = rnorm(n_obs, mean = 0, sd = sqrt(2)) 
  y_i = e_i 
  
 
  if (i == 1) {
    df<- tibble(
    x = x_i, 
    e = e_i, 
    y = y_i, 
    sim = i)
    }

  else {
  df <- rbind(df, tibble(
    x = x_i, 
    e = e_i, 
    y = y_i, 
    sim = i
  ))
  
  }
}

### a ### 


# Estimating simple regression #
 
estimates <- df %>%
  group_by(sim) %>%
  do(tidy(lm(y ~ x, data = .))) %>%
  ungroup()


# Plot beta_j values # 

df_beta_j <- estimates %>% 
  filter (term == "x")

# Plot # 

plot_1 <- ggplot(df_beta_j, aes(x = estimate)) + 
  geom_density(
    fill = "skyblue",
    color = "navy",
    linewidth = 1.2,
    alpha = 0.4
  ) +
  geom_vline(
    xintercept = mean(df_beta_j$estimate),
    linetype = "dashed",
    linewidth = 1
  ) +
  labs(
    x = expression(hat(beta)[1]),
    y = "Density"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  )

print(plot_1)

ggsave(
  filename = "beta_hist.pdf",
  plot     = plot_1,
  width    = 8,
  height   = 6
)


### b ###

# create a wide table with one beta0 and beta1 per sim

betas_wide <- estimates %>%
  select(sim, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  rename(beta0 = `(Intercept)`, beta1 = x) 

# join to original df and compute residuals
df_resid <- df %>%
  left_join(betas_wide, by = "sim") %>%      
  mutate(e_hat = y - (beta0 + beta1 * x))    

# Calculate covariance between X and e_hat # 

df_resid <- df_resid %>% 
  group_by(sim) %>% 
  mutate(
    cov_x_e = cov(x, e),
    cov_x_e_hat = cov(x, e_hat)
  )

# Save estimates in separate df # 

cov_results <- df_resid %>% 
  select(sim, beta0, beta1, cov_x_e, cov_x_e_hat) %>% 
  distinct()

# Plot beta_j and cov(x, e) # 

plot_2 <- ggplot(cov_results, aes(x = beta1, y = cov_x_e)) +
  geom_point(alpha = 0.5) +
  labs(
    x = expression(hat(beta)[1]),
    y = "Cov(x,e)"
  ) +
  theme_minimal()

print(plot_2)

ggsave(
  filename = "scatter.pdf",
  plot     = plot_2,
  width    = 8,
  height   = 6
)

#Real values of variance of beta # 

summary_beta1 <- betas_wide%>%
  summarise(
    mean_beta1 = mean(beta1), 
    var_beta1 = var(beta1), 
    sd_beta_1 = sd(beta1)
  )

kable(summary_beta1,
      format = "latex",
      booktabs = TRUE,
      digits = 3,
      caption = "Population Variance of Estimators")


## c ## 

# Solve for the variances # 

res_by_sim <- df %>%
  group_by(sim) %>%
  nest() %>%
  mutate(
    model    = map(data, ~ lm(y ~ x, data = .x)),
    vcov0    = map(model, vcov),                       # classical
    vcov_HC1 = map(model, ~ vcovHC(.x, type = "HC1")),
    vcov_HC2 = map(model, ~ vcovHC(.x, type = "HC2")),
    vcov_HC3 = map(model, ~ vcovHC(.x, type = "HC3")),
    
    # extract variance for the slope coefficient named "x"
    var_beta1_v0  = map_dbl(vcov0,    ~ .x["x","x"]),
    var_beta1_hc1 = map_dbl(vcov_HC1, ~ .x["x","x"]),
    var_beta1_hc2 = map_dbl(vcov_HC2, ~ .x["x","x"]),
    var_beta1_hc3 = map_dbl(vcov_HC3, ~ .x["x","x"])
  ) %>%
  ungroup() %>%
  mutate(
    se_beta1_v0  = sqrt(var_beta1_v0),
    se_beta1_hc1 = sqrt(var_beta1_hc1),
    se_beta1_hc2 = sqrt(var_beta1_hc2),
    se_beta1_hc3 = sqrt(var_beta1_hc3)
  ) %>%
  select(sim, model, starts_with("var_"), starts_with("se_"))

# Summarise stats# 

var_sum <- res_by_sim %>% 
  ungroup() %>%
  summarise(
      
    mean_v_0 = mean(var_beta1_v0), 
    mean_v_HC1 = mean(var_beta1_hc1),
    mean_v_HC2 = mean(var_beta1_hc2), 
    mean_v_HC3 = mean(var_beta1_hc3), 
    
    var_v_0 = var(var_beta1_v0), 
    var_v_HC1 = var(var_beta1_hc1),
    var_v_HC2 = var(var_beta1_hc2),
    var_v_HC3 = var(var_beta1_hc3),
    sd_v_0 = sd(var_beta1_v0), 
    
    sd_v_HC1 = sd(var_beta1_hc1), 
    sd_v_HC2 = sd(var_beta1_hc2),
    sd_v_HC3 = sd(var_beta1_hc3),
  )
#Make table look nmicer # 

table_clean <- var_sum %>%
  pivot_longer(
    everything(),
    names_to = c("stat", "estimator"),
    names_pattern = "(mean|var|sd)_v_(.*)"
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = value
  ) %>%
  rename(
    Mean = mean,
    Variance = var,
    SD = sd
  )

# Plot as overleaf# 

kable(table_clean,
      format = "latex",
      booktabs = TRUE,
      digits = 3,
      caption = "Comparison of Variance Estimators")



### d ### 

#Create new column for sample variance# 
 
k <- 2

df_s2 <- df_resid %>%
  group_by(sim, x) %>%
  summarise(
    n_group = n(),
    s2_e = sum(e_hat^2) / (n_group - k),
    .groups = "drop"
  )

# Combine  with estimator variances # 

df_res_var <- df_s2 %>% left_join(res_by_sim, by= c("sim"))

## Compute correlations ##

corr_results <- df_res_var %>%
  filter(is.finite(s2_e)) %>% 
  group_by(x) %>%
  summarise(
    corr_beta1_hat_0   = cor(s2_e, var_beta1_v0),
    corr_beta1_hat_HC2 = cor(s2_e, var_beta1_hc2),
    .groups = "drop"
  ) 


kable(corr_results,
      format = "latex",
      booktabs = TRUE,
      digits = 3,
      caption = "{Correlation $s^2_{\varepsilon}|X$ and Variance Estimators")

summary(df_s2$s2_e)


### e ### 

#Combine data# 
df_resid <- df_resid %>% 
  left_join(res_by_sim, by= c("sim"))

# Calculate t-stat

h_0 <- 0 

df_resid <- df_resid %>% 
  group_by(sim) %>% 
  mutate(t_stat_v0 = (beta1 - h_0)/se_beta1_v0, 
         t_stat_v_HC2 = (beta1 - h_0)/se_beta1_hc2)


# Plot the distribution over samples # 

df_t_stats <- df_resid %>% 
  select(sim, t_stat_v0, t_stat_v_HC2)%>%
  unique()

df_t_stats <- df_t_stats %>% 
  pivot_longer(
    cols = -sim,                
    names_to = "Estimator",
    values_to = "t_stat")

# Plot the stats # 

plot_3 <- ggplot(df_t_stats, aes(x = t_stat, color = Estimator)) +
  geom_density(size = 1) +
  stat_function(
    fun = dnorm,
    args = list(mean = 0, sd = 1),
    linetype = "dashed",
    color = "black"
  ) +
  theme_minimal()

print(plot_3)


ggsave(
  filename = "density_t_stat.pdf",
  plot     = plot_3,
  width    = 8,
  height   = 6
)



### e ###


# Redoing the generation with non-orthogonal errors# 

set.seed(123)
n_sims <- 10000
n_obs  <- 100

#Creating df_e1# 

for (i in 1:n_sims) {
  
  x_i <- rbinom(n_obs, 1, 0.1)
  
  e_i <- rnorm(n_obs,
               mean = 0,
               sd = ifelse(x_i == 0, sqrt(2), sqrt(1)))
  
  y_i <- e_i
  
  df_i <- tibble(
    x = x_i,
    e = e_i,
    y = y_i,
    sim = i
  )
  
  if (i == 1) {
    df_e_1 <- df_i
  } else {
    df_e_1 <- rbind(df_e_1, df_i)
  }
}


#Creating df_e2# 

for (i in 1:n_sims) {
  
  x_i <- rbinom(n_obs, 1, 0.1)
  
  e_i <- rnorm(n_obs,
               mean = 0,
               sd = ifelse(x_i == 0, sqrt(2), sqrt(4)))
  
  y_i <- e_i
  
  df_i <- tibble(
    x = x_i,
    e = e_i,
    y = y_i,
    sim = i
  )
  
  if (i == 1) {
    df_e_2 <- df_i
  } else {
    df_e_2 <- rbind(df_e_2, df_i)
  }
}


### Create function ## 

run_simulation_analysis <- function(df, label){

  estimates <- df %>%
    group_by(sim) %>%
    do(tidy(lm(y ~ x, data = .))) %>%
    ungroup()
  
  
  # Plot beta_j values # 
  
  df_beta_j <- estimates %>% 
    filter (term == "x")
  
  # Plot # 
  
  
  plot_1 <- ggplot(df_beta_j, aes(x = estimate)) + 
    geom_density(
      fill = "skyblue",
      color = "navy",
      linewidth = 1.2,
      alpha = 0.4
    ) +
    geom_vline(
      xintercept = mean(df_beta_j$estimate),
      linetype = "dashed",
      linewidth = 1
    ) +
    labs(
      x = expression(hat(beta)[1]),
      y = "Density"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.minor = element_blank()
    )
  
  ggsave(
    filename = paste0("beta_hist_", label, ".pdf"),
    plot     = plot_1,
    width    = 8,
    height   = 6
  )
  
  ## b ## 
  # create a wide table with one beta0 and beta1 per sim
  
  betas_wide <- estimates %>%
    select(sim, term, estimate) %>%
    pivot_wider(names_from = term, values_from = estimate) %>%
    rename(beta0 = `(Intercept)`, beta1 = x)   
  
  # join to original df and compute residuals
  df_resid <- df %>%
    left_join(betas_wide, by = "sim") %>%      
    mutate(e_hat = y - (beta0 + beta1 * x))   
  
  ## c ## 
  
  # Solve for the variances # 
  
  res_by_sim <- df %>%
    group_by(sim) %>%
    nest() %>%
    mutate(
      model    = map(data, ~ lm(y ~ x, data = .x)),
      vcov0    = map(model, vcov),                       # classical
      vcov_HC1 = map(model, ~ vcovHC(.x, type = "HC1")),
      vcov_HC2 = map(model, ~ vcovHC(.x, type = "HC2")),
      vcov_HC3 = map(model, ~ vcovHC(.x, type = "HC3")),
      
      # extract variance for the slope coefficient named "x"
      var_beta1_v0  = map_dbl(vcov0,    ~ .x["x","x"]),
      var_beta1_hc1 = map_dbl(vcov_HC1, ~ .x["x","x"]),
      var_beta1_hc2 = map_dbl(vcov_HC2, ~ .x["x","x"]),
      var_beta1_hc3 = map_dbl(vcov_HC3, ~ .x["x","x"])
    ) %>%
    ungroup() %>%
    mutate(
      se_beta1_v0  = sqrt(var_beta1_v0),
      se_beta1_hc1 = sqrt(var_beta1_hc1),
      se_beta1_hc2 = sqrt(var_beta1_hc2),
      se_beta1_hc3 = sqrt(var_beta1_hc3)
    ) %>%
    select(sim, model, starts_with("var_"), starts_with("se_"))
  
  # Summarise stats# 
  
  var_sum <- res_by_sim %>% 
    ungroup() %>%
    summarise(
      
      mean_v_0 = mean(var_beta1_v0), 
      mean_v_HC1 = mean(var_beta1_hc1),
      mean_v_HC2 = mean(var_beta1_hc2), 
      mean_v_HC3 = mean(var_beta1_hc3), 
      
      var_v_0 = var(var_beta1_v0), 
      var_v_HC1 = var(var_beta1_hc1),
      var_v_HC2 = var(var_beta1_hc2),
      var_v_HC3 = var(var_beta1_hc3),
      sd_v_0 = sd(var_beta1_v0), 
      
      sd_v_HC1 = sd(var_beta1_hc1), 
      sd_v_HC2 = sd(var_beta1_hc2),
      sd_v_HC3 = sd(var_beta1_hc3),
    )
  
  #Make table look nicer # 
  
  table_clean <- var_sum %>%
    pivot_longer(
      everything(),
      names_to = c("stat", "estimator"),
      names_pattern = "(mean|var|sd)_v_(.*)"
    ) %>%
    pivot_wider(
      names_from = stat,
      values_from = value
    ) %>%
    rename(
      Mean = mean,
      Variance = var,
      SD = sd
    )
  
  # Plot as overleaf# 
  
  kable(table_clean,
        format = "latex",
        booktabs = TRUE,
        digits = 3,
        caption = paste0("Comparison of Variance Estimators - Case", label))
  
  ### e ### 
  
  #Combine data# 
  df_resid <- df_resid %>% 
    left_join(res_by_sim, by= c("sim"))
  
  # Calculate t-stat
  
  h_0 <- 0 
  
  df_resid <- df_resid %>% 
    group_by(sim) %>% 
    mutate(t_stat_v0 = (beta1 - h_0)/se_beta1_v0, 
           t_stat_v_HC2 = (beta1 - h_0)/se_beta1_hc2)
  
  
  # Plot the distribution over samples # 
  
  df_t_stats <- df_resid %>% 
    select(sim, t_stat_v0, t_stat_v_HC2)%>%
    unique()
  
  df_t_stats <- df_t_stats %>% 
    pivot_longer(
      cols = -sim,                
      names_to = "Estimator",
      values_to = "t_stat")
  
  # Plot the stats # 
  
  plot_3 <- ggplot(df_t_stats, aes(x = t_stat, color = Estimator)) +
    geom_density(size = 1) +
    stat_function(
      fun = dnorm,
      args = list(mean = 0, sd = 1),
      linetype = "dashed",
      color = "black"
    ) +
    theme_minimal()
  
  print(plot_3)
  
  
  ggsave(
    filename = paste0("density_t_stat", label, ".pdf"),
    plot     = plot_3,
    width    = 8,
    height   = 6
  )
}

results_e1 <- run_simulation_analysis(df_e_1, "e1")
results_e2 <- run_simulation_analysis(df_e_2, "e2")


# Estimating simple regression #








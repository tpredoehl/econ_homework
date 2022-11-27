library(GRS.test)
library(tidyverse)
library(lubridate)
library(dplyr)
library(broom)
library(magrittr)
library(frenchdata)

#######################################################################
###
### DATA
###
#######################################################################

# load FF10 data from .Rda file
# load(file = "ff10_sample_1927_2022.Rda") # supplies ff_temp
# alternative load from FF website using FF package
ff_factors <- download_french_data('Fama/French 3 Factors')
ff_factors_monthly <- ff_factors$subsets$data[[1]]
names(ff_factors_monthly)[2] <- "Rm"

ff_size <- download_french_data('Portfolios Formed on Size')
targetSubset <- which(ff_size$subsets$name == "Value Weight Returns -- Monthly")
ff_10_size_monthly_vw <- ff_size$subsets$data[[targetSubset]]
names(ff_10_size_monthly_vw) <- gsub(pattern = " ", replacement = ".", x = names(ff_10_size_monthly_vw))
# # join data sets
decileCols <- c("Lo.10", "Dec.2","Dec.3", "Dec.4", "Dec.5", "Dec.6", "Dec.7", "Dec.8", "Dec.9", "Hi.10")
#decileCols <- c("Lo 10", "Dec 2","Dec 3", "Dec 4", "Dec 5", "Dec 6", "Dec 7", "Dec 8", "Dec 9", "Hi 10")
ff_raw <- left_join(x = ff_10_size_monthly_vw, y = ff_factors_monthly, by = "date") %>% filter(date > 192601 & date < 199901) %>% mutate(date = ym(date))

minusRF <- function(x) x-ff_raw$RF
N = 10    # small relative to T -> OK
T = 870   # large relative to N -> OK
K = 2
#######################################################################
###
### OLS
###
### multivariate OLS dplyr style - lm() uses sample variance (1/(T-1)) !!!
###
#######################################################################
ff_OLS <- ff_raw %>%
  select(all_of(decileCols), Rm) %>% 
  mutate_at(decileCols, minusRF) %>% 
  pivot_longer(cols = all_of(decileCols), names_to = "Decile", values_to = "Re") %>% 
  select(Decile, Re, Rm) %>% 
  nest(data = -Decile) %>% 
  mutate(
    ols = map(data, ~ lm(Re ~ Rm - 1, data = .x)), # run ols regression without intercept
    RE  = map(data, ~ mean(.x$Re)),                # extract E[Re_i]
    RM  = map(data, ~ mean(.x$Rm)),                # extract E[Rm]
    tidied = map(ols, tidy)
  ) %>% 
  unnest(c(tidied, RE, RM))
beta_ols <- ff_OLS %>% pull(estimate)

# excess market return
f_t <- ff_raw %>% pull(Rm)
f_t_hat <- mean(f_t); c(min(f_t), max(f_t))
V_f_t <- (1/(T)) * sum((f_t - mean(f_t))^2) # var(f_t) uses T-1
SE_f_t <- sqrt(V_f_t)

# excess asset return
r_it <- ff_raw %>% 
  mutate_at(decileCols, minusRF) %>% 
  select(all_of(decileCols)) %>% as.matrix()
mu_hat <- colMeans(r_it)

# OLS beta
X <- f_t             # [T x K], now WITH COLUMN OF ONE!!!
Sxx <- (1/T) * t(X)%*%X
Sxr <- (1/T) * t(X) %*% r_it
beta_ols <- t(solve(Sxx)%*%Sxz)

# OLS residuals
u_it <- r_it - f_t%*%t(beta_ols)
Vu2 <- (u_it)%*%t(u_it)                                         # [870x10 x 10x870 = TN x NT = T x T]
omega <- (1/T)*t(u_it)%*%(u_it) 

# OLS lambda estimate and SE
lambda_hat_ols <- solve(t(beta_ols)%*%beta_ols) %*% t(beta_ols) %*% mu_hat
BBB <- solve(t(beta_ols)%*%beta_ols) %*% t(beta_ols)
BBBt <- beta_ols %*% solve(t(beta_ols)%*%beta_ols)
V_lambda_hat_ols <- (1/T) * (V_f_t + BBB %*% omega %*% BBBt)
SE_lambda_hat_ols <- sqrt(V_lambda_hat_ols)
(lambda_hat_ols - 0)/sqrt(V_lambda_hat_ols)                     # significance? Yes, can reject H0, that lambda_hat_ols is equal to 0
(lambda_hat_ols - f_t_hat)/sqrt(V_lambda_hat_ols)               # equal to mean market return? No, cannot reject H0, that lambda_hat_ols is equal to mean market return
#######################################################################################################
###
### TimeSeries
###
### TS should be the average of the market, i.e. slope 0 -> E[Rm]/beta=1, which is 0.708
#######################################################################################################
lambda_hat_ts <- f_t_hat
u_ts <- mu_hat - lambda_hat_ts * beta_ols
omega_ts <- (1/T) * u_ts%*%t(u_ts)
V_lambda_hat_ts <- (1/T) * (V_f_t + BBB %*% omega_ts %*% BBBt) # using Ch 3, p. 59, (1.44) 
SE_lambda_hat_ts <- sqrt(V_lambda_hat_ts)
#######################################################################################################
### GLS Test
### lambda_hat_gls is identical with Cochrane, 2009, p. 282
lambda_hat_gls <- solve(t(beta_ols) %*% solve(omega) %*% beta_ols) %*% (t(beta_ols) %*% solve(omega) %*% mu_hat)
V_lambda_hat_gls <- (1/T) %*% (solve(t(beta_ols) %*% solve(omega) %*% beta_ols) + V_f_t)
SE_lambda_hat_gls <- sqrt(V_lambda_hat_gls)

# try an R package for GLS test
require(nlme)

r_it
gls(model = r_it ~ f_t - 1)

#######################################################################################################
#######################################################################################################
###
### Table 15.1
###
#######################################################################################################
TimeSeries <- c(lambda_hat_ts, SE_lambda_hat_ts)
CS_OLS <- c(lambda_hat_ols, SE_lambda_hat_ols)
CS_GLS <- c(lambda_hat_gls, SE_lambda_hat_gls)
bind_cols("TimeSeries" = TimeSeries, "OLS" = CS_OLS, "GLS" = CS_GLS)



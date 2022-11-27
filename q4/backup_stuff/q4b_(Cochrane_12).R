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

# excess asset return
r_it <- ff_raw %>%
  select(all_of(decileCols)) %>% 
  mutate_at(decileCols, minusRF) %>% 
  as.matrix()
E_ri <- colMeans(r_it)

# excess market return
f_t <- ff_raw %>% pull(Rm)
f_t_hat <- mean(f_t)
VAR_f_t <- (1/(T-1)) * sum((f_t - f_t_hat)^2)
SIG_f_t <- sqrt(VAR_f_t)    # sqrt(var(f_t)) # = 5.577758, uses sample var 1/(T-1), lecture note use pop. var (1/T)

# ts estimator
lambda_hat_ts <- mean(f_t)

# 12.3 Testing H0: all alpha = 0
GRS_test_stats <- GRS.test(ret.mat = r_it, factor.mat = f_t)
alpha_ols <- GRS_test_stats$coef[,"intercept"]
beta_ols <- GRS_test_stats$coef[,"Singlefactor"]
u_it <- GRS_test_stats$resid
omega <- (1/T) * t(u_it)%*%u_it
GRS_manual <- ((T-N-1)/N/(1 + (f_t_hat/SIG_f_t)^2)) %*% t(alpha_ols) %*% solve(omega) %*% alpha_ols # F ~ F_(N,T−N−1) = 0.8405352 = GRS test stat IIF alpha_ols is calculated !! see ff_OLS
row3 <- bind_cols(Type = "F-test", GRS = as.numeric(GRS_test_stats$GRS.stat), manual = as.numeric(GRS_manual), Diff = as.numeric(GRS_test_stats$GRS.stat - GRS_manual))
row3

require(AER)
ols_model <- lm(E_ri  ~ beta_ols - 1); summary(ols_model)
gls_model <- gls(E_ri ~ beta_ols - 1); summary(ols_model)
diag(vcov(object = ols_model))




NeweyWest(x = ols_model, lag = 3)
NeweyWest(x = gls_model, lag = 0)


data("Journals")
journals <- Journals[, c("subs", "price")]
journals$citeprice <- Journals$price/Journals$citations
journals$age <- 2022 - Journals$foundingyear
jour_lm <- lm(log(subs) ~ log(citeprice), data = journals)
# testing heteroskedasticity
bptest(jour_lm)


# OLS CS Regression, (12.11) Cochrane, p. 236
# run TS regression for each asset, THIS TIME WITHOUT intercept!!
ff_OLS_b <- ff_raw %>%
  select(all_of(decileCols), Rm) %>% 
  mutate_at(decileCols, minusRF) %>% 
  pivot_longer(cols = all_of(decileCols), names_to = "Decile", values_to = "Re") %>% 
  select(Decile, Re, Rm) %>% 
  nest(data = -Decile) %>% 
  mutate(
    ols = map(data, ~ lm(Re ~ Rm - 1, data = .x)), # NOW WITHOUT INTERCEPT
    RE  = map(data, ~ mean(.x$Re)),                # extract E[Re_i]
    RM  = map(data, ~ mean(.x$Rm)),                # extract E[Rm]
    tidied = map(ols, tidy)
  ) %>% 
  unnest(c(tidied, RE, RM))
beta_ols <- ff_OLS_b %>% filter(term == "Rm") %>% pull(estimate)
alpha_ols <- ff_OLS_b %>% filter(term == "(Intercept)") %>% pull(estimate)
beta <- cbind(alpha_ols, beta_ols)

lambda_hat_ols <- solve(t(beta_ols)%*%beta_ols)%*%t(beta_ols)%*%E_ri
a_i <- E_ri - lambda_hat_ols%*%beta_ols           # pricing error 
omega <- (1/T) * t(a_i)%*%a_i
# residual covariance matrix
e_it <- ff_OLS_b %>% 
  filter(term == "Rm") %>% 
  select(Decile, ols) %>% 
  mutate(e_i = map(ols, ~ residuals(.x))) %>% 
  unnest(e_i) %>% 
  select(-ols) %>% 
  mutate(id = rep(1:870, 10)) %>% 
  pivot_wider(names_from = Decile, values_from = e_i) %>% 
  select(-id) %>% 
  as.matrix()
E_hat <- (1/T) * t(e_it) %*% (e_it)                     # E[ee'] = 1/T Sum(ee')
E_f <- VAR_f_t <- (1/T)*t(f_t)%*%f_t
cov(f_t, f_t)
cov(e_it, e_it)
VAR_alpha_ols_b <- (1/T)*beta_ols*E_f%*%t(beta_ols) + E_hat

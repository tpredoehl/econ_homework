library(GRS.test)
library(tidyverse)
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
ff_raw <- left_join(x = ff_10_size_monthly_vw, y = ff_factors_monthly, by = "date") %>% filter(date > 192601 & date < 199901)

minusRF <- function(x) x-ff_raw$RF
N = 10    # small relative to T -> OK
T = 870   # large relative to N -> OK
K = 2
#######################################################################
###
### OLS
###
#######################################################################


#######################################################################
# multivariate OLS dplyr style - lm() uses sample variance (1/(T-1)) !!!
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
ff_OLS$estimate

TS_model <- summary(lm(ff_OLS$RE ~ ff_OLS$estimate - 1))
TS_estimate <- TS_model$coefficients[,"Estimate"]; TS_estimate
TS_SE <- TS_model$coefficients[,"Std. Error"]; TS_SE
#########################################################################################
# multivariate linear regression model, Ch 3, p. 11
# y_it, i=1..N, t=1..T
# N = 10 portfolios
# T = 870 monthly data points from 192607 to 199812
# systems of equations are stacked, Ch 3, p. 16:
#    y1 = data_Re_Lo.10, ..., y10 = data_Re_Hi.10
#    y = stack of all yi [NT x 1] = [8700 x 1]
#
#    X1 = date_Rm_Lo.10 = X2 = ... = X10 = data_Rm_Hi.10, 
#         these observations are common to all portfolios, ref (1.18) p.21
#    b1 = estimate_Lo.10, ..., b10 = estimate_Hi.10
#       = (X'X)^(-1)X'y

y <- ff_raw %>%
  select(all_of(decileCols), Rm) %>% 
  mutate_at(decileCols, minusRF) %>% 
  pivot_longer(cols = all_of(decileCols), names_to = "Decile", values_to = "Re") %>% 
  select(Decile, Re) %>% 
  arrange(match(Decile, decileCols)) %>% 
  pull(Re) %>% as.matrix()                                    # [NT x 1]

X <- ff_raw$Rm %>% as.matrix()                                # [T x 1], NO COLUMN OF ONES!!!
X_star <- diag(10)%x%X; dim(X_star)                           # [NT x NK = 8700 x 10]
ols_manual <- lm(y ~ X_star - 1)

# calculation of beta_hat Ch3, (1.19), p. 22
# direct via X_star
beta_hat <- solve(t(X_star)%*%X_star)%*%t(X_star) %*% y       # [10 x 1] beta_hat is identical to OLS
# indirect via X
XXX <- solve(t(X)%*%X)%*%t(X)
beta_hat2 <- diag(10) %x% XXX %*% y                           # [10 x 1] beta_hat is identical to OLSxx
################################################################################
###
### OLS: E[beta]
###
### manual vs smart OLS
bind_cols(Manual1 = beta_hat, manual2 = beta_hat2, dplyr = ff_OLS$estimate, Diff = beta_hat - ff_OLS$estimate)
### dplyr style E[beta] is -very much- identical with the manual beta_hat
###
################################################################################


# extract the var-cov matrix omega from the regression residuals (1.16), p. 17 & 24
# compare to ols dplyr style
# u: N x T, [10 x 870]
# omega: N x N, w_ij = 1/T sum(u_it * u_jt)
# Vu: NT x NT
u <- y - X_star %*% beta_hat                                  # [8700 x 1]
u <- t(matrix(data = u, nrow = 870, ncol = 10))               # [870 x 10]
omega <- (1/(T)) * u%*%t(u)                                   # [10 x 10]
Vu <- omega %x% diag(870)                                     # [8700 x 8700]
XXX_star <- solve(t(X_star)%*%X_star)%*%t(X_star)             # [10 x 8700]
var_beta <- XXX_star %*% solve(Vu) %*% t(XXX_star)            # [10 x 10]
se_beta <- (sqrt(var_beta))                                   # [10 x 10]
################################################################################
###
### OLS: V[beta]
###
### manual vs smart OLS
bind_cols(Manual = diag(se_beta), dplyr = ff_OLS$std.error, Diff = diag(se_beta)-ff_OLS$std.error)
### dplyr style Var(beta) is -almost- identical with the manual beta
### this proofs (1.27), p. 27
### if omega is calculated as the sample V(u) i.e. dividing by T-1 instead of T, 
### then the difference almost disappears.
### Consequently, it is shown that ols package in R (lm) uses the sample variance.
###
################################################################################

#######################################################################
###
### GRS
###
#######################################################################

#########################################################################################
# F-Test
# Asset excess return
z_it <- ff_raw %>%
  select(all_of(decileCols)) %>% 
  mutate_at(decileCols, minusRF)
# excess return of the risky assets
mu_i <- z_it %>% 
  summarise(across(everything(), mean))
var(z_it)

# market risk premium
z_mt <- ff_raw %>% 
  select(Rm)
mu_m <- z_mt %>% summarise(across(everything(), mean))
var_m <- var(z_mt)
sig_m <- sqrt(var_m)

# max theoretical sharpe ratio
S_rm <- mu_m/sig_m

# t-test manual (1.6), p. 7
sig_im <- cov(z_it, z_mt)
beta_i <- sig_im %x% solve(var(z_mt))
X <- bind_cols(One = 1, z_mt) %>% as.matrix()             # [T x K], now WITH COLUMN OF ONE!!!
Sxx <- (1/T) * t(X)%*%X %>% as.matrix(); class(Sxx)
Sxz <- (1/T) * t(X) %*% (z_it %>% as.matrix()); class(Sxz)
beta <- solve(Sxx)%*%Sxz; class(beta)
a_i <- beta[1,]
b_i <- beta[2,]                                           # same beta as by OLS, X WITHOUT ONES
u_it <- z_it - X%*%beta
u_it <- as.matrix(u_it)
# omega, p. 17
Vu2 <- (u_it)%*%t(u_it)                                   # [870x10 x 10x870 = TN x NT = T x T]
omega2 <- (1/T)*t(u_it)%*%(u_it)                          # identical to omega based on u, 
                                                          # based on X, based on Rm alone, 
                                                          # i.e. no column of ONEs!!!
w_i2 <- diag(omega2)
# the t-statistics to test H0: a_i = 0
eta1 <- as.matrix(a_i) * as.numeric(1/sqrt((w_i2/T)*(1+mu_m^2/sig_m^2)))
t.test(x = X, y = z_it[,1])                               # cannot reject H0 that a_i = 0,
                                                          # t-stat vastly different from eta1???
  
# F-Test, p.28 ff
f_t <- z_mt %>% as.matrix()                               # [870 x 1]
f_hat <- mu_m %>% as.numeric()                            # [1 x 1]
y_t <- z_it %>% as.matrix()
y_hat <- mu_i %>% as.numeric()                            # [10 x 1]
Eff <- (1/T)*t(f_t - f_hat)%*%(f_t - f_hat)               # [1x870 x 870x1]   -> [1x1]
                                                          # =(T-1)/T*cov(x = f_t, f_t)
Eff <- cov(x = f_t, f_t)
Eyf <- (1/T)*t(y_t - y_hat)%*%(f_t - f_hat)               # [870x10 x 1x870]  -> [1x1]
Eyy <- (1/T)*sum(t(y_t - y_hat)%*%(y_t - y_hat))          # [870x10 x 1x870]  -> [1x1]


# test if XX^-1[1,1] = 1 + f_hat' Eff^-1 f_hat -> CORRECT
solve(t(X)%*%X)[1,1] - (1/T) * (1+t(f_hat)%*%Eff^-1%*%f_hat)

# finite sample F-test for H0: a1=a2=...=aN = 0
eta_2 <- (T-N-(K-1))/N * (1+t(f_hat)%*%Eff^-1%*%f_hat)^-1 * (t(a_i)%*%solve(omega2)%*%a_i)

# GRS Test (finite sample)
GRS_test_stats <- GRS.test(ret.mat = y_t, factor.mat = f_t)

################################################################################
###
### GRS
### checking GRS stats against manual calcs
u_it - GRS_test_stats$resid
row1 <- bind_cols(Type = "ai", GRS = GRS_test_stats$coef[,1], manual = a_i, Diff = GRS_test_stats$coef[,1] - a_i)
row2 <- bind_cols(Type = "bi", GRS = GRS_test_stats$coef[,2], manual = b_i, Diff = GRS_test_stats$coef[,2] - b_i)
row3 <- bind_cols(Type = "F-test", GRS = as.numeric(GRS_test_stats$GRS.stat), manual = as.numeric(eta_2), Diff = as.numeric(GRS_test_stats$GRS.stat - eta_2))
bind_rows(row3, row2, row1)
### GRS.Test package uses Eff = cov(x = f_t, f_t), 
### which removes the bias by dividing (1/(T-1)) instead of (1/T) as per p. 26

# max Sharpe Ratio, p. 33
# ex-post maximum squared Sharpe Ratio achievable by a mean-variance investor investing in the K factors
theta_p2 <- t(f_hat)%*%Eff^-1%*%f_hat               
# ex-post maximum squared Sharpe Ratio achievable by a mean-variance investor investing in the K factors
# (assuming that they are portfolios!)  and the N test assets, which have returns yt
theta_star_p2 <- (t(a_i)%*%solve(omega2)%*%a_i) + theta_p2


#######################################################################################################
# CS test CAPM, p. 54ff
# mu_hat = beta lambda + e
# E[r_it] = beta E[lambda] + E[e], p. 59
r_it <- z_it
f_t_hat <- mean(f_t)
V_f_t <- (1/(T)) * sum((f_t - mean(f_t))^2) # var(f_t) uses T-1
SE_f_t <- sqrt(V_f_t)
beta_i
u_it
mu_hat <- as.numeric(r_it %>% summarise(across(everything(), mean)))
e_star <- beta_i * mean(f_t - as.numeric(mu_m)) + colMeans(u_it) # E[e_star] = 0
V_e_star <- E <- (1/T) * (beta_i %*% V_f_t %*% t(beta_i) + omega2)
lambda_hat_ols <- solve(t(beta_i)%*%beta_i) %*% t(beta_i) %*% as.numeric(mu_hat)
BBB <- solve(t(beta_i)%*%beta_i) %*% t(beta_i)
BBBt <- beta_i %*% solve(t(beta_i)%*%beta_i)
V_lambda_hat_ols <- (1/T) * (V_f_t + BBB %*% omega2 %*% BBBt)
SE_lambda_hat_ols <- sqrt(V_lambda_hat_ols)
(lambda_hat_ols - 0)/sqrt(V_lambda_hat_ols)                     # significance? Yes, can reject H0, that lambda_hat_ols is equal to 0
(lambda_hat_ols - f_t_hat)/sqrt(V_lambda_hat_ols)               # equal to mean market return? No, cannot reject H0, that lambda_hat_ols is equal to mean market return

#######################################################################################################
# GLS Test
lambda_hat_gls <- solve(t(beta_i) %*% solve(omega2) %*% beta_i) %*% t(beta_i) %*% solve(omega2) %*% mu_hat
# lambda_hat_gls is identical with Cochrane, 2009, p. 282
V_lambda_hat_gls <- (1/T) %*% (solve(t(beta_i) %*% solve(omega2) %*% beta_i) + V_f_t)
SE_lambda_hat_gls <- sqrt(V_lambda_hat_gls)
e_hat_gls <- mu_hat - beta_i%*%lambda_hat_gls
V_e_hat_gls <- (1/T) * (omega2 - beta_i %*%  solve(t(beta_i) %*% solve(omega2) %*% beta_i) %*% t(beta_i))


# Alternative Wald test for e_i = 0, p. 68 (see Cochrane, 2009, 12.2)
chisq_stat <- T %*% t(e_hat_gls) %*% solve(omega2) %*% e_hat_gls # is distributed Chi2(N-1)
pchisq(q = chisq_stat, df = N-1)

#######################################################################################################
# replicate table 15.1
# TS estimate -> mean market return in pct per month
# CS estimate -> lambda

# tbl 15.1    TS            CS-OLS            CS-GLS
# estimate    TS_estimate   lambda_hat_ols    lambda_hat_gls
# error iid   TS_SE         SE_lambda_hat_ols SE_lambda_hat_gls
# 0-lag
# 3-lag NW
# 24-lag


# TimeSeries E[z_it]  = lambda_hat_ts beta_i + E[u_it]
#            E[f_t]   = lambda_hat_ts
#             u_it    = z_it - lambda_hat_ts * beta_i
###############################################################
# stuff
###############################################################
# Alternative Wald test for e_i = 0, p. 68 (see Cochrane, 2009, 12.2)
chisq_stat <- T %*% t(e_hat_gls) %*% solve(omega_ols) %*% e_hat_gls # is distributed Chi2(N-1)
pchisq(q = chisq_stat, df = N-1)


# 12.3 Testing H0: all alpha = 0
GRS_test_stats <- GRS.test(ret.mat = r_it, factor.mat = f_t)
alpha_grs <- GRS_test_stats$coef[,"intercept"]
beta_grs <- GRS_test_stats$coef[,"Singlefactor"]
u_it <- GRS_test_stats$resid
omega <- (1/T) * t(u_it)%*%u_it
GRS_manual <- ((T-N-1)/N/(1 + (f_t_hat/SIG_f_t)^2)) %*% t(alpha_grs) %*% solve(omega) %*% alpha_grs # F ~ F_(N,T−N−1) = 0.8405352 = GRS test stat IIF alpha_ols is calculated !! see ff_OLS
GRS_Test_Comparison <- bind_cols(Type = "F-test", GRS = as.numeric(GRS_test_stats$GRS.stat), manual = as.numeric(GRS_manual), Diff = as.numeric(GRS_test_stats$GRS.stat - GRS_manual))
GRS_Test_Comparison

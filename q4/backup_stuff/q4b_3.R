#######################################################################
###
### SETUP
###
#######################################################################
require(sandwich)
library(GRS.test)
library(dplyr)
library(tidyr)
library(lubridate)
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
# function to deduct Rf from input vector
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
  f_t <- ff_raw %>% select(Rm) %>% as.matrix(); class(f_t)
  f_t_hat <- mean(f_t)
  VAR_f_t <- (1/(T-1)) * sum((f_t - f_t_hat)^2)
  SIG_f_t <- sqrt(VAR_f_t)    # sqrt(var(f_t)) # = 5.577758, uses sample var 1/(T-1), lecture note use pop. var (1/T)

###############################################################
# LAMBDA - TS Estimator
###############################################################
  # 1-pass: find alphas, betas
    TS_model_1 <- lm(r_it~f_t)                                          
    beta_TS <- t(TS_model_1$coefficients)
  # 2-pass: find lambdas by OLS
    TS_model_2 <- lm((E_ri)~beta_TS-1)               
    TS_model_2_s <- summary(TS_model_2)
    lambda_hat_ts <- TS_model_2_s$coefficients[2,"Estimate"]
    omega_TS <- (1/T)*t(TS_model_1$residuals)%*%TS_model_1$residuals
    V_lambda_hat_TS <- (1/T) * (solve(t(beta_TS) %*% solve(omega_TS) %*% beta_TS) + VAR_f_t)
    SE_lambda_hat_ts <- sqrt(V_lambda_hat_TS)[2,2]
  # iid
    omega_TS_iid <- mean(diag(omega_TS))*diag(1,nrow = N)
    V_lambda_hat_TS_iid <- (1/T)*(solve(t(beta_TS) %*% solve(omega_TS_iid) %*% beta_TS) + VAR_f_t)
    SE_lambda_hat_TS_iid <- sqrt(V_lambda_hat_TS_iid)[2,2]
  # 0-lag
    omega_TS_0lag <- diag(diag(omega_TS),nrow = N)
    V_lambda_hat_TS_0lag <- (1/T)*(solve(t(beta_TS) %*% solve(omega_TS_0lag) %*% beta_TS) + VAR_f_t)
    SE_lambda_hat_TS_0lag <- sqrt(V_lambda_hat_TS_0lag)[2,2]
  # 3NW-lag
    omega_TS_3NW <- NeweyWest(TS_model_1, lag = 0, sandwich = FALSE)
    f_Cols <- grep(pattern = ":f_t", x = colnames(omega_TS_3NW))
    omega_TS_3NW <- omega_TS_3NW[f_Cols,f_Cols]
    V_lambda_hat_TS_3lag <- (1/T)*(solve(t(beta_TS) %*% solve(omega_TS_3NW) %*% beta_TS) + VAR_f_t)
    SE_lambda_hat_TS_3lag <- sqrt(V_lambda_hat_TS_3lag)[2,2]
  # 24NW-lag
    omega_TS_24lag <- NeweyWest(TS_model_1, lag = 24)
    f_Cols <- grep(pattern = ":f_t", x = colnames(omega_TS_24lag))
    omega_TS_24lag <- omega_TS_24lag[f_Cols,f_Cols]
    V_lambda_hat_TS_24lag <- (1/T)*(solve(t(beta_TS) %*% solve(omega_TS_24lag) %*% beta_TS) + VAR_f_t)
    SE_lambda_hat_TS_24lag <- sqrt(V_lambda_hat_TS_24lag)[2,2]
  
###############################################################
# LAMBDA - OLS Estimator
###############################################################
  # 1-pass estimate betas for each portfolio
    OLS_model <- lm(r_it~f_t-1)                                             # set alpha=0, residuals capture pricing error
    omega_ols <- cov(OLS_model$residuals)                                   # (1.33) lecture notes 3, p. 44
    beta_ols <- t(OLS_model$coefficients)
    lambda_hat_ols <- solve(t(beta_ols)%*%beta_ols)%*%t(beta_ols) %*% E_ri  # (1.43) lecture notes 3, p. 59
    BBB <- solve(t(beta_ols)%*%beta_ols) %*% t(beta_ols)
    BBBt <- beta_ols %*% solve(t(beta_ols)%*%beta_ols)
    V_lambda_hat_ts <- (1/T)*(VAR_f_t+BBB%*%omega_ols%*%BBBt)               # (1.44) lecture notes 3, p. 59
    SE_lambda_hat_ts <- sqrt(V_lambda_hat_ts)                           
  # iid
    omega_ols_iid <- mean(diag(omega_ols))*diag(1,nrow = N)
    V_lambda_hat_ols_iid <- (1/T)*(VAR_f_t+BBB%*%omega_ols_iid%*%BBBt)      # (1.44) lecture notes 3, p. 59
    SE_lambda_hat_ols_iid <- sqrt(V_lambda_hat_ols_iid)
  # 0-lag
    omega_ols_0lag <- diag(diag(omega_ols),nrow = N)
    V_lambda_hat_ols_0lag <- (1/T)*(VAR_f_t+BBB%*%omega_ols_0lag%*%BBBt)
    SE_lambda_hat_ols_0lag <- sqrt(V_lambda_hat_ols_0lag)
  # 3NW-lag
    omega_ols_3NW <- NeweyWest(OLS_model, lag = 3, sandwich = TRUE)
    f_Cols <- grep(pattern = ":f_t", x = colnames(omega_ols_3NW))
    omega_ols_3NW <- omega_ols_3NW[f_Cols,f_Cols]
    V_lambda_hat_ols_3NW <- (1/T)*(VAR_f_t+BBB%*%omega_ols_3NW%*%BBBt)
    SE_lambda_hat_ols_3NW <- sqrt(V_lambda_hat_ols_3NW)
  # 24-NW-lag
    omega_ols_24lag <- vcovHC(OLS_model, lag = 24)
    f_Cols <- grep(pattern = ":f_t", x = colnames(omega_ols_24lag))
    omega_ols_24lag <- omega_ols_24lag[f_Cols,f_Cols]
    V_lambda_hat_ols_24lag <- (1/T)*(VAR_f_t+BBB%*%omega_ols_24lag%*%BBBt)
    SE_lambda_hat_ols_24lag <- sqrt(V_lambda_hat_ols_24lag)
###############################################################
# LAMBDA - GLS Estimator
###############################################################
  # GLS Test
  # 1-pass estimate betas for each portfolio
    GLS_model <- lm(r_it~f_t)                                         # set alpha=0, residuals capture pricing error
    beta_gls <- t(GLS_model$coefficients)
  # 2-pass: estimate lambda
    lambda_hat_gls_model <- lm(E_ri~beta_gls)
    u_gls <- r_it - cbind(1,f_t)%*%t(beta_gls)                        # extract model residuals
    omega_gls <- (1/T) * t(u_gls)%*%(u_gls) 
    lambda_hat_gls <- solve(t(beta_gls)%*%solve(omega_gls)%*%beta_gls)%*%t(beta_gls)%*%solve(omega_gls)%*%E_ri
    lambda_hat_gls <- as.numeric(lambda_hat_gls[2,1])
  # lambda_hat_gls is identical with Cochrane, 2009, p. 282
    V_lambda_hat_gls <- (1/T) * (solve(t(beta_gls) %*% solve(omega_gls) %*% beta_gls) + VAR_f_t)
    SE_lambda_hat_gls <- sqrt(V_lambda_hat_gls)[2,2]
  # iid
    omega_gls_iid <- mean(diag(omega_gls))*diag(1,nrow = N)
    V_lambda_hat_gls_iid <- (1/T)*(solve(t(beta_gls) %*% solve(omega_gls_iid) %*% beta_gls) + VAR_f_t)
    SE_lambda_hat_gls_iid <- sqrt(V_lambda_hat_gls_iid)[2,2]
  # 0-lag
    omega_gls_0lag <- diag(diag(omega_gls),nrow = N)
    V_lambda_hat_gls_0lag <- (1/T)*(solve(t(beta_gls) %*% solve(omega_gls_0lag) %*% beta_gls) + VAR_f_t)
    SE_lambda_hat_gls_0lag <- sqrt(V_lambda_hat_gls_0lag)[2,2]
  # 3NW-lag
    omega_gls_3NW <- NeweyWest(GLS_model, lag = 3)
    f_Cols <- grep(pattern = ":f_t", x = colnames(omega_gls_3NW))
    omega_gls_3NW <- omega_gls_3NW[f_Cols,f_Cols]
    V_lambda_hat_gls_3NW <- (1/T)*(solve(t(beta_gls) %*% solve(omega_gls_3NW) %*% beta_gls) + VAR_f_t)
    SE_lambda_hat_gls_3NW <- sqrt(V_lambda_hat_gls_3NW)[2,2]
  # 24NW-lag
    omega_gls_24lag <- NeweyWest(GLS_model, lag = 24)
    f_Cols <- grep(pattern = ":f_t", x = colnames(omega_gls_24lag))
    omega_gls_24lag <- omega_gls_24lag[f_Cols,f_Cols]
    V_lambda_hat_gls_24lag <- (1/T)*(solve(t(beta_gls) %*% solve(omega_gls_24lag) %*% beta_gls) + VAR_f_t)
    SE_lambda_hat_gls_24lag <- sqrt(V_lambda_hat_gls_24lag)[2,2]
###############################################################
# Table 15.1
###############################################################
  row_Estimate <- c(TS=lambda_hat_ts, OLS=lambda_hat_ols, GLS=lambda_hat_gls)
  row_iid <- c(TS=SE_lambda_hat_TS_iid, OLS=SE_lambda_hat_ols_iid, GLS=SE_lambda_hat_gls_iid)
  row_0lag <- c(TS=SE_lambda_hat_TS_0lag, OLS=SE_lambda_hat_ols_0lag, GLS=SE_lambda_hat_gls_0lag)
  row_3_NW <- c(TS=SE_lambda_hat_TS_0lag, OLS=SE_lambda_hat_ols_0lag, GLS=SE_lambda_hat_gls_0lag)
  row_24lag <- c(TS=SE_lambda_hat_TS_24lag, OLS=SE_lambda_hat_ols_24lag, GLS=SE_lambda_hat_gls_24lag)
  
  tbl_15_1 <- bind_rows(
    row_Estimate,
    row_iid,
    row_0lag,
    row_3_NW,
    row_24lag
  ) %>% 
    mutate(
      TS = round(TS,4),
      Title = c("Estimate", "iid", "0-lag", "3l-NW", "24l")
    ) %>% 
    select(Title, TS, OLS, GLS)
  tbl_15_1
  
  
  
  vcov(OLS_model, type = "const")
  u <- OLS_model$residuals
  cov(u)                                                                  # Zeileis (3): Vu <- (diag(1,T)-X%*%solve(t(X)%*%X)%*%t(X))%*%r_it, 
  sig2_hat_u <- 1/(N-1) * sum(diag(u)^2)
  phi <- sig2_hat_u*solve(t(X)%*%X)
  
  
  data("Investment")  
  fm.inv <- lm(RealInv ~ RealGNP + RealInt, data = Investment)
  summary(fm.inv)
  u <- fm.inv$residuals
  u%*%t(u)
  X <- as.matrix(cbind(1, as.numeric(Investment[,"RealGNP"]), as.numeric(Investment[,"RealInt"])))[-1,]
  n <- nrow(X)
  y_hat <- as.numeric(Investment[,"RealGNP"])[-1]
  Mx <- diag(n) - X%*%solve(t(X)%*%X)%*%t(X)
  u_hat <- Mx%*%y_hat
  sig2_hat_u <- 1/(n-3) * sum(u_hat^2)
  omega <- sig2_hat_u*diag(n)
  sig_hat <- sqrt(sig2_hat_u)
  phi_hat <- sig_hat*solve(t(X)%*%X)

  vcovHC(fm.inv, type = "const", sandwich = FALSE)

  
  
    NeweyWest(fm.inv, lag = 4, prewhite = FALSE)
  sqrt(0.0004229789)

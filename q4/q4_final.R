#######################################################################
###
### SETUP
###
#######################################################################
  library(tidyverse)
  require(sandwich)
  library(magick)
  library(GRS.test)
  library(lubridate)
  library(kableExtra)
  library(ggplot2)
  library(dplyr)
  library(broom)
  library(magrittr)
  library(frenchdata)
#######################################################################
###
### Helper Functions
###
#######################################################################
  get15_2 <- function(ols_model, type, uDist, finiteSample) {
  # ChiSq is the asymptotic test. 
  # The finite sample F-test requires multiplier (T-N-(K-1))/N instead of T
  # finiteSample = FALSE
  # type = "GLS"
  # uDist = "3NW"
  # Dimensions
  N                         <- ncol(ols_model$fitted.values)    # small relative to T -> OK
  T                         <- nrow(ols_model$residuals)        # large relative to N -> OK
  K                         <- ols_model$rank
  beta_vec                  <- t(ols_model$coefficients)
  alpha                     <- ols_model$coefficients["(Intercept)",]
  beta                      <- ols_model$coefficients["f_t",]
  omega                     <- cov(ols_model$residuals)
  # params
  f_t                       <- ols_model$model$f_t
  f_t_hat                   <- mean(f_t)
  V_f                       <- as.numeric(var(f_t))
  r_it                      <- ols_model$model$r_it
  # lambda
  if(type == "TS") lambda   <- mean(f_t)
  if(type == "OLS")lambda   <- solve(t(beta)%*%beta)%*%t(beta)%*%E_ri
  if(type == "GLS")lambda   <- (solve(t(beta_vec)%*%solve(omega)%*%beta_vec)%*%t(beta_vec)%*%solve(omega)%*%E_ri)["f_t",1]
  # omega
  if(uDist == "iid")   omega <- omega # mean(diag(omega))*diag(1,nrow = N) gives too low 
  if(uDist == "0lag")  omega <- omega # diag(diag(omega),nrow = N)
  if(uDist == "3NW") {
    omega                   <- NeweyWest(ols_model, lag = 3, prewhite = FALSE, sandwich = FALSE, adjust = TRUE)
    f_Cols                  <- grep(pattern = ":f_t", x = colnames(omega))
    omega                   <- omega[f_Cols,f_Cols]*as.numeric(solve(t(f_t)%*%f_t))
  }  
  if(uDist == "24lag") {
    omega                   <- NeweyWest(ols_model, lag = 24, prewhite = FALSE, sandwich = TRUE, adjust = TRUE)
    f_Cols                  <- grep(pattern = ":f_t", x = colnames(omega))
    omega                   <- omega[f_Cols,f_Cols]
  }
  
  # testing
  theta2_p    <- t(f_t_hat)%*%solve(V_f)%*%f_t_hat
  SHK         <- as.numeric((1+t(lambda%*%solve(V_f)%*%lambda)))
  aOa         <- (t(alpha)%*%solve(omega)%*%alpha)
  if (type == "TS" & finiteSample == TRUE) {
    eta       <- (T-N-(K-1))/N * (1+theta2_p)^-1 * aOa      # TS/SUR (finite sample) Ch3, p.29 (GRS Stat)
    eta_p     <- pf(q = eta, df1 = N, df2 = T-N-(K-1)) 
  }
  if (finiteSample == FALSE) {
    eta       <- T * (1+theta2_p)^-1 * aOa                  # TS/SUR (asymptotic) Ch3, p.31, (1.29)
    if(type == "GLS" | type == "OLS") N <- N-1              # comment Cochrane, p. 286: "the market premium is estimated from the cross section rather than from the market return"
    eta_p     <- pchisq(q = eta, df = N)
  }
  eta_SHK     <- SHK*eta
  eta_p_SHK   <- pchisq(q = eta_SHK, df = N-1)
  result      <- c(eta, eta_p, eta_SHK, eta_p_SHK)
  
  return(result)
}
  getNW <- function(L, X, e){
    X <- cbind(f_t)
    e <- u
    K = ncol(e)
    N = nrow(X)
    L = 3
    XSX <- matrix(data = 0, nrow = K, ncol = K)
    for (i in 0:L){
      wl <- 1-(i/(L+1))
      Xi <- X[1:(N-i),]
      Xj <- X[(1+i):N,]
      ei <- e[1:(N-i),]
      ej <- e[(1+i):N,]
      Wi <- (1/N)*(as.numeric(t(Xi)%*%(Xj))*t(ei)%*%(ej) + as.numeric(t(Xj)%*%(Xi))*t(ej)%*%(ei))
      W_temp <- wl*Wi
      XSX <- XSX+W_temp
    }
    XSX <- N/(N-K)*XSX
    XX <- as.numeric(solve(t(X)%*%X))
    Vb_HAC <- XX*XSX*XX
    Vb_HAC[1,1]
    Vb_NW <- NeweyWest(ols_model, lag = 3, prewhite = FALSE, sandwich = TRUE, adjust = FALSE)
    f_Cols <- grep(pattern = ":f_t", x = colnames(Vb_NW))
    Vb_NW <- Vb_NW[f_Cols,f_Cols]
    Vb_NW[1,1]
    
    Vb_NW[1,1]/Vb_HAC[1,1]
    
    # NO MATCH ! NO idea  ....
  }
  getTBL <- function(lambda, beta, type){
    iid   <- getSE(lambda, beta, omega_iid, V_f, type)
    lag0  <- getSE(lambda, beta, omega_0lag, V_f, type)
    NW3   <- getSE(lambda, beta, omega_3NW, V_f, type)
    lag24 <- getSE(lambda, beta, omega_24lag, V_f, type)
    title <- c("Estimate", "iid", "0-lag", "3 lags, NW", "24 lags")
    tbl_ <- t(bind_rows(Estimate=c(lambda, NA), iid=iid, lag0=lag0, NW3=NW3, lag24=lag24))
    tbl_lambda <- bind_cols(title, tbl_)
    colnames(tbl_lambda) <- c("Item", paste0(type,"_EstSE"), paste0(type, "_t"))
    return(tbl_lambda)
  }
  getSE <- function(lambda, beta_vec, omega, V_f, type){
    if(type == "ts" || type == "gls"){
      # omega <- omega_24lag
      V_ <- (1/T) * (solve(t(beta_vec) %*% solve(omega) %*% beta_vec) + V_f)
      if(!is.null(dim(beta_vec))) SE_ <- sqrt(V_)[2,2]
      if(is.null(dim(beta_vec))) SE_ <- sqrt(V_)  
    }
    if(type == "ols"){
      BBB <- solve(t(beta)%*%beta) %*% t(beta)
      BBBt <- beta %*% solve(t(beta)%*%beta)
      V_ <- (1/T)*(V_f+BBB%*%omega%*%BBBt)               # (1.44) lecture notes 3, p. 59
      SE_ <- sqrt(V_) 
    }
    t_ <- lambda/SE_
    return(c(SE_, t_))
  }
  getdmatrix <- function(N, X){
    d_11 <- diag(N)
    d_12 <- d_21 <- diag(N)*mean(X)
    d_22 <- as.numeric((1/(nrow(X)))*t(X)%*%X)*diag(N)
    d <- rbind(cbind(d_11, d_12), cbind(d_21, d_22))
    return(d)
  }
  getSMatrix <- function(X, e, L){                         
    N <- ncol(e)
    T <- nrow(X)
    K <- ncol(X)
    S_11 <- S_12 <- S_21 <- S_22 <- matrix(0,N,N)      # ref: Cochrane (12.7), p.234
    for (i in 0:L){                                 # set L=T-2 for full summation
      S_11_temp <- t(e[1:(T-i),])%*%e[(1+i):T,]
      S_12_temp <- t(mean(X[(1+i):(T)])*e[(1+i):(T),])%*%e[(1):(T-i),] # still assuming cov(f_t, e)=0
      S_21_temp <- t(e[1:(T-i),])%*%(e[(1+i):T,]*mean(X[(1):(T-i)])) # still assuming cov(f_t, e)=0
      S_22_temp <- t(mean(X[(1+i):(T)])*e[1:(T-i),])%*%(e[(1+i):T,]*mean(X[(1):(T-i)])) # still assuming cov(f_t, e)=0
      S_11 <- S_11 + (1/T)*S_11_temp
      S_12 <- S_12 + (1/T)*S_12_temp
      S_21 <- S_21 + (1/T)*S_21_temp
      S_22 <- S_22 + (1/T)*S_22_temp
    }
    S <- as.matrix(rbind(cbind(S_11, S_12), cbind(S_21, S_22)))
    return(S)
  }
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
# # join & filter data sets
  decileCols <- c("Lo.10", "Dec.2","Dec.3", "Dec.4", "Dec.5", "Dec.6", "Dec.7", "Dec.8", "Dec.9", "Hi.10")
  ff_raw <- left_join(x = ff_10_size_monthly_vw, y = ff_factors_monthly, by = "date") %>% filter(date > 192601 & date < 199901) %>% mutate(date = ym(date))
# function to deduct Rf from input vector
  minusRF <- function(x) x-ff_raw$RF
# excess asset return
  r_it <- ff_raw %>%
    select(all_of(decileCols)) %>% 
    mutate_at(decileCols, minusRF) %>% 
    as.matrix()
  E_ri <- colMeans(r_it)
# excess market return & variance of the factor
  f_t <- ff_raw %>% select(Rm) %>% as.matrix()
  f_t_hat <- colMeans(f_t)
  f_ <- f_t - f_t_hat
  sigma2_f <- t(f_)%*%(f_)/(T-1)
  V_f <- as.numeric(var(f_t))
#######################################################################
###
### 1st stage LRM
###
#######################################################################
  ols_model <- lm(r_it~f_t); ols_model_summary <- summary(ols_model)
# Dimensions
  N = ncol(ols_model$fitted.values)    # small relative to T -> OK
  T = nrow(ols_model$residuals)    # large relative to N -> OK
  K = ols_model$rank
# coefficients
  alpha <- coef(ols_model)["(Intercept)",]
  beta <- coef(ols_model)["f_t",]
  beta_vec <- t(coef(ols_model))
# var-cov matrix of (zero mean) residuals
  u <- ols_model$residuals
  omega <- (t(u)%*%u)/(T-K)

#######################################################################
###
### 4b. a replicate figures 15.1. 15.2
###
#######################################################################
  # estimate lambda by OLS and GLS
  lambda_hat_ols <- solve(t(beta)%*%beta)%*%t(beta) %*% E_ri  # (1.43) lecture notes 3, p. 59
  lambda_hat_gls <- (solve(t(beta_vec)%*%solve(omega)%*%beta_vec)%*%t(beta_vec)%*%solve(omega)%*%E_ri)[2,1]
  
  # QUESTION: why is gls using intercept and ols is not?
  
  plot_data <- bind_cols(E_ri=E_ri, beta=beta)
  plot_15_2 <- ggplot(plot_data, aes(x=beta, y=E_ri)) +
    guides(shape = guide_legend()) +
    geom_point(aes(shape = "Test assets"), size = 2) + 
    geom_point(aes(x=1, y=mean(f_t), shape = "Market"), size = 4) +
    geom_abline(aes(slope=lambda_hat_gls, intercept=0, colour="GLS"), linetype = "dashed", size=1) +
    geom_abline(aes(slope=lambda_hat_ols, intercept=0, colour="OLS"), size=1) +
    scale_x_continuous(limits = c(0,1.5)) +
    scale_y_continuous(limits = c(0,1.2)) +
    labs(
      x = "betas",
      y = "Mean mothly excess return %") +
    guides(
      shape=guide_legend(title="Portfolio"),
      colour=guide_legend(title="Regression")
    )
  
  plot_15_1 <- ggplot(plot_data, aes(x=beta, y=E_ri)) +
    guides(shape = guide_legend()) +
    geom_point(aes(shape = "Test assets"), size = 2) + 
    geom_point(aes(x=1, y=mean(f_t), shape = "Market"), size = 4) +
    geom_abline(aes(slope=mean(f_t), intercept=0, colour="TS"), linetype = "dashed", size=1) +
    scale_x_continuous(limits = c(0,1.5)) +
    scale_y_continuous(limits = c(0,1.2)) +
    labs(
      x = "betas",
      y = "Mean mothly excess return %") +
    guides(
      shape=guide_legend(title="Portfolio"),
      colour=guide_legend(title="Regression")
    )
  
  # par(mar = c(4, 4, .1, .1))
  # plot_15_1
  # plot_15_2
  ggsave(filename = "q4/replicated_15_1.png", plot = plot_15_1)
  ggsave(filename = "q4/replicated_15_2.png", plot = plot_15_2)
#######################################################################
###
### 4b. a replicate figures 15.1. 15.2
###
#######################################################################
  #######################################################################
  ###
  ### Estimates of lambda
  ###
  #######################################################################
  # estimate lambda
    lambda_hat_ols  <- solve(t(beta)%*%beta)%*%t(beta) %*% E_ri  # (1.43) lecture notes 3, p. 59
    TS_model        <- summary(lm((E_ri)~beta_vec-1))
    lambda_hat_ts   <- TS_model$coefficients[2,"Estimate"]
    lambda_hat_gls  <- (solve(t(beta_vec)%*%solve(omega)%*%beta_vec)%*%t(beta_vec)%*%solve(omega)%*%E_ri)[2,1]
  #######################################################################
  ###
  ### Estimates of SE of lambda -> Table 15.1
  ###
  #######################################################################
  # iid 
    omega_iid   <- mean(diag(omega))*diag(1,nrow = N)
  # 0-lag
    omega_0lag  <- diag(diag(omega),nrow = N)
  # 3NW-lag
  # sandwhich package outputs the whole sandwhich 
  # since we need to plug it into the V(b) = (1/T)(Vf + bbb'Obbb), the reults will look like
  # Var(b) = (1/T)(Vf + bbb' [(1/n X'X)^-1] [1/n phi] [(1/n X'X)^-1] bbb)
    omega_3NW   <- NeweyWest(ols_model, lag = 3, prewhite = FALSE, sandwich = TRUE, adjust = TRUE)
    f_Cols      <- grep(pattern = ":f_t", x = colnames(omega_3NW))
    omega_3NW   <- omega_3NW[f_Cols,f_Cols]
  # 24NW-lag
    omega_24lag <- NeweyWest(ols_model, lag = 24, prewhite = FALSE, sandwich = TRUE, adjust = TRUE)
    f_Cols      <- grep(pattern = ":f_t", x = colnames(omega_24lag))
    omega_24lag <- omega_24lag[f_Cols,f_Cols]
  # TS
  #getTBL(lambda_hat_ts, beta_vec, "ts")
    TS_iid      <- getSE(lambda_hat_ts, beta_vec, omega_iid, V_f, "ts")
    TS_0lag     <- getSE(lambda_hat_ts, beta_vec, omega_0lag, V_f, "ts")
    TS_3NW      <- getSE(lambda_hat_ts, beta_vec, omega_3NW, V_f, "ts")
    TS_24lag    <- getSE(lambda_hat_ts, beta_vec, omega_24lag, V_f, "ts")
  # OLS
  #getTBL(lambda_hat_ols, beta, "ols")
    OLS_iid     <- getSE(lambda_hat_ols, beta, omega_iid, V_f, "ols")
    OLS_0lag    <- getSE(lambda_hat_ols, beta, omega_0lag, V_f, "ols")
    OLS_3NW     <- getSE(lambda_hat_ols, beta, omega_3NW, V_f, "ols")
    OLS_24lag   <- getSE(lambda_hat_ols, beta, omega_24lag, V_f, "ols")
  # GLS
  #getTBL(lambda_hat_gls, beta_vec, "gls")
    GLS_iid     <- getSE(lambda_hat_gls, beta_vec, omega_iid, V_f, "gls")
    GLS_0lag    <- getSE(lambda_hat_gls, beta_vec, omega_0lag, V_f, "gls")
    GLS_3NW     <- getSE(lambda_hat_gls, beta_vec, omega_3NW, V_f, "gls")
    GLS_24lag   <- getSE(lambda_hat_gls, beta_vec, omega_24lag, V_f, "gls")
  # compile table
    title       <- c("Estimate", "iid", "0-lag", "3 lags, NW", "24 lags")
    tbl_ts      <- t(bind_rows(Estimate=c(lambda_hat_ts, NA), iid=TS_iid, lag0=TS_0lag, NW3=TS_3NW, lag24=TS_24lag))
    tbl_ols     <- t(bind_rows(Estimate=c(lambda_hat_ols, NA), iid=OLS_iid, lag0=OLS_0lag, NW3=OLS_3NW, lag24=OLS_24lag))
    tbl_gls     <- t(bind_rows(Estimate=c(lambda_hat_gls, NA), iid=GLS_iid, lag0=GLS_0lag, NW3=GLS_3NW, lag24=GLS_24lag))
    tbl_15_1    <- bind_cols(title, tbl_ts, tbl_ols, tbl_gls)
    colnames(tbl_15_1) <- c("Item", "TS_EstSE", "TS_t", "OLS_EstSE", "OLS_t", "GLS_EstSE", "GLS_t")

    tbl_15_1 <- kable(x = tbl_15_1, digits = c(0, 4, 2, 4, 2, 4, 2), align = "lrrrrrr", format = "latex")
    xtbl_15_1 <- xtable(tbl_15_1, align = "lrrrrrrr", digits = c(0, 0, 4, 2, 4, 2, 4, 2))
    print(xtbl_15_1, file = "tbl_15_1.tex")
  #######################################################################
  ###
  ### Table 15.2
  ### 0lag, 3NW, 24lag produce only crap ...
  ###
  #######################################################################
  # Table 15.2
  row_title <- c("iid", "GRS F", "0lag", "3lag-NW", "24lag")
  column_title <- c("TS_Chi2", "TS_pval", "TS_Chi2_SHK", "TS_pval_SHK", "CS_Chi2", "CS_pval", "CS_Chi2_SHK", "CS_pval_SHK")

  grs.stats <- GRS.test(ret.mat = r_it, factor.mat = f_t)
  grs_f     <- grs.stats$GRS.stat
  grs_pval  <- grs.stats$GRS.pval
  
  row_iid         <- c(get15_2(ols_model, type = "TS", uDist = "iid", finiteSample = FALSE), get15_2(ols_model, type = "GLS", uDist = "iid", finiteSample = FALSE))
  row_grs         <- c(grs_f, grs_pval, NA, NA, NA, NA, NA, NA)
  row_0lag        <- c(get15_2(ols_model, type = "TS", uDist = "0lag", finiteSample = FALSE), get15_2(ols_model, type = "GLS", uDist = "0lag", finiteSample = FALSE))
  row_3NW         <- c(get15_2(ols_model, type = "TS", uDist = "3NW", finiteSample = FALSE), get15_2(ols_model, type = "GLS", uDist = "3NW", finiteSample = FALSE))
  row_24lag       <- c(get15_2(ols_model, type = "TS", uDist = "24lag", finiteSample = FALSE), get15_2(ols_model, type = "GLS", uDist = "24lag", finiteSample = FALSE))
  
  tbl_15_2            <- rbind(row_iid, row_grs, row_0lag, row_3NW, row_24lag)
  colnames(tbl_15_2)  <- column_title
  rownames(tbl_15_2)  <- row_title
  tbl_15_2_ts         <- tbl_15_2[1:2,1:4]
  tbl_15_2_cs         <- tbl_15_2[1:2,5:8]
  
  xtbl_15_2_ts        <- xtable(tbl_15_2_ts, digits = c(0,2,2,2,2), align = "llrrr")
  xtbl_15_2_cs        <- xtable(tbl_15_2_cs, digits = c(0,2,2,2,2), align = "llrrr")
  print(xtbl_15_2_ts, file = "tbl_15_2_ts.tex")
  print(xtbl_15_2_cs, file = "tbl_15_2_cs.tex")
  
#######################################################################
###
### 4c. bootstrap
###
#######################################################################
  # step 1
    # the model: r_it = alpha_i + beta_i f_t + u_it
    # H0: alpha_i = 0
    # estimate all parameters by OLS (alpha, beta, omega)
    # calculate the original sample of residuals
    alpha           <- E_ri - mean(f_t)*beta
    u               <- r_it - f_t%*%beta
    omega           <- cov(u)
  # step 2
    # compute the test statistic
      tau <- (t(alpha)%*%solve(omega)%*%alpha)
  # step 3
      B <- 999
      tau_b <- rep(0,B)
      for (b in 1:B) {
        # 1 resample with replacement u_b of size T
          u_b <- matrix(data = 0, nrow = nrow(u), ncol = ncol(u))
          for (i in 1:ncol(u)) u_b[,i] <- sample(x = u[,i], size = T, replace = TRUE)
        # 2 generate new sample of r_it using original alpha, beta and u_b from step 2
          r_it_b <- alpha + f_t%*%beta + u_b
        # 3 re-estimate all parameters by OLS (alpha_b, beta_b, omega_b)
          ols_model_b <- lm(r_it_b~f_t); ols_model_b_summary <- summary(ols_model_b)
          beta_b      <- ols_model_b$coefficients[2,]
          alpha_b     <- ols_model_b$coefficients[1,]
          omega_b     <- cov(ols_model_b$residuals)
        # 4 re-compute the test statistic using alpha_b and omega_b
          tau_b[b]    <- (t(alpha_b)%*%solve(omega_b)%*%alpha_b)
      }
  # step 4
    # compute the bootstrapped p-values: p(tau_b) = 1/B sum(|tau_b|>=tau)
      p           <- (1/B) * sum(abs(tau_b)>rep(tau, B))
  # step 5
    # interpret p(t)
      # p(tau_b)>0.05, fail to reject H0
      # p(tau_b)<=0.05, reject H0
    
  
  
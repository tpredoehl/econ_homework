#######################################################################
###
### SETUP
###
#######################################################################
  library(tidyverse)
  require(sandwich)
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
  getNW <- function(L, X, e){
    X <- f_t
    e <- u
    K = ncol(e)
    N = nrow(X)
    XSX <- matrix(data = 0, nrow = K, ncol = K)
    for (i in 0:L){
      wl <- 1-(i/(L+1))
      Xi <- X[1:(N-i),]
      Xj <- X[(1+i):N,]
      ei <- e[1:(N-i),]
      ej <- e[(1+i):N,]
      
      Wi <- as.numeric(t(Xi)%*%Xj) * t(ei)%*%ej
      
      # Wi <- t(Xi)%*%(Xj*ei*ej) + t(Xj)%*%(Xi*ej*ei)
      W_temp <- wl*Wi
      XSX <- XSX+W_temp
    }
    XSX <- N/(N-K)*XSX
    XX <- as.numeric(solve(t(X)%*%X))
    Vb_HAC <- XX*XSX*XX
    Vb_HAC[1,1]
    Vb_NW <- NeweyWest(ols_model, lag = 3, prewhite = FALSE, sandwich = TRUE, adjust = TRUE)
    f_Cols <- grep(pattern = ":f_t", x = colnames(Vb_NW))
    Vb_NW <- Vb_NW[f_Cols,f_Cols]
    Vb_NW[1,1]
    # NO MATCH !
    Vb_HAC[1,1]-Vb_NW[1,1]
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
    if (type == "NW") {
      omega <- omega_3NW
      V_ <- (1/T)*(V_f + BBB%*%omega%*%BBBt)
      SE_ <- sqrt(V_)
    }
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
  ggsave(filename = "replicated_15_1.png", plot = plot_15_1)
  ggsave(filename = "replicated_15_2.png", plot = plot_15_2)
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
    lambda_hat_ols <- solve(t(beta)%*%beta)%*%t(beta) %*% E_ri  # (1.43) lecture notes 3, p. 59
    TS_model <- summary(lm((E_ri)~beta_vec-1))
    lambda_hat_ts <- TS_model$coefficients[2,"Estimate"]
    lambda_hat_gls <- (solve(t(beta_vec)%*%solve(omega)%*%beta_vec)%*%t(beta_vec)%*%solve(omega)%*%E_ri)[2,1]
  #######################################################################
  ###
  ### Estimates of SE of lambda -> Table 15.1
  ###
  #######################################################################
  # iid 
    omega_iid <- mean(diag(omega))*diag(1,nrow = N)
  # 0-lag
    omega_0lag <- diag(diag(omega),nrow = N)
  # 3NW-lag
  # sandwhich package outputs the whole sandwhich 
  # since we need to plug it into the V(b) = (1/T)(Vf + bbb'Obbb), the reults will look like
  # Var(b) = (1/T)(Vf + bbb' [(1/n X'X)^-1] [1/n phi] [(1/n X'X)^-1] bbb)
    omega_3NW <- NeweyWest(ols_model, lag = 3, prewhite = FALSE, sandwich = TRUE, adjust = TRUE)
    f_Cols <- grep(pattern = ":f_t", x = colnames(omega_3NW))
    omega_3NW <- omega_3NW[f_Cols,f_Cols]
  # 24NW-lag
    omega_24lag <- NeweyWest(ols_model, lag = 24, prewhite = FALSE, sandwich = TRUE, adjust = TRUE)
    f_Cols <- grep(pattern = ":f_t", x = colnames(omega_24lag))
    omega_24lag <- omega_24lag[f_Cols,f_Cols]
  # TS
  #getTBL(lambda_hat_ts, beta_vec, "ts")
    TS_iid <- getSE(lambda_hat_ts, beta_vec, omega_iid, V_f, "ts")
    TS_0lag <- getSE(lambda_hat_ts, beta_vec, omega_0lag, V_f, "ts")
    TS_3NW <- getSE(lambda_hat_ts, beta_vec, omega_3NW, V_f, "ts")
    TS_24lag <- getSE(lambda_hat_ts, beta_vec, omega_24lag, V_f, "ts")
  # OLS
  #getTBL(lambda_hat_ols, beta, "ols")
    OLS_iid <- getSE(lambda_hat_ols, beta, omega_iid, V_f, "ols")
    OLS_0lag <- getSE(lambda_hat_ols, beta, omega_0lag, V_f, "ols")
    OLS_3NW <- getSE(lambda_hat_ols, beta, omega_3NW, V_f, "ols")
    OLS_24lag <- getSE(lambda_hat_ols, beta, omega_24lag, V_f, "ols")
  # GLS
  #getTBL(lambda_hat_gls, beta_vec, "gls")
    GLS_iid <- getSE(lambda_hat_gls, beta_vec, omega_iid, V_f, "gls")
    GLS_0lag <- getSE(lambda_hat_gls, beta_vec, omega_0lag, V_f, "gls")
    GLS_3NW <- getSE(lambda_hat_gls, beta_vec, omega_3NW, V_f, "gls")
    GLS_24lag <- getSE(lambda_hat_gls, beta_vec, omega_24lag, V_f, "gls")
  # compile table
    title <- c("Estimate", "iid", "0-lag", "3 lags, NW", "24 lags")
    tbl_ts <- t(bind_rows(Estimate=c(lambda_hat_ts, NA), iid=TS_iid, lag0=TS_0lag, NW3=TS_3NW, lag24=TS_24lag))
    tbl_ols <- t(bind_rows(Estimate=c(lambda_hat_ols, NA), iid=OLS_iid, lag0=OLS_0lag, NW3=OLS_3NW, lag24=OLS_24lag))
    tbl_gls <- t(bind_rows(Estimate=c(lambda_hat_gls, NA), iid=GLS_iid, lag0=GLS_0lag, NW3=GLS_3NW, lag24=GLS_24lag))
    tbl_15_1 <- bind_cols(title, tbl_ts, tbl_ols, tbl_gls)
    colnames(tbl_15_1) <- c("Item", "TS_EstSE", "TS_t", "OLS_EstSE", "OLS_t", "GLS_EstSE", "GLS_t")
    
    opts <- options(knitr.kable.NA = "")
    kable(x = tbl_15_1, digits = c(0, 4, 2, 4, 2, 4, 2), align = "lrrrrrr", "simple")
  #######################################################################
  ###
  ### GRS Tests pricing error -> Table 15.2
  ### alpha takes the form of E_ri - beta lambda
  ###
  #######################################################################
  # critical value of F test
  c(qf(0.90, N, T-N-K), qf(0.95, N, T-N-K), qf(0.99, N, T-N-K))
  # GRS
  grs.stats <- GRS.test(ret.mat = r_it, factor.mat = f_t)
  grs_f <- grs.stats$GRS.stat; grs_f
  grs_pval <- grs.stats$GRS.pval; grs_pval
  # general test setup: (a-0)'V(a)^-1(a-0)'
  theta2_p    <- t(f_t_hat)%*%solve(V_f)%*%f_t_hat        # max Sharpe^2
  theta2_p2   <- f_t_hat^2/V_f
  var_a <- (1/T) * (1 + f_t_hat^2/V_f) * omega
  alpha%*%solve(var_a)%*%alpha
  # test statistics
  aOa <- (t(alpha)%*%solve(omega)%*%alpha)
  aOa_iid <- (t(alpha)%*%solve(omega_iid)%*%alpha)
  # Cochrane (12.3), p. 231     # =Cochrane (12.8), p.234
  T * (1+theta2_p)^-1 * aOa # =alpha%*%solve(var_a)%*%alpha
  T * (1+theta2_p)^-1 * aOa_iid
  # finite sample F-Test  for H0: a1=a2=...=aN = 0, Cochrane (12.4), p. 231
  # e are normal, uncorrelated and homoskedastic
  Xi_2 <- ((T-N-(K-1))/N)*aOa/(1+theta2_p); Xi_2
  Xi_p <- pf(Xi_2, N, T-N-(K-1)); Xi_p                  # pval .41 < .59?
  # finite sample F-test for H0: a1=a2=...=aN = 0
  theta2_p    <- t(f_t_hat)%*%solve(V_f)%*%f_t_hat
  aOa <- (t(alpha)%*%solve(omega)%*%alpha)
  eta2 <- (T-N-(K-1))/N * (1+theta2_p)^-1 * aOa; eta2
  eta2_p <- pf(eta2, N, T-N-(K-1)); eta2_p                # pval .41 < .59?
  
  # OLS Cochrane, p.237
  lambda_hat_ols <- solve(t(beta)%*%beta)%*%t(beta)%*%E_ri
  alpha_ols <- t(E_ri - lambda_hat_ols%*%beta)
  var_lambda_hat_ols <- (1/T)*(solve(t(beta)%*%beta)%*%(beta)%*%omega%*%beta%*%solve(t(beta)%*%beta)+V_f)
  SE_lambda_ols <- sqrt(var_lambda_hat_ols)
  var_a2 <- (1/T)*(V_f*beta%*%t(beta) + omega)
  Mb <- (diag(N)-beta%*%solve(t(beta)%*%beta)%*%t(beta))
  var_a_ols <- Mb%*%omega%*%t(Mb)
  T*t(alpha_ols)%*%solve(var_a2)%*%alpha_ols
  # t(alpha_ols)%*%solve(var_a_ols)%*%alpha_ols             # MbOMb not inversible ...?
  # Shanken correction, p. 240
  SHK <- as.numeric((1+t(lambda_hat_ols%*%solve(V_f)%*%lambda_hat_ols)))
  # t(alpha_ols)%*%solve(var_a_ols*SHK)%*%alpha_ols       # MbOMb not inversible ...?
  # t(alpha_ols)%*%solve(var_a2*SHK)%*%alpha_ols
  
  # GLS Cochrane, p. 238
  lambda_hat_gls <- (solve(t(beta_vec)%*%solve(omega)%*%beta_vec)%*%t(beta_vec)%*%solve(omega)%*%E_ri)[2,1]
  alpha_gls <- E_ri - beta*lambda_hat_gls
  var_lambda_hat_gls <- (1/T)*(solve(t(beta)%*%solve(omega)%*%beta)+V_f)
  SE_lambda_gls <- sqrt(var_lambda_hat_gls)
  var_a_gls <- (1/T)*(omega - beta%*%solve(t(beta)%*%solve(omega)%*%beta)%*%t(beta)) # (12.17), p. 238, = V(e) in (1.46) course notes, Ch 3, p. 64
  # T*t(alpha_gls)%*%solve(var_a_gls)%*%alpha_gls           # not invertible ...?
  T*t(alpha_gls)%*%solve(omega)%*%alpha_gls               # (1.52) Ch3, p. 68
  # Shanken correction, p. 240
  SHK <- as.numeric((1+t(lambda_hat_gls%*%solve(V_f)%*%lambda_hat_gls)))
  T*SHK*t(alpha_gls)%*%solve(omega)%*%alpha_gls           # (12.22), p. 240   
  
  # Table 15.2
  row_title <- c("iid", "GRS F", "0lag", "3lag-NW", "24lag")
  column_title <- c("TS_Chi2", "TS_pval", "TS_Chi2_SHK", "TS_pval_SHK", "CS_Chi2", "CS_pval", "CS_Chi2_SHK", "CS_pval_SHK")
  theta2_p  <- t(f_t_hat)%*%solve(V_f)%*%f_t_hat
  # TS - OK
  alpha           <- E_ri - mean(f_t)*beta
  u               <- r_it - f_t%*%beta
  omega           <- cov(u)
  aOa_iid         <- (t(alpha)%*%solve(omega)%*%alpha)
  eta_iid         <- T * (1+theta2_p)^-1 * aOa_iid; eta_iid
  eta_iid_p       <- pchisq(q = eta_iid, df = N); eta_iid_p
  eta_iid_SHK     <- SHK*eta_iid; eta_iid_SHK
  eta_iid_p_SHK   <- pchisq(q = eta_iid_SHK, df = N-1); eta_iid_p_SHK
  # CS - GLS
  lambda_hat_gls  <- (solve(t(beta_vec)%*%solve(omega)%*%beta_vec)%*%t(beta_vec)%*%solve(omega)%*%E_ri)[2,1]
  alpha_gls       <- E_ri - beta*lambda_hat_gls
  var_a_gls       <- (1/T)*(omega - beta%*%solve(t(beta)%*%solve(omega)%*%beta)%*%t(beta))
  aOa_CS_iid      <- T*(t(alpha_gls)%*%solve(omega)%*%alpha_gls); aOa_CS_iid
  eta_CS_iid_p    <- pchisq(q = aOa_CS_iid, df = N); eta_CS_iid_p
  eta_CS_iid_SHK  <- SHK*aOa_CS_iid; eta_CS_iid_SHK
  eta_CS_iid_p_SHK<- pchisq(q = eta_CS_iid_SHK, df = N-1); eta_CS_iid_p_SHK
  row_iid         <- c(eta_iid, eta_iid_p, eta_iid_SHK, eta_iid_p_SHK, aOa_CS_iid, eta_CS_iid_p, eta_CS_iid_SHK, eta_CS_iid_p_SHK)
  # put the table together...  
  tbl_15_2            <- rbind(row_iid, row_iid, row_iid, row_iid, row_iid)
  colnames(tbl_15_2)  <- column_title
  rownames(tbl_15_2)  <- row_title
  # GRS F row
  grs_row             <- c(grs_f, grs_pval, NA, NA, NA, NA, NA, NA)
  tbl_15_2[2,]        <- grs_row      
  # 0Lag row
  lag0_row            <- c(NA, NA, NA, NA, NA, NA, NA, NA)
  tbl_15_2[3,]        <- lag0_row     
  # 3NW row
  NW3_row             <- c(NA, NA, NA, NA, NA, NA, NA, NA)
  tbl_15_2[4,]        <- NW3_row     
  # lag24 row
  lag24_row           <- c(NA, NA, NA, NA, NA, NA, NA, NA)
  tbl_15_2[5,]        <- lag24_row     
  
  opts <- options(knitr.kable.NA = "")
  kable(x = tbl_15_2, digits = c(4,4,4,4,4,4,4,4), align = "lrrrrrrr", "simple")
  
  
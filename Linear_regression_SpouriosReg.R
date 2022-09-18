############################################################
#Lecture: Time Series
# new time series package:
install.packages('TSA')
library(TSA)

#
install.packages("car")
install.packages('AER')
install.packages("stargazer")
install.packages("tictoc")


# older time series packages already installed:
library(quantmod)  # allows to easily import data directly from downloading financial data from the internet, directly
# from some open sources, including Yahoo Finance, Google Finance, and the Federal
# Reserve Economic Data (FRED) of Federal Reserve Bank of St. Louis.
library(xts)
library(readr)
library(latex2exp) # to write latex formulas in graphs!
library(gridExtra) # multiple plots in one graph
library(summarytools)
library(qwraps2)
#library(normtest)
library(nortest)
library(moments)
library(xtable)
library(sm)
library(astsa)
#library(portes)
library(AER)
library(tseries)
library(car)

library(stargazer)
library(tictoc)


############################################################
# Simulate : White noise, Random Walk
set.seed(102) # fix "seed" in order to produce same random numbers every time next line of code is called!
sigma_ux <- 1 ; 

T <- 300      # T = sample size (length of TS to simulate)

# Simulate white Noise: u_{xt} ~ iidN(0,1), for t = 1, ..., T   
u = rnorm(T+1, mean = 0, sd = sigma_ux);  
# Loop to simulate observations of Random Walk : X_t = X_{t-1} + u_{xt} , 
X.RW      <- NA*matrix(0,T+1,1) # empty vector ( T X 1 )
X.RW[1]   <- 0 # initial value of RW set arbitrarily = 0 (unconditional mean)

# 'for' loop to generate all observations of RW form 2 to T
for (t in 2:(T+1)){ X.RW[t] <- X.RW[t-1] + u[t] }
X.RW <- X.RW[-1,] # delete initial observation which was set to 0    
########################################################
# plot simulated RW (X_t) and its innovations (u_xt)
# dev.off()
plot(X.RW, type='o', pch=19, col="blue", main=TeX('White Noise and Random Walk'), xlab='time', ylab='', cex = 0.5)
abline(0,0)
lines(u, type='o', pch=19, col="red", cex = 0.5)
legend("bottomleft", legend=c(TeX('RW: $X_t = X_{t-1} + u_{x,t}$   '), TeX('WN: $u_{x,t} \\sim iidN(0,1)$   ')), 
       pch=c(19,19), col=c("blue", "red"));
########################################################
########################################################
# Simulate AR(1) and plot acf
set.seed(102) # fix "seed" in order to produce same random numbers every time next line of code is called!
sigma_ux <- 1 ; 

T<- 200

u    <- rnorm(T+1, mean = 0, sd = sigma_ux) #WN
phi  <- 0.5                                 # AR 1 coefficient 


X.AR1.stat    <- NA*matrix(0,T,1) # empty vector ( T X 1 )
X.AR1.nonstat <- NA*matrix(0,T,1) # empty vector ( T X 1 )

X.AR1.nonstat[1] <- 0 # both processes start at t= 0 from the value 0
X.AR1.stat[1]    <- 0  
for (t in 2:T){
  X.AR1.nonstat[t] <- X.AR1.nonstat[t-1]      + u[t];
  X.AR1.stat[t]    <- phi*X.AR1.stat[t-1] + u[t] 
}

y.range <-range(as.matrix(c(X.AR1.stat,X.AR1.nonstat)),na.rm = TRUE, finite = FALSE) 

dev.off()
par(mfrow=c(2,2)) # 4 whiten noises
plot(X.AR1.nonstat, type='o', pch=19, col="blue", main=TeX('AR(1): $X_t = X_{t-1} + u_{x,t},  u_{x,t} \\sim iidN(0,1)$'), xlab='time', ylab='', cex = 0.5)
abline(0,0)
plot(X.AR1.stat, type='o', pch=19, col="red", main=TeX('AR(1): $X_t = 0.5 X_{t-1} + u_{x,t},  u_{x,t} \\sim iidN(0,1)$'), xlab='time', ylab='', cex = 0.5)
#abline(4,0)

# plot Autocorrelation function (acf)
lag.max.acf = 20;
acf(X.AR1.nonstat, main=TeX('ACF, AR(1) NON-stationary'), lag.max = lag.max.acf, xlab = "lag ") #, xlim = c(1,40)
acf(X.AR1.stat, main=TeX('ACF, AR(1) stationary'), lag.max = lag.max.acf, xlab = "lag ") # xlim = c(1,40)
###############################################
# Simulate two independent RW and run spurious regression

# function to simulate K independent AR process with :
# T       = n.f simulated dates (simulated sample size)
# phi_vec = [\phi_1, ..., phi_N]' = autoregressive coefficients for each AR process
# c_vec   = [c_1, ..., c_N] = constant in AR process
independent_AR_simulation<- function(phi_vec, c_vec, sigma_u_vec, T){

  N <- length(phi_vec) ; # number of AR(1) processes to simulate
  
  # create empty matrix T X N that will contain all simulated AR processes
  
  AR_sim_mat <- NA*matrix(0,T+1,N) # empty vector ( (T+1) X 1 )
   u_sim_mat <- NA*matrix(0,T+1,N) # empty vector ( (T+1) X 1 )
  
  # loop for different AR processes
  for (ind1 in seq(1,N)) {
    #########################
    # loop for different dates for one AR process

    c       <-       c_vec[ind1]
    phi     <-     phi_vec[ind1]
    sigma_u <- sigma_u_vec[ind1]

    # simulate vector of innovations zero mean, and sd given as input
#    set.seed(100+ind1) # fix "seed" in order to produce same random numbers every time next line of code is called!
    u    <- rnorm(T+1, mean = 0, sd = sigma_u) #WN
    u_sim_mat[,ind1] <- u

    # generate starting observation for the process: unconditional mean (stationary processes) or 0 (for non-stat processes only!)
    if (abs(phi)<1){
    AR_sim_mat[1,ind1] <- c/(1-phi) ; #start from unconditional value of the process
    } else{
    AR_sim_mat[1,ind1] <- 0 ; #start from unconditional value of the process
    }

    # simulate all other observations t=2 to T+1
    #######################
    for (t in seq(2,T+1)) {
      AR_sim_mat[t,ind1] <- phi*AR_sim_mat[t-1,ind1] + u[t];
      # print(c(t,ind1))
      
      }# end of t loop
    ########################

    } # end of ind1 loop
    #########################

  # Output of function
  # 1 = simulated AR process
  # 2 = simulated innovations
  
  return(list(AR_sim_mat = AR_sim_mat, u_sim_mat = u_sim_mat))
  } # end of function
########################
  
# RUN FUNCTION #############################################
phi_vec     = c(0.95, 0.95) 
c_vec       = c(0,0)
sigma_u_vec = c(1,1) 
T           = 60


output_temp <- independent_AR_simulation(phi_vec, c_vec, sigma_u_vec, T)
str(output_temp)

AR_sim <- output_temp$AR_sim_mat
u_sim  <- output_temp$u_sim_mat
########################################################
# plot simulated innovations of AR(1), that is u_{x,t} and $u_{y,t}
# dev.off()
dev.off()
par(mfrow=c(2,2)) # 4 pnales for plots
plot(u_sim[,1], type='o', pch=19, col="blue", main=TeX('AR innovation: $u_{x,t}$'), xlab='time', ylab='', cex = 0.5)
abline(0,0)
plot(u_sim[,2], type='o', pch=19, col="red", main=TeX('AR innovation: $u_{y,t}$'), xlab='time', ylab='', cex = 0.5)
abline(0,0)

# ACFs
lag.max.acf = 30
acf(u_sim[,1], main=TeX('ACF, AR: $x_t = \\phi_1 x_{t-1} + u_{x,t}$'), lag.max = lag.max.acf, xlab = "lag ") #, xlim = c(1,40)
acf(u_sim[,2], main=TeX('ACF, AR: $y_t = \\phi_2 y_{t-1} + u_{y,t}$'), lag.max = lag.max.acf, xlab = "lag ") # xlim = c(1,40)

# scatterplot + regression line
# dev.off()
par(mfrow=c(1,1)) # 4 panels for plots

#################################################
# scatterplot of u_yt vs. u_xt
#################################################
par(mfrow=c(1,1))
plot(x = u_sim[,1], y = u_sim[,2], 
     col = "blue", lwd =1, 
     xlab=TeX("$u_{x,t}$"), ylab=TeX("$u_{y,t}$"),
     main=TeX('$u_{y,t}$ vs. $u_{x,t}$'),   )

#dev.off()
reg1  <- lm(u_sim[,2] ~ u_sim[,1])
#res <-resid(reg1) ; 
# summary(reg1)

#stargazer(list(reg1,reg2,reg3f),type="text")
stargazer(list(reg1),type="text")
summary(reg1)
abline(reg1, col="red", lwd=2) # display regression line on scatterplot (y~x)
########################################################
# plot simulated AR(1)
# dev.off()
dev.off()
par(mfrow=c(2,2)) # 4 whiten noises
plot(AR_sim[,1], type='o', pch=19, col="blue", main=TeX('AR: $x_t = \\phi_1 x_{t-1} + u_{x,t}$'), xlab='time', ylab='', cex = 0.5)
abline(0,0)
plot(AR_sim[,2], type='o', pch=19, col="red", main=TeX('AR: $y_t = \\phi_2 y_{t-1} + u_{y,t}$'), xlab='time', ylab='', cex = 0.5)
abline(0,0)


# legend("bottomleft", legend=c(TeX('AR: $x_t = 0.5 x_{t-1} + u_{x,t}$'), TeX('AR: $y_t = 0.95 y_{t-1} + u_{y,t}$')), 
#        pch=c(19,19), col=c("blue", "red"));
lag.max.acf = 30
acf(AR_sim[,1], main=TeX('ACF, AR: $x_t = \\phi_1 x_{t-1} + u_{x,t}$'), lag.max = lag.max.acf, xlab = "lag ") #, xlim = c(1,40)
acf(AR_sim[,2], main=TeX('ACF, AR: $y_t = \\phi_2 y_{t-1} + u_{y,t}$'), lag.max = lag.max.acf, xlab = "lag ") # xlim = c(1,40)
#################################################
# scatterplot of y_t vs. x_t
#################################################
par(mfrow=c(1,1))
plot(x = AR_sim[,1], y = AR_sim[,2], 
     col = "blue", lwd =1, 
     xlab=TeX("$x_{t}$"), ylab=TeX("$y_{t}$"),
     main=TeX('$y_{t}$ vs. $x_{t}$'),   )

#dev.off()
#############################################################
# Run Spurious regression ###################################
#############################################################
reg2  <- lm(AR_sim[,2] ~ AR_sim[,1])
#res <-resid(reg1) ; 
# summary(reg1)

#stargazer(list(reg1,reg2,reg3f),type="text")
stargazer(list(reg2),type="text")
summary(reg2)
abline(reg2, col="red", lwd=2) # display regression line on scatterplot (y~x)

stargazer(list(reg1,reg2), type = "text")
# Magics: you can save tabel directly in latex !!!!!!!!!!!!!!!!!!!
# stargazer(list(reg1,reg2),
#           align=TRUE, type = "text", no.space = TRUE, 
#           title = "Table X", out = "texoutput/fit.tex")
###################################
###################################
# OLS with matrix algebra
y <- AR_sim[,2] # 
x <- AR_sim[,1] # 

T <- length(y)               # Sample size T = 60
ones_v <- matrix(1,nrow=T,1) # 'ones' vector (60 X 1) 
X <- cbind(ones_v, x )       # 'X' matrix (60 X 2) 
K <- ncol(X);                # K = 2, number of regressors

XX      <- t(X)%*%X          # X'X
# Note: %*% computes the product among conformable matrices -> "row times column" rule 

XX.inv  <- solve(t(X)%*%X)   # (X'X)^(-1) 
# Note: solve(A) computes the inverse of matrix A 

Xy      <- t(X)%*%y          # X'y

beta.OLS    <- XX.inv %*% Xy
beta.OLS.v2 <- ( solve(t(X)%*%X) ) %*% ((t(X)%*%y))

y.hat <- X %*% beta.OLS # fitted values 
e.hat <- y - y.hat      # residuals

SST <- sum((y -mean(y))^2)
SSR <- sum((y.hat -mean(y))^2)
SSE <- sum(e.hat^2)

MSE  <- SSR/(T-K)
RMSE <- sqrt(MSE)


V.hat.b <- SSE/(T-K) * ( solve(t(X)%*%X) )
V.hat.b
vcov(reg2) 

SE.hat.b1 <-sqrt(V.hat.b[1,1]) 
SE.hat.b2 <-sqrt(V.hat.b[2,2]) 

c(SE.hat.b1,SE.hat.b2)
summary(reg2)$coefficients[,2] # se from lm

# Display values on screen
XX
XX.inv
Xy 
beta.OLS
beta.OLS.v2
SST 
SSE 
SSR 

MSE  
RMSE 

R2     <-  1 - SSR/SST 
adj.R2 <-  1 - (SSR/(T-K))/(SST/(T-1)) 

R2
adj.R2
###############################################################
# We want to save regression coefficient
str(reg2)

reg2$coefficients

reg2$coefficients[1]
reg2$coefficients[2] # <- this is what we want to store!, the regression coefficient $\beta$


summary(reg_sim)

str(summary(reg_sim))
# R^2 :
summary(reg2)$r.squared
summary(reg2)$adj.r.squared

summary(reg2)$coefficients[,1] # coefficients
summary(reg2)$coefficients[,2] # se
summary(reg2)$coefficients[,3] # t-statistics
summary(reg2)$coefficients[,4] # p-values

#############################################################
# Simulate 200 Spurious regression ##########################
#############################################################
# N. of MC simulations
Nsim <- 2000

# empty matrix to store intercept, t-stat(intercept), beta, t-stat (beta), and R^2
MC_mat <- NA*matrix(0,Nsim,5)
colnames(MC_mat) <- c("intercept", "t_intercept", "beta", "t_beta", "R2")

# DGP parameters : as above #############
#phi_vec     = c(0.95, 0.95)
#phi_vec     = c(0.80, 0.80) 
phi_vec     = c(1, 1)
c_vec       = c(0,0)
sigma_u_vec = c(1,1) 
# T           = 120
T           = 100
#########################################
tic()
for(ind_sim in seq(1,Nsim)){
  # ind_sim = 1
  output_temp <- independent_AR_simulation(phi_vec, c_vec, sigma_u_vec, T)
  AR_sim <- output_temp$AR_sim_mat
  
  reg_sim  <- lm(AR_sim[,2] ~ AR_sim[,1])

  MC_mat[ind_sim, 1] <- summary(reg_sim)$coefficients[,1][1]  # intercept : value
  MC_mat[ind_sim, 2] <- summary(reg_sim)$coefficients[,3][1]  # intercept : t-stat
  
  MC_mat[ind_sim, 3] <- summary(reg_sim)$coefficients[,1][2]  # beta : value
  MC_mat[ind_sim, 4] <- summary(reg_sim)$coefficients[,3][2]  # beta : t-stat
  
  MC_mat[ind_sim, 5] <- summary(reg_sim)$r.squared            # R^2

}
toc()

head(MC_mat)

# Histogram of estimated t-statistics for beta vs. N(0,1)
# dev.off()
#postscript(file="outputtex/MC_spurious_reg1.eps",width=9,height=6,horizontal = FALSE, onefile=FALSE)
# par(mfrow=c(2,2))
par(mfrow=c(1,2))
data2plot = MC_mat[, 4]; # t-statistics for beta
hist_OUT <- hist(data2plot, freq = FALSE, breaks = 60, col="lightgreen",  xlab="", main=TeX('$\\hat{\\beta}$:  t-stat, MC distribution'), )
norm_y <-   dnorm(hist_OUT$mids, mean=0, sd=1);
lines(x=hist_OUT$mids, y=norm_y,col="red", lwd=2)
abline(v =  qnorm(1-(0.05/2)) , col="black", lwd=2, lty = 2)
abline(v = -qnorm(1-(0.05/2)) , col="black", lwd=2, lty = 2)



# Histogram of R^2
# par(mfrow=c(2,2))
data2plot = MC_mat[, 5]; # t-statistics for beta
hist_OUT <- hist(data2plot, freq = FALSE, breaks = 60, col="lightpink",  xlab="", main=TeX('$R^2$:  MC distribution'), )
abline(v = 0, col="blue", lwd=3)
#dev.off()


# Empirical size of (asymptotic) t-test:

# normal quantile for 5% two sided test:
alpha <- 0.05;
crit_value <- qnorm(1-(alpha/2)) ; crit_value

rejection_indic_vector <- abs(MC_mat[, 4])>crit_value
empircal_size <- sum(rejection_indic_vector)/length(rejection_indic_vector); empircal_size


################################################################################
# Bootstrap and resamping 
# Resampling

a <- c(1,2,3,4,5,6,7,8,9,100)
b<-sample(a,5,replace=TRUE);b

# install.packages("boot")
# library(boot)

#http://math.furman.edu/~dcs/courses/math47/R/library/tseries/html/tsbootstrap.html

install.packages("tseries")
library(tseries)

# same as sample
tsbootstrap(a, nb = 1, statistic = NULL, 
            m = 1, 
            b = 1,
            type = c("block"))

# block bootstrap : fixed block length b=4
tsbootstrap(a, nb = 5, statistic = NULL, 
            m = 1, 
            b = 4,
            type = c("block"))

# block bootstrap : varying block length, with mean block length b=3
tsbootstrap(a, nb = 5, statistic = NULL, 
            m = 1, 
            b = 3,
            type = c("stationary"))

#########################################################
# Estimation of AR(1) models + exogenous variables x_t

# OPTION 1: estimates AR  by least squares : rudimentary ...
y <- AR_sim[2];

?ar.ols


ar.ols(y) # automatic lag selection using AIC !!!!
ar.ols(y, aic=FALSE, order.max = 1) # force to fit AR(1)

model1

# OPTION 2: OLS by using lm, but we need to create new vectors of y_t = y_t and x_t =y_{t-1} of same length!
y_t   <- y[2:length(y)]
y_tm1 <- y[1:(length(y)-1)]

ar1.model <- lm(y_t ~ y_tm1) ; stargazer(list(ar1.model), type="text")

# OPTION 3 : Estimation of AR(1), .., AR(4) models using 'dyn' package
install.packages("dynlm")
library(dynlm)

y.ts <-ts(y) # y needs to be transformed in ``ts'' object ....!!!!!!

head(y.ts)
head(y)

fit.ar1  <- dynlm(y.ts ~ L(y.ts,1)   )
fit.ar2  <- dynlm(y.ts ~ L(y.ts,1) + L(y.ts,2) )
fit.ar3  <- dynlm(y.ts ~ L(y.ts,1) + L(y.ts,2) + L(y.ts,3)  )
fit.ar4  <- dynlm(y.ts ~ L(y.ts,1) + L(y.ts,2) + L(y.ts,3) + L(y.ts,4)  )
summary(fit.ar1)
stargazer(list(ar1.model, fit.ar1, fit.ar2, fit.ar3, fit.ar4), type="text")

xy.ts <-ts( cbind(u_sim[,1],AR_sim[,2]) )
x.ts  <- xy.ts[,1]
y.ts  <- xy.ts[,2]
head(xy.ts)

head(AR_sim)
head(u_sim)

fit.ar1       <- dynlm(y.ts ~ L(y.ts,1)             )
fit.ar1_xt    <- dynlm(y.ts ~ L(y.ts,1) + L(x.ts,0) )
fit.ar1_xt_v2 <- dynlm(y.ts ~ L(y.ts,1) + x.ts      )

stargazer(list(fit.ar1, fit.ar1_xt, fit.ar1_xt_v2), type="text")
########################################################
# END here for PHD #####################################
########################################################









# NOT FOR PHD ##################################################################
# NOT FOR PHD ##################################################################
# NOT FOR PHD ##################################################################
# NOT FOR PHD ##################################################################
# NOT FOR PHD ##################################################################
# NOT FOR PHD ##################################################################
# NOT FOR PHD ##################################################################
# NOT FOR PHD ##################################################################
# NOT FOR PHD ##################################################################
# NOT FOR PHD ##################################################################
# NOT FOR PHD ##################################################################
# NOT FOR PHD ##################################################################
# NOT FOR PHD ##################################################################
# NOT FOR PHD ##################################################################
# NOT FOR PHD ##################################################################
# NOT FOR PHD ##################################################################
# NOT FOR PHD ##################################################################
# NOT FOR PHD ##################################################################
########################################################
# Simulate white noises with different variance or distribution
set.seed(123) 
# Simulate white Noise: e_t ~ iidN(0,1), for t = 1, ..., T   
e1 <- rnorm(T, mean = 0, sd = 1) ;  
e2 <- sqrt(4)*e1

set.seed(123) 
df.t <- 3 
e3 <- ((sqrt(df.t/(df.t-2)))^(-1))*rt(T, df = df.t) # simulate from student's with 3 df and VARIANCE 1 !!!
e4 <- sqrt(4)*e3 # simulate from student's with 3 df and VARIANCE 4 !!!

# Plot of simulated time series
y.range <-range(as.matrix(c(e1,e2,e3,e4)),na.rm = TRUE, finite = FALSE) 

dev.off()
par(mfrow=c(2,2)) # 4 whiten noises
plot(e1, type='o', pch=19, col="blue", main=TeX('Gaussian WN: $\\epsilon_t \\sim iidN(0,1)$'), ylim = y.range, xlab='time', ylab='', cex = 0.5)
abline(0,0)
plot(e2, type='o', pch=19, col="red", main=TeX('Gaussian  WN: $\\epsilon_t \\sim iidN(0,4)$'), ylim = y.range, xlab='time', ylab='', cex = 0.5)
abline(0,0)
plot(e3, type='o', pch=19, col="blue", main=TeX('Students t WN: $\\epsilon_t \\sim iid t(3), var = 1$'), ylim = y.range, xlab='time', ylab='', cex = 0.5)
abline(0,0)
plot(e4, type='o', pch=19, col="red", main=TeX('Students t WN: $\\epsilon_t \\sim iid t(3), var = 4$'), ylim = y.range, xlab='time', ylab='', cex = 0.5)
abline(0,0)

# plot of the 4 ACF of white noises
lag.max.acf = 40;  lim.y.axes = c(-1,1)
dev.off(); par(mfrow=c(2,2)) # 4 ACF
Acf(e1, main=TeX('ACF, Gaussian WN: $\\epsilon_t \\sim iidN(0,1)$'), lag.max = lag.max.acf, xlab = "lag ", ylim=lim.y.axes, xlim = c(1,40))
Acf(e2, main=TeX('ACF, Gaussian WN: $\\epsilon_t \\sim iidN(0,4)$'), lag.max = lag.max.acf, xlab = "lag ", ylim=lim.y.axes, xlim = c(1,40))
Acf(e3, main=TeX('ACF, Students t WN: $\\epsilon_t \\sim iid t(3), var = 1$'), lag.max = lag.max.acf, xlab = "lag ", ylim=lim.y.axes, xlim = c(1,40))
Acf(e4, main=TeX('ACF, Students t WN: $\\epsilon_t \\sim iid t(3), var = 4$'), lag.max = lag.max.acf, xlab = "lag ", ylim=lim.y.axes, xlim = c(1,40))

############################################################
# Simulate : White noise, AR(1) and Random Walk
# Loop to generaterandow walk: X_t = X_{t-1} + e_t , 
X.RW.nodrift    <- NA*matrix(0,T,1) # empty vector ( T X 1 )
X.RW.drift      <- NA*matrix(0,T,1) # empty vector ( T X 1 )
X.RW.nodrift[1] <- 0 
X.RW.drift[1]   <- 0 

mu <- 0.1
# 'for' loop to generate all observations of RW form 2 to T
for (t in 2:T){ 
  X.RW.nodrift[t] <- 0  + X.RW.nodrift[t-1] + e1[t];
  X.RW.drift[t]   <- mu + X.RW.drift[t-1]   + e1[t]; }

y.range <-range(as.matrix(c(X.RW.nodrift,X.RW.drift)),na.rm = TRUE, finite = FALSE) 

dev.off()
par(mfrow=c(2,2)) # 4 whiten noises
plot(X.RW.nodrift, type='o', pch=19, col="blue", main=TeX('Gaussian RW: $X_t = X_{t-1} + \\epsilon_t,  \\epsilon_t \\sim iidN(0,1)$'), ylim = y.range, xlab='time', ylab='', cex = 0.5)
abline(0,0)
plot(X.RW.drift, type='o', pch=19, col="red", main=TeX('Gaussian RW: $X_t = 0.1 + X_{t-1} + \\epsilon_t,  \\epsilon_t \\sim iidN(0,1)$'), ylim = y.range, xlab='time', ylab='', cex = 0.5)
abline(0,0)

Acf(X.RW.nodrift, main=TeX('ACF, Gaussian RW: $X_t = X_{t-1} + \\epsilon_t,\\epsilon_t \\sim iidN(0,1)$'), lag.max = lag.max.acf, xlab = "lag ", ylim=lim.y.axes, xlim = c(1,40))
Acf(X.RW.drift, main=TeX('ACF, Gaussian RW: $X_t = 0.1 + X_{t-1} + \\epsilon_t,\\epsilon_t \\sim iidN(0,1)$'), lag.max = lag.max.acf, xlab = "lag ", ylim=lim.y.axes, xlim = c(1,40))
########################################################
# Simulate MA(1) and plot acf
T = 150

set.seed(123) 
e    <- rnorm(T, mean = 0, sd = 1) #WN
theta <- 0.6 # MA coefficient 
X.MA1.v1 <- NA*matrix(0,T,1) # empty vector ( T X 1 )

X.MA1.v1[1]  <- 0 
for (t in 2:T){X.MA1.v1[t] <-  theta*e[t-1] + e[t] }

set.seed(124) 
X.MA1.v2 <- arima.sim(model = list(order = c(0, 0, 1), ma = theta), n = T)

y.range <-range(as.matrix(c(X.MA1.v2,X.MA1.v2)),na.rm = TRUE, finite = FALSE) 

dev.off()
par(mfrow=c(2,2)) # 4 whiten noises
plot(X.MA1.v1, type='o', pch=19, col="blue", main=TeX('MA(1): $X_t = 0.6\\epsilon_{t-1} + \\epsilon_t,  \\epsilon_t \\sim iidN(0,1)$'), ylim = y.range, xlab='time', ylab='', cex = 0.5)
abline(0,0)
plot(X.MA1.v2, type='o', pch=19, col="red", main=TeX('MA(1): $X_t = 0.6\\epsilon_{t-1} +  \\epsilon_t,  \\epsilon_t \\sim iidN(0,1)$'), ylim = y.range, xlab='time', ylab='', cex = 0.5)
abline(0,0)

lim.y.axes = c(-0.2,0.6)
Acf(X.MA1.v1, main=TeX('ACF, MA(1): $X_t = 0.6\\epsilon_{t-1} + \\epsilon_t'), lag.max = lag.max.acf, xlab = "lag ", ylim=lim.y.axes, xlim = c(1,40))
Acf(X.MA1.v2, main=TeX('ACF, MA(1): $X_t = 0.6\\epsilon_{t-1} + \\epsilon_t'), lag.max = lag.max.acf, xlab = "lag ", ylim=lim.y.axes, xlim = c(1,40))

0.6/(1+0.6^2)
c(var(X.MA1.v1),var(X.MA1.v2))

acf.v1 <- Acf(X.MA1.v1, main=TeX('ACF, MA(1): $X_t = 0.6\\epsilon_{t-1} + \\epsilon_t'), lag.max = lag.max.acf, xlab = "lag ", ylim=lim.y.axes, xlim = c(1,40))
acf.v2 <- Acf(X.MA1.v2, main=TeX('ACF, MA(1): $X_t = 0.6\\epsilon_{t-1} + \\epsilon_t'), lag.max = lag.max.acf, xlab = "lag ", ylim=lim.y.axes, xlim = c(1,40))
########################################################
########################################################
# Simulate MA(1) and plot acf
T = 4000

set.seed(123) 
ma1.sim<-arima.sim(model=list(ma=c(.7)),n=T, sd=1) # MA(1)
set.seed(123) 
ma2.sim<-arima.sim(model=list(ma=c(.7, 0.5)),n=T, sd=1) # MA(2)
set.seed(123) 
ma3.sim<-arima.sim(model=list(ma=c(.7, 0.5, 0.3)),n=T, sd=1) # specifies standard deviation of Guausian innovations to be equal to 1
set.seed(123) 
ma4.sim<-arima.sim(model=list(ma=c(.7, 0.5, 0.3, -0.5)),n=T, sd=1) # specifies standard deviation of Guausian innovations to be equal to 1


dev.off()
lag.max.acf = 25
par(mfrow=c(2,2)) # 2 white noises
Acf(ma1.sim, main=TeX('ACF, MA(1): $X_t = 0.7\\epsilon_{t-1} + \\epsilon_t$'), lag.max = lag.max.acf, xlab = "lag ", xlim = c(1,lag.max.acf))
Acf(ma2.sim, main=TeX('ACF, MA(2): $X_t = 0.7\\epsilon_{t-1} + 0.5\\epsilon_{t-2} + 0.1X_{t-2} +\\epsilon_t$'), lag.max = lag.max.acf, xlab = "lag ", xlim = c(1,lag.max.acf))
Acf(ma3.sim, main=TeX('ACF, MA(3): $X_t = 0.7\\epsilon_{t-1} + 0.5\\epsilon_{t-2} + 0.3\\epsilon_{t-3} + \\epsilon_t$'), lag.max = lag.max.acf, xlab = "lag ", xlim = c(1,lag.max.acf))
Acf(ma4.sim, main=TeX('ACF, MA(4): $X_t = 0.7\\epsilon_{t-1} + 0.5\\epsilon_{t-2} + 0.3\\epsilon_{t-3} - 0.5\\epsilon_{t-4} +\\epsilon_t$'), lag.max = lag.max.acf, xlab = "lag ", xlim = c(1,lag.max.acf))

########################################################
# Simulate AR(1) and plot acf
T = 500

set.seed(123) 
e    <- rnorm(T, mean = 0, sd = sqrt(0.25)) #WN
phi  <- 0.5 # AR 1 coefficient 
X.AR1.stat    <- NA*matrix(0,T,1) # empty vector ( T X 1 )
X.AR1.nonstat <- NA*matrix(0,T,1) # empty vector ( T X 1 )

X.AR1.nonstat[1] <- 0
X.AR1.stat[1]    <- 4  
for (t in 2:T){
  X.AR1.nonstat[t] <- X.AR1.nonstat[t-1]      + e[t];
  X.AR1.stat[t]    <- 2 + phi*X.AR1.stat[t-1] + e[t] }

y.range <-range(as.matrix(c(X.AR1.stat,X.AR1.nonstat)),na.rm = TRUE, finite = FALSE) 

dev.off()
par(mfrow=c(2,2)) # 4 whiten noises
plot(X.AR1.nonstat, type='o', pch=19, col="blue", main=TeX('AR(1): $X_t = X_{t-1} + \\epsilon_t,  \\epsilon_t \\sim iidN(0,0.25)$'), xlab='time', ylab='', cex = 0.5)
abline(0,0)
plot(X.AR1.stat, type='o', pch=19, col="red", main=TeX('AR(1): $X_t = 2 + 0.5 X_{t-1} + \\epsilon_t,  \\epsilon_t \\sim iidN(0,0.25)$'), xlab='time', ylab='', cex = 0.5)
abline(4,0)

Acf(X.AR1.nonstat, main=TeX('ACF, AR(1) NON-stationary'), lag.max = lag.max.acf, xlab = "lag ", xlim = c(1,40))
Acf(X.AR1.stat, main=TeX('ACF, AR(1) stationary'), lag.max = lag.max.acf, xlab = "lag ", xlim = c(1,40))
########################################################
########################################################
# # Simulate AR(1) and AR(2) and plot acf
# T = 1000
# 
# set.seed(123) 
# ar2.sim<-arima.sim(model=list(ar=c(-.5, +0.2)),n=T, sd=1) # specifies standard deviation of Guausian innovations to be equal to 1
# 
# set.seed(123) 
# e       <- rnorm(T, mean = 0, sd = sqrt(1)) #WN
# ar3.sim <- NA*matrix(0,T,1) # empty vector ( T X 1 )
# 
# ar3.sim[1]<- 0
# ar3.sim[2]<- 0
# ar3.sim[3]<- 0
# 
# for (t in 4:T){
#   ar3.sim[t]    <- 0.2*ar3.sim[t-1] + 1.5*ar3.sim[t-3] + e[t] }
# 
# lag.max.acf = 25
# 
# dev.off()
# par(mfrow=c(2,2)) # 4 whiten noises
# plot(ar2.sim, type='o', pch=19, col="blue", main=TeX('AR(2) $X_t = -0.5X_{t-1} + 0.2X_{t-2} \\epsilon_t,  \\epsilon_t \\sim iidN(0,0.25)$'), xlab='time', ylab='', cex = 0.5)
# plot(ar3.sim, type='o', pch=19, col="red" , main=TeX('AR(1): $X_t = 0.2X_{t-1} + 1.5X_{t-3} +\\epsilon_t,  \\epsilon_t \\sim iidN(0,0.25)$'), xlab='time', ylab='', cex = 0.5)
# 
# Acf(ar2.sim, main=TeX('ACF, AR(2) stationary, $X_t = -0.5X_{t-1} + 0.2X_{t-2} \\epsilon_t$'), lag.max = lag.max.acf, xlab = "lag ", xlim = c(1,lag.max.acf))
# Acf(ar3.sim, main=TeX('ACF, AR(3) NON-stationary, $X_t = 0.2X_{t-1} + 1.5X_{t-3} +\\epsilon_t$'), lag.max = lag.max.acf, xlab = "lag ", xlim = c(1,lag.max.acf))
# 
# dev.off()
########################################################
# Simulate AR(1) and AR(2) and plot acf
T = 5000

set.seed(123) 
ar1.sim<-arima.sim(model=list(ar=c(.7)),n=T, sd=1) # specifies standard deviation of Guaussian innovations to be equal to 1
set.seed(123) 
ar2.sim<-arima.sim(model=list(ar=c(.7, .1)),n=T, sd=1) 

dev.off()
par(mfrow=c(2,2)) # 2 whiten noises
Acf(ar1.sim, main=TeX('ACF, AR(1) stationary, $X_t = 0.7X_{t-1} + \\epsilon_t$'), lag.max = lag.max.acf, xlab = "lag ", xlim = c(1,40))
Acf(ar2.sim, main=TeX('ACF, AR(2) stationary, $X_t = 0.7X_{t-1} + 0.1X_{t-2} +\\epsilon_t$'), lag.max = lag.max.acf, xlab = "lag ", xlim = c(1,40))
Pacf(ar1.sim, main=TeX('PACF, AR(1) stationary, $X_t = 0.7X_{t-1} + \\epsilon_t$'), lag.max = lag.max.acf, xlab = "lag ", xlim = c(1,40))
Pacf(ar2.sim, main=TeX('PACF, AR(2) stationary, $X_t = 0.7X_{t-1} + 0.1X_{t-2} +\\epsilon_t$'), lag.max = lag.max.acf, xlab = "lag ", xlim = c(1,40))
########################################################
# Simulate MA(1) and MA(2) and plot acf
T = 5000

set.seed(123) 
ma1.sim<-arima.sim(model=list(ma=c(.8)),n=T, sd=1) # specifies standard deviation of Gaussian innovations to be equal to 1
set.seed(123) 
ma2.sim<-arima.sim(model=list(ma=c(.8, .1)),n=T, sd=1) 

lag.max.acf = 25
dev.off()
par(mfrow=c(2,2)) # 2 whiten noises
Acf(ma1.sim, main=TeX('ACF, MA(1) stationary, $X_t = 0.8\\epsilon_{t-1} + \\epsilon_t$'), lag.max = lag.max.acf, xlab = "lag ", xlim = c(1,lag.max.acf))
Acf(ma2.sim, main=TeX('ACF, MA(2) stationary, $X_t = 0.8\\epsilon_{t-1} + 0.1\\epsilon_{t-2} +\\epsilon_t$'), lag.max = lag.max.acf, xlab = "lag ", xlim = c(1,lag.max.acf))
Pacf(ma1.sim, main=TeX('PACF, MA(1) stationary, $X_t = 0.8\\epsilon_{t-1} + \\epsilon_t$'), lag.max = lag.max.acf, xlab = "lag ", xlim = c(1,lag.max.acf))
Pacf(ma2.sim, main=TeX('PACF, MA(2) stationary, $X_t = 0.8\\epsilon_{t-1} + 0.1\\epsilon_{t-2} +\\epsilon_t$'), lag.max = lag.max.acf, xlab = "lag ", xlim = c(1,lag.max.acf))
########################################################
# Simulate ARMA(1,1) and ARMA(2,2) and plot acf
T = 5000

set.seed(123) 
arma11 <-arima.sim(model=list(order = c(1, 0, 1), ar=c(0.7), ma=c(.8)),n=T, sd=1) # specifies standard deviation of Guausian innovations to be equal to 1
set.seed(123) 
arma22  <-arima.sim(model=list(order = c(2, 0, 2), ar=c(0.7, 0.2), ma=c(.8, .1)),n=T, sd=1) 

lag.max.acf = 25
dev.off()
par(mfrow=c(2,2)) # 2 whiten noises
Acf(arma11, main=TeX('ACF, ARMA(1,1), $X_t = 0.7X_{t-1} + 0.8\\epsilon_{t-1} + \\epsilon_t$, T = 5000'), lag.max = lag.max.acf, xlab = "lag ", xlim = c(1,lag.max.acf), ylim = c(-0.5,1))
Acf(arma22, main=TeX('ACF, ARMA(2,2), $X_t = 0.7X_{t-1} + 0.2X_{t-2} + 0.8\\epsilon_{t-1} + 0.1\\epsilon_{t-2} +\\epsilon_t$,  T = 5000'), lag.max = lag.max.acf, xlab = "lag ", xlim = c(1,lag.max.acf), ylim = c(-0.5,1))
Pacf(arma11, main=TeX('PACF, ARMA(1,1), $X_t = 0.7X_{t-1} + 0.8\\epsilon_{t-1} + \\epsilon_t$,  T = 5000'), lag.max = lag.max.acf, xlab = "lag ", xlim = c(1,lag.max.acf), ylim = c(-0.5,1))
Pacf(arma22, main=TeX('PACF, ARMA(2,2), $X_t = 0.7X_{t-1} + 0.2X_{t-2} + 0.8\\epsilon_{t-1} + 0.1\\epsilon_{t-2} +\\epsilon_t$,  T = 5000'), lag.max = lag.max.acf, xlab = "lag ", xlim = c(1,lag.max.acf), ylim = c(-0.5,1))


########################################################
# Compute ACF and PACF on simulated data
getSymbols('INDPRO', src="FRED")
last_day_of_quarter  <- endpoints(INDPRO, on = "quarters")
ip.q        <- to.quarterly(INDPRO, OHLC = FALSE)
ss_dates    <- "197412/201812" 
ip.ss       <- INDPRO[ss_dates]
ip.ss.gr    <- diff(log(ip.q))
########################################################
# Study of one real time series: US IP 
# (Industrial Production, monthly data seasonally adjusted)
ip.gr       <- ip.ss.gr["197401/201812"]

dev.off
par(mfrow=c(2,2)) # 2 whiten noises
plot(ip.ss, main="IP Index (Q)")
plot(ip.gr, main="IP growth rates (Q)")
Acf(ip.gr, ylab="IP growth rates (Q): ACF")
Pacf(ip.gr,ylab="IP growth rates (Q): PACF")


acf.ip  <- Acf(ip.gr, lag.max = 10)
pacf.ip <- Pacf(ip.gr, lag.max = 10)
cbind(as.matrix(acf.ip$acf[2:11]), as.matrix(pacf.ip$acf))

pacf.ip$acf
########################################################
# Estimation OF AR(1) by MLE
arima(ip.gr,order=c(1,0,0), include.mean=T)
arima(ip.gr,order=c(1,0,1), include.mean=T)
arima(ip.gr,order=c(7,0,0), include.mean=T)


arima(ip.gr,order=c(1,0,0), include.mean=T)
arima(ip.gr,order=c(2,0,0), include.mean=T)
arima(ip.gr,order=c(6,0,0), include.mean=T)
arima(ip.gr,order=c(7,0,0), include.mean=T)
arima(ip.gr,order=c(8,0,0), include.mean=T)

# use least squares
# arima(ip.gr,order=c(1,0,0),method="CSS", include.mean=T)
# arima(ip.gr,order=c(7,0,0),method="CSS", include.mean=T)
# arima(ip.gr,order=c(8,0,0),method="CSS", include.mean=T)

# estimates AR  by least squares
ar.ols(ip.gr) # automatic lag seleciton using AIC !!!!

?ar.ols
ar.ols(ip.gr, aic = FALSE, order.max = 1)
arima(ip.gr,order=c(1,0,0), include.mean=T)

ar.ols.output <- ar.ols(ip.gr, aic = FALSE, order.max = 1)
ar.ols.output$asy.se.coef
# compare with mle

########################################################
my.data <- as.numeric(ip.gr)

my.max.lag      <- 10
lags.all        <- seq(1,my.max.lag,1)
my.acf          <- acf(my.data, lag.max = my.max.lag, plot = FALSE)
my.pacf         <- acf(my.data, lag.max = my.max.lag, type = "partial", plot = FALSE)
my.BoxPierce    <- BoxPierce(my.data, lags=lags.all)
#crit.value.5.BP <- qchisq(0.95,lags.all) 

my.table <- cbind(my.BoxPierce[,1],
                  as.numeric(my.acf$acf),
                  as.numeric(my.pacf$acf),
                  (1.96/sqrt(length(my.data)))*matrix(1,my.max.lag,1),
                  my.BoxPierce[,2],
                  my.BoxPierce[,4]) 
# change apperance to export more easily !
my.table.df <-as.data.frame(my.table)
names(my.table.df)  <- c("lag","acf","pacf","1.96/(T^0.5)","Box-Pierce stat","BP pval")
rownames(my.table.df) <-c()
options(scipen = 999)
a <- data.matrix(my.table.df)
round(a, digits = 3)

#xtable(b, type = "latex")
########################################################
# Estimation of AR(1), .., AR(4) models using 'dyn' package
install.packages("dyn")
library(dyn)
fit.ar1  <- dyn$lm(ip.gr ~ lag(ip.gr, 1))
fit.ar2  <- dyn$lm(ip.gr ~ lag(ip.gr, 1:2))
fit.ar3  <- dyn$lm(ip.gr ~ lag(ip.gr, 1:3))
fit.ar4  <- dyn$lm(ip.gr ~ lag(ip.gr, 1:4))
summary(fit.ar1)
summary(fit.ar2)
summary(fit.ar3)
summary(fit.ar4)
########################################################
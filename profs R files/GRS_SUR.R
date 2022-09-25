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

rm(list=ls()) # clear all variables from enviroment/workspace
# set working directory (wd)
mydirectory <- "C:/Users/mirco/Dropbox/PhD_course_Econometrics_byMircoRubin/R"
setwd(mydirectory); getwd() # set and show wd 


################################################################################
# Import 10 portfolios sorted on Size from Kenneth French website
# monthly data



# install.packages("frenchdata")
library(frenchdata)
# see : https://nareal.github.io/frenchdata/articles/basic_usage.html

browse_french_site()
list_ff_data <- get_french_data_list()

# 'Portfolios Formed on Size'
# browse_details_page(...)
list_ff_data$files_list$name

ff_3f         <- download_french_data('Fama/French 3 Factors')
ff_size_port  <- download_french_data('Portfolios Formed on Size')
ff_6_port     <- download_french_data('6 Portfolios Formed on Size and Book-to-Market (2 x 3)')
ff_25_port    <- download_french_data('25 Portfolios Formed on Size and Book-to-Market (5 x 5)')

browse_details_page(ff_3f)
browse_details_page(ff_size_port)
browse_details_page(ff_6_port)

##############################################################
ff_size_port$subsets$name

ff_10size_all <- ff_size_port$subsets$data[1] # Value Weight Returns -- Monthly

str(ff_10size_all)

ff_10size_all.df <- as.data.frame(ff_10size_all)

colnames(ff_10size_all.df)
keep_col       <-  c("date", "Lo.10",  "Dec.2",  "Dec.3",  "Dec.4", "Dec.5",  "Dec.6",  "Dec.7",  "Dec.8",  "Dec.9",  "Hi.10")
ff_10size.df   <- ff_10size_all.df[keep_col]
##############################################################
# # 6 FF portfolios 
# ff_6_port$subsets$name
# 
# ff_6_all <- ff_6_port$subsets$data[1] # Value Weight Returns -- Monthly
# # keep all columns
# ff_6_all.df <- as.data.frame(ff_6_all)
##############################################################
# # Fama and French 3 factors
ff_3f_all <- ff_3f$subsets$data[1] # Value Weight Returns -- Monthly

ff_3f_all.df <- as.data.frame(ff_3f_all)
ff_temp <-cbind(ff_10size.df,ff_3f_all.df)
# check visually that dates are aligned, and then delete redundant date column
ff_temp <- ff_temp[-12]

# note: size portfolios are long only: subtract risk free rate
for (i in seq(2,11)){
  ff_temp[,i] <- ff_temp[,i] -ff_temp[,15]  
}


save(ff_temp, file="ff10_sample_1927_2022.Rda")
load("ff10_sample_1927_2022.Rda")


# select subset of dates
ds <- subset(ff_temp, ff_temp$date>=192601 & ff_temp$date<=198212) 
# save(ds, file="ff10_sample_CLM_tab53.Rda")
# load("ff10_sample_CLM_tab53.Rda")



# select ff 6
#ff6_sample <- subset(ff_temp, ff_temp$date>=196501 & ff_temp$date<=199412) # sample in LM (book, 1996)


# test_assets <- ff_6_all.df
# install.packages("GRS.test")
library(GRS.test)

test_asstet_sel <- seq(2,11)  # selects 10 size-sorted portfolios
factor_sel      <- seq(12,12) # selects only market
#factor_sel      <- seq(12,14) # selects 3FF

Y <- as.matrix(ds[,test_asstet_sel]) # T X N
#Y <- as.matrix(ff6_sample) # T X N
X <- as.matrix(ds[,factor_sel])      # T X K

Ybar <- colMeans(Y) ; Xbar <- colMeans(X) ; 

T <- nrow(X)
K <- ncol(X)
N <- ncol(Y)

# create demean dataset
Y_tilde <- NA*matrix(0,T,N)
X_tilde <- NA*matrix(0,T,K)

for (t in seq(1,T)){
  Y_tilde[t,] <- Y[t,] - Ybar
  X_tilde[t,] <- X[t,] - Xbar
}


# OLS estimation of Beta matrix and alpha vector
B      <- (t(Y_tilde)%*%X_tilde) %*% solve(t(X_tilde)%*%X_tilde) ; B # attention : N X K
alpha  <- Ybar - B%*%Xbar; alpha


# just to check :
# reg1 <-lm(Y[,10]~X[,1]); summary(reg1)

# generates matrix of residuals

u_mat <- NA*matrix(0,T,N)

for (t in seq(1,T)){
  for (i in seq(1,N)){
  beta_vec_temp <- B[i,]; # K X 1 (A BIT COUNTERINTUITIVE....)
  u_mat[t,i] <- Y[t,i] - (alpha[i,1] + X[t,]%*%beta_vec_temp)
  }
}

# Sanity check for residuals: they should be 0 mean if everything went ok
colMeans(u_mat) # yes: residuals are zero mean (at machine precision)


# var-cov matrix of (zero mean) residuals
Omega <- (t(u_mat)%*%u_mat)/T
# var-cov matrix of factors

Sigma_f <- (t(X_tilde)%*%X_tilde)/T


theta2_p = t(Xbar)%*%solve(Sigma_f)%*%Xbar

# GRS test statistics
Xi_2 <- ((T-N-K)/N)*(t(alpha)%*%solve(Omega)%*%alpha)/(1+theta2_p)
#Xi_2 <- (T*(T - N - 1))/(N*(T - 2))*(t(alpha)%*%solve(Omega)%*%alpha)/(1+theta2_p)
Xi_2
# p-value of GRS test computed with F distribution
pf(Xi_2, N, T-N-K)

# critical value of F test
c(qf(0.90, N, T-N-K), qf(0.95, N, T-N-K), qf(0.99, N, T-N-K))


GRS_test_package <- GRS.test(Y, X)
GRS_test_package

# compare output of funciton with manually computed GRS test: same resutls!
c(GRS_test_package$GRS.stat,Xi_2)
# very similar to our Xi_2

GRS.test






# rep.row<-function(x,n){
#   matrix(rep(x,each=n),nrow=n)
# }
# rep.col<-function(x,n){
#   matrix(rep(x,each=n), ncol=n, byrow=TRUE)
# }
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
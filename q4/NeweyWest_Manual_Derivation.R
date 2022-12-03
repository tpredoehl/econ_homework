data("Investment")
y <- Investment[,"RealInv"][-1]
X <- cbind(1, RealGNP = Investment[,"RealGNP"], RealINT=Investment[,"RealInt"])[-1,]
fm.inv <- lm(y ~ X-1)

N <- nrow(X)
K <- fm.inv$rank
e <- fm.inv$residuals
Ve <- as.numeric((1/(N-K))*t(e)%*%e)

# Variance of the estimators
XXX <- solve(t(X)%*%X) %*% t(X) 
XXXt <- t(XXX)
Vb <- (XXX*Ve)%*%XXXt
Vb
Vb <- Ve*solve(t(X) %*% X)
Vb
vcov(fm.inv)

# Heteroskedasticity Robust
# https://www.real-statistics.com/multiple-regression/autocorrelation/newey-west-standard-errors/
# W0 <- X'X e^2

# lag truncator L
L = 3

# recursive summation
XSX <- matrix(data = 0, nrow = K, ncol = K)
for (i in 0:L){
  wl <- 1-(i/(L+1))
  if(i==0){
    W_temp <- wl*t(X)%*%(X*e^2)
  }
  if(i>0){
    Wi <- t(X[1:(N-i),])%*%(X[(1+i):N,]*e[1:(N-i)]*e[(1+i):N]) + t(X[(1+i):N,])%*%(X[1:(N-i),]*e[1:(N-i)]*e[(1+i):N])
    W_temp <- wl*Wi
  }
  XSX <- XSX+W_temp
}
XSX <- N/(N-K)*XSX
Vb_HAC <- solve(t(X)%*%X)%*%XSX%*%solve(t(X)%*%X)

Vb_HAC
Vb_NW <- NeweyWest(fm.inv, lag = 3, prewhite = FALSE, sandwich = TRUE, adjust = TRUE)
# MATCH !


SE_iid <- sqrt(diag(vcov(fm.inv)))
SE_HAC <- sqrt(diag(Vb_NW))

tblResult <- cbind(beta, SE_iid, beta/SE_iid, SE_HAC, beta/SE_HAC)
colnames(tblResult) <- c("Estimate", "iid", "t_iid", "3LNW", "t_3LNW")
tblResult

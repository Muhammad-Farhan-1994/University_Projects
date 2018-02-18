load("00732350.Rdata")
# Load in the data

#####Question 1#####

##### Part a)

myglm0 <- glm(y~.,family=gaussian,data=mydat)
# Fit the normal linear model

summary(myglm0)$deviance
# Find the deviance of the model
pchisq(summary(myglm0)$deviance,df.residual(myglm0),lower=FALSE)
# Calculate the p-value of obtaining a value at least as great as the deviance

Y <- as.matrix(mydat[1])
# Our observation Vector
Xmat <- cbind(rep(1,nrow(mydat)),mydat[,2:ncol(mydat)])
Xmat <- as.matrix(Xmat)
# Our design matrix

P <- Xmat%*%(solve(t(Xmat)%*%Xmat))%*%t(Xmat)
# Our hat matrix P

Pdiag <- rep(0,nrow(mydat))
for (i in 1:nrow(mydat)){
  Pdiag[i]=P[i,i]
}
# Extract the diagonal elements

res <- as.matrix(Y-Xmat%*%myglm0$coefficients)
# These are the residuals
RSS <- t(res)%*%res
# This is the residual sum of squares

sigmasqhat <- (RSS)/(nrow(mydat)-ncol(mydat))
# This calculates the unbiased estimator of sigma^2

sr <- (res)/(sqrt(sigmasqhat*(1-Pdiag)))
# This is the formula for the standardised residuals

par(mfrow=c(2,2))
plot(myglm0,which=1)
plot(myglm0,which=2)
hist(sr, n=25, freq=FALSE, xlab="Standardised Residuals", main="Histogram of standardised residuals")
plot(myglm0,which=3)
# Produce the plots

##### Part c) and d)

mylm0 <- lm(log(y+1)~., data=mydat)
# Fit the new linear model

NewModelRes <- as.matrix(log(Y+1)-Xmat%*%mylm0$coefficients)
# Work out the residuals for the new model
NewModelRSS <- t(NewModelRes)%*%NewModelRes
# Work out the new RSS
print(NewModelRSS)

library("MASS")
# For the Box-Cox transformation we need this library
par(mfrow=c(1,1))
boxcox(lm((y+1)~.,data=mydat),plotit=TRUE)
# Produce the plot

par(mfrow=c(1,2))
plot(mylm0, which=1)
plot(mylm0, which=3)
# Produce the plots which will be compared to myglm0

oldR2 <- 1-(RSS/sum((Y-mean(Y))^2))
# R^2 of myglm0 (old)

newR2 <- 1-(NewModelRSS/sum((Y-mean(Y))^2))
# R^2 of mylm0 (new)

##### Part e)  Zheng-Loh Method

Wald <- (coef(summary(mylm0))[,1])/(coef(summary(mylm0))[,2])
# Obtain all of the Wald statistics

OWald <- order(abs(Wald),decreasing=TRUE)[1:ncol(mydat)]
# Now order their absolute values

Xmat <- cbind(rep(1,nrow(mydat)),mydat[,2:ncol(mydat)])
Xmat <- as.matrix(Xmat)
# The original full design matrix X

Y <- mydat[,1]
Y <- as.matrix(Y)
# The response vector Y

V=as.vector(rep(0,ncol(mydat)))
# This will store our V quantities

for (j in 1:ncol(mydat)){
  
  newlm <- lm(log(Y+1)~Xmat[,OWald[1:j]]+0)
  # For each value of j we have to fit the model where the design matrix only contains the
  # the covariates corresponding to the j largest absolute Wald Statistics.
  # Because Owald contains the absolute Wald statistics in descending order
  # we can just use this vector to make sure we
  # extract the right columns of X for each j. e.g j=1 we only want the largest absolute 
  # Wald Statistic so we only want the OWald[1] column of X
  
  CurRes <- as.matrix(log(Y+1)-Xmat[,OWald[1:j]]%*%as.matrix(newlm$coefficients))
  # Work out the residuals for each fitted model
  CurRSS <- t(CurRes)%*%CurRes
  # Work out the RSS for each fitted model, RSS_j
  V[j]= CurRSS+j*((summary(mylm0)$sigma)^2)*log(nrow(mydat))
  # summary(mylm0)$sigma^2 is simply the unbiased estimate for sigma^2 from the full model 
}
end

print(V)
par(mfrow=c(1,1))
plot(V,xlab="j",ylab="V(j)")
# Plot V against j

jhat <- order(V)[1]
# Find jhat, the j corresponding to the smallest V
print(jhat)

finallm <- lm(log(Y+1)~Xmat[,OWald[1:jhat]]+0)
# The final model uses the jhat terms with the largest absolute Wald Statistics 
summary(finallm)$coefficients
# Find the fitted parameters of this final model


##### Question 2#####

myglm1 <- glm(y~., data=mydat, family=poisson)
# Fit the poisson GLM

# Part c)
interval <- 
  (coef(summary(myglm1))[1,1][1]-
     coef(summary(myglm1))[1,2][1]*c(-qnorm(0.975),qnorm(0.975)))
interval <- c(interval[2],interval[1])

print(interval)
#coef(summary(myglm1))[1,1][1] is betaHat_1 
#coef(summary(myglm1))[1,2][1] is the 
#standard error of betaHat_1

# Part d)
Dev0 <- deviance(glm(y~X1+X3,data=mydat, family=poisson))
# The deviance of the model in H_0
Dev1 <- deviance(glm(y~X1+X2+X3+X4,data=mydat, family=poisson))
# The deviance of the model in H_1

print(Dev0)
print(Dev1)

pchisq(Dev0-Dev1, 2,lower=FALSE)
# Work out the p-value of seeing something at least as large as Dev0-Dev1 in a
# chisquared with 5-3 degrees of freedom

# Part e)
newglm <- glm(y~X1+X2+X3+X4,data=mydat, family=poisson)
# Fit the glm detailed in part e

Xstar <- c(1,0.1,0.2,0.3,0.4)
# These are our new predictions

etahat <- t(Xstar)%*%coef(newglm)
print(etahat)
# The predicted linear response

muhat <- exp(etahat)
print (muhat)
# The predicted mean response

Xmat <- cbind(rep(1,nrow(mydat)),mydat[,2:5])
Xmat <- as.matrix(Xmat)
# Our design matrix X

Wvec <- newglm$weights
# The final weights of the fitting process

W <- matrix(0,nrow(mydat),nrow(mydat))
for (i in 1:nrow(mydat))
  W[i,i]=Wvec[i]
end
# Wvec is just the vector w_ii so we need to form the diagonal matrix W from Wvec

J <- as.matrix(t(Xmat)%*%W%*%Xmat)
# The matrix J
cov <- solve(J)
# Find the inverse. This is the variance covariance matrix of beta_hat

c=qnorm(0.955)

se <- sqrt(t(Xstar)%*%cov%*%Xstar)
interval <- c(exp(etahat-c*se),exp(etahat+c*se))
print (interval)
# Find the interval in the mean response scale

# Part f)
Xmat <- cbind(rep(1,nrow(mydat)),mydat[,2:ncol(mydat)])
Xmat <- as.matrix(Xmat)

Wvec <- myglm1$weights

W <- matrix(0,nrow(mydat),nrow(mydat))
for (i in 1:nrow(mydat))
  W[i,i]=Wvec[i]
end
# Once again we form the design matrix X and the weight matrix W (for myglm1)

P <- sqrt(W)%*%Xmat%*%solve(t(Xmat)%*%W%*%Xmat)%*%t(Xmat)%*%sqrt(W)
# This is the hat matrix P in the GLM framework

lev <- rep(0,nrow(mydat))
for (i in 1:nrow(mydat))
  lev[i]=P[i,i]
end
# The diagonal elements of P is the definition of the leverage of observation i
plot(lev,xlab='Index',ylab='Leverage')

Y <- mydat[,1]
# The response vector
eta <- Xmat%*%coef(myglm1)
# the linear predictor
mu <- exp(eta)
# The mean prediction
PR <- (Y-mu)/sqrt(mu)
# This is the definition of Pearsons Residuals
PSR <- PR/(sqrt(1-lev))
# Now standardise with the leverages
plot(PSR,xlab="Index",ylab="Pearson's standardised residuals")

# Part g)
order(lev,decreasing=TRUE)[1:5]
tail(sort(lev),5)
# 5 observations with highest leverages

order(PSR,decreasing=TRUE)[1:5]
tail(sort(PSR),5)
# 5 observations with highest residuals

plot(myglm1,which=5)
# Look at the Cook's distance plot

newdat <- rbind(mydat[1:78,],mydat[80:191,])
# remove the 79'th row from the data
newglm1 <- glm(y~., data=newdat, family=poisson)
# Refit the poisson model with the newdata

print(summary(myglm1)$deviance)
print(summary(newglm1)$deviance)
#Find the deviances of both models

par(mfrow=c(1,2))
plot(myglm1,which=3)
plot(newglm1,which=3)
# Compare the 2 diagnostic plots

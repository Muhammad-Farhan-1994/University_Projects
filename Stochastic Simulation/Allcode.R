"g" <- function(x)
{
  density <- function(x) (1/(x*(1-x)))*exp((-1/2)*((-2+log(x/(1-x)))^2))
  t <- 16
  # Set up a vector to hold f(1/16),f(2/16)...f(15/16)
  dvals <- rep(0,t+1)
  for (i in 2:t)
    dvals[i]=density((i-1)*1/t)
  end
  # Set f(0)=f(1)=0
  dvals[1]=0
  dvals[t+1]=0
  
  j=1
  while (x > j/t){
    # This will tell us where exactly x is (e.g between 3/16 and 4/16)
    j=j+1
  }
  # Now that we know where x is we can use the appropriate interpolant
  if (x<(t-1)/t)
    # If x< 15/16 then we will use one of the genuine fits of g(x)
    y <- ((dvals[j+1]-dvals[j])/(1/t))*x+(dvals[j+1]-(j/t)*(dvals[j+1]-dvals[j])/(1/t))
  else
    # If x is between 15/16 and 1 then we want g(x) to be f(15/16)
    y <- dvals[t]
  
  y
}


"Q1m1" <- function(n)
{
  rand <- rep(0,n)
  density <- function(x) (1/(x*(1-x)))*exp((-1/2)*((-2+log(x/(1-x)))^2))
  count <- 0
  k <- 1
  t <- 16
  # m=sup(f/g) can be calculated using optimize on R
  m <- 1.0158
  # Set up a vector to hold f(1/16),f(2/16)...f(15/16)
  dvals <- rep(0,t+1)
  for (i in 2:t)
    dvals[i]=density((i-1)*1/t)
  end
  # Set f(0)=f(1)=0
  dvals[1]=0
  dvals[t+1]=0
  
  # Set up a vector to hold the coefficient of x in g(x)
  # which is just the slope of the line
  lincoeffs <- rep(0,t)
  for (i in 1:t-1)
    lincoeffs[i]=(dvals[i+1]-dvals[i])/(1/t)
  end
  # The final coefficient equals 0
  lincoeffs[t]=0
  
  # Set up a vector to hold the constant terms in g(x)
  const <- rep(0,t)
  for (i in 1:t-1)
    const[i]=dvals[i+1]-(i/(t))*((dvals[i+1]-dvals[i])/(1/t))
  end
  # The final constant term is simply g(15/16) (=f(15/16))
  const[t]=dvals[t]
  
  # Calculate the area underneath each piece by splitting it up
  # into a triangle and a rectangle
  Area <- rep(0,t)
  for (i in 1:t-1)
    Area[i]=(1/t)*dvals[i]+(1/(2*t))*(dvals[i+1]-dvals[i])
  end
  # The last area is just a rectangle
  Area[t]=(1/t)*const[t]
  TArea <- sum(Area)
  
  # Set up a vector to hold the total area under the first piece
  # and then the second piece etc etc
  SumAreas <- rep(0,t)
  for (i in 1:t)
    SumAreas[i]=sum(Area[1:i])
  end
  
  
  Avec <- rep(0,t)
  Bvec <- rep(0,t)
  
  # The coefficient of x^2 in the cdf is half the coefficient
  # of x in g(x)
  for (i in 1:t-1)
    Avec[i]=(dvals[i+1]-dvals[i])/(2/t)
  end
  Avec[t]=0
  
  # The coefficient of x in the cdf is just the constant 
  # term in g(x)
  for (i in 1:t-1)
    Bvec[i]=(dvals[i+1]-(i/t)*((dvals[i+1]-dvals[i])/(1/t)))
  end
  Bvec[t]=const[t]
  
  # The constant of integration can be worked out by subbing in
  # the end points of each fit and equating the integral
  # to the area under the envelope up to that point
  Cvec <- rep(0,t)
  for (i in 1:t)
    Cvec[i]=sum(Area[1:i])-Avec[i]*((i/t)^2)-Bvec[i]*(i/t)
  end
  
  
  while (k <= n)
  {
    u1 <- TArea*runif(1)
    
    # Find out exactly where u1 falls
    j=1
    while (u1>SumAreas[j]){
      j=j+1
    }
    
    # If u1 <= are under the envelope between 0 and 15/16
    # we use inversion on a cdf which was a quadratic
    if(u1<SumAreas[t-1])
      x <- (-Bvec[j]+sqrt((Bvec[j]^2)-4*Avec[j]*(Cvec[j]-u1)))/(2*Avec[j])
    else
      # if the value of u1 is somewhere between the area under
      # the envelope between 15/16 and 1 then we have to use inversion
      # on a cdf that was linear
      x <- (u1-Cvec[t])/(Bvec[t])
    
    # Now we can sample from f(x)
    u2 <- runif(1)
    if (m*g(x)*u2 <= density(x)){rand[k]=x
                                      k <- k+1}
    
    else count <- count +1
  }
}

"shiftf" <- function(x)
{
  #Shifts the function as described in 1a)
  density <- function(x) (1/(x*(1-x)))*exp((-1/2)*((-2+log(x/(1-x)))^2))
  if(x<=0 || x>=1) y <- 0
  if(x>0 && x<=0.5) y <- density(x+0.5)
  if(x<1 && x>=0.5) y <- density(x-0.5)
  y
}

"rou" <- function(n)
{
  rand <- rep(0,n)
  # max value of the function
  max <- 13.49252
  # max of sqrt(f)
  a <- 3.673217
  #max of x*sqrt(f)
  c <- 1.662404
  i <- 1
  count <- 0
  while (i <= n)
  {
    u <- a*runif(1)
    v <- c*runif(1)
    if (u <= sqrt(shiftf((v/u)))) 
    {
      #Undo the shift
      if((v/u)<=0.5) rand[i] <- (v/u)+0.5
      else rand[i] <- (v/u)-0.5                      
      i <- i+1}
    else count <- count+1
  }
}


"HaMmc" <- function(n)
{
  phivec <- rep(0,n)
  #Set up a vector to store all of the 1s and 0s
  #from the indicator function
  ind <- rep(0,n)
  phi <- function(x) (1/(x*(1-x)))*exp((-1/2)*((-2+log(x/(1-x)))^2))
  max <- 13.49252
  u <- runif(n)
  v <- max*runif(n)
  
  
  for (i in 1:n)
    phivec[i]=phi(u[i])
  end
  
  for(i in 1:n)
    ind[i]= v[i]<=phivec[i]
  end
  
  theta <- (max/n)*(sum(ind))
  Var <- sqrt(2*pi)*(max-sqrt(2*pi))
}


"pf1MC" <- function(n)
{
  phivec <- rep(0,n)
  phi <- function(x) (1/(x*(1-x)))*exp((-1/2)*((-2+log(x/(1-x)))^2))
  unifs <- runif(n)
  #Simply enter the uniform samples into phi
  for (i in 1:n)
    phivec[i]=phi(unifs[i])
  end
  
  theta <- sum(phivec)/n
  UV <- (1/(n-1))*(sum((phivec-mean(phivec))^2))
}

"AntMC" <- function(n)
{
  phi1vec <- rep(0,n)
  phi2vec <- rep(0,n)
  phi <- function(x) (1/(x*(1-x)))*exp((-1/2)*((-2+log(x/(1-x)))^2))
  unifs <- runif(n)
  #This is for theta1 hat
  for (i in 1:n)
    phi1vec[i]=phi(unifs[i])
  end
  #This is for theta2 hat
  for (i in 1:n)
    phi2vec[i]=phi(1-unifs[i])
  end
  UV1 <- (1/(n-1))*(sum((phi1vec-mean(phi1vec))^2))
  UV2 <- (1/(n-1))*(sum((phi2vec-mean(phi2vec))^2))
  cat("Unbiased Variance theta 1= ",(1/n)*UV1, "\n")
  cat("Unbiased Variance theta 2= ",(1/n)*UV2, "\n")
  #Now we sum them to get theta*
  theta <- sum((phi1vec/2)+(phi2vec/2))/n
  cat("theta = ",theta, "\n")
  cat("Error", abs(theta-sqrt(2*pi)), "\n")
  UV=(1/(2))*(UV1)*(1+cor(phi1vec,phi2vec))
  cat("correlation = ",cor(phi1vec,phi2vec), "\n")
  cat("Unbiased Variance theta* = ",UV, "\n")
  
}

"gfunc" <- function(x)
{
  density <- function(x) (1/(x*(1-x)))*exp((-1/2)*((-2+log(x/(1-x)))^2))
  # Set up f(0.6)
  d <- density(0.6)
  # Set up f_max
  max <- 13.49252
  # Set up g as described
  if (x<=0.6) y <- d*x/0.6
  else if (x<=0.9) y <- ((max-d)/0.3)*(x-0.9)+max
  else y <- max
  y
}



"MC4" <- function(n)
{
  X <- rep(0,n)
  density <- function(x) (1/(x*(1-x)))*exp((-1/2)*((-2+log(x/(1-x)))^2))
  # Set up f(0.6)
  d <- density(0.6)
  # Set up f_max
  max <- 13.49252
  
  # The area under the first part of g
  A1 <- 0.3*d
  # The area under the second part of g
  A2 <- (0.3*d)+((max-d)*0.15)
  # The area under the third part of g
  A3 <- max*0.1
  # The Total area under g
  Area <- A1+A2+A3
  
  # This is needed for the CDF of g
  C2 <- (A1+A2)-((max-d)/0.6)*(0.9^2)-(max-3*(max-d))*(0.9)
  
  i  <- 1
  while (i <= n)
  {
    # Start sampling by inversion
    u1 <- (Area)*runif(1)
    # Depending on which part of the envelople we are under,
    # we will use a different form of the CDF for inversion
    if (u1<=A1) x <- sqrt(1.2*u1/d)
    else if(u1<= A2+A1) x <- ((2*max-3*d)+sqrt(((2*max-3*d)^2)-4*((max-d)/0.6)*(C2-(u1))))/((max-d)/0.3)
    else x <- (u1-(Area-max))/(max)
    X[i]=x
    i=i+1
  }
  # Subsitute each sample into phi
  phi <- rep(0,n)
  for (i in 1:n)
    phi[i]=(Area*density(X[i]))/gfunc(X[i])
  end
  
  theta <- sum(phi)/n
  var <- (1/(n-1))*(sum((phi-mean(phi))^2))
  cat("theta = ",theta, "\n")
  cat("Error=",abs(theta-sqrt(2*pi)),"\n")
  cat("variance Phi", var, "\n")
  cat("Variance Theta",var/n,"\n")
}

density <- function(t){
  d <- function(x) (1/sqrt(2*pi))*(1/(x*(1-x)))*exp((-1/2)*((-2+log(x/(1-x)))^2))
  
  if (t<=0 || t>=1) y <- 0
  else y <- d(t)
  y
}

"MHq1" <- function(n){
  X=rep(0.9,n)
  #Our variance for 20% acceptance
  v=0.3
  burnin=100
  gap=30
  d <- function(x) (1/sqrt(2*pi))*(1/(x*(1-x)))*exp((-1/2)*((-2+log(x/(1-x)))^2))
  for(i in 2:n)
  {
    #Go on the "random walk" notice how we select our Y
    #in keeping with the proposal density
    Y <- X[i-1]+rnorm(1,mean=0,sd=sqrt(v))
    alpha <- density(Y)/density(X[i-1])
    X[i]=X[i-1]+(Y-X[i-1])*(runif(1)<alpha)
  }
  X <- X[(burnin+1):n]
  
  Xthin <- rep(0,floor(length(X)/gap))
  for (i in 1:floor(length(X)/gap))
    Xthin[i]=X[gap*i]
  end
  
  accept_rate=length(unique(X))/(n)
  cat("Acceptance = ",accept_rate, "\n")
}


"MHq1" <- function(n){
  v=0.3
burnin=100
gap=30
#Take 30*n+100 so that after the burn in and
#thinning we are left with n
X=rep(0.9,gap*n+burnin)
d <- function(x) (1/sqrt(2*pi))*(1/(x*(1-x)))*exp((-1/2)*((-2+log(x/(1-x)))^2))
for(i in 2:(gap*n+burnin))
{
Y <- X[i-1]+rnorm(1,mean=0,sd=sqrt(v))
alpha <- density(Y)/density(X[i-1])
X[i]=X[i-1]+(Y-X[i-1])*(runif(1)<=alpha)
}
#Implement the burn in
X <- X[(burnin+1):(gap*n+burnin)]
#Now using thinning
Xthin <- X[seq(1,gap*n,gap)]
accept_rate=length(unique(X))/(gap*n+burnin)
}
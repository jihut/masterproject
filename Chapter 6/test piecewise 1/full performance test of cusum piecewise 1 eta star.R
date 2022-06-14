source("faster approx general population simulation.R")
source("test simulate in and out of control excess hazard.R")
source("Funksjoner.R")

vec.sim.excess.piecewise <- function(t.vec,chi.vec,x.matrix,beta.vec){ # in control simulation excess
  
  cum.baseline.excess.hazard.func <- function(t){
    cum.hazard.old <- 0
    if(t>= t.vec[length(t.vec)] ){
      
      cum.hazard.old <- exp(chi.vec[length(chi.vec)])*(t-t.vec[length(t.vec)])+as.numeric(exp(chi.vec[1:(length(chi.vec)-1)])%*%diff(t.vec))
    } else {
      for(i in 1:(length(t.vec)-1)){
        if(t.vec[i] <= t & t<t.vec[i+1]){
          cum.hazard.new <- cum.hazard.old + (t-t.vec[i])*exp(chi.vec[i])
          cum.hazard.old <- cum.hazard.new
          break
        } else {
          cum.hazard.new <- cum.hazard.old + (t.vec[i+1]-t.vec[i])*exp(chi.vec[i])
          cum.hazard.old <- cum.hazard.new
        }
      }
    }
    return(cum.hazard.old)
  }
  
  time.excess.sim <- rep(0,nrow(x.matrix))
  for(i in 1:nrow(x.matrix)){
    cum.excess.hazard.func <- function(t){
      as.numeric(exp(beta.vec%*%x.matrix[i,]))*cum.baseline.excess.hazard.func(t)
    }
    u<-runif(1)
    eq.solve <- function(t){
      1-exp(-cum.excess.hazard.func(t))-u
    }
    time.sim <- uniroot(eq.solve,lower=0,upper=t.vec[length(t.vec)]+10000000)$root
    time.excess.sim[i] <- time.sim
  }
  
  return(time.excess.sim)
  
}

vec.sim.excess.piecewise.proportional <- function(t.vec,chi.vec,x.matrix,beta.vec,rho){ # out of control simulation excess
  
  cum.baseline.excess.hazard.func <- function(t){
    cum.hazard.old <- 0
    if(t>= t.vec[length(t.vec)] ){
      
      cum.hazard.old <- exp(chi.vec[length(chi.vec)])*(t-t.vec[length(t.vec)])+as.numeric(exp(chi.vec[1:(length(chi.vec)-1)])%*%diff(t.vec))
    } else {
      for(i in 1:(length(t.vec)-1)){
        if(t.vec[i] <= t & t<t.vec[i+1]){
          cum.hazard.new <- cum.hazard.old + (t-t.vec[i])*exp(chi.vec[i])
          cum.hazard.old <- cum.hazard.new
          break
        } else {
          cum.hazard.new <- cum.hazard.old + (t.vec[i+1]-t.vec[i])*exp(chi.vec[i])
          cum.hazard.old <- cum.hazard.new
        }
      }
    }
    return(rho*cum.hazard.old)
  }
  
  time.excess.sim <- rep(0,nrow(x.matrix))
  for(i in 1:nrow(x.matrix)){
    cum.excess.hazard.func <- function(t){
      as.numeric(exp(beta.vec%*%x.matrix[i,]))*cum.baseline.excess.hazard.func(t)
    }
    u<-runif(1)
    eq.solve <- function(t){
      1-exp(-cum.excess.hazard.func(t))-u
    }
    time.sim <- uniroot(eq.solve,lower=0,upper=t.vec[length(t.vec)]+10000000)$root
    time.excess.sim[i] <- time.sim
  }
  
  return(time.excess.sim)
  
}

### First: Simulation to decide threshold c. Want to have a 5% probability of false alarm throughout the follow-up. Intensity=100.

set.seed(42)

max.R.0.75 <- rep(NA,1000) # store maximum value of each simulated CUSUM chart
max.R.1.25 <- rep(NA,1000) # store maximum value of each simulated CUSUM chart
max.R.1.50 <- rep(NA,1000) # store maximum value of each simulated CUSUM chart
start.year <- 2010
end.year <- 2020
chi.vec <- -c(7,6.75,6.5,6.25,6,5.75,5.5,5.75) # the chosen baseline parameters 
t.vec <- c(0,1,2,3,4,5,10,15) # the chosen splitting of follow-up interval 
beta.vec <- c(0.05,0.1,0.5)
t.grid <- seq(from=0,to=10,by=0.01)
for(i in 1:1000){
  lambda <- 100 # chosen arrival rate --> gives around 1000 observations per simulated data set. 
  arrival.sim.vec <- rexp(2000,rate=lambda)
  arrival.sim.vec <- cumsum(arrival.sim.vec)
  arrival.sim.vec <- arrival.sim.vec[arrival.sim.vec<=10]
  n.sim <- length(arrival.sim.vec)
  
  ## Simulate age
  age.sim.vec<- rnorm(n.sim,mean=70,sd=10) # simulate age from N(mean=70,sigma=10). 
  
  ## Simulate gender
  
  gender.sim.vec<- sample(0:1,size = n.sim,replace=T) # simulate gender (0=male, 1=female)
  
  ## Simulate treatment group 
  
  treatment.sim.vec <- sample(0:1,size = n.sim, replace=T) 
  
  ## Simulate population survival times and population censoring indicator 
  
  matrix.time.pop.sim <- vec.pop.sim.test(age.sim.vec,gender.sim.vec,start.year,end.year,arrival.sim.vec)
  time.pop.sim <- matrix.time.pop.sim[,"t.p"]
  delta.p<- matrix.time.pop.sim[,"delta.p"]
  
  ## Simulate excess survival times 
  
  x.matrix <- cbind(age.sim.vec,gender.sim.vec,treatment.sim.vec)
  time.excess.sim <- vec.sim.excess.piecewise(t.vec,chi.vec,x.matrix,beta.vec) # in-control the whole period
  
  time.censoring.sim<- cbind(rexp(n.sim,rate=0.001),end.year-start.year-arrival.sim.vec)
  time.censoring.sim <- apply(time.censoring.sim,1,FUN = min)
  
  ## Intermediate data set
  
  
  ## Pick out minimum of the three different survival times as the observed time. 
  ## Also, need to make a death indicator (delta.i) for all causes. 
  time.observed.sim <- pmin(time.pop.sim,time.excess.sim,time.censoring.sim)
  delta.i <- pmax(delta.p,as.numeric(time.excess.sim<time.pop.sim))*as.numeric(time.censoring.sim>time.observed.sim)
  
  cusum.curve.0.75 <- cusum_r.t.piecewise(start.year,age.sim.vec,gender.sim.vec,x.matrix,time.observed.sim,arrival.sim.vec,delta.i,chi.vec,t.vec,beta.vec,0.75,t.grid,pop.data.male,pop.data.female,2020)
  cusum.curve.1.25 <- cusum_r.t.piecewise(start.year,age.sim.vec,gender.sim.vec,x.matrix,time.observed.sim,arrival.sim.vec,delta.i,chi.vec,t.vec,beta.vec,1.25,t.grid,pop.data.male,pop.data.female,2020)
  cusum.curve.1.50 <- cusum_r.t.piecewise(start.year,age.sim.vec,gender.sim.vec,x.matrix,time.observed.sim,arrival.sim.vec,delta.i,chi.vec,t.vec,beta.vec,1.50,t.grid,pop.data.male,pop.data.female,2020)
  
  max.R.0.75[i] <- max(cusum.curve.0.75[,1])
  max.R.1.25[i] <- max(cusum.curve.1.25[,1])
  max.R.1.50[i] <- max(cusum.curve.1.50[,1])
  
}

quantiles.max.R.0.75 <- quantile(max.R.0.75,probs=c(0,0.9,0.95,0.99))
threshold.c.0.75 <- quantiles.max.R.0.75[3]
quantiles.max.R.1.25 <- quantile(max.R.1.25,probs=c(0,0.9,0.95,0.99))
threshold.c.1.25 <- quantiles.max.R.1.25[3]
quantiles.max.R.1.50 <- quantile(max.R.1.50,probs=c(0,0.9,0.95,0.99))
threshold.c.1.50 <- quantiles.max.R.1.50[3]

### Second: Testing in a case where eta^*=5 such that the jump to out-of control happens after 5 years of monitoring. 

max.test.star.0.75 <- rep(NA,1000)
prop.excess.deaths.star.0.75 <- rep(NA,1000)
max.test.star.1.25 <- rep(NA,1000)
prop.excess.deaths.star.1.25 <- rep(NA,1000)
max.test.star.1.50 <- rep(NA,1000)
prop.excess.deaths.star.1.50 <- rep(NA,1000)
set.seed(43)
for(i in 1:1000){
  lambda <- 100
  arrival.sim.vec.test <- rexp(2000,rate=lambda)
  arrival.sim.vec.test <- cumsum(arrival.sim.vec.test)
  arrival.sim.vec.test <- arrival.sim.vec.test[arrival.sim.vec.test<=10]
  n.sim.test <- length(arrival.sim.vec.test)
  
  ## Simulate age
  age.sim.vec.test <- rnorm(n.sim.test,mean=70,sd=10) # simulate age from N(mean=70,sigma=10). 
  
  ## Simulate gender
  
  gender.sim.vec.test <- sample(0:1,size = n.sim.test,replace=T) # simulate gender (0=male, 1=female)
  
  ## Simulate treatment group 
  
  treatment.sim.vec.test <- sample(0:1,size = n.sim.test,replace=T) 
  
  ## Simulate population survival times and population censoring indicator 
  
  matrix.time.pop.sim.test <- vec.pop.sim.test(age.sim.vec.test,gender.sim.vec.test,start.year,end.year,arrival.sim.vec.test)
  time.pop.sim.test <- matrix.time.pop.sim.test[,"t.p"]
  delta.p.test <- matrix.time.pop.sim.test[,"delta.p"]
  
  ## Simulate excess survival times 
  
  x.matrix.test<- cbind(age.sim.vec.test,gender.sim.vec.test,treatment.sim.vec.test)
  
  ### Patient who arrive during the first 5 years (from 2010-2015) --> still the usual excess hazard
  index.before.eta.star <- which(arrival.sim.vec.test<=5)
  index.after.eta.star <- -(index.before.eta.star)
  
  time.excess.sim.test.0.75 <- numeric(n.sim.test)
  time.excess.sim.test.0.75[index.before.eta.star] <- vec.sim.excess.piecewise(t.vec,chi.vec,x.matrix.test[index.before.eta.star,],beta.vec)
  time.excess.sim.test.0.75[index.after.eta.star] <- vec.sim.excess.piecewise.proportional(t.vec,chi.vec,x.matrix.test[index.after.eta.star,],beta.vec,0.75)
  
  time.excess.sim.test.1.25 <- numeric(n.sim.test)
  time.excess.sim.test.1.25[index.before.eta.star] <- vec.sim.excess.piecewise(t.vec,chi.vec,x.matrix.test[index.before.eta.star,],beta.vec)
  time.excess.sim.test.1.25[index.after.eta.star] <- vec.sim.excess.piecewise.proportional(t.vec,chi.vec,x.matrix.test[index.after.eta.star,],beta.vec,1.25)
  
  time.excess.sim.test.1.50 <- numeric(n.sim.test)
  time.excess.sim.test.1.50[index.before.eta.star] <- vec.sim.excess.piecewise(t.vec,chi.vec,x.matrix.test[index.before.eta.star,],beta.vec)
  time.excess.sim.test.1.50[index.after.eta.star] <- vec.sim.excess.piecewise.proportional(t.vec,chi.vec,x.matrix.test[index.after.eta.star,],beta.vec,1.50)
  
  
  time.censoring.sim.test <- cbind(rexp(n.sim.test,rate=0.001),end.year-start.year-arrival.sim.vec.test)
  time.censoring.sim.test <- apply(time.censoring.sim.test,1,FUN = min)
  
  
  ## Pick out minimum of the three different survival times as the observed time. 
  ## Also, need to make a death indicator (delta.i) for all causes. 
  time.observed.sim.test.0.75 <- pmin(time.pop.sim.test,time.excess.sim.test.0.75,time.censoring.sim.test)
  delta.i.test.0.75 <- pmax(delta.p.test,as.numeric(time.excess.sim.test.0.75<time.pop.sim.test))*as.numeric(time.censoring.sim.test>time.observed.sim.test.0.75)
  time.observed.sim.test.1.25 <- pmin(time.pop.sim.test,time.excess.sim.test.1.25,time.censoring.sim.test)
  delta.i.test.1.25 <- pmax(delta.p.test,as.numeric(time.excess.sim.test.1.25<time.pop.sim.test))*as.numeric(time.censoring.sim.test>time.observed.sim.test.1.25)
  time.observed.sim.test.1.50 <- pmin(time.pop.sim.test,time.excess.sim.test.1.50,time.censoring.sim.test)
  delta.i.test.1.50 <- pmax(delta.p.test,as.numeric(time.excess.sim.test.1.50<time.pop.sim.test))*as.numeric(time.censoring.sim.test>time.observed.sim.test.1.50)
  
  prop.excess.deaths.star.0.75[i] <- sum(time.excess.sim.test.0.75<=time.pop.sim.test & delta.i.test.0.75==1)/sum(delta.i.test.0.75==1)
  prop.excess.deaths.star.1.25[i] <- sum(time.excess.sim.test.1.25<=time.pop.sim.test & delta.i.test.1.25==1)/sum(delta.i.test.1.25==1)
  prop.excess.deaths.star.1.50[i] <- sum(time.excess.sim.test.1.50<=time.pop.sim.test & delta.i.test.1.50==1)/sum(delta.i.test.1.50==1)
  # eta^*=5
  t.grid <- seq(from=0,to=10,by=0.01)
  R.out.control.fast.0.75 <- cusum_r.t.piecewise(start.year,age.sim.vec.test,gender.sim.vec.test,x.matrix.test,time.observed.sim.test.0.75,arrival.sim.vec.test,delta.i.test.0.75,chi.vec,t.vec,beta.vec,0.75,t.grid,pop.data.male,pop.data.female,2020)
  max.test.star.0.75[i] <- max(R.out.control.fast.0.75[,1])
  R.out.control.fast.1.25 <- cusum_r.t.piecewise(start.year,age.sim.vec.test,gender.sim.vec.test,x.matrix.test,time.observed.sim.test.1.25,arrival.sim.vec.test,delta.i.test.1.25,chi.vec,t.vec,beta.vec,1.25,t.grid,pop.data.male,pop.data.female,2020)
  max.test.star.1.25[i] <- max(R.out.control.fast.1.25[,1])
  R.out.control.fast.1.50 <- cusum_r.t.piecewise(start.year,age.sim.vec.test,gender.sim.vec.test,x.matrix.test,time.observed.sim.test.1.50,arrival.sim.vec.test,delta.i.test.1.50,chi.vec,t.vec,beta.vec,1.50,t.grid,pop.data.male,pop.data.female,2020)
  max.test.star.1.50[i] <- max(R.out.control.fast.1.50[,1])
}

exceed.threshold.star.0.75 <- max.test.star.0.75>threshold.c.0.75
sum(exceed.threshold.star.0.75)
exceed.threshold.star.1.25 <- max.test.star.1.25>threshold.c.1.25
sum(exceed.threshold.star.1.25)
exceed.threshold.star.1.50 <- max.test.star.1.50>threshold.c.1.50
sum(exceed.threshold.star.1.50)

# Now, increase intensity from 100 to 1000

### First: Simulation to decide threshold c. Want to have a 5% probability of false alarm throughout the follow-up. Intensity=100.

set.seed(42)

max.R.0.75.larger.lambda <- rep(NA,100) # store maximum value of each simulated CUSUM chart
max.R.1.25.larger.lambda <- rep(NA,100) # store maximum value of each simulated CUSUM chart
max.R.1.50.larger.lambda <- rep(NA,100) # store maximum value of each simulated CUSUM chart
start.year <- 2010
end.year <- 2020
chi.vec <- -c(7,6.75,6.5,6.25,6,5.75,5.5,5.75) # the chosen baseline parameters 
t.vec <- c(0,1,2,3,4,5,10,15) # the chosen splitting of follow-up interval 
beta.vec <- c(0.05,0.1,0.5)
t.grid <- seq(from=0,to=10,by=0.01)
for(i in 1:100){
  lambda <- 1000 # chosen arrival rate --> gives around 10000 observations per simulated data set. 
  arrival.sim.vec <- rexp(20000,rate=lambda)
  arrival.sim.vec <- cumsum(arrival.sim.vec)
  arrival.sim.vec <- arrival.sim.vec[arrival.sim.vec<=10]
  n.sim <- length(arrival.sim.vec)
  
  ## Simulate age
  age.sim.vec<- rnorm(n.sim,mean=70,sd=10) # simulate age from N(mean=70,sigma=10). 
  
  ## Simulate gender
  
  gender.sim.vec<- sample(0:1,size = n.sim,replace=T) # simulate gender (0=male, 1=female)
  
  ## Simulate treatment group 
  
  treatment.sim.vec <- sample(0:1,size = n.sim, replace=T) 
  
  ## Simulate population survival times and population censoring indicator 
  
  matrix.time.pop.sim <- vec.pop.sim.test(age.sim.vec,gender.sim.vec,start.year,end.year,arrival.sim.vec)
  time.pop.sim <- matrix.time.pop.sim[,"t.p"]
  delta.p<- matrix.time.pop.sim[,"delta.p"]
  
  ## Simulate excess survival times 
  
  x.matrix <- cbind(age.sim.vec,gender.sim.vec,treatment.sim.vec)
  time.excess.sim <- vec.sim.excess.piecewise(t.vec,chi.vec,x.matrix,beta.vec) # in-control the whole period
  
  time.censoring.sim<- cbind(rexp(n.sim,rate=0.001),end.year-start.year-arrival.sim.vec)
  time.censoring.sim <- apply(time.censoring.sim,1,FUN = min)
  
  ## Intermediate data set
  
  
  ## Pick out minimum of the three different survival times as the observed time. 
  ## Also, need to make a death indicator (delta.i) for all causes. 
  time.observed.sim <- pmin(time.pop.sim,time.excess.sim,time.censoring.sim)
  delta.i <- pmax(delta.p,as.numeric(time.excess.sim<time.pop.sim))*as.numeric(time.censoring.sim>time.observed.sim)
  
  cusum.curve.0.75.larger.lambda <- cusum_r.t.piecewise(start.year,age.sim.vec,gender.sim.vec,x.matrix,time.observed.sim,arrival.sim.vec,delta.i,chi.vec,t.vec,beta.vec,0.75,t.grid,pop.data.male,pop.data.female,2020)
  cusum.curve.1.25.larger.lambda <- cusum_r.t.piecewise(start.year,age.sim.vec,gender.sim.vec,x.matrix,time.observed.sim,arrival.sim.vec,delta.i,chi.vec,t.vec,beta.vec,1.25,t.grid,pop.data.male,pop.data.female,2020)
  cusum.curve.1.50.larger.lambda <- cusum_r.t.piecewise(start.year,age.sim.vec,gender.sim.vec,x.matrix,time.observed.sim,arrival.sim.vec,delta.i,chi.vec,t.vec,beta.vec,1.50,t.grid,pop.data.male,pop.data.female,2020)
  
  max.R.0.75.larger.lambda[i] <- max(cusum.curve.0.75.larger.lambda[,1])
  max.R.1.25.larger.lambda[i] <- max(cusum.curve.1.25.larger.lambda[,1])
  max.R.1.50.larger.lambda[i] <- max(cusum.curve.1.50.larger.lambda[,1])
  
}

quantiles.max.R.0.75.larger.lambda <- quantile(max.R.0.75.larger.lambda,probs=c(0,0.9,0.95,0.99))
threshold.c.0.75.larger.lambda <- quantiles.max.R.0.75.larger.lambda[3]
quantiles.max.R.1.25.larger.lambda <- quantile(max.R.1.25.larger.lambda,probs=c(0,0.9,0.95,0.99))
threshold.c.1.25.larger.lambda <- quantiles.max.R.1.25.larger.lambda[3]
quantiles.max.R.1.50.larger.lambda <- quantile(max.R.1.50.larger.lambda,probs=c(0,0.9,0.95,0.99))
threshold.c.1.50.larger.lambda <- quantiles.max.R.1.50.larger.lambda[3]

# Very low for rho = 1.25, we try another seed

set.seed(50)

max.R.1.25.larger.lambda.new <- rep(NA,100) # store maximum value of each simulated CUSUM chart
start.year <- 2010
end.year <- 2020
chi.vec <- -c(7,6.75,6.5,6.25,6,5.75,5.5,5.75) # the chosen baseline parameters 
t.vec <- c(0,1,2,3,4,5,10,15) # the chosen splitting of follow-up interval 
beta.vec <- c(0.05,0.1,0.5)
t.grid <- seq(from=0,to=10,by=0.01)
for(i in 1:100){
  lambda <- 1000 # chosen arrival rate --> gives around 10000 observations per simulated data set. 
  arrival.sim.vec <- rexp(20000,rate=lambda)
  arrival.sim.vec <- cumsum(arrival.sim.vec)
  arrival.sim.vec <- arrival.sim.vec[arrival.sim.vec<=10]
  n.sim <- length(arrival.sim.vec)
  
  ## Simulate age
  age.sim.vec<- rnorm(n.sim,mean=70,sd=10) # simulate age from N(mean=70,sigma=10). 
  
  ## Simulate gender
  
  gender.sim.vec<- sample(0:1,size = n.sim,replace=T) # simulate gender (0=male, 1=female)
  
  ## Simulate treatment group 
  
  treatment.sim.vec <- sample(0:1,size = n.sim, replace=T) 
  
  ## Simulate population survival times and population censoring indicator 
  
  matrix.time.pop.sim <- vec.pop.sim.test(age.sim.vec,gender.sim.vec,start.year,end.year,arrival.sim.vec)
  time.pop.sim <- matrix.time.pop.sim[,"t.p"]
  delta.p<- matrix.time.pop.sim[,"delta.p"]
  
  ## Simulate excess survival times 
  
  x.matrix <- cbind(age.sim.vec,gender.sim.vec,treatment.sim.vec)
  time.excess.sim <- vec.sim.excess.piecewise(t.vec,chi.vec,x.matrix,beta.vec) # in-control the whole period
  
  time.censoring.sim<- cbind(rexp(n.sim,rate=0.001),end.year-start.year-arrival.sim.vec)
  time.censoring.sim <- apply(time.censoring.sim,1,FUN = min)
  
  ## Intermediate data set
  
  
  ## Pick out minimum of the three different survival times as the observed time. 
  ## Also, need to make a death indicator (delta.i) for all causes. 
  time.observed.sim <- pmin(time.pop.sim,time.excess.sim,time.censoring.sim)
  delta.i <- pmax(delta.p,as.numeric(time.excess.sim<time.pop.sim))*as.numeric(time.censoring.sim>time.observed.sim)
  
  cusum.curve.1.25.larger.lambda <- cusum_r.t.piecewise(start.year,age.sim.vec,gender.sim.vec,x.matrix,time.observed.sim,arrival.sim.vec,delta.i,chi.vec,t.vec,beta.vec,1.25,t.grid,pop.data.male,pop.data.female,2020)
  
  max.R.1.25.larger.lambda.new[i] <- max(cusum.curve.1.25.larger.lambda[,1])
  
}
quantiles.max.R.1.25.larger.lambda.new <- quantile(max.R.1.25.larger.lambda.new,probs=c(0,0.9,0.95,0.99))
threshold.c.1.25.larger.lambda.new <- quantiles.max.R.1.25.larger.lambda.new[3] # a bit more reasonable

### Second: Testing in a case where eta=5 such that the jump to out-of control happens after 5 years of monitoring. 

max.test.star.0.75.larger.lambda <- rep(NA,100)
prop.excess.deaths.star.0.75.larger.lambda <- rep(NA,100)
max.test.star.1.25.larger.lambda <- rep(NA,100)
prop.excess.deaths.star.1.25.larger.lambda <- rep(NA,100)
max.test.star.1.50.larger.lambda <- rep(NA,100)
prop.excess.deaths.star.1.50.larger.lambda <- rep(NA,100)
set.seed(43)
for(i in 1:100){
  lambda <- 1000
  arrival.sim.vec.test <- rexp(20000,rate=lambda)
  arrival.sim.vec.test <- cumsum(arrival.sim.vec.test)
  arrival.sim.vec.test <- arrival.sim.vec.test[arrival.sim.vec.test<=10]
  n.sim.test <- length(arrival.sim.vec.test)
  
  ## Simulate age
  age.sim.vec.test <- rnorm(n.sim.test,mean=70,sd=10) # simulate age from N(mean=70,sigma=10). 
  
  ## Simulate gender
  
  gender.sim.vec.test <- sample(0:1,size = n.sim.test,replace=T) # simulate gender (0=male, 1=female)
  
  ## Simulate treatment group 
  
  treatment.sim.vec.test <- sample(0:1,size = n.sim.test,replace=T) 
  
  ## Simulate population survival times and population censoring indicator 
  
  matrix.time.pop.sim.test <- vec.pop.sim.test(age.sim.vec.test,gender.sim.vec.test,start.year,end.year,arrival.sim.vec.test)
  time.pop.sim.test <- matrix.time.pop.sim.test[,"t.p"]
  delta.p.test <- matrix.time.pop.sim.test[,"delta.p"]
  
  ## Simulate excess survival times 
  
  x.matrix.test<- cbind(age.sim.vec.test,gender.sim.vec.test,treatment.sim.vec.test)
  
  ### First 5 year (from 2010-2015) --> still the usual excess hazard
  index.before.eta.star <- which(arrival.sim.vec.test<=5)
  index.after.eta.star <- -(index.before.eta.star)
  
  time.excess.sim.test.0.75 <- numeric(n.sim.test)
  time.excess.sim.test.0.75[index.before.eta.star] <- vec.sim.excess.piecewise(t.vec,chi.vec,x.matrix.test[index.before.eta.star,],beta.vec)
  time.excess.sim.test.0.75[index.after.eta.star] <- vec.sim.excess.piecewise.proportional(t.vec,chi.vec,x.matrix.test[index.after.eta.star,],beta.vec,0.75)
  
  time.excess.sim.test.1.25 <- numeric(n.sim.test)
  time.excess.sim.test.1.25[index.before.eta.star] <- vec.sim.excess.piecewise(t.vec,chi.vec,x.matrix.test[index.before.eta.star,],beta.vec)
  time.excess.sim.test.1.25[index.after.eta.star] <- vec.sim.excess.piecewise.proportional(t.vec,chi.vec,x.matrix.test[index.after.eta.star,],beta.vec,1.25)
  
  time.excess.sim.test.1.50 <- numeric(n.sim.test)
  time.excess.sim.test.1.50[index.before.eta.star] <- vec.sim.excess.piecewise(t.vec,chi.vec,x.matrix.test[index.before.eta.star,],beta.vec)
  time.excess.sim.test.1.50[index.after.eta.star] <- vec.sim.excess.piecewise.proportional(t.vec,chi.vec,x.matrix.test[index.after.eta.star,],beta.vec,1.50)
  
  
  time.censoring.sim.test <- cbind(rexp(n.sim.test,rate=0.001),end.year-start.year-arrival.sim.vec.test)
  time.censoring.sim.test <- apply(time.censoring.sim.test,1,FUN = min)
  
  
  ## Pick out minimum of the three different survival times as the observed time. 
  ## Also, need to make a death indicator (delta.i) for all causes. 
  time.observed.sim.test.0.75 <- pmin(time.pop.sim.test,time.excess.sim.test.0.75,time.censoring.sim.test)
  delta.i.test.0.75 <- pmax(delta.p.test,as.numeric(time.excess.sim.test.0.75<time.pop.sim.test))*as.numeric(time.censoring.sim.test>time.observed.sim.test.0.75)
  time.observed.sim.test.1.25 <- pmin(time.pop.sim.test,time.excess.sim.test.1.25,time.censoring.sim.test)
  delta.i.test.1.25 <- pmax(delta.p.test,as.numeric(time.excess.sim.test.1.25<time.pop.sim.test))*as.numeric(time.censoring.sim.test>time.observed.sim.test.1.25)
  time.observed.sim.test.1.50 <- pmin(time.pop.sim.test,time.excess.sim.test.1.50,time.censoring.sim.test)
  delta.i.test.1.50 <- pmax(delta.p.test,as.numeric(time.excess.sim.test.1.50<time.pop.sim.test))*as.numeric(time.censoring.sim.test>time.observed.sim.test.1.50)
  
  prop.excess.deaths.star.0.75.larger.lambda[i] <- sum(time.excess.sim.test.0.75<=time.pop.sim.test & delta.i.test.0.75==1)/sum(delta.i.test.0.75==1)
  prop.excess.deaths.star.1.25.larger.lambda[i] <- sum(time.excess.sim.test.1.25<=time.pop.sim.test & delta.i.test.1.25==1)/sum(delta.i.test.1.25==1)
  prop.excess.deaths.star.1.50.larger.lambda[i] <- sum(time.excess.sim.test.1.50<=time.pop.sim.test & delta.i.test.1.50==1)/sum(delta.i.test.1.50==1)
  # eta^*=5
  t.grid <- seq(from=0,to=10,by=0.01)
  R.out.control.fast.0.75 <- cusum_r.t.piecewise(start.year,age.sim.vec.test,gender.sim.vec.test,x.matrix.test,time.observed.sim.test.0.75,arrival.sim.vec.test,delta.i.test.0.75,chi.vec,t.vec,beta.vec,0.75,t.grid,pop.data.male,pop.data.female,2020)
  max.test.star.0.75.larger.lambda[i] <- max(R.out.control.fast.0.75[,1])
  R.out.control.fast.1.25 <- cusum_r.t.piecewise(start.year,age.sim.vec.test,gender.sim.vec.test,x.matrix.test,time.observed.sim.test.1.25,arrival.sim.vec.test,delta.i.test.1.25,chi.vec,t.vec,beta.vec,1.25,t.grid,pop.data.male,pop.data.female,2020)
  max.test.star.1.25.larger.lambda[i] <- max(R.out.control.fast.1.25[,1])
  R.out.control.fast.1.50 <- cusum_r.t.piecewise(start.year,age.sim.vec.test,gender.sim.vec.test,x.matrix.test,time.observed.sim.test.1.50,arrival.sim.vec.test,delta.i.test.1.50,chi.vec,t.vec,beta.vec,1.50,t.grid,pop.data.male,pop.data.female,2020)
  max.test.star.1.50.larger.lambda[i] <- max(R.out.control.fast.1.50[,1])
}

exceed.threshold.star.0.75.larger.lambda <- max.test.star.0.75.larger.lambda>threshold.c.0.75.larger.lambda
sum(exceed.threshold.star.0.75.larger.lambda)
exceed.threshold.star.1.25.larger.lambda <- max.test.star.1.25.larger.lambda>threshold.c.1.25.larger.lambda.new
sum(exceed.threshold.star.1.25.larger.lambda)
exceed.threshold.star.1.50.larger.lambda <- max.test.star.1.50.larger.lambda>threshold.c.1.50.larger.lambda
sum(exceed.threshold.star.1.50.larger.lambda)

# Now, decrease intensity from 100 to 25

### First: Simulation to decide threshold c. Want to have a 5% probability of false alarm throughout the follow-up. Intensity=100.

set.seed(42)

max.R.0.75.smaller.lambda <- rep(NA,1000) # store maximum value of each simulated CUSUM chart
max.R.1.25.smaller.lambda <- rep(NA,1000) # store maximum value of each simulated CUSUM chart
max.R.1.50.smaller.lambda <- rep(NA,1000) # store maximum value of each simulated CUSUM chart
start.year <- 2010
end.year <- 2020
chi.vec <- -c(7,6.75,6.5,6.25,6,5.75,5.5,5.75) # the chosen baseline parameters 
t.vec <- c(0,1,2,3,4,5,10,15) # the chosen splitting of follow-up interval 
beta.vec <- c(0.05,0.1,0.5)
t.grid <- seq(from=0,to=10,by=0.01)
for(i in 1:1000){
  lambda <- 1000 # chosen arrival rate --> gives around 10000 observations per simulated data set. 
  arrival.sim.vec <- rexp(20000,rate=lambda)
  arrival.sim.vec <- cumsum(arrival.sim.vec)
  arrival.sim.vec <- arrival.sim.vec[arrival.sim.vec<=10]
  n.sim <- length(arrival.sim.vec)
  
  ## Simulate age
  age.sim.vec<- rnorm(n.sim,mean=70,sd=10) # simulate age from N(mean=70,sigma=10). 
  
  ## Simulate gender
  
  gender.sim.vec<- sample(0:1,size = n.sim,replace=T) # simulate gender (0=male, 1=female)
  
  ## Simulate treatment group 
  
  treatment.sim.vec <- sample(0:1,size = n.sim, replace=T) 
  
  ## Simulate population survival times and population censoring indicator 
  
  matrix.time.pop.sim <- vec.pop.sim.test(age.sim.vec,gender.sim.vec,start.year,end.year,arrival.sim.vec)
  time.pop.sim <- matrix.time.pop.sim[,"t.p"]
  delta.p<- matrix.time.pop.sim[,"delta.p"]
  
  ## Simulate excess survival times 
  
  x.matrix <- cbind(age.sim.vec,gender.sim.vec,treatment.sim.vec)
  time.excess.sim <- vec.sim.excess.piecewise(t.vec,chi.vec,x.matrix,beta.vec) # in-control the whole period
  
  time.censoring.sim<- cbind(rexp(n.sim,rate=0.001),end.year-start.year-arrival.sim.vec)
  time.censoring.sim <- apply(time.censoring.sim,1,FUN = min)
  
  ## Intermediate data set
  
  
  ## Pick out minimum of the three different survival times as the observed time. 
  ## Also, need to make a death indicator (delta.i) for all causes. 
  time.observed.sim <- pmin(time.pop.sim,time.excess.sim,time.censoring.sim)
  delta.i <- pmax(delta.p,as.numeric(time.excess.sim<time.pop.sim))*as.numeric(time.censoring.sim>time.observed.sim)
  
  cusum.curve.0.75.smaller.lambda <- cusum_r.t.piecewise(start.year,age.sim.vec,gender.sim.vec,x.matrix,time.observed.sim,arrival.sim.vec,delta.i,chi.vec,t.vec,beta.vec,0.75,t.grid,pop.data.male,pop.data.female,2020)
  cusum.curve.1.25.smaller.lambda <- cusum_r.t.piecewise(start.year,age.sim.vec,gender.sim.vec,x.matrix,time.observed.sim,arrival.sim.vec,delta.i,chi.vec,t.vec,beta.vec,1.25,t.grid,pop.data.male,pop.data.female,2020)
  cusum.curve.1.50.smaller.lambda <- cusum_r.t.piecewise(start.year,age.sim.vec,gender.sim.vec,x.matrix,time.observed.sim,arrival.sim.vec,delta.i,chi.vec,t.vec,beta.vec,1.50,t.grid,pop.data.male,pop.data.female,2020)
  
  max.R.0.75.smaller.lambda[i] <- max(cusum.curve.0.75.smaller.lambda[,1])
  max.R.1.25.smaller.lambda[i] <- max(cusum.curve.1.25.smaller.lambda[,1])
  max.R.1.50.smaller.lambda[i] <- max(cusum.curve.1.50.smaller.lambda[,1])
  
}

quantiles.max.R.0.75.smaller.lambda <- quantile(max.R.0.75.smaller.lambda,probs=c(0,0.9,0.95,0.99))
threshold.c.0.75.smaller.lambda <- quantiles.max.R.0.75.smaller.lambda[3]
quantiles.max.R.1.25.smaller.lambda <- quantile(max.R.1.25.smaller.lambda,probs=c(0,0.9,0.95,0.99))
threshold.c.1.25.smaller.lambda <- quantiles.max.R.1.25.smaller.lambda[3]
quantiles.max.R.1.50.smaller.lambda <- quantile(max.R.1.50.smaller.lambda,probs=c(0,0.9,0.95,0.99))
threshold.c.1.50.smaller.lambda <- quantiles.max.R.1.50.smaller.lambda[3]

### Second: Testing in a case where eta=5 such that the jump to out-of control happens after 5 years of monitoring. 

max.test.star.0.75.smaller.lambda <- rep(NA,1000)
prop.excess.deaths.star.0.75.smaller.lambda <- rep(NA,1000)
max.test.star.1.25.smaller.lambda <- rep(NA,1000)
prop.excess.deaths.star.1.25.smaller.lambda <- rep(NA,1000)
max.test.star.1.50.smaller.lambda <- rep(NA,1000)
prop.excess.deaths.star.1.50.smaller.lambda <- rep(NA,1000)
set.seed(43)
for(i in 1:1000){
  lambda <- 25
  arrival.sim.vec.test <- rexp(20000,rate=lambda)
  arrival.sim.vec.test <- cumsum(arrival.sim.vec.test)
  arrival.sim.vec.test <- arrival.sim.vec.test[arrival.sim.vec.test<=10]
  n.sim.test <- length(arrival.sim.vec.test)
  
  ## Simulate age
  age.sim.vec.test <- rnorm(n.sim.test,mean=70,sd=10) # simulate age from N(mean=70,sigma=10). 
  
  ## Simulate gender
  
  gender.sim.vec.test <- sample(0:1,size = n.sim.test,replace=T) # simulate gender (0=male, 1=female)
  
  ## Simulate treatment group 
  
  treatment.sim.vec.test <- sample(0:1,size = n.sim.test,replace=T) 
  
  ## Simulate population survival times and population censoring indicator 
  
  matrix.time.pop.sim.test <- vec.pop.sim.test(age.sim.vec.test,gender.sim.vec.test,start.year,end.year,arrival.sim.vec.test)
  time.pop.sim.test <- matrix.time.pop.sim.test[,"t.p"]
  delta.p.test <- matrix.time.pop.sim.test[,"delta.p"]
  
  ## Simulate excess survival times 
  
  x.matrix.test<- cbind(age.sim.vec.test,gender.sim.vec.test,treatment.sim.vec.test)
  
  ### First 5 year (from 2010-2015) --> still the usual excess hazard
  index.before.eta.star <- which(arrival.sim.vec.test<=5)
  index.after.eta.star <- -(index.before.eta.star)
  
  time.excess.sim.test.0.75 <- numeric(n.sim.test)
  time.excess.sim.test.0.75[index.before.eta.star] <- vec.sim.excess.piecewise(t.vec,chi.vec,x.matrix.test[index.before.eta.star,],beta.vec)
  time.excess.sim.test.0.75[index.after.eta.star] <- vec.sim.excess.piecewise.proportional(t.vec,chi.vec,x.matrix.test[index.after.eta.star,],beta.vec,0.75)
  
  time.excess.sim.test.1.25 <- numeric(n.sim.test)
  time.excess.sim.test.1.25[index.before.eta.star] <- vec.sim.excess.piecewise(t.vec,chi.vec,x.matrix.test[index.before.eta.star,],beta.vec)
  time.excess.sim.test.1.25[index.after.eta.star] <- vec.sim.excess.piecewise.proportional(t.vec,chi.vec,x.matrix.test[index.after.eta.star,],beta.vec,1.25)
  
  time.excess.sim.test.1.50 <- numeric(n.sim.test)
  time.excess.sim.test.1.50[index.before.eta.star] <- vec.sim.excess.piecewise(t.vec,chi.vec,x.matrix.test[index.before.eta.star,],beta.vec)
  time.excess.sim.test.1.50[index.after.eta.star] <- vec.sim.excess.piecewise.proportional(t.vec,chi.vec,x.matrix.test[index.after.eta.star,],beta.vec,1.50)
  
  
  time.censoring.sim.test <- cbind(rexp(n.sim.test,rate=0.001),end.year-start.year-arrival.sim.vec.test)
  time.censoring.sim.test <- apply(time.censoring.sim.test,1,FUN = min)
  
  
  ## Pick out minimum of the three different survival times as the observed time. 
  ## Also, need to make a death indicator (delta.i) for all causes. 
  time.observed.sim.test.0.75 <- pmin(time.pop.sim.test,time.excess.sim.test.0.75,time.censoring.sim.test)
  delta.i.test.0.75 <- pmax(delta.p.test,as.numeric(time.excess.sim.test.0.75<time.pop.sim.test))*as.numeric(time.censoring.sim.test>time.observed.sim.test.0.75)
  time.observed.sim.test.1.25 <- pmin(time.pop.sim.test,time.excess.sim.test.1.25,time.censoring.sim.test)
  delta.i.test.1.25 <- pmax(delta.p.test,as.numeric(time.excess.sim.test.1.25<time.pop.sim.test))*as.numeric(time.censoring.sim.test>time.observed.sim.test.1.25)
  time.observed.sim.test.1.50 <- pmin(time.pop.sim.test,time.excess.sim.test.1.50,time.censoring.sim.test)
  delta.i.test.1.50 <- pmax(delta.p.test,as.numeric(time.excess.sim.test.1.50<time.pop.sim.test))*as.numeric(time.censoring.sim.test>time.observed.sim.test.1.50)
  
  prop.excess.deaths.star.0.75.smaller.lambda[i] <- sum(time.excess.sim.test.0.75<=time.pop.sim.test & delta.i.test.0.75==1)/sum(delta.i.test.0.75==1)
  prop.excess.deaths.star.1.25.smaller.lambda[i] <- sum(time.excess.sim.test.1.25<=time.pop.sim.test & delta.i.test.1.25==1)/sum(delta.i.test.1.25==1)
  prop.excess.deaths.star.1.50.smaller.lambda[i] <- sum(time.excess.sim.test.1.50<=time.pop.sim.test & delta.i.test.1.50==1)/sum(delta.i.test.1.50==1)
  # eta^*=5
  t.grid <- seq(from=0,to=10,by=0.01)
  R.out.control.fast.0.75 <- cusum_r.t.piecewise(start.year,age.sim.vec.test,gender.sim.vec.test,x.matrix.test,time.observed.sim.test.0.75,arrival.sim.vec.test,delta.i.test.0.75,chi.vec,t.vec,beta.vec,0.75,t.grid,pop.data.male,pop.data.female,2020)
  max.test.star.0.75.smaller.lambda[i] <- max(R.out.control.fast.0.75[,1])
  R.out.control.fast.1.25 <- cusum_r.t.piecewise(start.year,age.sim.vec.test,gender.sim.vec.test,x.matrix.test,time.observed.sim.test.1.25,arrival.sim.vec.test,delta.i.test.1.25,chi.vec,t.vec,beta.vec,1.25,t.grid,pop.data.male,pop.data.female,2020)
  max.test.star.1.25.smaller.lambda[i] <- max(R.out.control.fast.1.25[,1])
  R.out.control.fast.1.50 <- cusum_r.t.piecewise(start.year,age.sim.vec.test,gender.sim.vec.test,x.matrix.test,time.observed.sim.test.1.50,arrival.sim.vec.test,delta.i.test.1.50,chi.vec,t.vec,beta.vec,1.50,t.grid,pop.data.male,pop.data.female,2020)
  max.test.star.1.50.smaller.lambda[i] <- max(R.out.control.fast.1.50[,1])
}

exceed.threshold.star.0.75.smaller.lambda <- max.test.star.0.75.smaller.lambda>threshold.c.0.75.smaller.lambda
sum(exceed.threshold.star.0.75.smaller.lambda)
exceed.threshold.star.1.25.smaller.lambda <- max.test.star.1.25.smaller.lambda>threshold.c.1.25.smaller.lambda
sum(exceed.threshold.star.1.25.smaller.lambda)
exceed.threshold.star.1.50.smaller.lambda <- max.test.star.1.50.smaller.lambda>threshold.c.1.50.smaller.lambda
sum(exceed.threshold.star.1.50.smaller.lambda)


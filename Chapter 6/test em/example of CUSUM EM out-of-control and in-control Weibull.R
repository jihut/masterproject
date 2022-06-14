# Example of an application of CUSUM based on outputs from EM-model on simulated data 

source("Funksjoner.r")
source("final cusum em loess.R")
source("faster approx general population simulation.R")
library(relsurv)
library(dplyr)

lambda <- 0.75 # shape parameter
phi <- 0.005 # scale parameter
rho <- 1.25
beta.vec <- c(0.05,0.1,0.5) 
start.year.training <- 2000
end.year.training <- 2010

# Simulate from out-of-control (proportional alternative)
weibull.excess.one.sim.proportional <- function(lambda,phi,rho,beta.vec,x.vec){ # a function to simulate one excess time from the out-of-control based on one observation
  u <- runif(1)
  inner.prod <- as.numeric(beta.vec%*%x.vec)
  phi.i <- rho*phi*exp(inner.prod)
  t <- (-log(1-u)/phi.i)^(1/lambda) # inverse transform 
  return(t)
}

weibull.excess.sim.time.proportional <- function(lambda,phi,rho,beta.vec,x.matrix){ # lambda=shape, phi=scale of baseline, beta.vec=parameter vector, x.matrix=covariate matrix where each row corresponds to a single observation.
  time.vec <- rep(0,nrow(x.matrix))
  for(i in 1:nrow(x.matrix)){
    time.vec[i] <- weibull.excess.one.sim.proportional(lambda,phi,rho,beta.vec,x.matrix[i,])
  }
  return(time.vec)
}

# Training - used to fit an EM-based model

set.seed(420)

intensity <- 100 # arrival intensity
arrival.sim.vec.training <- rexp(2000,rate=intensity)
arrival.sim.vec.training <- cumsum(arrival.sim.vec.training)
arrival.sim.vec.training <- arrival.sim.vec.training[arrival.sim.vec.training<=10] # only consider arrivals in the period between 0 and 10 years. 
n.sim.training <- length(arrival.sim.vec.training)

## Simulate age
age.sim.vec.training <- rnorm(n.sim.training,mean=70,sd=10) # simulate age from N(mean=70,sigma=10). 

## Simulate gender

gender.sim.vec.training <- sample(0:1,size = n.sim.training,replace=T) # simulate gender (0=male, 1=female)

## Simulate treatment group 

treatment.sim.vec.training <- sample(0:1,size = n.sim.training,replace=T) # simulate treatment 

## Simulate population survival times and population censoring indicator 

matrix.time.pop.sim.training <- vec.pop.sim.test(age.sim.vec.training,gender.sim.vec.training,start.year.training,end.year.training,arrival.sim.vec.training)
time.pop.sim.training <- matrix.time.pop.sim.training[,"t.p"] # population times 
delta.p.training<- matrix.time.pop.sim.training[,"delta.p"] # population censoring indicator 

## Simulate excess survival times 

x.matrix.training <- cbind(age.sim.vec.training,gender.sim.vec.training,treatment.sim.vec.training)
time.excess.sim.training <- weibull.excess.sim.time(lambda,phi,beta.vec,x.matrix.training)

time.censoring.sim.training<- cbind(rexp(n.sim.training,rate=0.001),end.year.training-start.year.training-arrival.sim.vec.training) # simulate censoring times
time.censoring.sim.training <- apply(time.censoring.sim.training,1,FUN = min)

## Intermediate data set

inter.sim.dataset.training <- cbind(x.matrix.training,matrix.time.pop.sim.training,time.excess.sim.training,time.censoring.sim.training)

## Pick out minimum of the three different survival times as the observed time. 
## Also, need to make a death indicator (delta.i) for all causes. 
time.observed.sim.training <- pmin(time.pop.sim.training,time.excess.sim.training,time.censoring.sim.training)
delta.i.training <- pmax(delta.p.training,as.numeric(time.excess.sim.training<time.pop.sim.training))*as.numeric(time.censoring.sim.training>time.observed.sim.training)

# Final training data set

final.data.set.training <- cbind(x.matrix.training,start.year.training+arrival.sim.vec.training,time.observed.sim.training,delta.i.training)
colnames(final.data.set.training) <- c("age","gender","treatment","year","time.obs","delta.i")
final.data.set.training <- as.data.frame(final.data.set.training)

nortab <- transrate.hmd(male="NOR.mltper_1x1.txt",female="NOR.fltper_1x1.txt")
final.data.set.training$gender <- final.data.set.training$gender+1 # gender needs to take the value 1 or 2 in rsadd

test.fit.em.bwin.1 <- rsadd(Surv(time.obs*365.241,delta.i)~age+as.factor(gender)+as.factor(treatment),data=final.data.set.training, # fit an EM-based model for the baseline period
                            ratetable = nortab,method="EM",bwin=1,rmap = list(age=age*365.241,sex=gender,year=(year-1970)*365.241))


estimated.beta <- as.numeric(test.fit.em.bwin.1$coefficients) # estimated beta

time.vec.events <- test.fit.em.bwin.1$times # times to event from the baseline period (in days)
baseline.haz.vec.event.times <- test.fit.em.bwin.1$lambda0 # corresponding baseline excess hazard estimates (in days)
cumulative.baseline.haz.vec.event.times <- test.fit.em.bwin.1$Lambda0 # corresponding cumulative baseline excess hazard estimates (in days)

# Apply on test data to see how it works when eta=0, i.e. out of control right away 

start.year.test <- 2010
end.year.test <- 2020
intensity <- 100
arrival.sim.vec.test <- rexp(100000,rate=intensity)
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

matrix.time.pop.sim.test <- vec.pop.sim.test(age.sim.vec.test,gender.sim.vec.test,start.year.test,end.year.test,arrival.sim.vec.test)
time.pop.sim.test <- matrix.time.pop.sim.test[,"t.p"]
delta.p.test <- matrix.time.pop.sim.test[,"delta.p"]

## Simulate excess survival times 

x.matrix.test<- cbind(age.sim.vec.test,gender.sim.vec.test,treatment.sim.vec.test)

### Out of control right away on excess

time.excess.sim.test <- weibull.excess.sim.time.proportional(lambda,phi,1.25,beta.vec,x.matrix.test)

time.censoring.sim.test <- cbind(rexp(n.sim.test,rate=0.001),end.year.test-start.year.test-arrival.sim.vec.test)
time.censoring.sim.test <- apply(time.censoring.sim.test,1,FUN = min)

## Intermediate data set

inter.sim.dataset.test <- cbind(x.matrix.test,matrix.time.pop.sim.test,time.excess.sim.test,time.censoring.sim.test)

## Pick out minimum of the three different survival times as the observed time. 
## Also, need to make a death indicator (delta.i) for all causes. 
time.observed.sim.test <- pmin(time.pop.sim.test,time.excess.sim.test,time.censoring.sim.test)
delta.i.test <- pmax(delta.p.test,as.numeric(time.excess.sim.test<time.pop.sim.test))*as.numeric(time.censoring.sim.test>time.observed.sim.test)

# Final data set

final.data.set.test <- cbind(x.matrix.test,start.year.test+arrival.sim.vec.test,time.observed.sim.test,delta.i.test)
colnames(final.data.set.test) <- c("age","gender","treatment","year","time.obs","delta.i")
final.data.set.test <- as.data.frame(final.data.set.test)

t.grid <- seq(from=0,to=10,by=0.01) # a grid of time points to calculate the CUSUM chart

cusum.chart <- cusum_r.t.EM.loess(start.year.test,age.sim.vec.test,gender.sim.vec.test,x.matrix.test,time.observed.sim.test
                                          ,arrival.sim.vec.test,delta.i.test,estimated.beta,1.25,t.grid,baseline.haz.vec.event.times,
                                          cumulative.baseline.haz.vec.event.times,time.vec.events,span.value.lambda.0 = 0.5,span.value.Lambda.0 = 0.5,
                                          pop.data.male,pop.data.female,2020)
plot(t.grid,cusum.chart[,1],type="l",xlab="t (in years)",ylab=expression(paste(Psi,"(t)")),main=expression( paste("CUSUM chart for a Weibull baseline when ",eta,"=0 and ",rho,"=1.25")))

# Apply on test data to see how it works when the hazard stays in control throughout the follow-up 

arrival.sim.vec.test.2 <- rexp(100000,rate=intensity)
arrival.sim.vec.test.2 <- cumsum(arrival.sim.vec.test.2)
arrival.sim.vec.test.2 <- arrival.sim.vec.test.2[arrival.sim.vec.test.2<=10]
n.sim.test.2 <- length(arrival.sim.vec.test.2)

## Simulate age
age.sim.vec.test.2 <- rnorm(n.sim.test.2,mean=70,sd=10) # simulate age from N(mean=70,sigma=10). 

## Simulate gender

gender.sim.vec.test.2 <- sample(0:1,size = n.sim.test.2,replace=T) # simulate gender (0=male, 1=female)

## Simulate treatment group 

treatment.sim.vec.test.2 <- sample(0:1,size = n.sim.test.2,replace=T) 

## Simulate population survival times and population censoring indicator 

matrix.time.pop.sim.test.2 <- vec.pop.sim.test(age.sim.vec.test.2,gender.sim.vec.test.2,start.year.test,end.year.test,arrival.sim.vec.test.2)
time.pop.sim.test.2 <- matrix.time.pop.sim.test.2[,"t.p"]
delta.p.test.2 <- matrix.time.pop.sim.test.2[,"delta.p"]

## Simulate excess survival times 

x.matrix.test.2<- cbind(age.sim.vec.test.2,gender.sim.vec.test.2,treatment.sim.vec.test.2)

### Out of control right away on excess

time.excess.sim.test.2 <- weibull.excess.sim.time(lambda,phi,beta.vec,x.matrix.test.2)

time.censoring.sim.test.2 <- cbind(rexp(n.sim.test.2,rate=0.001),end.year.test-start.year.test-arrival.sim.vec.test.2)
time.censoring.sim.test.2 <- apply(time.censoring.sim.test.2,1,FUN = min)

## Intermediate data set

inter.sim.dataset.test.2 <- cbind(x.matrix.test.2,matrix.time.pop.sim.test.2,time.excess.sim.test.2,time.censoring.sim.test.2)

## Pick out minimum of the three different survival times as the observed time. 
## Also, need to make a death indicator (delta.i) for all causes. 
time.observed.sim.test.2 <- pmin(time.pop.sim.test.2,time.excess.sim.test.2,time.censoring.sim.test.2)
delta.i.test.2 <- pmax(delta.p.test.2,as.numeric(time.excess.sim.test.2<time.pop.sim.test.2))*as.numeric(time.censoring.sim.test.2>time.observed.sim.test.2)

# Final data set

final.data.set.test.2 <- cbind(x.matrix.test.2,start.year.test+arrival.sim.vec.test.2,time.observed.sim.test.2,delta.i.test.2)
colnames(final.data.set.test.2) <- c("age","gender","treatment","year","time.obs","delta.i")
final.data.set.test.2 <- as.data.frame(final.data.set.test.2)

t.grid <- seq(from=0,to=10,by=0.01)

cusum.chart.2 <- cusum_r.t.EM.loess(start.year.test,age.sim.vec.test.2,gender.sim.vec.test.2,x.matrix.test.2,time.observed.sim.test.2
                                          ,arrival.sim.vec.test.2,delta.i.test.2,estimated.beta,1.25,t.grid,baseline.haz.vec.event.times,
                                          cumulative.baseline.haz.vec.event.times,time.vec.events,span.value = 0.5,pop.data.male,pop.data.female,2020)
plot(t.grid,cusum.chart.2[,1],type="l",xlab="t (in years)",ylab=expression(paste(Psi,"(t)")),main=expression( paste("CUSUM chart for a Weibull baseline when in-control and ",rho,"=1.25")))


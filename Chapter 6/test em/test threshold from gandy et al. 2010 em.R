library(dplyr)
library(relsurv)
source("Funksjoner.r")
source("Nevents.r")
source("faster and approx general population simulation.R")

# Comments can be found from the piecewise case - similar idea here

lambda <- 0.75
phi <- 0.005
rho <- 1.25
start.year <- 2010
end.year <- 2020
beta.vec <- c(0.05,0.1,0.5)
max.events <- 100
max.R.em.1 <- numeric(1000)
n.times.events.stopping.em.1 <- numeric(1000)
number.em.1 <- numeric(1000)
max.R.em.2 <- numeric(1000)
n.times.events.stopping.em.2 <- numeric(1000)
number.em.2 <- numeric(1000)
t.grid <- seq(from=0,to=10,by=0.01)
t.grid.for.haz <- seq(from=0.01,to=10,by=0.01)
Lambda0 <- phi*t.grid.for.haz^lambda
lambda0 <- lambda*phi*t.grid.for.haz^(lambda-1)/365.241
t.grid.in.days <- t.grid.for.haz*365.241
threshold_c <- getthresholdPNevents(rho=1.25,Nmax=max.events, PN0=0.10) # threshold
set.seed(420)

for(i in 1:1000){
  intensity <- 100
  arrival.sim.vec <- rexp(2000,rate=intensity)
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
  time.excess.sim <- weibull.excess.sim.time(lambda,phi,beta.vec,x.matrix)
  
  time.censoring.sim<- cbind(rexp(n.sim,rate=0.001),end.year-start.year-arrival.sim.vec)
  time.censoring.sim <- apply(time.censoring.sim,1,FUN = min)
  
  ## Intermediate data set
  
  
  ## Pick out minimum of the three different survival times as the observed time. 
  ## Also, need to make a death indicator (delta.i) for all causes. 
  time.observed.sim <- pmin(time.pop.sim,time.excess.sim,time.censoring.sim)
  delta.i <- pmax(delta.p,as.numeric(time.excess.sim<time.pop.sim))*as.numeric(time.censoring.sim>time.observed.sim)
  
  cusum.curve <- cusum_r.t.EM.loess(start.year,age.sim.vec,gender.sim.vec,x.matrix,time.observed.sim
                                            ,arrival.sim.vec,delta.i,beta.vec,1.25,t.grid,lambda0,
                                            Lambda0,t.grid.in.days,span.value.lambda0 = 0.5,span.value.Lambda0=0.5,
                                            pop.data.male,pop.data.female,2020)
  max.R.em.1[i] <- max(cusum.curve[,1])
  index.stopping <- min(which(cusum.curve[,1]>=threshold_c))
  if(index.stopping==Inf){
    n.times.events.stopping.em.1[i] <- FALSE
  } else {
    n.events.stopping <- cusum.curve[index.stopping,2]
    number.em.1[i] <- n.events.stopping
    if(n.events.stopping<=max.events){
      n.times.events.stopping.em.1[i] <- T
    } else {
      n.times.events.stopping.em.1[i] <- F
    }
  }
  
  cusum.curve.2 <- cusum_r.t.EM.smoothing.spline(start.year,age.sim.vec,gender.sim.vec,x.matrix,time.observed.sim
                                    ,arrival.sim.vec,delta.i,beta.vec,1.25,t.grid,lambda0,
                                    Lambda0,t.grid.in.days,df.lambda0=5,df.Lambda0=5,pop.data.male,pop.data.female,2020)
  max.R.em.2[i] <- max(cusum.curve.2[,1])
  index.stopping.2 <- min(which(cusum.curve.2[,1]>=threshold_c))
  if(index.stopping.2==Inf){
    n.times.events.stopping.em.2[i] <- FALSE
  } else {
    n.events.stopping.2 <- cusum.curve.2[index.stopping.2,2]
    number.em.2[i] <- n.events.stopping.2
    if(n.events.stopping.2<=max.events){
      n.times.events.stopping.em.2[i] <- T
    } else {
      n.times.events.stopping.em.2[i] <- F
    }
  }
}
sum(n.times.events.stopping.em.2)/1000

max.R.em.1.young <- numeric(1000)
n.times.events.stopping.em.1.young <- numeric(1000)
number.em.1.young <- numeric(1000)
max.R.em.2.young <- numeric(1000)
n.times.events.stopping.em.2.young <- numeric(1000)
number.em.2.young <- numeric(1000)
for(i in 1:1000){
  intensity <- 100
  arrival.sim.vec <- rexp(2000,rate=intensity)
  arrival.sim.vec <- cumsum(arrival.sim.vec)
  arrival.sim.vec <- arrival.sim.vec[arrival.sim.vec<=10]
  n.sim <- length(arrival.sim.vec)
  
  ## Simulate age
  age.sim.vec<- rnorm(n.sim,mean=40,sd=5) # simulate age from N(mean=40,sigma=5) 
  
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
  time.excess.sim <- weibull.excess.sim.time(lambda,phi,beta.vec,x.matrix)
  
  time.censoring.sim<- cbind(rexp(n.sim,rate=0.001),end.year-start.year-arrival.sim.vec)
  time.censoring.sim <- apply(time.censoring.sim,1,FUN = min)
  
  ## Intermediate data set
  
  
  ## Pick out minimum of the three different survival times as the observed time. 
  ## Also, need to make a death indicator (delta.i) for all causes. 
  time.observed.sim <- pmin(time.pop.sim,time.excess.sim,time.censoring.sim)
  delta.i <- pmax(delta.p,as.numeric(time.excess.sim<time.pop.sim))*as.numeric(time.censoring.sim>time.observed.sim)
  
  cusum.curve <- cusum_r.t.EM.loess(start.year,age.sim.vec,gender.sim.vec,x.matrix,time.observed.sim
                                            ,arrival.sim.vec,delta.i,beta.vec,1.25,t.grid,lambda0,
                                            Lambda0,t.grid.in.days,span.value.lambda0 = 0.5,span.value.Lambda0=0,
                                            pop.data.male,pop.data.female,2020)
  max.R.em.1.young[i] <- max(cusum.curve[,1])
  index.stopping <- min(which(cusum.curve[,1]>=threshold_c))
  if(index.stopping==Inf){
    n.times.events.stopping.em.1.young[i] <- FALSE
  } else {
    n.events.stopping <- cusum.curve[index.stopping,2]
    number.em.1.young[i] <- n.events.stopping
    if(n.events.stopping<=max.events){
      n.times.events.stopping.em.1.young[i] <- T
    } else {
      n.times.events.stopping.em.1.young[i] <- F
    }
  }
  
  cusum.curve.2 <- cusum_r.t.EM.smoothing.spline(start.year,age.sim.vec,gender.sim.vec,x.matrix,time.observed.sim
                                                 ,arrival.sim.vec,delta.i,beta.vec,1.25,t.grid,lambda0,
                                                 Lambda0,t.grid.in.days,df.lambda0=5,df.Lambda0=5,pop.data.male,pop.data.female,2020)
  max.R.em.2.young[i] <- max(cusum.curve.2[,1])
  index.stopping.2 <- min(which(cusum.curve.2[,1]>=threshold_c))
  if(index.stopping.2==Inf){
    n.times.events.stopping.em.2.young[i] <- FALSE
  } else {
    n.events.stopping.2 <- cusum.curve.2[index.stopping.2,2]
    number.em.2.young[i] <- n.events.stopping.2
    if(n.events.stopping.2<=max.events){
      n.times.events.stopping.em.2.young[i] <- T
    } else {
      n.times.events.stopping.em.2.young[i] <- F
    }
  }
}

library(dplyr)
library(relsurv)
source("Funksjoner.r")
source("Nevents.r")
source("faster and approx general population simulation.R")

set.seed(42)

### First example ### 
max.R <- rep(NA,1000)
n.times.events.stopping <- rep(NA,1000)
number <- rep(NA,1000)
max.events <- 100
start.year <- 2010
end.year <- 2020
chi.vec <- -c(7,6.75,6.5,6.25,6,5.75,5.5,5.75)
t.vec <- c(0,1,2,3,4,5,10,15)
beta.vec <- c(0.05,0.1,0.5)
rho=1.25
t.grid <- seq(from=0,to=10,by=0.01)

for(i in 1:1000){
  lambda <- 100
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
  time.excess.sim <- vec.sim.excess.piecewise(t.vec,chi.vec,x.matrix,beta.vec)
  
  time.censoring.sim<- cbind(rexp(n.sim,rate=0.001),end.year-start.year-arrival.sim.vec)
  time.censoring.sim <- apply(time.censoring.sim,1,FUN = min)
  
  ## Intermediate data set
  
  
  ## Pick out minimum of the three different survival times as the observed time. 
  ## Also, need to make a death indicator (delta.i) for all causes. 
  time.observed.sim <- pmin(time.pop.sim,time.excess.sim,time.censoring.sim)
  delta.i <- pmax(delta.p,as.numeric(time.excess.sim<time.pop.sim))*as.numeric(time.censoring.sim>time.observed.sim)
  
  cusum.curve <- cusum_r.t.piecewise(start.year,age.sim.vec,gender.sim.vec,x.matrix,time.observed.sim,arrival.sim.vec,delta.i,chi.vec,t.vec,beta.vec,rho,t.grid,pop.data.male,pop.data.female,2020)
  max.R[i] <- max(cusum.curve[,1])
  index.stopping <- min(which(cusum.curve[,1]>=getthresholdPNevents(rho=1.25,Nmax=max.events, PN0=0.05))) # check which time step in the grid where the chart crosses the threshold
  if(index.stopping==Inf){ # if it does not cross at all 
    n.times.events.stopping[i] <- FALSE
  } else {
    n.events.stopping <- cusum.curve[index.stopping,2] # check number of events when the chart signals 
    number[i] <- n.events.stopping
    if(n.events.stopping<=max.events){ # if it actually signals when the number of events is less than Nmax defined 
      n.times.events.stopping[i] <- T
    } else {
      n.times.events.stopping[i] <- F
    }
  }
}

### Change now to a younger sample 

max.R.younger <- rep(NA,1000)
n.times.events.stopping.younger <- rep(NA,1000)
number.younger <- rep(NA,1000)

for(i in 1:1000){
  lambda <- 100
  arrival.sim.vec <- rexp(2000,rate=lambda)
  arrival.sim.vec <- cumsum(arrival.sim.vec)
  arrival.sim.vec <- arrival.sim.vec[arrival.sim.vec<=10]
  n.sim <- length(arrival.sim.vec)
  
  ## Simulate age
  age.sim.vec<- rnorm(n.sim,mean=50,sd=5) # simulate age from N(mean=70,sigma=10). 
  
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
  time.excess.sim <- vec.sim.excess.piecewise(t.vec,chi.vec,x.matrix,beta.vec)
  
  time.censoring.sim<- cbind(rexp(n.sim,rate=0.001),end.year-start.year-arrival.sim.vec)
  time.censoring.sim <- apply(time.censoring.sim,1,FUN = min)
  
  ## Intermediate data set
  
  
  ## Pick out minimum of the three different survival times as the observed time. 
  ## Also, need to make a death indicator (delta.i) for all causes. 
  time.observed.sim <- pmin(time.pop.sim,time.excess.sim,time.censoring.sim)
  delta.i <- pmax(delta.p,as.numeric(time.excess.sim<time.pop.sim))*as.numeric(time.censoring.sim>time.observed.sim)
  
  cusum.curve <- cusum_r.t.piecewise(start.year,age.sim.vec,gender.sim.vec,x.matrix,time.observed.sim,arrival.sim.vec,delta.i,chi.vec,t.vec,beta.vec,rho,t.grid,pop.data.male,pop.data.female,2020)
  max.R.younger[i] <- max(cusum.curve[,1])
  index.stopping <- min(which(cusum.curve[,1]>=getthresholdPNevents(rho=1.25,Nmax=max.events, PN0=0.05)))
  if(index.stopping==Inf){
    n.times.events.stopping.younger[i] <- FALSE
  } else {
    n.events.stopping <- cusum.curve[index.stopping,2]
    number.younger[i] <- n.events.stopping
    if(n.events.stopping<=max.events){
      n.times.events.stopping.younger[i] <- T
    } else {
      n.times.events.stopping.younger[i] <- F
    }
  }
}

# New parameter vector - chi_2

set.seed(42)

max.R.2 <- rep(NA,1000)
n.times.events.stopping.2 <- rep(NA,1000)
number.2 <- rep(NA,1000)
max.events <- 100
start.year <- 2010
end.year <- 2020
chi.vec.2 <- -c(7,6.75,6.5,6.25,6,5.75,5.5,5.75)+1
t.vec <- c(0,1,2,3,4,5,10,15)
beta.vec <- c(0.05,0.1,0.5)
rho=1.25
t.grid <- seq(from=0,to=10,by=0.01)

for(i in 1:1000){
  lambda <- 100
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
  time.excess.sim <- vec.sim.excess.piecewise(t.vec,chi.vec.2,x.matrix,beta.vec)
  
  time.censoring.sim<- cbind(rexp(n.sim,rate=0.001),end.year-start.year-arrival.sim.vec)
  time.censoring.sim <- apply(time.censoring.sim,1,FUN = min)
  
  ## Intermediate data set
  
  
  ## Pick out minimum of the three different survival times as the observed time. 
  ## Also, need to make a death indicator (delta.i) for all causes. 
  time.observed.sim <- pmin(time.pop.sim,time.excess.sim,time.censoring.sim)
  delta.i <- pmax(delta.p,as.numeric(time.excess.sim<time.pop.sim))*as.numeric(time.censoring.sim>time.observed.sim)
  
  cusum.curve <- cusum_r.t.piecewise(start.year,age.sim.vec,gender.sim.vec,x.matrix,time.observed.sim,arrival.sim.vec,delta.i,chi.vec.2,t.vec,beta.vec,rho,t.grid,pop.data.male,pop.data.female,2020)
  max.R.2[i] <- max(cusum.curve[,1])
  index.stopping <- min(which(cusum.curve[,1]>=getthresholdPNevents(rho=1.25,Nmax=max.events, PN0=0.05)))
  if(index.stopping==Inf){
    n.times.events.stopping.2[i] <- FALSE
  } else {
    n.events.stopping <- cusum.curve[index.stopping,2]
    number.2[i] <- n.events.stopping
    if(n.events.stopping<=max.events){
      n.times.events.stopping.2[i] <- T
    } else {
      n.times.events.stopping.2[i] <- F
    }
  }
}

### Change now to a younger sample 

max.R.younger.2 <- rep(NA,1000)
n.times.events.stopping.younger.2 <- rep(NA,1000)
number.younger.2 <- rep(NA,1000)

for(i in 1:1000){
  lambda <- 100
  arrival.sim.vec <- rexp(2000,rate=lambda)
  arrival.sim.vec <- cumsum(arrival.sim.vec)
  arrival.sim.vec <- arrival.sim.vec[arrival.sim.vec<=10]
  n.sim <- length(arrival.sim.vec)
  
  ## Simulate age
  age.sim.vec<- rnorm(n.sim,mean=50,sd=5) # simulate age from N(mean=70,sigma=10). 
  
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
  time.excess.sim <- vec.sim.excess.piecewise(t.vec,chi.vec.2,x.matrix,beta.vec)
  
  time.censoring.sim<- cbind(rexp(n.sim,rate=0.001),end.year-start.year-arrival.sim.vec)
  time.censoring.sim <- apply(time.censoring.sim,1,FUN = min)
  
  ## Intermediate data set
  
  
  ## Pick out minimum of the three different survival times as the observed time. 
  ## Also, need to make a death indicator (delta.i) for all causes. 
  time.observed.sim <- pmin(time.pop.sim,time.excess.sim,time.censoring.sim)
  delta.i <- pmax(delta.p,as.numeric(time.excess.sim<time.pop.sim))*as.numeric(time.censoring.sim>time.observed.sim)
  
  cusum.curve <- cusum_r.t.piecewise(start.year,age.sim.vec,gender.sim.vec,x.matrix,time.observed.sim,arrival.sim.vec,delta.i,chi.vec.2,t.vec,beta.vec,rho,t.grid,pop.data.male,pop.data.female,2020)
  max.R.younger.2[i] <- max(cusum.curve[,1])
  index.stopping <- min(which(cusum.curve[,1]>=getthresholdPNevents(rho=1.25,Nmax=max.events, PN0=0.05)))
  if(index.stopping==Inf){
    n.times.events.stopping.younger.2[i] <- FALSE
  } else {
    n.events.stopping <- cusum.curve[index.stopping,2]
    number.younger.2[i] <- n.events.stopping
    if(n.events.stopping<=max.events){
      n.times.events.stopping.younger.2[i] <- T
    } else {
      n.times.events.stopping.younger.2[i] <- F
    }
  }
}

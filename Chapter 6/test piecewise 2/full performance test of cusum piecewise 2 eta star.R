### Test of CUSUM performance for piecewise example nr.2 ### 

source("faster approx general population simulation.R")
source("Funksjoner.R")

vec.sim.excess.piecewise <- function(t.vec,chi.vec,x.matrix,beta.vec){ # in control simulation excess
  
  cum.baseline.excess.hazard.func <- function(t){
    if(length(t.vec)>1){
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
    } else {
      return(exp(chi.vec)*t)
    }
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
    if(length(t.vec)>1){
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
    } else {
      return(rho*exp(chi.vec)*t)
    }
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

### First: Simulation to decide threshold c. Want to have a 5% probability of false alarm throughout the follow-up. Rho=1.25 for this example and lambda=100.

set.seed(42)

max.R <- rep(NA,1000) # store maximum value of each simulated CUSUM chart
start.year <- 2010
end.year <- 2020
chi.vec <- -c(5,10,11,12) # the chosen baseline parameters 
t.vec <- c(0,0.25,5,10) # the chosen splitting of follow-up interval 
beta.vec <- c(0.05,0.1,0.5)
rho=1.25
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
  time.excess.sim <- vec.sim.excess.piecewise(t.vec,chi.vec,x.matrix,beta.vec)
  
  time.censoring.sim<- cbind(rexp(n.sim,rate=0.001),end.year-start.year-arrival.sim.vec)
  time.censoring.sim <- apply(time.censoring.sim,1,FUN = min)
  
  ## Intermediate data set
  
  
  ## Pick out minimum of the three different survival times as the observed time. 
  ## Also, need to make a death indicator (delta.i) for all causes. 
  time.observed.sim <- pmin(time.pop.sim,time.excess.sim,time.censoring.sim)
  delta.i <- pmax(delta.p,as.numeric(time.excess.sim<time.pop.sim))*as.numeric(time.censoring.sim>time.observed.sim)
  
  cusum.curve <- cusum_r.t.piecewise(start.year,age.sim.vec,gender.sim.vec,x.matrix,time.observed.sim,
                                 arrival.sim.vec, delta.i,chi.vec,t.vec,beta.vec,rho,t.grid,
                                 pop.data.male,pop.data.female,2020)
  max.R[i] <- max(cusum.curve[,1])
}

quantiles.max.R <- quantile(max.R,probs=c(0,0.9,0.95,0.99))
threshold.c <- quantiles.max.R[3] # upper 5% quantile under in control as threshold, i.e 5% of false alarm. 

### Second: Testing in a case where eta^*=5 such that only patients arriving after 5 years of monitoring will experience the out-of-control hazard.  

max.test.star <- rep(NA,1000)
prop.excess.deaths.star <- rep(NA,1000)
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
  
  time.excess.sim.test <- numeric(n.sim.test)
  index.before.eta <- which(arrival.sim.vec.test<=5)
  index.after.eta <- -(index.before.eta)
  time.excess.sim.test[index.before.eta] <- vec.sim.excess.piecewise(t.vec,chi.vec,x.matrix.test[index.before.eta,],beta.vec) # arrival before 5 years --> only experience in-control hazard
  time.excess.sim.test[index.after.eta] <- vec.sim.excess.piecewise.proportional(t.vec,chi.vec,x.matrix.test[index.after.eta,],beta.vec,rho)

  time.censoring.sim.test <- cbind(rexp(n.sim.test,rate=0.001),end.year-start.year-arrival.sim.vec.test)
  time.censoring.sim.test <- apply(time.censoring.sim.test,1,FUN = min)
  
  ## Intermediate data set
  
  inter.sim.dataset.test <- cbind(x.matrix.test,matrix.time.pop.sim.test,time.excess.sim.test,time.censoring.sim.test)
  
  ## Pick out minimum of the three different survival times as the observed time. 
  ## Also, need to make a death indicator (delta.i) for all causes. 
  time.observed.sim.test <- pmin(time.pop.sim.test,time.excess.sim.test,time.censoring.sim.test)
  delta.i.test <- pmax(delta.p.test,as.numeric(time.excess.sim.test<time.pop.sim.test))*as.numeric(time.censoring.sim.test>time.observed.sim.test)

  # Final data set
  
  prop.excess.deaths.star[i] <- sum(time.excess.sim.test<=time.pop.sim.test & delta.i.test==1)/sum(delta.i.test==1)

  
  t.grid <- seq(from=0,to=10,by=0.01)
  R.out.control.fast <- cusum_r.t.piecewise(start.year,age.sim.vec.test,gender.sim.vec.test,x.matrix.test,time.observed.sim.test,arrival.sim.vec.test,delta.i.test,chi.vec,t.vec,beta.vec,rho,t.grid,
                                            pop.data.male,pop.data.female,2020)
  max.test.star[i] <- max(R.out.control.fast[,1])

}

exceed.threshold.star <- max.test.star>threshold.c
sum(exceed.threshold.star)


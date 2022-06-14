### FULL TEST OF EDERER 2 WITH THE EXTREME CHOICE OF PARAMETER VALUES ### 

library(survival)
library(relsurv)
library(dplyr)
source("Funksjoner.R")
set.seed(420)

lambda = 0.5
phi=0.0005
beta.vec <- c(0.05,0.1)

age.sim.vec <- rnorm(10000,mean=70,sd=10) # Running 1000 simulations so that everything is a bit faster
age.sim.vec <- round(age.sim.vec)
n.sim <- length(age.sim.vec)
## Simulate gender

gender.sim.vec <- sample(0:1,size = n.sim,replace=T) # simulate gender (0=male, 1=female)

## Simulate population survival times and population censoring indicator 

matrix.time.pop.sim <- vec.pop.sim(age.sim.vec,rep(2010,n.sim),gender.sim.vec)
time.pop.sim <- matrix.time.pop.sim[,"t.p"]
delta.p <- matrix.time.pop.sim[,"delta.p"]

## Simulate excess survival times 

x.matrix <- cbind(age.sim.vec,gender.sim.vec)
time.excess.sim <- weibull.excess.sim.time(lambda,phi,beta.vec,x.matrix)

## Simulate interim censoring time - minimum of a random number from exp(rate=0.001) or maximum follow-up time which is almost 11 years in this case as we start out at the beginning of 2010 to end of 2020. 

time.censoring.sim <- cbind(rexp(n.sim,rate=0.001),rep(max(time.pop.sim),n.sim))
time.censoring.sim <- apply(time.censoring.sim,1,FUN = min)

## Intermediate data set

inter.sim.dataset <- cbind(x.matrix,matrix.time.pop.sim,time.excess.sim,time.censoring.sim)

## Pick out minimum of the three different survival times as the observed time. Also, need to make a death indicator (delta.i) for all causes. 

time.matrix.sim <- cbind(time.pop.sim,time.excess.sim,time.censoring.sim)
time.observed.sim <- apply(time.matrix.sim,1,FUN=min)
delta.i <- rep(0,n.sim)
for(i in 1:n.sim){
  if(time.observed.sim[i]==time.excess.sim[i]){
    delta.i[i] <- 1 # death indicator=1 if dead, 0 if censored
  }
  if(time.observed.sim[i]==time.pop.sim[i] & delta.p[i]==1){
    delta.i[i] <- 1
  }
  if(time.observed.sim[i]==time.pop.sim[i] & delta.p[i]==0){
    delta.i[i] <- 0
  }
  if(time.observed.sim[i]==time.censoring.sim[i]){
    delta.i[i] <- 0
  }
}

## Final standard data set with only observed survival times and death indicator delta.1 

final.data.set <- cbind(x.matrix,rep(2010,n.sim),time.observed.sim,delta.i)
colnames(final.data.set) <- c("age","gender","year","time.obs","delta.i")
final.data.set <- as.data.frame(final.data.set)


## Population part of Ederer 2 method 

numerator.func.ederer2.individual <- function(i){ # individual i
  force(i)
  if(final.data.set$gender[i]==0){
    life.table <- pop.data.male
  } else {
    life.table <- pop.data.female
  }
  haz.data <- haz.output(final.data.set$age[i],final.data.set$year[i],life.table)
  
  haz.func.single.int <- function(j){
    force(j)
    haz.value <- function(t){
      return((final.data.set$time.obs[i] >= t)*(j-1<=t)*(t<j)*haz.data[j,2]) # Y_i(u) lambda_{Pi} between two different years 
    }
  }
  list.haz.func.single <- (sapply(1:nrow(haz.data),haz.func.single.int))
  full.list.haz.func.single <- function(t){
    sum(unlist(sapply(list.haz.func.single,function(f) f(t)))) # sum over all yearly intervals to cover the whole follow-up
  }
}  

numerator.total.list <- sapply(1:nrow(final.data.set),numerator.func.ederer2.individual) # a list containing all the terms in numerator
total.numerator.ederer <- function(t){
  sum(unlist(sapply(numerator.total.list,function(f) f(t)))) # the whole sum in the numerator 
}

denominator.ederer.func <- function(i){
  force(i)
  indicator <- function(t)(final.data.set$time.obs[i] >= t)*1 # denominator for individual i
}

denominator.total.list <- sapply(1:nrow(final.data.set),denominator.ederer.func) 
total.denominator.ederer <- function(t){
  sum(unlist(sapply(denominator.total.list,function(f) f(t)))) # the whole sum in the denominator 
}

j.indicator <- function(t){
  indicator <- (total.denominator.ederer(t)>0)*1 # indicator function in case Y(t)=0
  return(indicator)
}
full.integrand <- function(t){
  ifelse(j.indicator(t)==0,0,total.numerator.ederer(t)/total.denominator.ederer(t)) # Define s.t J(t)/Y(t)=0 when Y(t)=0
}
vec.full.integrand <- Vectorize(full.integrand)

# Nelson-Aalen/Overall hazard part 
surv.fit.test <- summary(survfit(Surv(time.observed.sim,delta.i)~1))
nelson.aalen <- c(0,cumsum(surv.fit.test$n.event/surv.fit.test$n.risk))
time.vec <- c(0,surv.fit.test$time)

nelson.aalen.func.t <- function(t){ # The overall part at a time t 
  i <- max(which(time.vec <= t))
  return(nelson.aalen[i])
}

# Test

step.size.integral <- 0.01 # accuracy, smaller will be more accurate but longer time. 
t.grid.for.real.obs <- seq(from=0,to=11,by=0.1)
t.grid.for.int <- seq(from=0,to=11,by=step.size.integral)
# integrate.vec.full.integrand <- function(i){ # this takes too much time
#   integrate(vec.full.integrand,lower=t.grid.for.int[i],upper=t.grid.for.int[i+1])$value # integrate small intervals between each time from time grid for more precision
# }
# test.int <- c(0,sapply(1:(length(t.grid.for.int)-1),integrate.vec.full.integrand))
# cumsum.pop.test.est <- cumsum(test.int) # the integral from 0 to t for values in t.grid 
# nelson.aalen.est <- sapply(t.grid,nelson.aalen.func.t)
# obs.net.est <- nelson.aalen.est-cumsum.pop.test.est # Ederer 2 estimate 
real.obs.net <- weibull.obs.net.fast.ver(lambda,phi,beta.vec,x.matrix,2010,t.grid.for.real.obs) # "real" observable net 
plot(t.grid.for.real.obs,real.obs.net,type="l",ylim=c(0.9,1),xlab="time",ylab="survival")
# lines(t.grid.for.int,exp(-obs.net.est),col="red") # plot of what is supposed to be the estimated observable net survival curve

# Test with midpoint method for the integral of the population part 
midpoint.pop.haz.int <- function(i){ # element i from grid
  arg <- (t.grid.for.int[i]+t.grid.for.int[i+1])/2
  return(step.size.integral*full.integrand(arg))
}
test.int <- c(0,sapply(1:(length(t.grid.for.int)-1),midpoint.pop.haz.int))
cumsum.pop.test.est <- cumsum(test.int) # the integral of the population part from 0 to t for values in t.grid 
nelson.aalen.est <- sapply(t.grid.for.int,nelson.aalen.func.t)
obs.net.est <- nelson.aalen.est-cumsum.pop.test.est # Ederer 2 estimate
nortab <- transrate.hmd(male="NOR.mltper_1x1.txt",female="NOR.fltper_1x1.txt")
final.data.set$gender <- final.data.set$gender+1 
ederer2.estimate <- rs.surv(Surv(time.obs*365.241,delta.i)~1,data=final.data.set,ratetable = nortab,method="ederer2",add.times=c(0.5)*365.241,rmap = list(age=age*365.241,sex=gender,year=(year-1970)*365.241))
lines(ederer2.estimate$time/365.241,ederer2.estimate$surv,col="orange",lty=2)
lines(t.grid.for.int,exp(-obs.net.est),col="red")
legend("topleft",legend = c("Observable net","Ederer 2 from relsurv","Ederer 2 self-implemented"),
       col=c("black","orange","red"),lty=c(1,2,1))

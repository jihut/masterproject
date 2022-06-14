### Example of net survival vs observable net survival (Ederer 2 vs PP)

rm(list=ls()) 

library(dplyr)
library(relsurv)
library(popEpi)
library(Epi)
source("Funksjoner.r")
##### First case: Beta vector is the zero vector - lambda=1, phi=0.0025 #####
set.seed(42)
lambda <- 1
phi <- 0.0025
beta.vec <- c(0,0)
start.year <- 2000
end.year <- 2021

age.sim.vec <- rnorm(10000,mean=70,sd=10) # simulate age from N(mean=70,sigma=10). 
age.sim.vec <- round(age.sim.vec)
n.sim <- length(age.sim.vec)

## Simulate gender

gender.sim.vec <- sample(0:1,size = n.sim,replace=T) # simulate gender (0=male, 1=female)

## Simulate population survival times and population censoring indicator 

matrix.time.pop.sim <- vec.pop.sim(age.sim.vec,rep(start.year,n.sim),gender.sim.vec)
time.pop.sim <- matrix.time.pop.sim[,"t.p"]
delta.p <- matrix.time.pop.sim[,"delta.p"]

## Simulate excess survival times 

x.matrix <- cbind(age.sim.vec,gender.sim.vec)
time.excess.sim <- weibull.excess.sim.time(lambda,phi,beta.vec,x.matrix)

time.censoring.sim <- cbind(rexp(n.sim,rate=0.001),rep(max(time.pop.sim),n.sim))
time.censoring.sim <- apply(time.censoring.sim,1,FUN = min)

## Intermediate data set

inter.sim.dataset <- cbind(x.matrix,matrix.time.pop.sim,time.excess.sim,time.censoring.sim)

## Pick out minimum of the three different survival times as the observed time. 
## Also, need to make a death indicator (delta.i) for all causes. 
time.observed.sim <- pmin(time.pop.sim,time.excess.sim,time.censoring.sim)
delta.i <- pmax(delta.p,as.numeric(time.excess.sim<time.pop.sim))*as.numeric(time.censoring.sim>time.observed.sim)

# Final data set

final.data.set <- cbind(x.matrix,rep(start.year,n.sim),time.observed.sim,delta.i)
colnames(final.data.set) <- c("age","gender","year","time.obs","delta.i")
final.data.set <- as.data.frame(final.data.set)

sum(delta.i==1)
length(which(delta.i==1 & time.excess.sim<=time.pop.sim))
length(which(delta.i==1 & time.excess.sim>=time.pop.sim))
length(which(delta.i==1 & time.excess.sim<=time.pop.sim))/sum(delta.i==1)


# Calculate net survival 
t.grid <- seq(from=0,to=end.year-start.year,by=0.1)

net.survival.individual <- weibull.prop.excess.surv.test(lambda,phi,beta.vec,x.matrix,t.grid) # each row corresponds to an individual, i.e S_{Ei}. The columns correspond to the time steps defined by the time grid.
net.survival.group <- colMeans(net.survival.individual) # take the column mean to get the overall (mariginal) net survival of a group at the given value of time. 
plot(t.grid,net.survival.group,type="l",xlab="Time (in years)",ylab="Survival",lwd=2,main="Net vs Observable Net",ylim=c(0,800))

# Ederer 2 vs PP 

nortab <- transrate.hmd(male="NOR.mltper_1x1.txt",female="NOR.fltper_1x1.txt")
final.data.set$gender <- final.data.set$gender+1 
ederer2.estimate <- rs.surv(Surv(time.obs*365.241,delta.i)~1,data=final.data.set,ratetable = nortab,method="ederer2",rmap = list(age=age*365.241,sex=gender,year=mdy.date(1,1,year)))
lines(ederer2.estimate$time/365.241,ederer2.estimate$surv,col="orange",lty=2,lwd=2)
pp.estimate <- rs.surv(Surv(time.obs*365.241,delta.i)~1,data=final.data.set,ratetable = nortab,method="pohar-perme",rmap = list(age=age*365.241,sex=gender,year=mdy.date(1,1,year)))
lines(pp.estimate$time/365.241,pp.estimate$surv,col="green",lty=2,lwd=2)

legend("topright",legend = c("Net Survival","Ederer 2 with relsurv","PP with relsurv"),
       col=c("black","orange","green"),lty=c(1,2,2),lwd=2)

### Test against survtab ### 

final.data.set.survtab <- cbind(x.matrix,rep(start.year,n.sim),rep(start.year,n.sim)+time.observed.sim,delta.i)
colnames(final.data.set.survtab) <- c("AGE","sex","start.year","exit","delta.i")
final.data.set.survtab <- as.data.frame(final.data.set.survtab)


# relevant population table 
pop.male.table <- cbind.data.frame(rep(0,nrow(pop.data.male)),pop.data.male$Year,pop.data.male$Age,pop.data.male$mx)
colnames(pop.male.table) <- c("sex","year","age","haz")
pop.female.table <- cbind.data.frame(rep(1,nrow(pop.data.male)),pop.data.female$Year,pop.data.female$Age,pop.data.female$mx)
colnames(pop.female.table) <- c("sex","year","age","haz")
poptable.final <- rbind.data.frame(pop.male.table,pop.female.table)
poptable.final$age <- as.integer(poptable.final$age)
poptable.final$sex <- as.integer(poptable.final$sex)
poptable.final[is.na(poptable.final)] <- 110
x.sim <- Lexis(
  entry = list(FUT = 0,age=AGE,year=start.year), 
  exit = list(year = exit), 
  data = final.data.set.survtab,
  exit.status = delta.i,
  merge = TRUE
)

x.sim[,which(names(x.sim) %in% c("start.year", "AGE"))] <- NULL

st.e2 <- survtab(
  Surv(time = FUT, event = lex.Xst) ~ 1,
  data = x.sim,
  surv.type = "surv.rel",
  relsurv.method = "e2",
  breaks = list(FUT = seq(0, 21, 1/12)),
  pophaz = poptable.final
)

st.pp <- survtab(
  Surv(time = FUT, event = lex.Xst) ~ 1,
  data = x.sim,
  surv.type = "surv.rel",
  relsurv.method = "pp",
  breaks = list(FUT = seq(0, 21, 1/12)),
  pophaz = poptable.final
)

plot(t.grid,net.survival.group,type="l",xlab="Time (in years)",ylab="Survival",lwd=2,main="Net vs Observable Net",ylim=c(0.8,1.2))
lines(st.e2$Tstart,st.e2$r.e2,type="l",col="orange",lty=2,lwd=2)
lines(st.pp$Tstart,st.pp$r.pp,type="l",col="green",lty=2,lwd=2)


### Test against self-implemented function ### 

# J(t)

at.risk.indicator.individual <- function(i){
  force(i)
  indicator <- function(t)(final.data.set$time.obs[i] >= t)*1 # denominator for individual i
}
at.risk.list.func <- sapply(1:nrow(final.data.set),at.risk.indicator.individual) 
total.at.risk.func <- function(t){
  sum(unlist(sapply(at.risk.list.func,function(f) f(t)))) # the whole sum in the denominator 
}

j.indicator <- function(t){
  indicator <- (total.at.risk.func(t)>0)*1 # indicator function in case Y(t)=0
  return(indicator)
}

# Observed part of PP estimator

numerator.obs.pp.individual <- function(i){
  force(i)
  if(final.data.set$gender[i]==0){
    life.table <- pop.data.male
  } else {
    life.table <- pop.data.female
  }
  haz.data <- haz.output(final.data.set$age[i],final.data.set$year[i],life.table) # get out the hazard table for observation nr.i 
  numerator.func <- function(t){
    if(0<= t & t<=1){
      return(as.numeric(t == final.data.set$time.obs[i] & final.data.set$delta.i[i]==1)/exp(-(t-haz.data[1,1])*haz.data[1,2])) # basically N^w(u) - only non zero when t equal or larger than the corresponding uncensored time
    }
    if(t==nrow(haz.data)){
      return(as.numeric(t == final.data.set$time.obs[i] & final.data.set$delta.i[i]==1)/exp(-sum(haz.data[,2]))) # if t is a natural number
    } else {
      k <- max(which(c(1:nrow(haz.data)) <= t))
      return(as.numeric(t == final.data.set$time.obs[i] & final.data.set$delta.i[i]==1)/exp(-sum(haz.data[1:(k),2])-(t-haz.data[k+1,1])*haz.data[k+1,2])) # otherwise
    }
  }
}

numerator.obs.pp.total.list <- sapply(1:nrow(final.data.set),numerator.obs.pp.individual)   
total.numerator.obs.pp <- function(t){
  sum(unlist(sapply(numerator.obs.pp.total.list,function(f) f(t)))) # the whole sum in the numerator over all observations
}

denominator.obs.pp.individual <- function(i){ # similar as before, but for the denominator 
  force(i)
  if(final.data.set$gender[i]==0){
    life.table <- pop.data.male
  } else {
    life.table <- pop.data.female
  }
  haz.data <- haz.output(final.data.set$age[i],final.data.set$year[i],life.table)
  denominator.func <- function(t){
    if(0<= t & t<=1){
      return(as.numeric(final.data.set$time.obs[i] >= t)/exp(-(t-haz.data[1,1])*haz.data[1,2]))
    }
    if(t==nrow(haz.data)){
      return(as.numeric(final.data.set$time.obs[i] >= t)/exp(-sum(haz.data[,2])))
    } else {
      k <- max(which(c(1:nrow(haz.data)) <= t))
      return(as.numeric(final.data.set$time.obs[i] >= t)/exp(-sum(haz.data[1:(k),2])-(t-haz.data[k+1,1])*haz.data[k+1,2]))
    }
  }
}
denominator.obs.pp.total.list <- sapply(1:nrow(final.data.set),denominator.obs.pp.individual)  
total.denominator.obs.pp <- function(t){
  sum(unlist(sapply(denominator.obs.pp.total.list,function(f) f(t)))) # the whole sum in the denominator 
}

full.obs.function <- function(t){
  return(total.numerator.obs.pp(t)/total.denominator.obs.pp(t)) # the full integrand in the first term
}

index <- which(final.data.set$time.obs<21 & final.data.set$delta.i==1) # pick out which observations who experience an event 
event.time.obs <- final.data.set$time.obs[index]
sort.event.time.obs <- sort(event.time.obs) # sort after increasing times to event 
pp.nelson.aalen <- sapply(sort.event.time.obs,full.obs.function) # evaluate at all the times to event 
pp.cumsum.nelson.aalen <- c(0,cumsum(pp.nelson.aalen)) # sum them all over to get the observed part
time.vec <- c(0,sort.event.time.obs)

pp.nelson.aalen.func.t <- function(t){ 
  i <- max(which(time.vec <= t))
  return(pp.cumsum.nelson.aalen[i])
}

# Population part of PP method 
numerator.pop.pp.individual <- function(i){
  force(i)
  if(final.data.set$gender[i]==0){
    life.table <- pop.data.male
  } else {
    life.table <- pop.data.female
  }
  haz.data <- haz.output(final.data.set$age[i],final.data.set$year[i],life.table)
  haz.func.single.int <- function(j){ # population hazard function of an observation in a given interval
    force(j)
    haz.value <- function(t){
      return((final.data.set$time.obs[i] >= t)*(j-1<=t)*(t<j)*haz.data[j,2])
    }
  }
  list.haz.func.single <- (sapply(1:nrow(haz.data),haz.func.single.int))
  full.list.haz.func.single <- function(t){
    sum(unlist(sapply(list.haz.func.single,function(f) f(t)))) # sum all over all intervals to get the population hazard over the follow-up interval of the patient
  }
  numerator.func <- function(t){
    if(0<= t & t<=1){
      return(as.numeric(final.data.set$time.obs[i] >= t)*full.list.haz.func.single(t)/exp(-(t-haz.data[1,1])*haz.data[1,2])) # Y_i/S_{Pi}, the rest follows like before 
    }
    if(t==nrow(haz.data)){
      return(as.numeric(final.data.set$time.obs[i] >= t)*full.list.haz.func.single(t)/exp(-sum(haz.data[,2])))
    } else {
      k <- max(which(c(1:nrow(haz.data)) <= t))
      return(as.numeric(final.data.set$time.obs[i] >= t)*full.list.haz.func.single(t)/exp(-sum(haz.data[1:(k),2])-(t-haz.data[k+1,1])*haz.data[k+1,2]))
    }
  }
}
numerator.pop.pp.total.list <- sapply(1:nrow(final.data.set),numerator.pop.pp.individual)  
total.numerator.pop.pp <- function(t){
  sum(unlist(sapply(numerator.pop.pp.total.list,function(f) f(t)))) # the whole sum in the numerator 
}

full.integrand.pop.pp <- function(t){
  ifelse(j.indicator(t)==0,0,total.numerator.pop.pp(t)/total.denominator.obs.pp(t)) # integrand of the second term 
}

# Test

step.size.integral <- 0.1 # step size of midpoint as integrate function takes too much time now
t.grid.for.int <- seq(from=0,to=21,by=step.size.integral)
midpoint.pop.haz.int <- function(i){ # element i from grid, using the midpoint method to calculate the integral numerically
  arg <- (t.grid.for.int[i]+t.grid.for.int[i+1])/2
  return(step.size.integral*full.integrand.pop.pp(arg))
}
test.int <- c(0,sapply(1:(length(t.grid.for.int)-1),midpoint.pop.haz.int))
cumsum.pop.test.est <- cumsum(test.int) # the integral from 0 to t for values in t.grid of the population part 
pp.nelson.aalen.est <- sapply(t.grid.for.int,pp.nelson.aalen.func.t)
net.est <- pp.nelson.aalen.est-cumsum.pop.test.est # PP estimate 
lines(t.grid.for.int,exp(-net.est),col="red",lty=2,lwd=2)
legend("topright",legend = c("Net Survival","Ederer 2 with survtab","PP with survtab","PP self-implemented"),
       col=c("black","orange","green","red"),lty=c(1,2,2,2),lwd=2)

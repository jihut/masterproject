### FULL TEST OF PP ESTIMATOR WITH THE EXTREME CHOICE OF PARAMETER VALUES ### 
### Specifically for the simulated data set from the script "full ederer 2 implementation for the extreme case.R"

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
  haz.data <- haz.output(final.data.set$age[i],final.data.set$year[i],life.table)
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
  sum(unlist(sapply(numerator.obs.pp.total.list,function(f) f(t)))) # the whole sum in the numerator 
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
  sum(unlist(sapply(denominator.obs.pp.total.list,function(f) f(t)))) # the whole sum in the numerator 
}

full.obs.function <- function(t){
  return(total.numerator.obs.pp(t)/total.denominator.obs.pp(t)) # full integrand for the first term 
}

index <- which(final.data.set$time.obs<11 & final.data.set$delta.i==1) # pick out which observations who experience an event 
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
  ifelse(j.indicator(t)==0,0,total.numerator.pop.pp(t)/total.denominator.obs.pp(t))
}

# Test

step.size.integral <- 0.01
t.grid.for.int <- seq(from=0,to=11,by=step.size.integral)
midpoint.pop.haz.int <- function(i){ # element i from grid
  arg <- (t.grid.for.int[i]+t.grid.for.int[i+1])/2
  return(step.size.integral*full.integrand.pop.pp(arg))
}
test.int <- c(0,sapply(1:(length(t.grid.for.int)-1),midpoint.pop.haz.int))
cumsum.pop.test.est <- cumsum(test.int) # the integral from 0 to t for values in t.grid 
pp.nelson.aalen.est <- sapply(t.grid.for.int,pp.nelson.aalen.func.t)
net.est <- pp.nelson.aalen.est-cumsum.pop.test.est # PP estimate 
net.survival.individual <- weibull.prop.excess.surv.test(lambda,phi,beta.vec,x.matrix,t.grid.for.int) # each row corresponds to an individual, i.e S_{Ei}. The columns correspond to the time steps defined by the time grid.
net.survival.group <- colMeans(net.survival.individual) # take the column mean to get the overall (mariginal) net survival of a group at the given value of time. 
plot(t.grid.for.int,net.survival.group,type="l",ylim=c(0.8,1))
lines(t.grid.for.int,exp(-net.est),col="red")

# Ederer 2 function

# final.data.set = a simulated data set - needs to contain following columns with the given names: age, year and time.obs and delta.i

final.data.set 
numerator.func.ederer2.individual <- function(i){ # individual i
  force(i)
  if(final.data.set$gender[i]==0){
    life.table <- pop.data.male
  } else {
    life.table <- pop.data.female
  }
  haz.data <- haz.output(final.data.set$age[i],final.data.set$year[i],life.table)
  
  haz.func.single.int <- function(j){ # Y_i(u) lambda_{Pi} between two different years 
    force(j)
    haz.value <- function(t){
      return((final.data.set$time.obs[i] >= t)*(j-1<=t)*(t<j)*haz.data[j,2])
    }
  }
  list.haz.func.single <- (sapply(1:nrow(haz.data),haz.func.single.int)) 
  full.list.haz.func.single <- function(t){
    sum(unlist(sapply(list.haz.func.single,function(f) f(t)))) # sum over all yearly intervals to cover the whole follow-up
  }
}  

numerator.total.list <- sapply(1:nrow(final.data.set),numerator.func.ederer2.individual) 
total.numerator.ederer <- function(t){
  sum(unlist(sapply(numerator.total.list,function(f) f(t)))) # sum over all individuals to get the numerator 
}

denominator.ederer.func <- function(i){
  force(i)
  indicator <- function(t)(final.data.set$time.obs[i] >= t)*1
}

denominator.total.list <- sapply(1:nrow(final.data.set),denominator.ederer.func)
total.denominator.ederer <- function(t){
  sum(unlist(sapply(denominator.total.list,function(f) f(t)))) # sum over all individuals to get the denominator
}

j.indicator <- function(t){
  indicator <- (total.denominator.ederer(t)>0)*1 # J(t)
  return(indicator)
}
full.integrand <- function(t){
  ifelse(j.indicator(t)==0,0,total.numerator.ederer(t)/total.denominator.ederer(t)) # the integrand corresponding to the population part
}
vec.full.integrand <- Vectorize(full.integrand)

sort.obs.time <- c(0,sort(final.data.set$time.obs,decreasing = F)) # sort after observed time in the data set 

pop.integral.between.obs.times.func <- function(k){
  integrate(vec.full.integrand,lower=sort.obs.time[k],upper=sort.obs.time[k+1])$value # integral of population part between two consecutive observed times
}

pop.integral.between.obs.times <- cumsum(sapply(1:nrow(final.data.set),pop.integral.between.obs.times.func)) 

## Now need to use this to find at any given value of t by using pop.integral.between.obs.times

final.est.pop.cum.hazard <- function(t){
  i <- max(which(sort.obs.time <= t)) # find the index of the max ordered observed time less than t
  if(t==0){return(0)}
  if(i==1){
    return(integrate(vec.full.integrand,lower=0,upper=t)$value)
  } else {
    return(pop.integral.between.obs.times[i-1]+integrate(vec.full.integrand,lower=sort.obs.time[i],upper=t)$value) # integral of the second term in the Ederer 2 at a time t
  }
}

# 1. Function to simulate from population hazard

pop.data.male <- read.table("NOR.mltper_1x1.txt",skip=2,header=T) # Life table for male in Norway
pop.data.female <- read.table("NOR.fltper_1x1.txt",skip=2,header=T) # Life table for female in Norway

haz.output <- function(Age,Year,haz.data){
  
  # a function to simulate from population hazard when age and year change at the same time (in the beginning of a year)
  
  step <- which(haz.data$Year==(1950) & haz.data$Age==(42))-which(haz.data$Year==(1950-1) & haz.data$Age==(42-1)) # step size in table for observations less than 110 years old at a specific year
  step.110 <- which(haz.data$Year==(1950) & haz.data$Age==("110+"))-which(haz.data$Year==(1949) & haz.data$Age==("110+"))
  
  
  
  if(Age > 109){ # for life table from mortality.org if Age is 110+
    index <- which(haz.data$Year==(Year) & haz.data$Age==("110+")) # row index of hazard for time=0 in the original dataset.
    n <- floor((nrow(haz.data)-index)/step.110)
    
    haz.matrix <- matrix(0,nrow=n+1,ncol=2) # matrix to store the hazards for the individual
    haz.matrix[,1] <- c(0:n) # from t=0 to t=n Years
    
    for(i in 1:(n+1)){
      haz.matrix[i,2] <- haz.data[index+(i-1)*step.110,"mx"] # pick out the hazard from life table for Year i.
    }
    
  }
  
  
  if(Age <= 109){
    if(Age+haz.data[nrow(haz.data),"Year"]-Year > 110){
      index <- which(haz.data$Year==(Year) & haz.data$Age==(Age)) 
      year.diff <- haz.data[nrow(haz.data),"Year"]-Year
      age.diff <- 110 - Age # new amount of steps down the life table to get the next hazard when Age >109
      
      haz.matrix <- matrix(rep(0,2*(year.diff+1)),ncol=2)
      haz.matrix[,1] <- c(0:year.diff) 
      for(i in 1:(age.diff+1)){
        haz.matrix[i,2] <- haz.data[index+(i-1)*step,"mx"]
      }
      
      new.index <- index + age.diff*step
      for(i in (age.diff+2):(year.diff+1)){
        haz.matrix[i,2] <- haz.data[new.index+(i-(age.diff+1))*step.110,"mx"]
      }
    } 
    
    
    if(Age+haz.data[nrow(haz.data),"Year"]-Year <= 110) {
      
      index <- which(haz.data$Year==(Year) & haz.data$Age==(Age)) # row index of hazard for time=0.
      n <- as.numeric(floor((nrow(haz.data)-index)/step)) # lower limit of the last yearly follow up interval, i.e nrow(haz.data)=n+1 is the maximum follow up in years.
      
      haz.matrix <- matrix(0,nrow=n+1,ncol=2) # matrix to store the hazards for the individual
      haz.matrix[,1] <- c(0:n) # from t=0 to t=n Years
      for(i in 1:(n+1)){
        haz.matrix[i,2] <- haz.data[index+(i-1)*step,"mx"] # pick out the hazard from life table for Year i.
      }
    }
    
  }
  
  return(haz.matrix)
}

one.sim <- function(Age,Year,haz.data){
  haz.matrix <- haz.output(Age,Year,haz.data)
  n <- haz.matrix[nrow(haz.matrix),1]
  surv.matrix <- matrix(0,nrow=n+1,ncol=2) # matrix to store the survivor functions at
  surv.matrix[,1] <- c(0:n) # from t=0 to t=n Years
  surv.matrix[,2] <- exp(-cumsum(haz.matrix[,2])) # survivor functions.
  
  
  #inverse transform 
  u <- runif(1)
  if(u < surv.matrix[n+1,2]){
    t.p <- n+1
    delta.p <- 0
  }
  if(u> surv.matrix[1,2]){
    t.1 <- surv.matrix[1,1]
    t.2 <- surv.matrix[2,1]
    t.p <- t.1 + (u-1)*(t.2-t.1)/(surv.matrix[1,2]-1) # linear interpolation
    delta.p <- 1
  } else {
    for (j in 1:(n)){
      if(surv.matrix[j,2] >u & u>surv.matrix[j+1,2]){
        t.1 <- surv.matrix[j,1]+1
        t.2 <- surv.matrix[j,1]+2
        t.p <- t.1 + (u-surv.matrix[j,2])*(t.2-t.1)/(surv.matrix[j+1,2]-surv.matrix[j,2]) # linear interpolation
        delta.p <- 1
        break
      }
    }
  }
  return(c(t.p,delta.p))
}

vec.pop.sim <- function(Age,Year,Gender){ # Gender: 0=male, 1=female
  time.mat <- matrix(rep(0,2*length(Age)),ncol=2)
  colnames(time.mat) <- c("t.p","delta.p")
  for (i in 1:length(Age)){
    if(Gender[i]==0){
      haz.data <- pop.data.male  
      time.mat[i,] <- one.sim(Age[i],Year[i],haz.data)
    } else {
      haz.data <- pop.data.female  
      time.mat[i,] <- one.sim(Age[i],Year[i],haz.data)
    }
  }
  return(time.mat)
} 

# 2. Function to simulate from Weibull excess hazard

weibull.excess.one.sim <- function(lambda,phi,beta.vec,x.vec){
  u <- runif(1)
  inner.prod <- as.numeric(beta.vec%*%x.vec)
  phi.i <- phi*exp(inner.prod)
  t <- (-log(1-u)/phi.i)^(1/lambda) # inverse transform 
  return(t)
}

weibull.excess.sim.time <- function(lambda,phi,beta.vec,x.matrix){ # lambda=shape, phi=scale of baseline, beta.vec=parameter vector, x.matrix=covariate matrix where each row corresponds to a single observation.
  time.vec <- rep(0,nrow(x.matrix))
  for(i in 1:nrow(x.matrix)){
    time.vec[i] <- weibull.excess.one.sim(lambda,phi,beta.vec,x.matrix[i,])
  }
  return(time.vec)
}




# 4. Real net survival 

weibull.prop.excess.surv <- function(lambda,phi,beta.vec,x.vec,t){
  inner.prod <- as.numeric(exp(beta.vec%*%x.vec))
  hazard <- function(t){
    return(lambda*phi*t^(lambda-1)*inner.prod)
  }
  surv.func <- function(t){
    int.value <- integrate(hazard,lower=0,upper=t)$value
    return(exp(-int.value))
  }
  vec.surv.func <- Vectorize(surv.func)
}

weibull.prop.excess.surv.test <- function(lambda,phi,beta.vec,x.matrix,t.vec){
  surv.matrix <- matrix(nrow=nrow(x.matrix),ncol=length(t.vec))
  distinct.x.matrix <- distinct(as.data.frame(x.matrix)) # a matrix containing the different and distinct combinations of covariates
  distinct.x.matrix <- as.matrix(distinct.x.matrix)
  for(i in 1:nrow(distinct.x.matrix)){
    
    inner.prod <- as.numeric(exp(beta.vec%*%distinct.x.matrix[i,]))
    hazard <- function(t){
      return(lambda*phi*t^(lambda-1)*inner.prod)
    }
    
    surv.func <- function(t){
      int.value <- phi*inner.prod*t^lambda
      return(exp(-int.value))
    }
    
    vec.surv.func <- Vectorize(surv.func)
    surv.values <- vec.surv.func(t.vec)
    index <- which(x.matrix[,1]==distinct.x.matrix[i,1] & x.matrix[,2]==distinct.x.matrix[i,2]) # find which rows of the original covariate matrix with this combination of covariate values
    surv.matrix[index,] <- matrix(rep(surv.values,length(index)),nrow=length(index),byrow=T) # these rows will then have surv.values as the values of individual net survival at the different time steps defined by t.vec
    
  }
  return(surv.matrix)
}




# 5. Real "observable" net survival / cause-specific survival 

weibull.obs.net.fast.ver <- function(lambda,phi,beta.vec,x.matrix,Year,t.vec){ # Note that first column of x.matrix needs to be age, second column gender. 
  decoy.haz.matrix <- haz.output(x.matrix[1,1],Year,pop.data.male) # random test matrix to find out the amount of Yearly follow-up intervals
  
  distinct.x.matrix <- distinct(as.data.frame(x.matrix)) # a matrix containing the different and distinct combinations of covariates
  distinct.x.matrix <- as.matrix(distinct.x.matrix)
  
  numerator.list.time.step <- list() # an empty list storing the function in the numerator of the group observable net hazard for each at each Yearly interval
  denominator.list.time.step <- list()
  
  numerator.function.individual.1 <- function(n){ # numerator for t between 0 and 1 years
    force(n)
    if(distinct.x.matrix[n,2]==0){
      life.table <- pop.data.male
    } else {
      life.table <- pop.data.female
    }
    haz.data <- haz.output(distinct.x.matrix[n,1],Year,life.table)
    index <- which(x.matrix[,1] == distinct.x.matrix[n,1] & x.matrix[,2] == distinct.x.matrix[n,2])
    amount <- length(index) # amount of observations with this specific combination of covariate values
    # Individual net survival
    inner.prod <- as.numeric(exp(beta.vec%*%distinct.x.matrix[n,]))
    excess.hazard <- function(t){
      return(lambda*phi*t^(lambda-1)*inner.prod)
    }
    surv.excess.function <- function(t){
      int.value <- phi*inner.prod*t^lambda
      return(exp(-int.value))
    }
    vec.surv.excess.function <- Vectorize(surv.excess.function,"t")
    
    # Population individual survival
    cum.pop.hazard<- function(t){ 
      return((t-haz.data[1,1])*haz.data[1,2]) # cumulative population hazard between 0 and 1 when the population hazard is a step function
    }
    surv.pop.function <- function(t){
      return(exp(-cum.pop.hazard(t)))
    }
    vec.surv.pop.function <- Vectorize(surv.pop.function)
    
    # Overall individual survival 
    surv.overall.function<- function(t){ 
      return(vec.surv.pop.function(t)*vec.surv.excess.function(t))
    }
    vec.surv.overall.function <- Vectorize(surv.overall.function)
    
    # The numerator in the formula of observable net hazard for these x-values
    numerator.function <- function(t){
      return(amount*vec.surv.overall.function(t)*excess.hazard(t))
    }
    vec.numerator.function <- Vectorize(numerator.function)
  }
  
  denominator.function.individual.1 <- function(n){
    force(n)
    if(distinct.x.matrix[n,2]==0){
      life.table <- pop.data.male
    } else {
      life.table <- pop.data.female
    }
    haz.data <- haz.output(distinct.x.matrix[n,1],Year,life.table)
    index <- which(x.matrix[,1] == distinct.x.matrix[n,1] & x.matrix[,2] == distinct.x.matrix[n,2])
    amount <- length(index)
    # Individual net survival
    inner.prod <- as.numeric(exp(beta.vec%*%distinct.x.matrix[n,]))
    excess.hazard <- function(t){
      return(lambda*phi*t^(lambda-1)*inner.prod)
    }
    surv.excess.function <- function(t){
      int.value <- phi*inner.prod*t^lambda
      return(exp(-int.value))
    }
    vec.surv.excess.function <- Vectorize(surv.excess.function,"t")
    
    # Population individual survival
    cum.pop.hazard<- function(t){ 
      return((t-haz.data[1,1])*haz.data[1,2]) 
    }
    
    surv.pop.function <- function(t){
      return(exp(-cum.pop.hazard(t)))
    }
    vec.surv.pop.function <- Vectorize(surv.pop.function)
    
    # Overall individual survival = denominator as well
    surv.overall.function<- function(t){ 
      return(amount*vec.surv.pop.function(t)*vec.surv.excess.function(t))
    }
    vec.surv.overall.function <- Vectorize(surv.overall.function)
    
  }
  
  numerator.func.1 <- sapply(1:nrow(distinct.x.matrix),numerator.function.individual.1) # a list to store the numerators in the first interval 
  denominator.func.1 <- sapply(1:nrow(distinct.x.matrix),denominator.function.individual.1)
  
  numerator.list.time.step[[1]] <- numerator.func.1 # add this list to the "big" list 
  denominator.list.time.step[[1]] <- denominator.func.1
  
  numerator.function.individual <- function(n,k){ # the same as above for the remaining years 
    force(n)
    force(k)
    if(distinct.x.matrix[n,2]==0){
      life.table <- pop.data.male
    } else {
      life.table <- pop.data.female
    }
    haz.data <- haz.output(distinct.x.matrix[n,1],Year,life.table)
    index <- which(x.matrix[,1] == distinct.x.matrix[n,1] & x.matrix[,2] == distinct.x.matrix[n,2])
    amount <- length(index)
    # Individual net survival
    inner.prod <- as.numeric(exp(beta.vec%*%distinct.x.matrix[n,]))
    excess.hazard <- function(t){
      return(lambda*phi*t^(lambda-1)*inner.prod)
    }
    surv.excess.function <- function(t){
      int.value <- phi*inner.prod*t^lambda
      return(exp(-int.value))
    }
    vec.surv.excess.function <- Vectorize(surv.excess.function,"t")
    
    # Population individual survival
    cum.pop.hazard<- function(t){ 
      return(sum(haz.data[1:(k-1),2])+(t-haz.data[k,1])*haz.data[k,2]) 
    }
    surv.pop.function <- function(t){
      return(exp(-cum.pop.hazard(t)))
    }
    vec.surv.pop.function <- Vectorize(surv.pop.function)
    
    surv.overall.function<- function(t){ # vectorize both population and excess to get correct result with test function 
      return(vec.surv.pop.function(t)*vec.surv.excess.function(t))
    }
    vec.surv.overall.function <- Vectorize(surv.overall.function)
    
    numerator.function <- function(t){
      return(amount*vec.surv.overall.function(t)*excess.hazard(t))
    }
    vec.numerator.function <- Vectorize(numerator.function)
  }
  
  denominator.function.individual <- function(n,k){
    force(n)
    force(k)
    if(distinct.x.matrix[n,2]==0){
      life.table <- pop.data.male
    } else {
      life.table <- pop.data.female
    }
    haz.data <- haz.output(distinct.x.matrix[n,1],Year,life.table)
    index <- which(x.matrix[,1] == distinct.x.matrix[n,1] & x.matrix[,2] == distinct.x.matrix[n,2])
    amount <- length(index)
    # Individual net survival
    inner.prod <- as.numeric(exp(beta.vec%*%distinct.x.matrix[n,]))
    excess.hazard <- function(t){
      return(lambda*phi*t^(lambda-1)*inner.prod)
    }
    surv.excess.function <- function(t){
      int.value <- phi*inner.prod*t^lambda
      return(exp(-int.value))
    }
    vec.surv.excess.function <- Vectorize(surv.excess.function,"t")
    
    # Population individual survival
    cum.pop.hazard<- function(t){ 
      return(sum(haz.data[1:(k-1),2])+(t-haz.data[k,1])*haz.data[k,2]) 
    }
    surv.pop.function <- function(t){
      return(exp(-cum.pop.hazard(t)))
    }
    vec.surv.pop.function <- Vectorize(surv.pop.function)
    
    surv.overall.function<- function(t){ # vectorize both population and excess to get correct result with test function 
      return(amount*vec.surv.pop.function(t)*vec.surv.excess.function(t))
    }
    vec.surv.overall.function <- Vectorize(surv.overall.function)
  }
  
  for(i in 2:nrow(decoy.haz.matrix)){ # each time step 
    numerator.func <- list()
    denominator.func <- list()
    for(j in 1:nrow(distinct.x.matrix)){
      numerator.func[[j]] <- numerator.function.individual(j,i) # numerator in this time period for covariate pattern nr. j
      denominator.func[[j]] <- denominator.function.individual(j,i)
    }
    numerator.list.time.step[[i]] <- numerator.func
    denominator.list.time.step[[i]] <- denominator.func
  }
  
  relevant.integrand.time.interval.i <- function(i){
    relevant.numerator.list <- numerator.list.time.step[[i]] # take out the list of numerator functions i-th time interval 
    relevant.denominator.list <- denominator.list.time.step[[i]]
    total.numerator<- function(t){
      sum(unlist(lapply(relevant.numerator.list,function(f) f(t)))) # the numerator is the sum of all such functions (S_Oi * lambda_Ei)
    }
    vec.total.numerator <- Vectorize(total.numerator)
    
    total.denominator<- function(t){
      sum(unlist(lapply(relevant.denominator.list,function(f) f(t))))
    }
    vec.total.denominator <- Vectorize(total.denominator)
    
    integrand <- function(t){
      return(vec.total.numerator(t)/(vec.total.denominator(t))) # the observable net hazard in this interval itself 
    }
    vec.integrand <- Vectorize(integrand)
  }
  
  integral.between.time.intervals.func <- function(i){ # a function to integrate the relevant function between each time step 
    error.check <- try(integrate(relevant.integrand.time.interval.i(i),lower=i-1,upper=i),silent=T)
    loadError <- (is(error.check, 'try-error')|is(error.check,'error'))
    return(ifelse(loadError==T,999999,integrate(relevant.integrand.time.interval.i(i),lower=i-1,upper=i)$value))
  }
  integral.between.time.intervals <- sapply(1:nrow(decoy.haz.matrix),integral.between.time.intervals.func) # do this for all time step
  cumulative.integral.between.time.intervals <- cumsum(integral.between.time.intervals) # add them all to together 
  
  final.observable.net.hazard <- function(t){
    i <- floor(t) # check the value of t to the nearest year 
    if(t==0){return(0)}
    if(0<t & t<1){
      return(integrate(relevant.integrand.time.interval.i(i+1),lower=i,upper = t)$value)
    }
    if(i==t){
      return(cumulative.integral.between.time.intervals[i]) # if t is an integer, i.e. exactly a number of years
    } else {
      error.check <- try(integrate(relevant.integrand.time.interval.i(i+1),lower=i,upper=t),silent=T)
      loadError <- (is(error.check, 'try-error')|is(error.check,'error'))
      return(ifelse(loadError==T,999999,integrate(relevant.integrand.time.interval.i(i+1),lower=i,upper=t)$value)+cumulative.integral.between.time.intervals[i])
    }
  }
  
  return(exp(-sapply(t.vec,final.observable.net.hazard)))
}

ederer.2.func <- function(data.set,t.grid){ # data.set needs to contain the following columns with following names: "age","gender","year","time.obs","delta.i"
  ## Population part of Ederer 2 method  
  
  numerator.func.ederer2.individual <- function(i){ # individual i
    force(i)
    if(data.set$gender[i]==0){
      life.table <- pop.data.male
    } else {
      life.table <- pop.data.female
    }
    haz.data <- haz.output(data.set$age[i],data.set$year[i],life.table)
    
    haz.func.single.int <- function(j){
      force(j)
      haz.value <- function(t){
        return((data.set$time.obs[i] >= t)*(j-1<=t)*(t<j)*haz.data[j,2])
      }
    }
    list.haz.func.single <- (sapply(1:nrow(haz.data),haz.func.single.int))
    full.list.haz.func.single <- function(t){
      sum(unlist(sapply(list.haz.func.single,function(f) f(t))))
    }
  }  
  
  numerator.total.list <- sapply(1:nrow(data.set),numerator.func.ederer2.individual) # a list containing all the terms in numerator
  total.numerator.ederer <- function(t){
    sum(unlist(sapply(numerator.total.list,function(f) f(t)))) # the whole sum in the numerator 
  }
  
  denominator.ederer.func <- function(i){
    force(i)
    indicator <- function(t)(data.set$time.obs[i] >= t)*1 # denominator for individual i
  }
  
  denominator.total.list <- sapply(1:nrow(data.set),denominator.ederer.func) 
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
  
  # Nelson-Aalen part 
  surv.fit.test <- summary(survfit(Surv(data.set$time.obs,delta.i)~1))
  nelson.aalen <- c(0,cumsum(surv.fit.test$n.event/surv.fit.test$n.risk))
  time.vec <- c(0,surv.fit.test$time)
  
  nelson.aalen.func.t <- function(t){
    i <- max(which(time.vec <= t))
    return(nelson.aalen[i])
  }
  
  step.size <- t.grid[2]-t.grid[1]
  midpoint.pop.haz.int <- function(i){ # element i from grid
    arg <- (t.grid[i]+t.grid[i+1])/2
    return(step.size.integral*full.integrand(arg))
  }
  test.int <- c(0,sapply(1:(length(t.grid)-1),midpoint.pop.haz.int))
  cumsum.pop.test.est <- cumsum(test.int) # the integral from 0 to t for values in t.grid 
  nelson.aalen.est <- sapply(t.grid.for.int,nelson.aalen.func.t)
  obs.net.est <- nelson.aalen.est-cumsum.pop.test.est
  return(exp(-obs.net.est))
  
}

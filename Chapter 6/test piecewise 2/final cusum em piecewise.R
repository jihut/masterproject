# Final CUSUM for piecewise constant baseline excess hazard models #

# Function to calculate X_i(t) for any time t after origin=start.year

x_i.vec <- function(time.obs.vec,arrival.vec,u){ # X_i(t) - at-risk time for each observation after a time t with start.year as t=0.
  
  return(pmin(time.obs.vec,u-arrival.vec))
  
}

# Event indicator at a given time at risk.

delta_i.t <- function(time.obs.vec,arrival.vec,delta.i.vec,t){
  delta_i.vec.t <- as.numeric(time.obs.vec<=x_i.vec(time.obs.vec,arrival.vec,t) & delta.i.vec==1)
  return(delta_i.vec.t)
}

# Function to evaluate the excess hazard of an individual with a piecewise constant baseline.

# partition.t.vec = partition of monitoring period. E.g. if monitoring period is 10 years and the partition is 0-2.5,2.5-5,5-7.5,7.5-10 --> partition.t.vec = c(0,2.5,5,7.5) - no need for the end time, this should be specified as the last value of t.grid. 

excess.hazard.piecewise.func<- function(partition.t.vec,chi.vec,x.matrix,beta.vec,u){ #lambda_Ei
  
  baseline.vec <- exp(chi.vec)
  
  baseline.excess.hazard.func <- function(s){ # baseline part of lambda_Ei(t) when we have a piecewise constant baseline hazard model
    index <- findInterval(s,partition.t.vec)
    exp(chi.vec[index])
  }
  
  baseline.excess.hazard.func(u)*as.numeric(exp(beta.vec %*% t(x.matrix)))
  
}

# Function to evaluate cumulative excess hazard with a piecewise constant baseline

cumulative.excess.hazard.piecewise.func<- function(partition.t.vec,chi.vec,x.matrix,beta.vec,u){ #Lambda_Ei
  interval.length.partition.t.vec <- diff(partition.t.vec)
  num.interval.length.partition.t.vec <- length(interval.length.partition.t.vec)
  baseline.vec <- exp(chi.vec)
  cumulative.baseline.vec <- cumsum(baseline.vec[1:num.interval.length.partition.t.vec]*diff(partition.t.vec))
  new.cumulative.baseline.vec <- c(0,cumulative.baseline.vec)
  cumulative.baseline.excess.hazard.func <- function(s){
    index <- findInterval(s,t.vec)
    new.cumulative.baseline.vec[index]+(s-t.vec[index])*exp(chi.vec)[index]
  }
  cumulative.baseline.excess.hazard.func(u)*as.numeric(exp(beta.vec %*% t(x.matrix)))
}

# CUSUM function

### start.year = Calendar year when the monitoring starts
### age.vec = A vector containing age at arrival for the observations used in the monitoring
### gender.vec = A vector related to gender, 0=male and 1=female
### x.matrix = covariate matrix 
### time.obs.vec = observed survival times of the individuals
### arrival.vec = number of years until the patients arrive after start year
### delta.i.vec = censoring indicator when monitoring period ends, 0=censored and 1=death
### chi.vec = baseline parameters - usually estimated by one of the piecewise constant baseline models
### partition.t.vec = mentioned before
### beta.vec = estimated parameters from EM-based model related to the covariates. 
### rho = proportionality constant in the proportional alternative
### t.grid = a grid of points over the monitoring period - first element needs to be 0 (basically t=0, i.e. the time point corresponding to start.year)
### pop.data.male = life table for male from Human Mortality Database
### pop.data.female = life table for female from Human Mortality Database
### end.year.table = final year of the life tables 

cusum_r.t.piecewise <- function(start.year,age.vec,gender.vec,x.matrix,time.obs.vec,
                                arrival.vec, delta.i.vec,chi.vec,partition.t.vec,beta.vec,rho,t.grid,
                                pop.data.male,pop.data.female,end.year.table){
  new.pop.data.male <- pop.data.male
  new.pop.data.male$Age[new.pop.data.male$Age=="110+"] <- 110
  new.pop.data.male$Age <- as.integer(new.pop.data.male$Age)
  new.pop.data.female <- pop.data.female
  new.pop.data.female$Age[new.pop.data.female$Age=="110+"] <- 110
  new.pop.data.female$Age <- as.integer(new.pop.data.female$Age)
  
  is.integer0 <- function(x)
  {
    is.integer(x) && length(x) == 0L
  }
  
  R <- rep(NA,length(t.grid))
  pop.hazard.old <- c(NA)
  excess.hazard.old <- c(NA)
  index.pop.old <- c(NA)
  n.events <- numeric(length(t.grid))
  
  new.pop.data.male.inside.func <- new.pop.data.male[new.pop.data.male[,"Year"]>=floor(min(start.year+arrival.vec)),] # remove the part of life table which is not relevant
  new.pop.data.female.inside.func <- new.pop.data.female[new.pop.data.female[,"Year"]>=floor(min(start.year+arrival.vec)),]
  
  for(i in 1:length(t.grid)){
    index <- which(t.grid[i] >= arrival.vec) # Only consider observations who have arrived at time t
    if(is.integer0(index)==T){ # if no one has arrived
      r.t<- 0
      R[i] <- r.t
    } else {
      new.age.vec <- age.vec[index] # age of the patients who have arrived up to this time point
      new.gender.vec <- gender.vec[index]
      if(length(index)==1){ # covariate matrix corresponding to these patients
        new.x.matrix <- matrix(x.matrix[index,],nrow=1)
      } else {
        new.x.matrix <- x.matrix[index,]
      }
      new.time.obs.vec <- time.obs.vec[index]
      new.arrival.vec <- arrival.vec[index]
      new.delta.i.vec <- delta.i.vec[index]
      
      at.risk.time.vec <- x_i.vec(new.time.obs.vec,new.arrival.vec,t.grid[i]) # at risk time of the patients who have arrived at the specific time point
      
      ### POP HAZARD FUNCTION ###
      if(i==1){ # consider now patients who arrive between two consecutive time points on the grid and actually experience an event when the monitoring ends as only these contribute to the first term 
        index.pop <- which(t.grid[i] >= arrival.vec & delta.i.vec==1)
      } else {
        index.pop <- which(t.grid[i-1] < arrival.vec & arrival.vec <= t.grid[i] & delta.i.vec==1) 
      }
      if(is.integer0(index.pop)==F){
        pop.hazard.vec <- rep(NA,length(index.pop)) # a vector to store population hazards evaluated at the observed survival times for these observations
        new.age.vec.pop <- age.vec[index.pop] # age of these patients
        new.gender.vec.pop <- gender.vec[index.pop]
        new.time.obs.vec.pop <- time.obs.vec[index.pop]
        new.arrival.vec.pop <- arrival.vec[index.pop]
        new.age.floor.vec.pop <- floor(new.age.vec.pop+new.time.obs.vec.pop) # age at observed survival times 
        new.age.floor.vec.pop[new.age.floor.vec.pop>110] <- 110 # above 110 years old --> set to 110 
        new.year.floor.vec.pop <- floor(start.year+new.arrival.vec.pop+new.time.obs.vec.pop) # year at observed survival times 
        new.year.floor.vec.pop[new.year.floor.vec.pop>=(end.year.table+1)] <- end.year.table # set to the latest year that appear in the life tables
        if(length(index.pop)!=1){ # covariate values of these patients who arrive in this time interval 
          new.x.matrix.pop <- x.matrix[index.pop,]
        } else {
          new.x.matrix.pop <- matrix(x.matrix[index.pop,],nrow=1)
        }
        for(j in 1:length(index.pop)){ # among the patients who arrive in this time interval --> take out the population hazards and calculate the excess hazards at the observed survival times
          if(new.gender.vec.pop[j]==0){
            
            index.ind <- which(new.pop.data.male.inside.func[,"Age"]==new.age.floor.vec.pop[j] & new.pop.data.male.inside.func[,"Year"]==new.year.floor.vec.pop[j])
            
            pop.hazard.vec[j] <- new.pop.data.male.inside.func[index.ind,"mx"]
          } else {
            index.ind <- which(new.pop.data.female.inside.func[,"Age"]==new.age.floor.vec.pop[j] & new.pop.data.female.inside.func[,"Year"]==new.year.floor.vec.pop[j])
            
            pop.hazard.vec[j] <- new.pop.data.female.inside.func[index.ind,"mx"]
          } 
        }
        pop.hazard.old <- na.omit(c(pop.hazard.old,pop.hazard.vec)) # store it in a larger vector with the population hazards of the earlier arrived patients
        excess.hazard.vec <- excess.hazard.piecewise.func(partition.t.vec,chi.vec,new.x.matrix.pop,beta.vec,new.time.obs.vec.pop)
        excess.hazard.old <- na.omit(c(excess.hazard.old,excess.hazard.vec))
        index.pop.old <- na.omit(c(index.pop.old,index.pop))
      }
      ### ###
      
      h1.vec <- (pop.hazard.old+rho*excess.hazard.old)
      
      h0.vec <- (pop.hazard.old+excess.hazard.old)
      
      log.haz.ratio.vec <- log(h1.vec/h0.vec)
      
      diff.cumulative.haz.vec <- (rho-1)*cumulative.excess.hazard.piecewise.func(partition.t.vec,chi.vec,new.x.matrix,beta.vec,at.risk.time.vec)
      event.vec <- delta_i.t(time.obs.vec[index.pop.old],arrival.vec[index.pop.old],delta.i.vec[index.pop.old],t.grid[i]) # event indicator of patients that will experience an event in the end of monitoring period at this specific time point
      
      if(is.na(index.pop.old[1])==F){
        n.events[i] <- sum(event.vec)
        r.t <- sum(event.vec*log.haz.ratio.vec)-sum(diff.cumulative.haz.vec)
        R[i] <- r.t
      } else {
        n.events[i] <- 0
        r.t <- -sum(diff.cumulative.haz.vec)
        R[i] <- r.t
      }
    }
  }
  cbind(R-cummin(R),n.events)
}

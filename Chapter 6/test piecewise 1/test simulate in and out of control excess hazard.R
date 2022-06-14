vec.sim.excess.piecewise.in.or.out <- function(partition.t.vec,chi.vec,x.matrix,beta.vec,rho,arrival.time,eta){
  interval.length.partition.t.vec <- diff(partition.t.vec) # length of the bands except for the last one
  num.interval.length.partition.t.vec <- length(interval.length.partition.t.vec) # number of bands 
  baseline.vec <- exp(chi.vec) # baseline hazard vector - first value corresponds to the hazard during the first band etc. 
  cumulative.baseline.vec <- cumsum(baseline.vec[1:num.interval.length.partition.t.vec]*diff(partition.t.vec)) # cumulative hazard at a time defined by the bands
  new.cumulative.baseline.vec <- c(0,cumulative.baseline.vec) # the zero element corresponds to cumulative hazard at t=0
  in.cumulative.baseline.excess.hazard.func <- function(s){ # in-control 
    index <- findInterval(s,t.vec) # find the band which a time s is located
    new.cumulative.baseline.vec[index]+(s-t.vec[index])*exp(chi.vec)[index] # cumulative hazard at this specific time s
  }
  out.cumulative.baseline.excess.hazard.func <- function(s){ # out-of-control - otherwise, similar as above 
    index <- findInterval(s,t.vec)
    (new.cumulative.baseline.vec[index]+(s-t.vec[index])*exp(chi.vec)[index])*rho
  }
  
  time.excess.sim <- rep(0,nrow(x.matrix))
  for(i in 1:length(time.excess.sim)){
    final.cumulative.baseline.hazard.func <- function(t){
      if(eta>=arrival.time[i]){ # if a patient arrives before the shift at time eta 
        in.cumulative.baseline.excess.hazard.func(t)*(t<= (eta-arrival.time[i]))+(in.cumulative.baseline.excess.hazard.func(eta-arrival.time[i])+out.cumulative.baseline.excess.hazard.func(t)-out.cumulative.baseline.excess.hazard.func(eta-arrival.time[i]))*(t> (eta-arrival.time[i]))
      } else { # if a patient arrives after eta --> only experience the out-of-control state 
        out.cumulative.baseline.excess.hazard.func(t)
      }
    }
    final.cumulative.excess.hazard.func <- function(t){ # cumulative excess hazard function of this patient
      as.numeric(exp(beta.vec%*%x.matrix[i,]))*final.cumulative.baseline.hazard.func(t)
    }
    u<-runif(1)
    eq.solve <- function(t){
      1-exp(-final.cumulative.excess.hazard.func(t))-u # inverse transform
    }
    time.sim <- uniroot(eq.solve,lower=0,upper=t.vec[length(t.vec)]+10000000)$root
    time.excess.sim[i] <- time.sim
  }
  return(time.excess.sim)
}


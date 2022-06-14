vec.sim.excess.piecewise.in.or.out <- function(partition.t.vec,chi.vec,x.matrix,beta.vec,rho,arrival.time,eta){
  interval.length.partition.t.vec <- diff(partition.t.vec)
  num.interval.length.partition.t.vec <- length(interval.length.partition.t.vec)
  baseline.vec <- exp(chi.vec)
  cumulative.baseline.vec <- cumsum(baseline.vec[1:num.interval.length.partition.t.vec]*diff(partition.t.vec))
  new.cumulative.baseline.vec <- c(0,cumulative.baseline.vec)
  in.cumulative.baseline.excess.hazard.func <- function(s){
    index <- findInterval(s,t.vec)
    new.cumulative.baseline.vec[index]+(s-t.vec[index])*exp(chi.vec)[index]
  }
  out.cumulative.baseline.excess.hazard.func <- function(s){
    index <- findInterval(s,t.vec)
    (new.cumulative.baseline.vec[index]+(s-t.vec[index])*exp(chi.vec)[index])*rho
  }
  
  time.excess.sim <- rep(0,nrow(x.matrix))
  for(i in 1:length(time.excess.sim)){
    final.cumulative.baseline.hazard.func <- function(t){
      if(eta>=arrival.time[i]){
        in.cumulative.baseline.excess.hazard.func(t)*(t<= (eta-arrival.time[i]))+(in.cumulative.baseline.excess.hazard.func(eta-arrival.time[i])+out.cumulative.baseline.excess.hazard.func(t)-out.cumulative.baseline.excess.hazard.func(eta-arrival.time[i]))*(t> (eta-arrival.time[i]))
      } else {
        out.cumulative.baseline.excess.hazard.func(t)
      }
    }
    final.cumulative.excess.hazard.func <- function(t){
      as.numeric(exp(beta.vec%*%x.matrix[i,]))*final.cumulative.baseline.hazard.func(t)
    }
    u<-runif(1)
    eq.solve <- function(t){
      1-exp(-final.cumulative.excess.hazard.func(t))-u
    }
    time.sim <- uniroot(eq.solve,lower=0,upper=t.vec[length(t.vec)]+10000000)$root
    time.excess.sim[i] <- time.sim
  }
  return(time.excess.sim)
}


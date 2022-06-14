set.seed(43)
for(i in 1:4){ # number 4 is the example from Figure 6.1
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
  
  ### First 5 year (from 2010-2015) --> still the usual excess hazard
  
  time.excess.sim.test.1.25 <- vec.sim.excess.piecewise.in.or.out(t.vec,chi.vec,x.matrix.test,beta.vec,1.25,arrival.sim.vec.test,5)

  time.censoring.sim.test <- cbind(rexp(n.sim.test,rate=0.001),end.year-start.year-arrival.sim.vec.test)
  time.censoring.sim.test <- apply(time.censoring.sim.test,1,FUN = min)
  
  
  ## Pick out minimum of the three different survival times as the observed time. 
  ## Also, need to make a death indicator (delta.i) for all causes. 
  time.observed.sim.test.1.25 <- pmin(time.pop.sim.test,time.excess.sim.test.1.25,time.censoring.sim.test)
  delta.i.test.1.25 <- pmax(delta.p.test,as.numeric(time.excess.sim.test.1.25<time.pop.sim.test))*as.numeric(time.censoring.sim.test>time.observed.sim.test.1.25)

  # eta=5
  t.grid <- seq(from=0,to=10,by=0.01) # time grid to evaluate Psi(t)
  R.out.control.fast.1.25 <- cusum_r.t.piecewise(start.year,age.sim.vec.test,gender.sim.vec.test,x.matrix.test,time.observed.sim.test.1.25,arrival.sim.vec.test,delta.i.test.1.25,chi.vec,t.vec,beta.vec,1.25,t.grid,pop.data.male,pop.data.female,2020) # calculate the CUSUM chart
  max.test.1.25.individual <- max(R.out.control.fast.1.25[,1]) # first column to Psi(t)
  
  plot(t.grid,R.out.control.fast.1.25[,1],type="l",main=expression(paste("A signal from a CUSUM chart for a piecewise constant baseline when ",eta,"=5 and ",rho,"=1.25")),
       xlab="t (in years)",ylab=expression(paste(Psi,"(t)")))
  abline(h=threshold.c.1.25,col="red",lty=2)
  print(i)
}

rm(list=ls()) 

library(dplyr)
library(relsurv)
library(popEpi)
library(Epi)
source("Funksjoner.r")
# Example Nr.1
set.seed(42)
lambda <- 1
phi <- 0.05
beta.vec <- c(0.05,0.25,-0.25)
start.year <- 2000
end.year <- 2021

age.sim.vec <- rnorm(10000,mean=70,sd=10) # simulate age from N(mean=70,sigma=10). 
age.sim.vec <- round(age.sim.vec)
n.sim <- length(age.sim.vec)

start.year.sim.vec <- runif(10000,min=0,max=11)
start.year.sim.vec <- floor(start.year.sim.vec)
## Simulate gender

gender.sim.vec <- sample(0:1,size = n.sim,replace=T) # simulate gender (0=male, 1=female)

## Simulate population survival times and population censoring indicator 

matrix.time.pop.sim <- vec.pop.sim(age.sim.vec,rep(start.year,n.sim)+start.year.sim.vec,gender.sim.vec)
time.pop.sim <- matrix.time.pop.sim[,"t.p"]
delta.p <- matrix.time.pop.sim[,"delta.p"]

## Simulate excess survival times 

x.matrix <- cbind(age.sim.vec,gender.sim.vec,start.year.sim.vec)
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

final.data.set <- cbind(x.matrix,rep(start.year,n.sim)+start.year.sim.vec,time.observed.sim,delta.i)
colnames(final.data.set) <- c("age","gender","start.year.after.2000","start.year","time.obs","delta.i")
final.data.set <- as.data.frame(final.data.set)

sum(delta.i==1)
length(which(delta.i==1 & time.excess.sim<=time.pop.sim))
length(which(delta.i==1 & time.excess.sim>=time.pop.sim))
length(which(delta.i==1 & time.excess.sim<=time.pop.sim))/sum(delta.i==1) # proportion of excess events 

# Calculate net survival 
t.grid <- seq(from=0,to=end.year-start.year,by=0.1)

net.survival.individual <- weibull.prop.excess.surv.test(lambda,phi,beta.vec,x.matrix,t.grid) # each row corresponds to an individual, i.e S_{Ei}. The columns correspond to the time steps defined by the time grid.
net.survival.group <- colMeans(net.survival.individual) # take the column mean to get the overall (mariginal) net survival of a group at the given value of time. 
plot(t.grid,net.survival.group,type="l",xlab="Time (in years)",ylab="Survival",main="Net vs Observable Net",lwd=2,ylim=c(0,1.2))

# Calculate observable net survival 

# observable.net.surv.group <- weibull.obs.net.fast.ver(lambda,phi,beta.vec,x.matrix,start.year,t.grid) 
# lines(t.grid,observable.net.surv.group,col="red",lwd=2)
# legend("topright",legend = c("Net Survival","Observable Net Survival"),
#        col=c("black","red"),lty=1,lwd=2)
# 
# # Ederer 2 vs PP 
# 
# nortab <- transrate.hmd(male="NOR.mltper_1x1.txt",female="NOR.fltper_1x1.txt")
# final.data.set$gender <- final.data.set$gender+1 
# ederer2.estimate <- rs.surv(Surv(time.obs*365.241,delta.i)~1,data=final.data.set,ratetable = nortab,method="ederer2",rmap = list(age=age*365.241,sex=gender,year=mdy.date(1,1,start.year)))
# lines(ederer2.estimate$time/365.241,ederer2.estimate$surv,col="orange",lty=2,lwd=2)
# pp.estimate <- rs.surv(Surv(time.obs*365.241,delta.i)~1,data=final.data.set,ratetable = nortab,method="pohar-perme",rmap = list(age=age*365.241,sex=gender,year=mdy.date(1,1,start.year)))
# lines(pp.estimate$time/365.241,pp.estimate$surv,col="green",lty=2,lwd=2)

final.data.set.survtab <- cbind(x.matrix,rep(start.year,n.sim)+start.year.sim.vec,rep(start.year,n.sim)+start.year.sim.vec+time.observed.sim,delta.i)
colnames(final.data.set.survtab) <- c("AGE","sex","start.year.after.2000","start.year","exit","delta.i")
final.data.set.survtab <- as.data.frame(final.data.set.survtab)

# relevant population table for survtab 
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
lines(st.e2$Tstart,st.e2$r.e2,col="orange",lwd=2,lty=2)

st.pp <- survtab(
  Surv(time = FUT, event = lex.Xst) ~ 1,
  data = x.sim,
  surv.type = "surv.rel",
  relsurv.method = "pp",
  breaks = list(FUT = seq(0, 21, 1/12)),
  pophaz = poptable.final
)
lines(st.pp$Tstart,st.pp$r.pp,col="green",lwd=2,lty=2)

for(i in 2:5){ # to get Figure 5.9a
  age.sim.vec <- rnorm(10000,mean=70,sd=10) # simulate age from N(mean=70,sigma=10). 
  age.sim.vec <- round(age.sim.vec)
  n.sim <- length(age.sim.vec)
  
  start.year.sim.vec <- runif(10000,min=0,max=11)
  start.year.sim.vec <- floor(start.year.sim.vec)
  ## Simulate gender
  
  gender.sim.vec <- sample(0:1,size = n.sim,replace=T) # simulate gender (0=male, 1=female)
  
  ## Simulate population survival times and population censoring indicator 
  
  matrix.time.pop.sim <- vec.pop.sim(age.sim.vec,rep(start.year,n.sim)+start.year.sim.vec,gender.sim.vec)
  time.pop.sim <- matrix.time.pop.sim[,"t.p"]
  delta.p <- matrix.time.pop.sim[,"delta.p"]
  
  ## Simulate excess survival times 
  
  x.matrix <- cbind(age.sim.vec,gender.sim.vec,start.year.sim.vec)
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
  
  final.data.set <- cbind(x.matrix,rep(start.year,n.sim)+start.year.sim.vec,time.observed.sim,delta.i)
  colnames(final.data.set) <- c("age","gender","start.year.after.2000","start.year","time.obs","delta.i")
  final.data.set <- as.data.frame(final.data.set)
  
  sum(delta.i==1)
  length(which(delta.i==1 & time.excess.sim<=time.pop.sim))
  length(which(delta.i==1 & time.excess.sim>=time.pop.sim))
  length(which(delta.i==1 & time.excess.sim<=time.pop.sim))/sum(delta.i==1)
  
 
  final.data.set.survtab <- cbind(x.matrix,rep(start.year,n.sim)+start.year.sim.vec,rep(start.year,n.sim)+start.year.sim.vec+time.observed.sim,delta.i)
  colnames(final.data.set.survtab) <- c("AGE","sex","start.year.after.2000","start.year","exit","delta.i")
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
  lines(st.e2$Tstart,st.e2$r.e2,col="orange",lwd=2,lty=2)
  
  st.pp <- survtab(
    Surv(time = FUT, event = lex.Xst) ~ 1,
    data = x.sim,
    surv.type = "surv.rel",
    relsurv.method = "pp",
    breaks = list(FUT = seq(0, 21, 1/12)),
    pophaz = poptable.final
  )
  lines(st.pp$Tstart,st.pp$r.pp,col="green",lwd=2,lty=2)
}
legend("topright",legend = c("Net Survival","Ederer 2","PP"),
       col=c("black","orange","green"),lty=c(1,2,2),lwd=2)

# Example Nr.2
set.seed(42)
lambda <- 1
phi <- 0.05
beta.vec <- c(0.05,0.25,-0.5)
start.year <- 2000
end.year <- 2021

age.sim.vec <- rnorm(10000,mean=70,sd=10) # simulate age from N(mean=70,sigma=10). 
age.sim.vec <- round(age.sim.vec)
n.sim <- length(age.sim.vec)

start.year.sim.vec <- runif(10000,min=0,max=11)
start.year.sim.vec <- floor(start.year.sim.vec)
## Simulate gender

gender.sim.vec <- sample(0:1,size = n.sim,replace=T) # simulate gender (0=male, 1=female)

## Simulate population survival times and population censoring indicator 

matrix.time.pop.sim <- vec.pop.sim(age.sim.vec,rep(start.year,n.sim)+start.year.sim.vec,gender.sim.vec)
time.pop.sim <- matrix.time.pop.sim[,"t.p"]
delta.p <- matrix.time.pop.sim[,"delta.p"]

## Simulate excess survival times 

x.matrix <- cbind(age.sim.vec,gender.sim.vec,start.year.sim.vec)
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

final.data.set <- cbind(x.matrix,rep(start.year,n.sim)+start.year.sim.vec,time.observed.sim,delta.i)
colnames(final.data.set) <- c("age","gender","start.year.after.2000","start.year","time.obs","delta.i")
final.data.set <- as.data.frame(final.data.set)

sum(delta.i==1)
length(which(delta.i==1 & time.excess.sim<=time.pop.sim))
length(which(delta.i==1 & time.excess.sim>=time.pop.sim))
length(which(delta.i==1 & time.excess.sim<=time.pop.sim))/sum(delta.i==1) # proportion of excess events 

# Calculate net survival 
t.grid <- seq(from=0,to=end.year-start.year,by=0.1)

net.survival.individual <- weibull.prop.excess.surv.test(lambda,phi,beta.vec,x.matrix,t.grid) # each row corresponds to an individual, i.e S_{Ei}. The columns correspond to the time steps defined by the time grid.
net.survival.group <- colMeans(net.survival.individual) # take the column mean to get the overall (mariginal) net survival of a group at the given value of time. 
plot(t.grid,net.survival.group,type="l",xlab="Time (in years)",ylab="Survival",main="Net vs Observable Net",lwd=2,ylim=c(0,1.2))

# Calculate observable net survival 

# observable.net.surv.group <- weibull.obs.net.fast.ver(lambda,phi,beta.vec,x.matrix,start.year,t.grid) 
# lines(t.grid,observable.net.surv.group,col="red",lwd=2)
# legend("topright",legend = c("Net Survival","Observable Net Survival"),
#        col=c("black","red"),lty=1,lwd=2)
# 
# # Ederer 2 vs PP 
# 
# nortab <- transrate.hmd(male="NOR.mltper_1x1.txt",female="NOR.fltper_1x1.txt")
# final.data.set$gender <- final.data.set$gender+1 
# ederer2.estimate <- rs.surv(Surv(time.obs*365.241,delta.i)~1,data=final.data.set,ratetable = nortab,method="ederer2",rmap = list(age=age*365.241,sex=gender,year=mdy.date(1,1,start.year)))
# lines(ederer2.estimate$time/365.241,ederer2.estimate$surv,col="orange",lty=2,lwd=2)
# pp.estimate <- rs.surv(Surv(time.obs*365.241,delta.i)~1,data=final.data.set,ratetable = nortab,method="pohar-perme",rmap = list(age=age*365.241,sex=gender,year=mdy.date(1,1,start.year)))
# lines(pp.estimate$time/365.241,pp.estimate$surv,col="green",lty=2,lwd=2)

final.data.set.survtab <- cbind(x.matrix,rep(start.year,n.sim)+start.year.sim.vec,rep(start.year,n.sim)+start.year.sim.vec+time.observed.sim,delta.i)
colnames(final.data.set.survtab) <- c("AGE","sex","start.year.after.2000","start.year","exit","delta.i")
final.data.set.survtab <- as.data.frame(final.data.set.survtab)

# relevant population table for survtab 
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
lines(st.e2$Tstart,st.e2$r.e2,col="orange",lwd=2,lty=2)

st.pp <- survtab(
  Surv(time = FUT, event = lex.Xst) ~ 1,
  data = x.sim,
  surv.type = "surv.rel",
  relsurv.method = "pp",
  breaks = list(FUT = seq(0, 21, 1/12)),
  pophaz = poptable.final
)
lines(st.pp$Tstart,st.pp$r.pp,col="green",lwd=2,lty=2)

for(i in 2:5){ # to get Figure 5.9b
  age.sim.vec <- rnorm(10000,mean=70,sd=10) # simulate age from N(mean=70,sigma=10). 
  age.sim.vec <- round(age.sim.vec)
  n.sim <- length(age.sim.vec)
  
  start.year.sim.vec <- runif(10000,min=0,max=11)
  start.year.sim.vec <- floor(start.year.sim.vec)
  ## Simulate gender
  
  gender.sim.vec <- sample(0:1,size = n.sim,replace=T) # simulate gender (0=male, 1=female)
  
  ## Simulate population survival times and population censoring indicator 
  
  matrix.time.pop.sim <- vec.pop.sim(age.sim.vec,rep(start.year,n.sim)+start.year.sim.vec,gender.sim.vec)
  time.pop.sim <- matrix.time.pop.sim[,"t.p"]
  delta.p <- matrix.time.pop.sim[,"delta.p"]
  
  ## Simulate excess survival times 
  
  x.matrix <- cbind(age.sim.vec,gender.sim.vec,start.year.sim.vec)
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
  
  final.data.set <- cbind(x.matrix,rep(start.year,n.sim)+start.year.sim.vec,time.observed.sim,delta.i)
  colnames(final.data.set) <- c("age","gender","start.year.after.2000","start.year","time.obs","delta.i")
  final.data.set <- as.data.frame(final.data.set)
  
  sum(delta.i==1)
  length(which(delta.i==1 & time.excess.sim<=time.pop.sim))
  length(which(delta.i==1 & time.excess.sim>=time.pop.sim))
  length(which(delta.i==1 & time.excess.sim<=time.pop.sim))/sum(delta.i==1)
  
  
  final.data.set.survtab <- cbind(x.matrix,rep(start.year,n.sim)+start.year.sim.vec,rep(start.year,n.sim)+start.year.sim.vec+time.observed.sim,delta.i)
  colnames(final.data.set.survtab) <- c("AGE","sex","start.year.after.2000","start.year","exit","delta.i")
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
  lines(st.e2$Tstart,st.e2$r.e2,col="orange",lwd=2,lty=2)
  
  st.pp <- survtab(
    Surv(time = FUT, event = lex.Xst) ~ 1,
    data = x.sim,
    surv.type = "surv.rel",
    relsurv.method = "pp",
    breaks = list(FUT = seq(0, 21, 1/12)),
    pophaz = poptable.final
  )
  lines(st.pp$Tstart,st.pp$r.pp,col="green",lwd=2,lty=2)
}
legend("topright",legend = c("Net Survival","Ederer 2","PP"),
       col=c("black","orange","green"),lty=c(1,2,2),lwd=2)

# Example Nr.3

set.seed(42)
lambda <- 1
phi <- 0.05
beta.vec <- c(0.05,0.25,-1)
start.year <- 2000
end.year <- 2021

age.sim.vec <- rnorm(10000,mean=70,sd=10) # simulate age from N(mean=70,sigma=10). 
age.sim.vec <- round(age.sim.vec)
n.sim <- length(age.sim.vec)

start.year.sim.vec <- runif(10000,min=0,max=11)
start.year.sim.vec <- floor(start.year.sim.vec)
## Simulate gender

gender.sim.vec <- sample(0:1,size = n.sim,replace=T) # simulate gender (0=male, 1=female)

## Simulate population survival times and population censoring indicator 

matrix.time.pop.sim <- vec.pop.sim(age.sim.vec,rep(start.year,n.sim)+start.year.sim.vec,gender.sim.vec)
time.pop.sim <- matrix.time.pop.sim[,"t.p"]
delta.p <- matrix.time.pop.sim[,"delta.p"]

## Simulate excess survival times 

x.matrix <- cbind(age.sim.vec,gender.sim.vec,start.year.sim.vec)
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

final.data.set <- cbind(x.matrix,rep(start.year,n.sim)+start.year.sim.vec,time.observed.sim,delta.i)
colnames(final.data.set) <- c("age","gender","start.year.after.2000","start.year","time.obs","delta.i")
final.data.set <- as.data.frame(final.data.set)

sum(delta.i==1)
length(which(delta.i==1 & time.excess.sim<=time.pop.sim))
length(which(delta.i==1 & time.excess.sim>=time.pop.sim))
length(which(delta.i==1 & time.excess.sim<=time.pop.sim))/sum(delta.i==1)

# Calculate net survival 
t.grid <- seq(from=0,to=end.year-start.year,by=0.1)

net.survival.individual <- weibull.prop.excess.surv.test(lambda,phi,beta.vec,x.matrix,t.grid) # each row corresponds to an individual, i.e S_{Ei}. The columns correspond to the time steps defined by the time grid.
net.survival.group <- colMeans(net.survival.individual) # take the column mean to get the overall (mariginal) net survival of a group at the given value of time. 
plot(t.grid,net.survival.group,type="l",xlab="Time (in years)",ylab="Survival",main="Net vs Observable Net",lwd=2,ylim=c(0,1.2))

# Calculate observable net survival 

# observable.net.surv.group <- weibull.obs.net.fast.ver(lambda,phi,beta.vec,x.matrix,start.year,t.grid) 
# lines(t.grid,observable.net.surv.group,col="red",lwd=2)
# legend("topright",legend = c("Net Survival","Observable Net Survival"),
#        col=c("black","red"),lty=1,lwd=2)
# 
# # Ederer 2 vs PP 
# 
# nortab <- transrate.hmd(male="NOR.mltper_1x1.txt",female="NOR.fltper_1x1.txt")
# final.data.set$gender <- final.data.set$gender+1 
# ederer2.estimate <- rs.surv(Surv(time.obs*365.241,delta.i)~1,data=final.data.set,ratetable = nortab,method="ederer2",rmap = list(age=age*365.241,sex=gender,year=mdy.date(1,1,start.year)))
# lines(ederer2.estimate$time/365.241,ederer2.estimate$surv,col="orange",lty=2,lwd=2)
# pp.estimate <- rs.surv(Surv(time.obs*365.241,delta.i)~1,data=final.data.set,ratetable = nortab,method="pohar-perme",rmap = list(age=age*365.241,sex=gender,year=mdy.date(1,1,start.year)))
# lines(pp.estimate$time/365.241,pp.estimate$surv,col="green",lty=2,lwd=2)

final.data.set.survtab <- cbind(x.matrix,rep(start.year,n.sim)+start.year.sim.vec,rep(start.year,n.sim)+start.year.sim.vec+time.observed.sim,delta.i)
colnames(final.data.set.survtab) <- c("AGE","sex","start.year.after.2000","start.year","exit","delta.i")
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
lines(st.e2$Tstart,st.e2$r.e2,col="orange",lwd=2,lty=2)

st.pp <- survtab(
  Surv(time = FUT, event = lex.Xst) ~ 1,
  data = x.sim,
  surv.type = "surv.rel",
  relsurv.method = "pp",
  breaks = list(FUT = seq(0, 21, 1/12)),
  pophaz = poptable.final
)
lines(st.pp$Tstart,st.pp$r.pp,col="green",lwd=2,lty=2)

for(i in 2:5){ # to get Figure 5.10a
  age.sim.vec <- rnorm(10000,mean=70,sd=10) # simulate age from N(mean=70,sigma=10). 
  age.sim.vec <- round(age.sim.vec)
  n.sim <- length(age.sim.vec)
  
  start.year.sim.vec <- runif(10000,min=0,max=11)
  start.year.sim.vec <- floor(start.year.sim.vec)
  ## Simulate gender
  
  gender.sim.vec <- sample(0:1,size = n.sim,replace=T) # simulate gender (0=male, 1=female)
  
  ## Simulate population survival times and population censoring indicator 
  
  matrix.time.pop.sim <- vec.pop.sim(age.sim.vec,rep(start.year,n.sim)+start.year.sim.vec,gender.sim.vec)
  time.pop.sim <- matrix.time.pop.sim[,"t.p"]
  delta.p <- matrix.time.pop.sim[,"delta.p"]
  
  ## Simulate excess survival times 
  
  x.matrix <- cbind(age.sim.vec,gender.sim.vec,start.year.sim.vec)
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
  
  final.data.set <- cbind(x.matrix,rep(start.year,n.sim)+start.year.sim.vec,time.observed.sim,delta.i)
  colnames(final.data.set) <- c("age","gender","start.year.after.2000","start.year","time.obs","delta.i")
  final.data.set <- as.data.frame(final.data.set)
  
  sum(delta.i==1)
  length(which(delta.i==1 & time.excess.sim<=time.pop.sim))
  length(which(delta.i==1 & time.excess.sim>=time.pop.sim))
  length(which(delta.i==1 & time.excess.sim<=time.pop.sim))/sum(delta.i==1)
  
  
  final.data.set.survtab <- cbind(x.matrix,rep(start.year,n.sim)+start.year.sim.vec,rep(start.year,n.sim)+start.year.sim.vec+time.observed.sim,delta.i)
  colnames(final.data.set.survtab) <- c("AGE","sex","start.year.after.2000","start.year","exit","delta.i")
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
  lines(st.e2$Tstart,st.e2$r.e2,col="orange",lwd=2,lty=2)
  
  st.pp <- survtab(
    Surv(time = FUT, event = lex.Xst) ~ 1,
    data = x.sim,
    surv.type = "surv.rel",
    relsurv.method = "pp",
    breaks = list(FUT = seq(0, 21, 1/12)),
    pophaz = poptable.final
  )
  lines(st.pp$Tstart,st.pp$r.pp,col="green",lwd=2,lty=2)
}
legend("topright",legend = c("Net Survival","Ederer 2","PP"),
       col=c("black","orange","green"),lty=c(1,2,2),lwd=2)

# Example Nr.4

set.seed(42)
lambda <- 1
phi <- 0.005
beta.vec <- c(0.05,0.25,-0.25)
start.year <- 2000
end.year <- 2021

age.sim.vec <- rnorm(10000,mean=70,sd=10) # simulate age from N(mean=70,sigma=10). 
age.sim.vec <- round(age.sim.vec)
n.sim <- length(age.sim.vec)

start.year.sim.vec <- runif(10000,min=0,max=11)
start.year.sim.vec <- floor(start.year.sim.vec)
## Simulate gender

gender.sim.vec <- sample(0:1,size = n.sim,replace=T) # simulate gender (0=male, 1=female)

## Simulate population survival times and population censoring indicator 

matrix.time.pop.sim <- vec.pop.sim(age.sim.vec,rep(start.year,n.sim)+start.year.sim.vec,gender.sim.vec)
time.pop.sim <- matrix.time.pop.sim[,"t.p"]
delta.p <- matrix.time.pop.sim[,"delta.p"]

## Simulate excess survival times 

x.matrix <- cbind(age.sim.vec,gender.sim.vec,start.year.sim.vec)
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

final.data.set <- cbind(x.matrix,rep(start.year,n.sim)+start.year.sim.vec,time.observed.sim,delta.i)
colnames(final.data.set) <- c("age","gender","start.year.after.2000","start.year","time.obs","delta.i")
final.data.set <- as.data.frame(final.data.set)

sum(delta.i==1)
length(which(delta.i==1 & time.excess.sim<=time.pop.sim))
length(which(delta.i==1 & time.excess.sim>=time.pop.sim))
length(which(delta.i==1 & time.excess.sim<=time.pop.sim))/sum(delta.i==1)

# Calculate net survival 
t.grid <- seq(from=0,to=end.year-start.year,by=0.1)

net.survival.individual <- weibull.prop.excess.surv.test(lambda,phi,beta.vec,x.matrix,t.grid) # each row corresponds to an individual, i.e S_{Ei}. The columns correspond to the time steps defined by the time grid.
net.survival.group <- colMeans(net.survival.individual) # take the column mean to get the overall (mariginal) net survival of a group at the given value of time. 
plot(t.grid,net.survival.group,type="l",xlab="Time (in years)",ylab="Survival",main="Net vs Observable Net",lwd=2,ylim=c(0,1.2))

# Calculate observable net survival 

# observable.net.surv.group <- weibull.obs.net.fast.ver(lambda,phi,beta.vec,x.matrix,start.year,t.grid) 
# lines(t.grid,observable.net.surv.group,col="red",lwd=2)
# legend("topright",legend = c("Net Survival","Observable Net Survival"),
#        col=c("black","red"),lty=1,lwd=2)
# 
# # Ederer 2 vs PP 
# 
# nortab <- transrate.hmd(male="NOR.mltper_1x1.txt",female="NOR.fltper_1x1.txt")
# final.data.set$gender <- final.data.set$gender+1 
# ederer2.estimate <- rs.surv(Surv(time.obs*365.241,delta.i)~1,data=final.data.set,ratetable = nortab,method="ederer2",rmap = list(age=age*365.241,sex=gender,year=mdy.date(1,1,start.year)))
# lines(ederer2.estimate$time/365.241,ederer2.estimate$surv,col="orange",lty=2,lwd=2)
# pp.estimate <- rs.surv(Surv(time.obs*365.241,delta.i)~1,data=final.data.set,ratetable = nortab,method="pohar-perme",rmap = list(age=age*365.241,sex=gender,year=mdy.date(1,1,start.year)))
# lines(pp.estimate$time/365.241,pp.estimate$surv,col="green",lty=2,lwd=2)

final.data.set.survtab <- cbind(x.matrix,rep(start.year,n.sim)+start.year.sim.vec,rep(start.year,n.sim)+start.year.sim.vec+time.observed.sim,delta.i)
colnames(final.data.set.survtab) <- c("AGE","sex","start.year.after.2000","start.year","exit","delta.i")
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
lines(st.e2$Tstart,st.e2$r.e2,col="orange",lwd=2,lty=2)

st.pp <- survtab(
  Surv(time = FUT, event = lex.Xst) ~ 1,
  data = x.sim,
  surv.type = "surv.rel",
  relsurv.method = "pp",
  breaks = list(FUT = seq(0, 21, 1/12)),
  pophaz = poptable.final
)
lines(st.pp$Tstart,st.pp$r.pp,col="green",lwd=2,lty=2)

for(i in 2:5){ # to get Figure 5.10b
  age.sim.vec <- rnorm(10000,mean=70,sd=10) # simulate age from N(mean=70,sigma=10). 
  age.sim.vec <- round(age.sim.vec)
  n.sim <- length(age.sim.vec)
  
  start.year.sim.vec <- runif(10000,min=0,max=11)
  start.year.sim.vec <- floor(start.year.sim.vec)
  ## Simulate gender
  
  gender.sim.vec <- sample(0:1,size = n.sim,replace=T) # simulate gender (0=male, 1=female)
  
  ## Simulate population survival times and population censoring indicator 
  
  matrix.time.pop.sim <- vec.pop.sim(age.sim.vec,rep(start.year,n.sim)+start.year.sim.vec,gender.sim.vec)
  time.pop.sim <- matrix.time.pop.sim[,"t.p"]
  delta.p <- matrix.time.pop.sim[,"delta.p"]
  
  ## Simulate excess survival times 
  
  x.matrix <- cbind(age.sim.vec,gender.sim.vec,start.year.sim.vec)
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
  
  final.data.set <- cbind(x.matrix,rep(start.year,n.sim)+start.year.sim.vec,time.observed.sim,delta.i)
  colnames(final.data.set) <- c("age","gender","start.year.after.2000","start.year","time.obs","delta.i")
  final.data.set <- as.data.frame(final.data.set)
  
  sum(delta.i==1)
  length(which(delta.i==1 & time.excess.sim<=time.pop.sim))
  length(which(delta.i==1 & time.excess.sim>=time.pop.sim))
  length(which(delta.i==1 & time.excess.sim<=time.pop.sim))/sum(delta.i==1)
  
  
  final.data.set.survtab <- cbind(x.matrix,rep(start.year,n.sim)+start.year.sim.vec,rep(start.year,n.sim)+start.year.sim.vec+time.observed.sim,delta.i)
  colnames(final.data.set.survtab) <- c("AGE","sex","start.year.after.2000","start.year","exit","delta.i")
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
  lines(st.e2$Tstart,st.e2$r.e2,col="orange",lwd=2,lty=2)
  
  st.pp <- survtab(
    Surv(time = FUT, event = lex.Xst) ~ 1,
    data = x.sim,
    surv.type = "surv.rel",
    relsurv.method = "pp",
    breaks = list(FUT = seq(0, 21, 1/12)),
    pophaz = poptable.final
  )
  lines(st.pp$Tstart,st.pp$r.pp,col="green",lwd=2,lty=2)
}
legend("topright",legend = c("Net Survival","Ederer 2","PP"),
       col=c("black","orange","green"),lty=c(1,2,2),lwd=2)

# Example Nr.5

set.seed(42)
lambda <- 1
phi <- 0.0005
beta.vec <- c(0.05,0.25,0.25)
start.year <- 2000
end.year <- 2021

age.sim.vec <- rnorm(10000,mean=70,sd=10) # simulate age from N(mean=70,sigma=10). 
age.sim.vec <- round(age.sim.vec)
n.sim <- length(age.sim.vec)

start.year.sim.vec <- runif(10000,min=0,max=11)
start.year.sim.vec <- floor(start.year.sim.vec)
## Simulate gender

gender.sim.vec <- sample(0:1,size = n.sim,replace=T) # simulate gender (0=male, 1=female)

## Simulate population survival times and population censoring indicator 

matrix.time.pop.sim <- vec.pop.sim(age.sim.vec,rep(start.year,n.sim)+start.year.sim.vec,gender.sim.vec)
time.pop.sim <- matrix.time.pop.sim[,"t.p"]
delta.p <- matrix.time.pop.sim[,"delta.p"]

## Simulate excess survival times 

x.matrix <- cbind(age.sim.vec,gender.sim.vec,start.year.sim.vec)
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

final.data.set <- cbind(x.matrix,rep(start.year,n.sim)+start.year.sim.vec,time.observed.sim,delta.i)
colnames(final.data.set) <- c("age","gender","start.year.after.2000","start.year","time.obs","delta.i")
final.data.set <- as.data.frame(final.data.set)

sum(delta.i==1)
length(which(delta.i==1 & time.excess.sim<=time.pop.sim))
length(which(delta.i==1 & time.excess.sim>=time.pop.sim))
length(which(delta.i==1 & time.excess.sim<=time.pop.sim))/sum(delta.i==1)

# Calculate net survival 
t.grid <- seq(from=0,to=end.year-start.year,by=0.1)

net.survival.individual <- weibull.prop.excess.surv.test(lambda,phi,beta.vec,x.matrix,t.grid) # each row corresponds to an individual, i.e S_{Ei}. The columns correspond to the time steps defined by the time grid.
net.survival.group <- colMeans(net.survival.individual) # take the column mean to get the overall (mariginal) net survival of a group at the given value of time. 
plot(t.grid,net.survival.group,type="l",xlab="Time (in years)",ylab="Survival",main="Net vs Observable Net",lwd=2,ylim=c(0,1.2))

# Calculate observable net survival 

# observable.net.surv.group <- weibull.obs.net.fast.ver(lambda,phi,beta.vec,x.matrix,start.year,t.grid) 
# lines(t.grid,observable.net.surv.group,col="red",lwd=2)
# legend("topright",legend = c("Net Survival","Observable Net Survival"),
#        col=c("black","red"),lty=1,lwd=2)
# 
# # Ederer 2 vs PP 
# 
# nortab <- transrate.hmd(male="NOR.mltper_1x1.txt",female="NOR.fltper_1x1.txt")
# final.data.set$gender <- final.data.set$gender+1 
# ederer2.estimate <- rs.surv(Surv(time.obs*365.241,delta.i)~1,data=final.data.set,ratetable = nortab,method="ederer2",rmap = list(age=age*365.241,sex=gender,year=mdy.date(1,1,start.year)))
# lines(ederer2.estimate$time/365.241,ederer2.estimate$surv,col="orange",lty=2,lwd=2)
# pp.estimate <- rs.surv(Surv(time.obs*365.241,delta.i)~1,data=final.data.set,ratetable = nortab,method="pohar-perme",rmap = list(age=age*365.241,sex=gender,year=mdy.date(1,1,start.year)))
# lines(pp.estimate$time/365.241,pp.estimate$surv,col="green",lty=2,lwd=2)

final.data.set.survtab <- cbind(x.matrix,rep(start.year,n.sim)+start.year.sim.vec,rep(start.year,n.sim)+start.year.sim.vec+time.observed.sim,delta.i)
colnames(final.data.set.survtab) <- c("AGE","sex","start.year.after.2000","start.year","exit","delta.i")
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
lines(st.e2$Tstart,st.e2$r.e2,col="orange",lwd=2,lty=2)

st.pp <- survtab(
  Surv(time = FUT, event = lex.Xst) ~ 1,
  data = x.sim,
  surv.type = "surv.rel",
  relsurv.method = "pp",
  breaks = list(FUT = seq(0, 21, 1/12)),
  pophaz = poptable.final
)
lines(st.pp$Tstart,st.pp$r.pp,col="green",lwd=2,lty=2)

for(i in 2:5){ # to get Figure 5.11a
  age.sim.vec <- rnorm(10000,mean=70,sd=10) # simulate age from N(mean=70,sigma=10). 
  age.sim.vec <- round(age.sim.vec)
  n.sim <- length(age.sim.vec)
  
  start.year.sim.vec <- runif(10000,min=0,max=11)
  start.year.sim.vec <- floor(start.year.sim.vec)
  ## Simulate gender
  
  gender.sim.vec <- sample(0:1,size = n.sim,replace=T) # simulate gender (0=male, 1=female)
  
  ## Simulate population survival times and population censoring indicator 
  
  matrix.time.pop.sim <- vec.pop.sim(age.sim.vec,rep(start.year,n.sim)+start.year.sim.vec,gender.sim.vec)
  time.pop.sim <- matrix.time.pop.sim[,"t.p"]
  delta.p <- matrix.time.pop.sim[,"delta.p"]
  
  ## Simulate excess survival times 
  
  x.matrix <- cbind(age.sim.vec,gender.sim.vec,start.year.sim.vec)
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
  
  final.data.set <- cbind(x.matrix,rep(start.year,n.sim)+start.year.sim.vec,time.observed.sim,delta.i)
  colnames(final.data.set) <- c("age","gender","start.year.after.2000","start.year","time.obs","delta.i")
  final.data.set <- as.data.frame(final.data.set)
  
  sum(delta.i==1)
  length(which(delta.i==1 & time.excess.sim<=time.pop.sim))
  length(which(delta.i==1 & time.excess.sim>=time.pop.sim))
  length(which(delta.i==1 & time.excess.sim<=time.pop.sim))/sum(delta.i==1)
  
  
  final.data.set.survtab <- cbind(x.matrix,rep(start.year,n.sim)+start.year.sim.vec,rep(start.year,n.sim)+start.year.sim.vec+time.observed.sim,delta.i)
  colnames(final.data.set.survtab) <- c("AGE","sex","start.year.after.2000","start.year","exit","delta.i")
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
  lines(st.e2$Tstart,st.e2$r.e2,col="orange",lwd=2,lty=2)
  
  st.pp <- survtab(
    Surv(time = FUT, event = lex.Xst) ~ 1,
    data = x.sim,
    surv.type = "surv.rel",
    relsurv.method = "pp",
    breaks = list(FUT = seq(0, 21, 1/12)),
    pophaz = poptable.final
  )
  lines(st.pp$Tstart,st.pp$r.pp,col="green",lwd=2,lty=2)
}
legend("topright",legend = c("Net Survival","Ederer 2","PP"),
       col=c("red","orange","green"),lty=c(1,2,2),lwd=2)

# Example Nr.6

set.seed(42)
lambda <- 1
phi <- 0.001
beta.vec <- c(0.05,0.25,0)
start.year <- 2000
end.year <- 2021

age.sim.vec <- rnorm(10000,mean=70,sd=10) # simulate age from N(mean=70,sigma=10). 
age.sim.vec <- round(age.sim.vec)
n.sim <- length(age.sim.vec)

start.year.sim.vec <- runif(10000,min=0,max=11)
start.year.sim.vec <- floor(start.year.sim.vec)
## Simulate gender

gender.sim.vec <- sample(0:1,size = n.sim,replace=T) # simulate gender (0=male, 1=female)

## Simulate population survival times and population censoring indicator 

matrix.time.pop.sim <- vec.pop.sim(age.sim.vec,rep(start.year,n.sim)+start.year.sim.vec,gender.sim.vec)
time.pop.sim <- matrix.time.pop.sim[,"t.p"]
delta.p <- matrix.time.pop.sim[,"delta.p"]

## Simulate excess survival times 

x.matrix <- cbind(age.sim.vec,gender.sim.vec,start.year.sim.vec)
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

final.data.set <- cbind(x.matrix,rep(start.year,n.sim)+start.year.sim.vec,time.observed.sim,delta.i)
colnames(final.data.set) <- c("age","gender","start.year.after.2000","start.year","time.obs","delta.i")
final.data.set <- as.data.frame(final.data.set)

sum(delta.i==1)
length(which(delta.i==1 & time.excess.sim<=time.pop.sim))
length(which(delta.i==1 & time.excess.sim>=time.pop.sim))
length(which(delta.i==1 & time.excess.sim<=time.pop.sim))/sum(delta.i==1)

# Calculate net survival 
t.grid <- seq(from=0,to=end.year-start.year,by=0.1)

net.survival.individual <- weibull.prop.excess.surv.test(lambda,phi,beta.vec,x.matrix,t.grid) # each row corresponds to an individual, i.e S_{Ei}. The columns correspond to the time steps defined by the time grid.
net.survival.group <- colMeans(net.survival.individual) # take the column mean to get the overall (mariginal) net survival of a group at the given value of time. 
plot(t.grid,net.survival.group,type="l",xlab="Time (in years)",ylab="Survival",main="Net vs Observable Net",lwd=2,ylim=c(0,1.2))

# Calculate observable net survival 

# observable.net.surv.group <- weibull.obs.net.fast.ver(lambda,phi,beta.vec,x.matrix,start.year,t.grid) 
# lines(t.grid,observable.net.surv.group,col="red",lwd=2)
# legend("topright",legend = c("Net Survival","Observable Net Survival"),
#        col=c("black","red"),lty=1,lwd=2)
# 
# # Ederer 2 vs PP 
# 
# nortab <- transrate.hmd(male="NOR.mltper_1x1.txt",female="NOR.fltper_1x1.txt")
# final.data.set$gender <- final.data.set$gender+1 
# ederer2.estimate <- rs.surv(Surv(time.obs*365.241,delta.i)~1,data=final.data.set,ratetable = nortab,method="ederer2",rmap = list(age=age*365.241,sex=gender,year=mdy.date(1,1,start.year)))
# lines(ederer2.estimate$time/365.241,ederer2.estimate$surv,col="orange",lty=2,lwd=2)
# pp.estimate <- rs.surv(Surv(time.obs*365.241,delta.i)~1,data=final.data.set,ratetable = nortab,method="pohar-perme",rmap = list(age=age*365.241,sex=gender,year=mdy.date(1,1,start.year)))
# lines(pp.estimate$time/365.241,pp.estimate$surv,col="green",lty=2,lwd=2)

final.data.set.survtab <- cbind(x.matrix,rep(start.year,n.sim)+start.year.sim.vec,rep(start.year,n.sim)+start.year.sim.vec+time.observed.sim,delta.i)
colnames(final.data.set.survtab) <- c("AGE","sex","start.year.after.2000","start.year","exit","delta.i")
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
lines(st.e2$Tstart,st.e2$r.e2,col="orange",lwd=2,lty=2)

st.pp <- survtab(
  Surv(time = FUT, event = lex.Xst) ~ 1,
  data = x.sim,
  surv.type = "surv.rel",
  relsurv.method = "pp",
  breaks = list(FUT = seq(0, 21, 1/12)),
  pophaz = poptable.final
)
lines(st.pp$Tstart,st.pp$r.pp,col="green",lwd=2,lty=2)

for(i in 2:5){ # to get Figure 5.11b
  age.sim.vec <- rnorm(10000,mean=70,sd=10) # simulate age from N(mean=70,sigma=10). 
  age.sim.vec <- round(age.sim.vec)
  n.sim <- length(age.sim.vec)
  
  start.year.sim.vec <- runif(10000,min=0,max=11)
  start.year.sim.vec <- floor(start.year.sim.vec)
  ## Simulate gender
  
  gender.sim.vec <- sample(0:1,size = n.sim,replace=T) # simulate gender (0=male, 1=female)
  
  ## Simulate population survival times and population censoring indicator 
  
  matrix.time.pop.sim <- vec.pop.sim(age.sim.vec,rep(start.year,n.sim)+start.year.sim.vec,gender.sim.vec)
  time.pop.sim <- matrix.time.pop.sim[,"t.p"]
  delta.p <- matrix.time.pop.sim[,"delta.p"]
  
  ## Simulate excess survival times 
  
  x.matrix <- cbind(age.sim.vec,gender.sim.vec,start.year.sim.vec)
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
  
  final.data.set <- cbind(x.matrix,rep(start.year,n.sim)+start.year.sim.vec,time.observed.sim,delta.i)
  colnames(final.data.set) <- c("age","gender","start.year.after.2000","start.year","time.obs","delta.i")
  final.data.set <- as.data.frame(final.data.set)
  
  sum(delta.i==1)
  length(which(delta.i==1 & time.excess.sim<=time.pop.sim))
  length(which(delta.i==1 & time.excess.sim>=time.pop.sim))
  length(which(delta.i==1 & time.excess.sim<=time.pop.sim))/sum(delta.i==1)
  
  
  final.data.set.survtab <- cbind(x.matrix,rep(start.year,n.sim)+start.year.sim.vec,rep(start.year,n.sim)+start.year.sim.vec+time.observed.sim,delta.i)
  colnames(final.data.set.survtab) <- c("AGE","sex","start.year.after.2000","start.year","exit","delta.i")
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
  lines(st.e2$Tstart,st.e2$r.e2,col="orange",lwd=2,lty=2)
  
  st.pp <- survtab(
    Surv(time = FUT, event = lex.Xst) ~ 1,
    data = x.sim,
    surv.type = "surv.rel",
    relsurv.method = "pp",
    breaks = list(FUT = seq(0, 21, 1/12)),
    pophaz = poptable.final
  )
  lines(st.pp$Tstart,st.pp$r.pp,col="green",lwd=2,lty=2)
}
legend("topright",legend = c("Net Survival","Ederer 2","PP"),
       col=c("red","orange","green"),lty=c(1,2,2),lwd=2)




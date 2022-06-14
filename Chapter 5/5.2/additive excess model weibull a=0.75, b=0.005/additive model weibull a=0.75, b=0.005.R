library(relsurv)
library(dplyr)
source("Funksjoner.r")

lambda <- 0.75
phi <- 0.005
beta.vec <- c(0.05,0.1,0.5)
start.year <- 2000
end.year <- 2021
n.sim.models <- 250
t.grid <- seq(from=0,to=end.year-start.year,by=0.1)

beta.age.em.default <- rep(0,n.sim.models)
beta.age.em.bwin.0 <- rep(0,n.sim.models)
beta.age.em.bwin.1 <- rep(0,n.sim.models)
beta.age.em.bwin.100 <- rep(0,n.sim.models)
beta.age.poisson.1 <- rep(0,n.sim.models)
beta.age.poisson.2 <- rep(0,n.sim.models)
beta.age.ml.1 <- rep(0,n.sim.models)
beta.age.ml.2 <- rep(0,n.sim.models)

beta.gender.em.default <- rep(0,n.sim.models)
beta.gender.em.bwin.0 <- rep(0,n.sim.models)
beta.gender.em.bwin.1 <- rep(0,n.sim.models)
beta.gender.em.bwin.100 <- rep(0,n.sim.models)
beta.gender.poisson.1 <- rep(0,n.sim.models)
beta.gender.poisson.2 <- rep(0,n.sim.models)
beta.gender.ml.1 <- rep(0,n.sim.models)
beta.gender.ml.2 <- rep(0,n.sim.models)

beta.treatment.em.default <- rep(0,n.sim.models)
beta.treatment.em.bwin.0 <- rep(0,n.sim.models)
beta.treatment.em.bwin.1 <- rep(0,n.sim.models)
beta.treatment.em.bwin.100 <- rep(0,n.sim.models)
beta.treatment.poisson.1 <- rep(0,n.sim.models)
beta.treatment.poisson.2 <- rep(0,n.sim.models)
beta.treatment.ml.1 <- rep(0,n.sim.models)
beta.treatment.ml.2 <- rep(0,n.sim.models)

br.max.test.em.default.age <- rep(0,n.sim.models)
br.max.test.em.bwin.0.age <- rep(0,n.sim.models)
br.max.test.em.bwin.1.age <- rep(0,n.sim.models)
br.max.test.em.bwin.100.age <- rep(0,n.sim.models)
br.max.test.poisson.1.age <- rep(0,n.sim.models)
br.max.test.poisson.2.age <- rep(0,n.sim.models)
br.max.test.ml.1.age <- rep(0,n.sim.models)
br.max.test.ml.2.age <- rep(0,n.sim.models)

br.max.test.em.default.gender <- rep(0,n.sim.models)
br.max.test.em.bwin.0.gender <- rep(0,n.sim.models)
br.max.test.em.bwin.1.gender <- rep(0,n.sim.models)
br.max.test.em.bwin.100.gender <- rep(0,n.sim.models)
br.max.test.poisson.1.gender <- rep(0,n.sim.models)
br.max.test.poisson.2.gender <- rep(0,n.sim.models)
br.max.test.ml.1.gender <- rep(0,n.sim.models)
br.max.test.ml.2.gender <- rep(0,n.sim.models)

br.max.test.em.default.treatment <- rep(0,n.sim.models)
br.max.test.em.bwin.0.treatment <- rep(0,n.sim.models)
br.max.test.em.bwin.1.treatment <- rep(0,n.sim.models)
br.max.test.em.bwin.100.treatment <- rep(0,n.sim.models)
br.max.test.poisson.1.treatment <- rep(0,n.sim.models)
br.max.test.poisson.2.treatment <- rep(0,n.sim.models)
br.max.test.ml.1.treatment <- rep(0,n.sim.models)
br.max.test.ml.2.treatment <- rep(0,n.sim.models)

br.max.weighted.test.em.default.age <- rep(0,n.sim.models)
br.max.weighted.test.em.bwin.0.age <- rep(0,n.sim.models)
br.max.weighted.test.em.bwin.1.age <- rep(0,n.sim.models)
br.max.weighted.test.em.bwin.100.age <- rep(0,n.sim.models)
br.max.weighted.test.poisson.1.age <- rep(0,n.sim.models)
br.max.weighted.test.poisson.2.age <- rep(0,n.sim.models)
br.max.weighted.test.ml.1.age <- rep(0,n.sim.models)
br.max.weighted.test.ml.2.age <- rep(0,n.sim.models)

br.max.weighted.test.em.default.gender <- rep(0,n.sim.models)
br.max.weighted.test.em.bwin.0.gender <- rep(0,n.sim.models)
br.max.weighted.test.em.bwin.1.gender <- rep(0,n.sim.models)
br.max.weighted.test.em.bwin.100.gender <- rep(0,n.sim.models)
br.max.weighted.test.poisson.1.gender <- rep(0,n.sim.models)
br.max.weighted.test.poisson.2.gender <- rep(0,n.sim.models)
br.max.weighted.test.ml.1.gender <- rep(0,n.sim.models)
br.max.weighted.test.ml.2.gender <- rep(0,n.sim.models)

br.max.weighted.test.em.default.treatment <- rep(0,n.sim.models)
br.max.weighted.test.em.bwin.0.treatment <- rep(0,n.sim.models)
br.max.weighted.test.em.bwin.1.treatment <- rep(0,n.sim.models)
br.max.weighted.test.em.bwin.100.treatment <- rep(0,n.sim.models)
br.max.weighted.test.poisson.1.treatment <- rep(0,n.sim.models)
br.max.weighted.test.poisson.2.treatment <- rep(0,n.sim.models)
br.max.weighted.test.ml.1.treatment <- rep(0,n.sim.models)
br.max.weighted.test.ml.2.treatment <- rep(0,n.sim.models)

cvm.test.em.default.age <- rep(0,n.sim.models)
cvm.test.em.bwin.0.age <- rep(0,n.sim.models)
cvm.test.em.bwin.1.age <- rep(0,n.sim.models)
cvm.test.em.bwin.100.age <- rep(0,n.sim.models)
cvm.test.poisson.1.age <- rep(0,n.sim.models)
cvm.test.poisson.2.age <- rep(0,n.sim.models)
cvm.test.ml.1.age <- rep(0,n.sim.models)
cvm.test.ml.2.age <- rep(0,n.sim.models)

cvm.test.em.default.gender <- rep(0,n.sim.models)
cvm.test.em.bwin.0.gender <- rep(0,n.sim.models)
cvm.test.em.bwin.1.gender <- rep(0,n.sim.models)
cvm.test.em.bwin.100.gender <- rep(0,n.sim.models)
cvm.test.poisson.1.gender <- rep(0,n.sim.models)
cvm.test.poisson.2.gender <- rep(0,n.sim.models)
cvm.test.ml.1.gender <- rep(0,n.sim.models)
cvm.test.ml.2.gender <- rep(0,n.sim.models)

cvm.test.em.default.treatment <- rep(0,n.sim.models)
cvm.test.em.bwin.0.treatment <- rep(0,n.sim.models)
cvm.test.em.bwin.1.treatment <- rep(0,n.sim.models)
cvm.test.em.bwin.100.treatment <- rep(0,n.sim.models)
cvm.test.poisson.1.treatment <- rep(0,n.sim.models)
cvm.test.poisson.2.treatment <- rep(0,n.sim.models)
cvm.test.ml.1.treatment <- rep(0,n.sim.models)
cvm.test.ml.2.treatment <- rep(0,n.sim.models)

min.diff.time.obs <- rep(0,n.sim.models)
prop.excess.deaths <- rep(0,n.sim.models)

Lambda0.em.default <- list()
lambda0.em.default <- list()

Lambda0.em.bwin.0 <- list()
lambda0.em.bwin.0 <- list()

Lambda0.em.bwin.1 <- list()
lambda0.em.bwin.1 <- list()

Lambda0.em.bwin.100 <- list()
lambda0.em.bwin.100 <- list()
set.seed(4242)
for(i in 1:n.sim.models){
  age.sim.vec <- rnorm(1000,mean=70,sd=10) # simulate age from N(mean=70,sigma=10). 
  age.sim.vec <- round(age.sim.vec)
  n.sim <- length(age.sim.vec)
  
  ## Simulate gender
  
  gender.sim.vec <- sample(0:1,size = n.sim,replace=T) # simulate gender (0=male, 1=female)
  
  ## Simulate treatment group 
  
  treatment.sim.vec <- sample(0:1,size = n.sim,replace=T) 
  
  ## Simulate population survival times and population censoring indicator 
  
  matrix.time.pop.sim <- vec.pop.sim(age.sim.vec,rep(start.year,n.sim),gender.sim.vec)
  time.pop.sim <- matrix.time.pop.sim[,"t.p"]
  delta.p <- matrix.time.pop.sim[,"delta.p"]
  
  ## Simulate excess survival times - Weibull with baseline shape=1, baseline scale = 0.0025, beta_1=effect of age=0.05, beta_2=effect of gender=0.1 
  
  x.matrix <- cbind(age.sim.vec,gender.sim.vec,treatment.sim.vec)
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
  colnames(final.data.set) <- c("age","gender","treatment","year","time.obs","delta.i")
  final.data.set <- as.data.frame(final.data.set)
  
  nortab <- transrate.hmd(male="NOR.mltper_1x1.txt",female="NOR.fltper_1x1.txt")
  final.data.set$gender <- final.data.set$gender+1 
  
  prop.excess.deaths[i] <- length(which(time.excess.sim<=time.pop.sim & delta.i==1))/sum(delta.i==1)
  
  # EM-based model? 
  
  test.fit.em.default <- rsadd(Surv(time.obs*365.241,delta.i)~age+as.factor(gender)+as.factor(treatment),data=final.data.set,
                               ratetable = nortab,method="EM",bwin=-1,rmap = list(age=age*365.241,sex=gender,year=mdy.date(1,1,year)))
  
  test.fit.em.bwin.0 <- rsadd(Surv(time.obs*365.241,delta.i)~age+as.factor(gender)+as.factor(treatment),data=final.data.set,
                              ratetable = nortab,method="EM",bwin=0,rmap = list(age=age*365.241,sex=gender,year=mdy.date(1,1,year)))
  
  test.fit.em.bwin.1 <- rsadd(Surv(time.obs*365.241,delta.i)~age+as.factor(gender)+as.factor(treatment),data=final.data.set,
                              ratetable = nortab,method="EM",bwin=1,rmap = list(age=age*365.241,sex=gender,year=mdy.date(1,1,year)))
  
  test.fit.em.bwin.100 <- rsadd(Surv(time.obs*365.241,delta.i)~age+as.factor(gender)+as.factor(treatment),data=final.data.set,
                                ratetable = nortab,method="EM",bwin=100,rmap = list(age=age*365.241,sex=gender,year=mdy.date(1,1,year)))
  
  # Poisson GLM? 
  
  test.fit.poisson.1 <- rsadd(Surv(time.obs*365.241,delta.i)~age+as.factor(gender)+as.factor(treatment),data=final.data.set,
                              ratetable = nortab,int=c(0,2,5,10,15,21),method="glm.poi",rmap = list(age=age*365.241,sex=gender,year=mdy.date(1,1,year)))
  test.fit.poisson.2 <- rsadd(Surv(time.obs*365.241,delta.i)~age+as.factor(gender)+as.factor(treatment),data=final.data.set,
                              ratetable = nortab,int=c(0,seq(from=1,to=10,by=1),15,21),method="glm.poi",rmap = list(age=age*365.241,sex=gender,year=mdy.date(1,1,year)))
  
  # Est?ve full maximum likelihood
  
  test.fit.ml.1 <- rsadd(Surv(time.obs*365.241,delta.i)~age+as.factor(gender)+as.factor(treatment),data=final.data.set,
                         ratetable = nortab,int=c(0,2,5,10,15,21),method="max.lik",rmap = list(age=age*365.241,sex=gender,year=mdy.date(1,1,year)))
  test.fit.ml.2 <- rsadd(Surv(time.obs*365.241,delta.i)~age+as.factor(gender)+as.factor(treatment),data=final.data.set,
                         ratetable = nortab,int=c(0,seq(from=1,to=10,by=1),15,21),method="max.lik",rmap = list(age=age*365.241,sex=gender,year=mdy.date(1,1,year)))
  
  # Proportional hazard tests
  standard.res.em.default <- residuals.rsadd(test.fit.em.default)
  standard.res.em.bwin.0 <- residuals.rsadd(test.fit.em.bwin.0)
  standard.res.em.bwin.1 <- residuals.rsadd(test.fit.em.bwin.1)
  standard.res.em.bwin.100 <- residuals.rsadd(test.fit.em.bwin.100)
  standard.res.poisson.1 <- residuals.rsadd(test.fit.poisson.1)
  standard.res.poisson.2 <- residuals.rsadd(test.fit.poisson.2)
  standard.res.ml.1 <- residuals.rsadd(test.fit.ml.1)
  standard.res.ml.2 <- residuals.rsadd(test.fit.ml.2)
  
  br.max.test.em.default <- rs.br(test.fit.em.default,standard.res.em.default)
  br.max.test.em.bwin.0 <- rs.br(test.fit.em.bwin.0,standard.res.em.bwin.0)
  br.max.test.em.bwin.1 <- rs.br(test.fit.em.bwin.1,standard.res.em.bwin.1)
  br.max.test.em.bwin.100 <- rs.br(test.fit.em.bwin.100,standard.res.em.bwin.100)
  br.max.test.poisson.1 <- rs.br(test.fit.poisson.1,standard.res.poisson.1)
  br.max.test.poisson.2 <- rs.br(test.fit.poisson.2,standard.res.poisson.2)
  br.max.test.ml.1 <- rs.br(test.fit.ml.1,standard.res.ml.1)
  br.max.test.ml.2 <- rs.br(test.fit.ml.2,standard.res.ml.2)
  
  br.max.weighted.test.em.default <- rs.br(test.fit.em.default,rho=1,standard.res.em.default)
  br.max.weighted.test.em.bwin.0 <- rs.br(test.fit.em.bwin.0,rho=1,standard.res.em.bwin.0)
  br.max.weighted.test.em.bwin.1 <- rs.br(test.fit.em.bwin.1,rho=1,standard.res.em.bwin.1)
  br.max.weighted.test.em.bwin.100 <- rs.br(test.fit.em.bwin.100,rho=1,standard.res.em.bwin.100)
  br.max.weighted.test.poisson.1 <- rs.br(test.fit.poisson.1,rho=1,standard.res.poisson.1)
  br.max.weighted.test.poisson.2 <- rs.br(test.fit.poisson.2,rho=1,standard.res.poisson.2)
  br.max.weighted.test.ml.1 <- rs.br(test.fit.ml.1,rho=1,standard.res.ml.1)
  br.max.weighted.test.ml.2 <- rs.br(test.fit.ml.2,rho=1,standard.res.ml.2)
  
  cvm.test.em.default <- rs.br(test.fit.em.default,test="cvm",standard.res.em.default)
  cvm.test.em.bwin.0 <- rs.br(test.fit.em.bwin.0,test="cvm",standard.res.em.bwin.0)
  cvm.test.em.bwin.1 <- rs.br(test.fit.em.bwin.1,test="cvm",standard.res.em.bwin.1)
  cvm.test.em.bwin.100 <- rs.br(test.fit.em.bwin.100,test="cvm",standard.res.em.bwin.100)
  cvm.test.poisson.1<- rs.br(test.fit.poisson.1,test="cvm",standard.res.poisson.1)
  cvm.test.poisson.2<- rs.br(test.fit.poisson.2,test="cvm",standard.res.poisson.2)
  cvm.test.ml.1<- rs.br(test.fit.ml.1,test="cvm",standard.res.ml.1)
  cvm.test.ml.2<- rs.br(test.fit.ml.2,test="cvm",standard.res.ml.2)
  
  br.max.test.em.default.age[i] <- br.max.test.em.default$table[1,"p"]<0.05 # alpha=0.05
  br.max.test.em.bwin.0.age[i] <- br.max.test.em.bwin.0$table[1,"p"]<0.05
  br.max.test.em.bwin.1.age[i] <- br.max.test.em.bwin.1$table[1,"p"]<0.05
  br.max.test.em.bwin.100.age[i] <- br.max.test.em.bwin.100$table[1,"p"]<0.05
  br.max.test.poisson.1.age[i] <- br.max.test.poisson.1$table[1,"p"]<0.05
  br.max.test.poisson.2.age[i] <- br.max.test.poisson.2$table[1,"p"]<0.05
  br.max.test.ml.1.age[i] <- br.max.test.ml.1$table[1,"p"]<0.05
  br.max.test.ml.2.age[i] <- br.max.test.ml.2$table[1,"p"]<0.05
  
  br.max.test.em.default.gender[i] <- br.max.test.em.default$table[2,"p"]<0.05 # alpha=0.05
  br.max.test.em.bwin.0.gender[i] <- br.max.test.em.bwin.0$table[2,"p"]<0.05
  br.max.test.em.bwin.1.gender[i] <- br.max.test.em.bwin.1$table[2,"p"]<0.05
  br.max.test.em.bwin.100.gender[i] <- br.max.test.em.bwin.100$table[2,"p"]<0.05
  br.max.test.poisson.1.gender[i] <- br.max.test.poisson.1$table[2,"p"]<0.05
  br.max.test.poisson.2.gender[i] <- br.max.test.poisson.2$table[2,"p"]<0.05
  br.max.test.ml.1.gender[i] <- br.max.test.ml.1$table[2,"p"]<0.05
  br.max.test.ml.2.gender[i] <- br.max.test.ml.2$table[2,"p"]<0.05
  
  br.max.test.em.default.treatment[i] <- br.max.test.em.default$table[3,"p"]<0.05 # alpha=0.05
  br.max.test.em.bwin.0.treatment[i] <- br.max.test.em.bwin.0$table[3,"p"]<0.05
  br.max.test.em.bwin.1.treatment[i] <- br.max.test.em.bwin.1$table[3,"p"]<0.05
  br.max.test.em.bwin.100.treatment[i] <- br.max.test.em.bwin.100$table[3,"p"]<0.05
  br.max.test.poisson.1.treatment[i] <- br.max.test.poisson.1$table[3,"p"]<0.05
  br.max.test.poisson.2.treatment[i] <- br.max.test.poisson.2$table[3,"p"]<0.05
  br.max.test.ml.1.treatment[i] <- br.max.test.ml.1$table[3,"p"]<0.05
  br.max.test.ml.2.treatment[i] <- br.max.test.ml.2$table[3,"p"]<0.05
  
  br.max.weighted.test.em.default.age[i] <- br.max.weighted.test.em.default$table[1,"p"]<0.05 # alpha=0.05
  br.max.weighted.test.em.bwin.0.age[i] <- br.max.weighted.test.em.bwin.0$table[1,"p"]<0.05
  br.max.weighted.test.em.bwin.1.age[i] <- br.max.weighted.test.em.bwin.1$table[1,"p"]<0.05
  br.max.weighted.test.em.bwin.100.age[i] <- br.max.weighted.test.em.bwin.100$table[1,"p"]<0.05
  br.max.weighted.test.poisson.1.age[i] <- br.max.weighted.test.poisson.1$table[1,"p"]<0.05
  br.max.weighted.test.poisson.2.age[i] <- br.max.weighted.test.poisson.2$table[1,"p"]<0.05
  br.max.weighted.test.ml.1.age[i] <- br.max.weighted.test.ml.1$table[1,"p"]<0.05
  br.max.weighted.test.ml.2.age[i] <- br.max.weighted.test.ml.2$table[1,"p"]<0.05
  
  br.max.weighted.test.em.default.gender[i] <- br.max.weighted.test.em.default$table[2,"p"]<0.05 # alpha=0.05
  br.max.weighted.test.em.bwin.0.gender[i] <- br.max.weighted.test.em.bwin.0$table[2,"p"]<0.05
  br.max.weighted.test.em.bwin.1.gender[i] <- br.max.weighted.test.em.bwin.1$table[2,"p"]<0.05
  br.max.weighted.test.em.bwin.100.gender[i] <- br.max.weighted.test.em.bwin.100$table[2,"p"]<0.05
  br.max.weighted.test.poisson.1.gender[i] <- br.max.weighted.test.poisson.1$table[2,"p"]<0.05
  br.max.weighted.test.poisson.2.gender[i] <- br.max.weighted.test.poisson.2$table[2,"p"]<0.05
  br.max.weighted.test.ml.1.gender[i] <- br.max.weighted.test.ml.1$table[2,"p"]<0.05
  br.max.weighted.test.ml.2.gender[i] <- br.max.weighted.test.ml.2$table[2,"p"]<0.05
  
  br.max.weighted.test.em.default.treatment[i] <- br.max.weighted.test.em.default$table[3,"p"]<0.05 # alpha=0.05
  br.max.weighted.test.em.bwin.0.treatment[i] <- br.max.weighted.test.em.bwin.0$table[3,"p"]<0.05
  br.max.weighted.test.em.bwin.1.treatment[i] <- br.max.weighted.test.em.bwin.1$table[3,"p"]<0.05
  br.max.weighted.test.em.bwin.100.treatment[i] <- br.max.weighted.test.em.bwin.100$table[3,"p"]<0.05
  br.max.weighted.test.poisson.1.treatment[i] <- br.max.weighted.test.poisson.1$table[3,"p"]<0.05
  br.max.weighted.test.poisson.2.treatment[i] <- br.max.weighted.test.poisson.2$table[3,"p"]<0.05
  br.max.weighted.test.ml.1.treatment[i] <- br.max.weighted.test.ml.1$table[3,"p"]<0.05
  br.max.weighted.test.ml.2.treatment[i] <- br.max.weighted.test.ml.2$table[3,"p"]<0.05
  
  cvm.test.em.default.age[i] <- cvm.test.em.default$table[1,"p"]<0.05 # alpha=0.05
  cvm.test.em.bwin.0.age[i] <- cvm.test.em.bwin.0$table[1,"p"]<0.05
  cvm.test.em.bwin.1.age[i] <- cvm.test.em.bwin.1$table[1,"p"]<0.05
  cvm.test.em.bwin.100.age[i] <- cvm.test.em.bwin.100$table[1,"p"]<0.05
  cvm.test.poisson.1.age[i] <- cvm.test.poisson.1$table[1,"p"]<0.05
  cvm.test.poisson.2.age[i] <- cvm.test.poisson.2$table[1,"p"]<0.05
  cvm.test.ml.1.age[i] <- cvm.test.ml.1$table[1,"p"]<0.05
  cvm.test.ml.2.age[i] <- cvm.test.ml.2$table[1,"p"]<0.05
  
  cvm.test.em.default.gender[i] <- cvm.test.em.default$table[2,"p"]<0.05 # alpha=0.05
  cvm.test.em.bwin.0.gender[i] <- cvm.test.em.bwin.0$table[2,"p"]<0.05
  cvm.test.em.bwin.1.gender[i] <- cvm.test.em.bwin.1$table[2,"p"]<0.05
  cvm.test.em.bwin.100.gender[i] <- cvm.test.em.bwin.100$table[2,"p"]<0.05
  cvm.test.poisson.1.gender[i] <- cvm.test.poisson.1$table[2,"p"]<0.05
  cvm.test.poisson.2.gender[i] <- cvm.test.poisson.2$table[2,"p"]<0.05
  cvm.test.ml.1.gender[i] <- cvm.test.ml.1$table[2,"p"]<0.05
  cvm.test.ml.2.gender[i] <- cvm.test.ml.2$table[2,"p"]<0.05
  
  cvm.test.em.default.treatment[i] <- cvm.test.em.default$table[3,"p"]<0.05 # alpha=0.05
  cvm.test.em.bwin.0.treatment[i] <- cvm.test.em.bwin.0$table[3,"p"]<0.05
  cvm.test.em.bwin.1.treatment[i] <- cvm.test.em.bwin.1$table[3,"p"]<0.05
  cvm.test.em.bwin.100.treatment[i] <- cvm.test.em.bwin.100$table[3,"p"]<0.05
  cvm.test.poisson.1.treatment[i] <- cvm.test.poisson.1$table[3,"p"]<0.05
  cvm.test.poisson.2.treatment[i] <- cvm.test.poisson.2$table[3,"p"]<0.05
  cvm.test.ml.1.treatment[i] <- cvm.test.ml.1$table[3,"p"]<0.05
  cvm.test.ml.2.treatment[i] <- cvm.test.ml.2$table[3,"p"]<0.05
  
  beta.age.em.default[i] <- as.numeric(test.fit.em.default$coefficients[1])
  beta.age.em.bwin.0[i] <- as.numeric(test.fit.em.bwin.0$coefficients[1])
  beta.age.em.bwin.1[i] <- as.numeric(test.fit.em.bwin.1$coefficients[1])
  beta.age.em.bwin.100[i] <- as.numeric(test.fit.em.bwin.100$coefficients[1])
  beta.age.poisson.1[i] <- as.numeric(test.fit.poisson.1$coefficients[1])
  beta.age.poisson.2[i] <- as.numeric(test.fit.poisson.2$coefficients[1])
  beta.age.ml.1[i] <- as.numeric(test.fit.ml.1$coefficients[1])
  beta.age.ml.2[i] <- as.numeric(test.fit.ml.2$coefficients[1])
  
  beta.gender.em.default[i] <- as.numeric(test.fit.em.default$coefficients[2])
  beta.gender.em.bwin.0[i] <- as.numeric(test.fit.em.bwin.0$coefficients[2])
  beta.gender.em.bwin.1[i] <- as.numeric(test.fit.em.bwin.1$coefficients[2])
  beta.gender.em.bwin.100[i] <- as.numeric(test.fit.em.bwin.100$coefficients[2])
  beta.gender.poisson.1[i] <- as.numeric(test.fit.poisson.1$coefficients[2])
  beta.gender.poisson.2[i] <- as.numeric(test.fit.poisson.2$coefficients[2])
  beta.gender.ml.1[i] <- as.numeric(test.fit.ml.1$coefficients[2])
  beta.gender.ml.2[i] <- as.numeric(test.fit.ml.2$coefficients[2])
  
  beta.treatment.em.default[i] <- as.numeric(test.fit.em.default$coefficients[3])
  beta.treatment.em.bwin.0[i] <- as.numeric(test.fit.em.bwin.0$coefficients[3])
  beta.treatment.em.bwin.1[i] <- as.numeric(test.fit.em.bwin.1$coefficients[3])
  beta.treatment.em.bwin.100[i] <- as.numeric(test.fit.em.bwin.100$coefficients[3])
  beta.treatment.poisson.1[i] <- as.numeric(test.fit.poisson.1$coefficients[3])
  beta.treatment.poisson.2[i] <- as.numeric(test.fit.poisson.2$coefficients[3])
  beta.treatment.ml.1[i] <- as.numeric(test.fit.ml.1$coefficients[3])
  beta.treatment.ml.2[i] <- as.numeric(test.fit.ml.2$coefficients[3])
  
  Lambda0.times.default <- test.fit.em.default$times
  Lambda0.values.default <- test.fit.em.default$Lambda0
  Lambda0.em.default[[i]] <- cbind(Lambda0.times.default/365.241,Lambda0.values.default)
  lambda0.values.default <- test.fit.em.default$lambda0
  lambda0.em.default[[i]] <- cbind(Lambda0.times.default/365.241,lambda0.values.default*365.241)
  
  Lambda0.times.bwin.0 <- test.fit.em.bwin.0$times
  Lambda0.values.bwin.0 <- test.fit.em.bwin.0$Lambda0
  Lambda0.em.bwin.0[[i]] <- cbind(Lambda0.times.bwin.0/365.241,Lambda0.values.bwin.0)
  lambda0.values.bwin.0 <- test.fit.em.bwin.0$lambda0
  lambda0.em.bwin.0[[i]] <- cbind(Lambda0.times.bwin.0/365.241,lambda0.values.bwin.0*365.241)
  
  Lambda0.times.bwin.1 <- test.fit.em.bwin.1$times
  Lambda0.values.bwin.1 <- test.fit.em.bwin.1$Lambda0
  Lambda0.em.bwin.1[[i]] <- cbind(Lambda0.times.bwin.1/365.241,Lambda0.values.bwin.1)
  lambda0.values.bwin.1 <- test.fit.em.bwin.1$lambda0
  lambda0.em.bwin.1[[i]] <- cbind(Lambda0.times.bwin.1/365.241,lambda0.values.bwin.1*365.241)
  
  Lambda0.times.bwin.100 <- test.fit.em.bwin.100$times
  Lambda0.values.bwin.100 <- test.fit.em.bwin.100$Lambda0
  Lambda0.em.bwin.100[[i]] <- cbind(Lambda0.times.bwin.100/365.241,Lambda0.values.bwin.100)
  lambda0.values.bwin.100 <- test.fit.em.bwin.100$lambda0
  lambda0.em.bwin.100[[i]] <- cbind(Lambda0.times.bwin.100/365.241,lambda0.values.bwin.100*365.241)
  
  min.diff.time.obs[i] <- min(diff(sort(time.observed.sim[time.observed.sim<21])))
}

par(mfrow=c(2,4))
hist(beta.age.em.default,xlab = expression(paste(hat(beta)["age"])),
     col="forestgreen",main=expression(paste("EM with automatic choice of bandwidth")))
abline(v=0.05,col="red",lwd=3)
hist(beta.age.em.bwin.0,xlab = expression(paste(hat(beta)["age"])),
     col="gold",main=expression(paste("EM with no smoothing")))
abline(v=0.05,col="red",lwd=3)
hist(beta.age.em.bwin.1,xlab = expression(paste(hat(beta)["age"])),
     col="dodgerblue",main=expression(paste("EM with bwin=1")))
abline(v=0.05,col="red",lwd=3)
hist(beta.age.em.bwin.100,xlab = expression(paste(hat(beta)["age"])),
     col="lightpink",main=expression(paste("EM with bwin=100")))
abline(v=0.05,col="red",lwd=3)
hist(beta.age.poisson.1,xlab = expression(paste(hat(beta)["age"])),
     col="mediumpurple",main=expression(paste("Poisson (Partition 1)")))
abline(v=0.05,col="red",lwd=3)
hist(beta.age.poisson.2,xlab = expression(paste(hat(beta)["age"])),
     col="wheat",main=expression(paste("Poisson (Partition 2)")))
abline(v=0.05,col="red",lwd=3)
hist(beta.age.ml.1,xlab = expression(paste(hat(beta)["age"])),
     col="gainsboro",main=expression(paste("ML (Partition 1)")))
abline(v=0.05,col="red",lwd=3)
hist(beta.age.ml.2,xlab = expression(paste(hat(beta)["age"])),
     col="darkturquoise",main=expression(paste("ML (Partition 2)")))
abline(v=0.05,col="red",lwd=3)
c(mean(beta.age.em.default),mean(beta.age.em.bwin.0),mean(beta.age.em.bwin.1),mean(beta.age.em.bwin.100),mean(beta.age.poisson.1),mean(beta.age.poisson.2),mean(beta.age.ml.1),mean(beta.age.ml.2))
c(sd(beta.age.em.default),sd(beta.age.em.bwin.0),sd(beta.age.em.bwin.1),sd(beta.age.em.bwin.100),sd(beta.age.poisson.1),sd(beta.age.poisson.2),sd(beta.age.ml.1),sd(beta.age.ml.2))


hist(beta.gender.em.default,xlab = expression(paste(hat(beta)["gender"])),
     col="forestgreen",main=expression(paste("EM with automatic choice of bandwidth")))
abline(v=0.1,col="red",lwd=3)
hist(beta.gender.em.bwin.0,xlab = expression(paste(hat(beta)["gender"])),
     col="gold",main=expression(paste("EM with no smoothing")))
abline(v=0.1,col="red",lwd=3)
hist(beta.gender.em.bwin.1,xlab = expression(paste(hat(beta)["gender"])),
     col="dodgerblue",main=expression(paste("EM with bwin=1")))
abline(v=0.1,col="red",lwd=3)
hist(beta.gender.em.bwin.100,xlab = expression(paste(hat(beta)["gender"])),
     col="lightpink",main=expression(paste("EM with bwin=100")))
abline(v=0.1,col="red",lwd=3)
hist(beta.gender.poisson.1,xlab = expression(paste(hat(beta)["gender"])),
     col="mediumpurple",main=expression(paste("Poisson (Partition 1)")))
abline(v=0.1,col="red",lwd=3)
hist(beta.gender.poisson.2,xlab = expression(paste(hat(beta)["gender"])),
     col="wheat",main=expression(paste("Poisson (Partition 2)")))
abline(v=0.1,col="red",lwd=3)
hist(beta.gender.ml.1,xlab = expression(paste(hat(beta)["gender"])),
     col="gainsboro",main=expression(paste("ML (Partition 1)")))
abline(v=0.1,col="red",lwd=3)
hist(beta.gender.ml.2,xlab = expression(paste(hat(beta)["gender"])),
     col="darkturquoise",main=expression(paste("ML (Partition 2)")))
abline(v=0.1,col="red",lwd=3)
c(mean(beta.gender.em.default),mean(beta.gender.em.bwin.0),mean(beta.gender.em.bwin.1),mean(beta.gender.em.bwin.100),mean(beta.gender.poisson.1),mean(beta.gender.poisson.2),mean(beta.gender.ml.1),mean(beta.gender.ml.2))
c(sd(beta.gender.em.default),sd(beta.gender.em.bwin.0),sd(beta.gender.em.bwin.1),sd(beta.gender.em.bwin.100),sd(beta.gender.poisson.1),sd(beta.gender.poisson.2),sd(beta.gender.ml.1),sd(beta.gender.ml.2))

hist(beta.treatment.em.default,xlab = expression(paste(hat(beta)["treatment"])),
     col="forestgreen",main=expression(paste("EM with automatic choice of bandwidth")))
abline(v=0.5,col="red",lwd=3)
hist(beta.treatment.em.bwin.0,xlab = expression(paste(hat(beta)["treatment"])),
     col="gold",main=expression(paste("EM with no smoothing")))
abline(v=0.5,col="red",lwd=3)
hist(beta.treatment.em.bwin.1,xlab = expression(paste(hat(beta)["treatment"])),
     col="dodgerblue",main=expression(paste("EM with bwin=1")))
abline(v=0.5,col="red",lwd=3)
hist(beta.treatment.em.bwin.100,xlab = expression(paste(hat(beta)["treatment"])),
     col="lightpink",main=expression(paste("EM with bwin=100")))
abline(v=0.5,col="red",lwd=3)
hist(beta.treatment.poisson.1,xlab = expression(paste(hat(beta)["treatment"])),
     col="mediumpurple",main=expression(paste("Poisson (Partition 1)")))
abline(v=0.5,col="red",lwd=3)
hist(beta.treatment.poisson.2,xlab = expression(paste(hat(beta)["treatment"])),
     col="wheat",main=expression(paste("Poisson (Partition 2)")))
abline(v=0.5,col="red",lwd=3)
hist(beta.treatment.ml.1,xlab = expression(paste(hat(beta)["treatment"])),
     col="gainsboro",main=expression(paste("ML (Partition 1)")))
abline(v=0.5,col="red",lwd=3)
hist(beta.treatment.ml.2,xlab = expression(paste(hat(beta)["treatment"])),
     col="darkturquoise",main=expression(paste("ML (Partition 2)")))
abline(v=0.5,col="red",lwd=3)
c(mean(beta.treatment.em.default),mean(beta.treatment.em.bwin.0),mean(beta.treatment.em.bwin.1),mean(beta.treatment.em.bwin.100),mean(beta.treatment.poisson.1),mean(beta.treatment.poisson.2),mean(beta.treatment.ml.1),mean(beta.treatment.ml.2))
c(sd(beta.treatment.em.default),sd(beta.treatment.em.bwin.0),sd(beta.treatment.em.bwin.1),sd(beta.treatment.em.bwin.100),sd(beta.treatment.poisson.1),sd(beta.treatment.poisson.2),sd(beta.treatment.ml.1),sd(beta.treatment.ml.2))

# Plots of baseline cumulative hazard/hazard of default EM 
plot(Lambda0.em.default[[1]][,1],Lambda0.em.default[[1]][,2],type="l",xlab="t (in years)",
     ,ylab="",main="EM with automatic choice of bandwidth")
title(ylab=expression(paste(hat(Lambda)[0](t))),line=2.5)
for(i in 2:n.sim.models){
  lines(Lambda0.em.default[[i]][,1],Lambda0.em.default[[i]][,2],col=i)
}
lines(t.grid,phi*t.grid^lambda,col="red",lwd=3)

plot(lowess(lambda0.em.bwin.0[[1]][,1],lambda0.em.bwin.0[[1]][,2],f=0.15),type="l",ylim=c(0,0.025),
     xlab="t (in years)",ylab="",main="EM with automatic choice of bandwidth")
title(ylab=expression(paste(hat(lambda)[0](t))),line=2.5)
for(i in 2:n.sim.models){
  lines(lowess(lambda0.em.default[[i]][,1],lambda0.em.default[[i]][,2],f=0.5),col=i)
}
lines(t.grid,lambda*phi*t.grid^(lambda-1),col="red",lwd=3)

# Plots of baseline cumulative hazard/hazard of EM with no smoothing 
plot(Lambda0.em.bwin.0[[1]][,1],Lambda0.em.bwin.0[[1]][,2],type="l",xlab="t (in years)",
     ,ylab="",main="EM with no smoothing")
title(ylab=expression(paste(hat(Lambda)[0](t))),line=2.5)
for(i in 2:n.sim.models){
  lines(Lambda0.em.bwin.0[[i]][,1],Lambda0.em.bwin.0[[i]][,2],col=i)
}
lines(t.grid,phi*t.grid^lambda,col="red",lwd=3)

plot(lowess(lambda0.em.bwin.0[[1]][,1],lambda0.em.bwin.0[[1]][,2],f=0.15),type="l",ylim=c(0,0.025),
     xlab="t (in years)",ylab="",main="EM with no smoothing")
title(ylab=expression(paste(hat(lambda)[0](t))),line=2.5)
for(i in 2:n.sim.models){
  lines(lowess(lambda0.em.bwin.0[[i]][,1],lambda0.em.bwin.0[[i]][,2],f=0.15),col=i)
}
lines(t.grid,lambda*phi*t.grid^(lambda-1),col="red",lwd=3)

# Plots of baseline cumulative hazard/hazard of EM with bwin=1 
plot(Lambda0.em.bwin.1[[1]][,1],Lambda0.em.bwin.1[[1]][,2],type="l",xlab="t (in years)",
     ,ylab="",main="EM with bwin=1")
title(ylab=expression(paste(hat(Lambda)[0](t))),line=2.5)
for(i in 2:n.sim.models){
  lines(Lambda0.em.bwin.1[[i]][,1],Lambda0.em.bwin.1[[i]][,2],col=i)
}
lines(t.grid,phi*t.grid^lambda,col="red",lwd=3)

plot(lowess(lambda0.em.bwin.1[[1]][,1],lambda0.em.bwin.1[[1]][,2],f=0.15),type="l",ylim=c(0,0.025),
     xlab="t (in years)",ylab="",main="EM with bwin=1")
title(ylab=expression(paste(hat(lambda)[0](t))),line=2.5)
for(i in 2:n.sim.models){
  lines(lowess(lambda0.em.bwin.1[[i]][,1],lambda0.em.bwin.1[[i]][,2],f=0.15),col=i)
}
lines(t.grid,lambda*phi*t.grid^(lambda-1),col="red",lwd=3)

# Plots of baseline cumulative hazard/hazard of EM with bwin=100
plot(Lambda0.em.bwin.100[[1]][,1],Lambda0.em.bwin.100[[1]][,2],type="l",xlab="t (in years)",
     ,ylab="",main="EM with bwin=100")
title(ylab=expression(paste(hat(Lambda)[0](t))),line=2.5)
for(i in 2:n.sim.models){
  lines(Lambda0.em.bwin.100[[i]][,1],Lambda0.em.bwin.100[[i]][,2],col=i)
}
lines(t.grid,phi*t.grid^lambda,col="red",lwd=3)

plot(lowess(lambda0.em.bwin.1[[1]][,1],lambda0.em.bwin.1[[1]][,2],f=0.15),type="l",ylim=c(0,0.025),
     xlab="t (in years)",ylab="",main="EM with bwin=100")
title(ylab=expression(paste(hat(lambda)[0](t))),line=2.5)
for(i in 2:n.sim.models){
  lines(lowess(lambda0.em.bwin.100[[i]][,1],lambda0.em.bwin.100[[i]][,2],f=0.15),col=i)
}
lines(t.grid,lambda*phi*t.grid^(lambda-1),col="red",lwd=3)

# Proportional excess hazard tests summary 

# Normal maximum value of Brownian bridge
prop.rej.br.max.em.default.age <- sum(br.max.test.em.default.age)/n.sim.models
prop.rej.br.max.em.bwin.0.age <- sum(br.max.test.em.bwin.0.age)/n.sim.models
prop.rej.br.max.em.bwin.1.age <-sum(br.max.test.em.bwin.1.age)/n.sim.models
prop.rej.br.max.em.bwin.100.age <-sum(br.max.test.em.bwin.100.age)/n.sim.models
prop.rej.br.max.poisson.1.age <-sum(br.max.test.poisson.1.age)/n.sim.models
prop.rej.br.max.poisson.2.age <-sum(br.max.test.poisson.2.age)/n.sim.models
prop.rej.br.max.ml.1.age <-sum(br.max.test.ml.1.age)/n.sim.models
prop.rej.br.max.ml.2.age <-sum(br.max.test.ml.2.age)/n.sim.models
c(prop.rej.br.max.em.default.age,prop.rej.br.max.em.bwin.0.age,prop.rej.br.max.em.bwin.1.age
  ,prop.rej.br.max.em.bwin.100.age,prop.rej.br.max.poisson.1.age,prop.rej.br.max.poisson.2.age,
  prop.rej.br.max.ml.1.age,prop.rej.br.max.ml.2.age)

prop.rej.br.max.em.default.gender <- sum(br.max.test.em.default.gender)/n.sim.models
prop.rej.br.max.em.bwin.0.gender <- sum(br.max.test.em.bwin.0.gender)/n.sim.models
prop.rej.br.max.em.bwin.1.gender <-sum(br.max.test.em.bwin.1.gender)/n.sim.models
prop.rej.br.max.em.bwin.100.gender <-sum(br.max.test.em.bwin.100.gender)/n.sim.models
prop.rej.br.max.poisson.1.gender <-sum(br.max.test.poisson.1.gender)/n.sim.models
prop.rej.br.max.poisson.2.gender <-sum(br.max.test.poisson.2.gender)/n.sim.models
prop.rej.br.max.ml.1.gender <-sum(br.max.test.ml.1.gender)/n.sim.models
prop.rej.br.max.ml.2.gender <-sum(br.max.test.ml.2.gender)/n.sim.models
c(prop.rej.br.max.em.default.gender,prop.rej.br.max.em.bwin.0.gender,prop.rej.br.max.em.bwin.1.gender
  ,prop.rej.br.max.em.bwin.100.gender,prop.rej.br.max.poisson.1.gender,prop.rej.br.max.poisson.2.gender,
  prop.rej.br.max.ml.1.gender,prop.rej.br.max.ml.2.gender)

prop.rej.br.max.em.default.treatment <- sum(br.max.test.em.default.treatment)/n.sim.models
prop.rej.br.max.em.bwin.0.treatment <- sum(br.max.test.em.bwin.0.treatment)/n.sim.models
prop.rej.br.max.em.bwin.1.treatment <-sum(br.max.test.em.bwin.1.treatment)/n.sim.models
prop.rej.br.max.em.bwin.100.treatment <-sum(br.max.test.em.bwin.100.treatment)/n.sim.models
prop.rej.br.max.poisson.1.treatment <-sum(br.max.test.poisson.1.treatment)/n.sim.models
prop.rej.br.max.poisson.2.treatment <-sum(br.max.test.poisson.2.treatment)/n.sim.models
prop.rej.br.max.ml.1.treatment <-sum(br.max.test.ml.1.treatment)/n.sim.models
prop.rej.br.max.ml.2.treatment <-sum(br.max.test.ml.2.treatment)/n.sim.models
c(prop.rej.br.max.em.default.treatment,prop.rej.br.max.em.bwin.0.treatment,prop.rej.br.max.em.bwin.1.treatment
  ,prop.rej.br.max.em.bwin.100.treatment,prop.rej.br.max.poisson.1.treatment,prop.rej.br.max.poisson.2.treatment,
  prop.rej.br.max.ml.1.treatment,prop.rej.br.max.ml.2.treatment)

# Weighted version of maximum value of Brownian bridge
prop.rej.br.max.weighted.em.default.age <- sum(br.max.weighted.test.em.default.age)/n.sim.models
prop.rej.br.max.weighted.em.bwin.0.age <- sum(br.max.weighted.test.em.bwin.0.age)/n.sim.models
prop.rej.br.max.weighted.em.bwin.1.age <-sum(br.max.weighted.test.em.bwin.1.age)/n.sim.models
prop.rej.br.max.weighted.em.bwin.100.age <-sum(br.max.weighted.test.em.bwin.100.age)/n.sim.models
prop.rej.br.max.weighted.poisson.1.age <-sum(br.max.weighted.test.poisson.1.age)/n.sim.models
prop.rej.br.max.weighted.poisson.2.age <-sum(br.max.weighted.test.poisson.2.age)/n.sim.models
prop.rej.br.max.weighted.ml.1.age <-sum(br.max.weighted.test.ml.1.age)/n.sim.models
prop.rej.br.max.weighted.ml.2.age <-sum(br.max.weighted.test.ml.2.age)/n.sim.models
c(prop.rej.br.max.weighted.em.default.age,prop.rej.br.max.weighted.em.bwin.0.age,prop.rej.br.max.weighted.em.bwin.1.age
  ,prop.rej.br.max.weighted.em.bwin.100.age,prop.rej.br.max.weighted.poisson.1.age,prop.rej.br.max.weighted.poisson.2.age,
  prop.rej.br.max.weighted.ml.1.age,prop.rej.br.max.weighted.ml.2.age)

prop.rej.br.max.weighted.em.default.gender <- sum(br.max.weighted.test.em.default.gender)/n.sim.models
prop.rej.br.max.weighted.em.bwin.0.gender <- sum(br.max.weighted.test.em.bwin.0.gender)/n.sim.models
prop.rej.br.max.weighted.em.bwin.1.gender <-sum(br.max.weighted.test.em.bwin.1.gender)/n.sim.models
prop.rej.br.max.weighted.em.bwin.100.gender <-sum(br.max.weighted.test.em.bwin.100.gender)/n.sim.models
prop.rej.br.max.weighted.poisson.1.gender <-sum(br.max.weighted.test.poisson.1.gender)/n.sim.models
prop.rej.br.max.weighted.poisson.2.gender <-sum(br.max.weighted.test.poisson.2.gender)/n.sim.models
prop.rej.br.max.weighted.ml.1.gender <-sum(br.max.weighted.test.ml.1.gender)/n.sim.models
prop.rej.br.max.weighted.ml.2.gender <-sum(br.max.weighted.test.ml.2.gender)/n.sim.models
c(prop.rej.br.max.weighted.em.default.gender,prop.rej.br.max.weighted.em.bwin.0.gender,prop.rej.br.max.weighted.em.bwin.1.gender
  ,prop.rej.br.max.weighted.em.bwin.100.gender,prop.rej.br.max.weighted.poisson.1.gender,prop.rej.br.max.weighted.poisson.2.gender,
  prop.rej.br.max.weighted.ml.1.gender,prop.rej.br.max.weighted.ml.2.gender)

prop.rej.br.max.weighted.em.default.treatment <- sum(br.max.weighted.test.em.default.treatment)/n.sim.models
prop.rej.br.max.weighted.em.bwin.0.treatment <- sum(br.max.weighted.test.em.bwin.0.treatment)/n.sim.models
prop.rej.br.max.weighted.em.bwin.1.treatment <-sum(br.max.weighted.test.em.bwin.1.treatment)/n.sim.models
prop.rej.br.max.weighted.em.bwin.100.treatment <-sum(br.max.weighted.test.em.bwin.100.treatment)/n.sim.models
prop.rej.br.max.weighted.poisson.1.treatment <-sum(br.max.weighted.test.poisson.1.treatment)/n.sim.models
prop.rej.br.max.weighted.poisson.2.treatment <-sum(br.max.weighted.test.poisson.2.treatment)/n.sim.models
prop.rej.br.max.weighted.ml.1.treatment <-sum(br.max.weighted.test.ml.1.treatment)/n.sim.models
prop.rej.br.max.weighted.ml.2.treatment <-sum(br.max.weighted.test.ml.2.treatment)/n.sim.models
c(prop.rej.br.max.weighted.em.default.treatment,prop.rej.br.max.weighted.em.bwin.0.treatment,prop.rej.br.max.weighted.em.bwin.1.treatment
  ,prop.rej.br.max.weighted.em.bwin.100.treatment,prop.rej.br.max.weighted.poisson.1.treatment,prop.rej.br.max.weighted.poisson.2.treatment,
  prop.rej.br.max.weighted.ml.1.treatment,prop.rej.br.max.weighted.ml.2.treatment)

# Cramer von Mises 
prop.rej.cvm.em.default.age <- sum(cvm.test.em.default.age)/n.sim.models
prop.rej.cvm.em.bwin.0.age <- sum(cvm.test.em.bwin.0.age)/n.sim.models
prop.rej.cvm.em.bwin.1.age <-sum(cvm.test.em.bwin.1.age)/n.sim.models
prop.rej.cvm.em.bwin.100.age <-sum(cvm.test.em.bwin.100.age)/n.sim.models
prop.rej.cvm.poisson.1.age <-sum(cvm.test.poisson.1.age)/n.sim.models
prop.rej.cvm.poisson.2.age <-sum(cvm.test.poisson.2.age)/n.sim.models
prop.rej.cvm.ml.1.age <-sum(cvm.test.ml.1.age)/n.sim.models
prop.rej.cvm.ml.2.age <-sum(cvm.test.ml.2.age)/n.sim.models
c(prop.rej.cvm.em.default.age,prop.rej.cvm.em.bwin.0.age,prop.rej.cvm.em.bwin.1.age
  ,prop.rej.cvm.em.bwin.100.age,prop.rej.cvm.poisson.1.age,prop.rej.cvm.poisson.2.age,
  prop.rej.cvm.ml.1.age,prop.rej.cvm.ml.2.age)

prop.rej.cvm.em.default.gender <- sum(cvm.test.em.default.gender)/n.sim.models
prop.rej.cvm.em.bwin.0.gender <- sum(cvm.test.em.bwin.0.gender)/n.sim.models
prop.rej.cvm.em.bwin.1.gender <-sum(cvm.test.em.bwin.1.gender)/n.sim.models
prop.rej.cvm.em.bwin.100.gender <-sum(cvm.test.em.bwin.100.gender)/n.sim.models
prop.rej.cvm.poisson.1.gender <-sum(cvm.test.poisson.1.gender)/n.sim.models
prop.rej.cvm.poisson.2.gender <-sum(cvm.test.poisson.2.gender)/n.sim.models
prop.rej.cvm.ml.1.gender <-sum(cvm.test.ml.1.gender)/n.sim.models
prop.rej.cvm.ml.2.gender <-sum(cvm.test.ml.2.gender)/n.sim.models
c(prop.rej.cvm.em.default.gender,prop.rej.cvm.em.bwin.0.gender,prop.rej.cvm.em.bwin.1.gender
  ,prop.rej.cvm.em.bwin.100.gender,prop.rej.cvm.poisson.1.gender,prop.rej.cvm.poisson.2.gender,
  prop.rej.cvm.ml.1.gender,prop.rej.cvm.ml.2.gender)

prop.rej.cvm.em.default.treatment <- sum(cvm.test.em.default.treatment)/n.sim.models
prop.rej.cvm.em.bwin.0.treatment <- sum(cvm.test.em.bwin.0.treatment)/n.sim.models
prop.rej.cvm.em.bwin.1.treatment <-sum(cvm.test.em.bwin.1.treatment)/n.sim.models
prop.rej.cvm.em.bwin.100.treatment <-sum(cvm.test.em.bwin.100.treatment)/n.sim.models
prop.rej.cvm.poisson.1.treatment <-sum(cvm.test.poisson.1.treatment)/n.sim.models
prop.rej.cvm.poisson.2.treatment <-sum(cvm.test.poisson.2.treatment)/n.sim.models
prop.rej.cvm.ml.1.treatment <-sum(cvm.test.ml.1.treatment)/n.sim.models
prop.rej.cvm.ml.2.treatment <-sum(cvm.test.ml.2.treatment)/n.sim.models
c(prop.rej.cvm.em.default.treatment,prop.rej.cvm.em.bwin.0.treatment,prop.rej.cvm.em.bwin.1.treatment
  ,prop.rej.cvm.em.bwin.100.treatment,prop.rej.cvm.poisson.1.treatment,prop.rej.cvm.poisson.2.treatment,
  prop.rej.cvm.ml.1.treatment,prop.rej.cvm.ml.2.treatment)

# Check why the weighted version does not work in this case

which(is.na(br.max.weighted.test.em.bwin.0.age))
min.diff.time.obs[which(is.na(br.max.weighted.test.em.bwin.0.age))]
which(!(min.diff.time.obs>10e-8))

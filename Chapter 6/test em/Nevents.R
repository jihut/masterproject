##############################################################################
## These functions implement the method described in the Appendix of
## Gandy, Kvaloy, Bottle & Zhou (2010), Biometrika.  The code has not
## been extensively tested - so it should be treated with caution.
## Please let me know of any error/problems that you encounter.
##
## Please cite the underlying paper if you use the method in your own work!
##
## Copyright (C)2010 Axel Gandy
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
###########################################################################

##
##Computes E(N(\tau)) for given c
##
getENevents <- function(rho=0.8,lambda=1,c,S0=0){
  if (S0>=c) stop("headstart X0 must be less than threshold c")
  if (rho==1) stop("rho cannot equal 1")

  c <- c/abs(log(rho))
  S0 <- S0/abs(log(rho))
  a <- abs(log(rho)/(rho-1))
  if (rho>1)
    getEPNeventsup(a=a,lambda=lambda,c=c,S0=S0)
  else
    getEPNeventsdown(a=a,lambda=lambda,c=c,S0=S0)
}

##
##Computes P(N(\tau)<=NMax) 
##
getPNevents <- function(rho=0.8,lambda=1,c,S0=0,Nmax){
  if (S0>=c) stop("headstart X0 must be less than threshold c")
  if (rho==1) stop("rho cannot equal 1")

  c <- c/abs(log(rho))
  S0 <- S0/abs(log(rho))
  a <- abs(log(rho)/(rho-1))
  if (rho>1)
    getEPNeventsup(a=a,lambda=lambda,c=c,S0=S0,Nmax=Nmax)
  else
    getEPNeventsdown(a=a,lambda=lambda,c=c,S0=S0,Nmax=Nmax)
}

getEPNeventsup <- function(a,lambda=1,c,S0=0,Nmax){
  if (c<0) stop("c must be positive")
  if (S0>=c) stop("headstart S0 must be less than threshold c")
  if (a<=0) stop("a must be positive")
  if (!missing(Nmax)&&Nmax<0) stop("Nmax must be nonnegative")

  if (missing(Nmax)){
    if (c<= 1) return(1) #the first jump will hit the thresholds
  }else{
    if (c<=1) return(Nmax>=1)
  }
  
  k <- ceiling(c-1);
  xi <-a*(c-k)
  eta <- a-xi
  main <- getdefslargerone(k=k,lambda=lambda,xi=xi,eta=eta)   
  nu <- ceiling(c-S0)

  if (!missing(Nmax)&&nu>Nmax) return(0)

  t0 <- a*(S0-(c-nu))
  if (S0>c-k){
    pstart <- c(rep(0,k-nu),(lambda*t0)^(0:(nu-1))*exp(-lambda*t0)/factorial(0:(nu-1)))
  }else{
    start <- getdefslargerone(k=k,lambda=lambda,xi=a*S0,eta)
    pstart <- start$Q[1,]
  }
  if (missing(Nmax))
    nu+pstart%*%solve(diag(rep(1,k))-main$Q,rep(1,k))
  else{
    Nmax <- floor(Nmax)
    1-pstart%*%matrix.power(main$Q,Nmax-nu)%*%rep(1,k)
  }
}

getdefslargerone <- function(k,lambda,xi,eta){
  if (k==0) stop("k==0 not implemented")
  Q <- matrix(0,nrow=k,ncol=k)
  a <- xi+eta
  pa <- (lambda*a)^(0:k)*exp(-lambda*a)/factorial(0:k)
  peta <- (lambda*eta)^(0:k)*exp(-lambda*eta)/factorial(0:k)
  p0xi <- exp(-lambda*xi)
  Q[1,] <- pa[-1]+p0xi*(peta[-(k+1)]-peta[-1])
  if (k>1){
    for (i in 2:k)
      for (j in (i-1):k)
        Q[i,j] <- pa[j-i+2]
  }
  list(Q=Q)
}

getEPNeventsdown <- function(a,lambda=1,c,S0=0,Nmax){
  if (c<0) stop("c must be positive")
  if (S0>=c) stop("headstart S0 must be less than threshold c")
  if (a<=0) stop("a must be positive")
  if (!missing(Nmax)&&Nmax<0) stop("Nmax must be nonnegative")
  
  k <- ceiling(c-1);
  xi <-a*(c-k)
  eta <- a-xi
  main <- getdefssmallerone(k=k,lambda=lambda,xi=xi,eta=eta)
  nu <- ceiling(c-S0)
  
  t0 <- a*(S0-(c-nu))  
  
  start <- getdefssmallerone(k=k,lambda=lambda,xi=min(a-t0,xi),eta=max(a-t0-xi,0))
  
  pstart <- ifelse((1:(k+1))==(k+2-nu),1,0)%*%start$Q
  
  if (missing(Nmax))
    pstart%*%matrix.power(main$Q,nu-1)%*%solve(diag(rep(1,k+1))-main$Q,rep(1,k+1))
  else {
    Nmax <- floor(Nmax)
    1-pstart%*%matrix.power(main$Q,nu-1+Nmax)%*%rep(1,k+1)
  }
}

getdefssmallerone <- function(k,lambda,xi,eta){
  QA <- matrix(0,nrow=k+1,ncol=k+1)
  pxi <- (lambda*xi)^(0:k)*exp(-lambda*xi)/factorial(0:k)
  if (k>=1){
    for (j in 0:(k-1)){
      for (i in 2:(k+1-j))
        QA[i+j,i] <- pxi[j+2]
    }
    for (i in 1:k)
      QA[i,i+1] <- pxi[1]
  }
  QA[,1] <- 1-rowSums(QA)-c(rep(0,k),pxi[1])
  
  QB <- matrix(0,nrow=k+1,ncol=k+1)
  peta <- (lambda*eta)^(0:(k+1))*exp(-lambda*eta)/factorial(0:(k+1))
  if (k>=1){
    for (j in 0:(k-1)){
      for (i in 2:(k+1-j))
        QB[i+j,i] <- peta[j+1]
    }
  }
  QB[,1] <- 1-rowSums(QB)
  
  list(Q=QB%*%QA)
}  

matrix.power <- function(mat, n){
  if (n==1) return(mat)
  result <- diag(1,ncol(mat))
  while (n > 0) {
    if (n %% 2 != 0){
      result <- result %*% mat
      n <- n - 1
    }
    mat <- mat%*%mat
    n <- n/2
  }
  result
}

##
##compute thresholds to get a desired in-control perfromance.
##
getthresholdENevents <- function(rho,headstartfrac=0,EN0=1000){
  uniroot(function(x) getENevents(rho=rho,lambda=1,c=x,S0=x*headstartfrac)-EN0,c(1,10))$root
}

getthresholdPNevents <- function(rho,headstartfrac=0,Nmax=100,PN0=0.01){
  f <- function(x) getPNevents(rho=rho,lambda=1,c=x,Nmax,S0=x*headstartfrac)-PN0
  cup <- 10;
  while (f(cup)>0){
    cup <- cup*1.25
    if (cup>1e5) stop("cup getting too large")
  }
  clow <- 1;
  while (f(clow)<0){
    clow <- clow*0.75
    if (clow<1e-5) stop("clow getting too small")
  }
  uniroot(f,c(clow,cup))$root
}

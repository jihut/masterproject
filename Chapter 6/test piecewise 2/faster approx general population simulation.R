round2 = function(x, n) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5 + sqrt(.Machine$double.eps)
  z = trunc(z)
  z = z/10^n
  z*posneg
}

haz.output.male <- function(age,start.year,end.year,arrival.time,n=8){ # n - "precision" --> can get error if smaller n
  
  age <- round2(age,n)
  arrival.time <- round2(arrival.time,n)
  
  step.age <- round2(ceiling(age)-age,n)
  step.year <- round2(ceiling(start.year+arrival.time)-(start.year+arrival.time),n)
  max.follow.up <- round2((end.year-start.year-arrival.time),n)
  
  if(max.follow.up< round2(step.age+1,n)){
    change.age.time.vec <- c(step.age)
  } else {
    change.age.time.vec <- c(step.age,(seq(from=round2(step.age+1,n),to=max.follow.up,by=1)))
  }
  
  if(max.follow.up< round2(step.year+1,n)){
    change.year.time.vec <- c(step.year)
  } else {
    change.year.time.vec <- c(step.year,(seq(from=round2(step.year+1,n),to=max.follow.up,by=1)))
  }
  
  ### If step.age is equal to step.year ###
  if(step.age==step.year){
    if(step.age==0){
      store.matrix.equal <- matrix(NA,ncol=2,nrow = length(change.age.time.vec)-1)
      store.matrix.equal[,1] <- change.age.time.vec[1:(length(change.age.time.vec)-1)]
      age.throughout.vec <- floor(round2(store.matrix.equal[,1]+age,n))
      if(age.throughout.vec[1]<110){
        index <- which(pop.data.male$Age==age.throughout.vec[1] & pop.data.male$Year==floor(round2(start.year+arrival.time,n)))
      } else {
        index <- which(pop.data.male$Age=="110+" & pop.data.male$Year==floor(round2(start.year+arrival.time,n)))
      }
      index.age.over.110 <- index+(sum(age.throughout.vec<110))*112
      if(sum(age.throughout.vec<=110)!=0){
        for(i in 1:sum(age.throughout.vec<=110)){
          store.matrix.equal[i,2] <- pop.data.male[index+(i-1)*112,"mx"]
        }
      }
      if(sum(age.throughout.vec>110)!=0){
        if(sum(age.throughout.vec<=110)!=0){
          for(i in 1:sum(age.throughout.vec>110)){
            store.matrix.equal[i+sum(age.throughout.vec<=110),2] <- pop.data.male[index.age.over.110+(i)*111,"mx"]
          }
        } else {
          for(i in 1:sum(age.throughout.vec>110)){
            store.matrix.equal[i+sum(age.throughout.vec<=110),2] <- pop.data.male[index.age.over.110+(i-1)*111,"mx"]
          }
        }
      }
    } else {
      if(length(change.age.time.vec)==1){
        store.matrix.equal <- matrix(NA,ncol=2,nrow=1)
        store.matrix.equal[1,1] <- 0
        age.floor <- floor(age)
        if(age.floor<110){
          index <- which(pop.data.male$Age==age.floor & pop.data.male$Year==floor(round2(start.year+arrival.time,n)))
        } else {
          index <- which(pop.data.male$Age=="110+" & pop.data.male$Year==floor(round2(start.year+arrival.time,n)))
        }
        store.matrix.equal[1,2] <- pop.data.male[index,"mx"]
      } else {
        store.matrix.equal <- matrix(NA,ncol=2,nrow=length(change.age.time.vec))
        store.matrix.equal[,1] <- c(0,change.age.time.vec[1:(length(change.age.time.vec)-1)])
        #store.matrix.equal[,2] <- vec.pop.hazard.func.in.sim(age,gender,arrival.time,start.year,store.matrix.equal[,1])
        age.throughout.vec <- floor(round2(store.matrix.equal[,1]+age,n))
        if(age.throughout.vec[1]<110){
          index <- which(pop.data.male$Age==age.throughout.vec[1] & pop.data.male$Year==floor(round2(start.year+arrival.time,n)))
        } else {
          index <- which(pop.data.male$Age=="110+" & pop.data.male$Year==floor(round2(start.year+arrival.time,n)))
        }
        index.age.over.110 <- index+(sum(age.throughout.vec<110))*112
        if(sum(age.throughout.vec<=110)!=0){
          for(i in 1:sum(age.throughout.vec<=110)){
            store.matrix.equal[i,2] <- pop.data.male[index+(i-1)*112,"mx"]
          }
        }
        if(sum(age.throughout.vec>110)!=0){
          if(sum(age.throughout.vec<=110)!=0){
            for(i in 1:sum(age.throughout.vec>110)){
              store.matrix.equal[i+sum(age.throughout.vec<=110),2] <- pop.data.male[index.age.over.110+(i)*111,"mx"]
            }
          } else {
            for(i in 1:sum(age.throughout.vec>110)){
              store.matrix.equal[i+sum(age.throughout.vec<=110),2] <- pop.data.male[index.age.over.110+(i-1)*111,"mx"]
            }
          }
        }
      }
    }
    return(store.matrix.equal)
  }
  
  ### If step.age is different from step.year ### 
  
  if(step.age!=step.year){
    
    if(step.age < step.year){
      small.time.step <- step.age
      small.time.vec <- change.age.time.vec
      large.time.step <- step.year
      large.time.vec <- change.year.time.vec[1:(length(change.year.time.vec)-1)]
      
      if(small.time.step!=0){
        sort.time.vec <- sort(c(0,small.time.vec,large.time.vec))
        small.time.floor.vec.plus.age <- floor(round2(small.time.vec+age,n))
        year.after.first.age.change <- floor(round2(start.year+arrival.time+small.time.vec[1],n))
      } else {
        sort.time.vec <- sort(c(small.time.vec,large.time.vec))
        small.time.floor.vec.plus.age <- floor(round2(small.time.vec+age,n))[-1]
        year.after.first.age.change <- floor(round2(start.year+arrival.time+small.time.vec[2],n))
      }
      
      if(length(change.year.time.vec)==1){
        sort.time.vec <- sort.time.vec[-length(sort.time.vec)]
      }
      
      store.matrix <- matrix(NA,nrow=length(sort.time.vec),ncol=2)
      store.matrix[,1] <- sort.time.vec[1:(length(sort.time.vec))]
      #store.matrix[,2] <- sapply(store.matrix[,1],pop.hazard.func.in.sim,age=age,gender=gender,arrival.time=arrival.time, start.year=start.year)
      age.floor <- floor(age)
      
      large.time.floor.vec.plus.age <- floor(round2(large.time.vec+age,n))
      year.after.first.year.change <- floor(round2(start.year+arrival.time+large.time.vec[1],n))
      if(year.after.first.year.change>=2021){
        year.after.first.year.change <- 2020
      }
      
      if(age.floor<110){
        index <- which(pop.data.male$Age==age.floor & pop.data.male$Year==floor(round2(start.year+arrival.time,n)))
      } else {
        index <- which(pop.data.male$Age=="110+" & pop.data.male$Year==floor(round2(start.year+arrival.time,n)))
      }
      
      if(length(sort.time.vec)==1){
        store.matrix[1,2] <- pop.data.male[index,"mx"]
      } else {
        
        if(small.time.floor.vec.plus.age[1]<110){
          index.small <- which(pop.data.male$Age==small.time.floor.vec.plus.age[1] & pop.data.male$Year==year.after.first.age.change)
        } else {
          index.small <- which(pop.data.male$Age=="110+" & pop.data.male$Year==year.after.first.age.change)
        }      
        if(large.time.floor.vec.plus.age[1]<110){
          index.large <- which(pop.data.male$Age==large.time.floor.vec.plus.age[1] & pop.data.male$Year==year.after.first.year.change)
        } else {
          index.large <- which(pop.data.male$Age=="110+" & pop.data.male$Year==year.after.first.year.change)
        }            
        index.age.over.110.small <- index.small+(sum(small.time.floor.vec.plus.age<110))*112
        index.age.over.110.large <- index.large+(sum(large.time.floor.vec.plus.age<110))*112
        
        num.small.time.vec.age.less.110 <- sum(small.time.floor.vec.plus.age<=110)
        num.large.time.vec.age.less.110 <- sum(large.time.floor.vec.plus.age<=110)
        num.small.time.vec.age.larger.110 <- sum(small.time.floor.vec.plus.age>110)
        num.large.time.vec.age.larger.110 <- sum(large.time.floor.vec.plus.age>110)
        
        store.matrix[1,2] <- pop.data.male[index,"mx"]
        if(num.small.time.vec.age.less.110!=0){
          if(small.time.step!=0){
            for(i in 1:num.small.time.vec.age.less.110){
              store.matrix[2*i,2] <- pop.data.male[index.small+(i-1)*112,"mx"]
            }
          } else {
            for(i in 1:num.small.time.vec.age.less.110){
              store.matrix[2*i+1,2] <- pop.data.male[index.small+(i-1)*112,"mx"]
            }
          }
        }
        if(length(change.year.time.vec)!=1 & num.large.time.vec.age.less.110!=0){ 
          if(small.time.step!=0){
            for(i in 1:num.large.time.vec.age.less.110){
              store.matrix[2*i+1,2] <- pop.data.male[index.large+(i-1)*112,"mx"]
            }
          } else {
            for(i in 1:num.large.time.vec.age.less.110){
              store.matrix[2*i,2] <- pop.data.male[index.large+(i-1)*112,"mx"]
            }
          }
        }
        
        if(num.small.time.vec.age.larger.110!=0){
          if(num.small.time.vec.age.less.110!=0){
            if(small.time.step!=0){
              for(i in 1:num.small.time.vec.age.larger.110){
                store.matrix[2*(i+num.small.time.vec.age.less.110),2] <- pop.data.male[index.age.over.110.small+(i)*111,"mx"]
              }
            } else {
              for(i in 1:num.small.time.vec.age.larger.110){
                store.matrix[2*(i+num.small.time.vec.age.less.110)+1,2] <- pop.data.male[index.age.over.110.small+(i)*111,"mx"]
              }
            }
          } else {
            if(small.time.step!=0){
              for(i in 1:num.small.time.vec.age.larger.110){
                store.matrix[2*(i+num.small.time.vec.age.less.110),2] <- pop.data.male[index.age.over.110.small+(i-1)*111,"mx"]
              }
            } else {
              for(i in 1:num.small.time.vec.age.larger.110){
                store.matrix[2*(i+num.small.time.vec.age.less.110)+1,2] <- pop.data.male[index.age.over.110.small+(i-1)*111,"mx"]
              }
            }
          }
        }
        
        if(num.large.time.vec.age.larger.110!=0 & length(change.year.time.vec)!=1){
          if(num.large.time.vec.age.less.110!=0){
            if(small.time.step!=0){
              for(i in 1:num.large.time.vec.age.larger.110){
                store.matrix[2*(i+num.large.time.vec.age.less.110)+1,2] <- pop.data.male[index.age.over.110.large+(i)*111,"mx"]
              }
            } else {
              for(i in 1:num.large.time.vec.age.larger.110){
                store.matrix[2*(i+num.large.time.vec.age.less.110),2] <- pop.data.male[index.age.over.110.large+(i)*111,"mx"]
              }
            }
          } else {
            if(small.time.step!=0){
              for(i in 1:num.large.time.vec.age.larger.110){
                store.matrix[2*(i+num.large.time.vec.age.less.110)+1,2] <- pop.data.male[index.age.over.110.large+(i-1)*111,"mx"]
              }
            } else {
              for(i in 1:num.large.time.vec.age.larger.110){
                store.matrix[2*(i+num.large.time.vec.age.less.110),2] <- pop.data.male[index.age.over.110.large+(i-1)*111,"mx"]
              }
            }
          }
        }
      }
    } else {
      small.time.step <- step.year
      small.time.vec <- change.year.time.vec[1:(length(change.year.time.vec)-1)]
      large.time.step <- step.age
      large.time.vec <- change.age.time.vec
      
      if(small.time.step!=0){
        sort.time.vec <- sort(c(0,small.time.vec,large.time.vec))
        small.time.floor.vec.plus.age <- floor(round2(small.time.vec+age,n))
        year.after.first.year.change <- floor(round2(start.year+arrival.time+small.time.vec[1],n))
      } else {
        sort.time.vec <- sort(c(small.time.vec,large.time.vec))
        small.time.floor.vec.plus.age <- floor(round2(small.time.vec+age,n))[-1]
        year.after.first.year.change <- floor(round2(start.year+arrival.time+small.time.vec[2],n))
      }
      
      if(length(change.year.time.vec)==1){
        sort.time.vec <- sort.time.vec[-((length(sort.time.vec)-1):length(sort.time.vec))]
      }
      
      store.matrix <- matrix(NA,nrow=length(sort.time.vec),ncol=2)
      store.matrix[,1] <- sort.time.vec[1:(length(sort.time.vec))]
      age.floor <- floor(age)
      large.time.floor.vec.plus.age <- floor(round2(large.time.vec+age,n))
      year.after.first.age.change <- floor(round2(start.year+arrival.time+large.time.vec[1],n))
      
      
      
      if(age.floor<110){
        index <- which(pop.data.male$Age==age.floor & pop.data.male$Year==floor(round2(start.year+arrival.time,n)))
      } else {
        index <- which(pop.data.male$Age=="110+" & pop.data.male$Year==floor(round2(start.year+arrival.time,n)))
      }
      
      if(length(sort.time.vec)==1){
        store.matrix[1,2] <- pop.data.male[index,"mx"]
      } else {
        
        # if(floor(round2(small.time.vec[1]+age,10))<110){
        #   index.small <- which(pop.data.male$Age==floor(round2(small.time.vec[1]+age,10)) & pop.data.male$Year==floor(round2(start.year+arrival.time+small.time.vec[1],10)))
        # } else {
        #   index.small <- which(pop.data.male$Age=="110+" & pop.data.male$Year==floor(round2(start.year+arrival.time+small.time.vec[1],10)))
        # }
        # 
        # if(floor(round2(large.time.vec[1]+age,10))<110){
        #   index.large <- which(pop.data.male$Age==floor(round2(large.time.vec[1]+age,10)) & pop.data.male$Year==floor(round2(start.year+arrival.time+large.time.vec[1],10)))
        # } else {
        #   index.large <- which(pop.data.male$Age=="110+" & pop.data.male$Year==floor(round2(start.year+arrival.time+large.time.vec[1],10)))
        # }      
        
        if(length(small.time.vec)!=1){
          if(small.time.floor.vec.plus.age[1]<110){
            index.small <- which(pop.data.male$Age==small.time.floor.vec.plus.age[1] & pop.data.male$Year==year.after.first.year.change)
          } else {
            index.small <- which(pop.data.male$Age=="110+" & pop.data.male$Year==year.after.first.year.change)
          }      
          if(large.time.floor.vec.plus.age[1]<110){
            index.large <- which(pop.data.male$Age==large.time.floor.vec.plus.age[1] & pop.data.male$Year==year.after.first.age.change)
          } else {
            index.large <- which(pop.data.male$Age=="110+" & pop.data.male$Year==year.after.first.age.change)
          }   
          
          if(small.time.step!=0){
            index.age.over.110.small <- index.small+(sum(floor(round2(small.time.vec+age,n))<110))*112
          } else {
            index.age.over.110.small <- index.small+(sum(floor(round2(small.time.vec+age,n))[-1]<110))*112
          }
          
          index.age.over.110.large <- index.large+(sum(floor(round2(large.time.vec+age,n))<110))*112
          
          num.small.time.vec.age.less.110 <- sum(floor(round2(small.time.vec+age,n))<=110)
          num.large.time.vec.age.less.110 <- sum(floor(round2(large.time.vec+age,n))<=110)
          num.small.time.vec.age.larger.110 <- sum(floor(round2(small.time.vec+age,n))>110)
          num.large.time.vec.age.larger.110 <- sum(floor(round2(large.time.vec+age,n))>110)
          
          store.matrix[1,2] <- pop.data.male[index,"mx"]
          # if(num.small.time.vec.age.less.110!=0){
          #   for(i in 1:num.small.time.vec.age.less.110){
          #     store.matrix[2*i,2] <- pop.data.male[index.small+(i-1)*112,"mx"]
          #   }
          # }
          # 
          # if(num.large.time.vec.age.less.110!=0){
          #   for(i in 1:num.large.time.vec.age.less.110){
          #     store.matrix[2*i+1,2] <- pop.data.male[index.large+(i-1)*112,"mx"]
          #   }
          # }  
          
          # if(num.small.time.vec.age.larger.110!=0){ # FIX THIS
          #   for(i in 1:num.small.time.vec.age.larger.110){
          #     store.matrix[2*(i+num.small.time.vec.age.less.110),2] <- pop.data.male[index.age.over.110.small+(i)*111,"mx"]
          #   }
          # }
          
          if(num.small.time.vec.age.less.110!=0 & age.floor!=110){
            if(small.time.step!=0){
              for(i in 1:num.small.time.vec.age.less.110){
                store.matrix[2*i,2] <- pop.data.male[index.small+(i-1)*112,"mx"]
              }
            } else {
              for(i in 1:(num.small.time.vec.age.less.110-1)){
                store.matrix[2*i+1,2] <- pop.data.male[index.small+(i-1)*112,"mx"]
              }
            }
          }
          if(length(change.year.time.vec)!=1 & num.large.time.vec.age.less.110!=0){ 
            if(small.time.step!=0){
              for(i in 1:num.large.time.vec.age.less.110){
                store.matrix[2*i+1,2] <- pop.data.male[index.large+(i-1)*112,"mx"]
              }
            } else {
              for(i in 1:num.large.time.vec.age.less.110){
                store.matrix[2*i,2] <- pop.data.male[index.large+(i-1)*112,"mx"]
              }
            }
          }
          
          if(num.small.time.vec.age.larger.110!=0){ # TEST 
            if(num.small.time.vec.age.less.110!=0 & age.floor!=110){
              if(small.time.step!=0){
                for(i in 1:num.small.time.vec.age.larger.110){
                  store.matrix[2*(i+num.small.time.vec.age.less.110),2] <- pop.data.male[index.age.over.110.small+(i)*111,"mx"]
                }
              } else {
                for(i in 1:num.small.time.vec.age.larger.110){
                  store.matrix[2*(i+num.small.time.vec.age.less.110-1)+1,2] <- pop.data.male[index.age.over.110.small+(i)*111,"mx"]
                }
              } 
            } else if(num.small.time.vec.age.less.110==0 & age.floor<110) {
              if(small.time.step!=0){
                for(i in 1:num.small.time.vec.age.larger.110){
                  store.matrix[2*(i+num.small.time.vec.age.less.110),2] <- pop.data.male[index.age.over.110.small+(i-1)*111,"mx"]
                }
              } else {
                for(i in 1:num.small.time.vec.age.larger.110){
                  store.matrix[2*(i+num.small.time.vec.age.less.110-1)+1,2] <- pop.data.male[index.age.over.110.small+(i-1)*111,"mx"] # fix here
                }
                
              }
            } else if (num.small.time.vec.age.less.110==0 & age.floor>110) {
              if(small.time.step!=0){
                for(i in 1:(num.small.time.vec.age.larger.110)){
                  store.matrix[2*(i+num.small.time.vec.age.less.110),2] <- pop.data.male[index.age.over.110.small+(i-1)*111,"mx"]
                }
              } else {
                for(i in 1:(num.small.time.vec.age.larger.110-1)){
                  store.matrix[2*(i+num.small.time.vec.age.less.110)+1,2] <- pop.data.male[index.age.over.110.small+(i-1)*111,"mx"]
                }
              }
            } else {
              if(small.time.step!=0){
                for(i in 1:(num.small.time.vec.age.larger.110+1)){
                  store.matrix[2*(i+num.small.time.vec.age.less.110-1),2] <- pop.data.male[index.age.over.110.small+(i-1)*111,"mx"]
                }
              } else {
                for(i in 1:(num.small.time.vec.age.larger.110)){
                  store.matrix[2*(i+num.small.time.vec.age.less.110-1)+1,2] <- pop.data.male[index.age.over.110.small+(i-1)*111,"mx"]
                }
              }
            }
          }
          
          # if(num.large.time.vec.age.larger.110!=0){
          #   for(i in 1:num.large.time.vec.age.larger.110){
          #     store.matrix[2*(i+num.large.time.vec.age.less.110)+1,2] <- pop.data.male[index.age.over.110.large+(i)*111,"mx"]
          #   }
          # }
          
          if(num.large.time.vec.age.larger.110!=0 & length(change.year.time.vec)!=1){
            if(num.large.time.vec.age.less.110!=0){
              if(small.time.step!=0){
                for(i in 1:num.large.time.vec.age.larger.110){
                  store.matrix[2*(i+num.large.time.vec.age.less.110)+1,2] <- pop.data.male[index.age.over.110.large+(i)*111,"mx"]
                }
              } else {
                for(i in 1:num.large.time.vec.age.larger.110){
                  store.matrix[2*(i+num.large.time.vec.age.less.110),2] <- pop.data.male[index.age.over.110.large+(i)*111,"mx"]
                }
              }
            } else {
              if(small.time.step!=0){
                for(i in 1:num.large.time.vec.age.larger.110){
                  store.matrix[2*(i+num.large.time.vec.age.less.110)+1,2] <- pop.data.male[index.age.over.110.large+(i-1)*111,"mx"]
                }
              } else {
                for(i in 1:num.large.time.vec.age.larger.110){
                  store.matrix[2*(i+num.large.time.vec.age.less.110),2] <- pop.data.male[index.age.over.110.large+(i-1)*111,"mx"]
                }
              }
            }
          }
        } else {
          if(small.time.step==0){
            store.matrix[1,2] <- pop.data.male[index,"mx"]
            if(large.time.floor.vec.plus.age[1]<110){
              index.large <- which(pop.data.male$Age==large.time.floor.vec.plus.age[1] & pop.data.male$Year==year.after.first.age.change)
            } else {
              index.large <- which(pop.data.male$Age=="110+" & pop.data.male$Year==year.after.first.age.change)
            }
            store.matrix[2,2] <- pop.data.male[index.large,"mx"]
          } else {
            store.matrix[1,2] <- pop.data.male[index,"mx"]
            if(small.time.floor.vec.plus.age[1]<110){
              index.small <- which(pop.data.male$Age==small.time.floor.vec.plus.age[1] & pop.data.male$Year==year.after.first.year.change)
            } else {
              index.small <- which(pop.data.male$Age=="110+" & pop.data.male$Year==year.after.first.age.change)
            }      
            if(large.time.floor.vec.plus.age[1]<110){
              index.large <- which(pop.data.male$Age==large.time.floor.vec.plus.age[1] & pop.data.male$Year==year.after.first.age.change)
            } else {
              index.large <- which(pop.data.male$Age=="110+" & pop.data.male$Year==year.after.first.age.change)
            } 
            store.matrix[2,2] <- pop.data.male[index.small,"mx"]
            store.matrix[3,2] <- pop.data.male[index.large,"mx"]
            
          }
        }
      }
    }
    return(store.matrix)
  }
}

haz.output.female <- function(age,start.year,end.year,arrival.time,n=8){ # n - "precision" --> can get error if smaller n
  
  age <- round2(age,n)
  arrival.time <- round2(arrival.time,n)
  
  step.age <- round2(ceiling(age)-age,n)
  step.year <- round2(ceiling(start.year+arrival.time)-(start.year+arrival.time),n)
  max.follow.up <- round2((end.year-start.year-arrival.time),n)
  
  if(max.follow.up< round2(step.age+1,n)){
    change.age.time.vec <- c(step.age)
  } else {
    change.age.time.vec <- c(step.age,(seq(from=round2(step.age+1,n),to=max.follow.up,by=1)))
  }
  
  if(max.follow.up< round2(step.year+1,n)){
    change.year.time.vec <- c(step.year)
  } else {
    change.year.time.vec <- c(step.year,(seq(from=round2(step.year+1,n),to=max.follow.up,by=1)))
  }
  
  ### If step.age is equal to step.year ###
  if(step.age==step.year){
    if(step.age==0){
      store.matrix.equal <- matrix(NA,ncol=2,nrow = length(change.age.time.vec)-1)
      store.matrix.equal[,1] <- change.age.time.vec[1:(length(change.age.time.vec)-1)]
      age.throughout.vec <- floor(round2(store.matrix.equal[,1]+age,n))
      if(age.throughout.vec[1]<110){
        index <- which(pop.data.female$Age==age.throughout.vec[1] & pop.data.female$Year==floor(round2(start.year+arrival.time,n)))
      } else {
        index <- which(pop.data.female$Age=="110+" & pop.data.female$Year==floor(round2(start.year+arrival.time,n)))
      }
      index.age.over.110 <- index+(sum(age.throughout.vec<110))*112
      if(sum(age.throughout.vec<=110)!=0){
        for(i in 1:sum(age.throughout.vec<=110)){
          store.matrix.equal[i,2] <- pop.data.female[index+(i-1)*112,"mx"]
        }
      }
      if(sum(age.throughout.vec>110)!=0){
        if(sum(age.throughout.vec<=110)!=0){
          for(i in 1:sum(age.throughout.vec>110)){
            store.matrix.equal[i+sum(age.throughout.vec<=110),2] <- pop.data.female[index.age.over.110+(i)*111,"mx"]
          }
        } else {
          for(i in 1:sum(age.throughout.vec>110)){
            store.matrix.equal[i+sum(age.throughout.vec<=110),2] <- pop.data.female[index.age.over.110+(i-1)*111,"mx"]
          }
        }
      }
    } else {
      if(length(change.age.time.vec)==1){
        store.matrix.equal <- matrix(NA,ncol=2,nrow=1)
        store.matrix.equal[1,1] <- 0
        age.floor <- floor(age)
        if(age.floor<110){
          index <- which(pop.data.female$Age==age.floor & pop.data.female$Year==floor(round2(start.year+arrival.time,n)))
        } else {
          index <- which(pop.data.female$Age=="110+" & pop.data.female$Year==floor(round2(start.year+arrival.time,n)))
        }
        store.matrix.equal[1,2] <- pop.data.female[index,"mx"]
      } else {
        store.matrix.equal <- matrix(NA,ncol=2,nrow=length(change.age.time.vec))
        store.matrix.equal[,1] <- c(0,change.age.time.vec[1:(length(change.age.time.vec)-1)])
        #store.matrix.equal[,2] <- vec.pop.hazard.func.in.sim(age,gender,arrival.time,start.year,store.matrix.equal[,1])
        age.throughout.vec <- floor(round2(store.matrix.equal[,1]+age,n))
        if(age.throughout.vec[1]<110){
          index <- which(pop.data.female$Age==age.throughout.vec[1] & pop.data.female$Year==floor(round2(start.year+arrival.time,n)))
        } else {
          index <- which(pop.data.female$Age=="110+" & pop.data.female$Year==floor(round2(start.year+arrival.time,n)))
        }
        index.age.over.110 <- index+(sum(age.throughout.vec<110))*112
        if(sum(age.throughout.vec<=110)!=0){
          for(i in 1:sum(age.throughout.vec<=110)){
            store.matrix.equal[i,2] <- pop.data.female[index+(i-1)*112,"mx"]
          }
        }
        if(sum(age.throughout.vec>110)!=0){
          if(sum(age.throughout.vec<=110)!=0){
            for(i in 1:sum(age.throughout.vec>110)){
              store.matrix.equal[i+sum(age.throughout.vec<=110),2] <- pop.data.female[index.age.over.110+(i)*111,"mx"]
            }
          } else {
            for(i in 1:sum(age.throughout.vec>110)){
              store.matrix.equal[i+sum(age.throughout.vec<=110),2] <- pop.data.female[index.age.over.110+(i-1)*111,"mx"]
            }
          }
        }
      }
    }
    return(store.matrix.equal)
  }
  
  ### If step.age is different from step.year ### 
  
  if(step.age!=step.year){
    
    if(step.age < step.year){
      small.time.step <- step.age
      small.time.vec <- change.age.time.vec
      large.time.step <- step.year
      large.time.vec <- change.year.time.vec[1:(length(change.year.time.vec)-1)]
      
      if(small.time.step!=0){
        sort.time.vec <- sort(c(0,small.time.vec,large.time.vec))
        small.time.floor.vec.plus.age <- floor(round2(small.time.vec+age,n))
        year.after.first.age.change <- floor(round2(start.year+arrival.time+small.time.vec[1],n))
      } else {
        sort.time.vec <- sort(c(small.time.vec,large.time.vec))
        small.time.floor.vec.plus.age <- floor(round2(small.time.vec+age,n))[-1]
        year.after.first.age.change <- floor(round2(start.year+arrival.time+small.time.vec[2],n))
      }
      
      if(length(change.year.time.vec)==1){
        sort.time.vec <- sort.time.vec[-length(sort.time.vec)]
      }
      
      store.matrix <- matrix(NA,nrow=length(sort.time.vec),ncol=2)
      store.matrix[,1] <- sort.time.vec[1:(length(sort.time.vec))]
      #store.matrix[,2] <- sapply(store.matrix[,1],pop.hazard.func.in.sim,age=age,gender=gender,arrival.time=arrival.time, start.year=start.year)
      age.floor <- floor(age)
      
      large.time.floor.vec.plus.age <- floor(round2(large.time.vec+age,n))
      year.after.first.year.change <- floor(round2(start.year+arrival.time+large.time.vec[1],n))
      if(year.after.first.year.change>=2021){
        year.after.first.year.change <- 2020
      }
      
      if(age.floor<110){
        index <- which(pop.data.female$Age==age.floor & pop.data.female$Year==floor(round2(start.year+arrival.time,n)))
      } else {
        index <- which(pop.data.female$Age=="110+" & pop.data.female$Year==floor(round2(start.year+arrival.time,n)))
      }
      
      if(length(sort.time.vec)==1){
        store.matrix[1,2] <- pop.data.female[index,"mx"]
      } else {
        
        if(small.time.floor.vec.plus.age[1]<110){
          index.small <- which(pop.data.female$Age==small.time.floor.vec.plus.age[1] & pop.data.female$Year==year.after.first.age.change)
        } else {
          index.small <- which(pop.data.female$Age=="110+" & pop.data.female$Year==year.after.first.age.change)
        }      
        if(large.time.floor.vec.plus.age[1]<110){
          index.large <- which(pop.data.female$Age==large.time.floor.vec.plus.age[1] & pop.data.female$Year==year.after.first.year.change)
        } else {
          index.large <- which(pop.data.female$Age=="110+" & pop.data.female$Year==year.after.first.year.change)
        }            
        index.age.over.110.small <- index.small+(sum(small.time.floor.vec.plus.age<110))*112
        index.age.over.110.large <- index.large+(sum(large.time.floor.vec.plus.age<110))*112
        
        num.small.time.vec.age.less.110 <- sum(small.time.floor.vec.plus.age<=110)
        num.large.time.vec.age.less.110 <- sum(large.time.floor.vec.plus.age<=110)
        num.small.time.vec.age.larger.110 <- sum(small.time.floor.vec.plus.age>110)
        num.large.time.vec.age.larger.110 <- sum(large.time.floor.vec.plus.age>110)
        
        store.matrix[1,2] <- pop.data.female[index,"mx"]
        if(num.small.time.vec.age.less.110!=0){
          if(small.time.step!=0){
            for(i in 1:num.small.time.vec.age.less.110){
              store.matrix[2*i,2] <- pop.data.female[index.small+(i-1)*112,"mx"]
            }
          } else {
            for(i in 1:num.small.time.vec.age.less.110){
              store.matrix[2*i+1,2] <- pop.data.female[index.small+(i-1)*112,"mx"]
            }
          }
        }
        if(length(change.year.time.vec)!=1 & num.large.time.vec.age.less.110!=0){ 
          if(small.time.step!=0){
            for(i in 1:num.large.time.vec.age.less.110){
              store.matrix[2*i+1,2] <- pop.data.female[index.large+(i-1)*112,"mx"]
            }
          } else {
            for(i in 1:num.large.time.vec.age.less.110){
              store.matrix[2*i,2] <- pop.data.female[index.large+(i-1)*112,"mx"]
            }
          }
        }
        
        if(num.small.time.vec.age.larger.110!=0){
          if(num.small.time.vec.age.less.110!=0){
            if(small.time.step!=0){
              for(i in 1:num.small.time.vec.age.larger.110){
                store.matrix[2*(i+num.small.time.vec.age.less.110),2] <- pop.data.female[index.age.over.110.small+(i)*111,"mx"]
              }
            } else {
              for(i in 1:num.small.time.vec.age.larger.110){
                store.matrix[2*(i+num.small.time.vec.age.less.110)+1,2] <- pop.data.female[index.age.over.110.small+(i)*111,"mx"]
              }
            }
          } else {
            if(small.time.step!=0){
              for(i in 1:num.small.time.vec.age.larger.110){
                store.matrix[2*(i+num.small.time.vec.age.less.110),2] <- pop.data.female[index.age.over.110.small+(i-1)*111,"mx"]
              }
            } else {
              for(i in 1:num.small.time.vec.age.larger.110){
                store.matrix[2*(i+num.small.time.vec.age.less.110)+1,2] <- pop.data.female[index.age.over.110.small+(i-1)*111,"mx"]
              }
            }
          }
        }
        
        if(num.large.time.vec.age.larger.110!=0 & length(change.year.time.vec)!=1){
          if(num.large.time.vec.age.less.110!=0){
            if(small.time.step!=0){
              for(i in 1:num.large.time.vec.age.larger.110){
                store.matrix[2*(i+num.large.time.vec.age.less.110)+1,2] <- pop.data.female[index.age.over.110.large+(i)*111,"mx"]
              }
            } else {
              for(i in 1:num.large.time.vec.age.larger.110){
                store.matrix[2*(i+num.large.time.vec.age.less.110),2] <- pop.data.female[index.age.over.110.large+(i)*111,"mx"]
              }
            }
          } else {
            if(small.time.step!=0){
              for(i in 1:num.large.time.vec.age.larger.110){
                store.matrix[2*(i+num.large.time.vec.age.less.110)+1,2] <- pop.data.female[index.age.over.110.large+(i-1)*111,"mx"]
              }
            } else {
              for(i in 1:num.large.time.vec.age.larger.110){
                store.matrix[2*(i+num.large.time.vec.age.less.110),2] <- pop.data.female[index.age.over.110.large+(i-1)*111,"mx"]
              }
            }
          }
        }
      }
    } else {
      small.time.step <- step.year
      small.time.vec <- change.year.time.vec[1:(length(change.year.time.vec)-1)]
      large.time.step <- step.age
      large.time.vec <- change.age.time.vec
      
      if(small.time.step!=0){
        sort.time.vec <- sort(c(0,small.time.vec,large.time.vec))
        small.time.floor.vec.plus.age <- floor(round2(small.time.vec+age,n))
        year.after.first.year.change <- floor(round2(start.year+arrival.time+small.time.vec[1],n))
      } else {
        sort.time.vec <- sort(c(small.time.vec,large.time.vec))
        small.time.floor.vec.plus.age <- floor(round2(small.time.vec+age,n))[-1]
        year.after.first.year.change <- floor(round2(start.year+arrival.time+small.time.vec[2],n))
      }
      
      if(length(change.year.time.vec)==1){
        sort.time.vec <- sort.time.vec[-((length(sort.time.vec)-1):length(sort.time.vec))]
      }
      
      store.matrix <- matrix(NA,nrow=length(sort.time.vec),ncol=2)
      store.matrix[,1] <- sort.time.vec[1:(length(sort.time.vec))]
      age.floor <- floor(age)
      large.time.floor.vec.plus.age <- floor(round2(large.time.vec+age,n))
      year.after.first.age.change <- floor(round2(start.year+arrival.time+large.time.vec[1],n))
      
      
      
      if(age.floor<110){
        index <- which(pop.data.female$Age==age.floor & pop.data.female$Year==floor(round2(start.year+arrival.time,n)))
      } else {
        index <- which(pop.data.female$Age=="110+" & pop.data.female$Year==floor(round2(start.year+arrival.time,n)))
      }
      
      if(length(sort.time.vec)==1){
        store.matrix[1,2] <- pop.data.female[index,"mx"]
      } else {
        
        # if(floor(round2(small.time.vec[1]+age,10))<110){
        #   index.small <- which(pop.data.female$Age==floor(round2(small.time.vec[1]+age,10)) & pop.data.female$Year==floor(round2(start.year+arrival.time+small.time.vec[1],10)))
        # } else {
        #   index.small <- which(pop.data.female$Age=="110+" & pop.data.female$Year==floor(round2(start.year+arrival.time+small.time.vec[1],10)))
        # }
        # 
        # if(floor(round2(large.time.vec[1]+age,10))<110){
        #   index.large <- which(pop.data.female$Age==floor(round2(large.time.vec[1]+age,10)) & pop.data.female$Year==floor(round2(start.year+arrival.time+large.time.vec[1],10)))
        # } else {
        #   index.large <- which(pop.data.female$Age=="110+" & pop.data.female$Year==floor(round2(start.year+arrival.time+large.time.vec[1],10)))
        # }      
        
        if(length(small.time.vec)!=1){
          if(small.time.floor.vec.plus.age[1]<110){
            index.small <- which(pop.data.female$Age==small.time.floor.vec.plus.age[1] & pop.data.female$Year==year.after.first.year.change)
          } else {
            index.small <- which(pop.data.female$Age=="110+" & pop.data.female$Year==year.after.first.year.change)
          }      
          if(large.time.floor.vec.plus.age[1]<110){
            index.large <- which(pop.data.female$Age==large.time.floor.vec.plus.age[1] & pop.data.female$Year==year.after.first.age.change)
          } else {
            index.large <- which(pop.data.female$Age=="110+" & pop.data.female$Year==year.after.first.age.change)
          }   
          
          if(small.time.step!=0){
            index.age.over.110.small <- index.small+(sum(floor(round2(small.time.vec+age,n))<110))*112
          } else {
            index.age.over.110.small <- index.small+(sum(floor(round2(small.time.vec+age,n))[-1]<110))*112
          }
          
          index.age.over.110.large <- index.large+(sum(floor(round2(large.time.vec+age,n))<110))*112
          
          num.small.time.vec.age.less.110 <- sum(floor(round2(small.time.vec+age,n))<=110)
          num.large.time.vec.age.less.110 <- sum(floor(round2(large.time.vec+age,n))<=110)
          num.small.time.vec.age.larger.110 <- sum(floor(round2(small.time.vec+age,n))>110)
          num.large.time.vec.age.larger.110 <- sum(floor(round2(large.time.vec+age,n))>110)
          
          store.matrix[1,2] <- pop.data.female[index,"mx"]
          # if(num.small.time.vec.age.less.110!=0){
          #   for(i in 1:num.small.time.vec.age.less.110){
          #     store.matrix[2*i,2] <- pop.data.female[index.small+(i-1)*112,"mx"]
          #   }
          # }
          # 
          # if(num.large.time.vec.age.less.110!=0){
          #   for(i in 1:num.large.time.vec.age.less.110){
          #     store.matrix[2*i+1,2] <- pop.data.female[index.large+(i-1)*112,"mx"]
          #   }
          # }  
          
          # if(num.small.time.vec.age.larger.110!=0){ # FIX THIS
          #   for(i in 1:num.small.time.vec.age.larger.110){
          #     store.matrix[2*(i+num.small.time.vec.age.less.110),2] <- pop.data.female[index.age.over.110.small+(i)*111,"mx"]
          #   }
          # }
          
          if(num.small.time.vec.age.less.110!=0 & age.floor!=110){
            if(small.time.step!=0){
              for(i in 1:num.small.time.vec.age.less.110){
                store.matrix[2*i,2] <- pop.data.female[index.small+(i-1)*112,"mx"]
              }
            } else {
              for(i in 1:(num.small.time.vec.age.less.110-1)){
                store.matrix[2*i+1,2] <- pop.data.female[index.small+(i-1)*112,"mx"]
              }
            }
          }
          if(length(change.year.time.vec)!=1 & num.large.time.vec.age.less.110!=0){ 
            if(small.time.step!=0){
              for(i in 1:num.large.time.vec.age.less.110){
                store.matrix[2*i+1,2] <- pop.data.female[index.large+(i-1)*112,"mx"]
              }
            } else {
              for(i in 1:num.large.time.vec.age.less.110){
                store.matrix[2*i,2] <- pop.data.female[index.large+(i-1)*112,"mx"]
              }
            }
          }
          
          if(num.small.time.vec.age.larger.110!=0){ # TEST 
            if(num.small.time.vec.age.less.110!=0 & age.floor!=110){
              if(small.time.step!=0){
                for(i in 1:num.small.time.vec.age.larger.110){
                  store.matrix[2*(i+num.small.time.vec.age.less.110),2] <- pop.data.female[index.age.over.110.small+(i)*111,"mx"]
                }
              } else {
                for(i in 1:num.small.time.vec.age.larger.110){
                  store.matrix[2*(i+num.small.time.vec.age.less.110-1)+1,2] <- pop.data.female[index.age.over.110.small+(i)*111,"mx"]
                }
              } 
            } else if(num.small.time.vec.age.less.110==0 & age.floor<110) {
              if(small.time.step!=0){
                for(i in 1:num.small.time.vec.age.larger.110){
                  store.matrix[2*(i+num.small.time.vec.age.less.110),2] <- pop.data.female[index.age.over.110.small+(i-1)*111,"mx"]
                }
              } else {
                for(i in 1:num.small.time.vec.age.larger.110){
                  store.matrix[2*(i+num.small.time.vec.age.less.110-1)+1,2] <- pop.data.female[index.age.over.110.small+(i-1)*111,"mx"] # fix here
                }
                
              }
            } else if (num.small.time.vec.age.less.110==0 & age.floor>110) {
              if(small.time.step!=0){
                for(i in 1:(num.small.time.vec.age.larger.110)){
                  store.matrix[2*(i+num.small.time.vec.age.less.110),2] <- pop.data.female[index.age.over.110.small+(i-1)*111,"mx"]
                }
              } else {
                for(i in 1:(num.small.time.vec.age.larger.110-1)){
                  store.matrix[2*(i+num.small.time.vec.age.less.110)+1,2] <- pop.data.female[index.age.over.110.small+(i-1)*111,"mx"]
                }
              }
            } else {
              if(small.time.step!=0){
                for(i in 1:(num.small.time.vec.age.larger.110+1)){
                  store.matrix[2*(i+num.small.time.vec.age.less.110-1),2] <- pop.data.female[index.age.over.110.small+(i-1)*111,"mx"]
                }
              } else {
                for(i in 1:(num.small.time.vec.age.larger.110)){
                  store.matrix[2*(i+num.small.time.vec.age.less.110-1)+1,2] <- pop.data.female[index.age.over.110.small+(i-1)*111,"mx"]
                }
              }
            }
          }
          
          # if(num.large.time.vec.age.larger.110!=0){
          #   for(i in 1:num.large.time.vec.age.larger.110){
          #     store.matrix[2*(i+num.large.time.vec.age.less.110)+1,2] <- pop.data.female[index.age.over.110.large+(i)*111,"mx"]
          #   }
          # }
          
          if(num.large.time.vec.age.larger.110!=0 & length(change.year.time.vec)!=1){
            if(num.large.time.vec.age.less.110!=0){
              if(small.time.step!=0){
                for(i in 1:num.large.time.vec.age.larger.110){
                  store.matrix[2*(i+num.large.time.vec.age.less.110)+1,2] <- pop.data.female[index.age.over.110.large+(i)*111,"mx"]
                }
              } else {
                for(i in 1:num.large.time.vec.age.larger.110){
                  store.matrix[2*(i+num.large.time.vec.age.less.110),2] <- pop.data.female[index.age.over.110.large+(i)*111,"mx"]
                }
              }
            } else {
              if(small.time.step!=0){
                for(i in 1:num.large.time.vec.age.larger.110){
                  store.matrix[2*(i+num.large.time.vec.age.less.110)+1,2] <- pop.data.female[index.age.over.110.large+(i-1)*111,"mx"]
                }
              } else {
                for(i in 1:num.large.time.vec.age.larger.110){
                  store.matrix[2*(i+num.large.time.vec.age.less.110),2] <- pop.data.female[index.age.over.110.large+(i-1)*111,"mx"]
                }
              }
            }
          }
        } else {
          if(small.time.step==0){
            store.matrix[1,2] <- pop.data.female[index,"mx"]
            if(large.time.floor.vec.plus.age[1]<110){
              index.large <- which(pop.data.female$Age==large.time.floor.vec.plus.age[1] & pop.data.female$Year==year.after.first.age.change)
            } else {
              index.large <- which(pop.data.female$Age=="110+" & pop.data.female$Year==year.after.first.age.change)
            }
            store.matrix[2,2] <- pop.data.female[index.large,"mx"]
          } else {
            store.matrix[1,2] <- pop.data.female[index,"mx"]
            if(small.time.floor.vec.plus.age[1]<110){
              index.small <- which(pop.data.female$Age==small.time.floor.vec.plus.age[1] & pop.data.female$Year==year.after.first.year.change)
            } else {
              index.small <- which(pop.data.female$Age=="110+" & pop.data.female$Year==year.after.first.age.change)
            }      
            if(large.time.floor.vec.plus.age[1]<110){
              index.large <- which(pop.data.female$Age==large.time.floor.vec.plus.age[1] & pop.data.female$Year==year.after.first.age.change)
            } else {
              index.large <- which(pop.data.female$Age=="110+" & pop.data.female$Year==year.after.first.age.change)
            } 
            store.matrix[2,2] <- pop.data.female[index.small,"mx"]
            store.matrix[3,2] <- pop.data.female[index.large,"mx"]
            
          }
        }
      }
    }
    return(store.matrix)
  }
}

one.sim.test <- function(age,gender,start.year,end.year,arrival.time){
  if(gender==0){
    haz.matrix <- haz.output.male(age,start.year,end.year,arrival.time)
  } else {
    haz.matrix <- haz.output.female(age,start.year,end.year,arrival.time)
  }
  
  if(nrow(haz.matrix)==1){
    cum.hazard.func <- function(t){
      cum.hazard.old <- t*haz.matrix[1,2]
      return(cum.hazard.old)
    }
  } else {
    
    cum.hazard.func <- function(t){
      cum.hazard.old <- 0
      
      for(i in 1:(nrow(haz.matrix)-1)){
        if(haz.matrix[i,1] <= t & t<haz.matrix[i+1,1]){
          cum.hazard.new <- cum.hazard.old + (t-haz.matrix[i,1])*haz.matrix[i,2]
          cum.hazard.old <- cum.hazard.new
          break
        } else {
          cum.hazard.old <- cum.hazard.old + (haz.matrix[i+1,1]-haz.matrix[i,1])*haz.matrix[i,2]
        }
      }
      if(t >= haz.matrix[nrow(haz.matrix),1]){
        cum.hazard.old <- cum.hazard.old+(t-haz.matrix[nrow(haz.matrix),1])*haz.matrix[nrow(haz.matrix),2]
      }
      return(cum.hazard.old)
    }
  }
  
  u<-runif(1)
  eq.solve <- function(t){
    1-exp(-cum.hazard.func(t))-u
  }
  time.sim <- uniroot(eq.solve,lower=0,upper=nrow(haz.matrix)+10000000)$root
  if(time.sim >= end.year-arrival.time-start.year){
    t <- end.year-arrival.time-start.year
    delta.p <- 0
  } else {
    t<- time.sim
    delta.p <- 1
  }
  return(c(t,delta.p))
}

vec.pop.sim.test <- function(age,gender,start.year,end.year,arrival.time){
  time.mat <- matrix(NA,ncol=2,nrow=length(age))
  colnames(time.mat) <- c("t.p","delta.p")
  for (i in 1:length(age)){
    time.mat[i,] <- one.sim.test(age[i],gender[i],start.year,end.year,arrival.time[i])
  }
  return(time.mat)
}

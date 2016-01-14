## Graph
# List of neighbors
Neighbors=vector("list",9)
Neighbors[[1]]=c(6)
Neighbors[[2]]=c(6)
Neighbors[[3]]=c(7)
Neighbors[[4]]=c(7)
Neighbors[[5]]=c(9)
Neighbors[[6]]=c(1,2,8)
Neighbors[[7]]=c(3,4,8)
Neighbors[[8]]=c(6,7,9)
Neighbors[[9]]=c(5,8)
Neighbors
# Evidence nodes
evidence=c(1,2,3,4,5)
f=9 #arbitrary root

### Definition of Potentials
## Dual cliques
def_M = function (m){
  # M1
  M1=matrix(c(m[1],m[4],m[2],m[3],
              m[4],m[1],m[3],m[2],
              m[2],m[3],m[1],m[4],
              m[3],m[2],m[4],m[1]),4)

  # M2
  M2=matrix(c(m[1],3*m[4],4*m[2],3*m[3],
              3*m[4],m[1],3*m[3],4*m[2],
              4*m[2],3*m[3],m[1],3*m[4],
              3*m[3],4*m[2],3*m[4],m[1]),4)
  # Potentials
  M = array(dim=c(9,9,4,4))
  M[1,6,,]=M1
  M[2,6,,]=M1
  M[3,7,,]=M1
  M[4,7,,]=M1
  M[6,8,,]=M1
  M[7,8,,]=M1
  M[8,9,,]=M1
  M[5,9,,]=M2
  
  M[6,1,,]=M1
  M[6,2,,]=M1
  M[7,3,,]=M1
  M[7,4,,]=M1
  M[8,6,,]=M1
  M[8,7,,]=M1
  M[9,8,,]=M1
  M[9,5,,]=M2
  return(M)
}


# Single Potentials
psi=matrix(rep(1,9*4),ncol=9)
psi[,9]=c(8,9,9,8)



sum_product = function(Neighbors,psiE,evidence,x,M,f){
  psiE=psi
  psiE[,evidence]=x*psi[,evidence]
  # Messages
  mess=array(data=0,dim=c(9,9,4))
  #Marginals
  marginals=vector("list",9)
  
  # Send message
  send_message=function(j,i){
    # Compute product of mkj
    product_mkj=c(1,1,1,1)
    for (k in Neighbors[[j]]){
      if (k != i){
        product_mkj=product_mkj*mess[k,j,]
      }
    }
    # Compute mji
    mji=c(0,0,0,0)
    possible_xj=diag(4) # sum over all values
    for (j0 in 1:4){
      xj=possible_xj[,j0]
      # Message to be sent
      #mji=mji+ ( M[j,i,,] %*% (psiE[,j]) )* product_mkj*xj
      psiE_xj=psiE[j0,j] # real
      mkj_xj=product_mkj[j0] # real
      psi_xi_xj=M[j,i,,j0] # vector
      mji=mji+ psiE_xj*psi_xi_xj*mkj_xj
    }
    mess[j,i,] <<- mji
  }
  
  # Collect
  collect = function (i,j){
    for (k in Neighbors[[j]]){
      if ( k !=i ){
        collect(j,k)
      }
    }
    send_message(j,i)
  }
  
  # Distribute
  distribute = function (i,j){
    send_message(i,j)
    for (k in Neighbors[[j]]){
      if ( k !=i ){
        distribute(j,k)
      }
    }
  }
  
  # Marginal
  compute_marginal = function(i){
    product_mji=1
    for (j in Neighbors[[i]]){
      product_mji=product_mji*mess[j,i,]
    }
    marginals[[i]] <<- psiE[,i]*product_mji
  }
  
  # Algorithm
  for (e in Neighbors[[f]]){
    collect(f,e)
  }
  for (e in Neighbors[[f]]){
    distribute(f,e)
  }
  for (i in 1:9){
    compute_marginal(i)
  }
  return(lapply(marginals,function (x) x/sum(x)))
}

######
#par(mfrow=c(1,2))
############
## Exemple 1
############

m=c(8,3,2,1)
M=def_M(m)
x=matrix(c(0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1,0),ncol=5)
s1=sum_product(Neighbors,psiE,evidence,x,M,f)
barplot(s1[[9]],names.arg =c("A","C","T","G"),ylab="probability",ylim = c(0,.5),
        main= "First example")

###########
# Exemple 2
###########

m=c(7,4,3,2)
M=def_M(m)
x=matrix(c(0,1,0,0,1,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0),ncol=5)
s2=sum_product(Neighbors,psiE,evidence,x,M,f)
barplot(s2[[9]],names.arg =c("A","C","T","G"),ylab="probability",ylim = c(0,0.5),
        main= "Second example")
s1[[9]]
s2[[9]]
##
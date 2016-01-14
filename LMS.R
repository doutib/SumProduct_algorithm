
# Load data
setwd('~/Documents/STAT241')
data=read.table('./hw2.data',col.names = c('x1','x2','y'))
X=as.matrix(data[,1:2])
Y=as.matrix(data[,3])

# solve(A,b) solve b = Ax
theta_star = solve(t(X)%*%X,t(X) %*% Y)
row.names(theta_star)=NULL

# Eigen vectors / values
eigen_values = eigen(t(X)%*%X)$values
eigen_vectors = eigen(t(X)%*%X)$vectors
lambda_max=max(eigen_values)

# Cost function with vector input theta
J = function (theta){
  a= Y - X %*% theta
  return( (t(a)%*%a) [1,1])
}
J_opt=J(theta_star)
J_opt
eigen_vectors
eigen_values

# J when arguments are sequences of numbers
J_seq = function(theta1_seq,theta2_seq){
  theta_seq=matrix(c(theta1_seq,theta2_seq),nrow=2,byrow=T)
  return(apply(theta_seq,2,J))
}

# Plot
plot_J=function(npoints=100,nlevels=20,xlim=c(theta_star[1]-1,theta_star[1]+1),
                ylim=c(theta_star[2]-1,theta_star[2]+1),
                title="",col="black",eigen=T){
  x=seq(xlim[1],xlim[2],length.out = npoints)
  y=seq(ylim[1],ylim[2],length.out = npoints)
  z=outer(x,y,function(x,y) J_seq(x,y))
  contour(x,y,z,main= title,levels= seq(J_opt,J_opt*20,length.out = nlevels),col=col)
  if (eigen)
    abline(-eigen_vectors[1,2]*theta_star[1]+theta_star[2],eigen_vectors[1,2],col="darkgray")
    abline(-eigen_vectors[2,2]*theta_star[1]+theta_star[2],eigen_vectors[2,2],col="darkgray")
  points(theta_star[1],theta_star[2],col="red",cex=2,pch=20)
}
plot_J(title="Cost function J")

## LMS
theta_t = function(t,rho,theta0,Y,X,init_res){
  if (t==0)
    return(list(theta=theta0,res=init_res))
  else
    i=sample(1:length(Y),1)
    yn=Y[i]
    xn=as.matrix(X[i,],nrow=2)
    # Precedent values obtained by induction
    old=theta_t(t-1,rho,theta0,Y,X,init_res)
    old_theta=old$theta
    old_res=old$res
    # Update theta
    new_theta=old_theta+rho*(yn - (t(old_theta) %*% xn)[1,1]) * xn
    # Update res
    new_res=old_res
    new_res[,t+1]=new_theta
    return(list(theta=new_theta,res=new_res))
}

lms=function(t,rho,theta0,Y,X){
  init_res=matrix(nrow=2,ncol=t+1)
  init_res[,1]=theta0
  return(theta_t(t,rho,theta0,Y,X,init_res)$res)
}

lms_plot=function(t,rho,theta0,Y,X,it=1000){
  plot_J(xlim=theta_star[1]+c(-.2,.2),ylim=theta_star[2]+c(-.2,.2),nlevels=80,col="darkgray",eigen=F)
  l=lms(it,rho,theta0,Y,X)
  lines(l[1,],l[2,],type="l",xlab="x",ylab="y")
  points(theta_star[1],theta_star[2],col="red",pch=19)  
  points(l[1,length(l)/2],l[2,length(l)/2],col="green",pch=19)
  legend( "bottomright",legend=c("Theta_LMS","Theta_star"),col=c("green","red"),pch=c(16,16) )
}

theta0=theta_star+.2
rho1=2/max(eigen_values)
rho2=1/(2*max(eigen_values))
rho3=1/(8*max(eigen_values))
set.seed(4)
it=200
par(mfrow=c(1,3),oma = c(3, 3, 7, 3))
lms_plot(t,rho1,theta0,Y,X,it,title="rho1")
lms_plot(t,rho2,theta0,Y,X,it,title="rho2")
lms_plot(t,rho3,theta0,Y,X,it,title="rho3")
title("LMS algorithm, 200 itérations",outer=T,lwd=3,cex=3)
set.seed(15)
it=1000
par(mfrow=c(1,3),oma = c(3, 3, 7, 3))
lms_plot(t,rho1,theta0,Y,X,it,title="rho1")
lms_plot(t,rho2,theta0,Y,X,it,title="rho2")
lms_plot(t,rho3,theta0,Y,X,it,title="rho3")
title("LMS algorithm, 1000 itérations",outer=T,lwd=3,cex=3)

  
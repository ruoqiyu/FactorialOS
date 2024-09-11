factorial_weight<-function(X,z,y,o=1,variance=TRUE,model='interaction',lambda0=NULL,alpha=1,tau=0.5,c=0.2,tol=10^(-7),max.iter=10000){
  #X: basis function that we want to balance
  #z: treatment combination: each column is a treatment +1 or -1
  #o: non-negligible order of interaction effects
  
  if (!(model %in%c('additive','interaction'))){
    stop("model has to be either additive or interaction.")
  } 

  nz=ncol(z)

  
  zfactor=factor(z[,1])
  if (nz>1){
    for (zk in 2:nz) zfactor=zfactor:factor(z[,zk])
  }
  ztable=table(zfactor)
  if (any(ztable==0)) result=factorial_weight_incomplete(X,z,y,o,variance,model,lambda0,alpha,tau,c,tol,max.iter)
  else result=factorial_weight_full(X,z,y,o,variance,model,lambda0,alpha,tau,c,tol,max.iter)
  result
}

imbalance_measure<-function(X,z,o=1,G=NULL,weight=NULL,standardized=T){
  #X: basis function that we want to balance
  #z: treatment combination
  nz=ncol(z)
  n=nrow(X)
  stopifnot(o>0)
  if (is.null(G)){
    G0=expand.grid(rep(list(c(-1,1)), nz))
    G0=G0[,nz:1,drop=FALSE]
    G=c()
    for (oo in 1:o){
      ix=combn(1:nz, oo)
      for (oix in 1:ncol(ix)){
        gixi=apply(G0[,ix[,oix],drop=FALSE],1,prod)
        G=cbind(G,gixi)
      }
    }
  }
  nc=ncol(G)
  
  if (is.data.frame(X)) X=data.matrix(X)
  if(!(is.vector(z)|ncol(z)==1)){
    zz=interaction(z[,1])
    for (k in 2:ncol(z)){
      zz=interaction(zz,z[,k])
    }
  }else{
    zz=interaction(z)
  }
  w=ade4::acm.disjonctif(matrix(zz,ncol=1))
  w=data.matrix(w)
  if (ncol(w)!=nrow(G)) stop("some group is empty")
  
  if (is.null(weight)){
    weight=rep(1,n)
  }
  
  res=c()
  A0=t(t(G)%*%t(w))
  # sdev=numeric(ncol(X))
  for (k in 1:nc){
    AX=A0[,k]*X
    weightk=weight
    weightk[A0[,k]==1]=weightk[A0[,k]==1]/sum(weightk[A0[,k]==1])
    weightk[A0[,k]==-1]=weightk[A0[,k]==-1]/sum(weightk[A0[,k]==-1])
    res=rbind(res,apply(AX,2,function(x) weighted.mean(x,weightk)))
    # if (standardized){
    #   sdev=sqrt((apply(X[which(A0[,k]==1),],2,var)+apply(X[which(A0[,k]==-1),],2,var))/2)
    #   res[k,]=res[k,]/sdev
    # }
  }
  if (standardized){
    sdev=apply(X,2,sd)
    for (k in 1:nc) res[k,]=res[k,]/sdev
  }
  res
}

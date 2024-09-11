factorial_weight_incomplete<-function(X,z,y,o=1,variance=TRUE,model='interaction',lambda0=NULL,alpha=1,tau=0.5,c=0.2,tol=10^(-7),max.iter=10000){
  #X: basis function that we want to balance
  #z: treatment combination: each column is a treatment +1 or -1
  #o: non-negligible order of interaction effects
  
  if (!(model %in%c('additive','interaction'))){
    stop("model has to be either additive or interaction.")
  } 
  
  if (is.vector(y)){
    ny=1
    y=matrix(y,ncol=1)
  }else{
    ny=ncol(y)
  }
  nz=ncol(z)
  n=nrow(X)
  stopifnot(o>0 & o<nz)
  G0=expand.grid(rep(list(c(-1,1)), nz))
  G0=G0[,nz:1,drop=FALSE]
  G=c()
  effect_index=list()
  effect_ix=0
  for (oo in 1:o){
    ix=combn(1:nz, oo)
    for (oix in 1:ncol(ix)){
      effect_ix=effect_ix+1
      effect_index[[effect_ix]]=ix[,oix]
      gixi=apply(G0[,ix[,oix],drop=FALSE],1,prod)
      G=cbind(G,gixi)
    }
  }
  nc=ncol(G)
  Gp=cbind(rep(1,2^nz),G)
  Gm=c()
  for (oo in (o+1):nz){
    ix=combn(1:nz, oo)
    for (oix in 1:ncol(ix)){
      effect_ix=effect_ix+1
      effect_index[[effect_ix]]=ix[,oix]
      gixi=apply(G0[,ix[,oix],drop=FALSE],1,prod)
      Gm=cbind(Gm,gixi)
    }
  }
  
  z.factor=factor(z[,1])
  if (nz>1){
    for (kz in 2:nz){
      z.factor=z.factor:factor(z[,kz])
    }
  }
  tb=table(z.factor)
  observe.ix=as.numeric(which(tb>0))
  unobserve.ix=as.numeric(which(tb==0))
  Gop=Gp[observe.ix,,drop=FALSE]
  Gup=Gp[unobserve.ix,,drop=FALSE]
  Gom=Gm[observe.ix,,drop=FALSE]
  Gum=Gm[unobserve.ix,,drop=FALSE]
  Gnew=t(t(Gop)-t(Gup)%*%(MASS::ginv(t(Gum))%*%t(Gom)))
  
  if (is.data.frame(X)) X=data.matrix(X)
  X=cbind(rep(1,n),X)
  
  #make treatments +1, -1
  for (zi in 1:nz){
    zvaluei=unique(z[,zi])
    stopifnot(length(zvaluei)<=2)
    z[which(z[,zi]==max(zvaluei)),zi]=1
    z[which(z[,zi]==min(zvaluei)),zi]=-1
  }
  
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
  
  Gplus=Gnew
  Gplus[Gplus<0]=0
  Gminus=-Gnew
  Gminus[Gminus<0]=0
  
  zfull=w%*%Gnew
  Aplus=w%*%Gplus
  Aminus=w%*%Gminus
  
  Bm=c()
  bm=c()
  
  if (model=='additive'){
    for (k in 1:nc){
      if (k==1){
        Bm=cbind(Bm,X)
        bm=cbind(bm,2*X)
      }
      Bm=cbind(Bm,Aplus[,k]*X)
      bm=cbind(bm,X)
      for (j in 1:nc){
        CR1=length(effect_index[[k]])==length(effect_index[[j]])
        CR2=(length(effect_index[[k]])!=length(effect_index[[j]])) & 
          (!all(effect_index[[k]])%in%length(effect_index[[j]]))
        if (CR1 | CR2){
          ind_obs=zfull[,j]
          ind_des=Gnew[,j]
          if (k==1){
            Bm=cbind(Bm,ind_obs)
            # bm=cbind(bm,sum(ind_des)/(2^(nz-1)))
            bm=cbind(bm,0)
          }
          Bm=cbind(Bm,Aplus[,k]*ind_obs)
          bm=cbind(bm,sum(Gplus[,k]*ind_des)/(2^(nz-1)))
        }
      }
    }
  }else{
    for (k in 1:nc){
      for (j in 1:nc){
        CR1=length(effect_index[[k]])==length(effect_index[[j]])
        CR2=(length(effect_index[[k]])!=length(effect_index[[j]])) & 
          (!all(effect_index[[k]])%in%length(effect_index[[j]]))
        print(CR1 | CR2)
        if (CR1 | CR2){
          ind_obs=zfull[,j]
          ind_des=Gnew[,j]
          if (k==1){
            Bm=cbind(Bm,ind_obs*X)
            # bm=cbind(bm,sum(ind_des)*X/2^(nz-1))
            bm=cbind(bm,matrix(0,nrow=nrow(X),ncol=ncol(X)))
          }
          Bm=cbind(Bm,(Aplus[,k]*ind_obs*X))
          bm=cbind(bm,sum(Gplus[,k]*ind_des)*X/2^(nz-1))
        }
      }
    }
  }
  qr_res=qr(Bm)
  reduced_Bm=Bm[, qr_res$pivot[1:qr_res$rank]]
  reduced_bm=bm[, qr_res$pivot[1:qr_res$rank]]
  Bm=t(reduced_Bm)
  bm=t(reduced_bm)
  if (is.null(lambda0)) lambda0=-2*MASS::ginv(Bm%*%t(Bm))%*%apply(bm,1,sum)
  lambda_res=gradient_descent(Bm,bm,lambda0,alpha,tau,c,tol,max.iter)
  lambda=lambda_res$lambda
  product=matrix(matrix(lambda,nrow=1)%*%Bm,ncol=1)
  gamma=product
  gamma[which(gamma<0)]=0
  weight=(gamma-product)/2
  if (variance){
    pr=apply(Bm,2,function(v) sum(lambda*v))
    pr_numeric=as.numeric(pr<0)
    psi_gradient=sapply(1:ncol(Bm),function(i) -Bm[,i]*pr[i]*pr_numeric[i]/2-bm[,i])
    Bm_ind=sapply(1:n,function(i) Bm[,i]*pr_numeric[i])
    L_second=MASS::ginv(Bm%*%t(Bm_ind)/n)
  }
  eff=matrix(nrow=nc,ncol=ny)
  # aug_eff=matrix(nrow=nc,ncol=ny)
  var_est=matrix(nrow=nc,ncol=ny)
  for (ei in 1:ny){
    yei=y[,ei]
    for (ec in 1:nc){
      Aec=zfull[,ec]
      yei_weighted=weight*yei*Aec
      mean_weighted=mean(yei_weighted)
      eff[ec,ei]=mean_weighted
      if (variance){
        Bm_ay=apply(Bm_ind,1,function(v) v*Aec*yei)
        L_first=matrix(apply(Bm_ay,2,mean),nrow=1)
        L=matrix(c(L_first%*%L_second,-1),ncol=1)
        eta=rbind(psi_gradient,matrix(yei_weighted-mean_weighted,nrow=1))
        var_est[ec,ei]=mean((t(eta)%*%L)^2)
      }
    }
  }
  list(effects=eff,variance=var_est,weights=weight,G=Gnew,lambda=lambda,gamma=gamma,gr=lambda_res$g,iteration=lambda_res$iter,ee=lambda_res$ee,ss=lambda_res$ss,B=Bm,b=bm)
}

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

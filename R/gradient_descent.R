gradient_descent<-function(B,b,lambda0,alpha=1,tau=0.5,c=0.2,tol=0.0001,max.iter=100){
  
  gradient<-function(B,b,lambda){
    pr=apply(B,2,function(v) sum(lambda*v))
    psi_gradient=sapply(1:ncol(B),function(i) B[,i]*pr[i]*as.numeric(pr[i]<0)/2+b[,i])
    if (nrow(B)==1) mean(psi_gradient)
    else apply(psi_gradient,1,mean)
  }
  
  # backtracking line search
  # Armijo used 1⁄2 for both tau and c in Armijo (1966).
  backtracking <- function(B,b,lambda,g,alpha=1,tau=0.5,c=0.2){
    
    f<-function(B,b,lambda){
      pr=apply(B,2,function(v) sum(lambda*v))
      psi=sapply(1:ncol(B),function(i) (pr[i]^2)*as.numeric(pr[i]<0)/4+sum(b[,i]*lambda))
      mean(psi)
    }
    
    t=c*sum(g^2)
    fv=f(B,b,lambda)
    TF=(fv-f(B,b,lambda-alpha*g))<(alpha*t)
    while(TF){
      alpha=alpha*tau
      TF=(fv-f(B,b,lambda-alpha*g))<(alpha*t)
    }
    alpha
  }
  
  g=gradient(B,b,lambda0)
  error=sqrt(mean(g^2))
  TF=error>=tol
  iter=0
  ss=c()
  ee=error
  while(TF){
    iter=iter+1
    step=backtracking(B,b,lambda0,g,alpha,tau,c)
    ss=c(ss,step)
    lambda=lambda0-step*g
    g=gradient(B,b,lambda)
    error=sqrt(mean(g^2))
    ee=c(ee,error)
    TF=(error>=tol)&(iter<max.iter)
    lambda0=lambda
  }
  list(lambda=lambda0,g=g,error=error,iter=iter,ss=ss,ee=ee)
}


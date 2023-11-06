l.deltar.original<-function(data,KMR.params,family='gaussian'){ #original log likelihood, meaning not including prior log likelihood
  y=data$y
  X=data$X
  beta=KMR.params$beta
  tau=KMR.params$tau
  K=KMR.params$K
  phi=KMR.params$phi
  if (family=='gaussian') return (dmvnorm(as.vector(y-X%*%beta),as.vector(rep(0,length(y))),diag(phi)+tau*K,log=TRUE))
} 

l.deltar.original.diff<-function(data,KMR.params,KMR.params.star,family='gaussian',PC.num=10,int.num=100000,int.h.stat=FALSE){ 
  #KMR.params.star-KMR.params
  y=data$y
  X=data$X
  beta=KMR.params$beta
  tau=KMR.params$tau
  phi=KMR.params$phi
  K=KMR.params$K
  K.star=KMR.params.star$K
  
  if (int.h.stat==TRUE){
    # if (family=='gaussian'){
    #   V.temp=diag(length(y))+tau/phi[1]*K
    #   V.star.temp=diag(length(y))+tau/phi[1]*K.star
    #   return(
    #     1/2*log(det(V.temp))-1/2*log(det(V.star.temp))+(
    #       -phi[1]^(-1)*1/2*t(y-X%*%beta)%*%(solve(V.star.temp)-solve(V.temp))%*%(y-X%*%beta)
    #     )
    #   )
    # }
    
    if (family=='gaussian'){
      b<-function(theta) theta^2/2
      b.p<-function(theta){
        return (theta)
      }
      b.pp<-function(theta) return (matrix(rep(1,length(theta)),ncol=1))
      g<-function(mu) mu
      g.p<-function(mu) return (matrix(rep(1,length(theta)),ncol=1))
    }
    
    if (family=='poisson'){
      b<-function(theta) exp(theta)
      b.p<-function(theta){
        return (exp(theta))
      } 
      b.pp<-function(theta) return (exp(theta))
      g<-function(mu) log(mu)
      g.p<-function(mu) return (1/mu)
    }
    
    
    if (family=='binomial-probit'){
      b<-function(theta) log(1+exp(theta))
      b.p<-function(theta){
        exp.theta=exp(theta)
        return (exp.theta/(1+exp.theta))
      }
      b.pp<-function(theta) exp(theta)/(1+exp(theta))^2
      g<-function(mu) qnorm(mu)
      g.p<-function(mu) 1/dnorm(qnorm(mu))
    }
    
    # #Option 1: this is using Iterative reweighted least sequare form - 
    # #assuming g(Y) \approx g(\mu)+g'(\mu)(Y-\mu) ~ normal(X\beta+h,\phi b''(\theta)g'(\mu)^2)
    # # then integrate h out - when family == gaussian, this method is equivalent to closed form gaussian
    eta=X%*%KMR.params$beta+KMR.params$A%*%KMR.params$h
    theta=eta
    mu=b.p(theta)
    y.tilde=eta+(y-mu)*g.p(mu)
    V.temp=diag(as.vector(b.pp(theta)*g.p(mu)^2))+KMR.params$tau/KMR.params$phi[1]*KMR.params$K
    
    eta.star=X%*%KMR.params.star$beta+KMR.params.star$A%*%KMR.params.star$h
    theta.star=eta.star
    mu.star=b.p(theta.star)
    y.tilde.star=eta.star+(y-mu.star)*g.p(mu.star)
    V.temp.star=diag(as.vector(b.pp(theta.star)*g.p(mu.star)^2))+KMR.params.star$tau/KMR.params.star$phi[1]*KMR.params.star$K
    return(
      (-1/2*log(det(V.temp.star))-1/2*1/KMR.params.star$phi[1]*t(y.tilde.star-X%*%KMR.params.star$beta)%*%(solve(V.temp.star))%*%(y.tilde.star-X%*%KMR.params.star$beta))
      -(-1/2*log(det(V.temp))-1/2*1/KMR.params$phi[1]*t(y.tilde-X%*%KMR.params$beta)%*%(solve(V.temp))%*%(y.tilde-X%*%KMR.params$beta))
    )
    
    # #Option 2: this is numerical integration approximation, seems do not work well
    # eig=KMR.params$eig
    # eig.star=KMR.params.star$eig
    # tau.star=KMR.params.star$tau
    # beta.star=KMR.params.star$beta
    # L=max(sum(eig[1:PC.num]>1),sum(eig.star[1:PC.num]>1) )
    # A=KMR.params$A[,1:L]
    # A.star=KMR.params.star$A[,1:L]
    # h.mat=matrix(rnorm(L*int.num,0,sd=1),ncol=int.num,nrow=L)
    # A.h=A%*%h.mat*sqrt(tau)
    # A.h.star=A.star%*%h.mat*sqrt(tau.star)
    # rm(h.mat)
    # X.beta=X%*%beta
    # X.beta.star=X%*%beta.star
    # if (family != 'binomial'){
    #   phi=KMR.params$phi
    #   lk=exp(colSums((A.h*y-b(matrix(X.beta,nrow=nrow(A.h),ncol=ncol(A.h))+A.h))/as.vector(phi))) #need notice!!20221226
    #   lk.star=exp(colSums((A.h.star*y-b(matrix(X.beta.star,nrow=nrow(A.h.star),ncol=ncol(A.h.star))+A.h.star))/as.vector(phi)))
    #   return (log(mean(lk.star))-log(mean(lk))+sum((X%*%(beta.star-beta))*(y/phi))) #here we assume phi=phi.star=rep(1,length(y))
    # }
    # if (family=='binomial'){
    #   exp.theta.ratio.star=1/(1+exp(matrix(X.beta.star,nrow=nrow(A.h.star),ncol=ncol(A.h.star))+A.h.star))
    #   exp.theta.ratio=1/(1+exp(matrix(X.beta,nrow=nrow(A.h),ncol=ncol(A.h))+A.h))
    #   exp.theta.ratio.star[as.logical(y),]=1-exp.theta.ratio.star[as.logical(y),]
    #   exp.theta.ratio[as.logical(y),]=1-exp.theta.ratio[as.logical(y),]
    #   return (log(mean(colProds(exp.theta.ratio.star)))-log(mean(colProds(exp.theta.ratio))))
    # }
  }
  if (int.h.stat==FALSE){
    h=KMR.params$h
    A=KMR.params$A
    h.star=KMR.params.star$h
    beta.star=KMR.params.star$beta
    A.star=KMR.params.star$A
    if (family == 'gaussian'){
      phi=KMR.params.star$phi
      return (-sum((y-X%*%beta.star-A.star%*%h.star)^2/phi[1])/2+sum((y-X%*%beta-A%*%h)^2/phi[1])/2)
    }
    if (family=='binomial'){
      exp.theta.ratio=1/(1+exp(X%*%beta+A%*%h))
      exp.theta.star.ratio=1/(1+exp(X%*%beta.star+A.star%*%h.star))
      exp.theta.ratio[as.logical(y)]=1-exp.theta.ratio[as.logical(y)]
      exp.theta.star.ratio[as.logical(y)]=1-exp.theta.star.ratio[as.logical(y)]
      return (log(prod(exp.theta.star.ratio))-log(prod(exp.theta.ratio)))
    }
    
    if (family=='binomial-probit'){
      exp.theta.ratio=pnorm(X%*%beta+A%*%h)
      exp.theta.star.ratio=pnorm(X%*%beta.star+A.star%*%h.star)
      exp.theta.ratio[as.logical(1-y)]=1-exp.theta.ratio[as.logical(1-y)]
      exp.theta.star.ratio[as.logical(1-y)]=1-exp.theta.star.ratio[as.logical(1-y)]
      return (log(prod(exp.theta.star.ratio))-log(prod(exp.theta.ratio)))
    }
    
    if (family=='poisson'){
      log.lambda=(X%*%beta+A%*%h)
      log.lambda.star=(X%*%beta.star+A.star%*%h.star)
      return (sum((log.lambda.star-log.lambda)*y+(exp(log.lambda)-exp(log.lambda.star))))
    }
    
    if (!(family %in% c('binomial','poisson','binomial-probit'))) stop('Need to specify the likelihood function')
  }
}

l.r.prior<-function(m.choose,KMR.params,prior.params,r.prior.dist='invuniform'){
  #f(r[m]|delta[m]=1) used in move2
  r=KMR.params$r
  if (r.prior.dist=='invuniform'){
    a.r.prior = prior.params$a.r.if_invunif.prior
    b.r.prior = prior.params$b.r.if_invunif.prior
    return (ifelse(1/b.r.prior <= r[m.choose] & r[m.choose] <= 1/a.r.prior, -2*log(r[m.choose]) - log(b.r.prior - a.r.prior), log(0)))
  }
  if (r.prior.dist=='uniform'){
    a.r.prior = prior.params$a.r.if_unif.prior
    b.r.prior = prior.params$b.r.if_unif.prior
    return (ifelse(a.r.prior <= r[m.choose] & r[m.choose] <= b.r.prior,  - log(b.r.prior - a.r.prior), log(0)))
  }
}

l.r.in_deltar.prior<-function(m.choose,KMR.params,prior.params,r.prior.dist='invuniform'){
  #f(r[m]|delta[m])
  r=KMR.params$r
  delta=KMR.params$delta
  if (r.prior.dist=='invuniform'){
    a.r.prior = prior.params$a.r.if_invunif.prior
    b.r.prior = prior.params$b.r.if_invunif.prior
    return (ifelse (delta[m.choose]==TRUE,
                    ifelse(1/b.r.prior <= r[m.choose] & r[m.choose] <= 1/a.r.prior, -2*log(r[m.choose]) - log(b.r.prior - a.r.prior), log(0)),
                    log(1))
    )
  }
  if (r.prior.dist=='uniform'){
    a.r.prior = prior.params$a.r.if_unif.prior
    b.r.prior = prior.params$b.r.if_unif.prior
    return (ifelse (delta[m.choose]==TRUE,
                    ifelse(a.r.prior <= r[m.choose] & r[m.choose] <= b.r.prior, - log(b.r.prior - a.r.prior), log(0)),
                    log(1))
    )
  }
}


l.delta.prior<-function(KMR.params,prior.params){
  #f(delta); after integrate p out, p=prob(delta[i]=1)
  delta=KMR.params$delta
  a.delta.prior=prior.params$a.delta.prior
  b.delta.prior=prior.params$b.delta.prior
  return (-lgamma(a.delta.prior+sum(delta))-lgamma(length(delta)-sum(delta)+b.delta.prior)) #?? 20230804
}

# # updating r solely; currently no use
# 
# generate.r.in_r.prop <- function(m.choose,KMR.params,prop.params,r.prior.dist='invuniform'){
#   r=KMR.params$r
#   r.sd.prop           = prop.params$r.in_r.sd.prop
#   if (r.prior.dist=='invuniform'){
#     a.r.prop = prop.params$a.r.in_r.if_invunif.prop
#     b.r.prop = prop.params$b.r.in_r.if_invunif.prop
#     return(rtruncnorm(1, a = 1/b.r.prop, b = 1/a.r.prop, mean = r[m.choose], sd = r.sd.prop))
#   }
# }
# 
# l.r.in_r.prop<-function(m.choose,current.KMR.params,after.KMR.params,prop.params,r.prior.dist='invuniform'){
#   r.current=current.KMR.params$r
#   r.after=after.KMR.params$r
#   if (r.prior.dist=='invuniform'){
#     r.sd.prop           = prop.params$r.in_r.sd.prop
#     a.r.prop = prop.params$a.r.in_r.if_invunif.prop
#     b.r.prop = prop.params$b.r.in_r.if_invunif.prop
#     return(log(dtruncnorm(r.after[m.choose], a = 1/b.r.prop, b = 1/a.r.prop, mean = r.current[m.choose], sd = r.sd.prop)))
#   }
# }

#move 1 in updating delta,r 
generate.r.in_deltar.move1.prop <- function(m.choose,prop.params,r.prior.dist='invuniform'){
  if (r.prior.dist=='invuniform'){
    a.r.prop = prop.params$a.r.in_deltar.if_invunif.move1.prop
    b.r.prop = prop.params$b.r.in_deltar.if_invunif.move1.prop
    return(1/runif(1,a.r.prop,b.r.prop))
  }
  if (r.prior.dist=='uniform'){
    a.r.prop = prop.params$a.r.in_deltar.if_unif.move1.prop
    b.r.prop = prop.params$b.r.in_deltar.if_unif.move1.prop
    return(runif(1,a.r.prop,b.r.prop))
  }
}


l.deltar.move1.prop<-function(m.choose,current.KMR.params,after.KMR.params,prop.params,r.prior.dist='invuniform'){ #log prob of proposal function  delta=KMR.params$delta
  r.current=current.KMR.params$r
  r.after=after.KMR.params$r
  delta.current=current.KMR.params$delta
  delta.after=after.KMR.params$delta

  if (r.prior.dist=='invuniform'){ 
    a.r.prop=prop.params$a.r.in_deltar.if_invunif.move1.prop
    b.r.prop=prop.params$b.r.in_deltar.if_invunif.move1.prop
    if (sum(delta.current)==0 & delta.current[m.choose]==0) res=ifelse(1/b.r.prop <= r.after[m.choose] & r.after[m.choose] <= 1/a.r.prop, -2*log(r.after[m.choose]) - log(b.r.prop - a.r.prop), log(0)) -log(length(delta.current))
    if (sum(delta.current)!=0 & delta.current[m.choose]==0) res=ifelse(1/b.r.prop <= r.after[m.choose] & r.after[m.choose] <= 1/a.r.prop, -2*log(r.after[m.choose]) - log(b.r.prop - a.r.prop), log(0)) -log(2)-log(length(delta.current))
  }
  if (r.prior.dist=='uniform'){ 
    a.r.prop=prop.params$a.r.in_deltar.if_unif.move1.prop
    b.r.prop=prop.params$b.r.in_deltar.if_unif.move1.prop
    if (sum(delta.current)==0 & delta.current[m.choose]==0) res=ifelse(a.r.prop <= r.after[m.choose] & r.after[m.choose] <= b.r.prop, - log(b.r.prop - a.r.prop), log(0)) -log(length(delta.current))
    if (sum(delta.current)!=0 & delta.current[m.choose]==0) res=ifelse(a.r.prop <= r.after[m.choose] & r.after[m.choose] <= b.r.prop, - log(b.r.prop - a.r.prop), log(0)) -log(2)-log(length(delta.current))
  }
  if (delta.current[m.choose]==1) res=-log(2)-log(length(delta.current))
  return (res)
}

#move 2 in updating delta,r 
generate.r.in_deltar.move2.prop <- function(m.choose,KMR.params,prop.params,r.prior.dist='invuniform'){
  r=KMR.params$r
  if (r.prior.dist=='invuniform'){
    r.sd.prop           = prop.params$r.in_deltar.sd.move2.prop
    a.r.prop = prop.params$a.r.in_deltar.if_invunif.move2.prop
    b.r.prop = prop.params$b.r.in_deltar.if_invunif.move2.prop
    return(rtruncnorm(1, a = 1/b.r.prop, b = 1/a.r.prop, mean = r[m.choose], sd = r.sd.prop))
  }
  if (r.prior.dist=='uniform'){
    r.sd.prop           = prop.params$r.in_deltar.sd.move2.prop
    a.r.prop = prop.params$a.r.in_deltar.if_unif.move2.prop
    b.r.prop = prop.params$b.r.in_deltar.if_unif.move2.prop
    return(rtruncnorm(1, a = a.r.prop, b = b.r.prop, mean = r[m.choose], sd = r.sd.prop))
  }
}


l.r.in_deltar.move2.prop<-function(m.choose,current.KMR.params,after.KMR.params,prop.params,r.prior.dist='invuniform'){
  r.current=current.KMR.params$r
  r.after=after.KMR.params$r
  r.sd.prop           = prop.params$r.in_deltar.sd.move2.prop
  if (r.prior.dist=='invuniform'){
    a.r.prop = prop.params$a.r.in_deltar.if_invunif.move2.prop
    b.r.prop = prop.params$b.r.in_deltar.if_invunif.move2.prop
    return(log(dtruncnorm(r.after[m.choose], a = 1/b.r.prop, b = 1/a.r.prop, mean = r.current[m.choose], sd = r.sd.prop)))
  }
  if (r.prior.dist=='uniform'){
    a.r.prop = prop.params$a.r.in_deltar.if_unif.move2.prop
    b.r.prop = prop.params$b.r.in_deltar.if_unif.move2.prop
    return(log(dtruncnorm(r.after[m.choose], a = a.r.prop, b = b.r.prop, mean = r.current[m.choose], sd = r.sd.prop)))
  }
}

#When r and delta are changed, K and A must also be updated
update.KMR.params<-function(KMR.params,data){
  res=KMR.params
  delta=KMR.params$delta
  r=KMR.params$r
  Z.diff.precal=data$Z.diff.precal
  if (!all.equal(r==0,delta==FALSE)) stop('when delta[i]=0, r[i] must equals 0, vice versa.')
  if (sum(delta)==0){
    K=matrix(1,nrow=nrow(KMR.params$K),ncol=ncol(KMR.params$K))
  }else if (sum(delta)==1){
    K=exp(-r[delta]*Z.diff.precal[,,delta])
  }else{
    r.Z.diff.precal.temp=array(NA,c(nrow(KMR.params$K),ncol(KMR.params$K),sum(delta)))
    count=0
    for (i in which(delta==TRUE)){
      count=count+1
      r.Z.diff.precal.temp[,,count]=Z.diff.precal[,,i]*r[i]
    }
    K=exp(-rowSums(r.Z.diff.precal.temp,dims=2))
  }
  temp=eigen(K)
  temp$values[(temp$values)<0]=0
  eig=temp$values
  A=temp$vectors%*%diag(sqrt(temp$values))
  res$K=K
  res$A=A
  res$eig=eig
  return (res)
}

#this is to compare whether the approx and true llk diff are close
l.deltar.original.diff.closedform<-function(data,KMR.params,KMR.params.star,family='gaussian'){ 
  #KMR.params.star-KMR.params
  y=data$y
  X=data$X
  beta=KMR.params$beta
  tau=KMR.params$tau
  K=KMR.params$K
  K.star=KMR.params.star$K
  if (family=='gaussian'){
    V.temp=diag(length(y))+tau/phi[1]*K
    V.star.temp=diag(length(y))+tau/phi[1]*K.star
    return(
      1/2*log(det(V.temp))-1/2*log(det(V.star.temp))+(
        -phi[1]^(-1)*1/2*t(y-X%*%beta)%*%(solve(V.star.temp)-solve(V.temp))%*%(y-X%*%beta)
      )
    )
  }
}


update.expfam.params<-function(data,KMR.params,family){
  if (family=='binomial'){
    b<-function(theta) log(1+exp(theta))
    b.p<-function(theta){
      exp.theta=exp(theta)
      return (exp.theta/(1+exp.theta))
    }
    b.pp<-function(theta) exp(theta)/(1+exp(theta))^2
    g<-function(mu) log(mu/(1-mu))
    g.p<-function(mu) 1/mu/(1-mu)
  }
  
  if (family=='gaussian'){
    b<-function(theta) theta^2/2
    b.p<-function(theta){
      return (theta)
    } 
    b.pp<-function(theta) return (matrix(rep(1,length(theta)),ncol=1))
    g<-function(mu) mu
    g.p<-function(mu) return (matrix(rep(1,length(theta)),ncol=1))
  }
  
  if (family=='poisson'){
    b<-function(theta) exp(theta)
    b.p<-function(theta){
      return (exp(theta))
    } 
    b.pp<-function(theta) return (exp(theta))
    g<-function(mu) log(mu)
    g.p<-function(mu) return (1/mu)
  }
  
  
  if (family=='binomial-probit'){
    b<-function(theta) log(1+exp(theta))
    b.p<-function(theta){
      exp.theta=exp(theta)
      return (exp.theta/(1+exp.theta))
    }
    b.pp<-function(theta) exp(theta)/(1+exp(theta))^2
    g<-function(mu) qnorm(mu)
    g.p<-function(mu) 1/dnorm(qnorm(mu))
  }
  
  #input parameters from data and KMR.params
  y=data$y
  X=data$X
  Z=data$Z
  Z.diff.precal=data$Z.diff.precal
  beta=KMR.params$beta
  K=KMR.params$K
  A=KMR.params$A
  eig=KMR.params$eig
  tau=KMR.params$tau
  r=KMR.params$r
  delta=KMR.params$delta
  h=KMR.params$h
  phi=KMR.params$phi
  
  eta=X%*%beta+A%*%h
  theta=eta
  mu=b.p(theta)
  y.tilde=eta+(y-mu)*g.p(mu)
  W.diag=as.vector(1/b.pp(theta)/phi/g.p(mu)^2)
  
  return (list(eta=eta,theta=theta,mu=mu,y.tilde=y.tilde,W.diag=W.diag))
}

update.beta<-function(data,KMR.params,prior.params,family){
  
  if (family=='binomial'){
    b<-function(theta) log(1+exp(theta))
    b.p<-function(theta){
      exp.theta=exp(theta)
      return (exp.theta/(1+exp.theta))
    }
    b.pp<-function(theta) exp(theta)/(1+exp(theta))^2
    g<-function(mu) log(mu/(1-mu))
    g.p<-function(mu) 1/mu/(1-mu)
  }
  
  if (family=='gaussian'){
    b<-function(theta) theta^2/2
    b.p<-function(theta){
      return (theta)
    } 
    b.pp<-function(theta) return (matrix(rep(1,length(theta)),ncol=1))
    g<-function(mu) mu
    g.p<-function(mu) return (matrix(rep(1,length(theta)),ncol=1))
  }
  
  if (family=='poisson'){
    b<-function(theta) exp(theta)
    b.p<-function(theta){
      return (exp(theta))
    } 
    b.pp<-function(theta) return (exp(theta))
    g<-function(mu) log(mu)
    g.p<-function(mu) return (1/mu)
  }
  
  
  if (family=='binomial-probit'){
    b<-function(theta) log(1+exp(theta))
    b.p<-function(theta){
      exp.theta=exp(theta)
      return (exp.theta/(1+exp.theta))
    }
    b.pp<-function(theta) exp(theta)/(1+exp(theta))^2
    g<-function(mu) qnorm(mu)
    g.p<-function(mu) 1/dnorm(qnorm(mu))
  }
  
  y=data$y
  X=data$X
  Z=data$Z
  Z.diff.precal=data$Z.diff.precal
  beta=KMR.params$beta
  K=KMR.params$K
  A=KMR.params$A
  eig=KMR.params$eig
  tau=KMR.params$tau
  r=KMR.params$r
  delta=KMR.params$delta
  h=KMR.params$h
  phi=KMR.params$phi
  lambda=KMR.params$lambda
  expfam.params=update.expfam.params(data=data,KMR.params=KMR.params,family=family)
  eta=expfam.params$eta
  theta=expfam.params$theta
  mu=expfam.params$mu
  y.tilde=expfam.params$y.tilde
  W.diag=expfam.params$W.diag
  mu.beta.prior=prior.params$mu.beta.prior
  V.beta.prior=prior.params$V.beta.prior
  
  
  
  C.beta=solve(solve(V.beta.prior)+t(X)%*%(W.diag*X))
  m.beta=C.beta%*%(solve(V.beta.prior)%*%mu.beta.prior+t(X)%*%(W.diag*(y.tilde-A%*%h)))
  log.pi.beta=sum((y*theta-b(theta))/phi)-t(beta-mu.beta.prior)%*%solve(V.beta.prior)%*%(beta-mu.beta.prior)/2
  
  beta.star=matrix(rmvnorm(1,m.beta,C.beta),ncol=1)
  eta.star=X%*%beta.star+A%*%h
  theta.star=eta.star
  mu.star=b.p(theta.star)
  y.tilde.star=eta.star+(y-mu.star)*g.p(mu.star)
  W.diag.star=(as.vector(1/b.pp(theta.star)/phi/g.p(mu.star)^2))
  C.beta.star=solve(solve(V.beta.prior)+t(X)%*%(W.diag.star*X))
  m.beta.star=C.beta.star%*%(solve(V.beta.prior)%*%mu.beta.prior+t(X)%*%(W.diag.star*(y.tilde.star-A%*%h)))
  log.pi.beta.star=sum((y*theta.star-b(theta.star))/phi)-t(beta-mu.beta.prior)%*%solve(V.beta.prior)%*%(beta-mu.beta.prior)/2
  
  log.q.tminus1.star.beta=dmvnorm(as.vector(beta),as.vector(m.beta.star),C.beta.star,log=TRUE)
  log.q.star.tminus1.beta=dmvnorm(as.vector(beta.star),as.vector(m.beta),C.beta,log=TRUE)
  
  log.gamma.prob.beta=log.q.tminus1.star.beta-log.q.star.tminus1.beta+log.pi.beta.star-log.pi.beta
  gamma.prob.beta=exp(log.gamma.prob.beta)
  
  beta.na.flag=FALSE
  beta.acpt.flag=FALSE
  if (is.nan(gamma.prob.beta)){
    beta=beta
    beta.na.flag=TRUE
  }else{
    if (runif(1)<=gamma.prob.beta){
      beta=beta.star
      beta.acpt.flag=TRUE
    } else {
      beta=beta
    }
  }
  KMR.params$beta=beta
  # print(c('beta: ',beta))
  return (list(KMR.params=KMR.params,beta.na.flag=beta.na.flag,beta.acpt.flag=beta.acpt.flag))
  
  
  # if (family=='gaussian'){
  #   V.temp=diag(length(y))+lambda*K
  #   beta.gaussian.covariance.temp=phi[1]*solve(t(X)%*%solve(V.temp)%*%X+solve(V.beta.prior))
  #   beta.gaussian.mean.temp=beta.gaussian.covariance.temp%*%t(t(y)%*%solve(phi[1]*V.temp)%*%X+t(mu.beta.prior)%*%V.beta.prior)
  #   beta=matrix(rmvnorm(1,mean=beta.gaussian.mean.temp,sigma=beta.gaussian.covariance.temp),nrow=ncol(X),ncol=1)
  #   KMR.params$beta=beta
  #   return (list(KMR.params=KMR.params,beta.na.flag=FALSE,beta.acpt.flag=TRUE))
  # }
  
}

update.h<-function(data,KMR.params,family){
  
  if (family=='binomial'){
    b<-function(theta) log(1+exp(theta))
    b.p<-function(theta){
      exp.theta=exp(theta)
      return (exp.theta/(1+exp.theta))
    }
    b.pp<-function(theta) exp(theta)/(1+exp(theta))^2
    g<-function(mu) log(mu/(1-mu))
    g.p<-function(mu) 1/mu/(1-mu)
  }
  
  if (family=='gaussian'){
    b<-function(theta) theta^2/2
    b.p<-function(theta){
      return (theta)
    } 
    b.pp<-function(theta) return (matrix(rep(1,length(theta)),ncol=1))
    g<-function(mu) mu
    g.p<-function(mu) return (matrix(rep(1,length(theta)),ncol=1))
  }
  
  if (family=='poisson'){
    b<-function(theta) exp(theta)
    b.p<-function(theta){
      return (exp(theta))
    } 
    b.pp<-function(theta) return (exp(theta))
    g<-function(mu) log(mu)
    g.p<-function(mu) return (1/mu)
  }
  
  
  if (family=='binomial-probit'){
    b<-function(theta) log(1+exp(theta))
    b.p<-function(theta){
      exp.theta=exp(theta)
      return (exp.theta/(1+exp.theta))
    }
    b.pp<-function(theta) exp(theta)/(1+exp(theta))^2
    g<-function(mu) qnorm(mu)
    g.p<-function(mu) 1/dnorm(qnorm(mu))
  }
  
  y=data$y
  X=data$X
  Z=data$Z
  Z.diff.precal=data$Z.diff.precal
  beta=KMR.params$beta
  K=KMR.params$K
  A=KMR.params$A
  eig=KMR.params$eig
  tau=KMR.params$tau
  r=KMR.params$r
  delta=KMR.params$delta
  h=KMR.params$h
  phi=KMR.params$phi
  expfam.params=update.expfam.params(data=data,KMR.params=KMR.params,family=family)
  eta=expfam.params$eta
  theta=expfam.params$theta
  mu=expfam.params$mu
  y.tilde=expfam.params$y.tilde
  W.diag=expfam.params$W.diag
  
  
  
  h.acpt.flag.vec=rep(FALSE,length(y))
  h.na.flag.vec=rep(FALSE,length(y))
  L=length(y) #which(eig==0)[1]-1 #
  for (i in 1:L){  #update h_i  
    eta=X%*%beta+A%*%h
    theta=eta
    mu=b.p(theta)
    y.tilde=eta+(y-mu)*g.p(mu)
    W.diag=as.vector(1/b.pp(theta)/phi/g.p(mu)^2)
    C.h=(tau^(-1)+sum(A[,i]^2*W.diag))^(-1)
    m.h=C.h*sum(A[,i]*(y.tilde-X%*%beta-A[,-i]%*%h[-i])*W.diag)
    log.pi.h=(sum((y*theta-b(theta))/phi)-tau^(-1)*h[i]^2/2)
    
    h.star=h
    
    h.star[i]=rnorm(1,m.h,sqrt(C.h))
    eta.star=X%*%beta+A%*%h.star
    theta.star=eta.star
    mu.star=b.p(theta.star)
    y.tilde.star=eta.star+(y-mu.star)*g.p(mu.star)
    W.diag.star=as.vector(1/b.pp(theta.star)/phi/g.p(mu.star)^2)
    C.h.star=(tau^(-1)+sum(A[,i]^2*W.diag.star))^(-1)
    m.h.star=C.h.star*sum(A[,i]*(y.tilde.star-X%*%beta-A[,-i]%*%h.star[-i])*W.diag.star)
    log.pi.h.star=(sum((y*theta.star-b(theta.star))/phi)-tau^(-1)*h.star[i]^2/2)
    
    log.q.tminus1.star.h=dnorm(h[i],m.h.star,sqrt(C.h.star),log=TRUE)
    log.q.star.tminus1.h=dnorm(h.star[i],m.h,sqrt(C.h),log=TRUE)
    
    log.gamma.prob.h=log.q.tminus1.star.h-log.q.star.tminus1.h+(log.pi.h.star-log.pi.h)
    gamma.prob.h=exp(log.gamma.prob.h) 
    
    if (is.nan(gamma.prob.h)){
      h[i]=h[i]
      h.na.flag.vec[i]=TRUE
    }else{
      if (runif(1)<=gamma.prob.h){
        h[i]=h.star[i]
        h.acpt.flag.vec[i]=TRUE
      } else {
        h[i]=h[i]
      }
    }
  }
  KMR.params$h=h
  return (list(KMR.params=KMR.params,h.na.flag.vec=h.na.flag.vec,h.acpt.flag.vec=h.acpt.flag.vec))
  
  
}

update.tau<-function(data,KMR.params,prop.params,prior.params,family){
  y=data$y
  X=data$X
  Z=data$Z
  Z.diff.precal=data$Z.diff.precal
  beta=KMR.params$beta
  K=KMR.params$K
  A=KMR.params$A
  eig=KMR.params$eig
  tau=KMR.params$tau
  r=KMR.params$r
  delta=KMR.params$delta
  h=KMR.params$h
  phi=KMR.params$phi
  lambda.sd.prop=prop.params$lambda.sd.prop
  a.lambda.prior=prior.params$a.lambda.prior
  b.lambda.prior=prior.params$b.lambda.prior
  
  a.tau.prior=prior.params$a.tau.prior
  b.tau.prior=prior.params$b.tau.prior
  
  tau.na.flag=FALSE
  tau.acpt.flag=FALSE # only valid in gaussian
  
  tau=rinvgamma(1,length(y)/2+a.tau.prior,(t(h)%*%h)/2+b.tau.prior)
  KMR.params$tau=tau
  KMR.params$lambda=tau/phi[1] #no use for family != gaussian, just keep lambda updated
  # print(c('tau: ',tau))
  return (list(KMR.params=KMR.params))
  
}

update.deltar<-function(data,KMR.params,prop.params,prior.params,r.prior.dist,family,PC.num,int.num,int.h.stat){
  y=data$y
  X=data$X
  Z=data$Z
  Z.diff.precal=data$Z.diff.precal
  beta=KMR.params$beta
  K=KMR.params$K
  A=KMR.params$A
  eig=KMR.params$eig
  tau=KMR.params$tau
  r=KMR.params$r
  delta=KMR.params$delta
  h=KMR.params$h
  phi=KMR.params$phi
  
  KMR.params.star=KMR.params
  delta.star=delta
  r.star=r
  move.type=ifelse(sum(delta)==0,1,sample(c(1,2),1))
  if (move.type==1){
    m.choose=ifelse(length(delta)>1,sample(1:length(delta),1),1)
    if (delta[m.choose]==TRUE){
      delta.star[m.choose]=FALSE
      r.star[m.choose]=0
      KMR.params.star$delta=delta.star
      KMR.params.star$r=r.star
      KMR.params.star=update.KMR.params(KMR.params=KMR.params.star,data=data)
    }
    if (delta[m.choose]==FALSE){
      r.star[m.choose]=generate.r.in_deltar.move1.prop(m.choose=m.choose,prop.params=prop.params,r.prior.dist=r.prior.dist)
      delta.star[m.choose]=TRUE
      KMR.params.star$delta=delta.star
      KMR.params.star$r=r.star
      KMR.params.star=update.KMR.params(KMR.params=KMR.params.star,data=data)
    }
    acpt.rate.temp=exp(
      l.deltar.original.diff(data=data,KMR.params=KMR.params,KMR.params.star=KMR.params.star,family=family,PC.num=PC.num,int.num=int.num,int.h.stat=int.h.stat)
      +l.r.in_deltar.prior(m.choose=m.choose,KMR.params=KMR.params.star,prior.params=prior.params,r.prior.dist=r.prior.dist)
      -l.r.in_deltar.prior(m.choose=m.choose,KMR.params=KMR.params,prior.params=prior.params,r.prior.dist=r.prior.dist)
      +l.delta.prior(KMR.params.star,prior.params)-l.delta.prior(KMR.params,prior.params)
      +l.deltar.move1.prop(m.choose=m.choose,current.KMR.params=KMR.params.star,after.KMR.params=KMR.params,prop.params=prop.params,r.prior.dist=r.prior.dist)
      -l.deltar.move1.prop(m.choose=m.choose,current.KMR.params=KMR.params,after.KMR.params=KMR.params.star,prop.params=prop.params,r.prior.dist=r.prior.dist)
    )
    # test.record[t]=l.deltar.original.diff(data=data,KMR.params=KMR.params,KMR.params.star=KMR.params.star,family=family)-
    # l.deltar.original.diff.closedform(data=data,KMR.params=KMR.params,KMR.params.star=KMR.params.star,family=family)
    deltar.acpt.flag=FALSE
    deltar.na.flag=FALSE
    if (!is.nan(acpt.rate.temp)){
      if(runif(1)<acpt.rate.temp){
        delta[m.choose]=delta.star[m.choose]
        r[m.choose]=r.star[m.choose]
        deltar.acpt.flag=TRUE
      }
    }else{
      deltar.na.flag=TRUE
    }
  }
  if (move.type==2){
    m.choose=ifelse(sum(delta)>1,sample(which(delta==TRUE),1),which(delta==TRUE))
    r.star[m.choose]=generate.r.in_deltar.move2.prop(m.choose=m.choose,KMR.params=KMR.params,prop.params=prop.params,r.prior.dist=r.prior.dist)
    KMR.params.star$r=r.star
    KMR.params.star=update.KMR.params(KMR.params=KMR.params.star,data=data)
    acpt.rate.temp=exp(
      l.deltar.original.diff(data=data,KMR.params=KMR.params,KMR.params.star=KMR.params.star,family=family,PC.num=PC.num,int.num=int.num,int.h.stat=int.h.stat)
      +l.r.prior(m.choose=m.choose,KMR.params=KMR.params.star,prior.params=prior.params,r.prior.dist=r.prior.dist)
      -l.r.prior(m.choose=m.choose,KMR.params=KMR.params,prior.params=prior.params,r.prior.dist=r.prior.dist)
      +l.r.in_deltar.move2.prop(m.choose=m.choose,current.KMR.params=KMR.params.star,after.KMR.params=KMR.params,prop.params=prop.params,r.prior.dist=r.prior.dist)
      -l.r.in_deltar.move2.prop(m.choose=m.choose,current.KMR.params=KMR.params,after.KMR.params=KMR.params.star,prop.params=prop.params,r.prior.dist=r.prior.dist)
    )
    
    # test.record[t]=l.deltar.original.diff(data=data,KMR.params=KMR.params,KMR.params.star=KMR.params.star,family=family)-
    #   l.deltar.original.diff.closedform(data=data,KMR.params=KMR.params,KMR.params.star=KMR.params.star,family=family)
    deltar.acpt.flag=FALSE
    deltar.na.flag=FALSE
    if (!is.nan(acpt.rate.temp)){
      if(runif(1)<acpt.rate.temp){
        r[m.choose]=r.star[m.choose]
        deltar.acpt.flag=TRUE
      }
    }else{
      deltar.na.flag=TRUE  
    }
  }
  # print(c('r: ',r))
  KMR.params$r=r
  KMR.params$delta=delta
  KMR.params=update.KMR.params(KMR.params=KMR.params,data=data)
  return (list(KMR.params=KMR.params,deltar.acpt.flag=deltar.acpt.flag,deltar.na.flag=deltar.na.flag,move.type=move.type))
  
}


update.deltar.fixdelta<-function(data,KMR.params,prop.params,prior.params,r.prior.dist,family,PC.num,int.num,int.h.stat){
  y=data$y
  X=data$X
  Z=data$Z
  Z.diff.precal=data$Z.diff.precal
  beta=KMR.params$beta
  K=KMR.params$K
  A=KMR.params$A
  eig=KMR.params$eig
  tau=KMR.params$tau
  r=KMR.params$r
  delta=KMR.params$delta
  h=KMR.params$h
  phi=KMR.params$phi
  
  KMR.params.star=KMR.params
  delta.star=delta
  r.star=r
  move.type=2
  if (move.type==1){
    m.choose=ifelse(length(delta)>1,sample(1:length(delta),1),1)
    if (delta[m.choose]==TRUE){
      delta.star[m.choose]=FALSE
      r.star[m.choose]=0
      KMR.params.star$delta=delta.star
      KMR.params.star$r=r.star
      KMR.params.star=update.KMR.params(KMR.params=KMR.params.star,data=data)
    }
    if (delta[m.choose]==FALSE){
      r.star[m.choose]=generate.r.in_deltar.move1.prop(m.choose=m.choose,prop.params=prop.params,r.prior.dist=r.prior.dist)
      delta.star[m.choose]=TRUE
      KMR.params.star$delta=delta.star
      KMR.params.star$r=r.star
      KMR.params.star=update.KMR.params(KMR.params=KMR.params.star,data=data)
    }
    acpt.rate.temp=exp(
      l.deltar.original.diff(data=data,KMR.params=KMR.params,KMR.params.star=KMR.params.star,family=family,PC.num=PC.num,int.num=int.num,int.h.stat=int.h.stat)
      +l.r.in_deltar.prior(m.choose=m.choose,KMR.params=KMR.params.star,prior.params=prior.params,r.prior.dist=r.prior.dist)
      -l.r.in_deltar.prior(m.choose=m.choose,KMR.params=KMR.params,prior.params=prior.params,r.prior.dist=r.prior.dist)
      +l.delta.prior(KMR.params.star,prior.params)-l.delta.prior(KMR.params,prior.params)
      +l.deltar.move1.prop(m.choose=m.choose,current.KMR.params=KMR.params.star,after.KMR.params=KMR.params,prop.params=prop.params,r.prior.dist=r.prior.dist)
      -l.deltar.move1.prop(m.choose=m.choose,current.KMR.params=KMR.params,after.KMR.params=KMR.params.star,prop.params=prop.params,r.prior.dist=r.prior.dist)
    )
    # test.record[t]=l.deltar.original.diff(data=data,KMR.params=KMR.params,KMR.params.star=KMR.params.star,family=family)-
    # l.deltar.original.diff.closedform(data=data,KMR.params=KMR.params,KMR.params.star=KMR.params.star,family=family)
    deltar.acpt.flag=FALSE
    deltar.na.flag=FALSE
    if (!is.nan(acpt.rate.temp)){
      if(runif(1)<acpt.rate.temp){
        delta[m.choose]=delta.star[m.choose]
        r[m.choose]=r.star[m.choose]
        deltar.acpt.flag=TRUE
      }
    }else{
      deltar.na.flag=TRUE
    }
  }
  if (move.type==2){
    m.choose=ifelse(sum(delta)>1,sample(which(delta==TRUE),1),which(delta==TRUE))
    r.star[m.choose]=generate.r.in_deltar.move2.prop(m.choose=m.choose,KMR.params=KMR.params,prop.params=prop.params,r.prior.dist=r.prior.dist)
    KMR.params.star$r=r.star
    KMR.params.star=update.KMR.params(KMR.params=KMR.params.star,data=data)
    acpt.rate.temp=exp(
      l.deltar.original.diff(data=data,KMR.params=KMR.params,KMR.params.star=KMR.params.star,family=family,PC.num=PC.num,int.num=int.num,int.h.stat=int.h.stat)
      +l.r.prior(m.choose=m.choose,KMR.params=KMR.params.star,prior.params=prior.params,r.prior.dist=r.prior.dist)
      -l.r.prior(m.choose=m.choose,KMR.params=KMR.params,prior.params=prior.params,r.prior.dist=r.prior.dist)
      +l.r.in_deltar.move2.prop(m.choose=m.choose,current.KMR.params=KMR.params.star,after.KMR.params=KMR.params,prop.params=prop.params,r.prior.dist=r.prior.dist)
      -l.r.in_deltar.move2.prop(m.choose=m.choose,current.KMR.params=KMR.params,after.KMR.params=KMR.params.star,prop.params=prop.params,r.prior.dist=r.prior.dist)
    )
    
    # test.record[t]=l.deltar.original.diff(data=data,KMR.params=KMR.params,KMR.params.star=KMR.params.star,family=family)-
    #   l.deltar.original.diff.closedform(data=data,KMR.params=KMR.params,KMR.params.star=KMR.params.star,family=family)
    deltar.acpt.flag=FALSE
    deltar.na.flag=FALSE
    if (!is.nan(acpt.rate.temp)){
      if(runif(1)<acpt.rate.temp){
        r[m.choose]=r.star[m.choose]
        deltar.acpt.flag=TRUE
      }
    }else{
      deltar.na.flag=TRUE  
    }
  }
  # print(c('r: ',r))
  KMR.params$r=r
  KMR.params$delta=delta
  KMR.params=update.KMR.params(KMR.params=KMR.params,data=data)
  return (list(KMR.params=KMR.params,deltar.acpt.flag=deltar.acpt.flag,deltar.na.flag=deltar.na.flag,move.type=move.type))
  
}


update.sigmasq<-function(data,KMR.params,prior.params,family){
  if (family != 'gaussian') {stop('only for gaussian distribution')}
  y=data$y
  X=data$X
  Z=data$Z
  Z.diff.precal=data$Z.diff.precal
  beta=KMR.params$beta
  K=KMR.params$K
  A=KMR.params$A
  eig=KMR.params$eig
  tau=KMR.params$tau
  r=KMR.params$r
  delta=KMR.params$delta
  h=KMR.params$h
  phi=KMR.params$phi
  lambda=KMR.params$lambda
  a.sigmasq.prior=prior.params$a.sigmasq.prior
  b.sigmasq.prior=prior.params$b.sigmasq.prior
  
  V.temp=diag(length(y))+lambda*K
  phi=rep(rinvgamma(1,length(y)/2+a.sigmasq.prior,t(y-X%*%beta)%*%solve(V.temp)%*%(y-X%*%beta)/2+b.sigmasq.prior),length(y))
  KMR.params$phi=phi
  return (list(KMR.params=KMR.params))
}

update.lambda<-function(data,KMR.params,prop.params,prior.params,family){
  if (family !=  'gaussian') stop('update.lambda only applies for gaussian')
  y=data$y
  X=data$X
  Z=data$Z
  Z.diff.precal=data$Z.diff.precal
  beta=KMR.params$beta
  K=KMR.params$K
  A=KMR.params$A
  eig=KMR.params$eig
  tau=KMR.params$tau
  r=KMR.params$r
  delta=KMR.params$delta
  h=KMR.params$h
  phi=KMR.params$phi
  lambda=KMR.params$lambda
  lambda.sd.prop=prop.params$lambda.sd.prop
  a.lambda.prior=prior.params$a.lambda.prior
  b.lambda.prior=prior.params$b.lambda.prior
  
  
  lambda.na.flag=FALSE
  lambda.acpt.flag=FALSE # only valid in gaussian
  
  
  # update lambda
  V.temp=diag(length(y))+lambda*K
  lambda.star=rgamma(1,shape=lambda^2/lambda.sd.prop^2,rate=lambda/lambda.sd.prop^2)
  V.star.temp=diag(length(y))+lambda.star*K
  acpt.rate.temp=sqrt(det(V.temp)/det(V.star.temp))*exp( 
    -phi[1]^(-1)*1/2*t(y-X%*%beta)%*%(solve(V.star.temp)-solve(V.temp))%*%(y-X%*%beta)+               #original log-likelihood diff
      dgamma(lambda.star,a.lambda.prior,b.lambda.prior,log=TRUE)-                                     #prior log-likelihood diff
      dgamma(lambda,a.lambda.prior,b.lambda.prior,log=TRUE)+
      dgamma(lambda,shape=lambda.star^2/lambda.sd.prop^2,rate=lambda.star/lambda.sd.prop^2,log=TRUE)- #proposal log-likelihood diff
      dgamma(lambda.star,shape=lambda^2/lambda.sd.prop^2,rate=lambda/lambda.sd.prop^2,log=TRUE)
  )
  if (!is.nan(acpt.rate.temp)){
    if(runif(1)<acpt.rate.temp){
      lambda.acpt.flag=TRUE
      lambda=lambda.star
    }
  }else{
    lambda.na.flag=TRUE
  }
  tau=lambda*phi[1]
  KMR.params$lambda=lambda
  KMR.params$tau=tau
  return (list(KMR.params=KMR.params,lambda.acpt.flag=lambda.acpt.flag,lambda.na.flag=lambda.na.flag))
}


update.deltarh<-function(data,KMR.params,prop.params,prior.params,r.prior.dist,family,PC.num,int.num,int.h.stat){
  y=data$y
  X=data$X
  Z=data$Z
  Z.diff.precal=data$Z.diff.precal
  beta=KMR.params$beta
  K=KMR.params$K
  A=KMR.params$A
  eig=KMR.params$eig
  tau=KMR.params$tau
  r=KMR.params$r
  delta=KMR.params$delta
  h=KMR.params$h
  phi=KMR.params$phi
  
  KMR.params.star=KMR.params
  delta.star=delta
  r.star=r
  move.type=ifelse(sum(delta)==0,1,sample(c(1,2),1))
  if (move.type==1){
    m.choose=ifelse(length(delta)>1,sample(1:length(delta),1),1)
    if (delta[m.choose]==TRUE){
      delta.star[m.choose]=FALSE
      r.star[m.choose]=0
      KMR.params.star$delta=delta.star
      KMR.params.star$r=r.star
      KMR.params.star=update.KMR.params(KMR.params=KMR.params.star,data=data)
      KMR.params.star$h=generate.h.in_detarh.prop(data,KMR.params,family)
    }
    if (delta[m.choose]==FALSE){
      r.star[m.choose]=generate.r.in_deltar.move1.prop(m.choose=m.choose,prop.params=prop.params,r.prior.dist=r.prior.dist)
      delta.star[m.choose]=TRUE
      KMR.params.star$delta=delta.star
      KMR.params.star$r=r.star
      KMR.params.star=update.KMR.params(KMR.params=KMR.params.star,data=data)
      KMR.params.star$h=generate.h.in_detarh.prop(data,KMR.params,family)
    }
    acpt.rate.temp=exp(
      l.deltar.original.diff(data=data,KMR.params=KMR.params,KMR.params.star=KMR.params.star,family=family,PC.num=PC.num,int.num=int.num,int.h.stat=int.h.stat)
      +l.r.in_deltar.prior(m.choose=m.choose,KMR.params=KMR.params.star,prior.params=prior.params,r.prior.dist=r.prior.dist)
      -l.r.in_deltar.prior(m.choose=m.choose,KMR.params=KMR.params,prior.params=prior.params,r.prior.dist=r.prior.dist)
      +l.delta.prior(KMR.params.star,prior.params)-l.delta.prior(KMR.params,prior.params)
      +l.deltar.move1.prop(m.choose=m.choose,current.KMR.params=KMR.params.star,after.KMR.params=KMR.params,prop.params=prop.params,r.prior.dist=r.prior.dist)
      -l.deltar.move1.prop(m.choose=m.choose,current.KMR.params=KMR.params,after.KMR.params=KMR.params.star,prop.params=prop.params,r.prior.dist=r.prior.dist)
      +sum(dnorm(KMR.params.star$h,mean=0,sd=sqrt(tau),log=TRUE)-dnorm(KMR.params$h,mean=0,sd=sqrt(tau),log=TRUE)) #prior for h
      +l.h.in_deltarh.prop(data=data,current.KMR.params=KMR.params.star,after.KMR.params=KMR.params,family)
      -l.h.in_deltarh.prop(data=data,current.KMR.params=KMR.params,after.KMR.params=KMR.params.star,family)
    )
    
  }
  if (move.type==2){
    m.choose=ifelse(sum(delta)>1,sample(which(delta==TRUE),1),which(delta==TRUE))
    r.star[m.choose]=generate.r.in_deltar.move2.prop(m.choose=m.choose,KMR.params=KMR.params,prop.params=prop.params,r.prior.dist=r.prior.dist)
    KMR.params.star$r=r.star
    KMR.params.star=update.h(data,KMR.params.star,family)$KMR.params #
    KMR.params.star$h=generate.h.in_detarh.prop(data,KMR.params,family)
    acpt.rate.temp=exp(
      l.deltar.original.diff(data=data,KMR.params=KMR.params,KMR.params.star=KMR.params.star,family=family,PC.num=PC.num,int.num=int.num,int.h.stat=int.h.stat)
      +l.r.prior(m.choose=m.choose,KMR.params=KMR.params.star,prior.params=prior.params,r.prior.dist=r.prior.dist)
      -l.r.prior(m.choose=m.choose,KMR.params=KMR.params,prior.params=prior.params,r.prior.dist=r.prior.dist)
      +l.r.in_deltar.move2.prop(m.choose=m.choose,current.KMR.params=KMR.params.star,after.KMR.params=KMR.params,prop.params=prop.params,r.prior.dist=r.prior.dist)
      -l.r.in_deltar.move2.prop(m.choose=m.choose,current.KMR.params=KMR.params,after.KMR.params=KMR.params.star,prop.params=prop.params,r.prior.dist=r.prior.dist)
      +sum(dnorm(KMR.params.star$h,mean=0,sd=sqrt(tau),log=TRUE)-dnorm(KMR.params$h,mean=0,sd=sqrt(tau),log=TRUE)) #prior for h
      +l.h.in_deltarh.prop(data=data,current.KMR.params=KMR.params.star,after.KMR.params=KMR.params,family)
      -l.h.in_deltarh.prop(data=data,current.KMR.params=KMR.params,after.KMR.params=KMR.params.star,family)
    )
  }
  
  
  deltar.acpt.flag=FALSE
  deltar.na.flag=FALSE
  if (!is.nan(acpt.rate.temp)){
    if(runif(1)<acpt.rate.temp){
      KMR.params=KMR.params.star
      deltar.acpt.flag=TRUE
    }
  }else{
    deltar.na.flag=TRUE  
  }
  
  return (list(KMR.params=KMR.params,deltar.acpt.flag=deltar.acpt.flag,deltar.na.flag=deltar.na.flag,move.type=move.type))
  
}


generate.h.in_detarh.prop<-function(data,KMR.params,family){
  
  if (family=='binomial'){
    b<-function(theta) log(1+exp(theta))
    b.p<-function(theta){
      exp.theta=exp(theta)
      return (exp.theta/(1+exp.theta))
    }
    b.pp<-function(theta) exp(theta)/(1+exp(theta))^2
    g<-function(mu) log(mu/(1-mu))
    g.p<-function(mu) 1/mu/(1-mu)
  }
  
  if (family=='gaussian'){
    b<-function(theta) theta^2/2
    b.p<-function(theta){
      return (theta)
    } 
    b.pp<-function(theta) return (matrix(rep(1,length(theta)),ncol=1))
    g<-function(mu) mu
    g.p<-function(mu) return (matrix(rep(1,length(theta)),ncol=1))
  }
  
  if (family=='poisson'){
    b<-function(theta) exp(theta)
    b.p<-function(theta){
      return (exp(theta))
    } 
    b.pp<-function(theta) return (exp(theta))
    g<-function(mu) log(mu)
    g.p<-function(mu) return (1/mu)
  }
  
  
  if (family=='binomial-probit'){
    b<-function(theta) log(1+exp(theta))
    b.p<-function(theta){
      exp.theta=exp(theta)
      return (exp.theta/(1+exp.theta))
    }
    b.pp<-function(theta) exp(theta)/(1+exp(theta))^2
    g<-function(mu) qnorm(mu)
    g.p<-function(mu) 1/dnorm(qnorm(mu))
  }
  
  y=data$y
  X=data$X
  Z=data$Z
  Z.diff.precal=data$Z.diff.precal
  beta=KMR.params$beta
  K=KMR.params$K
  A=KMR.params$A
  eig=KMR.params$eig
  tau=KMR.params$tau
  r=KMR.params$r
  delta=KMR.params$delta
  h=KMR.params$h
  phi=KMR.params$phi
  # expfam.params=update.expfam.params(data=data,KMR.params=KMR.params,family=family)
  # eta=expfam.params$eta
  # theta=expfam.params$theta
  # mu=expfam.params$mu
  # y.tilde=expfam.params$y.tilde
  # W.diag=expfam.params$W.diag
  
  
  
  eta=X%*%beta+A%*%h
  theta=eta
  mu=b.p(theta)
  y.tilde=eta+(y-mu)*g.p(mu)
  W.diag=as.vector(1/b.pp(theta)/phi/g.p(mu)^2)
  C.h=(t(W.diag*A))%*%A
  diag(C.h)=diag(C.h)+1/tau
  C.h=solve(C.h)
  m.h=C.h%*%t(A)%*%((y.tilde-X%*%beta)*W.diag)
  return (matrix(rmvnorm(1,mean=m.h,sigma=C.h),ncol=1,nrow=length(y)))
}


l.h.in_deltarh.prop<-function(data,current.KMR.params,after.KMR.params,family){
  
  if (family=='binomial'){
    b<-function(theta) log(1+exp(theta))
    b.p<-function(theta){
      exp.theta=exp(theta)
      return (exp.theta/(1+exp.theta))
    }
    b.pp<-function(theta) exp(theta)/(1+exp(theta))^2
    g<-function(mu) log(mu/(1-mu))
    g.p<-function(mu) 1/mu/(1-mu)
  }
  
  if (family=='gaussian'){
    b<-function(theta) theta^2/2
    b.p<-function(theta){
      return (theta)
    } 
    b.pp<-function(theta) return (matrix(rep(1,length(theta)),ncol=1))
    g<-function(mu) mu
    g.p<-function(mu) return (matrix(rep(1,length(theta)),ncol=1))
  }
  
  if (family=='poisson'){
    b<-function(theta) exp(theta)
    b.p<-function(theta){
      return (exp(theta))
    } 
    b.pp<-function(theta) return (exp(theta))
    g<-function(mu) log(mu)
    g.p<-function(mu) return (1/mu)
  }
  
  
  if (family=='binomial-probit'){
    b<-function(theta) log(1+exp(theta))
    b.p<-function(theta){
      exp.theta=exp(theta)
      return (exp.theta/(1+exp.theta))
    }
    b.pp<-function(theta) exp(theta)/(1+exp(theta))^2
    g<-function(mu) qnorm(mu)
    g.p<-function(mu) 1/dnorm(qnorm(mu))
  }
  
  y=data$y
  X=data$X
  Z=data$Z
  Z.diff.precal=data$Z.diff.precal
  beta=current.KMR.params$beta
  K=current.KMR.params$K
  A=current.KMR.params$A
  eig=current.KMR.params$eig
  tau=current.KMR.params$tau
  r=current.KMR.params$r
  delta=current.KMR.params$delta
  h=current.KMR.params$h
  phi=current.KMR.params$phi
  
  eta=X%*%beta+A%*%h
  theta=eta
  mu=b.p(theta)
  y.tilde=eta+(y-mu)*g.p(mu)
  W.diag=as.vector(1/b.pp(theta)/phi/g.p(mu)^2)
  C.h=(t(W.diag*A))%*%A
  diag(C.h)=diag(C.h)+1/tau
  C.h=solve(C.h)
  m.h=C.h%*%t(A)%*%((y.tilde-X%*%beta)*W.diag)
  
  return (dmvnorm(as.vector(after.KMR.params$h),mean=as.vector(m.h),sigma=C.h,log=TRUE))
  
}



update.ALLh<-function(data,KMR.params,family){
  
  tau=KMR.params$tau
  KMR.params.star=KMR.params
  KMR.params.star$h=generate.h.in_detarh.prop(data,KMR.params,family)
  
  gamma.prob.h=exp(sum(dnorm(KMR.params.star$h,mean=0,sd=sqrt(tau),log=TRUE)-dnorm(KMR.params$h,mean=0,sd=sqrt(tau),log=TRUE)) #prior for h
                   +l.h.in_deltarh.prop(data=data,current.KMR.params=KMR.params.star,after.KMR.params=KMR.params,family)
                   -l.h.in_deltarh.prop(data=data,current.KMR.params=KMR.params,after.KMR.params=KMR.params.star,family))
  
  h.na.flag.vec=FALSE
  h.acpt.flag.vec=FALSE
  
  if (is.nan(gamma.prob.h)){
    h.na.flag.vec=TRUE
  }else{
    if (runif(1)<=gamma.prob.h){
      h.acpt.flag.vec=TRUE
    } 
  }
  return (list(KMR.params=KMR.params,h.na.flag.vec=h.na.flag.vec,h.acpt.flag.vec=h.acpt.flag.vec))
  
}

#sep means update delta, r, and separately h[i]; initial goal is to improve the success rate
update.deltarhSep<-function(data,KMR.params,prop.params,prior.params,r.prior.dist,family,PC.num,int.num,int.h.stat){
  # y=data$y
  # X=data$X
  # Z=data$Z
  # Z.diff.precal=data$Z.diff.precal

  deltar.acpt.flag.vec=rep(FALSE,length(y))
  deltar.na.flag.vec=rep(FALSE,length(y))  
  
  for (i in 1:length(y)){

    
    KMR.params.star=KMR.params
    delta.star=KMR.params.star$delta
    r.star=KMR.params.star$r
    delta=KMR.params$delta
    move.type=ifelse(sum(delta)==0,1,sample(c(1,2),1))
    if (move.type==1){
      m.choose=ifelse(length(delta)>1,sample(1:length(delta),1),1)
      if (delta[m.choose]==TRUE){
        delta.star[m.choose]=FALSE
        r.star[m.choose]=0
        KMR.params.star$delta=delta.star
        KMR.params.star$r=r.star
        KMR.params.star=update.KMR.params(KMR.params=KMR.params.star,data=data)
        KMR.params.star$h=generate.hSep.in_detarh.prop(data,KMR.params,family,i)
      }
      if (delta[m.choose]==FALSE){
        r.star[m.choose]=generate.r.in_deltar.move1.prop(m.choose=m.choose,prop.params=prop.params,r.prior.dist=r.prior.dist)
        delta.star[m.choose]=TRUE
        KMR.params.star$delta=delta.star
        KMR.params.star$r=r.star
        KMR.params.star=update.KMR.params(KMR.params=KMR.params.star,data=data)
        KMR.params.star$h=generate.hSep.in_detarh.prop(data,KMR.params,family,i)
      }
      acpt.rate.temp=exp(
        l.deltar.original.diff(data=data,KMR.params=KMR.params,KMR.params.star=KMR.params.star,family=family,PC.num=PC.num,int.num=int.num,int.h.stat=int.h.stat)
        +l.r.in_deltar.prior(m.choose=m.choose,KMR.params=KMR.params.star,prior.params=prior.params,r.prior.dist=r.prior.dist)
        -l.r.in_deltar.prior(m.choose=m.choose,KMR.params=KMR.params,prior.params=prior.params,r.prior.dist=r.prior.dist)
        +l.delta.prior(KMR.params.star,prior.params)-l.delta.prior(KMR.params,prior.params)
        +l.deltar.move1.prop(m.choose=m.choose,current.KMR.params=KMR.params.star,after.KMR.params=KMR.params,prop.params=prop.params,r.prior.dist=r.prior.dist)
        -l.deltar.move1.prop(m.choose=m.choose,current.KMR.params=KMR.params,after.KMR.params=KMR.params.star,prop.params=prop.params,r.prior.dist=r.prior.dist)
        +dnorm(KMR.params.star$h[i],mean=0,sd=sqrt(tau),log=TRUE)-dnorm(KMR.params$h[i],mean=0,sd=sqrt(tau),log=TRUE) #prior for h
        +l.h.in_deltarhSep.prop(data=data,current.KMR.params=KMR.params.star,after.KMR.params=KMR.params,family=family,index=i)
        -l.h.in_deltarhSep.prop(data=data,current.KMR.params=KMR.params,after.KMR.params=KMR.params.star,family=family,index=i)
      )
      
    }
    if (move.type==2){
      m.choose=ifelse(sum(delta)>1,sample(which(delta==TRUE),1),which(delta==TRUE))
      r.star[m.choose]=generate.r.in_deltar.move2.prop(m.choose=m.choose,KMR.params=KMR.params,prop.params=prop.params,r.prior.dist=r.prior.dist)
      KMR.params.star$r=r.star
      KMR.params.star$h=generate.hSep.in_detarh.prop(data,KMR.params,family,i)
      acpt.rate.temp=exp(
        l.deltar.original.diff(data=data,KMR.params=KMR.params,KMR.params.star=KMR.params.star,family=family,PC.num=PC.num,int.num=int.num,int.h.stat=int.h.stat)
        +l.r.prior(m.choose=m.choose,KMR.params=KMR.params.star,prior.params=prior.params,r.prior.dist=r.prior.dist)
        -l.r.prior(m.choose=m.choose,KMR.params=KMR.params,prior.params=prior.params,r.prior.dist=r.prior.dist)
        +l.r.in_deltar.move2.prop(m.choose=m.choose,current.KMR.params=KMR.params.star,after.KMR.params=KMR.params,prop.params=prop.params,r.prior.dist=r.prior.dist)
        -l.r.in_deltar.move2.prop(m.choose=m.choose,current.KMR.params=KMR.params,after.KMR.params=KMR.params.star,prop.params=prop.params,r.prior.dist=r.prior.dist)
        +dnorm(KMR.params.star$h[i],mean=0,sd=sqrt(tau),log=TRUE)-dnorm(KMR.params$h[i],mean=0,sd=sqrt(tau),log=TRUE) #prior for h
        +l.h.in_deltarhSep.prop(data=data,current.KMR.params=KMR.params.star,after.KMR.params=KMR.params,family=family,index=i)
        -l.h.in_deltarhSep.prop(data=data,current.KMR.params=KMR.params,after.KMR.params=KMR.params.star,family=family,index=i)
      )
    }
    
    

    if (!is.nan(acpt.rate.temp)){
      if(runif(1)<acpt.rate.temp){
        KMR.params=KMR.params.star
        deltar.acpt.flag.vec[i]=TRUE
      }
    }else{
      deltar.na.flag.vec[i]=TRUE  
    }

  }

  return (list(KMR.params=KMR.params,deltar.acpt.flag.vec=deltar.acpt.flag.vec,deltar.na.flag.vec=deltar.na.flag.vec,move.type=move.type))
  
}


generate.hSep.in_detarh.prop<-function(data,KMR.params,family,index){ #update the index'th h
  
  if (family=='binomial'){
    b<-function(theta) log(1+exp(theta))
    b.p<-function(theta){
      exp.theta=exp(theta)
      return (exp.theta/(1+exp.theta))
    }
    b.pp<-function(theta) exp(theta)/(1+exp(theta))^2
    g<-function(mu) log(mu/(1-mu))
    g.p<-function(mu) 1/mu/(1-mu)
  }
  
  if (family=='gaussian'){
    b<-function(theta) theta^2/2
    b.p<-function(theta){
      return (theta)
    } 
    b.pp<-function(theta) return (matrix(rep(1,length(theta)),ncol=1))
    g<-function(mu) mu
    g.p<-function(mu) return (matrix(rep(1,length(theta)),ncol=1))
  }
  
  if (family=='poisson'){
    b<-function(theta) exp(theta)
    b.p<-function(theta){
      return (exp(theta))
    } 
    b.pp<-function(theta) return (exp(theta))
    g<-function(mu) log(mu)
    g.p<-function(mu) return (1/mu)
  }
  
  
  if (family=='binomial-probit'){
    b<-function(theta) log(1+exp(theta))
    b.p<-function(theta){
      exp.theta=exp(theta)
      return (exp.theta/(1+exp.theta))
    }
    b.pp<-function(theta) exp(theta)/(1+exp(theta))^2
    g<-function(mu) qnorm(mu)
    g.p<-function(mu) 1/dnorm(qnorm(mu))
  }
  
  y=data$y
  X=data$X
  Z=data$Z
  Z.diff.precal=data$Z.diff.precal
  beta=KMR.params$beta
  K=KMR.params$K
  A=KMR.params$A
  eig=KMR.params$eig
  tau=KMR.params$tau
  r=KMR.params$r
  delta=KMR.params$delta
  h=KMR.params$h
  phi=KMR.params$phi

  
  i=index
  
  eta=X%*%beta+A%*%h
  theta=eta
  mu=b.p(theta)
  y.tilde=eta+(y-mu)*g.p(mu)
  W.diag=as.vector(1/b.pp(theta)/phi/g.p(mu)^2)
  C.h=(tau^(-1)+sum(A[,i]^2*W.diag))^(-1)
  m.h=C.h*sum(A[,i]*(y.tilde-X%*%beta-A[,-i]%*%h[-i])*W.diag)
  
  h.star=h
  
  h.star[i]=rnorm(1,m.h,sqrt(C.h))
  return (h.star)
}


l.h.in_deltarhSep.prop<-function(data,current.KMR.params,after.KMR.params,family,index){
  
  if (family=='binomial'){
    b<-function(theta) log(1+exp(theta))
    b.p<-function(theta){
      exp.theta=exp(theta)
      return (exp.theta/(1+exp.theta))
    }
    b.pp<-function(theta) exp(theta)/(1+exp(theta))^2
    g<-function(mu) log(mu/(1-mu))
    g.p<-function(mu) 1/mu/(1-mu)
  }
  
  if (family=='gaussian'){
    b<-function(theta) theta^2/2
    b.p<-function(theta){
      return (theta)
    } 
    b.pp<-function(theta) return (matrix(rep(1,length(theta)),ncol=1))
    g<-function(mu) mu
    g.p<-function(mu) return (matrix(rep(1,length(theta)),ncol=1))
  }
  
  if (family=='poisson'){
    b<-function(theta) exp(theta)
    b.p<-function(theta){
      return (exp(theta))
    } 
    b.pp<-function(theta) return (exp(theta))
    g<-function(mu) log(mu)
    g.p<-function(mu) return (1/mu)
  }
  
  
  if (family=='binomial-probit'){
    b<-function(theta) log(1+exp(theta))
    b.p<-function(theta){
      exp.theta=exp(theta)
      return (exp.theta/(1+exp.theta))
    }
    b.pp<-function(theta) exp(theta)/(1+exp(theta))^2
    g<-function(mu) qnorm(mu)
    g.p<-function(mu) 1/dnorm(qnorm(mu))
  }
  
  y=data$y
  X=data$X
  Z=data$Z
  Z.diff.precal=data$Z.diff.precal
  beta=current.KMR.params$beta
  K=current.KMR.params$K
  A=current.KMR.params$A
  eig=current.KMR.params$eig
  tau=current.KMR.params$tau
  r=current.KMR.params$r
  delta=current.KMR.params$delta
  h=current.KMR.params$h
  phi=current.KMR.params$phi
  
  i=index
  
  eta=X%*%beta+A%*%h
  theta=eta
  mu=b.p(theta)
  y.tilde=eta+(y-mu)*g.p(mu)
  W.diag=as.vector(1/b.pp(theta)/phi/g.p(mu)^2)
  C.h=(tau^(-1)+sum(A[,i]^2*W.diag))^(-1)
  m.h=C.h*sum(A[,i]*(y.tilde-X%*%beta-A[,-i]%*%h[-i])*W.diag)
  
  
  return (dnorm(after.KMR.params$h[i],m.h,sqrt(C.h),log=TRUE))
  
}

#check if we can get the correct delta and r given true h,
# just let KMR.params$h=h.true and  KMR.params$tau=tau.true
# note that h.true~N(0,tau K). Some versions of the data generation
# makes h.true~N(0,tau*I)
update.deltar.trueh<-function(data,KMR.params,prop.params,prior.params,r.prior.dist,family,PC.num,int.num,int.h.stat){
  y=data$y
  X=data$X
  Z=data$Z
  Z.diff.precal=data$Z.diff.precal
  beta=KMR.params$beta
  K=KMR.params$K
  A=KMR.params$A
  eig=KMR.params$eig
  tau=KMR.params$tau
  r=KMR.params$r
  delta=KMR.params$delta
  h=KMR.params$h
  phi=KMR.params$phi
  
  KMR.params.star=KMR.params
  delta.star=delta
  r.star=r
  move.type=ifelse(sum(delta)==0,1,sample(c(1,2),1))
  if (move.type==1){
    m.choose=ifelse(length(delta)>1,sample(1:length(delta),1),1)
    if (delta[m.choose]==TRUE){
      delta.star[m.choose]=FALSE
      r.star[m.choose]=0
      KMR.params.star$delta=delta.star
      KMR.params.star$r=r.star
      KMR.params.star=update.KMR.params(KMR.params=KMR.params.star,data=data)
    }
    if (delta[m.choose]==FALSE){
      r.star[m.choose]=generate.r.in_deltar.move1.prop(m.choose=m.choose,prop.params=prop.params,r.prior.dist=r.prior.dist)
      delta.star[m.choose]=TRUE
      KMR.params.star$delta=delta.star
      KMR.params.star$r=r.star
      KMR.params.star=update.KMR.params(KMR.params=KMR.params.star,data=data)
    }
    
    # V.temp=diag(length(y))+tau/phi[1]*KMR.params$K
    # V.star.temp=diag(length(y))+tau/phi[1]*KMR.params.star$K
    # 1/2*log(det(V.temp))-1/2*log(det(V.star.temp))+(
    #   -phi[1]^(-1)*1/2*t(y-X%*%beta)%*%(solve(V.star.temp)-solve(V.temp))%*%(y-X%*%beta)
    # )
    
    acpt.rate.temp=exp(

      dmvnorm(as.vector(h),mean=as.vector(rep(0,length(y))),sigma=KMR.params.star$tau*KMR.params.star$K,log=TRUE)-dmvnorm(as.vector(h),mean=as.vector(rep(0,length(y))),sigma=KMR.params$tau*KMR.params$K,log=TRUE)
      +l.r.in_deltar.prior(m.choose=m.choose,KMR.params=KMR.params.star,prior.params=prior.params,r.prior.dist=r.prior.dist)
      -l.r.in_deltar.prior(m.choose=m.choose,KMR.params=KMR.params,prior.params=prior.params,r.prior.dist=r.prior.dist)
      +l.delta.prior(KMR.params.star,prior.params)-l.delta.prior(KMR.params,prior.params)
      +l.deltar.move1.prop(m.choose=m.choose,current.KMR.params=KMR.params.star,after.KMR.params=KMR.params,prop.params=prop.params,r.prior.dist=r.prior.dist)
      -l.deltar.move1.prop(m.choose=m.choose,current.KMR.params=KMR.params,after.KMR.params=KMR.params.star,prop.params=prop.params,r.prior.dist=r.prior.dist)
    )
    # test.record[t]=l.deltar.original.diff(data=data,KMR.params=KMR.params,KMR.params.star=KMR.params.star,family=family)-
    # l.deltar.original.diff.closedform(data=data,KMR.params=KMR.params,KMR.params.star=KMR.params.star,family=family)
    deltar.acpt.flag=FALSE
    deltar.na.flag=FALSE
    if (!is.nan(acpt.rate.temp)){
      if(runif(1)<acpt.rate.temp){
        delta[m.choose]=delta.star[m.choose]
        r[m.choose]=r.star[m.choose]
        deltar.acpt.flag=TRUE
      }
    }else{
      deltar.na.flag=TRUE
    }
  }
  if (move.type==2){
    m.choose=ifelse(sum(delta)>1,sample(which(delta==TRUE),1),which(delta==TRUE))
    r.star[m.choose]=generate.r.in_deltar.move2.prop(m.choose=m.choose,KMR.params=KMR.params,prop.params=prop.params,r.prior.dist=r.prior.dist)
    KMR.params.star$r=r.star
    KMR.params.star=update.KMR.params(KMR.params=KMR.params.star,data=data)
    
    # V.temp=diag(length(y))+tau/phi[1]*KMR.params$K
    # V.star.temp=diag(length(y))+tau/phi[1]*KMR.params.star$K
    # 1/2*log(det(V.temp))-1/2*log(det(V.star.temp))+(
    #   -phi[1]^(-1)*1/2*t(y-X%*%beta)%*%(solve(V.star.temp)-solve(V.temp))%*%(y-X%*%beta)
    # )
    
    
    acpt.rate.temp=exp(
      dmvnorm(as.vector(h),mean=as.vector(rep(0,length(y))),sigma=KMR.params.star$tau*KMR.params.star$K,log=TRUE)-dmvnorm(as.vector(h),mean=as.vector(rep(0,length(y))),sigma=KMR.params$tau*KMR.params$K,log=TRUE)
      +l.r.prior(m.choose=m.choose,KMR.params=KMR.params.star,prior.params=prior.params,r.prior.dist=r.prior.dist)
      -l.r.prior(m.choose=m.choose,KMR.params=KMR.params,prior.params=prior.params,r.prior.dist=r.prior.dist)
      +l.r.in_deltar.move2.prop(m.choose=m.choose,current.KMR.params=KMR.params.star,after.KMR.params=KMR.params,prop.params=prop.params,r.prior.dist=r.prior.dist)
      -l.r.in_deltar.move2.prop(m.choose=m.choose,current.KMR.params=KMR.params,after.KMR.params=KMR.params.star,prop.params=prop.params,r.prior.dist=r.prior.dist)
    )
    
    # test.record[t]=l.deltar.original.diff(data=data,KMR.params=KMR.params,KMR.params.star=KMR.params.star,family=family)-
    #   l.deltar.original.diff.closedform(data=data,KMR.params=KMR.params,KMR.params.star=KMR.params.star,family=family)
    deltar.acpt.flag=FALSE
    deltar.na.flag=FALSE
    if (!is.nan(acpt.rate.temp)){
      if(runif(1)<acpt.rate.temp){
        r[m.choose]=r.star[m.choose]
        deltar.acpt.flag=TRUE
      }
    }else{
      deltar.na.flag=TRUE  
    }
  }
  # print(c('r: ',r))
  KMR.params$r=r
  KMR.params$delta=delta
  KMR.params=update.KMR.params(KMR.params=KMR.params,data=data)
  return (list(KMR.params=KMR.params,deltar.acpt.flag=deltar.acpt.flag,deltar.na.flag=deltar.na.flag,move.type=move.type))
  
}


flip<-function(delta){ # flip a random number of statuses in a TRUE or FALSE vector
  m.choose=sample(1:length(delta),sample(length(delta),1))
  delta[m.choose]=!delta[m.choose] ###
  return (delta)
}

############
#move 1 in updating delta,r 
generate.r.in_deltar.move1.prop.MultiJump <- function(prop.params,r.prior.dist='invuniform'){
  if (r.prior.dist=='invuniform'){
    a.r.prop = prop.params$a.r.in_deltar.if_invunif.move1.prop
    b.r.prop = prop.params$b.r.in_deltar.if_invunif.move1.prop
    return(1/runif(1,a.r.prop,b.r.prop))
  }
  if (r.prior.dist=='uniform'){
    a.r.prop = prop.params$a.r.in_deltar.if_unif.move1.prop
    b.r.prop = prop.params$b.r.in_deltar.if_unif.move1.prop
    return(runif(1,a.r.prop,b.r.prop))
  }
}


update.deltar.MultiJump<-function(data,KMR.params,prop.params,prior.params,r.prior.dist,family,PC.num,int.num,int.h.stat){
  y=data$y
  X=data$X
  Z=data$Z
  Z.diff.precal=data$Z.diff.precal
  beta=KMR.params$beta
  K=KMR.params$K
  A=KMR.params$A
  eig=KMR.params$eig
  tau=KMR.params$tau
  r=KMR.params$r
  delta=KMR.params$delta
  h=KMR.params$h
  phi=KMR.params$phi
  
  KMR.params.star=KMR.params
  delta.star=delta
  r.star=r
  move.type=ifelse(sum(delta)==0,1,sample(c(1,2),1))
  if (move.type==1){
    delta.star=flip(delta)
    
    #for those from true to false
    index.temp=(delta.star-delta)==-1
    r.star[index.temp]=0
    #for those from false to true
    index.temp=(delta.star-delta)==1
    for (i in which(index.temp==TRUE)){
      r.star[i]=generate.r.in_deltar.move1.prop.MultiJump(prop.params=prop.params,r.prior.dist=r.prior.dist)
    }
    KMR.params.star$delta=delta.star
    KMR.params.star$r=r.star
    KMR.params.star=update.KMR.params(KMR.params=KMR.params.star,data=data)
    
    #move 1 is only for variable selection, no need to sample posterior distribution
    #just need to find the delta to make likelihood largest
    acpt.rate.temp=exp(
      l.deltar.original.diff(data=data,KMR.params=KMR.params,KMR.params.star=KMR.params.star,family=family,PC.num=PC.num,int.num=int.num,int.h.stat=int.h.stat)
      # +l.r.in_deltar.prior(m.choose=m.choose,KMR.params=KMR.params.star,prior.params=prior.params,r.prior.dist=r.prior.dist)
      # -l.r.in_deltar.prior(m.choose=m.choose,KMR.params=KMR.params,prior.params=prior.params,r.prior.dist=r.prior.dist)
      # +l.delta.prior(KMR.params.star,prior.params)-l.delta.prior(KMR.params,prior.params)
      # +l.deltar.move1.prop(m.choose=m.choose,current.KMR.params=KMR.params.star,after.KMR.params=KMR.params,prop.params=prop.params,r.prior.dist=r.prior.dist)
      # -l.deltar.move1.prop(m.choose=m.choose,current.KMR.params=KMR.params,after.KMR.params=KMR.params.star,prop.params=prop.params,r.prior.dist=r.prior.dist)
    )
    deltar.acpt.flag=FALSE
    deltar.na.flag=FALSE
    if (!is.nan(acpt.rate.temp)){
      if(1<acpt.rate.temp){
        delta=delta.star
        r=r.star
        deltar.acpt.flag=TRUE
      }
    }else{
      deltar.na.flag=TRUE
    }
  }
  if (move.type==2){
    m.choose=ifelse(sum(delta)>1,sample(which(delta==TRUE),1),which(delta==TRUE))
    r.star[m.choose]=generate.r.in_deltar.move2.prop(m.choose=m.choose,KMR.params=KMR.params,prop.params=prop.params,r.prior.dist=r.prior.dist)
    KMR.params.star$r=r.star
    KMR.params.star=update.KMR.params(KMR.params=KMR.params.star,data=data)
    acpt.rate.temp=exp(
      l.deltar.original.diff(data=data,KMR.params=KMR.params,KMR.params.star=KMR.params.star,family=family,PC.num=PC.num,int.num=int.num,int.h.stat=int.h.stat)
      +l.r.prior(m.choose=m.choose,KMR.params=KMR.params.star,prior.params=prior.params,r.prior.dist=r.prior.dist)
      -l.r.prior(m.choose=m.choose,KMR.params=KMR.params,prior.params=prior.params,r.prior.dist=r.prior.dist)
      +l.r.in_deltar.move2.prop(m.choose=m.choose,current.KMR.params=KMR.params.star,after.KMR.params=KMR.params,prop.params=prop.params,r.prior.dist=r.prior.dist)
      -l.r.in_deltar.move2.prop(m.choose=m.choose,current.KMR.params=KMR.params,after.KMR.params=KMR.params.star,prop.params=prop.params,r.prior.dist=r.prior.dist)
    )
    
    # test.record[t]=l.deltar.original.diff(data=data,KMR.params=KMR.params,KMR.params.star=KMR.params.star,family=family)-
    #   l.deltar.original.diff.closedform(data=data,KMR.params=KMR.params,KMR.params.star=KMR.params.star,family=family)
    deltar.acpt.flag=FALSE
    deltar.na.flag=FALSE
    if (!is.nan(acpt.rate.temp)){
      if(runif(1)<acpt.rate.temp){
        r[m.choose]=r.star[m.choose]
        deltar.acpt.flag=TRUE
      }
    }else{
      deltar.na.flag=TRUE  
    }
  }
  # print(c('r: ',r))
  KMR.params$r=r
  KMR.params$delta=delta
  KMR.params=update.KMR.params(KMR.params=KMR.params,data=data)
  return (list(KMR.params=KMR.params,deltar.acpt.flag=deltar.acpt.flag,deltar.na.flag=deltar.na.flag,move.type=move.type))
  
}

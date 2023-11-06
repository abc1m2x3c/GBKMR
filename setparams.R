set.start.param<-function(data,family,delta.start){
  y=data$y
  X=data$X
  Z=data$Z
  
  n=length(y)   # sample size
  p=ncol(Z)     # number of candidate Z's
  p0=ncol(X)   # number of X variables

  #pre calculate and store some global variables to speed up the program
  # Z.diff.precal is an n by n by p array that stores the pair-wise difference of Z variables between subject i and j
  Z.diff.precal=array(NA,c(n,n,p))
  for (i in 1:n){
    for (j in 1:i){
      Z.diff.precal[i,j,]=(Z[i,]-Z[j,])^2
      Z.diff.precal[j,i,]=Z.diff.precal[i,j,]
    }
  }
  
  #starting values of beta, K, delta, r, tau, h
  if (family == 'gaussian'){
    if (p0==1){
      start.beta.lm.fit=summary(lm(y~1))
    }else{
      start.beta.lm.fit=summary(lm(y~X[,2:(p0+1)]))
    }
  }
  if (family %in% c('binomial','binomial-probit')){
    if (p0==1){
      start.beta.lm.fit=summary(glm(y~1,family='binomial'))
    }else{
      start.beta.lm.fit=summary(glm(y~X[,2:(p0+1)],family='binomial'))
    }
  }
  if (family == 'poisson'){
    if (p0==1){
      start.beta.lm.fit=summary(glm(y~1,family='poisson'))
    }else{
      start.beta.lm.fit=summary(glm(y~X[,2:(p0+1)],family='poisson'))
    }
  }
  beta.start=start.beta.lm.fit$coefficients[1:p0,1]
  if (family !='gaussian'){
    phi.start=matrix(1,nrow=n,ncol=1) #phi
  }else{
    phi.start=matrix(start.beta.lm.fit$sigma^2,nrow=n,ncol=1) #phi
  }
  

  # if (var.select==FALSE){
  #   delta.start=delta.true # prediction , know delta
  # }
  r.start=rep(0,p)
  r.start[delta.start]=1 # prediction , know delta
  if (sum(delta.start)==0){
    K.start=matrix(1,nrow=n,ncol=n)
  }else if (sum(delta.start)==1){
    K.start=exp(-r.start[delta.start]*Z.diff.precal[,,delta.start])
  }else{
    r.Z.diff.precal.temp=array(NA,c(n,n,sum(delta.start)))
    count=0
    for (i in which(delta.start==TRUE)){
      count=count+1
      r.Z.diff.precal.temp[,,count]=Z.diff.precal[,,i]*r.start[i]
    }
    K.start=exp(-rowSums(r.Z.diff.precal.temp,dims=2))
  }
  temp=eigen(K.start)
  temp$values[(temp$values)<0]=0
  eig.start=temp$values
  A.start=temp$vectors%*%diag(sqrt(temp$values))
  tau.start=10*phi.start[1]
  rm(temp)
  h.start=matrix(0,nrow=n,ncol=1)
  lambda.start=tau.start/phi.start[1]
  
  return (list(beta.start=beta.start,tau.start=tau.start,delta.start=delta.start,
               r.start=r.start,K.start=K.start,A.start=A.start,eig.start=eig.start,
               phi.start=phi.start,h.start=h.start,lambda.start=lambda.start))
  
}

set.prior.param<-function(data){
  y=data$y
  X=data$X
  Z=data$Z
  
  n=length(y)   # sample size
  p=ncol(Z)     # number of candidate Z's
  p0=ncol(X)   # number of X variables
  
  # beta ~ N(mu.beta.prior,V.beta.prior), where mu.beta.prior=0, a p by 1 matrix; V.beta.prior=diag(10^6,p)
  mu.beta.prior=matrix(0,nrow=p0,ncol=1)
  if (p0==1){
    V.beta.prior=10^8
  }else{
    V.beta.prior=diag(rep(10^8,p0))
  }
  a.sigmasq.prior=0.001
  b.sigmasq.prior=0.001
  a.tau.prior=0.001
  b.tau.prior=0.001
  a.lambda.prior=1 # for Gaussian case
  b.lambda.prior=0.1
  # r.prior shared in both move1 and move2
  a.r.if_invunif.prior=0#1/r.true-0.001 # for prior of r
  b.r.if_invunif.prior=100#1/r.true+0.001
  a.r.if_unif.prior=0.01#1/r.true-0.001 # for prior of r
  b.r.if_unif.prior=2  #1/r.true+0.001
  a.r.if_gamma.prior=1
  b.r.if_gamma.prior=0.2
  a.delta.prior=1
  b.delta.prior=1
  
  return(list(mu.beta.prior=mu.beta.prior,V.beta.prior=V.beta.prior,
              a.sigmasq.prior=a.sigmasq.prior,b.sigmasq.prior=b.sigmasq.prior,
              a.tau.prior=a.tau.prior,b.tau.prior=b.tau.prior,
              a.lambda.prior=a.lambda.prior,b.lambda.prior=b.lambda.prior,
              a.r.if_invunif.prior=a.r.if_invunif.prior,b.r.if_invunif.prior=b.r.if_invunif.prior,
              a.r.if_unif.prior=a.r.if_unif.prior,b.r.if_unif.prior=b.r.if_unif.prior,
              a.r.if_gamma.prior=a.r.if_gamma.prior,b.r.if_gamma.prior=b.r.if_gamma.prior,
              a.delta.prior=a.delta.prior,b.delta.prior=b.delta.prior))
  
}


set.prop.param<-function(){
  # Proposal parameters for MH
  lambda.sd.prop=10
  r.in_r.sd.prop=0.1
  a.r.in_r.if_invunif.prop=0 # parameters in solely updating r step, currently no use
  b.r.in_r.if_invunif.prop=100 # currently no use
  r.in_deltar.sd.move1.prop=2  # currently no use
  a.r.in_deltar.if_invunif.move1.prop=0#1/r.true-0.001 # parameters in updating r step
  b.r.in_deltar.if_invunif.move1.prop=100#1/r.true+0.001
  a.r.in_deltar.if_unif.move1.prop=0.001
  b.r.in_deltar.if_unif.move1.prop=2
  r.in_deltar.sd.move2.prop=0.2
  a.r.in_deltar.if_invunif.move2.prop=a.r.in_deltar.if_invunif.move1.prop # parameters in updating r step
  b.r.in_deltar.if_invunif.move2.prop=b.r.in_deltar.if_invunif.move1.prop
  a.r.in_deltar.if_unif.move2.prop=a.r.in_deltar.if_unif.move1.prop
  b.r.in_deltar.if_unif.move2.prop=b.r.in_deltar.if_unif.move1.prop
  return (list(lambda.sd.prop=lambda.sd.prop,r.in_r.sd.prop=r.in_r.sd.prop,
               a.r.in_r.if_invunif.prop=a.r.in_r.if_invunif.prop,
               b.r.in_r.if_invunif.prop=b.r.in_r.if_invunif.prop,
               r.in_deltar.sd.move1.prop=r.in_deltar.sd.move1.prop,
               a.r.in_deltar.if_invunif.move1.prop=a.r.in_deltar.if_invunif.move1.prop,
               b.r.in_deltar.if_invunif.move1.prop=b.r.in_deltar.if_invunif.move1.prop,
               a.r.in_deltar.if_unif.move1.prop=a.r.in_deltar.if_unif.move1.prop,
               b.r.in_deltar.if_unif.move1.prop=b.r.in_deltar.if_unif.move1.prop,
               a.r.in_deltar.if_unif.move2.prop=a.r.in_deltar.if_unif.move1.prop,
               b.r.in_deltar.if_unif.move2.prop=b.r.in_deltar.if_unif.move1.prop,
               r.in_deltar.sd.move2.prop=r.in_deltar.sd.move2.prop,
               a.r.in_deltar.if_invunif.move2.prop=a.r.in_deltar.if_invunif.move2.prop,
               b.r.in_deltar.if_invunif.move2.prop=b.r.in_deltar.if_invunif.move2.prop))
}



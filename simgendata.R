########################################
#data generation begins
########################################
##input
# n: sample size
# p: dimension of Z
# p0: dimension of X
# beta.true: true value for beta
# family: type of outcome ('gaussian','binomial','poisson','binomial-probit')
# pattern: how X\beta+h(Z) is linked to E(Y|X,Z) in the simulation
# 

##output:
# y: outcome
# X: X
# Z: Z
# 
sim.gen.data<-function(n,p,p0,beta.true,family,pattern){
  X=cbind(1,matrix(rnorm(n*(p0-1),0,2),nrow=n,ncol=p0-1))
  Z=matrix(NA,nrow=n,ncol=p)
  for (i in 1:p){
    Z[,i]=runif(n,-1,1)
    delta.true=c(rep(TRUE,3),rep(FALSE,p-3))
  }
  delta.true=c(rep(TRUE,3),rep(FALSE,p-3))
  
  #pattern determines how X and Z are linked to E[Y|X,Z]
  if (pattern == 'nonlinear'){
    eta.true=X%*%beta.true+2*sin(Z[,delta.true][,1]*Z[,delta.true][,2]*pi)+Z[,delta.true][,3]
  }
  if (pattern =='linear'){
    eta.true=X%*%beta.true+Z[,delta.true][,1]+1.5*Z[,delta.true][,2]+2*Z[,delta.true][,3]
  }
  if (pattern =='quadratic'){
    eta.true=X%*%beta.true+0.75*(Z[,delta.true][,1]+Z[,delta.true][,2])^2+Z[,delta.true][,3]
  }
  
  # the distribution of the outcome
  if (family=='gaussian'){
    sigma.true=0.1
    phi.true=rep(sigma.true^2,n)
    y=rnorm(n,eta.true,sigma.true)
  }
  if (family=='binomial'){
    sigma.true=1
    phi.true=rep(sigma.true^2,n)
    y=rbinom(n,1,exp(eta.true)/(1+exp(eta.true)))
  }
  if (family=='poisson')  {
    sigma.true=1
    phi.true=rep(sigma.true^2,n)
    y=rpois(n,exp(eta.true))
  }
  if (family=='binomial-probit'){
    sigma.true=1
    phi.true=rep(sigma.true^2,n)
    y=rbinom(n,1,pnorm(eta.true))
    # ## generate using latent normal representation
    # eps <- rnorm(n)
    # ystar <- X * beta.true + h + eps
    # y <- ifelse(ystar > 0, 1, 0)
    # datp <- list(n = n, M = M, beta.true = beta.true, Z = Z, h = h, X = cbind(X), y = y, eps = eps, ystar = ystar)
  }
  
  
  
  data=list(y=y,X=X,Z=Z)
  
  return (data)
}

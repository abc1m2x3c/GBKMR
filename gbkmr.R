source('src.R')
gbkmr.varselect<-function(data,family,start.params=start.params,prior.params=prior.params,prop.params=prop.params,T=10000,multijump=TRUE,r.prior.dist='invuniform'){
  y=data$y
  X=data$X
  Z=data$Z
  
  n=length(y)   # sample size
  p=ncol(Z)     # number of candidate Z's
  p0=ncol(X)   # number of X variables
  
  if (family=='gaussian'){ 
    int.h.stat=TRUE
  }else{
    int.h.stat=FALSE
  }
  
  #pre calculate and store some global variables to speed up the program
  # Z.diff.precal is an n by n by p array that stores the pair-wise difference of Z variables between subject i and j
  Z.diff.precal=array(NA,c(n,n,p))
  for (i in 1:n){
    for (j in 1:i){
      Z.diff.precal[i,j,]=(Z[i,]-Z[j,])^2
      Z.diff.precal[j,i,]=Z.diff.precal[i,j,]
    }
  }
  data=list(Z=data$Z,X=data$X,y=data$y,Z.diff.precal=Z.diff.precal) #add Z.diff.precal to backward compatible old codes
  
  
  
  # assign starting values
  beta=start.params$beta.start
  tau=start.params$tau.start
  delta=start.params$delta.start
  r=start.params$r.start
  K=start.params$K.start
  A=start.params$A.start
  eig=start.params$eig.start
  phi=start.params$phi.start
  h=start.params$h.start
  lambda=start.params$lambda.start
  KMR.params=list(beta=beta,K=K,A=A,eig=eig,tau=tau,r=r,delta=delta,h=h,phi=phi,lambda=lambda) #lambda is not used when family != gaussian
  
  
  # variables recording each iteration's results
  beta.acpt.record=c()
  beta.na.record=c()
  h.acpt.record=c()
  h.na.record=c()
  tau.acpt.record=c()
  tau.na.record=c()
  move.record=c()
  deltar.na.record=c()
  deltar.acpt.record=c()
  delta.record=c()
  delta.floppos.record=c()
  r.record=c()
  beta.record=c()
  tau.record=c()
  h.new.record=c()
  acpt.deltar.count=0
  for (t in 1:T){
    # if (t%%1==0) (print(t))
    #update beta
    res=update.beta(data=data,KMR.params=KMR.params,prior.params=prior.params,family=family)
    KMR.params=res$KMR.params
    if (t%%10==0){
      beta.record=cbind(beta.record,KMR.params$beta)
      beta.acpt.record=c(beta.acpt.record,res$beta.acpt.flag)
      beta.na.record=c(beta.na.record,res$beta.na.flag)
    }
    # KMR.params$beta=beta.true
    
    #update h
    res=update.h(data=data,KMR.params=KMR.params,family=family)
    KMR.params=res$KMR.params
    if (t%%10==0){
      h.acpt.record=c(h.acpt.record,sum(res$h.acpt.flag.vec)/length(KMR.params$h))
      h.na.record=c(h.na.record,sum(res$h.na.flag.vec)/length(KMR.params$h))
    }
    # KMR.params$h=h.true
    
    #update sigmasq (phi) when gaussian
    if (family=='gaussian'){
      res=update.sigmasq(data,KMR.params,prior.params,family)
      KMR.params=res$KMR.params
      # KMR.params$phi=phi.true
    }
    
    #update tau
    res=update.tau(data=data,KMR.params=KMR.params,prop.params=prop.params,prior.params=prior.params,family=family)
    KMR.params=res$KMR.params
    if (t%%10==0){
      tau.record=c(tau.record,KMR.params$tau)
    }
    # KMR.params$tau=tau.true
    
    #update r and delta together
    # if (var.select==TRUE){
    if (multijump==TRUE){
      res=update.deltar.MultiJump(data=data,KMR.params=KMR.params,prop.params=prop.params,prior.params=prior.params,r.prior.dist=r.prior.dist,family=family,PC.num=PC.num,int.num=int.num,int.h.stat=int.h.stat)  
    }else{
      res=update.deltar(data=data,KMR.params=KMR.params,prop.params=prop.params,prior.params=prior.params,r.prior.dist=r.prior.dist,family=family,PC.num=PC.num,int.num=int.num,int.h.stat=int.h.stat)
    }
    
    # 
    # res=update.deltarhSep(data=data,KMR.params=KMR.params,prop.params=prop.params,prior.params=prior.params,r.prior.dist=r.prior.dist,family=family,PC.num=PC.num,int.num=int.num,int.h.stat=int.h.stat)
    # res=update.deltar.trueh(data=data,KMR.params=KMR.params,prop.params=prop.params,prior.params=prior.params,r.prior.dist=r.prior.dist,family=family,PC.num=PC.num,int.num=int.num,int.h.stat=int.h.stat)
    KMR.params=res$KMR.params
    # }else{
    #   res=update.deltar.fixdelta(data=data,KMR.params=KMR.params,prop.params=prop.params,prior.params=prior.params,r.prior.dist=r.prior.dist,family=family,PC.num=PC.num,int.num=int.num,int.h.stat=int.h.stat)
    #   KMR.params=res$KMR.params
    # }
    if (t%%1==0){
      deltar.acpt.record=c(deltar.acpt.record,res$deltar.acpt.flag)
      if (res$deltar.acpt.flag==TRUE){
        acpt.deltar.count=acpt.deltar.count+1
      } 
      if (res$deltar.acpt.flag==TRUE & res$move.type==1){
        delta.record=cbind(delta.record,KMR.params$delta)
        delta.floppos.record=c(delta.floppos.record,t)
      }
      
      deltar.na.record=c(deltar.na.record,res$deltar.na.flag)
      
      move.record=c(move.record,res$move.type)
      r.record=cbind(r.record,KMR.params$r)
    }
    if (t %%10 ==0){
      print(paste('Iteration ',t,' out of ',T,' completed',sep=''))
      # print(KMR.params$r)
    }
  }
  
  ##recover the delta status in each iteration, and calculate the pip
  
  valid.index=(T/2+1):T #index after burn in process
  temp=delta.floppos.record
  temp.delta=delta.record
  #recover the delta track matrix
  if (temp[1]!=1){
    temp=c(1,temp)
    temp.delta=cbind(start.params$delta.start,temp.delta)
  }
  temp=c(temp,valid.index[length(valid.index)]+1)
  temp.record=c()
  for (j in 1:ncol(temp.delta)){
    temp.record=cbind(temp.record,matrix(rep(temp.delta[,j],diff(temp)[j]),nrow=p))
  }
  
  if (length(valid.index)>1){
    pip=rowMeans(temp.record[,valid.index])
  }
  if (length(valid.index)==1) pip=temp.record[,valid.index]
  
  
  return (list(beta.acpt.record=beta.acpt.record,
               beta.na.record=beta.na.record,
               h.acpt.record=h.acpt.record,
               h.na.record=h.na.record,
               move.record=move.record,
               deltar.na.record=deltar.na.record,
               deltar.acpt.record=deltar.acpt.record,
               delta.record=delta.record,
               delta.floppos.record=delta.floppos.record,
               r.record=r.record,
               beta.record=beta.record,
               tau.record=tau.record,
               pip=pip))
}


gbkmr.pred<-function(data,family,Z.new,start.params,prior.params,prop.params,T,r.prior.dist='invuniform'){
  y=data$y
  X=data$X
  Z=data$Z
  
  n=length(y)   # sample size
  p=ncol(Z)     # number of candidate Z's
  p0=ncol(X)   # number of X variables
  
  if (family=='gaussian'){ 
    int.h.stat=TRUE
  }else{
    int.h.stat=FALSE
  }
  
  #pre calculate and store some global variables to speed up the program
  # Z.diff.precal is an n by n by p array that stores the pair-wise difference of Z variables between subject i and j
  Z.diff.precal=array(NA,c(n,n,p))
  for (i in 1:n){
    for (j in 1:i){
      Z.diff.precal[i,j,]=(Z[i,]-Z[j,])^2
      Z.diff.precal[j,i,]=Z.diff.precal[i,j,]
    }
  }
  data=list(Z=data$Z,X=data$X,y=data$y,Z.diff.precal=Z.diff.precal) #add Z.diff.precal to backward compatible old codes
  
  Z.diff.ON.precal=array(NA,c(n+nrow(Z.new),n+nrow(Z.new),p)) #NO = new and original
  Z.temp=rbind(Z,Z.new)
  for (i in 1:(n+nrow(Z.new))){
    for (j in 1:i){
      Z.diff.ON.precal[i,j,]=(Z.temp[i,]-Z.temp[j,])^2
      Z.diff.ON.precal[j,i,]=Z.diff.ON.precal[i,j,]
    }
  }
  rm(Z.temp)
  
  # assign starting values
  beta=start.params$beta.start
  tau=start.params$tau.start
  delta=start.params$delta.start
  r=start.params$r.start
  K=start.params$K.start
  A=start.params$A.start
  eig=start.params$eig.start
  phi=start.params$phi.start
  h=start.params$h.start
  lambda=start.params$lambda.start
  KMR.params=list(beta=beta,K=K,A=A,eig=eig,tau=tau,r=r,delta=delta,h=h,phi=phi,lambda=lambda) #lambda is not used when family != gaussian
  
  
  # variables recording each iteration's results
  beta.acpt.record=c()
  beta.na.record=c()
  h.acpt.record=c()
  h.na.record=c()
  tau.acpt.record=c()
  tau.na.record=c()
  move.record=c()
  deltar.na.record=c()
  deltar.acpt.record=c()
  delta.record=c()
  delta.floppos.record=c()
  r.record=c()
  beta.record=c()
  tau.record=c()
  h.new.record=c()
  acpt.deltar.count=0
  for (t in 1:T){
    # if (t%%1==0) (print(t))
    #update beta
    res=update.beta(data=data,KMR.params=KMR.params,prior.params=prior.params,family=family)
    KMR.params=res$KMR.params
    if (t%%10==0){
      beta.record=cbind(beta.record,KMR.params$beta)
      beta.acpt.record=c(beta.acpt.record,res$beta.acpt.flag)
      beta.na.record=c(beta.na.record,res$beta.na.flag)
    }
    # KMR.params$beta=beta.true
    
    #update h
    res=update.h(data=data,KMR.params=KMR.params,family=family)
    KMR.params=res$KMR.params
    if (t%%10==0){
      h.acpt.record=c(h.acpt.record,sum(res$h.acpt.flag.vec)/length(KMR.params$h))
      h.na.record=c(h.na.record,sum(res$h.na.flag.vec)/length(KMR.params$h))
    }
    # KMR.params$h=h.true
    
    #update sigmasq (phi) when gaussian
    if (family=='gaussian'){
      res=update.sigmasq(data,KMR.params,prior.params,family)
      KMR.params=res$KMR.params
      # KMR.params$phi=phi.true
    }
    
    #update tau
    res=update.tau(data=data,KMR.params=KMR.params,prop.params=prop.params,prior.params=prior.params,family=family)
    KMR.params=res$KMR.params
    if (t%%10==0){
      tau.record=c(tau.record,KMR.params$tau)
    }
    # KMR.params$tau=tau.true
    
    #update r and delta together
    # if (var.select==TRUE){
    #   if (multijump==TRUE){
    #     res=update.deltar.MultiJump(data=data,KMR.params=KMR.params,prop.params=prop.params,prior.params=prior.params,r.prior.dist=r.prior.dist,family=family,PC.num=PC.num,int.num=int.num,int.h.stat=int.h.stat)  
    #   }else{
    #     res=update.deltar(data=data,KMR.params=KMR.params,prop.params=prop.params,prior.params=prior.params,r.prior.dist=r.prior.dist,family=family,PC.num=PC.num,int.num=int.num,int.h.stat=int.h.stat)
    #   }
    #   
    #   # 
    #   # res=update.deltarhSep(data=data,KMR.params=KMR.params,prop.params=prop.params,prior.params=prior.params,r.prior.dist=r.prior.dist,family=family,PC.num=PC.num,int.num=int.num,int.h.stat=int.h.stat)
    #   # res=update.deltar.trueh(data=data,KMR.params=KMR.params,prop.params=prop.params,prior.params=prior.params,r.prior.dist=r.prior.dist,family=family,PC.num=PC.num,int.num=int.num,int.h.stat=int.h.stat)
    #   KMR.params=res$KMR.params
    # }else{
    res=update.deltar.fixdelta(data=data,KMR.params=KMR.params,prop.params=prop.params,prior.params=prior.params,r.prior.dist=r.prior.dist,family=family,PC.num=PC.num,int.num=int.num,int.h.stat=int.h.stat)
    KMR.params=res$KMR.params
    # }
    if (t%%10==0){
      # print(t)
      # print(KMR.params$r)
      deltar.acpt.record=c(deltar.acpt.record,res$deltar.acpt.flag)
      if (res$deltar.acpt.flag==TRUE){
        acpt.deltar.count=acpt.deltar.count+1
        # print(acpt.deltar.count)
        # print(KMR.params$r)
      } 
      if (res$deltar.acpt.flag==TRUE & res$move.type==1){
        delta.record=cbind(delta.record,KMR.params$delta)
        delta.floppos.record=c(delta.floppos.record,t)
      }
      
      deltar.na.record=c(deltar.na.record,res$deltar.na.flag)
      
      move.record=c(move.record,res$move.type)
      r.record=cbind(r.record,KMR.params$r)
    }
    
    # if (pred==TRUE){
    if (sum(KMR.params$delta)==0){
      K.ON=matrix(1,nrow=n+nrow(Z.new),ncol=n+nrow(Z.new)) 
    }else if (sum(KMR.params$delta)==1){
      K.ON=exp(-KMR.params$r[KMR.params$delta]*Z.diff.ON.precal[,,KMR.params$delta])
    }else{
      
      #reuse  r.Z.diff.precal.temp
      r.Z.diff.precal.temp=array(NA,c(n+nrow(Z.new),n+nrow(Z.new),sum(KMR.params$delta)))
      count=0
      for (i in which(KMR.params$delta==TRUE)){
        count=count+1
        r.Z.diff.precal.temp[,,count]=Z.diff.ON.precal[,,i]*KMR.params$r[i]
      }
      K.ON=exp(-rowSums(r.Z.diff.precal.temp,dims=2))
    }
    
    
    if (sum(KMR.params$delta)==0){
      h.new=rep(0,nrow(Z.new))
    }else{
      K.11=K.ON[(n+1):(n+nrow(Z.new)),(n+1):(n+nrow(Z.new))]
      K.22=K.ON[1:n,1:n]
      K.12=K.ON[(n+1):(n+nrow(Z.new)),1:n]
      K.21=K.ON[1:n,(n+1):(n+nrow(Z.new))]
      if (family == 'gaussian'){
        cov.h.new=KMR.params$tau*(K.11-KMR.params$tau*K.12%*%solve(KMR.params$tau*K.22+KMR.params$phi[1]*diag(n))%*%K.21)
        mu.h.new=KMR.params$tau*K.12%*%solve(KMR.params$tau*K.22+KMR.params$phi[1]*diag(n))%*%(y-X%*%KMR.params$beta)
      }else{
        temp=eigen(KMR.params$K)
        cov.h.new=KMR.params$tau*(K.11-K.12%*% temp$vectors %*% (1/temp$values*t(temp$vectors))  %*%K.21)
        mu.h.new=K.12%*%temp$vectors %*%(sqrt(1/temp$values)*KMR.params$h)
      }
      if (isSymmetric.matrix(cov.h.new)==FALSE){
        cov.h.new[lower.tri(cov.h.new)]=t(cov.h.new)[lower.tri(cov.h.new)]
      } 
      h.new=t(rmvnorm(1,mu.h.new,cov.h.new))
    }
    
    if (t%%10==0){ 
      h.new.record=cbind(h.new.record,h.new)
      print(paste('Iteration ',t,' out of ',T,' completed',sep=''))
    }
    # }
  }
  
  ###
  valid.index=(T/20+1):(T/10)
  h.est=rowMeans(h.new.record[,valid.index])
  h.est=h.est-mean(h.est)
  
  
  return (list(beta.acpt.record=beta.acpt.record,
               beta.na.record=beta.na.record,
               h.acpt.record=h.acpt.record,
               h.na.record=h.acpt.record,
               move.record=move.record,
               deltar.na.record=deltar.na.record,
               deltar.acpt.record=deltar.acpt.record,
               delta.record=delta.record,
               delta.floppos.record=delta.floppos.record,
               r.record=r.record,
               beta.record=beta.record,
               tau.record=tau.record,
               h.new.record=h.new.record,
               h.est=h.est))
  
}
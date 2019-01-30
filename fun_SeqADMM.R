#######################################
#                                     #
# Optimization using ADMM:            #
# Iteratively solve for projection    #
# matrix H.                           #
#                                     #
#######################################

seqADMM_0330=function( S, ndim, PrevPi, PrevPi_d, alpha, lambda, p_m, maxiter,eps,verbose=T){
  p=ncol(S)
  if (alpha==0 & lambda==0){
    S=deflate(S, PrevPi)
    D=eigen(realsym(S))$values
    id=order(D,decreasing=T) #sort eigen value from large to small
    V=eigen(realsym(S))$vectors[,id][,1:ndim]
    projH=V%*%t(V)
  }
  else {
    #starting value
    eps=ndim*eps
    tau=0.1*max(abs(S)) #search step size
    tauStep=2
    
    y0=matrix(0,p,p)
    w0=matrix(0,p,p)
    
    niter=0
    maxnorm=eps+1
    
    while((niter<maxiter)&&(maxnorm>eps)){
      #update h
      h=FantopeProj(y0-w0+S/tau, PrevPi, PrevPi_d, ndim)
      #update y
      
      y1=matrix(0,p,p)
      r=length(p_m)
      for(ell in 1:r){
        for(m in 1:r){
          locRow = ((ell-1)*(p/r)+1):(ell*p/r)
          locCol = ((m-1)*(p/r)+1):(m*p/r)
          yellm = SoftThreshold(h[locRow,locCol] + w0[locRow,locCol],lambda/tau)
          norm_yellm = norm(yellm, type = "F")
          if( norm_yellm > alpha*((p/r)^1)/tau)  {
            y1[locRow,locCol] = (norm_yellm - alpha*((p/r)^1)/tau)*yellm/norm_yellm
          }
        }
      }
      
      
      #update w
      w1=w0+h-y1
      #stop criterion
      normr1=(norm(h-y1,"F"))^2
      norms1=(norm(tau*(y0-y1),"F"))^2
      maxnorm=max(normr1,norms1)
      niter=niter+1
      
      
      #update
      y0=y1
      if (normr1>100*norms1){
        tau=tau*tauStep
        w0=w1/tauStep
      } else if (norms1>100*normr1){
        tau=tau/tauStep
        w0=w1*tauStep
      } else {w0=w1}
    }
    
    if (verbose){
      #display warning message
      if (niter<maxiter){
        warning(paste("seqADMM has converged after",niter,"iterations"))
      } else {
        warning("seqADMM could not converge")
      }
    }
    
    projH=y1
  }
  return(projH)
}


#######################################
#                                     #
# First step in ADMM: get h           #
#                                     #
#######################################

GetTheta=function(v,ndim){
  if ( (v[ndim]-v[ndim+1])>1 || (v[ndim]-v[ndim+1])==1 ){
    theta=v[ndim]-1
    return(theta)
    break
  }
  p=length(v)
  v1=1:(p+1)
  v1[1:p]=v
  v1[p+1]=v[p]-ndim/p
  ddnew=0
  fnew=0
  dnew=max(ndim-2,0)
  while (fnew<ndim){
    f=fnew
    dd=ddnew
    d=dnew
    dnew=dnew+1
    theta=v1[dnew]
    ddnew=which((v1-theta)<1)[1]
    fnew=(ddnew-1)+sum(v1[ddnew:dnew])-(dnew-ddnew+1)*theta
  }
  if (fnew==ndim){
    return(theta)
    break
  }
  theta=v1[d]
  m0=min(1-(v1[dd]-theta),theta-v1[d+1])
  while ((f+(d-dd+1)*m0)<ndim){
    f=f+(d-dd+1)*m0
    dd=dd+1
    theta=theta-m0
    m0=min(1-(v1[dd]-theta),theta-v1[d+1])
  }
  theta=theta-(ndim-f)/(d-dd+1)
}


FantopeProj=function(mat,mat0,d,ndim){
  #If we have previous projection matrix pi, we can find its orthogonal complement basis U, and multiply U to mat to get new mat.
  if (!is.null(mat0)){
    p=dim(mat0)[1]
    D=eigen(diag(p)-mat0)$values
    U=eigen(diag(p)-mat0)$vectors
    id=order(D,decreasing=T) #sort eigen value from large to small
    U=U[,id]
    U=U[,1:(p-d)]
    mat=t(U)%*%mat%*%U
  }
  #Decompose mat as in spectral decomposition, reform eigenvalue, then form mat back.
  mat=Re((mat+t(mat))/2)
  D=eigen(mat)$values
  V=eigen(mat)$vectors
  id=order(D,decreasing = T) #sort eigen value from large to small
  D=D[id]
  V=V[,id]
  theta=GetTheta(D,ndim)
  new.values=pmin(pmax(D-theta,0),1)
  newmat=V%*%diag(new.values)%*%t(V)
  #If we have previous projection matrix pi, we again multiply newmat by U.
  if (!is.null(mat0)){
    newmat=U%*%newmat%*%t(U)
  }
  return(newmat)
}

#######################################
#                                     #
# Second step in ADMM: get y          #
#                                     #
#######################################

SoftThreshold=function(x,lambda){
  newvalue=sign(x)*pmax(abs(x)-lambda,0) #pmax is pointwise max
  return(newvalue)
}
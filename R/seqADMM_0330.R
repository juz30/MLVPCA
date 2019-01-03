


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


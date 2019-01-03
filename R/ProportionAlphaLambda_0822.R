


ProportionAlphaLambda_0822=function(S,xcov,k,PrevPi,PrevPi_d,alphaSeq,lambdaSeq, p_m, op, totV){
  # for center X
  nalpha=length(alphaSeq)
  nlambda=length(lambdaSeq)
  FVE=matrix(0, nrow = nalpha, ncol = nlambda)

  for (i in 1:nalpha){
    for (j in 1:nlambda){
      alpha=alphaSeq[i]
      lambda=lambdaSeq[j]
      projH=seqADMM_0330(S,k,PrevPi,PrevPi_d, alpha, lambda, p_m, op$maxiter,op$eps)

      eigV=eigen(realsym(projH))$vectors[,1] #this is pxk

      FVE[i,j]=sum(diag(t(eigV)%*%xcov%*%eigV))/totV
    }
  }

  prop=FVE/FVE[1,1]
  print(prop)
  #Method 1:
  #Among candidate who can explain enough variance,choose the combination of alpha+lambda that yields the smallest rFVE to localize more.
  #I=which(prop==min(prop[prop>op$alphaproportion]),arr.ind = T)

  #Method 2.1:
  #If multiple largest combo of lambda+alpha, choose the one with largest lambda,
  #because alpha doesn't have to be the largest to achive maximum influence.
  # if (is.null(dim(eligible1))){ #if there is only one row in eligible1
  #   I=eligible1
  # }else{
  #   I=tail(eligible1,1)
  # }

  #Method 2.2:
  #If multiple largest combo of lambda+alpha, choose the one with largest alpha,
  #because we would like to see more variate being shrunken to 0.
  if (is.null(dim(eligible1))){ #if there is only one row in eligible1
    I=eligible1
  }else{
    I=eligible1[which.max(eligible1[1,]),] #choose the combo with largest alpha
  }

  cat("Selected index",I)

  alpha=alphaSeq[I[1]]
  lambda=lambdaSeq[I[2]]

  return(c(alpha,lambda))

}

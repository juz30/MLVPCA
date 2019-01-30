####################################
#                                  #
#  Generate a sequence of rho1     #
#                                  #
####################################

GenerateRho1=function(S,m,nsol,minv,maxv,skew){
  p=ncol(S)
  if (is.null(minv)){
    minv=0
  }
  if (is.null(maxv)){
    l1=eigen(realsym(S))$values[1]
    maxv=l1*p/m
  }
  if (is.null(skew)){
    skew=1
  }
  step=1/(nsol-1)
  lx=(1-exp(skew*seq(0,1,by=step)))/(1-exp(skew))
  rho1Seq=(1-lx)*minv+lx*maxv
  return(rho1Seq)
}

############################################
#                                          #
#   Select rho1 through cross-validation   #
#                                          #
############################################


SplitRho1_list=function(X,op){
  
  
  if (op$model=="2WayNested"|op$model=="2WayNested.spatial"){
    
    L1=op$levels[1]
    L2=op$levels[2]
    nsplit=floor(L1/op$nfold)
    # starting value
    Vsplit.w=0
    Vsplit.z=0
    i1.w=0
    i1.z=0
    nrho1=length(op$rho1Seq.list$rho1Seq.w) # WE REQUEST SEQUENCE OF Z AND W  ARE THE SAME.
    
    # To get CV <S,H> for each of the 10 candidate rho, using parallel computing or not
    
    if (option$parallel==T){
      
      # parallel computing
      no_cores=op$no_cores
      cl=makeCluster(no_cores,type="FORK")
      
      EachElement=function(x){
        rho1.w=op$rho1Seq.list$rho1Seq.w[x]
        rho1.z=op$rho1Seq.list$rho1Seq.z[x]
        V.w=0
        V.z=0
        # For each rho, do cv nfold times, and aggregate test tr(t(H)*S) to get V. Then compare V for every rho
        # Split data by L1 (subject), if has more levels, split the data by the uppest level (to make the data structure complete within train and test)
        for (ifold in 1:op$nfold){
          idtest=((ifold-1)*nsplit*L2+1):(ifold*nsplit*L2)
          Xtrain=X[-idtest,]
          Xtest=X[idtest,]
          G.train=MultilevelS(Xtrain,L1-nsplit,L2,op)
          G.test=MultilevelS(Xtest,nsplit,L2,op)
          
          # Within subject
          K.train.w=G.train$G.w-rho1.w*op$SmoothD
          projH.w=seqADMM_0330(K.train.w, 1, op$PrevPi, op$PrevPi_d,0,0)
          K.test.w=G.test$G.w
          V.w=V.w+sum(projH.w*K.test.w)
          # Between subject
          K.train.z=G.train$G.z-rho1.z*op$SmoothD
          projH.z=seqADMM_0330(K.train.z, 1, op$PrevPi, op$PrevPi_d,0,0)
          K.test.z=G.test$G.z
          V.z=V.z+sum(projH.z*K.test.z)
        }
        cat("current total sum",V.w, V.z, "\t")
        return(c(V.w,V.z))
      }
      
      Vmatrix=matrix(parSapply(cl, 1:nrho1, EachElement),byrow = F, ncol=nrho1, nrow=2)
      stopCluster(cl)
      
      i1.w=which.max(Vmatrix[1,])
      i1.z=which.max(Vmatrix[2,])
      
    } else{
      
      for (i in 1:nrho1){
        
        rho1.w=op$rho1Seq.list$rho1Seq.w[i]
        rho1.z=op$rho1Seq.list$rho1Seq.z[i]
        V.w=0
        V.z=0
        # For each rho, do cv nfold times, and aggregate test tr(t(H)*S) to get V. Then compare V for every rho
        # Split data by L1 (subject), if has more levels, split the data by the uppest level (to make the data structure complete within train and test)
        for (ifold in 1:op$nfold){
          idtest=((ifold-1)*nsplit*L2+1):(ifold*nsplit*L2)
          Xtrain=X[-idtest,]
          Xtest=X[idtest,]
          G.train=MultilevelS(Xtrain,L1-nsplit,L2,op)
          G.test=MultilevelS(Xtest,nsplit,L2,op)
          
          # Within subject
          K.train.w=G.train$G.w-rho1.w*op$SmoothD
          projH.w=seqADMM_0330(K.train.w, 1, op$PrevPi, op$PrevPi_d,0,0)
          K.test.w=G.test$G.w
          V.w=V.w+sum(projH.w*K.test.w)
          # Between subject
          K.train.z=G.train$G.z-rho1.z*op$SmoothD
          projH.z=seqADMM_0330(K.train.z, 1, op$PrevPi, op$PrevPi_d,0,0)
          K.test.z=G.test$G.z
          V.z=V.z+sum(projH.z*K.test.z)
        }
        cat("current total sum",V.w, V.z, "\t")
        
        # If the V for the current rho is larger than the V for the previous rho, we replace Vsplit and il (to track down V and i).
        if (V.w>=Vsplit.w){
          Vsplit.w=V.w
          i1.w=i
        }
        if (V.z>=Vsplit.z){
          Vsplit.z=V.z
          i1.z=i
        }
      }
      
      
    }
    
    # Select best rho based on i1.w and i1.z
    
    if (i1.w<nrho1){
      rho_1.w=op$rho1Seq.list$rho1Seq.w[i1.w+1]
    }else {
      rho_1.w=op$rho1Seq.list$rho1Seq.w[i1.w]
    }
    
    if (i1.z<nrho1){
      rho_1.z=op$rho1Seq.list$rho1Seq.z[i1.z+1]
    }else{
      rho_1.z=op$rho1Seq.list$rho1Seq.z[i1.z]
    }
    
    # rho_1.w=op$rho1Seq.list$rho1Seq.w[i1.w]
    # rho_1.z=op$rho1Seq.list$rho1Seq.z[i1.z]
    
    return(rho_1=list("rho_1.w"=rho_1.w, "rho_1.z"=rho_1.z))
    
  }
  
  
}

############################################
#                                          #
# Generate sequences  of alpha and lambda  #
#                                          #
############################################


GenerateRho2=function(S,nsol,minv,maxv,skew){
  if (is.null(minv)){
    minv=0
  }
  if (is.null(maxv)){
    p=nrow(S)
    v=abs(S)+diag(rep(NaN,p))
    maxv=quantile(as.vector(v),0.95, na.rm = T)
  }
  if (maxv<10){
    maxv=maxv*4 ## this is newly added, 4/26
  }
  
  if (is.null(skew)){
    skew=1
  }
  
  step0=1/(nsol-1)
  lx=(1-exp(skew*seq(0,1,by=step0)))/(1-exp(skew))
  rho2Seq=(1-lx)*minv+lx*maxv
  return(rho2Seq)
}

#######################################################
#                                                     #
#  Select alpha and lambda through cross-validation   #
#                                                     #
#######################################################


CV_2WayNested_1020=function(x_c, k, PrevPi.w,PrevPi.z, PrevPi_d, alphaSeq.w, lambdaSeq.w, alphaSeq.z, lambdaSeq.z, p_m, option){
  # for center x
  nalpha=length(alphaSeq.z) # WE REQUEST SEQUENCE OF Z AND W  ARE THE SAME.
  nlambda=length(lambdaSeq.z)
  n=nrow(x_c)
  p=ncol(x_c)
  L1=option$levels[1]
  L2=option$levels[2]
  nfold=option$nfold
  nsplit=floor(L1/nfold)
  
  #alpha0.w=0
  alpha0.w=alphaSeq.w[2]
  lambda0.w=0
  #alpha0.z=0
  alpha0.z=alphaSeq.z[2]
  lambda0.z=0
  cont.w=1
  cont.z=1
  maxiter=20
  niter=1
  
  Each=function(i){
    if (select=="lambda"){
      alpha.w=alpha0.w
      lambda.w=lambdaSeq.w[i]
      alpha.z=alpha0.z
      lambda.z=lambdaSeq.z[i]
    }else if(select=="alpha"){
      alpha.w=alphaSeq.w[i]
      lambda.w=lambda0.w
      alpha.z=alphaSeq.z[i]
      lambda.z=lambda0.z
    }
    #cat("in function lambda.z=",lambda.z,"\t")
    #cat("in function lambda.w=",lambda.w,"\t")
    #cat("in function alpha.z=",alpha.z,"\t")
    #cat("in function alpha.w=",alpha.w,"\t")
    
    V.w=0
    V.z=0
    for (ifold in 1:nfold){
      idtest=((ifold-1)*nsplit*L2+1):(ifold*nsplit*L2)
      Xtrain=x_c[-idtest,]
      Xtest=x_c[idtest,]
      G.train=MultilevelS(Xtrain,L1-nsplit,L2,option)
      G.test=MultilevelS(Xtest,nsplit,L2,option)
      
      # within subject
      if (cont.w>0){
        K.train.w=G.train$G.w-option$rho_1.list$rho_1.w*option$SmoothD
        projH.w=seqADMM_0330(K.train.w, k, PrevPi.w, PrevPi_d, alpha.w, lambda.w, p_m, option$maxiter, option$eps)
        K.test.w=G.test$G.w
        V.w=V.w+sum(projH.w*K.test.w)
        #cat("fold i ",sum(projH.w*K.test.w))
      }
      
      # between subject
      if (cont.z>0){
        K.train.z=G.train$G.z-option$rho_1.list$rho_1.z*option$SmoothD
        projH.z=seqADMM_0330(K.train.z, k, PrevPi.z, PrevPi_d, alpha.z, lambda.z, p_m, option$maxiter, option$eps)
        K.test.z=G.test$G.z
        V.z=V.z+sum(projH.z*K.test.z)
        
      }
    }
    Vsplit1.w=V.w/nfold
    Vsplit1.z=V.z/nfold
    return(c(Vsplit1.z,Vsplit1.w))
  }
  
  
  
  if (is.null(option$alpha.list)&is.null(option$lambda.list)){
    ############### if both alpha and lambda need to be selected
    while ((cont.w>0 | cont.z>0)&(niter<maxiter)){
      if (niter %% 2==1){ # Now select lambda
        select="lambda"
        if (option$parallel==T){
          # parallel computing
          no_cores=option$no_cores
          cl=makeCluster(no_cores,type="FORK")
          Vsplit1=parSapply(cl, 1:nlambda, Each)
          stopCluster(cl)
        } else{
          Vsplit1=sapply(1:nlambda, Each)
        }
        
        I.z=which(Vsplit1[1,]==max(Vsplit1[1,]))
        if (cont.z>0){
          lambda1.z=lambdaSeq.z[I.z]
          cat("lambda.z=",lambda1.z,"maxCV.z=",Vsplit1[1,I.z],"\t")
        }else {
          lambda1.z=lambda0.z # now z is already converged, don't update lambda.z
        }
        
        I.w=which(Vsplit1[2,]==max(Vsplit1[2,]))
        if (cont.w>0){
          lambda1.w=lambdaSeq.w[I.w]
          cat("lambda.w=",lambda1.w,"maxCV.w=",Vsplit1[2,I.w],"\t")
        }else{
          lambda1.w=lambda0.w # now w is already converged, don't update lambda.w
        }
        
        ### stop criteria
        if ( lambda1.w==lambda0.w & niter>1 ){
          cont.w=0
          cat(" w Converged")
        }
        if ( lambda1.z==lambda0.z & niter>1 ){
          cont.z=0
          cat(" z Converged")
        }
        ### update lambda0
        lambda0.w=lambda1.w
        lambda0.z=lambda1.z
        niter=niter+1
        
      }else if (niter %% 2==0){  # Now select alpha
        select="alpha"
        if (option$parallel==T){
          no_cores=option$no_cores
          cl=makeCluster(no_cores,type="FORK")
          Vsplit2=parSapply(cl, 1:nalpha, Each)
          stopCluster(cl)
        }else{
          Vsplit2=sapply(1:nalpha,Each)
        }
        
        I.z=which(Vsplit2[1,]==max(Vsplit2[1,]))
        if (cont.z>0){
          alpha1.z=alphaSeq.z[I.z]
          cat("alpha.z=",alpha1.z,"maxCV.z=",Vsplit2[1,I.z],"\t")
        }else{
          alpha1.z=alpha0.z  # now z is already converged, don't update alpha.z
        }
        
        I.w=which(Vsplit2[2,]==max(Vsplit2[2,]))
        if (cont.w>0){
          alpha1.w=alphaSeq.w[I.w]
          cat("alpha.w=",alpha1.w,"maxCV.w=",Vsplit2[2,I.w],"\t")
        }else{
          alpha1.w=alpha0.w  # now w is already converged, don't update alpha.w
        }
        
        ### stop criteria
        if ( alpha1.w==alpha0.w & niter>1 ){
          cont.w=0
          cat(" w Converged")
        }
        if ( alpha1.z==alpha0.z & niter>1 ){
          cont.z=0
          cat(" z Converged")
        }
        ### update alpha0
        alpha0.w=alpha1.w
        alpha0.z=alpha1.z
        niter=niter+1
      }
    }
    
  } else if (!is.null(option$alpha.list)&is.null(option$lambda.list)){ # if only select lambda
    ########################################## if only select lambda
    select="lambda"
    ki=PrevPi_d+1
    alpha0.w=option$alpha.list$alpha.w[ki]
    alpha0.z=option$alpha.list$alpha.z[ki]
    if (option$parallel==T){
      # parallel computing
      no_cores=option$no_cores
      cl=makeCluster(no_cores,type="FORK")
      Vsplit1=parSapply(cl, 1:nlambda, Each)
      stopCluster(cl)
    } else{
      Vsplit1=sapply(1:nlambda, Each)
    }
    I.z=which(Vsplit1[1,]==max(Vsplit1[1,]))
    lambda1.z=lambdaSeq.z[I.z]
    cat("lambda.z=",lambda1.z,"\t")
    I.w=which(Vsplit1[2,]==max(Vsplit1[2,]))
    lambda1.w=lambdaSeq.w[I.w]
    cat("lambda.w=",lambda1.w,"\t")
    
    alpha1.z=alpha0.z
    alpha1.w=alpha0.w
    
  }else if (is.null(option$alpha.list)&!is.null(option$lambda.list)){
    ########################################## if only select alpha
    select="alpha"
    ki=PrevPi_d+1
    lambda0.w=option$lambda.list$lambda.w[ki]
    lambda0.z=option$lambda.list$lambda.z[ki]
    if (option$parallel==T){
      # parallel computing
      no_cores=option$no_cores
      cl=makeCluster(no_cores,type="FORK")
      Vsplit2=parSapply(cl, 1:nalpha, Each)
      stopCluster(cl)
    } else{
      Vsplit2=sapply(1:nalpha, Each)
    }
    I.z=which(Vsplit2[1,]==max(Vsplit2[1,]))
    alpha1.z=alphaSeq.z[I.z]
    cat("alpha.z=",alpha1.z,"\t")
    I.w=which(Vsplit2[2,]==max(Vsplit2[2,]))
    alpha1.w=alphaSeq.w[I.w]
    cat("alpha.w=",alpha1.w,"\t")
    
    lambda1.z=lambda0.z
    lambda1.w=lambda0.w
  }
  
  
  
  return(list('alpha1.w'=alpha1.w,'lambda1.w'=lambda1.w,'alpha1.z'=alpha1.z,'lambda1.z'=lambda1.z))
  
}


######################################################
#                                                    #
#    Select alpha and lambda through FVE method      #
#                                                    #
######################################################


ProportionAlphaLambda_0829=function(S,xcov,k,PrevPi,PrevPi_d,alphaSeq,lambdaSeq, p_m, op, totV){
  # for center X
  nalpha=length(alphaSeq)
  nlambda=length(lambdaSeq)
  FVE=matrix(0, nrow = nalpha, ncol = nlambda)
  
  
  # parallel computing
  no_cores=op$no_cores
  cl=makeCluster(no_cores,type="FORK")
  
  EachElement=function(x){
    index=cbind(rep(1:nalpha,rep(nlambda,nalpha)), rep(1:nlambda,nalpha))
    i=index[x,][1]
    j=index[x,][2]
    alpha=alphaSeq[i]
    lambda=lambdaSeq[j]
    projH=seqADMM_0330(S,k,PrevPi,PrevPi_d, alpha, lambda, p_m, op$maxiter,op$eps)
    
    eigV=eigen(realsym(projH))$vectors[,1] #this is pxk
    
    FVE[i,j]=sum(diag(t(eigV)%*%xcov%*%eigV))/totV
    
    return(FVE[i,j])
  }
  
  FVE=matrix(parSapply(cl, 1:(nalpha*nlambda), EachElement),byrow = T, ncol=nlambda, nrow=nalpha)
  
  stopCluster(cl)
  
  prop=FVE/FVE[1,1]
  print(prop)
  
  #Method 1:
  #Among candidate who can explain enough variance,choose the combination of alpha+lambda that yields the smallest rFVE to localize more.
  #I=which(prop==min(prop[prop>op$alphaproportion]),arr.ind = T)
  
  #Method 2:
  #Among candidate who can explain enough variance,choose the largest alpha+lambda to localize more.
  eligible=which(prop>=op$alphaproportion,arr.ind = T)
  eligible1=eligible[which(apply(eligible,1,sum)==max(apply(eligible,1,sum))),]
  
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


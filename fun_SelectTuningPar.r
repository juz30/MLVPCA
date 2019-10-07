#' @import parallel
#' @importFrom stats quantile

####################################
#                                  #
#  Generate a sequence of gamma     #
#                                  #
####################################

Generategamma=function(S,m,nsol,minv,maxv,skew){
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
  gammaSeq=(1-lx)*minv+lx*maxv
  return(gammaSeq)
}


############################################
#                                          #
#   Select gamma through cross-validation   #
#                                          #
############################################

Splitgamma_list=function(X, option){


    gammaSeq.list <- option$gammaSeq.list # is the candidate rho sequence to perform CV
    L1 <- option$L1 # is used to split data for CV
    L2 <- option$L2 # is used to split data for CV
    nfold <- option$nfold # nfold CV
    SmoothD <- option$SmoothD # SmoothD can be used for each splited dataset (same dimension)
    parallel <- option$parallel
    no_cores <- option$no_cores


    nsplit=floor(L1/nfold)
    # starting value
    Vsplit.w=0
    Vsplit.z=0
    i1.w=0
    i1.z=0
    ngamma=length(gammaSeq.list$gammaSeq.w) # WE REQUEST SEQUENCE OF Z AND W  ARE THE SAME.

    # For each of the 10 candidate rho perform CV and get CV <S,H>

    EachElement=function(x){
      gamma.w=gammaSeq.list$gammaSeq.w[x]
      gamma.z=gammaSeq.list$gammaSeq.z[x]
      V.w=0
      V.z=0
      # For each rho, do cv nfold times, and aggregate test tr(t(H)*S) to get V. Then compare V for every rho
      # Split data by L1 (subject), if has more levels, split the data by the uppest level (to make the data structure complete within train and test)
      for (ifold in 1:nfold){
        idtest=((ifold-1)*nsplit*L2+1):(ifold*nsplit*L2)
        Xtrain=X[-idtest,]
        Xtest=X[idtest,]
        G.train=MultilevelS(Xtrain,L1-nsplit,L2, option)
        G.test=MultilevelS(Xtest,nsplit,L2, option)

        # Within subject
        K.train.w=G.train$G.w-gamma.w*SmoothD
        projH.w=seqADMM_0330(K.train.w, 1, option$PrevPi, option$PrevPi_d,0,0) # gamma is selected based on PC1, PrevPi=NULL
        K.test.w=G.test$G.w
        V.w=V.w+sum(projH.w*K.test.w)
        # Between subject
        K.train.z=G.train$G.z-gamma.z*SmoothD
        projH.z=seqADMM_0330(K.train.z, 1, option$PrevPi, option$PrevPi_d,0,0)
        K.test.z=G.test$G.z
        V.z=V.z+sum(projH.z*K.test.z)
      }
      return(c(V.w,V.z))
    }

    # Use parallel computing or sapply to get CV <S,H> for 10 candidate rhos.

    if (parallel==T){

      # parallel computing
      cl=makeCluster(no_cores,type="FORK")

      Vmatrix=matrix(parSapply(cl, 1:ngamma, EachElement),byrow = F, ncol=ngamma, nrow=2)
      stopCluster(cl)

      i1.w=which.max(Vmatrix[1,])
      i1.z=which.max(Vmatrix[2,])

    } else{

      Vmatrix=matrix(sapply(1:ngamma, EachElement),byrow = F, ncol=ngamma, nrow=2)
      i1.w=which.max(Vmatrix[1,])
      i1.z=which.max(Vmatrix[2,])

      }


    # Select best rho based on largest CV <S,H> ( V.w and V.z ).

    if (i1.w<ngamma){
      gamma.w=gammaSeq.list$gammaSeq.w[i1.w+1]
    }else {
      gamma.w=gammaSeq.list$gammaSeq.w[i1.w]
    }

    if (i1.z<ngamma){
      gamma.z=gammaSeq.list$gammaSeq.z[i1.z+1]
    }else{
      gamma.z=gammaSeq.list$gammaSeq.z[i1.z]
    }


    return(gamma=list("gamma.w"=gamma.w, "gamma.z"=gamma.z))


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

CV_2WayNested_1020=function(x_c, alphaSeq.w, lambdaSeq.w, alphaSeq.z, lambdaSeq.z,
                            Fantope_d, PrevPi.w, PrevPi.z, PrevPi_d, p_m, option){
  # for center x
  nalpha=length(alphaSeq.z) # WE REQUEST SEQUENCE OF Z AND W  ARE THE SAME.
  nlambda=length(lambdaSeq.z)
  n=nrow(x_c)
  p=ncol(x_c)
  L1=option$L1
  L2=option$L2
  nfold=option$nfold
  nsplit=floor(L1/nfold)


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
        K.train.w=G.train$G.w-option$gamma.list$gamma.w*option$SmoothD
        projH.w=seqADMM_0330(K.train.w, Fantope_d, PrevPi.w, PrevPi_d, alpha.w, lambda.w, p_m, option$maxiter, option$eps)
        K.test.w=G.test$G.w
        V.w=V.w+sum(projH.w*K.test.w)
      }

      # between subject
      if (cont.z>0){
        K.train.z=G.train$G.z-option$gamma.list$gamma.z*option$SmoothD
        projH.z=seqADMM_0330(K.train.z, Fantope_d, PrevPi.z, PrevPi_d, alpha.z, lambda.z, p_m, option$maxiter, option$eps)
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

        I.z=which.max(Vsplit1[1,])
        if (cont.z>0){
          lambda1.z=lambdaSeq.z[I.z]
          #cat("lambda.z=",lambda1.z,"maxCV.z=",Vsplit1[1,I.z],"\t")
        }else {
          lambda1.z=lambda0.z # now z is already converged, don't update lambda.z
        }

        I.w=which.max(Vsplit1[2,])
        if (cont.w>0){
          lambda1.w=lambdaSeq.w[I.w]
          #cat("lambda.w=",lambda1.w,"maxCV.w=",Vsplit1[2,I.w],"\t")
        }else{
          lambda1.w=lambda0.w # now w is already converged, don't update lambda.w
        }

        ### stop criteria
        if ( lambda1.w==lambda0.w & niter>1 & cont.w>0 ){
          cont.w=0
          #cat(" Second level Converged on lambda = ", lambda1.w)
        }
        if ( lambda1.z==lambda0.z & niter>1 & cont.z>0 ){
          cont.z=0
          #cat(" First level Converged on lambda = ", lambda1.z)
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

        I.z=which.max(Vsplit2[1,])
        if (cont.z>0){
          alpha1.z=alphaSeq.z[I.z]
          #cat("alpha.z=",alpha1.z,"maxCV.z=",Vsplit2[1,I.z],"\t")
        }else{
          alpha1.z=alpha0.z  # now z is already converged, don't update alpha.z
        }

        I.w=which.max(Vsplit2[2,])
        if (cont.w>0){
          alpha1.w=alphaSeq.w[I.w]
          #cat("alpha.w=",alpha1.w,"maxCV.w=",Vsplit2[2,I.w],"\t")
        }else{
          alpha1.w=alpha0.w  # now w is already converged, don't update alpha.w
        }

        ### stop criteria
        if ( alpha1.w==alpha0.w & niter>1 & cont.w>0 ){
          cont.w=0
          #cat(" Second level Converged on alpha = ", alpha1.w)
        }
        if ( alpha1.z==alpha0.z & niter>1 & cont.z>0 ){
          cont.z=0
          #cat(" First level Converged on alpha = ", alpha1.z)
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
    if (!is.null(option$k)){
      alpha0.w=option$alpha.list$alpha.w[ki]
      alpha0.z=option$alpha.list$alpha.z[ki]
    }else{
      alpha0.w=option$alpha.list$alpha.w
      alpha0.z=option$alpha.list$alpha.z
    }
    cont.w=1
    cont.z=1
    if (option$parallel==T){
      # parallel computing
      no_cores=option$no_cores
      cl=makeCluster(no_cores,type="FORK")
      Vsplit1=parSapply(cl, 1:nlambda, Each)
      stopCluster(cl)
    } else{
      Vsplit1=sapply(1:nlambda, Each)
    }
    I.z=which.max(Vsplit1[1,])
    lambda1.z=lambdaSeq.z[I.z]
    #cat("lambda.z=",lambda1.z,"\t")
    I.w=which.max(Vsplit1[2,])
    lambda1.w=lambdaSeq.w[I.w]
    #cat("lambda.w=",lambda1.w,"\t")

    alpha1.z=alpha0.z
    alpha1.w=alpha0.w

  }else if (is.null(option$alpha.list)&!is.null(option$lambda.list)){
    ########################################## if only select alpha
    select="alpha"
    ki=PrevPi_d+1
    if (!is.null(option$k)){
      lambda0.w=option$lambda.list$lambda.w[ki]
      lambda0.z=option$lambda.list$lambda.z[ki]
    }else{
      lambda0.w=option$lambda.list$lambda.w
      lambda0.z=option$lambda.list$lambda.z
    }
    cont.w=1
    cont.z=1
    if (option$parallel==T){
      # parallel computing
      no_cores=option$no_cores
      cl=makeCluster(no_cores,type="FORK")
      Vsplit2=parSapply(cl, 1:nalpha, Each)
      stopCluster(cl)
    } else{
      Vsplit2=sapply(1:nalpha, Each)
    }
    I.z=which.max(Vsplit2[1,])
    alpha1.z=alphaSeq.z[I.z]
    #cat("alpha.z=",alpha1.z,"\t")
    I.w=which.max(Vsplit2[2,])
    alpha1.w=alphaSeq.w[I.w]
    #cat("alpha.w=",alpha1.w,"\t")

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

ProportionAlphaLambda=function(K,G,alphaSeq,lambdaSeq,totV,
                               Fantope_d,PrevPi,PrevPi_d, p_m, option, select){

  # If both alpha and lambda need to be selected

  if (is.null(option$alpha.list)& is.null(option$lambda.list)){

    nalpha=length(alphaSeq)
    nlambda=length(lambdaSeq)
    FVE=matrix(0, nrow = nalpha, ncol = nlambda)

    EachElement=function(x){
      index=cbind(rep(1:nalpha,rep(nlambda,nalpha)), rep(1:nlambda,nalpha))
      i=index[x,][1]
      j=index[x,][2]
      alpha=alphaSeq[i]
      lambda=lambdaSeq[j]
      projH=seqADMM_0330(K,Fantope_d,PrevPi,PrevPi_d, alpha, lambda, p_m, option$maxiter,option$eps)

      eigV=eigen(realsym(projH))$vectors[,1] # this is pxk

      FVE[i,j]=sum(diag(t(eigV)%*%G%*%eigV))/totV  # FVE (fraction of variance explained at ki-th PC)
                                                      # use raw cov instead of smoothed S as in Chen2015

      return(FVE[i,j])
    }

    if (option$parallel==T){
      # parallel computing
      no_cores=option$no_cores
      cl=makeCluster(no_cores,type="FORK")
      FVE=matrix(parSapply(cl, 1:(nalpha*nlambda), EachElement),byrow = T, ncol=nlambda, nrow=nalpha)
      stopCluster(cl)

    } else {
      FVE=matrix(sapply(1:(nalpha*nlambda), EachElement),byrow = T, ncol=nlambda, nrow=nalpha)
    }

    prop=FVE/FVE[1,1]
    #print(prop)

    #Method 1:
    #Among candidate who can explain enough variance,choose the combination of alpha+lambda that yields the smallest rFVE to localize more.
    #I=which(prop==min(prop[prop>op$rFVEproportion]),arr.ind = T)

    #Method 2:
    #Among candidate who can explain enough variance,choose the largest alpha+lambda to localize more.
    eligible=which(prop>=option$rFVEproportion,arr.ind = T)
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
    #This is for variates that are almost shrunken to 0, we want them to be 0 for simplicity.
    if (is.null(dim(eligible1))){ #if there is only one row in eligible1
      I=eligible1
    }else{
      I=eligible1[which.max(eligible1[1,]),] #choose the combo with largest alpha
    }


    alpha=alphaSeq[I[1]]
    lambda=lambdaSeq[I[2]]


  } else if (!is.null(option$alpha.list)& is.null(option$lambda.list)){

    # If only need to select lambda
    ki=PrevPi_d+1
    if (select=="w"){
      if (!is.null(option$k)){alpha=option$alpha.list$alpha.w[ki]}else{alpha=option$alpha.list$alpha.w}
    } else {
      if (!is.null(option$k)){alpha=option$alpha.list$alpha.z[ki]}else{alpha=option$alpha.list$alpha.z}
    }

    nlambda=length(lambdaSeq)
    FVE=rep(0, nlambda)

    EachLambda=function(x){
      lambda=lambdaSeq[x]
      projH=seqADMM_0330(K,Fantope_d,PrevPi,PrevPi_d, alpha, lambda, p_m, option$maxiter,option$eps)
      eigV=eigen(realsym(projH))$vectors[,1] #this is pxk
      FVE[x]=sum(diag(t(eigV)%*%G%*%eigV))/totV
      return(FVE[x])
    }

    if (option$parallel==T){
      # parallel computing
      no_cores=option$no_cores
      cl=makeCluster(no_cores,type="FORK")
      FVE=parSapply(cl, 1:nlambda, EachLambda)
      stopCluster(cl)
    } else {
      FVE=sapply(1:nlambda, EachLambda)
    }

    prop=FVE/FVE[1]
    #print(prop)
    I=max(which(prop>=option$rFVEproportion))
    lambda=lambdaSeq[I]

  } else if (is.null(option$alpha.list)& !is.null(option$lambda.list)){

    # If only need to select alpha
    ki=PrevPi_d+1
    if (select=="w"){
      if (!is.null(option$k)){lambda=option$lambda.list$lambda.w[ki]} else {lambda=option$lambda.list$lambda.w}
    } else {
      if (!is.null(option$k)){lambda=option$lambda.list$lambda.z[ki]} else {lambda=option$lambda.list$lambda.z}
    }


    nalpha=length(alphaSeq)
    FVE=rep(0, nalpha)

    EachAlpha=function(x){
      alpha=alphaSeq[x]
      projH=seqADMM_0330(K,Fantope_d,PrevPi,PrevPi_d, alpha, lambda, p_m, option$maxiter,option$eps)
      eigV=eigen(realsym(projH))$vectors[,1] #this is pxk
      FVE[x]=sum(diag(t(eigV)%*%G%*%eigV))/totV
      return(FVE[x])
    }

    if (option$parallel==T){
      # parallel computing
      no_cores=option$no_cores
      cl=makeCluster(no_cores,type="FORK")
      FVE=parSapply(cl, 1:nalpha, EachAlpha)
      stopCluster(cl)
    } else {
      FVE=sapply(1:nalpha, EachAlpha)
    }

    prop=FVE/FVE[1]
    #print(prop)
    I=max(which(prop>=option$rFVEproportion))
    alpha=alphaSeq[I]

  }


  return(list("alpha"=alpha,"lambda"=lambda))

}


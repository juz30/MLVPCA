

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

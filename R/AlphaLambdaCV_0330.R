


AlphaLambdaCV_0330=function(x_c, k, PrevPi, PrevPi_d, alpha, lambda, p_m, option){
  # for center x

  n=nrow(x_c)
  p=ncol(x_c)
  nsplit=floor(n/option$nfold)

  if (option$SmoothD=="2Diff"){
    option$SmoothD=getSmoothD(p)
  }

  nfold=option$nfold


    V=0
    for (ifold in 1:nfold){
      idtest=((ifold-1)*nsplit+1):(ifold*nsplit)
      Xtrain=x_c[-idtest,]
      Xtest=x_c[idtest,]
      Xcov1=EmpiricalS(Xtrain)
      Xcov2=EmpiricalS(Xtest)
      Strain=Xcov1-option$rho_1*option$SmoothD
      Stest=Xcov2
      projH=seqADMM_0330(Strain, k, PrevPi, PrevPi_d, alpha, lambda, p_m, option$maxiter, option$eps)
      V=V+sum(projH*Stest)
    }
    Vsplit=V/nfold

  return(Vsplit)
}

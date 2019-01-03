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


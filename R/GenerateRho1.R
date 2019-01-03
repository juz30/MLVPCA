
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


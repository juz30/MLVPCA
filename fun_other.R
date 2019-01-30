realsym=function(A){
  B=Re((A+t(A))/2)
  return(B)
}

library(Matrix)

getSmoothD_0610=function(p,m){
  Delta=matrix(0,p/m-2,p/m)
  for (i in 1:(p/m-2)){
    Delta[i,i]=0.5
    Delta[i,i+2]=0.5
    Delta[i,i+1]=-1
  }
  D0=t(Delta)%*%Delta
  
  D0_list=lapply(1:m,function(x){D0})
  D=as.matrix(do.call(bdiag,D0_list))
  return(D)
}


deflate=function(S, PrevPi){
  if (!is.null(PrevPi)){
    p=dim(S)[1]
    I=diag(p)
    S=(I-PrevPi)%*%S%*%(I-PrevPi)
    #S=S-PrevPi%*%S%*%PrevPi
  }
  return(S)
}
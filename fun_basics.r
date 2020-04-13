#' @import Matrix

#################
#    realsym    #
#################

realsym=function(A){
  B=Re((A+t(A))/2)
  return(B)
}


##################
#  get smooth D  #
##################

#library(Matrix)

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


##################
#     deflate    #
##################

deflate=function(S, PrevPi){
  if (!is.null(PrevPi)){
    p=dim(S)[1]
    I=diag(p)
    S=(I-PrevPi)%*%S%*%(I-PrevPi)
    #S=S-PrevPi%*%S%*%PrevPi
  }
  return(S)
}


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


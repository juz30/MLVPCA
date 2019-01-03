
FantopeProj=function(mat,mat0,d,ndim){
  #If we have previous projection matrix pi, we can find its orthogonal complement basis U, and multiply U to mat to get new mat.
  if (!is.null(mat0)){
    p=dim(mat0)[1]
    D=eigen(diag(p)-mat0)$values
    U=eigen(diag(p)-mat0)$vectors
    id=order(D,decreasing=T) #sort eigen value from large to small
    U=U[,id]
    U=U[,1:(p-d)]
    mat=t(U)%*%mat%*%U
  }
  #Decompose mat as in spectral decomposition, reform eigenvalue, then form mat back.
  mat=Re((mat+t(mat))/2)
  D=eigen(mat)$values
  V=eigen(mat)$vectors
  id=order(D,decreasing = T) #sort eigen value from large to small
  D=D[id]
  V=V[,id]
  theta=GetTheta(D,ndim)
  new.values=pmin(pmax(D-theta,0),1)
  newmat=V%*%diag(new.values)%*%t(V)
  #If we have previous projection matrix pi, we again multiply newmat by U.
  if (!is.null(mat0)){
    newmat=U%*%newmat%*%t(U)
  }
  return(newmat)
}




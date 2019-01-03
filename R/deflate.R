deflate=function(S, PrevPi){
    if (!is.null(PrevPi)){
      p=dim(S)[1]
      I=diag(p)
      S=(I-PrevPi)%*%S%*%(I-PrevPi)
      #S=S-PrevPi%*%S%*%PrevPi
    }
  return(S)
}
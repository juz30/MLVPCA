CombineTime=function(time){
  m=length(time)
  start=rep(0,m)
  t_x=NULL
  for (i in 1:m){
    if (i>1){
      start[i]=t_x[length(t_x)] 
    }
    t_x=c(t_x,time[[i]]+start[i])
  }
  return(t_x)
}
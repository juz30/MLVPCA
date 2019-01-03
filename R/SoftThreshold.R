SoftThreshold=function(x,lambda){
  newvalue=sign(x)*pmax(abs(x)-lambda,0) #pmax is pointwise max
  return(newvalue)
}
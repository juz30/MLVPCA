
VectorSoftThreshold=function(Soft,lambda,p_m){
 y1=NULL
 m=length(p_m)
 for (i in 1:m){
   y_i=NULL
   for (j in 1:m){
     b=block(Soft,i,j,p_m)
     if (sqrt(sum(b^2))==0){
       y_ij=matrix(0,nrow=p_m[i],ncol=p_m[j])
     }else{
     y_ij=max(1-lambda*p_m[i]/sqrt(sum(b^2)),0)*b
     }
     y_i=cbind(y_i,y_ij)
   }
   y1=rbind(y1,y_i)
 }
  return(y1)
}


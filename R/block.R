block=function(h,i,j,p_m,shape="square"){
  if (shape=="square"){
    if (i==1){a=1}
    else{a=sum(p_m[1:i-1])+1}
    if (j==1){b=1}
    else{b=sum(p_m[1:j-1])+1}
    c=sum(p_m[1:i])
    d=sum(p_m[1:j])
    return(h[a:c,b:d])  
  }
  else if(shape=="long"){ # split rows into m parts with full columns in each
    if (i==1){a=1}
    else{a=sum(p_m[1:i-1])+1}
    c=sum(p_m[1:i])
    return(h[a:c,1:ncol(h)])
  }
  else if(shape=="wide"){ # split columns into m parts with full rows in each
      if (j==1){b=1}
      else{b=sum(p_m[1:j-1])+1}
      d=sum(p_m[1:j])
      return(h[1:nrow(h),b:d])
  }
}


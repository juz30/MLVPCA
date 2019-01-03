MultilevelS=function(Y,L1,L2,option){
  if (option$model=="2WayNested"){
    
    # Model:Y_ij(t)=mu_j(t)+Z_i(t)+W_ij(t)+e_ij(t)
    I_index=rep(1:L1,each=L2)
    n=length(I_index)

    I=max(I_index)
    n_I0=table(I_index)

    k=sum(n_I0^2)
    Y_I=rowsum(Y,I_index)

    Hw=(t(Y)%*%diag(n_I0[I_index])%*%Y-t(Y_I)%*%(Y_I))*2/(k-n)
    Hz=(n*t(Y)%*%Y-colSums(Y)%*%t(colSums(Y))-(k-n)/2*Hw)*2/(n^2-k)

    Gw=Hw/2
    Gz=(Hz-Hw)/2
    corr.rho=diag(L2)
    c=1
    
    return(list("G.w"=Gw,"G.z"=Gz, "corr.rho"=corr.rho, "c"=c))
  }
  
  
  
  if (option$model=="2WayNested.spatial"){
    n <- nrow(Y) # number of subjects
    p <- ncol(Y) # number of total time points
    
    if (is.null(option$corr.rho)){
      ### 1.Get correlation constant c
      # create function cross.semivar
      cross.semivar=function(Y,j,u){
        Y_diff=Y[seq(from=j,to=n,by=L2),]-Y[seq(from=u,to=n,by=L2),]
        cross.semivar=t(Y_diff)%*%Y_diff/(2*L1)
        return(cross.semivar)
      }
      # get sum of cross.semivar for each of the 14x14 cell
      h.w.sum=matrix(0, nrow = L2, ncol = L2)
      for (j in 1:L2){
        for (u in 1:L2){
          h.w.sum[j,u]=sum(cross.semivar(Y,j,u)*(1-diag(p))) # remove the diagonal
        }
      }
      # calculate h.w.far
      if (is.null(option$corr.far.quantile)){
        option$corr.far.quantile=0.85
      }
      upper15= quantile(h.w.sum[which(h.w.sum>0)],probs=option$corr.far.quantile) #take upper 15%
      far.index=which(h.w.sum>upper15,arr.ind = T)
      h.w.far.sum=matrix(0, nrow = p,ncol = p)
      tot=rep(0,dim(far.index)[1])
      for (i in 1:dim(far.index)[1]){
        j=far.index[i,1]
        u=far.index[i,2]
        h.w.far.sum=h.w.far.sum+cross.semivar(Y,j,u)
        tot[i]=sum(cross.semivar(Y,j,u))
      }
      h.w.far=h.w.far.sum/dim(far.index)[1]
      h.w.far.sum=sum(h.w.far*(1-diag(p)))
      # correlation between different electrode effect
      corr.rho=(h.w.far.sum-h.w.sum)/h.w.far.sum
      corr.rho[far.index]=0 # far pairs are forced to be 0.
    
      c=(L2-sum(corr.rho)/L2)/(L2-1)

    }else{
      corr.rho=option$corr.rho
      c=(L2-sum(corr.rho)/L2)/(L2-1)
    }
    
    
    ### 2.Get level specific covariance matrix
    
    # Model:Y_ij(t)=mu_j(t)+Z_i(t)+W_ij(t)+e_ij(t)
    I_index=rep(1:L1,each=L2)
    n=length(I_index)
    
    I=max(I_index)
    n_I0=table(I_index)
    
    k=sum(n_I0^2)
    Y_I=rowsum(Y,I_index)
    
    Hw=(t(Y)%*%diag(n_I0[I_index])%*%Y-t(Y_I)%*%(Y_I))*2/(k-n)
    Hz=(n*t(Y)%*%Y-colSums(Y)%*%t(colSums(Y))-(k-n)/2*Hw)*2/(n^2-k)
    
    Gw=Hw/(2*c)
    Gz=Hz/2-Hw/(2*c)
    
    return(list("G.w"=Gw,"G.z"=Gz, "corr.rho"=corr.rho, "c"=c))
    
  }
  
}




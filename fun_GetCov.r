#' @importFrom stats quantile

MultilevelS=function(Y, L1, L2, option){

  model <- option$model # different model use different method
  corr.rho <- option$corr.rho # if specified then don't need to estimate rho
  corr.uncorrpct <- option$corr.uncorrpct # is used when estimating rho

  if (model=="1Way"){
    n=nrow(Y)
    xcov=t(Y)%*%Y/n
    G.z=xcov

    p=ncol(Y)
    G.w=matrix(0,ncol=p,nrow=p)
    return(list("G.w"=G.w,"G.z"=G.z, "corr.rho"=NULL, "c"=NULL))

  } else if (model=="2WayNested"){
    n <- nrow(Y) # number of subjects
    p <- ncol(Y) # number of total time points

    if (is.null(corr.rho)){
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
      upper= quantile(h.w.sum[which(h.w.sum>0)],probs=1-corr.uncorrpct) #take upper corr.uncorrpct% quantile
      far.index=which(h.w.sum>upper,arr.ind = T)
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
      c=(L2-sum(corr.rho)/L2)/(L2-1)
      h.w.sum=NULL
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

    return(list("G.w"=Gw,"G.z"=Gz, "corr.rho"=corr.rho, "c"=c, "h.w.sum"=h.w.sum))

  }

}




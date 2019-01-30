#' Perform Multilevel localized-variate functional principal component analysis (MLVPCA)
#'
#' About the methodology:
#' Multilevel localized-variate functional principal component analysis (MLVPCA) extracts interpretable basis functions that account for intra- and inter-cluster variability in a clustered multivariate process. The methodology is motivated by popular neuroscience experiments where participant brain activity is recorded using modalities such as EEG or fMRI, summarized as power within multiple time-varying frequency bands in multiple brain regions. The basis functions found by MLVFPCA can be both localized within a variate (e.g. nonzero only within a subinterval of each frequency band) and sparse among variates (e.g. zero across an entire frequency band). The sparsity is achieved by a rank-one based convex optimization with matrix L1 and block Frobenius norm based penalties. To jointly model data across repeated measures (e.g. electrodes), the functional variability is decomposed into subject level and within subject level (e.g. electrode level) variability.
#'
#' Output:
#' predx: L1XL2 by p*m matrix, is the predicted functional data for each observation within each subject; \cr
#' pred.margin: L1XL2 by p*m matrix, is the predicted marginal functional data for each observation; \cr
#' phi.w: p*m by k matrix, each column is the within subject level eigenvector for 1 to kth principal component; \cr
#' phi.z: p*m by k matrix, each column is the subject level eigenvector for 1 to kth principal component; \cr
#' PCscore.w: L1*L2 by k matrix, contains within subject level principal component score; \cr
#' PCscore.z:  L1*L2 by k matrix, contains subject level principal component score; \cr
#' FVE.w: k-element vector, is the fraction of variance explained by each PC of the within subject level; \cr
#' FVE.z: k-element vector, is the fraction of variance explained by each PC of the subject level; \cr
#' eigValue.w: k-element vector, contains within subject level eigenvalues; \cr
#' eigValue.w: k-element vector, contains subject level eigenvalues. \cr
#'
#'
#' @name MLVPCA
#' @aliases MLVPCA
#' @param data a dataset where each row is individual time series that concatenates m length-p vectors. Each column contains L2 repeated measures nested within L1 subjects. The number of columns is pXm; and the number of rows is L1*L2.
#' @param m number of variates
#' @param L1 number of subjects
#' @param L2 the number of repeated measures within each subject
#' @param option K: the number of principal components to estimate in each level, default is 3; Rho_1.list: roughness parameter, if specified as NULL cross-validation will be used to select rho. default is NULL; Alpha.list: penalty parameter that controls the degree of among variate sparsity, if specified as NULL cross-validation or FVE will be used to select alpha, default is NULL; Lambda.list: penalty parameter that controls the degree of within variate sparsity, if specified as NULL cross-validation or FVE will be used to select alpha, default is NULL; alphaproportion: the fraction of variance one would like to sacrifice in return of interpretability using FVE method to select penalty parameters. If specified as a number, FVE method will be used; if specified as NULL, cross-validation will be used, default is NULL.
#' @return See Details
#' @import parallel
#' @import pracma
#' @import magic
#' @import Matrix
#' @export

source("fun_GetCovarianceEst.r")
source("fun_SelectTuningPar.r")
source("fun_SeqADMM.r")
source("fun_other.r")

library(pracma)
library(magic)

MLVPCA=function (data, m, L1, L2, k=3, model="2WayNested.spatial",
                                 rho_1.list = NULL, rho1Seq.list = NULL, FVE_threshold=0.85, alphaproportion = NULL,
                                 alpha.list = NULL,lambda.list = NULL, alphaSeq.list = NULL, lambdaSeq.list = NULL,
                                 parallel = F, no_cores=1, maxiter=100, nfold=5, nsol=10){
  
  option = list(SmoothD = "2Diff", rho_1 = NULL, rho1Seq = NULL, nsol = nsol, FVE_M = 20, PrevPi = NULL, PrevPi_d = NULL,
                maxiter = maxiter, eps = 0.01, verbose = FALSE, k = k, alpha = NULL, lambda = NULL, alphaSeq = NULL, lambdaSeq = NULL, FVE_threshold = FVE_threshold,
                nfold = nfold, model = model, levels = c(L1 = L1, L2 = L2), rho1Seq.list = rho1Seq.list, rho_1.list = rho_1.list, alpha.list = alpha.list,lambda.list = lambda.list, 
                alphaSeq.list = alphaSeq.list, lambdaSeq.list = lambdaSeq.list, parallel = parallel, no_cores = no_cores, alphaproportion = alphaproportion, c = NULL,
                corr.rho = NULL, corr.far.quantile = 0.7)
  
  x <- as.matrix(data)
  n <- nrow(x) # number of subjects
  p <- ncol(x) # number of total time points
  p_m <- rep(p/m,m) # vector of number of time points for each variable
  t_x <- (1:p)/p
  
  option <<- option
  option$levels[1] <- L1
  option$levels[2] <- L2
  
  ### get S=xcov-D*Rho_1
  if (option$SmoothD=="2Diff"){
    option$SmoothD <- getSmoothD_0610(p,m)
  }
  # center data
  eta=matrix(0, L2, p)
  for(j in 1:L2) {
    eta[j,] <- apply(x[ (0:(L1-1)*L2) + j,  ], 2, mean)
  }
  x_c=matrix(0, nrow=n, ncol=p)
  for(j in 1:L2) {
    x_c[ (0:(L1-1)*L2) + j,  ] <- t ( t( x[ (0:(L1-1)*L2) + j,  ] ) - (eta[j,]) )
  }
  
  # get level specific covariance
  G <- MultilevelS(x_c,L1,L2,option)
  G.w <- G$G.w
  G.z <- G$G.z
  option$c <- G$c
  option$corr.rho <- G$corr.rho
  
  if (is.null(option$rho_1.list)){
    if (is.null(option$rho1Seq.list)){
      option$rho1Seq.list=list("rho1Seq.w"=GenerateRho1(G.w,m,option$nsol,NULL,NULL,NULL),
                               "rho1Seq.z"=GenerateRho1(G.z,m,option$nsol,NULL,NULL,NULL))       # get candidate rho1
    }
    option$rho_1.list <- SplitRho1_list(x_c,option)                                       # get rho1
  }
  
  cat("Rho1Seq.z",option$rho1Seq.list$rho1Seq.z,"\t",
      "Rho1Seq.w",option$rho1Seq.list$rho1Seq.w,"\t")
  cat("rho1.z=",option$rho_1.list$rho_1.z, "rho1.w=",option$rho_1.list$rho_1.w)
  
  K.w <- G.w-option$rho_1.list$rho_1.w*option$SmoothD                                     # get S=xcov-rho1*D
  K.z <- G.z-option$rho_1.list$rho_1.z*option$SmoothD
  
  
  ### get totV
  d.w <- eigen(realsym(K.w))$values
  d.w <- Re(d.w)
  totV.w <- sum(d.w[d.w>0])
  
  d.z <- eigen(realsym(K.z))$values
  d.z <- Re(d.z)
  totV.z <- sum(d.z[d.z>0])
  
  ### get alpha and lambda and enter into the algorithm
  vec.w <- matrix(0,p,p)
  alpha.w <- rep(0,p)
  lambda.w <- rep(0,p)
  FVE.w <- rep(0,p)
  PrevPi.w <- NULL
  alphaSeq.w <- list()
  lambdaSeq.w <- list()
  
  vec.z <- matrix(0,p,p)
  alpha.z <- rep(0,p)
  lambda.z <- rep(0,p)
  FVE.z <- rep(0,p)
  PrevPi.z <- NULL
  alphaSeq.z <- list()
  lambdaSeq.z <- list()
  
  ki <- 0
  cont <- 1
  
  while (cont>0 && ki<option$FVE_M){
    ki <- ki+1
    cat("PC",ki,"start")
    
    if(!is.null(option$alpha.list) & !is.null(option$lambda.list)){
      if (!is.null(option$k)){
        alpha.w[ki] <- option$alpha.list$alpha.w[ki]
        lambda.w[ki] <- option$lambda.list$lambda.w[ki]
        alpha.z[ki] <- option$alpha.list$alpha.z[ki]
        lambda.z[ki] <- option$lambda.list$lambda.z[ki]
      } else {
        alpha.w[ki] <- option$alpha.list$alpha.w
        lambda.w[ki] <- option$lambda.list$lambda.w
        alpha.z[ki] <- option$alpha.list$alpha.z
        lambda.z[ki] <- option$lambda.list$lambda.z
        
      }
      
    } else if (is.null(option$alpha.list) | is.null(option$lambda.list)){
      
      # get alpha sequence
      Gnow.w=deflate(G.w,PrevPi.w) # deflate, make the range narrower.
      alphaSeq.w[[ki]] <- GenerateRho2(Gnow.w,option$nsol,NULL,NULL,NULL)
      Gnow.z=deflate(G.z,PrevPi.z) # deflate, make the range narrower.
      alphaSeq.z[[ki]] <- GenerateRho2(Gnow.z,option$nsol,NULL,NULL,NULL)
      
      alphaSeqnow.w <- alphaSeq.w[[ki]]
      alphaSeqnow.z <- alphaSeq.z[[ki]]
      cat("alphaSeqnow.w",alphaSeqnow.w,"\t",
          "alphaSeqnow.z",alphaSeqnow.z,"\t")
      
      # get lambda sequence
      lambdaSeq.w[[ki]] <- alphaSeq.w[[ki]]
      lambdaSeq.z[[ki]] <- alphaSeq.z[[ki]]
      
      lambdaSeqnow.w <- lambdaSeq.w[[ki]]
      lambdaSeqnow.z <- lambdaSeq.z[[ki]]
      
      # get alpha and lambda
      if (!is.null(option$alphaproportion)){               # NOT ABLE TO CHOOSE ONE PARAMETER USING FVE, NEED THIS!
        
        ### Choose alpha and lambda by balancing FVE
        if (option$parallel==T){
          FVEchoice.w <- ProportionAlphaLambda_0829(K.w, G.w, 1, PrevPi.w, (ki-1), alphaSeqnow.w, lambdaSeqnow.w, p_m, option, totV.w)
          FVEchoice.z <- ProportionAlphaLambda_0829(K.z, G.z, 1, PrevPi.z, (ki-1), alphaSeqnow.z, lambdaSeqnow.z, p_m, option, totV.z)
        } else {
          FVEchoice.w <- ProportionAlphaLambda_0822(K.w, G.w, 1, PrevPi.w, (ki-1), alphaSeqnow.w, lambdaSeqnow.w, p_m, option, totV.w)
          FVEchoice.z <- ProportionAlphaLambda_0822(K.z, G.z, 1, PrevPi.z, (ki-1), alphaSeqnow.z, lambdaSeqnow.z, p_m, option, totV.z)
        }
        
        alpha.w[ki] <- FVEchoice.w[1]
        lambda.w[ki] <- FVEchoice.w[2]
        alpha.z[ki] <- FVEchoice.z[1]
        lambda.z[ki] <- FVEchoice.z[2]
        cat("FVEchoice.w=",FVEchoice.w, "and FVEchoice.z=",FVEchoice.z,"for PC",ki,"\t")# haven't consider other combinations
      } else {
        CVchoice <- CV_2WayNested_1020(x_c, 1, PrevPi.w,PrevPi.z, (ki-1), alphaSeqnow.w, lambdaSeqnow.w,alphaSeqnow.z, lambdaSeqnow.z, p_m, option)
        alpha.w[ki] <- CVchoice$alpha1.w
        lambda.w[ki] <- CVchoice$lambda1.w
        alpha.z[ki] <- CVchoice$alpha1.z
        lambda.z[ki] <- CVchoice$lambda1.z
        cat(" CVchoic.W=",alpha.w[ki],lambda.w[ki],
            " CVchoic.z=",alpha.z[ki],lambda.z[ki],
            "for PC",ki,"\t")# haven't consider other combinations
      }
      
      
    }
    
    
    # w
    projH.w <- seqADMM_0330(K.w,1,PrevPi.w, (ki-1), alpha.w[ki], lambda.w[ki], p_m, option$maxiter, option$eps, option$verbose)     # get H matrix (case1)
    # z
    projH.z <- seqADMM_0330(K.z,1,PrevPi.z, (ki-1), alpha.z[ki], lambda.z[ki], p_m, option$maxiter, option$eps, option$verbose)     # get H matrix (case1)
    
    
    # update parameters
    vec.w[,ki] <- eigen(realsym(projH.w))$vectors[,1]                       # update vec (eigenvector in the ki-th PC)
    if (vec.w[,ki][which.max(abs(vec.w[,ki]))]<0){
      vec.w[,ki] <- -vec.w[,ki]                                             # Decide the sign for eigenvector (rule: the largest abs needs to be positive)
    }
    
    vec.z[,ki] <- eigen(realsym(projH.z))$vectors[,1]                       # update vec (eigenvector in the ki-th PC)
    if (vec.z[,ki][which.max(abs(vec.z[,ki]))]<0){
      vec.z[,ki] <- -vec.z[,ki]                                             # Decide the sign for eigenvector (rule: the largest abs needs to be positive)
    }
    
    
    FVE.w[ki] <- sum(diag(t(vec.w[,ki])%*%G.w%*%vec.w[,ki]))/totV.w            # update FVE (fraction of variance explained at ki-th PC)
    FVE.z[ki] <- sum(diag(t(vec.z[,ki])%*%G.z%*%vec.z[,ki]))/totV.z            # update FVE (fraction of variance explained at ki-th PC)
    #FVE[ki]=sum(diag(t(vec[,ki])%*%S%*%vec[,ki]))/totV            # this won't work, give negative values. WHY?
    
    if (is.null(option$k)){                                          # update cont, stop criterion
      ####################################################################
      # Should the w part and z part be written toghether or seperate?   #
      # If written together, they cannot stop seperately, but if written #
      # seperately, I still want them to be estimated toghether for CV   #
      # For now I let them stop together when the total variance explaine#
      # exceed threshold.                                                #
      ####################################################################
      FVE0 <- sum(FVE.w[1:ki]*totV.w+FVE.z[1:ki]*totV.z)/(totV.w+totV.z)
      if (FVE0> option$FVE_threshold){
        cont <- 0
      }
    } else if (ki==option$k){
      cont <- 0
    }
    PrevPi.w <- vec.w[,1:ki]%*%t(vec.w[,1:ki])                                 # update PrevPi (projection matrix Pi)
    PrevPi.z <- vec.z[,1:ki]%*%t(vec.z[,1:ki])                                 # update PrevPi (projection matrix Pi)
    
  }
  
  
  ### Estimate PC scores
  # prepare parameters
  k.w <- ki  # k.w IS NOT A GOOD NAME BECAUSE WE ALREADY HAVE K.w. CHANGE NAME!
  k.z <- ki
  
  # prepare small modules
  eigVec.w <- vec.w[,1:k.w] # T X 4
  eigVec.z <- vec.z[,1:k.z]
  h <- (range(t_x)[2]-range(t_x)[1])/(p-1)
  phi.w <- (eigVec.w/sqrt(h)) # THINK ABOUT WHY, MAKE SURE IT IS NECESSARY.
  phi.z <- (eigVec.z/sqrt(h))
  eigValue.w <- diag(t(eigVec.w)%*%G.w%*%eigVec.w)*(h)
  eigValue.z <- diag(t(eigVec.z)%*%G.z%*%eigVec.z)*(h)
  r1 <- sum(eigen(K.w)$values>0)
  noise1 <- (sum(eigen(G.w)$values)-sum(eigen(K.w)$values[1:r1]))/p*option$c
  if (option$model=='2WayNested.spatial'){
    r2 <- sum(eigen(K.z)$values>0)
    noise2 <- (sum(eigen(G.z)$values)-sum(eigen(K.z)$values[1:r2]))/p/(1-1/option$c)
  }else{
    noise2 <- NA
  }
  
  
  # prepare large modules
  Gi1 <- diag(eigValue.z)
  Gi2 <- kronecker(diag(eigValue.w),option$corr.rho)
  Gi <- adiag(Gi1,Gi2)
  Zi1 <- kronecker(rep(1,L2),phi.z)
  Zi2 <- do.call('cbind',lapply(1:k.w,function(x){kronecker(diag(1,L2),phi.w[,x])}))
  Zi <- cbind(Zi1,Zi2)
  Ri <- diag(L2*p)*noise1
  Vi <- Zi%*%Gi%*%t(Zi)+Ri
  GZVi <- Gi%*%t(Zi)%*%solve(Vi)
  
  
  # Calculate PC score for each subject, and predict each subject's time series
  PCscore <- matrix(NA, nrow = n, ncol = (k.z+k.w) )
  predx <- matrix(NA,nrow = n, ncol = p)
  predx.margin <- matrix(NA,nrow = n, ncol = p)
  for (i in 1:L1){
    xi_c <- as.vector(t(x_c[1:L2+L2*(i-1),]))
    PCscorei.long <- GZVi%*%xi_c
    # reshape PCscore
    PCscorei <- cbind(kronecker(rep(1,L2),t(PCscorei.long[1:k.z])),
                      matrix(PCscorei.long[-(1:k.z)],nrow=L2))
    PCscore[1:L2+L2*(i-1),] <- PCscorei
    # predict individual time series
    predx[1:L2+L2*(i-1),] <- eta + t(cbind(phi.z,phi.w) %*% t(PCscorei))
    mu <- apply(eta,2,mean)
    predx.margin[1:L2+L2*(i-1),] <-  kronecker(rep(1,L2),t(mu)) + t(phi.z %*% t(PCscorei[,1:k.z]))
  }
  PCscore.z <- PCscore[,1:k.z]
  PCscore.w <- PCscore[,-(1:k.z)]
  
  # Estimate corr.rho from PCscore (THIS NEED TO BE CHANGED ACCORDING TO #PC ! )
  PC1score <- t(matrix(PCscore.w[,1],nrow = L2))
  PC2score <- t(matrix(PCscore.w[,2],nrow=L2))
  PC3score <- t(matrix(PCscore.w[,3],nrow=L2))
  corr.rho.PCscore <- (cor(PC1score)*eigValue.w[1]+cor(PC2score)*eigValue.w[2]+cor(PC3score)*eigValue.w[3])/(sum(eigValue.w))
  option$corr.rho.PCscore <- corr.rho.PCscore
  
  ### Prepare output
  FVE.w <- FVE.w[1:k.w]
  FVE.z <- FVE.z[1:k.z]
  totVprop.w <- totV.w/(totV.w+totV.z)
  totVprop.z <- totV.z/(totV.w+totV.z)
  option$k <- list("k.w"=k.w,"k.z"=k.z)
  option$alpha.list <- list("alpha.w"=alpha.w[1:k.w],"alpha.z"=alpha.z[1:k.z])
  option$lambda.list <- list("lambda.w"=lambda.w[1:k.w],"lambda.z"=lambda.z[1:k.z])
  option$alphaSeq.list <- list("alphaSeq.w"=alphaSeq.w,"alphaSeq.z"=alphaSeq.z)
  option$lambdaSeq.list <- list("lambdaSeq.w"=lambdaSeq.w,"lambdaSeq.z"=lambdaSeq.z)
  
  res=list('predx'=predx, 'predx.margin'=predx.margin,'PCscore.w'=PCscore.w,'PCscore.z'=PCscore.z,
           'FVE.w'=FVE.w, 'eigValue.w'=eigValue.w,'phi.w'=phi.w,'eigVec.w'=eigVec.w,'totV.w'=totV.w,"totVprop.w"=totVprop.w,
           'FVE.z'=FVE.z, 'eigValue.z'=eigValue.z,'phi.z'=phi.z,'eigVec.z'=eigVec.z,'totV.z'=totV.z,"totVprop.z"=totVprop.z,
           'option'=option, "noise1"=noise1, "noise2"=noise2)
  
  return(res)
  
}

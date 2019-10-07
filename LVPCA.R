#' Localized sparse-variate functional principal component analysis
#'
#' This function performs localized sparse-variate functional principal component analysis
#' for multilevel multivariate functional data
#'
#' @name LVPCA
#' @aliases LVPCA
#' @param data A dataset where in each row m vectorized variates are concatenated, with values evaluated at p time points per variate.
#'  Each column contains \code{L2} replicates nested within \code{L1} subjects. The dimension is \code{L1}\code{L2} by \code{m}\code{p}.
#' @param m The number of variates.
#' @param L1 The number of subjects.
#' @param L2 The number of replicate within each subject. The default is \code{1}.
#' @param k The number of principal components to select for subject level and replicate-within-subject level. The default is \code{NULL}.
#' @param FVE_threshold If \code{k} is not specified, a threshold of fraction of variance explained (FVE) to decide the number of principal components selected.
#' The default is \code{0.85}.
#' @param rFVEproportion when selecting tuning parameters alpha and lambda, how much fraction of variance explained relative to no penalty (rFVE)
#' to gaurantee while encouraging more localization. If not specified cross-validation will be used to select alpha and lambda. The default is \code{NULL}.
#' @param model \code{"1Way"} for 1-way functional model or \code{"2WayNested"} for 2-way nested functional model. The default is \code{"1Way"}.
#' @param correlation If \code{model = "2WayNested"}, whether replicates-within-subject deviations are assumed to be correlated. The default is \code{TRUE}.
#' @param corr.rho If \code{model = "2WayNested"} and \code{correlation = TURE}, the correlation matrix of replicates-within-subject deviations. If not specified
#' an unstructured sparse correlation will be assumed and estimated. The default is \code{NULL}.
#' @param corr.uncorrpct If \code{corr.rho} is not specified, the percentage of replicate-within-subject-deviation combinations that are assumed to be uncorrelated.
#' Denoted as \eqn{\delta} in the manuscript section 4.4. The default is \code{0.2}.
#' @param gamma.list The smoothing parameter gamma. For \code{model = "2WayNested"}, specify using
#' the form \code{gamma.list = list("gamma.z" = c,}\code{ "gamma.w" = c)} for subject level and
#' replicate-within-subject level; while for \code{model = "1Way"}, specify using the form \code{gamma.list = c},
#' where \code{c} is a number if \code{k} is not specified, and is a vector of k numbers if \code{k} is specified. The default is NULL.
#' @param gammaSeq.list The candidate sequence of smoothing parameter gamma to select from. For
#' \code{model = "2WayNested"}, specify \code{gammaSeq.list = list("gammaSeq.z" = c,} \cr \code{ "gammaSeq.w" = c)}
#' for subject level and replicate-within-subject level; while for \code{model = "1Way"},
#' specify using the form \code{gammaSeq.list = c}, where \code{c} is a vector if \code{k}
#' is not specified, and is a list of \code{k} vectors if \code{k} is specified. The default is NULL.
#' @param lambda.list The tuning parameter lambda. For \code{model = "2WayNested"}, specify using
#'  the form \code{lambda.list = list("lambda.z" = c, "lambda.w" = c)} for subject level and
#'  replicate-within-subject level; while for \code{model = "1Way"}, specify using the form \code{lambda.list = c},
#' where \code{c} is a number if \code{k} is not specified, and is a vector of \code{k} numbers if \code{k} is specified. The default is NULL.
#' @param alpha.list The tuning parameter alpha. For \code{model = "2WayNested"}, specify using
#' the for \code{alpha.list = list("alpha.z" = c, "alpha.w" = c)} for subject level and
#' replicate-within-subject level; while for \code{model = "1Way"}, specify using the form \code{alpha.list = c},
#' where \code{c} is a number if \code{k} is not specified, and is a vector of \code{k} numbers if \code{k} is specified. The default is NULL.
#' @param lambdaSeq.list The candidate sequence of tuning parameter lambda to select from. For \code{model =} \cr
#' \code{"2WayNested"}, specify \code{lambdaSeq.list = list("lambdaSeq.z" = c, "lambdaSeq.w" = c)}
#' for subject level and replicate-within-subject level; while for \code{model = "1Way"},
#' specify using the form \code{lambdaSeq.list = c}, where \code{c} is a vector if \code{k}
#' is not specified, and is a list of \code{k} vectors if \code{k} is specified. The default is NULL.
#' @param alphaSeq.list The candidate sequence of tuning parameter alpha to select from. For \code{model =} \cr
#' \code{"2WayNested"}, specify \code{alphaSeq.list = list("alphaSeq.z" = c, "alphaSeq.w" = c)}
#' for subject level and replicate-within-subject level; while for \code{model = "1Way"},
#' specify using the form \code{alphaSeq.list = c}, where \code{c} is a vector if \code{k}
#' is not specified, and is a list of \code{k} vectors if \code{k} is specified. The default is NULL.
#' @param SmoothD Method to smooth eigenvector estiamtes. The default is \code{"2Diff"}, where the sum of squared second difference acoss time is used as roughess matrix.
#' @param nsol The number of candidate parameters to choose from in tuning parameter selection. The default is \code{10}.
#' @param nfold The number of folds to split in cross-validation for tuning parameter selection. The default is \code{5}.
#' @param FVE_k The number of principal components to select at most. The default is \code{10}.
#' @param maxiter The maximum number of iterations allowed when performing ADMM. The default is \code{100}.
#' @param eps Convergence tolerance when performing ADMM. The default is \code{0.01}.
#' @param verbose Whether to show warning messages in ADMM. The default is \code{FALSE}.
#' @param parallel Whether to use parallel computing. The default is \code{FALSE}.
#' @param no_cores If use parallel computing, the number of cores to use. The default is \code{10}.
#' @import magic
#' @import Matrix
#' @import caTools
#' @import parallel
#' @export
#'
#' @details
#' When the argument \code{model = "2WayNested"}, a two-way functional ANOVA model is used: \eqn{Y_{ij}(t)=\mu(t)+\eta_j(t)+Z_i(t)+W_{ij}(t)+\epsilon_{ij}(t)}
#' where \eqn{\mu(t)} and \eqn{\eta_j(t)} are fixed effects of the overall mean and the replicate-specific shift from the overall mean. \eqn{Z_i(t)} and \eqn{W_{ij}(t)}
#' are the random effects of subject-level deviation and replicate-within-subject-level deviation. In some of the return values and arguments, \code{.z} and \code{.w} are
#' used as suffix to represent the subject-level and replicate-within-subject-level values such as eigenvectors, eigenvalues and principal component scores.
#'
#' @return
#' Below items are avaliable for both \code{"1Way"} and \code{"2WayNested"} model types.
#' If \code{model = "2WayNested"}, items have suffix \code{.z} or \code{.w} to indicate subject level or replicate-within-subject level.
#' \item{phi}{\code{m}\code{p} by \code{k} matrix. The estimated eigenfunctions.}
#' \item{PCscore}{\code{L1}\code{L2} by \code{k} matrix. The estimated principal component scores.}
#' \item{predx}{\code{L1}\code{L2} by \code{m}\code{p} matrix. The predicted functional data.}
#' \item{eigValue}{Length \code{k} vector. The estimated eigenvalues.}
#' \item{FVE}{Length \code{k} vector. The fraction of (level-specific) variance explained by each principal component.}
#' \item{option}{A list containing updated argument paramters such as the selected tuning paramaters, the estimated correlation matrix across the
#' replicate-within-subject-level deviations, and the estimated noise terms.
#'
#' \code{F.hat} is a \code{L2} by \code{L2} matrix containing \eqn{\hat{F}_{jk}} same as in the manuscript section 4.4. Looking for changing point
#' among the off-diagnal values of F.hat may help decide \eqn{\delta}. See Online Supplement Section 2.2 for details.}
#' Below items are only avaliable when argument \code{model = "2WayNested"}
#' \item{pred.margin}{\code{L1}\code{L2} by \code{m}\code{p} matrix. The predicted marginal functional data for each subject.}
#' \item{totVprop.z}{Numeric value. The proportion of total variance attributable to the between subject effects.}
#' \item{totVprop.w}{Numeric value. The proportion of total variance attributable to the within subject effects.}
#' \item{FVE.tot}{Numeric value. The fraction of total variance explained by principal components from both levels.}
#' @examples
#'
#' data("data.2levels")
#' # I. Analyze multilevel multivariate functional data
#' result1=LVPCA(data.2levels,m=3,L1=100,L2=5, FVE_threshold = 0.9, model="2WayNested")
#' # User define tuning parameters
#' result1.1=LVPCA(data.2levels,m=3,L1=100,L2=5, FVE_threshold = 0.9, model="2WayNested",
#' alpha.list = list("alpha.z"=0, "alpha.w"=0),lambda.list = list("lambda.z"=0, "lambda.w"=0))
#' # II.  Analyze single-level multivariate functional data
#' result2=LVPCA(data.1level,m=3,L1=100,L2=1, k = 3, model="1Way")
#' # User define tuning parameters
#' result2.2=LVPCA(data.1level,m=3,L1=100, L2=1, k = 3, model="1Way", alpha.list = c(0,0,0))
#' # III.  Analyze multilevel univariate functional data
#' result3=LVPCA(data.2levels,m=1,L1=100,L2=5, k=3, model="2WayNested")



# source("fun_basics.r")
# source("fun_GetCov.r")
# source("fun_SelectTuningPar.r")
# source("fun_SeqADMM.r")
# require(magic)
# library(magic) # used in BLUP
# require(caTools)
# library(caTools) # used in predx for 1Way model
# require(parallel)
# library(parallel)


LVPCA=function(data,m,L1,L2=1,model="1Way",k=NULL,FVE_threshold=0.85,
                                correlation=T,  rFVEproportion=NULL,  # model basic par, have to specify
                                parallel=F, no_cores=10, # computing basic, better specify
                                corr.rho=NULL, corr.uncorrpct=0.2,        # correlation par
                                gammaSeq.list=NULL, gamma.list=NULL, lambda.list=NULL, alpha.list=NULL,    # tuning par for multilevel models
                                lambdaSeq.list=NULL, alphaSeq.list=NULL,
                                SmoothD="2Diff", nsol=10, nfold=5, FVE_k=10, maxiter=100, eps=0.01, verbose=F # good to trust default

                                ){

  x <- as.matrix(data)
  n <- nrow(x) # number of observations
  p <- ncol(x) # number of total time points
  p_m <- rep(p/m,m) # vector of number of time points for each variable
  t_x <- (1:p)/p
  option <- list("m"=m,"L1"=L1,"L2"=L2,"model"=model, "k"=k, # model basic par
                 "rFVEproportion"=rFVEproportion, "FVE_threshold"=FVE_threshold, # model basic par
                 "parallel"=parallel, "no_cores"=no_cores, # computing basic
                 "corr.rho"=corr.rho, "corr.uncorrpct"=corr.uncorrpct,        # correlation par
                 "gammaSeq.list"=gammaSeq.list, "gamma.list"=gamma.list,"lambda.list"=lambda.list, "alpha.list"=alpha.list,    # tuning par
                 "lambdaSeq.list"=lambdaSeq.list, "alphaSeq.list"=alphaSeq.list,
                 "SmoothD"=SmoothD, "nsol"=nsol, "nfold"=nfold, "FVE_k"=FVE_k, "maxiter"=maxiter, "eps"=eps, # good to trust default
                 "PrevPi"=NULL,"PrevPi_d"=NULL)
  # For 1Way models, transform their tuning parameters to the 2way models format.
  if (model=="1Way"){
    if (!is.null(option$gammaSeq.list)){option$gammaSeq.list=list("gammaSeq.w"=rep(0,10),"gammaSeq.z"=option$gammaSeq.list)}
    if (!is.null(option$gamma.list)){option$gamma.list=list("gamma.w"=0, "gamma.z"=option$gamma.list)}
    if (!is.null(option$alphaSeq.list)){option$alphaSeq.list=list("alphaSeq.w"=option$alphaSeq.list,"alphaSeq.z"=option$alphaSeq.list)}
    if (!is.null(option$lambdaSeq.list)){option$lambdaSeq.list=list("lambdaSeq.w"=option$lambdaSeq.list,"lambdaSeq.z"=option$lambdaSeq.list)}
    if (!is.null(option$alpha.list)){option$alpha.list=list("alpha.w"=option$alpha.list,"alpha.z"=option$alpha.list)}
    if (!is.null(option$lambda.list)){option$lambda.list=list("lambda.w"=option$lambda.list,"lambda.z"=option$lambda.list)}
  }
  # If univariate, then only perform localized FPCA by restrict alpha=0
  if (m==1){
    if (is.null(option$k)){
      option$alpha.list=list("alpha.w"=0,"alpha.z"=0)
    }else{
      option$alpha.list=list("alpha.w"=rep(0,option$k),"alpha.z"=rep(0,option$k))
    }
  }
  # If don't consider correlation, will force correlation matrix to be diagonal matrix
  if (correlation==FALSE){
    option$corr.rho=diag(L2)
  }


  ### get S=xcov-D*Rho_1
  if (option$SmoothD=="2Diff"){
    option$SmoothD <- getSmoothD_0610(p,m)
  }
  # center data
  eta <- matrix(0, L2, p)
  for(j in 1:L2) {
    eta[j,] <- apply(x[ (0:(L1-1)*L2) + j,  ], 2, mean)
  }
  x_c <- matrix(0, nrow=n, ncol=p)
  for(j in 1:L2) {
    x_c[ (0:(L1-1)*L2) + j,  ] <- t ( t( x[ (0:(L1-1)*L2) + j,  ] ) - (eta[j,]) )
  }

  # get level specific covariance
  G <- MultilevelS(x_c, L1, L2, option)
  G.w <- G$G.w
  G.z <- G$G.z
  option$c <- G$c
  option$corr.rho <- G$corr.rho
  option$F.hat <- G$h.w.sum

  if (is.null(option$gamma.list)){
    if (is.null(option$gammaSeq.list)){
      option$gammaSeq.list <- list("gammaSeq.w"=Generategamma(G.w,m,nsol,NULL,NULL,NULL),
                               "gammaSeq.z"=Generategamma(G.z,m,nsol,NULL,NULL,NULL))     # get candidate gamma
    }
      option$gamma.list <- Splitgamma_list(x_c,option)                                       # get gamma
  }

  #cat("First level rho candidates = ",option$gammaSeq.list$gammaSeq.z,"\t","Second level rho candidates = ",option$gammaSeq.list$gammaSeq.w,"\n")
  #cat("First level selected rho = ",option$gamma.list$gamma.z, "\t","Second level selected rho = ",option$gamma.list$gamma.w, "\n")

  K.w <- G.w-option$gamma.list$gamma.w*option$SmoothD                                     # get S=xcov-gamma*D
  K.z <- G.z-option$gamma.list$gamma.z*option$SmoothD


  ### get totV
  d.w=eigen(realsym(K.w))$values
  d.w=Re(d.w)
  totV.w=sum(d.w[d.w>0])

  d.z=eigen(realsym(K.z))$values
  d.z=Re(d.z)
  totV.z=sum(d.z[d.z>0])


  ### get alpha and lambda and enter into the algorithm
  vec.w=matrix(0,p,p)
  alpha.w=rep(0,p)
  lambda.w=rep(0,p)
  FVE.w=rep(0,p)
  PrevPi.w=NULL
  alphaSeq.w=list()
  lambdaSeq.w=list()

  vec.z=matrix(0,p,p)
  alpha.z=rep(0,p)
  lambda.z=rep(0,p)
  FVE.z=rep(0,p)
  PrevPi.z=NULL
  alphaSeq.z=list()
  lambdaSeq.z=list()

  ki=0
  cont=1

  while (cont>0 && ki<option$FVE_k){
    ki = ki+1
    #cat("PC",ki,"start")

    if(!is.null(option$alpha.list) & !is.null(option$lambda.list)){
      if (!is.null(option$k)){
        alpha.w[ki]=option$alpha.list$alpha.w[ki]
        lambda.w[ki]=option$lambda.list$lambda.w[ki]
        alpha.z[ki]=option$alpha.list$alpha.z[ki]
        lambda.z[ki]=option$lambda.list$lambda.z[ki]
      } else {
        alpha.w[ki]=option$alpha.list$alpha.w
        lambda.w[ki]=option$lambda.list$lambda.w
        alpha.z[ki]=option$alpha.list$alpha.z
        lambda.z[ki]=option$lambda.list$lambda.z
      }

    } else if (is.null(option$alpha.list) | is.null(option$lambda.list)){

      # get alpha sequence

      if (!is.null(option$alphaSeq.list)){
        if (!is.null(option$k)){
          alphaSeq.w[[ki]]=option$alphaSeq.list$alphaSeq.w[[ki]]
          alphaSeq.z[[ki]]=option$alphaSeq.list$alphaSeq.z[[ki]]
        } else {
          alphaSeq.w[[ki]]=option$alphaSeq.list$alphaSeq.w
          alphaSeq.z[[ki]]=option$alphaSeq.list$alphaSeq.z
        }
      } else {
        Gnow.w=deflate(G.w,PrevPi.w) # deflate, make the range narrower.
        alphaSeq.w[[ki]]=GenerateRho2(Gnow.w,option$nsol,NULL,NULL,NULL)
        Gnow.z=deflate(G.z,PrevPi.z) # deflate, make the range narrower.
        alphaSeq.z[[ki]]=GenerateRho2(Gnow.z,option$nsol,NULL,NULL,NULL)
      }

      alphaSeqnow.w=alphaSeq.w[[ki]]
      alphaSeqnow.z=alphaSeq.z[[ki]]

      # get lambda sequence

      if (!is.null(option$lambdaSeq.list)){
        if (!is.null(option$k)){
          lambdaSeq.w[[ki]]=option$lambdaSeq.list$lambdaSeq.w[[ki]]
          lambdaSeq.z[[ki]]=option$lambdaSeq.list$lambdaSeq.z[[ki]]
        } else {
          lambdaSeq.w[[ki]]=option$lambdaSeq.list$lambdaSeq.w
          lambdaSeq.z[[ki]]=option$lambdaSeq.list$lambdaSeq.z
        }
      } else {
        lambdaSeq.w[[ki]]=alphaSeq.w[[ki]]
        lambdaSeq.z[[ki]]=alphaSeq.z[[ki]]
      }

      lambdaSeqnow.w=lambdaSeq.w[[ki]]
      lambdaSeqnow.z=lambdaSeq.z[[ki]]

      #cat("Z level alpha and lambda candidates = ", alphaSeqnow.z,"\t","W level alpha and lambda candidates = ",alphaSeqnow.w,"\n")


      # get alpha and lambda

      if (!is.null(option$rFVEproportion)){
        if (model=="1Way"){
          FVEchoice.w=list("alpha"=0,"lambda"=0)
        } else{
        FVEchoice.w=ProportionAlphaLambda(K.w, G.w, alphaSeqnow.w, lambdaSeqnow.w, totV.w,
                                          Fantope_d=1, PrevPi=PrevPi.w, PrevPi_d=(ki-1), p_m, option, select="w")
        }
        FVEchoice.z=ProportionAlphaLambda(K.z, G.z, alphaSeqnow.z, lambdaSeqnow.z, totV.z,
                                          Fantope_d=1, PrevPi=PrevPi.z, PrevPi_d=(ki-1), p_m, option, select="z")
        alpha.w[ki]=FVEchoice.w$alpha
        lambda.w[ki]=FVEchoice.w$lambda
        alpha.z[ki]=FVEchoice.z$alpha
        lambda.z[ki]=FVEchoice.z$lambda
        if (model=="1Way"){
          cat("Choice by FVE: ", "alpha=", alpha.z[ki], "lambda=", lambda.z[ki], "for PC",ki,"\n")
        }else{
          cat("Z level choice by FVE: ", "alpha=", alpha.z[ki], "lambda=", lambda.z[ki], "\t",
              "W level choice by FVE: ", "alpha=", alpha.w[ki], "lambda=", lambda.w[ki],"for PC",ki,"\n")
        }

      } else {
        CVchoice=CV_2WayNested_1020(x_c,alphaSeqnow.w, lambdaSeqnow.w, alphaSeqnow.z, lambdaSeqnow.z,
                                    Fantope_d=1, PrevPi.w=PrevPi.w, PrevPi.z=PrevPi.z, PrevPi_d=(ki-1), p_m, option)
        # Using CV w level and z level tuning pars are doing cross-validation together to save time
        alpha.w[ki]=CVchoice$alpha1.w
        lambda.w[ki]=CVchoice$lambda1.w
        alpha.z[ki]=CVchoice$alpha1.z
        lambda.z[ki]=CVchoice$lambda1.z
        if (model=="1Way"){
          cat(" Choice by CV: ", "alpha=", alpha.z[ki], "lambda=", lambda.z[ki],"for PC",ki,"\n")
        }else{
          cat("Z level choice by CV: ", "alpha=", alpha.z[ki], "lambda=", lambda.z[ki], "\t",
            "W level choice by CV: ", "alpha=", alpha.w[ki], "lambda=", lambda.w[ki],
            "for PC",ki,"\n")
        }

      }

    }


    # w
    projH.w=seqADMM_0330(K.w,1,PrevPi.w, (ki-1), alpha.w[ki], lambda.w[ki], p_m, option$maxiter, option$eps, verbose)     # get H matrix (case1)
    # z
    projH.z=seqADMM_0330(K.z,1,PrevPi.z, (ki-1), alpha.z[ki], lambda.z[ki], p_m, option$maxiter, option$eps, verbose)     # get H matrix (case1)


    # update parameters
    vec.w[,ki]=eigen(realsym(projH.w))$vectors[,1]                       # update vec (eigenvector in the ki-th PC)
    if (vec.w[,ki][which.max(abs(vec.w[,ki]))]<0){
      vec.w[,ki]=-vec.w[,ki]                                             # Decide the sign for eigenvector (rule: the largest abs needs to be positive)
    }

    vec.z[,ki]=eigen(realsym(projH.z))$vectors[,1]                       # update vec (eigenvector in the ki-th PC)
    if (vec.z[,ki][which.max(abs(vec.z[,ki]))]<0){
      vec.z[,ki]=-vec.z[,ki]                                             # Decide the sign for eigenvector (rule: the largest abs needs to be positive)
    }

    if (model=="1Way"){
      FVE.w[ki]=0
    } else {
      FVE.w[ki]=sum(diag(t(vec.w[,ki])%*%G.w%*%vec.w[,ki]))/totV.w            # update FVE (fraction of variance explained at ki-th PC)
    }
      FVE.z[ki]=sum(diag(t(vec.z[,ki])%*%G.z%*%vec.z[,ki]))/totV.z            # update FVE (fraction of variance explained at ki-th PC)
                                                                              # use raw covariance G intead of K as in Chen2015
      if (FVE.z[ki]<0){
        ki.z.pos=max(which(FVE.z[1:ki]>0))                                    # Since G.z is consistent up to Sigma-(1-1/c)*noise, it may not be positive definite
        FVE.z[ki]=0                                                           # which means eigenvalue/FVE[ki] can be<0. When that occurs, force FVE[ki]=0, and toss PCs from then.
        }else{ki.z.pos=ki}                                                    # When the loop ends, instead of keep ki for z level, keep ki.z.pos for G level.


    if (is.null(option$k)){                                                 # update cont, stop criterion
      FVE.tot=sum(FVE.w[1:ki]*totV.w+FVE.z[1:ki]*totV.z)/(totV.w+totV.z)
      if (FVE.tot > option$FVE_threshold){
        cont=0
      }
    } else if (ki==option$k){
      cont=0
    }
    PrevPi.w=vec.w[,1:ki]%*%t(vec.w[,1:ki])                                 # update PrevPi (projection matrix Pi)
    PrevPi.z=vec.z[,1:ki]%*%t(vec.z[,1:ki])                                 # update PrevPi (projection matrix Pi)

  }

  ### Estimate PC scores


  if (model=="2WayNested" ){

  # prepare parameters
  k.w=ki
  k.z=ki.z.pos

  # prepare small modules
  eigVec.w=vec.w[,1:k.w] # T X 4
  eigVec.z=vec.z[,1:k.z]
  h=(range(t_x)[2]-range(t_x)[1])/(p-1)
  phi.w=(eigVec.w/sqrt(h)) # THINK ABOUT WHY, MAKE SURE IT IS NECESSARY.
  phi.z=(eigVec.z/sqrt(h))
  eigValue.w=diag(t(eigVec.w)%*%G.w%*%eigVec.w)*(h)
  eigValue.z=diag(t(eigVec.z)%*%G.z%*%eigVec.z)*(h)

  npositive=sum(eigen(K.w)$values>0)
  noise=(sum(eigen(G.w)$values)-sum(eigen(K.w)$values[1:npositive]))/p*option$c


  # prepare large modules
  Gi1=diag(eigValue.z)
  Gi2=kronecker(diag(eigValue.w),option$corr.rho)
  Gi=adiag(Gi1,Gi2)
  Zi1=kronecker(rep(1,L2),phi.z)
  Zi2=do.call('cbind',lapply(1:k.w,function(x){kronecker(diag(1,L2),phi.w[,x])}))
  Zi=cbind(Zi1,Zi2)
  Ri=diag(L2*p)*noise
  Vi=Zi%*%Gi%*%t(Zi)+Ri
  GZVi=Gi%*%t(Zi)%*%solve(Vi)


  # Calculate PC score for each subject, and predict each subject's time series
  PCscore=matrix(NA, nrow = n, ncol = (k.z+k.w) )
  predx=matrix(NA,nrow = n, ncol = p)
  predx.margin=matrix(NA,nrow = n, ncol = p)
  for (i in 1:L1){
    xi_c=as.vector(t(x_c[1:L2+L2*(i-1),]))
    PCscorei.long=GZVi%*%xi_c
    # reshape PCscore
    PCscorei=cbind(kronecker(rep(1,L2),t(PCscorei.long[1:k.z])),
                   matrix(PCscorei.long[-(1:k.z)],nrow=L2))
    PCscore[1:L2+L2*(i-1),]=PCscorei
    # predict individual time series
    predx[1:L2+L2*(i-1),]=eta + t(cbind(phi.z,phi.w) %*% t(PCscorei))
    mu=apply(eta,2,mean)
    predx.margin[1:L2+L2*(i-1),]= kronecker(rep(1,L2),t(mu)) + t(phi.z %*% t(PCscorei[,1:k.z]))
  }
  PCscore.z=PCscore[,1:k.z]
  PCscore.w=PCscore[,-(1:k.z)]


  ### Prepare output
  FVE.w=FVE.w[1:k.w]
  FVE.z=FVE.z[1:k.z]
  totVprop.w=totV.w/(totV.w+totV.z)
  totVprop.z=totV.z/(totV.w+totV.z)
  FVE.tot=(sum(FVE.w*totV.w)+sum(FVE.z*totV.z))/(totV.w+totV.z)
  option$k=list("k.w"=k.w,"k.z"=k.z)
  option$alpha.list=list("alpha.w"=alpha.w[1:k.w],"alpha.z"=alpha.z[1:k.z])
  option$lambda.list=list("lambda.w"=lambda.w[1:k.w],"lambda.z"=lambda.z[1:k.z])
  option$alphaSeq.list=list("alphaSeq.w"=alphaSeq.w,"alphaSeq.z"=alphaSeq.z)
  option$lambdaSeq.list=list("lambdaSeq.w"=lambdaSeq.w,"lambdaSeq.z"=lambdaSeq.z)
  option$noise=noise
  option$G.w=G.w
  option$G.z=G.z
  option$K.w=K.w
  option$K.z=K.z


  res=list('predx'=predx, 'predx.margin'=predx.margin,'PCscore.w'=PCscore.w,'PCscore.z'=PCscore.z,
           'FVE.w'=FVE.w, 'eigValue.w'=eigValue.w, 'phi.w'=phi.w, "totVprop.w"=totVprop.w,
           'FVE.z'=FVE.z, 'eigValue.z'=eigValue.z, 'phi.z'=phi.z, "totVprop.z"=totVprop.z, "FVE.tot"=FVE.tot,
           'option'=option)


  } else if (model=="1Way"){

    k=ki
    eigVec=vec.z[,1:k]
    projH=eigVec%*%t(eigVec)
    FVE=FVE.z[1:k]
    totV=totV.z
    eigValue=FVE*totV
    option$k=k
    option$alpha=alpha.z[1:k]
    option$lambda=lambda.z[1:k]
    option$alphaSeq=alphaSeq.z
    option$lambdaSeq=lambdaSeq.z


    h=(range(t_x)[2]-range(t_x)[1])/(p-1)
    phi=(eigVec/sqrt(h))
    PCscore=matrix(0,n,k)
    for (i in 1:k){
      prod=x_c*matrix(rep(phi[,i],n),nrow=n,byrow = T)
      PCscore[,i]=apply(prod,1,function(x) trapz(t_x,x))
    }
    predx=matrix(rep(eta,n),nrow=n,byrow = T)+PCscore%*%t(phi)
    res=list('predx'=predx,'PCscore'=PCscore,'phi'=phi,'FVE'=FVE,'eigValue'=eigValue,'option'=option)

  }

    class(res)="LVPCA"
    return(res)


}

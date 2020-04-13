## usethis namespace: start
#' @useDynLib LVPCA, .registration = TRUE
## usethis namespace: end
NULL
## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
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
#' @param corr_rho If \code{model = "2WayNested"} and \code{correlation = TURE}, the correlation matrix of replicates-within-subject deviations. If not specified
#' an unstructured sparse correlation will be assumed and estimated. The default is \code{NULL}.
#' @param corr_uncorrpct If \code{corr_rho} is not specified, the percentage of replicate-within-subject-deviation combinations that are assumed to be uncorrelated.
#' Denoted as \eqn{\delta} in the manuscript section 4.4. The default is \code{0.2}.
#' @param gamma_list The smoothing parameter gamma. For \code{model = "2WayNested"}, specify using
#' the form \code{gamma_list = list("gamma_z" = c,}\code{ "gamma_w" = c)} for subject level and
#' replicate-within-subject level; while for \code{model = "1Way"}, specify using the form \code{gamma_list = c},
#' where \code{c} is a number if \code{k} is not specified, and is a vector of k numbers if \code{k} is specified. The default is NULL.
#' @param gammaSeq_list The candidate sequence of smoothing parameter gamma to select from. For
#' \code{model = "2WayNested"}, specify \code{gammaSeq_list = list("gammaSeq_z" = c,} \cr \code{ "gammaSeq_w" = c)}
#' for subject level and replicate-within-subject level; while for \code{model = "1Way"},
#' specify using the form \code{gammaSeq_list = c}, where \code{c} is a vector if \code{k}
#' is not specified, and is a list of \code{k} vectors if \code{k} is specified. The default is NULL.
#' @param lambda_list The tuning parameter lambda. For \code{model = "2WayNested"}, specify using
#'  the form \code{lambda_list = list("lambda_z" = c, "lambda_w" = c)} for subject level and
#'  replicate-within-subject level; while for \code{model = "1Way"}, specify using the form \code{lambda_list = c},
#' where \code{c} is a number if \code{k} is not specified, and is a vector of \code{k} numbers if \code{k} is specified. The default is NULL.
#' @param alpha_list The tuning parameter alpha. For \code{model = "2WayNested"}, specify using
#' the for \code{alpha_list = list("alpha_z" = c, "alpha_w" = c)} for subject level and
#' replicate-within-subject level; while for \code{model = "1Way"}, specify using the form \code{alpha_list = c},
#' where \code{c} is a number if \code{k} is not specified, and is a vector of \code{k} numbers if \code{k} is specified. The default is NULL.
#' @param lambdaSeq_list The candidate sequence of tuning parameter lambda to select from. For \code{model =} \cr
#' \code{"2WayNested"}, specify \code{lambdaSeq_list = list("lambdaSeq_z" = c, "lambdaSeq_w" = c)}
#' for subject level and replicate-within-subject level; while for \code{model = "1Way"},
#' specify using the form \code{lambdaSeq_list = c}, where \code{c} is a vector if \code{k}
#' is not specified, and is a list of \code{k} vectors if \code{k} is specified. The default is NULL.
#' @param alphaSeq_list The candidate sequence of tuning parameter alpha to select from. For \code{model =} \cr
#' \code{"2WayNested"}, specify \code{alphaSeq_list = list("alphaSeq_z" = c, "alphaSeq_w" = c)}
#' for subject level and replicate-within-subject level; while for \code{model = "1Way"},
#' specify using the form \code{alphaSeq_list = c}, where \code{c} is a vector if \code{k}
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
#' @import Rcpp
#'
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
#' If \code{model = "2WayNested"}, items have suffix \code{_z} or \code{_w} to indicate subject level or replicate-within-subject level.
#' \item{phi}{\code{m}\code{p} by \code{k} matrix. The estimated eigenfunctions.}
#' \item{PCscore}{\code{L1}\code{L2} by \code{k} matrix. The estimated principal component scores.}
#' \item{predx}{\code{L1}\code{L2} by \code{m}\code{p} matrix. The predicted functional data.}
#' \item{eigValue}{Length \code{k} vector. The estimated eigenvalues.}
#' \item{FVE}{Length \code{k} vector. The fraction of (level-specific) variance explained by each principal component.}
#' \item{option}{A list containing updated argument paramters such as the selected tuning paramaters, the estimated correlation matrix across the
#' replicate-within-subject-level deviations, and the estimated noise terms.
#'
#' \code{F_hat} is a \code{L2} by \code{L2} matrix containing \eqn{\hat{F}_{jk}} same as in the manuscript section 4.4. Looking for changing point
#' among the off-diagnal values of F_hat may help decide \eqn{\delta}. See Online Supplement Section 2.2 for details.}
#' Below items are only avaliable when argument \code{model = "2WayNested"}
#' \item{pred_margin}{\code{L1}\code{L2} by \code{m}\code{p} matrix. The predicted marginal functional data for each subject.}
#' \item{totVprop_z}{Numeric value. The proportion of total variance attributable to the between subject effects.}
#' \item{totVprop_w}{Numeric value. The proportion of total variance attributable to the within subject effects.}
#' \item{FVE_tot}{Numeric value. The fraction of total variance explained by principal components from both levels.}
#' @examples
#'
#' data("data.2levels")
#' # I. Analyze multilevel multivariate functional data
#' result1=LVPCA(data.2levels,m=3,L1=100,L2=5, FVE_threshold = 0.9, model="2WayNested")
#' # User define tuning parameters
#' result1.1=LVPCA(data.2levels,m=3,L1=100,L2=5, FVE_threshold = 0.9, model="2WayNested",
#' alpha_list = list("alpha_z"=0, "alpha_w"=0),lambda_list = list("lambda_z"=0, "lambda_w"=0))
#' # II.  Analyze single-level multivariate functional data
#' result2=LVPCA(data.1level,m=3,L1=100,L2=1, k = 3, model="1Way")
#' # User define tuning parameters
#' result2.2=LVPCA(data.1level,m=3,L1=100, L2=1, k = 3, model="1Way", alpha_list = c(0,0,0))
#' # III.  Analyze multilevel univariate functional data
#' result3=LVPCA(data.2levels,m=1,L1=100,L2=5, k=3, model="2WayNested")


#source("fun_basics.r")
#library(magic)
#library(Matrix)
#library(caTools)
#library(Rcpp)
#sourceCpp("fun_GetCov_c.cpp")
#sourceCpp("fun_SeqADMM_c.cpp")
#sourceCpp("fun_SelectTuningPar_c.cpp")
#sourceCpp("fun_GetCov_SeqADMM_SelectTuningPar_c.cpp")


LVPCA=function(data,m,L1,L2=1,model="1Way",k=NULL,
              rFVEproportion=NULL, FVE_threshold=0.85, # model basic par, have to specify
              correlation=T,  corr_rho=NULL, corr_uncorrpct=0.2, # correlation par
              gammaSeq_list=NULL, gamma_list=NULL,
              lambda_list=NULL, alpha_list=NULL,lambdaSeq_list=NULL, alphaSeq_list=NULL, # tuning par for multilevel models
              SmoothD="2Diff", nsol=10, nfold=5, FVE_k=10, maxiter=100, maxiter_cv=20, eps=0.01, verbose=F # good to trust default
                                ){

  x <- as.matrix(data)
  n <- nrow(x) # number of observations
  p <- ncol(x) # number of total time points
  p_m <- rep(p/m,m) # vector of number of time points for each variable
  t_x <- (1:p)/p
  option <- list("m"=m,"L1"=L1,"L2"=L2,"model"=model, "k"=k, # model basic par
                 "rFVEproportion"=rFVEproportion, "FVE_threshold"=FVE_threshold, # model basic par
                 "corr_rho"=corr_rho, "corr_uncorrpct"=corr_uncorrpct,        # correlation par
                 "gammaSeq_list"=gammaSeq_list, "gamma_list"=gamma_list,"lambda_list"=lambda_list, "alpha_list"=alpha_list,    # tuning par
                 "lambdaSeq_list"=lambdaSeq_list, "alphaSeq_list"=alphaSeq_list,
                 "SmoothD"=SmoothD, "nsol"=nsol, "nfold"=nfold, "FVE_k"=FVE_k, "maxiter"=maxiter, "maxiter_cv"=maxiter_cv, "eps"=eps, # good to trust default
                 "PrevPi"=NULL,"PrevPi_d"=NULL)
  # For 1Way models, transform their tuning parameters to the 2way models format.
  if (model=="1Way"){
    if (!is.null(option$gammaSeq_list)){option$gammaSeq_list=list("gammaSeq_w"=rep(0,10),"gammaSeq_z"=option$gammaSeq.list)}
    if (!is.null(option$gamma_list)){option$gamma_list=list("gamma_w"= 0, "gamma_z"=option$gamma_list)}
    if (!is.null(option$alphaSeq_list)){option$alphaSeq_list=list("alphaSeq_w"=option$alphaSeq_list,"alphaSeq_z"=option$alphaSeq_list)}
    if (!is.null(option$lambdaSeq_list)){option$lambdaSeq_list=list("lambdaSeq_w"=option$lambdaSeq_list,"lambdaSeq_z"=option$lambdaSeq_list)}
    if (!is.null(option$alpha_list)){option$alpha_list=list("alpha_w"=option$alpha_list,"alpha_z"=option$alpha_list)}
    if (!is.null(option$lambda_list)){option$lambda_list=list("lambda_w"=option$lambda_list,"lambda_z"=option$lambda_list)}
  }
  # If univariate, then only perform localized FPCA by restrict alpha=0
  if (m==1){
    if (is.null(option$k)){
      option$alpha_list=list("alpha_w"=0,"alpha_z"=0)
    }else{
      option$alpha_list=list("alpha_w"=rep(0,option$k),"alpha_z"=rep(0,option$k))
    }
  }
  # If don't consider correlation, will force correlation matrix to be diagonal matrix
  if (correlation==FALSE){
    option$corr_rho=diag(L2)
  }


  ##########################################
  #                LVPCA                   #
  ##########################################


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
  G <- MultilevelS_c(x_c, L1, L2, option)
  G_w <- G$G_w
  G_z <- G$G_z
  option$c <- G$c
  if (model=="2WayNested"){
    option$corr_rho <- G$corr_rho
  }
  option$F_hat <- G$h_w_sum

  if (is.null(option$gamma_list)){
    if (is.null(option$gammaSeq_list)){
      option$gammaSeq_list <- list("gammaSeq_w"=Generategamma(G_w,m,nsol,NULL,NULL,NULL),
                                   "gammaSeq_z"=Generategamma(G_z,m,nsol,NULL,NULL,NULL))     # get candidate gamma
    }
    option$gamma_list <- CV_Gamma_c(x_c,option)                                       # get gamma
  }

  #cat("First level rho candidates = ",option$gammaSeq.list$gammaSeq.z,"\t","Second level rho candidates = ",option$gammaSeq.list$gammaSeq.w,"\n")
  #cat("First level selected rho = ",option$gamma.list$gamma.z, "\t","Second level selected rho = ",option$gamma.list$gamma.w, "\n")

  K_w <- G_w-option$gamma_list$gamma_w*option$SmoothD                                     # get S=xcov-gamma*D
  K_z <- G_z-option$gamma_list$gamma_z*option$SmoothD


  ### get totV
  d_w=eigen(realsym(K_w))$values
  d_w=Re(d_w)
  totV_w=sum(d_w[d_w>0])

  d_z=eigen(realsym(K_z))$values
  d_z=Re(d_z)
  totV_z=sum(d_z[d_z>0])


  ### get alpha and lambda and enter into the algorithm
  vec_w=matrix(0,p,p)
  alpha_w=rep(0,p)
  lambda_w=rep(0,p)
  FVE_w=rep(0,p)
  PrevPi_w=NULL
  alphaSeq_w=list()
  lambdaSeq_w=list()

  vec_z=matrix(0,p,p)
  alpha_z=rep(0,p)
  lambda_z=rep(0,p)
  FVE_z=rep(0,p)
  PrevPi_z=NULL
  alphaSeq_z=list()
  lambdaSeq_z=list()

  ki=0
  cont=1

  while (cont>0 && ki<option$FVE_k){
    ki = ki+1
    #cat("PC",ki,"start")

    if(!is.null(option$alpha_list) & !is.null(option$lambda_list)){
      if (!is.null(option$k)){
        alpha_w[ki]=option$alpha_list$alpha_w[ki]
        lambda_w[ki]=option$lambda_list$lambda_w[ki]
        alpha_z[ki]=option$alpha_list$alpha_z[ki]
        lambda_z[ki]=option$lambda_list$lambda_z[ki]
      } else {
        alpha_w[ki]=option$alpha_list$alpha_w
        lambda_w[ki]=option$lambda_list$lambda_w
        alpha_z[ki]=option$alpha_list$alpha_z
        lambda_z[ki]=option$lambda_list$lambda_z
      }

    } else if (is.null(option$alpha_list) | is.null(option$lambda_list)){

      # get alpha sequence

      if (!is.null(option$alphaSeq_list)){
        if (!is.null(option$k)){
          alphaSeq_w[[ki]]=option$alphaSeq_list$alphaSeq_w[[ki]]
          alphaSeq_z[[ki]]=option$alphaSeq_list$alphaSeq_z[[ki]]
        } else {
          alphaSeq_w[[ki]]=option$alphaSeq_list$alphaSeq_w
          alphaSeq_z[[ki]]=option$alphaSeq_list$alphaSeq_z
        }
      } else {
        Gnow_w=deflate(G_w,PrevPi_w) # deflate, make the range narrower.
        alphaSeq_w[[ki]]=GenerateRho2(Gnow_w,option$nsol,NULL,NULL,NULL)
        Gnow_z=deflate(G_z,PrevPi_z) # deflate, make the range narrower.
        alphaSeq_z[[ki]]=GenerateRho2(Gnow_z,option$nsol,NULL,NULL,NULL)
      }

      alphaSeqnow_w=alphaSeq_w[[ki]]
      alphaSeqnow_z=alphaSeq_z[[ki]]

      # get lambda sequence

      if (!is.null(option$lambdaSeq_list)){
        if (!is.null(option$k)){
          lambdaSeq_w[[ki]]=option$lambdaSeq_list$lambdaSeq_w[[ki]]
          lambdaSeq_z[[ki]]=option$lambdaSeq_list$lambdaSeq_z[[ki]]
        } else {
          lambdaSeq_w[[ki]]=option$lambdaSeq_list$lambdaSeq_w
          lambdaSeq_z[[ki]]=option$lambdaSeq_list$lambdaSeq_z
        }
      } else {
        lambdaSeq_w[[ki]]=alphaSeq_w[[ki]]
        lambdaSeq_z[[ki]]=alphaSeq_z[[ki]]
      }

      lambdaSeqnow_w=lambdaSeq_w[[ki]]
      lambdaSeqnow_z=lambdaSeq_z[[ki]]

      #cat("Z level alpha and lambda candidates = ", alphaSeqnow.z,"\t","W level alpha and lambda candidates = ",alphaSeqnow.w,"\n")


      # get alpha and lambda

      if (!is.null(option$rFVEproportion)){
        if (model=="1Way"){
          FVEchoice_w=list("alpha1"=0,"lambda1"=0)
        } else{
          FVEchoice_w=FVE_AlphaLambda_c(K_w, G_w, alphaSeqnow_w, lambdaSeqnow_w, totV_w,
                                        Fantope_d=1, PrevPi_d=(ki-1), option, select="w", PrevPi=PrevPi_w)
        }
        FVEchoice_z=FVE_AlphaLambda_c(K_z, G_z, alphaSeqnow_z, lambdaSeqnow_z, totV_z,
                                      Fantope_d=1, PrevPi_d=(ki-1), option, select="z", PrevPi=PrevPi_z)
        alpha_w[ki]=FVEchoice_w$alpha1
        lambda_w[ki]=FVEchoice_w$lambda1
        alpha_z[ki]=FVEchoice_z$alpha1
        lambda_z[ki]=FVEchoice_z$lambda1
        if (model=="1Way"){
          cat("Choice by FVE: ", "alpha=", alpha_z[ki], "lambda=", lambda_z[ki], "for PC",ki,"\n")
        }else{
          cat("Z level choice by FVE: ", "alpha=", alpha_z[ki], "lambda=", lambda_z[ki], "\t",
              "W level choice by FVE: ", "alpha=", alpha_w[ki], "lambda=", lambda_w[ki],"for PC",ki,"\n")
        }

      } else {
        CVchoice=CV_AlphaLambda_c(x_c,alphaSeqnow_w, lambdaSeqnow_w, alphaSeqnow_z, lambdaSeqnow_z,
                                  Fantope_d=1, PrevPi_d=(ki-1), option, PrevPi_w=PrevPi_w, PrevPi_z=PrevPi_z )
        # Using CV w level and z level tuning pars are doing cross-validation together to save time
        alpha_w[ki]=CVchoice$alpha1_w
        lambda_w[ki]=CVchoice$lambda1_w
        alpha_z[ki]=CVchoice$alpha1_z
        lambda_z[ki]=CVchoice$lambda1_z
        if (model=="1Way"){
          cat(" Choice by CV: ", "alpha=", alpha_z[ki], "lambda=", lambda_z[ki],"for PC",ki,"\n")
        }else{
          cat("Z level choice by CV: ", "alpha=", alpha_z[ki], "lambda=", lambda_z[ki], "\t",
              "W level choice by CV: ", "alpha=", alpha_w[ki], "lambda=", lambda_w[ki],
              "for PC",ki,"\n")
        }

      }

    }

    if (model=="1Way"){
      # w
      projH_w=K_w   # get H matrix (case1)
      # z
      projH_z=seqADMM_c(K_z,1, (ki-1), alpha_z[ki], lambda_z[ki], option, PrevPi_z, verbose)     # get H matrix (case1)

    } else if (model=="2WayNested"){
      # w
      projH_w=seqADMM_c(K_w,1, (ki-1), alpha_w[ki], lambda_w[ki], option, PrevPi_w, verbose)     # get H matrix (case1)
      # z
      projH_z=seqADMM_c(K_z,1, (ki-1), alpha_z[ki], lambda_z[ki], option, PrevPi_z, verbose)     # get H matrix (case1)
    }

    # update parameters
    vec_w[,ki]=eigen(realsym(projH_w))$vectors[,1]                       # update vec (eigenvector in the ki-th PC)
    if (vec_w[,ki][which.max(abs(vec_w[,ki]))]<0){
      vec_w[,ki]=-vec_w[,ki]                                             # Decide the sign for eigenvector (rule: the largest abs needs to be positive)
    }

    vec_z[,ki]=eigen(realsym(projH_z))$vectors[,1]                       # update vec (eigenvector in the ki-th PC)
    if (vec_z[,ki][which.max(abs(vec_z[,ki]))]<0){
      vec_z[,ki]=-vec_z[,ki]                                             # Decide the sign for eigenvector (rule: the largest abs needs to be positive)
    }

    if (model=="1Way"){
      FVE_w[ki]=0
    } else {
      FVE_w[ki]=sum(diag(t(vec_w[,ki])%*%G_w%*%vec_w[,ki]))/totV_w            # update FVE (fraction of variance explained at ki-th PC)
    }
    FVE_z[ki]=sum(diag(t(vec_z[,ki])%*%G_z%*%vec_z[,ki]))/totV_z            # update FVE (fraction of variance explained at ki-th PC)
    # use raw covariance G intead of K as in Chen2015
    if (FVE_z[ki]<0){
      ki_z_pos=max(which(FVE_z[1:ki]>0))                                    # Since G.z is consistent up to Sigma-(1-1/c)*noise, it may not be positive definite
      FVE_z[ki]=0                                                           # which means eigenvalue/FVE[ki] can be<0. When that occurs, force FVE[ki]=0, and toss PCs from then.
    }else{ki_z_pos=ki}                                                    # When the loop ends, instead of keep ki for z level, keep ki.z.pos for z level.


    if (is.null(option$k)){                                                 # update cont, stop criterion
      FVE_tot=sum(FVE_w[1:ki]*totV_w+FVE_z[1:ki]*totV_z)/(totV_w+totV_z)
      if (FVE_tot > option$FVE_threshold){
        cont=0
      }
    } else if (ki==option$k){
      cont=0
    }
    PrevPi_w=vec_w[,1:ki]%*%t(vec_w[,1:ki])                                 # update PrevPi (projection matrix Pi)
    PrevPi_z=vec_z[,1:ki]%*%t(vec_z[,1:ki])                                 # update PrevPi (projection matrix Pi)

  }

  ### Estimate PC scores


  if (model=="2WayNested" ){

    # prepare parameters
    k_w=ki
    k_z=ki_z_pos

    # prepare small modules
    eigVec_w=vec_w[,1:k_w] # T X 4
    eigVec_z=vec_z[,1:k_z]
    h=(range(t_x)[2]-range(t_x)[1])/(p-1)
    phi_w=(eigVec_w/sqrt(h)) # THINK ABOUT WHY, MAKE SURE IT IS NECESSARY.
    phi_z=(eigVec_z/sqrt(h))
    eigValue_w=diag(t(eigVec_w)%*%G_w%*%eigVec_w)*(h)
    eigValue_z=diag(t(eigVec_z)%*%G_z%*%eigVec_z)*(h)

    npositive=sum(eigen(K_w)$values>0)
    noise=(sum(eigen(G_w)$values)-sum(eigen(K_w)$values[1:npositive]))/p*option$c


    # prepare large modules
    Gi1=diag(eigValue_z)
    Gi2=kronecker(diag(eigValue_w),option$corr_rho)
    Gi=adiag(Gi1,Gi2)
    Zi1=kronecker(rep(1,L2),phi_z)
    Zi2=do.call('cbind',lapply(1:k_w,function(x){kronecker(diag(1,L2),phi_w[,x])}))
    Zi=cbind(Zi1,Zi2)
    Ri=diag(L2*p)*noise
    Vi=Zi%*%Gi%*%t(Zi)+Ri
    GZVi=Gi%*%t(Zi)%*%solve(Vi)


    # Calculate PC score for each subject, and predict each subject's time series
    PCscore=matrix(NA, nrow = n, ncol = (k_z+k_w) )
    predx=matrix(NA,nrow = n, ncol = p)
    predx_margin=matrix(NA,nrow = n, ncol = p)
    for (i in 1:L1){
      xi_c=as.vector(t(x_c[1:L2+L2*(i-1),]))
      PCscorei_long=GZVi%*%xi_c
      # reshape PCscore
      PCscorei=cbind(kronecker(rep(1,L2),t(PCscorei_long[1:k_z])),
                     matrix(PCscorei_long[-(1:k_z)],nrow=L2))
      PCscore[1:L2+L2*(i-1),]=PCscorei
      # predict individual time series
      predx[1:L2+L2*(i-1),]=eta + t(cbind(phi_z,phi_w) %*% t(PCscorei))
      mu=apply(eta,2,mean)
      predx_margin[1:L2+L2*(i-1),]= kronecker(rep(1,L2),t(mu)) + t(phi_z %*% t(PCscorei[,1:k_z]))
    }
    PCscore_z=PCscore[,1:k_z]
    PCscore_w=PCscore[,-(1:k_z)]


    ### Prepare output
    FVE_w=FVE_w[1:k_w]
    FVE_z=FVE_z[1:k_z]
    totVprop_w=totV_w/(totV_w+totV_z)
    totVprop_z=totV_z/(totV_w+totV_z)
    FVE_tot=(sum(FVE_w*totV_w)+sum(FVE_z*totV_z))/(totV_w+totV_z)
    option$k=list("k_w"=k_w,"k_z"=k_z)
    option$alpha_list=list("alpha_w"=alpha_w[1:k_w],"alpha_z"=alpha_z[1:k_z])
    option$lambda_list=list("lambda_w"=lambda_w[1:k_w],"lambda_z"=lambda_z[1:k_z])
    option$alphaSeq_list=list("alphaSeq_w"=alphaSeq_w,"alphaSeq_z"=alphaSeq_z)
    option$lambdaSeq_list=list("lambdaSeq_w"=lambdaSeq_w,"lambdaSeq_z"=lambdaSeq_z)
    option$noise=noise
    option$G_w=G_w
    option$G_z=G_z
    option$K_w=K_w
    option$K_z=K_z


    res=list('predx'=predx, 'predx_margin'=predx_margin,'PCscore_w'=PCscore_w,'PCscore_z'=PCscore_z,
             'FVE_w'=FVE_w, 'eigValue_w'=eigValue_w, 'phi_w'=phi_w, "totVprop_w"=totVprop_w,
             'FVE_z'=FVE_z, 'eigValue_z'=eigValue_z, 'phi_z'=phi_z, "totVprop_z"=totVprop_z, "FVE_tot"=FVE_tot,
             'option'=option)


  } else if (model=="1Way"){

    k=ki
    eigVec=vec_z[,1:k]
    projH=eigVec%*%t(eigVec)
    FVE=FVE_z[1:k]
    totV=totV_z
    eigValue=FVE*totV
    option$k=k
    option$alpha=alpha_z[1:k]
    option$lambda=lambda_z[1:k]
    option$alphaSeq=alphaSeq_z
    option$lambdaSeq=lambdaSeq_z


    h=(range(t_x)[2]-range(t_x)[1])/(p-1)
    phi=as.matrix(eigVec/sqrt(h))
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

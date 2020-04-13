#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


//////////////////////////////////////////////////////////////////////////////////////////
//                                                                                      //
//                                   fun_GetCov_c                                       //
//                                                                                      //
//////////////////////////////////////////////////////////////////////////////////////////



//create function cross.semivar

namespace semi{

arma::mat cross_semivar(arma::mat Y, int j, int u, int L2){
  int n = Y.n_rows;
  int L1 = n/L2;
  arma::uvec rowloc_j = arma::regspace<arma::uvec>(j,L2,n-1);//seq(from=j,to=n,by=L2)
  arma::uvec rowloc_u = arma::regspace<arma::uvec>(u,L2,n-1);//seq(from=u,to=n,by=L2)
  arma::mat Y_diff = Y.rows(rowloc_j)-Y.rows(rowloc_u);
  arma::mat cross_semivar = Y_diff.t()*Y_diff/(2*L1);
  return cross_semivar;
}

}

//[[Rcpp::export]]
Rcpp::List MultilevelS_c(arma::mat Y, int L1, int L2, Rcpp::List option){

  std::string model = option("model");
  Rcpp::Nullable<Rcpp::NumericMatrix> corr_rho0 = option("corr_rho");
  double corr_uncorrpct = option("corr_uncorrpct");
  int n = Y.n_rows;
  int p = Y.n_cols;

  if (model=="1Way"){
    arma::mat xcov = Y.t()*Y/n;
    arma::mat G_z = xcov;
    arma::mat G_w = arma::zeros(p,p);
    Rcpp::Nullable<double> c = R_NilValue;
    corr_rho0 = R_NilValue;

    Rcpp::List S_list = Rcpp::List::create(Rcpp::Named("G_w")=G_w, Rcpp::Named("G_z")=G_z, Rcpp::Named("corr_rho")=corr_rho0, Rcpp::Named("c")=c);

    return S_list;

  } else if (model=="2WayNested"){

    arma::mat corr_rho;
    arma::mat h_w_sum = arma::zeros(L2,L2);

    //1.Get correlation matrix and c
    if (corr_rho0.isNull()){
      //get sum of cross.semivar for each of the L2xL2 cell
      for (int j=0; j<L2; ++j){
        for (int u=0; u<L2; ++u){
          arma::mat h_w_ju = semi::cross_semivar(Y,j,u,L2);
          h_w_ju.diag().zeros(); // #remove the diagonal
          h_w_sum(j,u) = accu(h_w_ju); //sum(h_w_ju)
        }
      }
      //calculate h_w_sum_far
      arma::vec h_w_sum_vec = vectorise(h_w_sum);
      arma::vec upperv = quantile(h_w_sum_vec.elem(find(h_w_sum_vec>0)), arma::vec {1-corr_uncorrpct}); //quantile(h.w.sum[which(h.w.sum>0)],probs=0.8)
      double upper = upperv(0);//this 2 steps take upper corr.uncorrpct% quantile. Note that the quantile can be a little different from in R because of algorithm difference
      double h_w_sum_far = mean( h_w_sum.elem(find(h_w_sum > upper)) ) ;//mean(h_w_sum[which(h_w_sum>upper)])

      //correlation between different electrode effect
      corr_rho = (h_w_sum_far-h_w_sum)/h_w_sum_far;
      corr_rho.elem(find(h_w_sum > upper)).zeros(); // far pairs are forced to be 0.
    } else {
      corr_rho = Rcpp::as<arma::mat>(corr_rho0);
    }
    double c = (L2-accu(corr_rho)/L2)/(L2-1); //c=(L2-sum(corr.rho)/L2)/(L2-1)

    //2.Get level specific covariance matrix
    //Model:Y_ij(t)=mu_j(t)+Z_i(t)+W_ij(t)+e_ij(t). Note that this part is re-written to be more efficient than in R;
    arma::mat B = kron(arma::eye(L1,L1),L2*arma::eye(L2,L2));
    arma::mat E = kron(arma::eye(L1,L1),arma::ones(L2).t());

    arma::mat Hw = (Y.t()*(B-E.t()*E)*Y)*2/(n*L2-n);
    arma::mat Hz = (Y.t()*(n*arma::eye(n,n)-arma::ones(n)*arma::ones(n).t() -B+E.t()*E)*Y)*2/(n*n-n*L2);

    arma::mat G_w = Hw/(2*c);
    arma::mat G_z = Hz/2-Hw/(2*c);

    Rcpp::List S_list = Rcpp::List::create(Rcpp::Named("G_w")=G_w, Rcpp::Named("G_z")=G_z, Rcpp::Named("corr_rho")=corr_rho, Rcpp::Named("c")=c, Rcpp::Named("h_w_sum")=h_w_sum);

    return S_list;
  }

}

/***R

# library(LVPCA)
# Y=data.1level
# option=list("model"="2WayNested","corr_rho"=NULL,"corr_uncorrpct"=0.2)
# a=MultilevelS_c(Y,L1=20,L2=5,option)

*/






//////////////////////////////////////////////////////////////////////////////////////////
//                                                                                      //
//                                   fun_SeqADMM_c                                      //
//                                                                                      //
//////////////////////////////////////////////////////////////////////////////////////////


////////////////////
//     deflate    //
////////////////////

//[[Rcpp::export]]

arma::mat deflate_c(arma::mat S,
              Rcpp::Nullable<Rcpp::NumericMatrix> PrevPi=R_NilValue){
  if (PrevPi.isNotNull()){
    arma::mat PrevPi0 = Rcpp::as<arma::mat>(PrevPi);
    int p = S.n_rows;
    S = (arma::eye(p,p)-PrevPi0)*S*(arma::eye(p,p)-PrevPi0);
  }
  return S;
}

/*** R

# mat1 <- matrix(1:9,nrow=3)
# mat0 <- c(1,1,0)%*%t(c(1,1,0))
# deflate_c(mat1,mat0)

*/


////////////////////////////////////////
//                                    //
//     First step in ADMM: get h      //
//                                    //
////////////////////////////////////////

// [[Rcpp::export]]
double GetTheta_c(arma::vec v,
                  int ndim){
  if ( (v(ndim-1)-v(ndim))>=1) {
    double theta = v(ndim-1)-1;//theta=v[ndim]-1
    return theta;
  } else {
    int p = v.n_elem; //p=length(v)
    arma::vec v1 = arma::linspace(1,p+1,p+1); //v1=1:(p+1)
    v1.subvec(0,p-1) = v;//v1[1:p]=v;
    v1.subvec(p,p) = v(p-1) - 1.0*ndim/p; //v1[p+1]=v[p]-ndim/p;
    double theta=0;
    double fnew = 0.0;
    int ddnew = 0;
    int dnew = max(arma::vec {(ndim-2)*1.0,0.0}); // dnew=max(ndim-2,0);
    double f = fnew;
    int dd = ddnew;
    int d = dnew;
    while (fnew < ndim){
      f = fnew;
      dd = ddnew;
      d = dnew;
      dnew += 1;//dnew = dnew +1
      theta = v1(dnew-1);//Don't re-define double theta, because the outside theta will not be changed;
      arma::uvec loc1 = find((v1 - theta)<1);//which((v1-theta)<1)
      ddnew = loc1(0)+1;
      fnew = (ddnew-1)*1.0+sum(v1.subvec(ddnew-1,dnew-1))-(dnew-ddnew+1)*theta;
      //Rcpp::Rcout << "dnew= " << dnew << "ddnew="<< ddnew <<"fnew="<< fnew << "theta="<<theta<<"\n";
    }
    if (fnew==ndim){
      return theta;
    }else{
      theta = v1(d-1);
      double m0 = min(arma::vec{1-(v1(dd-1)-theta),theta-v1(d)});
      while ((f+(d-dd+1)*m0) <ndim){
        f = f+(d-dd+1)*m0;
        dd += 1;
        theta -= m0;
        m0 = min(arma::vec{1-(v1(dd-1)-theta),theta-v1(d)});
        //Rcpp::Rcout << "f= " << f << "dd="<< dd <<"theta="<< theta << "m0="<< m0 <<"\n";
      }
      theta = theta - (ndim-f)/(d-dd+1) ;
      return theta;
    }
  }
}

/*** R
#GetTheta_c(c(2,1.5,0.75,0.1),1)
*/


//[[Rcpp::export]]
arma::mat FantopeProj_c(arma::mat mat1, int ndim, int d,
                  Rcpp::Nullable<Rcpp::NumericMatrix> mat0 = R_NilValue
){
  //If we have previous projection matrix pi, we can find its orthogonal complement basis U, and multiply U to arma::mat to get new mat.
  arma::mat U;
  if (mat0.isNotNull()){
    arma::mat mat00 = Rcpp::as<arma::mat>(mat0);
    int p = mat00.n_rows; //p=dim(mat0)[1]
    arma::vec D;
    eig_sym(D, U, arma::eye(p,p)-mat00); //U=eigen(diag(p)-mat0)$vectors;D=eigen(diag(p)-mat0)$values
    U = reverse(U,1); //id=order(D,decreasing=T); U=U[,id]
    U = U.cols(0,p-d-1); //U=U[,1:(p-d)]
    mat1=U.t()*mat1*U; //mat=t(U)%*%mat%*%U
  }
  //Decompose arma::mat as in spectral decomposition, reform eigenvalue, then form arma::mat back.
  mat1 = (mat1+mat1.t())/2; //mat=Re((mat+t(mat))/2)
  arma::mat V; arma::vec D;
  eig_sym(D, V, mat1);//D=eigen(mat)$values; V=eigen(mat)$vectors;
  V = reverse(V,1); D = reverse(D);//id=order(D,decreasing = T); D=D[id];V=V[,id];
  double theta = GetTheta_c(D,ndim);
  arma::vec newvalues = min(max(D-theta,arma::zeros(D.n_elem)),arma::ones(D.n_elem)); //new.values=pmin(pmax(D-theta,0),1)
  arma::mat newmat=V*diagmat(newvalues)*V.t(); //newmat=V%*%diag(new.values)%*%t(V)
  //return newmat;
  if (mat0.isNotNull()){
    newmat=U*newmat*U.t();
  }
  return newmat;
}

/*** R

# mat1 <- matrix(1:9,nrow=3)
# mat0 <- c(1,1,0)%*%t(c(1,1,0))
# FantopeProj_c(mat1=mat1,ndim=1,d=1,mat0 = mat0)

*/


////////////////////////////////////////
//                                    //
//    Second step in ADMM: get y      //
//                                    //
////////////////////////////////////////


// [[Rcpp::export]]
arma::mat SoftThreshold_c(arma::mat x, double lambda){
  int n=x.n_rows;
  arma::mat newvalue=sign(x)%max(abs(x)-lambda,arma::zeros(n,n));
  return newvalue;
}


////////////////////////////////////////
//                                    //
// Optimization using ADMM:           //
// Iteratively solve for projection   //
// matrix H.                          //
//                                    //
////////////////////////////////////////

//[[Rcpp::export]]
arma::mat seqADMM_c(arma::mat S, int ndim, int PrevPi_d, double alpha, double lambda, Rcpp::List option,
              Rcpp::Nullable<Rcpp::NumericMatrix> PrevPi=R_NilValue,
              bool verbose=false
){
  int p = S.n_cols;
  int m = option("m");
  int maxiter = option("maxiter");
  double eps = option("eps");

  if ( (alpha==0) & (lambda==0)){
    S = deflate_c(S, PrevPi);
    S = (S + S.t())/2; //realsym(S)
    arma::vec D; arma::mat V;
    eig_sym(D,V,S); //D=eigen(realsym(S))$values
    V = reverse(V,1); //id=order(D,decreasing=T)
    V = V.cols(0, ndim-1); //V=eigen(realsym(S))$vectors[,id][,1:ndim]
    arma::mat projH = V*V.t();
    return projH;
  }else {
    //starting value
    eps = ndim*eps;
    double tau = 0.1*abs(S).max(); //tau=0.1*max(abs(S)) ;search step size
    int tauStep = 2;

    arma::mat y0 = arma::zeros(p,p); //y0=matrix(0,p,p)
    arma::mat w0 = arma::zeros(p,p); //w0=matrix(0,p,p)

    int niter = 0;
    double maxnorm = eps+1;

    while((niter<maxiter)&(maxnorm>eps)){
      //update h
      arma::mat h = FantopeProj_c(y0-w0+S/tau, ndim, PrevPi_d, PrevPi);
      //update y
      arma::mat y1=arma::zeros(p,p);
      for(int rowi = 0; rowi<m; ++rowi){ //rowi in 1:m
        for(int colj = 0; colj<m; ++colj){ //colj in 1:m
          arma::uvec locRow = arma::regspace<arma::uvec>( (rowi*(p/m)), ((rowi+1)*(p/m)-1) ); // locRow=((rowi-1)*(p/m)+1):(rowi*(p/m))
          arma::uvec locCol = arma::regspace<arma::uvec>( (colj*(p/m)), ((colj+1)*(p/m)-1) ); // locCol=((colj-1)*(p/m)+1):(colj*(p/m))
          arma::mat y_ij = SoftThreshold_c(h(locRow, locCol)+ w0(locRow, locCol),lambda/tau); //h[locRow,locCol]
          double normy_ij = norm(y_ij, "fro"); //normy_ij=norm(y_ij,"F")
          if( normy_ij > alpha*(p/m)/tau)  {
            y1(locRow,locCol) = (normy_ij - alpha*(p/m)/tau)*y_ij/normy_ij;
          }
        }
      }
      //update w
      arma::mat w1 = w0+h-y1;

      //stop criterion
      double normr1 = (norm(h-y1,"fro"))*(norm(h-y1,"fro")); //normr1 = (norm(h-y1,"fro"))^2
      double norms1 = (norm(tau*(y0-y1),"fro"))*(norm(tau*(y0-y1),"fro"));
      maxnorm = max(arma::vec{normr1,norms1}); //maxnorm = max(normr1,norms1)
      niter += 1; //niter=niter+1

      //update
      y0 = y1;
      if (normr1>100*norms1){
        tau=tau*tauStep;
        w0=w1/tauStep;
      } else if (norms1>100*normr1){
        tau=tau/tauStep;
        w0=w1*tauStep;
      } else {
        w0=w1;
      }
    }

    if (verbose){
      //display warning message
      if (niter<maxiter){
        Rcpp::Rcout << "Warning messege:\n" << "seqADMM has converged after " << niter << " iterations" << "\n";
      } else {
        Rcpp::Rcout << "Warning messege:\n" << "seqADMM could not converge";
      }
    }

    arma::mat projH=y0;
    return projH;
  }

}

/*** R

# mat1 <- matrix(1:9,nrow=3); mat1<-(mat1+t(mat1))/2;
# mat0 <- c(1,1,0)%*%t(c(1,1,0))
# option <- list("m"=1,"maxiter"=100,"eps"=0.01)
# seqADMM_c(S=mat1,ndim=1,PrevPi_d = 0, alpha=0, lambda=0.1, option,
#           PrevPi = NULL)

*/






//////////////////////////////////////////////////////////////////////////////////////////
//                                                                                      //
//                             fun_SelectTuningPar_c                                    //
//                                                                                      //
//////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////
//                                          //
//   Select gamma through cross-validation  //
//                                          //
//////////////////////////////////////////////


// [[Rcpp::export]]

Rcpp::List CV_Gamma_c(arma::mat X, Rcpp::List option){

  Rcpp::List gammaSeq_list = option("gammaSeq_list");
  int ngamma = option("nsol");
  int L1 = option("L1");
  int L2 = option("L2");
  int nfold = option("nfold");
  arma::mat SmoothD  = option("SmoothD");

  // starting value
  arma::vec V_w = arma::zeros(ngamma);
  arma::vec V_z = arma::zeros(ngamma);
  int i1_w = 0;
  int i1_z = 0;
  arma::vec gamma_w_vec = gammaSeq_list("gammaSeq_w");
  arma::vec gamma_z_vec = gammaSeq_list("gammaSeq_z");
  int nsplit = floor(L1/nfold);

  // for each of the 10 candidate gamma perform CV and get CV <S,H>;
  for (int i=0; i<ngamma; ++i){

    double gamma_w = gamma_w_vec(i);
    double gamma_z = gamma_z_vec(i);
    // For each rho, do CV nfold times, and aggregate test tr(t(H)S) to get V. Then compare V for every rho
    // Split data by L1 (subject), if has more levels, split the data by the upper level to preserve data structure
    for (int ifold=0; ifold<nfold; ++ifold){
      arma::uvec idtest = arma::regspace<arma::uvec>(ifold*nsplit*L2,(ifold+1)*nsplit*L2-1);
      arma::mat Xtest = X.rows(idtest); // Xtest=X[idtest,]
      arma::mat Xtrain = X; Xtrain.shed_rows(idtest); // Xtrain=X[-idtest,]
      Rcpp::List G_train = MultilevelS_c(Xtrain, L1-nsplit, L2, option);
      Rcpp::List G_test = MultilevelS_c(Xtest, nsplit, L2, option);
      // Within subject
      arma::mat G_train_w = G_train("G_w");
      arma::mat K_train_w = G_train_w - gamma_w*SmoothD;
      arma::mat projH_w = seqADMM_c(K_train_w, 1, 0, 0, 0, option);
      arma::mat K_test_w = G_test("G_w");
      V_w(i) += accu(projH_w%K_test_w); // V_w=V_w+sum(projH_w*K_test_w)
      // Between subject
      arma::mat G_train_z = G_train("G_z");
      arma::mat K_train_z = G_train_z - gamma_z*SmoothD;
      arma::mat projH_z = seqADMM_c(K_train_z, 1, 0, 0, 0, option);
      arma::mat K_test_z = G_test("G_z");
      V_z(i) += accu(projH_z%K_test_z); //
    }
  }

  // Select bese rho based on largest CV <S,H>;
  i1_w = V_w.index_max(); // i1_w=which.max(V_w)
  i1_z = V_z.index_max(); // i1_z=which.max(V_z)

  double gamma_w_best;
  double gamma_z_best;
  if (i1_w< (ngamma-1)){
    gamma_w_best = gamma_w_vec(i1_w+1);
  }else{
    gamma_w_best = gamma_w_vec(i1_w);
  }
  if (i1_z< (ngamma-1)){
    gamma_z_best = gamma_z_vec(i1_z+1);
  }else{
    gamma_z_best = gamma_z_vec(i1_z);
  }

  Rcpp::List gamma_list = Rcpp::List::create(Rcpp::Named("gamma_w")=gamma_w_best, Rcpp::Named("gamma_z")=gamma_z_best);
  return gamma_list;

}


/*** R
# library(Matrix)
# source("fun_basics.r")
# library(LVPCA)
# X=data.1level
# option=list("L1"=20,"L2"=5,"nfold"=5,"SmoothD"=getSmoothD_0610(90,1),"nsol"=10,
#             "m"=1,"maxiter"=100,"eps"=0.01,
#             "model"="2WayNested","corr_rho"=NULL,"corr_uncorrpct"=0.2)
#
# source("fun_GetCov.r")
# source("fun_SelectTuningPar.r")
# optiona=list("model"="2WayNested","corr.rho"=NULL,"corr.uncorrpct"=0.2)
# a=MultilevelS(X,L1=20,L2=5,optiona)
# option$gammaSeq_list <- list("gammaSeq_w"=Generategamma(a$G.w,1,10,NULL,NULL,NULL),
#                              "gammaSeq_z"=Generategamma(a$G.z,1,10,NULL,NULL,NULL))
#
# gamma_c=CV_Gamma_c(X,option)

*/



//////////////////////////////////////////////////////
//                                                  //
// Select alpha and lambda through cross-validation //
//                                                  //
//////////////////////////////////////////////////////


namespace cv{

arma::vec each(arma::mat x_c, double alpha_w, double lambda_w, double alpha_z, double lambda_z,int cont_w, int cont_z,
         int Fantope_d, int PrevPi_d, Rcpp::List option,
         Rcpp::Nullable<Rcpp::NumericMatrix> PrevPi_w=R_NilValue,
         Rcpp::Nullable<Rcpp::NumericMatrix> PrevPi_z=R_NilValue){
  int L1 = option("L1");
  int L2 = option("L2");
  int nfold = option("nfold");
  int nsplit = floor(L1/nfold);
  Rcpp::List gamma_list = option("gamma_list");
  double gamma_w = gamma_list("gamma_w");
  double gamma_z = gamma_list("gamma_z");
  arma::mat SmoothD  = option("SmoothD");

  double V_w = 0;
  double V_z = 0;
  double Vsplit1_w = 0;
  double Vsplit1_z = 0;

  for (int ifold=0; ifold<nfold; ++ifold){
    arma::uvec idtest = arma::regspace<arma::uvec>(ifold*nsplit*L2,(ifold+1)*nsplit*L2-1);
    arma::mat Xtest = x_c.rows(idtest);
    arma::mat Xtrain = x_c; Xtrain.shed_rows(idtest);
    Rcpp::List G_train = MultilevelS_c(Xtrain, L1-nsplit, L2, option);
    Rcpp::List G_test = MultilevelS_c(Xtest, nsplit, L2, option);
    // within subject
    if (cont_w>0){
      arma::mat G_train_w = G_train("G_w");
      arma::mat K_train_w = G_train_w - gamma_w*SmoothD;
      arma::mat projH_w = seqADMM_c(K_train_w, Fantope_d, PrevPi_d, alpha_w, lambda_w, option,
                              PrevPi_w);
      arma::mat K_test_w = G_test("G_w");
      V_w += accu(projH_w%K_test_w);
    }
    // between subject
    if (cont_z>0){
      arma::mat G_train_z = G_train("G_z");
      arma::mat K_train_z = G_train_z - gamma_z*SmoothD;
      arma::mat projH_z = seqADMM_c(K_train_z, Fantope_d, PrevPi_d, alpha_z, lambda_z, option,
                              PrevPi_z);
      arma::mat K_test_z = G_test("G_z");
      V_z += accu(projH_z%K_test_z);
    }
  }

  Vsplit1_w = V_w/nfold;
  Vsplit1_z = V_z/nfold;
  arma::vec Vsplit = {Vsplit1_z,Vsplit1_w};
  return Vsplit;

}

}





// [[Rcpp::export]]
Rcpp::List CV_AlphaLambda_c(arma::mat x_c, arma::vec alphaSeq_w, arma::vec lambdaSeq_w, arma::vec alphaSeq_z, arma::vec lambdaSeq_z,
                      int Fantope_d, int PrevPi_d, Rcpp::List option,
                      Rcpp::Nullable<Rcpp::NumericMatrix> PrevPi_w=R_NilValue,
                      Rcpp::Nullable<Rcpp::NumericMatrix> PrevPi_z=R_NilValue){

  Rcpp::Nullable<Rcpp::List> alpha_list = option("alpha_list");
  Rcpp::Nullable<Rcpp::List> lambda_list = option("lambda_list");
  int nalpha = alphaSeq_z.n_elem; // we request sequence of z and w are the same
  int nlambda = lambdaSeq_z.n_elem;
  double alpha1_w;
  double lambda1_w;
  double alpha1_z;
  double lambda1_z;

  //////////////////////////////////////////////////
  // If both alpha and lambda need to be selected //
  //////////////////////////////////////////////////

  if ((alpha_list.isNull()) & (lambda_list.isNull())){
    // start value
    double alpha0_w = alphaSeq_w(1);
    double lambda0_w = 0;
    double alpha0_z = alphaSeq_z(1);
    double lambda0_z = 0;
    int cont_w = 1;
    int cont_z = 1;
    int maxiter_cv = option("maxiter_cv");
    int niter = 1;

    while( ((cont_w>0)|(cont_z>0))&(niter<maxiter_cv) ){


      if (niter % 2 == 1){

        /////////////////////////
        // first select lambda //
        /////////////////////////

        // cross-validation
        arma::mat Vsplit1 = arma::zeros(2,nlambda);
        for (int i=0; i<nlambda; ++i){
          Vsplit1.col(i) = cv::each(x_c, alpha0_w, lambdaSeq_w(i), alpha0_z, lambdaSeq_z(i), cont_w, cont_z,
                      Fantope_d, PrevPi_d, option, PrevPi_w, PrevPi_z);
        }
        // choose the best lambda.z
        int I_z = Vsplit1.row(0).index_max();
        if (cont_z>0){
          lambda1_z = lambdaSeq_z[I_z];
          //Rcpp::Rcout<<"lambda_z = "<<lambda1_z<<" maxCV_z = "<<Vsplit1(0,I_z)<<"\t";
          //Rcpp::Rcout<<"\n"<<Vsplit1<<"\n";
        }else{
          lambda1_z = lambda0_z; // now z is already converged, don't update lambda.z
        }
        // choose the best lambda.w
        int I_w = Vsplit1.row(1).index_max();
        if (cont_w>0){
          lambda1_w = lambdaSeq_w[I_w];
          //Rcpp::Rcout<<"lambda_w = "<<lambda1_w<<" maxCV_w = "<<Vsplit1(1,I_w)<<"\t";
        }else{
          lambda1_w = lambda0_w; // now w is already converged, don't update lambda.w
        }
        //stop criteria
        if ((lambda1_w==lambda0_w)&(niter>1)&(cont_w>0)){
          cont_w = 0;
          //Rcpp::Rcout<<"Second level converged on lambda = "<<lambda1_w<<"\t";
        }
        if ((lambda1_z==lambda0_z)&(niter>1)&(cont_z>0)){
          cont_z = 0;
          //Rcpp::Rcout<<"First level converged on lambda = "<<lambda1_z<<"\t";
        }
        //update lambda0
        lambda0_w = lambda1_w;
        lambda0_z = lambda1_z;
        niter += 1;

      } else if (niter % 2 == 0){

        /////////////////////////
        //  then select alpha  //
        /////////////////////////

        // cross-validation
        arma::mat Vsplit2 = arma::zeros(2,nalpha);
        for (int j=0; j<nalpha; ++j){
          Vsplit2.col(j) = cv::each(x_c, alphaSeq_w(j), lambda0_w, alphaSeq_z(j), lambda0_z, cont_w, cont_z,
                      Fantope_d, PrevPi_d, option, PrevPi_w, PrevPi_z);
        }
        // choose the best alpha.z
        int J_z = Vsplit2.row(0).index_max();
        if (cont_z>0){
          alpha1_z = alphaSeq_z[J_z];
          //Rcpp::Rcout<<"alpha_z = "<<alpha1_z<<" maxCV_z = "<<Vsplit2(0,J_z)<<"\t";
        }else{
          alpha1_z = alpha0_z; // now z is already converged, don't update lambda.z
        }
        // choose the best alpha.w
        int J_w = Vsplit2.row(1).index_max();
        if (cont_w>0){
          alpha1_w = alphaSeq_w[J_w];
          //Rcpp::Rcout<<"alpha_w = "<<alpha1_w<<" maxCV_w = "<<Vsplit2(1,J_w)<<"\t";
        }else{
          alpha1_w = alpha0_w; // now w is already converged, don't update lambda.w
        }
        //stop criteria
        if ((alpha1_w==alpha0_w)&(niter>1)&(cont_w>0)){
          cont_w = 0;
          //Rcpp::Rcout<<"Second level converged on alpha = "<<alpha1_w<<"\t";
        }
        if ((alpha1_z==alpha0_z)&(niter>1)&(cont_z>0)){
          cont_z = 0;
          //Rcpp::Rcout<<"First level converged on alpha = "<<alpha1_w<<"\t";
        }
        //update alpha0
        alpha0_w = alpha1_w;
        alpha0_z = alpha1_z;
        niter += 1;
      }
    }

  } else if (alpha_list.isNotNull() & lambda_list.isNull()){

    ////////////////////////////
    // If only select lambda  //
    ////////////////////////////

    Rcpp::Nullable<int> k = option("k");
    int ki = PrevPi_d+1;
    Rcpp::List alpha0_list = Rcpp::as<Rcpp::List>(alpha_list);
    double alpha0_w;
    double alpha0_z;
    if (k.isNotNull()){
      arma::vec alpha0_w_vec = alpha0_list("alpha_w");
      alpha0_w = alpha0_w_vec(ki-1);
      arma::vec alpha0_z_vec = alpha0_list("alpha_z");
      alpha0_z = alpha0_z_vec(ki-1);
    } else{
      alpha0_w = alpha0_list("alpha_w");
      alpha0_z = alpha0_list("alpha_z");
    }
    int cont_w = 1;
    int cont_z = 1;

    // cross-validation
    arma::mat Vsplit1 = arma::zeros(2,nlambda);
    for (int i=0; i<nlambda; ++i){
      Vsplit1.col(i) = cv::each(x_c, alpha0_w, lambdaSeq_w(i), alpha0_z, lambdaSeq_z(i), cont_w, cont_z,
                  Fantope_d, PrevPi_d, option, PrevPi_w, PrevPi_z);
    }
    // choose the best lambda.z
    int I_z = Vsplit1.row(0).index_max();
    lambda1_z = lambdaSeq_z[I_z];
    //Rcpp::Rcout<<"lambda_z="<<lambda1_z<<"\t";
    // choose the best lambda.w
    int I_w = Vsplit1.row(1).index_max();
    lambda1_w = lambdaSeq_w[I_w];
    //Rcpp::Rcout<<"lambda_w="<<lambda1_w<<"\t";

    alpha1_z = alpha0_z;
    alpha1_w = alpha0_w;

  } else if (alpha_list.isNull() & lambda_list.isNotNull()){

    ////////////////////////////
    // If only select alpha   //
    ////////////////////////////

    Rcpp::Nullable<int> k = option("k");
    int ki = PrevPi_d+1;
    Rcpp::List lambda0_list = Rcpp::as<Rcpp::List>(lambda_list);
    double lambda0_w;
    double lambda0_z;
    if (k.isNotNull()){
      arma::vec lambda0_w_vec = lambda0_list("lambda_w");
      lambda0_w = lambda0_w_vec(ki-1);
      arma::vec lambda0_z_vec = lambda0_list("lambda_z");
      lambda0_z = lambda0_z_vec(ki-1);
    } else{
      lambda0_w = lambda0_list("lambda_w");
      lambda0_z = lambda0_list("lambda_z");
    }
    int cont_w = 1;
    int cont_z = 1;

    // cross-validation
    arma::mat Vsplit2 = arma::zeros(2,nalpha);
    for (int j=0; j<nlambda; ++j){
      Vsplit2.col(j) = cv::each(x_c, alphaSeq_w(j), lambda0_w, alphaSeq_z(j), lambda0_z, cont_w, cont_z,
                  Fantope_d, PrevPi_d, option, PrevPi_w, PrevPi_z);
    }
    // choose the best lambda.z
    int J_z = Vsplit2.row(0).index_max();
    alpha1_z = alphaSeq_z[J_z];
    //Rcpp::Rcout<<"alpha_z="<<alpha1_z<<"\t";
    // choose the best lambda.w
    int J_w = Vsplit2.row(1).index_max();
    alpha1_w = alphaSeq_w[J_w];
    //Rcpp::Rcout<<"alpha_w="<<alpha1_w<<"\t";

    lambda1_z = lambda0_z;
    lambda1_w = lambda0_w;

  }

  Rcpp::List AlphaLambda_list = Rcpp::List::create(Rcpp::Named("alpha1_w")=alpha1_w, Rcpp::Named("lambda1_w")=lambda1_w,
                                       Rcpp::Named("alpha1_z")=alpha1_z, Rcpp::Named("lambda1_z")=lambda1_z);
  return AlphaLambda_list;
}






/*** R
# library(Matrix)
# source("fun_basics.r")
# library(LVPCA)
#   x_c=data.1level
#   option=list("L1"=20,"L2"=5,"nfold"=5,"SmoothD"=getSmoothD_0610(90,1),"nsol"=10,
#               "m"=1,"maxiter"=100,"eps"=0.01,"maxiter_cv"=20,
#               "model"="2WayNested","corr_rho"=NULL,"corr_uncorrpct"=0.2,
#               "alpha_list"=NULL,"lambda_list"=NULL,
#               "gamma_list"=list("gamma_w"=1038,"gamma_z"=170))
#
# source("fun_GetCov.r")
# source("fun_SelectTuningPar.r")
#   optiona=list("model"="2WayNested","corr.rho"=NULL,"corr.uncorrpct"=0.2)
#   a=MultilevelS(x_c,L1=20,L2=5,optiona)
#   alphaSeq_w=GenerateRho2(a$G.w,10,NULL,NULL,NULL)
#   alphaSeq_z=GenerateRho2(a$G.z,10,NULL,NULL,NULL)
#   lambdaSeq_w=alphaSeq_w
#   lambdaSeq_z=alphaSeq_z
#
#   alphalambda_c=CV_AlphaLambda_c(x_c,alphaSeq_w,lambdaSeq_w,alphaSeq_z,lambdaSeq_z,
#                                  1,0,option,NULL,NULL)

*/



//////////////////////////////////////////////////////
//                                                  //
//    Select alpha and lambda through FVE method    //
//                                                  //
//////////////////////////////////////////////////////


// [[Rcpp::export]]
Rcpp::List FVE_AlphaLambda_c(arma::mat K, arma::mat G, arma::vec alphaSeq, arma::vec lambdaSeq, double totV,
                       int Fantope_d, int PrevPi_d, Rcpp::List option, std::string select,
                       Rcpp::Nullable<Rcpp::NumericMatrix> PrevPi=R_NilValue){

  Rcpp::Nullable<Rcpp::List> alpha_list = option("alpha_list");
  Rcpp::Nullable<Rcpp::List> lambda_list = option("lambda_list");
  Rcpp::Nullable<int> k = option("k");

  double rFVEproportion = option("rFVEproportion");
  int nalpha = alphaSeq.n_elem;
  int nlambda = lambdaSeq.n_elem;
  double alpha1;
  double lambda1;

  if (alpha_list.isNull() & lambda_list.isNull()){

    /////////////////////////////////////
    // If select both alpha and lambda //
    /////////////////////////////////////

    arma::mat FVE = arma::zeros(nalpha,nlambda);
    for (int i=0; i<nalpha; ++i){
      for (int j=0; j<nlambda; ++j){
        double alpha = alphaSeq(i);
        double lambda = lambdaSeq(j);
        arma::mat projH = seqADMM_c(K, Fantope_d, PrevPi_d, alpha, lambda, option, PrevPi);
        arma::vec D; arma::mat eigV;
        eig_sym(D, eigV, (projH+projH.t())/2); //eigen(realsym(projH))$vectors[,1]
        eigV = reverse(eigV,1);
        arma::vec eigV1 = eigV.col(0);
        FVE(i,j) = accu(eigV1.t()*G*eigV1)/totV; //use raw cov instead of smoothed S as in Chen2015
      }
    }
    arma::mat prop = FVE/FVE(0,0);
    //Rcpp::Rcout<<"prop = "<<prop;
    //Among candidates who can explain enough variance, choose the largest rank of alpha+lambda to localize more
    arma::uvec eligible_vec = find(prop>=rFVEproportion);
    arma::mat eligible_mat = arma::zeros(eligible_vec.n_elem,2);
    for (int i=0; i<eligible_vec.n_elem; ++i){
      int loc = eligible_vec(i);
      eligible_mat(i,0) = loc%nalpha; //eligible_mat=which(prop>=rFVEproportion,arr.ind=T)
      eligible_mat(i,1) = ceil(loc/nalpha);
    }
    arma::mat eligible1 = eligible_mat.rows(find(sum(eligible_mat,1)==sum(eligible_mat,1).max()));
    //If multiple largest combo of rank lambda+alpha (e.g. 2+4, 4+2), choose the one with largest alpha
    arma::rowvec I;
    if (eligible1.n_rows==1){
      I = eligible1;
    } else {
      I = eligible1.row(eligible1.col(0).index_max()); //I=eligible1[which.max(eligible1[1,]),]
    }

    alpha1 = alphaSeq(I(0));
    lambda1 = lambdaSeq(I(1));

  } else if (alpha_list.isNotNull() & lambda_list.isNull()){

    /////////////////////////////////////
    //       If only select lambda     //
    /////////////////////////////////////

    Rcpp::List alpha0_list = Rcpp::as<Rcpp::List>(alpha_list);
    int ki = PrevPi_d+1;
    double alpha;
    if (select=="w"){
      if (k.isNotNull()){
        arma::vec alpha_w = alpha0_list("alpha_w");
        alpha = alpha_w(ki-1);
      } else {
        alpha = alpha0_list("alpha_w");
      }
    } else if (select=="z"){
      if (k.isNotNull()){
        arma::vec alpha_z = alpha0_list("alpha_z");
        alpha = alpha_z(ki-1);
      } else {
        alpha = alpha0_list("alpha_z");
      }
    }

    arma::vec FVE = arma::zeros(nlambda);
    for (int i=0;i<nlambda;++i){
      double lambda = lambdaSeq(i);
      arma::mat projH = seqADMM_c(K, Fantope_d, PrevPi_d, alpha, lambda, option, PrevPi);        arma::vec D; arma::mat eigV;
      eig_sym(D, eigV, (projH+projH.t())/2); //eigen(realsym(projH))$vectors[,1]
      eigV = reverse(eigV,1);
      arma::vec eigV1 = eigV.col(0);
      FVE(i) = accu(eigV1.t()*G*eigV1)/totV; //use raw cov instead of smoothed S as in Chen2015
    }
    arma::vec prop = FVE/FVE(0);
    //Rcpp::Rcout<<"prop = "<<prop;
    int I = max(find(prop>=rFVEproportion)); //I=max(which(prop>=option$rFVEproportion))

    lambda1 = lambdaSeq(I);
    alpha1 = alpha;

  } else if (alpha_list.isNull() & lambda_list.isNotNull()){

    /////////////////////////////////////
    //       If only select alpha      //
    /////////////////////////////////////

    Rcpp::List lambda0_list = Rcpp::as<Rcpp::List>(lambda_list);
    int ki = PrevPi_d+1;
    double lambda;
    if (select=="w"){
      if (k.isNotNull()){
        arma::vec lambda_w = lambda0_list("lambda_w");
        lambda = lambda_w(ki-1);
      } else {
        lambda = lambda0_list("lambda_w");
      }
    } else if (select=="z"){
      if (k.isNotNull()){
        arma::vec lambda_z = lambda0_list("lambda_z");
        lambda = lambda_z(ki-1);
      } else {
        lambda = lambda0_list("lambda_z");
      }
    }

    arma::vec FVE = arma::zeros(nalpha);
    for (int i=0;i<nalpha;++i){
      double alpha = alphaSeq(i);
      arma::mat projH = seqADMM_c(K, Fantope_d, PrevPi_d, alpha, lambda, option, PrevPi);
      arma::vec D; arma::mat eigV;
      eig_sym(D, eigV, (projH+projH.t())/2); //eigen(realsym(projH))$vectors[,1]
      eigV = reverse(eigV,1);
      arma::vec eigV1 = eigV.col(0);
      FVE(i) = accu(eigV1.t()*G*eigV1)/totV; //use raw cov instead of smoothed S as in Chen2015
    }
    arma::vec prop = FVE/FVE(0);
    //Rcpp::Rcout<<"prop = "<<prop;
    int I = max(find(prop>=rFVEproportion)); //I=max(which(prop>=option$rFVEproportion))

    alpha1 = alphaSeq(I);
    lambda1 = lambda;

  }

  Rcpp::List AlphaLambda_list = Rcpp::List::create(Rcpp::Named("alpha1")=alpha1, Rcpp::Named("lambda1")=lambda1);
  return AlphaLambda_list;


}

/*** R
# library(Matrix)
# source("fun_basics.r")
#   library(LVPCA)
#   x_c=data.1level
#   option_c=list("L1"=20,"L2"=5,"nfold"=5,"SmoothD"=getSmoothD_0610(90,1),"nsol"=10,
#                 "m"=1,"maxiter"=100,"eps"=0.01,"maxiter_cv"=20,
#                 "model"="2WayNested","corr_rho"=NULL,"corr_uncorrpct"=0.2,
#                 "alpha_list"=NULL,"lambda_list"=list("lambda_w"=c(0,0,0),"lambda_z"=c(0,0,0)),k=3,
#                 "gamma_list"=list("gamma_w"=1038,"gamma_z"=170),
#                 "rFVEproportion"=0.98)
#
#   source("fun_GetCov.r")
#   source("fun_SelectTuningPar.r")
#   optiona=list("model"="2WayNested","corr.rho"=NULL,"corr.uncorrpct"=0.2)
#   a=MultilevelS(x_c,L1=20,L2=5,optiona)
#   alphaSeq_w=GenerateRho2(a$G.w,10,NULL,NULL,NULL)
#   alphaSeq_z=GenerateRho2(a$G.z,10,NULL,NULL,NULL)
#   lambdaSeq_w=alphaSeq_w
#   lambdaSeq_z=alphaSeq_z
#
#   G.w=a$G.w
#   K.w=G.w-option$gamma.list$gamma.w*option$SmoothD
#   d.w=eigen((K.w+t(K.w))/2)$values
#   totV.w=sum(d.w[d.w>0])
#
#   G.z=a$G.z
#   K.z=G.z-option$gamma.list$gamma.z*option$SmoothD
#   d.z=eigen((K.z+t(K.z))/2)$values
#   totV.z=sum(d.z[d.z>0])
#
#   alphalambda_c=FVE_AlphaLambda_c(K.w, G.w, alphaSeq_w,lambdaSeq_w, totV.w,
#                                   1, 0, option_c, "w", NULL)
#   alphalambda_c



*/


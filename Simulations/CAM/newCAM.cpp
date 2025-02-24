/*
 This code file is based on the original code file developed by
 
 Francesco Denti, Federico Camerlenghi, Michele Guindani, and Antonietta Mira (2021).
 A common atom model for the Bayesian nonparametric analysis of nested data. Journal of the American Statistical Association.
 
 This paper (Denti et al., 2021) proposed a common atoms model to deal with one-dimensional data. We modified the code to be applied to 2-dimensional data, assuming that each 2-dimensional data vector y_{ij} is distributed as a 2-dimensional normal distribution \prod_{p=1}^{2} N(y_{pi,j} | \theta_{p}). We then compare it with our TreeTC model in the simulation study. The modified functions are "Update_Cij" and "Update_theta".
 */


#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;


///////////////////////////////////////////////////// Useful
// [[Rcpp::export]]
arma::mat PSM(arma::mat inds){
  int nsim = inds.n_rows;
  int n= inds.n_cols;
  arma::mat PSM(n,n), D(n,n); 
  D.eye(); PSM.zeros();
  
  for(int i=0; i<n; i++){
    for(int j=i+1; j<n; j++){
    arma::colvec Z = (inds.col(i)-inds.col(j));
    arma::uvec success = find(Z==0);
    PSM(i,j) = success.n_elem;
    PSM(j,i) = PSM(i,j);
    }
  }
  
  return(PSM/nsim + D);
  
}


double rt_cpp(double nu, double lambda){
  double TAA;
  TAA = R::rnorm(lambda,1.0) / exp(.5*log( R::rchisq(nu)/nu ));
  return(TAA);
}

// [[Rcpp::export]]
arma::colvec SB_given_u2(arma::colvec V) {
  
  int N = V.size();
  arma::colvec pi(N), mV2(N);
  arma::colvec mV = 1-V;
  mV2=arma::shift(cumprod(mV),+1);
  mV2(0) =1;
  return(V%mV2);
}

/////////////////////////////////////////////////////////// Update pi
// [[Rcpp::export]]
arma::colvec Update_Distributional_Sticks(arma::colvec zj, 
                        int NN_z, double alpha){
  arma::colvec v_z(NN_z);
  for(int j=0; j<(NN_z); j++){
    v_z[j] = R::rbeta(1 + accu(zj == (j+1)), alpha + accu(zj > (j+1)) );
  }
  //v_z[NN_z-1] =1.; memento: this is wrong
  return v_z; 
}

/////////////////////////////////////////////////////////// Update Omega
// [[Rcpp::export]]
arma::mat Update_Observational_Sticks(arma::colvec cij, 
                       arma::colvec zj_pg,
                       int NN_c, 
                       int NN_z, 
                       double beta){
  arma::colvec zj_pg_cpp = zj_pg -1 ,  cij_cpp = cij -1;
  
  arma::mat v_omega(NN_c, NN_z);
  
  v_omega.fill(0); v_omega.fill(0);
  
  for(int jj=0; jj<NN_z; jj++){
    for(int ii=0; ii<NN_c;  ii++){
      v_omega(ii,jj) = 
          R::rbeta( 1 + accu(zj_pg_cpp == jj && cij_cpp == ii),
                 beta + accu(zj_pg_cpp == jj && cij_cpp >  ii));
    }
   // v_omega(NN_c-1,jj) = 1.; // new line
   // omega.col(jj)      = SB_given_u2(v_omega.col(jj));
  }
  return(v_omega); 
}



// [[Rcpp::export]]
arma::mat Update_omega(arma::colvec cij, 
                       arma::colvec zj_pg,
                       int NN_c, 
                       int NN_z, 
                       double beta){
  arma::colvec zj_pg_cpp = zj_pg -1 ,  cij_cpp = cij -1;
  
  arma::mat v_omega(NN_c, NN_z), omega(NN_c, NN_z);
  
  v_omega.fill(0); v_omega.fill(0);
  
  for(int jj=0; jj<NN_z; jj++){
    for(int ii=0; ii<NN_c;  ii++){
      v_omega(ii,jj) = 
        R::rbeta( 1 + accu(zj_pg_cpp == jj && cij_cpp == ii),
                  beta + accu(zj_pg_cpp == jj && cij_cpp >  ii));
    }
    // v_omega(NN_c-1,jj) = 1.; // new line
    omega.col(jj)      = SB_given_u2(v_omega.col(jj));
  }
  return(omega); 
}


//////////////////////////////////////////////////////////// Update Zj
// [[Rcpp::export]]
arma::colvec Update_Zj_v2(arma::colvec Uj,    // uij collapsed
                                     arma::colvec xi_z, 
                                     arma::colvec xi_c,
                                     arma::colvec pi_z,
                                     arma::colvec cij,
                                     arma::mat omega,
                                     arma::colvec y_group,
                                     int NN_z, int J){
  arma::colvec newZj(J), p_zj_k(NN_z),
  possible_label = arma::linspace<arma::vec>(1, NN_z, NN_z);
  
  for(int q=0; q<J; q++){
    arma::uvec      ind = find(y_group==(q+1));
    arma::colvec    subCij = cij.elem(ind);
    arma::uvec indCij   = arma::conv_to<arma::uvec>::from(subCij-1);
    arma::mat subOmega  = omega.rows(indCij);
    
    for(int k=0; k<NN_z; k++){
      p_zj_k[k] =   log( Uj[q] < xi_z[k] ) - log(xi_z[k]) + log(pi_z[k]) +
        accu(log( subOmega.col(k))) ;
    }
    arma::colvec pp = exp(p_zj_k-max(p_zj_k));
    newZj[q] =   RcppArmadillo::sample(possible_label, 1, 1, pp)[0];      
  }
  return newZj;
}



/*
 We modified this function.
 The input data y_obser is a N*2 matrix (arma::mat)
 */
////////////////////////////////////////////////////////////////////// // Update Cij
// [[Rcpp::export]]
arma::mat Update_Cij( arma::mat y_obser,
                               arma::colvec Uij,
                               arma::colvec xi_c,
                               arma::mat omega,
                               arma::colvec zj_pg,
                               arma::mat theta,
                               int N, int NN_c){
  arma::colvec possible_label = arma::linspace<arma::vec>(1, NN_c, NN_c);
  arma::colvec IND(N);
  arma::colvec p(NN_c);
  
  for(int i=0; i<N; i++){
    for(int k=0; k<NN_c; k++){
      p[k] =
        log(xi_c[k] > Uij[i]) +
        log(omega(k,zj_pg[i]-1)) +
        R::dnorm(y_obser(i,0), theta(k,0), sqrt(theta(k,1)), 1) +  /* 1 means log=TRUE in R */
        R::dnorm(y_obser(i,1), theta(k,2), sqrt(theta(k,3)), 1) -
        log(xi_c[k]);
    }
    
    if(arma::is_finite(max(p))){
      arma::colvec pp = exp(p-max(p));
      IND(i) = RcppArmadillo::sample(possible_label, 1, 1, pp)[0];
    }else{
      IND(i) = RcppArmadillo::sample(possible_label, 1, 1)[0];
    }
    
  }
  
  return(IND);
}



/*
 We modified this function.
 The input data y_obser is a n*2 matrix (arma::mat), and theta is a NN_c*4 matrix
 */
////////////////////////////////////////////////////////////////////  Update theta
// [[Rcpp::export]]
arma::mat Update_theta(arma::mat y_obser,  /* before modified: arma::colvec y_obser */
                   arma::colvec cij,
                   double a0, double b0, double k0, double m0,
                   int NN_c, int J){
  arma::colvec cij_cpp = cij -1;
  arma::colvec y_eachDim;  /* new added */
  arma::mat theta(NN_c,4); /* before modified: theta(NN_c,2) */
  
  double ybar_i, ss_i;
  double astar, bstar, mustar, kstar;
  
  for(int i = 0; i<NN_c; i++ ){
    arma::uvec   ind = find(cij_cpp==i);
    int          n_i = ind.n_elem;
    
    for(int p = 0; p<2; p++){
        y_eachDim = y_obser.col(p).as_col();  /* Mat.col() returns arma::mat */
        
        /* arma::uvec   ind = find(cij_cpp==i); */
        arma::colvec YYY = y_eachDim.elem(ind);
        
        /* int          n_i = YYY.n_elem; */
        if(n_i > 0) {  ybar_i = mean(YYY);} else {ybar_i = 0;}
        
        if(n_i > 1) {  ss_i   = accu( pow((YYY-ybar_i),2) );} else {ss_i = 0;}
        
          astar = (a0 + n_i / 2);
          bstar = b0 + .5 * ss_i + ((k0 * n_i) * (ybar_i - m0)* (ybar_i - m0))  / ( 2 * (k0+n_i) );
          theta(i,p * 2 + 1) = 1 / rgamma(1, astar, 1/bstar)[0];
          mustar     = (k0*m0+ybar_i*n_i)/(k0+n_i);
          kstar      = k0 + n_i;
          theta(i,p * 2) = rt_cpp(2*astar, 0.0) * sqrt(bstar/(kstar*astar))+mustar;
        }
    }
    
  return(theta);
}


// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
  //
  //   http://www.rcpp.org/
  //   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
  //
  // log f(y|X,s) ~ g_1(y,x^T\beta(s)) + g_2(y,s)
  // [[Rcpp::export]]
arma::rowvec dummy(arma::vec x, arma::vec y){
  arma::mat A(3,5,fill::zeros);
  for(int i=0;i<5;i++)
    A.col(i) = x(span(1,3)) * y(i)/10;
  double a = A.max();
  arma::mat B = (exp(A.rows(span(0,1))) - exp(A.rows(span(1,2))))/(A.rows(span(0,1)) - A.rows(span(1,2)));
  return (vectorise(B.t()-a*0)).t();
}

// [[Rcpp::export]]
arma::rowvec colsums_cpp(arma::mat A){
  int k=A.n_cols;
  arma::rowvec c(k,fill::zeros);
  for(int i=0;i<k;i++)
    c(i) = sum(A.col(i));
  return(c);
}

//' hBasis function.
//' @param x A numeric vector.
//' @param knot A numeric vector.
//' @return A matrix.
//' @export
// [[Rcpp::export]]
arma::mat hBasis(arma::vec x, arma::vec knot){
  double knot_length = knot.n_rows;
  double delta_N =  (1/(knot_length - 1));
  double x_length = x.n_rows;
  arma::mat A(x_length, knot_length, fill::zeros);
  for(int i=0;i<x_length;i++){
    if(x(i)!=1){
      int a = floor(x(i) * (knot_length-1));
      A(i,a) = 1 - (x(i)-knot(a))/delta_N;
      A(i,a+1) = 1 - (knot(a+1)-x(i))/delta_N;
    }
    else
      A(i,(knot_length-1)) = 1;
  }
  return A ;
}
// [[Rcpp::export]]
arma::mat PhiBasis(arma::vec x, arma::vec knot){
  double knot_length = knot.n_rows;
  double delta_N =  (1/(knot_length - 1));
  double x_length = x.n_rows;
  arma::mat A(x_length, knot_length, fill::zeros);
  A.fill(delta_N);
  for(int i=0;i<x_length;i++){
    double a = floor(x[i] * (knot_length-1));
    if(a<(knot_length-2)){
      A(i,a) = ((delta_N + knot[a])*(x[i]-knot[a]) - 0.5 * (pow(x[i],2)-pow(knot[a],2)))/delta_N + delta_N/2;
      A(i,a+1) = delta_N/2 - ((delta_N - knot[a+1])*(knot[a+1]-x[i]) + 0.5*(pow(knot[a+1],2)-pow(x[i],2)))/delta_N;
      A(i,span(a+2,(knot_length-1))).fill(0);
    }
    else {
      A(i,a) = ((delta_N + knot[a])*(x[i]-knot[a]) - 0.5 * (pow(x[i],2)-pow(knot[a],2)))/delta_N + delta_N/2;
      A(i,a+1) = delta_N/2 - ((delta_N - knot[a+1])*(knot[a+1]-x[i]) + 0.5*(pow(knot[a+1],2)-pow(x[i],2)))/delta_N;
    }
  }
  A.col(0) = A.col(0) - delta_N/2;
  return A;
}

// arma::mat X(200,5,fill::randu); 
// arma::mat y(200,10,fill::randu);
// double n = X.n_rows;
// double K = y.n_cols;
// double p = X.n_cols; 
// arma::vec test_knots = linspace(0,1,101);
// arma::vec test_knots2 = linspace(0,1,51);
// arma::vec location_grid = linspace(0,K-1,K-1)/(K-1);
// arma::mat h_loc = hBasis(location_grid,test_knots2);
// arma::vec y_new = vectorise(y);
// arma::mat h_y = hBasis(y_new,test_knots);
// arma::mat X_new = kron(diagmat(linspace(1,1,K)),X);
// double intsize = 1/(test_knots.n_rows-1);

// [[Rcpp::export]]
List logcond_likelihood_beta(arma::mat xi, arma::mat beta, arma::mat tau, 
                             arma::mat h_loc, arma::mat X_new, arma::mat h_y, 
                             arma::vec test_knots){
  double knot_length = xi.n_rows;
  double intsize = 1/(knot_length-1);
  // double p = beta.n_rows;
  double K = beta.n_cols;
  double n = double(X_new.n_rows)/K;
  arma::mat beta_tilde = normalise(beta);
  arma::vec beta_tilde_aug = vectorise(beta_tilde,0);
  arma::vec Z = (X_new*beta_tilde_aug + 1)/2;
  //arma::vec g2 = h_loc*tau;
  arma::mat H_Z = hBasis(Z,test_knots);
  arma::mat xi_tilde = xi*H_Z.t() + repelem(tau*h_loc.t(),1,n);
  double a_max = xi_tilde.max();
  arma::mat A1 = xi_tilde.rows(span(0,knot_length-2)) - a_max;
  arma::mat A2 = xi_tilde.rows(span(1,knot_length-1)) - a_max;
  double cons = sum(log(colsums_cpp((exp(A2)-exp(A1))/((A2-A1)/intsize)))) + n*K*a_max;
  
  if(exp(cons) <=0){
    return List::create(Named("L") =-1e9 , _["H_Z.val"] = H_Z);
  }
  else{
    double temp = dot(vectorise(h_y,1),vectorise(xi_tilde,0));
        return List::create(Named("L") = temp - cons,
                        _["H_Z.val"] = H_Z);
  }
}

// [[Rcpp::export]]
List logcond_likelihood(arma::mat xi, arma::mat beta, arma::mat tau, 
                        arma::mat h_loc, arma::mat X_new,
                        arma::mat h_y, arma::vec test_knots,
                        arma::mat H_Z){
  double knot_length = xi.n_rows;
  double intsize = 1/(knot_length-1);
  // double p = beta.n_rows;
  double K = beta.n_cols;
  double n = double(X_new.n_rows)/K;
  // arma::mat beta_tilde = normalise(beta);
  // arma::vec beta_tilde_aug = vectorise(beta_tilde,0);
  // arma::vec Z = (X_new*beta_tilde_aug + 1)/2;
  // arma::mat H_Z = PhiBasis(Z,test_knots);
  // arma::mat xi_tilde = xi*H_Z.t();
  // arma::vec top(n*K,fill::zeros);
  // for(int i=0;i<n*K;i++)
    // top(i) = h_y.row(i)*xi_tilde.col(i);
   //arma::vec g2 = h_loc*tau;
   arma::mat xi_tilde = xi*H_Z.t() + repelem(tau*h_loc.t(),1,n);
   double a_max = xi_tilde.max();
   arma::mat A1 = xi_tilde.rows(span(0,knot_length-2)) - a_max;
   arma::mat A2 = xi_tilde.rows(span(1,knot_length-1)) - a_max;
   double cons = sum(log(colsums_cpp((exp(A2)-exp(A1))/((A2-A1)/intsize)))) + n*K*a_max;
   
    if(exp(cons) <=0){
      return List::create(Named("L") =-1e9 , _["H_Z.val"] = H_Z);
    }
    else{
      double temp = dot(vectorise(h_y,1),vectorise(xi_tilde,0));
      return List::create(Named("L") = temp - cons,
                          _["H_Z.val"] = H_Z);
}
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
  
  # /*** R
  # timesTwo(42)
  # */
  

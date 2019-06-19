// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

//' Ordinary Least Squares
//' 
//' Fits the standard OLS model.
//' 
//' @param y Numeric vector.
//' @param X Numeric matrix.
//' 
//' @return List containing the following:
//' \item{Beta}{Regression coefficient.}
//' \item{V}{Outcome variance.}
//' \item{Ibb}{Information matrix for beta.}
//' \item{Resid}{Outcome residuals.}
//' 
// [[Rcpp::export]]

SEXP fitOLS(const Eigen::Map<Eigen::VectorXd> y, const Eigen::Map<Eigen::MatrixXd> X){
  // Observations
  const int n = y.size();
  // Estimated parameters
  const int p = X.cols();
  // Information
  const Eigen::MatrixXd A = X.transpose()*X;
  // Estimate beta
  const Eigen::VectorXd b = A.ldlt().solve(X.transpose()*y);
  // Calculate residuals
  const Eigen::VectorXd eps = (y-X*b);
  // Scale
  const double qf = (eps.transpose()*eps);
  const double v = qf/(n-p);
  // Information
  const Eigen::MatrixXd Ibb = A/v;
  return Rcpp::List::create(Rcpp::Named("Beta")=b,Rcpp::Named("V")=v,Rcpp::Named("Ibb")=Ibb,Rcpp::Named("Resid")=eps);
}
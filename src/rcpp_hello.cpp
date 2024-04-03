#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
int rtheta(double lambda1, double lambda2, NumericVector X){
  int n = 112;
  int theta = 0;
  NumericVector prob(n - 1);
  NumericVector Xfac = factorial(X);
  for (int k = 1; k < n; k++) {
    double p_theta = 1.0;
    for (int t = 0; t < k; t++) {
      p_theta *= pow(lambda1, X(t)) * exp(-lambda1) / Xfac(t);
    }
    for (int t = k; t < n; t++) {
      p_theta *= pow(lambda2, X(t)) * exp(-lambda2) / Xfac(t);
    }
    prob(k - 1) = p_theta;
  }
  prob = prob / sum(prob);
  IntegerVector indices = seq(1, n - 1);
  theta = sample(indices, 1, false, prob)(0);
  return(theta);
}

// [[Rcpp::export]]
NumericMatrix gibbs(NumericVector X, int burnin, int iter) {
  double alpha = 1;
  int theta = 56;
  double lambda1 = 0.0;
  double lambda2 = 0.0;
  int n = 112;
  NumericMatrix res(iter, 4);

  for (int i = 0; i < burnin; i++) {
    double x1 = sum(X[Rcpp::Range(0, theta - 1)]);
    double x2 = sum(X[Rcpp::Range(theta, n - 1)]);
    lambda1 = R::rgamma(3 + x1, 1 / (theta + alpha));
    lambda2 = R::rgamma(3 + x2, 1 / (n - theta + alpha));
    alpha = R::rgamma(16, 1 / (10 + lambda1 + lambda2));
    theta = rtheta(lambda1,lambda2,X);
  }

  for (int i = 0; i < iter; i++) {
    double x1 = sum(X[Rcpp::Range(0, theta - 1)]);
    double x2 = sum(X[Rcpp::Range(theta, n - 1)]);
    lambda1 = R::rgamma(3 + x1, 1 / (theta + alpha));
    lambda2 = R::rgamma(3 + x2, 1 / (n - theta + alpha));
    alpha = R::rgamma(16, 1 / (19 + lambda1 + lambda2));
    theta = rtheta(lambda1,lambda2,X);
    res(i, 0) = lambda1;
    res(i, 1) = lambda2;
    res(i, 2) = alpha;
    res(i, 3) = theta;
  }

  return (res);
}

// [[Rcpp::export]]
double halfnormal(double x, double sigma, double theta, double s) {
  double halfnormal_part = exp(-pow(x, 2.0) / (2 * pow(sigma, 2.0)));
  double poisson_part = pow(x, s) * exp(-x * theta);
  return (poisson_part * halfnormal_part);
}

// [[Rcpp::export]]
double MH_ratio(double current, double candidate, double target, double theta, double s) {
  double ratio = (halfnormal(candidate, target, theta, s) * R::dgamma(current, 3, 1/5, 0)) /
    (halfnormal(current, target, theta, s) * R::dgamma(candidate, 3, 1/5, 0));
  return (ratio);
}




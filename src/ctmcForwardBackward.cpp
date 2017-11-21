#include<Rcpp.h>
using namespace std;
using namespace Rcpp;

void forwardLoop2(NumericVector& ts, double theta, NumericMatrix& obs, NumericMatrix& alpha)
{
  double T = (double)ts.size(); // col
  double ns = (double)obs.nrow(); // row
  
  for (int i = 0; i < ns; ++i)
    alpha(i, 0) = 1 / ns * obs(i, 0);
  
  for (int t = 1; t < T; ++t) {
    vector<double> tmp1, tmp2;
    double tmp = 1 - exp(-(ts[t] - ts[t-1]) * theta);
    double p0 = 1 - (ns-1)/ns * tmp;
    double p1 = 1/ns * tmp;
    for (int j = 0; j < ns; ++j){
      tmp1.push_back(alpha(j, t-1) * p0);
      tmp2.push_back(alpha(j, t-1) * p1);
    }
    double sum1 = 0, sum2 = 0;
    for (int i = 0; i < (int)tmp1.size(); ++i) {
      sum1 += tmp1[i];
      sum2 += tmp2[i];
    }
    for (int i = 0; i < ns; i++)
      alpha(i, t) = (sum2 - tmp2[i] + tmp1[i]) * obs(i, t);
  }
}

void backwardLoop(NumericVector& ts, double theta, NumericMatrix& obs, NumericMatrix& beta)
{
  double T = (double)ts.size(); // col
  double ns = (double)obs.nrow(); // row
  
  for (int i= 0; i < ns; ++i)
    beta(i, T-1) = 1;
  
  for (int t = T-2; t >= 0; --t) {
    double tmp = 1 - exp(-(ts[t + 1] - ts[t]) * theta);
    double p0 = 1 - (ns-1)/ns * tmp;
    double p1 = 1/ns * tmp;
    vector<double> tmp1, tmp2;
    for (int j = 0; j < ns; ++j) {
      tmp1.push_back(p0 * beta(j, t + 1) * obs(j, t + 1));
      tmp2.push_back(p1 * beta(j, t + 1) * obs(j, t + 1));
    }
    double sum1 = 0, sum2 = 0;
    for (int i = 0; i < (int)tmp1.size(); ++i) {
      sum1 += tmp1[i];
      sum2 += tmp2[i];
    }
    for (int i = 0; i < ns; ++i)
      beta(i, t) = sum2 - tmp2[i] + tmp1[i];
  }
}

// Do not forget to add documentation for your package using roxygen2
// [[Rcpp::export]]
NumericMatrix ctmcForwardBackward(NumericVector ts, double theta, NumericMatrix obs) {
  int m = (int)ts.size();
  int n = (int)obs.nrow();
  
  if ( obs.ncol() != m )
    stop("The input matrix does not conform to the other parameters");
  
  NumericMatrix condProb(n,m);
  
  NumericMatrix alpha(n, m);
  NumericMatrix beta(n, m);
  
  forwardLoop2(ts, theta, obs, alpha);
  backwardLoop(ts, theta, obs, beta);
  
  for (int t = 0; t < m; ++t) {
    double sum = 0;
    
    for (int i = 0; i < n; ++i)
      sum += (alpha(i, t) * beta(i, t));
    for (int i = 0; i < n; ++i)
      condProb(i, t) = alpha(i, t) * beta(i, t) / sum;
  }
  
  return condProb;
}
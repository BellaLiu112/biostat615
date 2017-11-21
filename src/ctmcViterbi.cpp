#include<Rcpp.h>
using namespace std;
using namespace Rcpp;

//[[Rcpp::export]]
IntegerVector ctmcViterbi(NumericVector ts, double theta, NumericMatrix obs)
{
  double m = (double)ts.size(); // col
  double n = (double)obs.nrow(); // row
  
  if ( obs.ncol() != m )
    stop("The input matrix does not conform to the other parameters");
  
  NumericMatrix delta(n, m+1);
  IntegerMatrix phi(n ,m);
  
  NumericVector tm;
  tm.push_back(0);
  for (int i = 0; i < ts.size(); ++i) {
    tm.push_back(ts[i]);
  }
  for (int i = 0; i < n; ++i){
    delta(i, 0) = 1/n;
  }
  vector<double> tmp1, tmp2;
  for (int t = 1; t < m + 1; ++t) {
    double tmp = 1 - exp(-(tm[t] - tm[t-1]) * theta);
    double p0 = 1 - (n-1)/n * tmp;
    double p1 = 1/n * tmp;
    tmp1.clear();
    tmp2.clear();
    
    for (int j = 0; j < n; ++j) {
      tmp1.push_back(p0 * delta(j, t-1));
      tmp2.push_back(p1 * delta(j, t-1));
    }
    
    // find the largest value index in tmp2
    int index = 0;
    double maxTmp2 = tmp2[0];
    for (int i = 0; i < tmp2.size(); ++i) {
      if (tmp2[i] > maxTmp2) {
        index = i;
        maxTmp2 = tmp2[i];
      }
    }
    if (p0 > p1) {
      for (int j = 0; j < n; ++j) {
        delta(j, t) = (tmp1[j] > maxTmp2 ? tmp1[j] : maxTmp2) * obs(j, t-1);
        phi(j, t-1) = tmp1[j] > maxTmp2 ? j : index;
      }
    } else {
      for (int j = 0; j < n; ++j) {
        if (tmp1[j] > maxTmp2) {
          delta(j, t) = tmp1[j] * obs(j, t-1);
          phi(j, t-1) = j;
        } else {
          if (maxTmp2 != tmp2[j]) {
            delta(j, t) = maxTmp2 * obs(j, t-1);
            phi(j, t-1) = index;
          } else {
            double leftMax = tmp2[0], rightMax = tmp2[j + 1];
            int leftIndex = 0, rightIndex = j + 1;
            for (int i = 0; i < j; ++i) {
              if (tmp2[i] > leftMax) {
                leftMax = tmp2[i];
                leftIndex = i;
              }
            }
            for (int i = j + 1; i < tmp2.size(); ++i) {
              if (tmp2[i] > rightMax) {
                rightMax = tmp2[i];
                rightIndex = i;
              }
            }
            double secondMax = leftMax > rightMax ? leftMax : rightMax;
            int secondIndex = leftMax > rightMax ? leftIndex : rightIndex;
            if (tmp2[j] > secondMax) {
              delta(j, t) = tmp2[j] * obs(j, t-1);
              phi(j, t-1) = j;
            } else {
              delta(j, t) = secondMax * obs(j, t-1);
              phi(j, t-1) = secondIndex;
            }
          }
        }
      }
    }
  }
  IntegerVector viterbiPath(m, 0);
  double ml = -1;
  for (int i = 0; i < n; ++i) {
    if (delta(i, m-1) > ml) {
      viterbiPath[m-1] = i;
      ml = delta(i, m-1);
    }
  }
  for (int i = m-1; i > 0; i--)
    viterbiPath[i-1] = phi(viterbiPath[i], i);
  
  return viterbiPath;
}


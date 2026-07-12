// src/sample_entropy.cpp — Sample Entropy (SampEn) core
//
// Direct C++ port of the counting core in mse.c (Costa, PhysioNet),
// validated to reproduce the compiled reference binary's output exactly
// (to its own displayed precision) on synthetic test data. See
// inst/COPYRIGHTS.
//
// Original copyright notice (preserved per GPL requirements):
//   mse.c   M. Costa   1 August 2004, last revised 4 August 2004 (GBM)
//   Copyright (C) 2004 Madalena Costa
//   Licensed under the GNU General Public License, version 2, or (at
//   your option) any later version. See LICENSE.md / inst/COPYRIGHTS.
//   Original method: Costa M, Goldberger AL, Peng CK. Multiscale entropy
//   analysis of complex physiologic time series. Phys Rev Lett
//   2002;89:068102.

#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Returns cont[1..m+1] as a length-(m+1) vector: out[0] = cont[1] (count of
// length-1 template matches, i.e. pairs (i,l), i<l, with |y[i]-y[l]| <=
// tolerance), ..., out[m] = cont[m+1] (length-(m+1) template matches).
// Self-matches (i == l) are excluded, matching the original algorithm.

//' @keywords internal
// [[Rcpp::export]]
NumericVector sample_entropy_counts_cpp(NumericVector y, int m, double tolerance) {
  int N = y.size();
  int nlin_m = N - m;

  NumericVector out(m + 1, 0.0);
  if (nlin_m < 2) return out;

  std::vector<long long> cont(m + 2, 0);

  for (int i = 0; i < nlin_m; ++i) {
    for (int l = i + 1; l < nlin_m; ++l) {
      int k = 0;
      while (k < m && std::fabs(y[i + k] - y[l + k]) <= tolerance) {
        ++k;
        cont[k]++;
      }
      if (k == m && std::fabs(y[i + m] - y[l + m]) <= tolerance) {
        cont[m + 1]++;
      }
    }
  }

  for (int k = 1; k <= m + 1; ++k) out[k - 1] = (double)cont[k];
  return out;
}

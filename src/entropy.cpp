// src/entropy.cpp — Entropy measures
//
// Ported from Lucas França's own mrpheus package (src/staging_features.cpp,
// perm_entropy_cpp), which was itself validated for exact parity against
// the antropy Python library (antropy.perm_entropy) as part of mrpheus's
// AASM staging feature pipeline. Re-validated here directly against
// antropy 0.2.2 on synthetic test data (exact match to displayed
// precision). Same author, no licensing conflict porting between own
// packages; see inst/COPYRIGHTS.
//
// Reference: Bandt C, Pompe B. Permutation entropy: a natural complexity
// measure for time series. Phys Rev Lett 2002;88:174102.

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

//' @keywords internal
// [[Rcpp::export]]
double perm_entropy_raw_cpp(NumericVector x, int order = 3, int delay = 1) {
  int N = x.size();
  int n = N - (order - 1) * delay;
  if (n <= 0) return NA_REAL;

  std::vector<int> fact(order, 1);
  for (int i = 1; i < order; i++) fact[i] = fact[i - 1] * i;

  int n_perms = fact[order - 1] * order;
  std::vector<int> counts(n_perms, 0);
  std::vector<int> argsort(order);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < order; j++) argsort[j] = j;
    std::stable_sort(argsort.begin(), argsort.end(),
      [&](int a, int b) { return x[i + a * delay] < x[i + b * delay]; });
    int hash = 0;
    for (int j = 0; j < order - 1; j++) {
      int cnt = 0;
      for (int k = j + 1; k < order; k++) if (argsort[k] < argsort[j]) cnt++;
      hash += cnt * fact[order - 1 - j];
    }
    counts[hash]++;
  }

  double h = 0.0;
  for (int k = 0; k < n_perms; k++) {
    if (counts[k] > 0) {
      double p = (double)counts[k] / (double)n;
      h -= p * std::log(p);
    }
  }
  return h; // raw Shannon entropy in nats; R wrapper normalises if requested
}

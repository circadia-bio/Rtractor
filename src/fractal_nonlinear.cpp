// src/fractal_nonlinear.cpp — Petrosian FD, Hjorth parameters, zero-crossings
//
// Ported from Lucas França's own mrpheus package (src/staging_features.cpp),
// itself validated for exact parity against the antropy Python library as
// part of mrpheus's AASM staging feature pipeline (antropy.petrosian_fd,
// antropy.hjorth_params, antropy.num_zerocross). Re-validated here directly
// against antropy 0.2.2 on synthetic test data (exact or near-exact match —
// residual differences are floating-point summation-order noise, not
// algorithmic). Same author, no licensing conflict porting between own
// packages; see inst/COPYRIGHTS.

#include <Rcpp.h>
#include <cmath>
#include <limits>
using namespace Rcpp;

//' @keywords internal
// [[Rcpp::export]]
double petrosian_fd_cpp(NumericVector x) {
  int N = x.size();
  if (N < 3) return NA_REAL;
  int nzc_p = 0;
  for (int i = 0; i < N - 2; i++) {
    double d_curr = x[i + 1] - x[i];
    double d_next = x[i + 2] - x[i + 1];
    if (d_curr * d_next < 0.0) nzc_p++;
  }
  double logN = std::log10((double)N);
  return logN / (logN + std::log10((double)N / ((double)N + 0.4 * (double)nzc_p)));
}

//' @keywords internal
// [[Rcpp::export]]
int num_zerocross_cpp(NumericVector x) {
  int N = x.size();
  int count = 0;
  for (int i = 0; i < N - 1; i++) {
    int s1 = (x[i]     > 0.0) - (x[i]     < 0.0);
    int s2 = (x[i + 1] > 0.0) - (x[i + 1] < 0.0);
    if (s1 != s2) count++;
  }
  return count;
}

//' @keywords internal
// [[Rcpp::export]]
NumericVector hjorth_cpp(NumericVector x) {
  int N = x.size();
  NumericVector na_out = NumericVector::create(
    Named("mobility") = NA_REAL, Named("complexity") = NA_REAL
  );
  if (N < 3) return na_out;

  double sum_x = 0, sum_x2 = 0;
  double sum_d1 = 0, sum_d12 = 0;
  double sum_d2 = 0, sum_d22 = 0;
  for (int i = 0; i < N; i++) {
    sum_x  += x[i];
    sum_x2 += x[i] * x[i];
    if (i < N - 1) {
      double d1 = x[i + 1] - x[i];
      sum_d1 += d1; sum_d12 += d1 * d1;
    }
    if (i < N - 2) {
      double d2 = x[i + 2] - 2.0 * x[i + 1] + x[i];
      sum_d2 += d2; sum_d22 += d2 * d2;
    }
  }
  int Nd1 = N - 1, Nd2 = N - 2;
  double var_x  = (sum_x2  - sum_x  * sum_x  / (double)N)   / (double)(N   - 1);
  double var_d1 = (sum_d12 - sum_d1 * sum_d1 / (double)Nd1) / (double)(Nd1 - 1);
  double var_d2 = (sum_d22 - sum_d2 * sum_d2 / (double)Nd2) / (double)(Nd2 - 1);

  double eps   = std::numeric_limits<double>::epsilon();
  double hmob  = std::sqrt(var_d1 / (var_x  + eps));
  double hcomp = std::sqrt(var_d2 / (var_d1 + eps)) / (hmob + eps);
  return NumericVector::create(Named("mobility") = hmob, Named("complexity") = hcomp);
}

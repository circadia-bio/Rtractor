// src/higuchi.cpp — Higuchi Fractal Dimension
//
// Clean-room C++ reimplementation of Higuchi's (1988) algorithm, developed
// from the published method rather than translated from any licensed
// source. A MATLAB reference implementation (Higuchi_FD.m, J. Monge
// Alvarez, Univ. of Valladolid) was consulted for validation only; see
// inst/COPYRIGHTS. Validated against a transliteration of that reference
// on synthetic test data (max abs. difference ~1e-11).
//
// Reference: Higuchi T. Approach to an irregular time series on the basis
// of the fractal theory. Physica D 1988;31(2):277-283.

#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

//' @keywords internal
// [[Rcpp::export]]
NumericVector higuchi_length_cpp(NumericVector serie, int kmax) {
  int N = serie.size();
  NumericVector L(kmax);

  for (int k = 1; k <= kmax; ++k) {
    double Lk_sum = 0.0;
    for (int m = 1; m <= k; ++m) {
      int limit = (N - m) / k;
      if (limit < 1) stop("Series too short for the requested k_max.");
      int count = limit + 1;

      std::vector<double> aux(count);
      for (int idx = 0; idx < count; ++idx) {
        int i = m + idx * k;
        aux[idx] = serie[i - 1];
      }

      double R = (double)(N - 1) / (double)(limit * k);
      double Lm = 0.0;
      for (int idx = 0; idx < count - 1; ++idx) Lm += std::fabs(aux[idx + 1] - aux[idx]);
      Lm = (Lm * R) / k;
      Lk_sum += Lm;
    }
    L[k - 1] = Lk_sum / k;
  }
  return L;
}

// src/fractal_multifractal.cpp — Multifractal analysis
//
// Clean-room C++ reimplementations from published algorithms (no license
// headers on the reference MATLAB sources consulted during development;
// see inst/COPYRIGHTS). Both cores validated against Python
// transliterations of the originals on synthetic test data (exact match
// to displayed precision).

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// ---- MFDMA: multifractal detrending moving average ------------------------
//
// Reference: Gu GF, Zhou WX. Detrending moving average algorithm for
// multifractals. Phys Rev E 2010;82:011136.
//
// Computes, for each scale n, the segment-level RMS fluctuations after
// detrending the profile y with a moving average of window n. Downstream
// steps (q-th order fluctuation function Fq, scaling exponent tau(q),
// singularity spectrum alpha/f via Legendre transform) are done in R
// (R/fractal.R, mfdma()) since they're small regression loops over q,
// not performance-critical.

//' @keywords internal
// [[Rcpp::export]]
List mfdma_fluctuations_cpp(NumericVector y, IntegerVector scales, double theta) {
  int M = y.size();
  List F(scales.size());

  for (int si = 0; si < scales.size(); ++si) {
    int n = scales[si];
    int len_y1 = M - n + 1;
    if (len_y1 < 1) {
      F[si] = NumericVector(0);
      continue;
    }

    std::vector<double> y1(len_y1);
    double window_sum = 0.0;
    for (int k = 0; k < n; ++k) window_sum += y[k];
    y1[0] = window_sum / n;
    for (int j = 1; j < len_y1; ++j) {
      window_sum += y[j + n - 1] - y[j - 1];
      y1[j] = window_sum / n;
    }

    int offset = std::max(1, (int)std::floor(n * (1.0 - theta))); // 1-indexed, as in the original
    std::vector<double> e(len_y1);
    for (int k = 0; k < len_y1; ++k) {
      int idx = offset - 1 + k; // 0-indexed
      e[k] = y[idx] - y1[k];
    }

    int n_seg = len_y1 / n;
    NumericVector Fi(n_seg);
    for (int k = 0; k < n_seg; ++k) {
      double sumsq = 0.0;
      for (int idx = k * n; idx < (k + 1) * n; ++idx) sumsq += e[idx] * e[idx];
      Fi[k] = std::sqrt(sumsq / n);
    }
    F[si] = Fi;
  }

  return F;
}

// ---- Chhabra-Jensen: multifractal spectrum via box-counting ---------------
//
// Reference: Chhabra A, Jensen RV. Direct determination of the f(alpha)
// singularity spectrum. Phys Rev Lett 1989;62:1327-1330.
//
// Computes, for each (q, scale) pair, the Ma/Mf/Md moments used to recover
// alpha(q), f(alpha(q)), and the generalised dimension D(q) via linear
// regression against -log10(2^scale) (done in R, R/fractal.R,
// chhabra_jensen()).

//' @keywords internal
// [[Rcpp::export]]
List chhabra_jensen_moments_cpp(NumericVector timeseries, NumericVector q_values,
                                 IntegerVector scales) {
  int L = timeseries.size();
  int nq = q_values.size();
  int ns = scales.size();

  NumericMatrix Ma(nq, ns), Mf(nq, ns), Md(nq, ns);

  double total_sum = 0.0;
  for (int i = 0; i < L; ++i) total_sum += timeseries[i];

  for (int j = 0; j < ns; ++j) {
    int window = 1 << scales[j];
    if (L % window != 0) {
      stop("Length of `x` is not divisible by 2^scale for scale = %d.", scales[j]);
    }
    int nrows = L / window;

    std::vector<double> p(window);
    for (int k = 0; k < window; ++k) {
      double s = 0.0;
      for (int r = 0; r < nrows; ++r) s += timeseries[k * nrows + r];
      p[k] = s / total_sum;
    }

    for (int i = 0; i < nq; ++i) {
      double q = q_values[i];
      double Nor = 0.0;
      for (int k = 0; k < window; ++k) Nor += std::pow(p[k], q);

      if (q > 0 && q <= 1) {
        double s = 0.0;
        for (int k = 0; k < window; ++k) s += p[k] * std::log10(p[k]);
        Md(i, j) = s / Nor;
      } else {
        Md(i, j) = std::log10(Nor);
      }

      double ma = 0.0, mf = 0.0;
      for (int k = 0; k < window; ++k) {
        double mu = std::pow(p[k], q) / Nor;
        ma += mu * std::log10(p[k]);
        mf += mu * std::log10(mu);
      }
      Ma(i, j) = ma;
      Mf(i, j) = mf;
    }
  }

  return List::create(_["Ma"] = Ma, _["Mf"] = Mf, _["Md"] = Md);
}

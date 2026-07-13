// src/sda.cpp — Standard Deviation Analysis (SDA) fluctuation function
//
// Windowed Root-Mean-Square deviation from the local mean, no detrending
// (Russ 1994). Ported from Lucas Franca's own legacy `mrug` C++ tool
// (calcM()/localM() code path), which was used to compute the body-sway
// Hurst exponents in:
//   da Silva Costa I, Gamundi A, Miranda JGV, Franca LGS, De Santana CN,
//   Montoya P. Altered functional performance in patients with
//   fibromyalgia. Front Hum Neurosci 2017;11:14.
//
// Unlike every other file in this table, this is NOT a faithful port --
// three bugs found in the available copy of `mrug` were fixed rather than
// reproduced, since the exact code snapshot used for the 2017 paper could
// not be identified. See inst/COPYRIGHTS for the full write-up:
//   1. this->k / local k member-shadowing bug (undefined first-window size
//      at every scale) -- fixed by using a single, correctly-scoped scale
//      value throughout.
//   2. read() eof()-loop off-by-one (duplicated final row) -- not
//      applicable here; R handles I/O, this file only computes fluctuations.
//   3. windowing always started at index 1, dropping the series' first
//      sample -- fixed by windowing from index 0.
//
// The geometric scale-stepping scheme itself (growth factor 1.1, forced
// +1 increments below window size 20) is preserved as-is (R/fractal.R,
// .sda_scales()) -- that's a deliberate sampling design, not a bug.

#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' @keywords internal
// [[Rcpp::export]]
NumericVector sda_fluctuation_cpp(NumericVector x, IntegerVector scales) {
  int N = x.size();
  int n_scales = scales.size();
  NumericVector F(n_scales);

  for (int si = 0; si < n_scales; ++si) {
    int k = scales[si];
    if (k < 2 || k > N) {
      F[si] = NA_REAL;
      continue;
    }

    // Running sums for the current window, updated incrementally as the
    // window slides (O(N) per scale rather than O(N*k)) -- needed since
    // these signals can run to tens of thousands of samples (the original
    // mrug.h sized its buffers for 60005 points).
    double sum_y = 0.0, sum_y2 = 0.0;
    for (int j = 0; j < k; ++j) {
      sum_y  += x[j];
      sum_y2 += x[j] * x[j];
    }

    double sum_s = 0.0;
    int ih = 0;
    int n_windows = N - k + 1;

    for (int i = 0; i < n_windows; ++i) {
      if (i > 0) {
        double out = x[i - 1];
        double in  = x[i + k - 1];
        sum_y  += in - out;
        sum_y2 += in * in - out * out;
      }
      // sum((y - mean)^2) == sum_y2 - sum_y^2/k, the standard identity.
      // Clamped at 0 to guard against benign floating-point cancellation
      // noise (which would otherwise occasionally produce a tiny negative
      // value for a near-constant window).
      double dh = sum_y2 - (sum_y * sum_y) / (double)k;
      if (dh < 0.0) dh = 0.0;
      if (dh > 0.0) {
        ih++;
        sum_s += std::sqrt(dh / (double)k);
      }
    }
    F[si] = (ih > 0) ? (sum_s / (double)ih) : NA_REAL;
  }
  return F;
}

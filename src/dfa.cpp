// src/dfa.cpp — Detrended Fluctuation Analysis
//
// This is a direct C++ port of dfa.c (v4.9), preserving its algorithm and
// numerical behaviour exactly (validated against the original compiled
// binary on synthetic test data — identical output to printed precision).
// See inst/COPYRIGHTS for full attribution and citation requirements.
//
// Original copyright notice (preserved per GPL requirements):
//   dfa.c   J. Mietus, C-K Peng, and G. Moody   8 February 2001
//   Copyright (C) 2001-2005 Joe Mietus, C-K Peng, and George B. Moody
//   Licensed under the GNU General Public License, version 2, or (at
//   your option) any later version. See LICENSE / inst/COPYRIGHTS.
//   Original method: Peng et al., Phys Rev E 1994;49:1685-1689.

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// Fits a polynomial of degree (nfit - 1) to y[0..boxsize-1] against
// x = 1..boxsize via OLS normal equations (Gauss-Jordan with partial
// pivoting), matching the original's polyfit(). Returns the residual
// sum of squares (chisq).
static double dfa_polyfit(const double* y, long boxsize, int nfit) {
  std::vector<double> A(nfit * nfit, 0.0);
  std::vector<double> b(nfit, 0.0);

  for (long i = 1; i <= boxsize; ++i) {
    std::vector<double> xi(nfit);
    xi[0] = 1.0;
    for (int j = 1; j < nfit; ++j) xi[j] = xi[j - 1] * (double)i;

    double yi = y[i - 1];
    for (int j = 0; j < nfit; ++j) {
      b[j] += xi[j] * yi;
      for (int k = 0; k < nfit; ++k) A[j * nfit + k] += xi[j] * xi[k];
    }
  }

  std::vector<double> M(A), rhs(b);
  for (int col = 0; col < nfit; ++col) {
    int piv = col;
    double best = std::fabs(M[col * nfit + col]);
    for (int r = col + 1; r < nfit; ++r) {
      double v = std::fabs(M[r * nfit + col]);
      if (v > best) { best = v; piv = r; }
    }
    if (piv != col) {
      for (int c = 0; c < nfit; ++c) std::swap(M[col * nfit + c], M[piv * nfit + c]);
      std::swap(rhs[col], rhs[piv]);
    }
    double pivval = M[col * nfit + col];
    if (pivval == 0.0) stop("Singular matrix encountered in DFA polynomial fit.");
    for (int c = 0; c < nfit; ++c) M[col * nfit + c] /= pivval;
    rhs[col] /= pivval;
    for (int r = 0; r < nfit; ++r) {
      if (r == col) continue;
      double factor = M[r * nfit + col];
      if (factor == 0.0) continue;
      for (int c = 0; c < nfit; ++c) M[r * nfit + c] -= factor * M[col * nfit + c];
      rhs[r] -= factor * rhs[col];
    }
  }
  std::vector<double>& beta = rhs;

  double chisq = 0.0;
  for (long i = 1; i <= boxsize; ++i) {
    double p = 1.0, fitted = 0.0;
    for (int j = 0; j < nfit; ++j) { fitted += beta[j] * p; p *= (double)i; }
    double resid = fitted - y[i - 1];
    chisq += resid * resid;
  }
  return chisq;
}

//' @keywords internal
// [[Rcpp::export]]
List dfa_cpp(NumericVector x, int nfit = 2, int minbox = 0, int maxbox = 0,
             bool integrate = true, bool sliding_window = false) {
  long npts = x.size();

  std::vector<double> seq(npts);
  if (integrate) {
    double acc = 0.0;
    for (long i = 0; i < npts; ++i) { acc += x[i]; seq[i] = acc; }
  } else {
    for (long i = 0; i < npts; ++i) seq[i] = x[i];
  }

  long lo_box = minbox, hi_box = maxbox;
  if (lo_box < 2 * nfit) lo_box = 2 * nfit;
  if (hi_box == 0 || hi_box > npts / 4) hi_box = npts / 4;
  if (lo_box > hi_box) {
    std::swap(lo_box, hi_box);
    if (lo_box < 2 * nfit) lo_box = 2 * nfit;
  }

  double boxratio = std::pow(2.0, 1.0 / 8.0);
  int rslen = (int)(std::log10((double)hi_box / (double)lo_box) / std::log10(boxratio) + 1.5);
  std::vector<long> rs(rslen + 2, 0);
  int n = 2;
  rs[1] = lo_box;
  for (int ir = 1; n <= rslen && rs[n - 1] < hi_box; ++ir) {
    long rw = (long)(lo_box * std::pow(boxratio, ir) + 0.5);
    if (rw > rs[n - 1]) { rs[n] = rw; ++n; }
  }
  --n;
  if (rs[n] > hi_box) --n;
  int nr = n;

  if (nr < 1) stop("Series too short for the requested box-size range.");

  NumericVector n_out(nr);
  NumericVector F_out(nr);

  for (int i = 1; i <= nr; ++i) {
    long boxsize = rs[i];
    long inc, stat;
    if (sliding_window) {
      inc = 1;
      stat = (npts - boxsize + 1) * boxsize;
    } else {
      inc = boxsize;
      stat = (npts / boxsize) * boxsize;
    }
    double mse = 0.0;
    for (long j = 0; j <= npts - boxsize; j += inc)
      mse += dfa_polyfit(&seq[j], boxsize, nfit);
    mse /= (double)stat;
    n_out[i - 1] = (double)boxsize;
    F_out[i - 1] = std::sqrt(mse);
  }

  return List::create(_["n"] = n_out, _["F"] = F_out);
}

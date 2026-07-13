// src/embed.cpp — False Nearest Neighbours (embedding dimension estimation)
//
// Clean-room C++ implementation of Kennel, Brown & Abarbanel's (1992)
// false-nearest-neighbours method for choosing an embedding dimension.
// No reference implementation was specified for this port (Lucas had no
// preference) -- implemented directly from the published algorithm. See
// inst/COPYRIGHTS.
//
// Reference: Kennel MB, Brown R, Abarbanel HDI. Determining embedding
// dimension for phase-space reconstruction using a geometrical
// construction. Phys Rev A 1992;45:3403-3411.

#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' @keywords internal
// [[Rcpp::export]]
NumericVector fnn_fraction_cpp(NumericVector x, int tau, IntegerVector dims,
                                double rtol, double atol, int theiler) {
  int N = x.size();
  double Ra = Rcpp::sd(x);  // attractor-size proxy for criterion 2
  int n_dims = dims.size();
  NumericVector frac(n_dims);

  for (int di = 0; di < n_dims; ++di) {
    int m = dims[di];
    int n_vec = N - m * tau;  // need coordinate i + m*tau to exist
    if (n_vec < 2) { frac[di] = NA_REAL; continue; }

    int false_count = 0, valid_count = 0;

    for (int i = 0; i < n_vec; ++i) {
      double best_d2 = -1.0;
      int best_j = -1;
      for (int j = 0; j < n_vec; ++j) {
        if (j == i || std::abs(j - i) <= theiler) continue;
        double d2 = 0.0;
        for (int d = 0; d < m; ++d) {
          double diff = x[i + d * tau] - x[j + d * tau];
          d2 += diff * diff;
        }
        if (best_j < 0 || d2 < best_d2) { best_d2 = d2; best_j = j; }
      }
      if (best_j < 0) continue;

      double Rm = std::sqrt(best_d2);
      double extra = x[i + m * tau] - x[best_j + m * tau];
      valid_count++;

      bool is_false = false;
      if (Rm > 0.0) {
        if (std::fabs(extra) / Rm > rtol) is_false = true;
      } else if (extra != 0.0) {
        is_false = true;  // coincident in m dims but diverge at m+1
      }
      double Rm1 = std::sqrt(best_d2 + extra * extra);
      if (Rm1 / Ra > atol) is_false = true;

      if (is_false) false_count++;
    }
    frac[di] = (valid_count > 0) ? ((double)false_count / (double)valid_count) : NA_REAL;
  }
  return frac;
}

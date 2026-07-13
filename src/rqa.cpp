// src/rqa.cpp — recurrence matrix construction and line-length scanning
//
// Clean-room C++ implementation of standard recurrence quantification
// analysis (RQA) primitives. No reference implementation was specified for
// this port (Lucas had no preference) -- implemented directly from the
// widely-used formulas in the Marwan et al. (2007) review, which are
// consistent across the major RQA toolboxes (Marwan's CRP Toolbox, pyRQA,
// R's crqa/nonlinearTseries). See inst/COPYRIGHTS.
//
// Reference: Marwan N, Romano MC, Thiel M, Kurths J. Recurrence plots for
// the analysis of complex systems. Phys Rep 2007;438(5-6):237-329.

#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

//' @keywords internal
// [[Rcpp::export]]
NumericMatrix distance_matrix_cpp(NumericMatrix traj, int norm_code) {
  // norm_code: 1 = euclidean, 2 = maximum (Chebyshev), 3 = manhattan
  int N = traj.nrow();
  int m = traj.ncol();
  NumericMatrix D(N, N);

  for (int i = 0; i < N; ++i) {
    for (int j = i; j < N; ++j) {
      double d = 0.0;
      if (norm_code == 2) {
        double mx = 0.0;
        for (int k = 0; k < m; ++k) {
          double diff = std::fabs(traj(i, k) - traj(j, k));
          if (diff > mx) mx = diff;
        }
        d = mx;
      } else if (norm_code == 3) {
        double s = 0.0;
        for (int k = 0; k < m; ++k) s += std::fabs(traj(i, k) - traj(j, k));
        d = s;
      } else {
        double s = 0.0;
        for (int k = 0; k < m; ++k) {
          double diff = traj(i, k) - traj(j, k);
          s += diff * diff;
        }
        d = std::sqrt(s);
      }
      D(i, j) = d;
      D(j, i) = d;
    }
  }
  return D;
}

//' @keywords internal
// [[Rcpp::export]]
List rqa_line_stats_cpp(LogicalMatrix R, int theiler_window) {
  int N = R.nrow();

  // Mask the Theiler band (main diagonal +/- theiler_window) out of every
  // count below -- standard practice, since points there are trivially
  // recurrent (temporal, not dynamical, proximity).
  std::vector<std::vector<bool>> M(N, std::vector<bool>(N));
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      M[i][j] = R(i, j) && (std::abs(i - j) > theiler_window);
    }
  }

  long long total_points = 0;
  int n_excluded = 0;
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      if (M[i][j]) total_points++;
      if (std::abs(i - j) <= theiler_window) n_excluded++;
    }
  }

  std::vector<int> diag_lengths;
  for (int k = -(N - 1); k <= N - 1; ++k) {
    int i0 = std::max(0, -k);
    int j0 = i0 + k;
    int run = 0;
    for (int i = i0, j = j0; i < N && j < N; ++i, ++j) {
      if (M[i][j]) {
        run++;
      } else {
        if (run > 0) diag_lengths.push_back(run);
        run = 0;
      }
    }
    if (run > 0) diag_lengths.push_back(run);
  }

  std::vector<int> vert_lengths;
  for (int j = 0; j < N; ++j) {
    int run = 0;
    for (int i = 0; i < N; ++i) {
      if (M[i][j]) {
        run++;
      } else {
        if (run > 0) vert_lengths.push_back(run);
        run = 0;
      }
    }
    if (run > 0) vert_lengths.push_back(run);
  }

  return List::create(
    Named("total_points") = (double)total_points,
    Named("n_excluded")   = n_excluded,
    Named("diag_lengths") = wrap(diag_lengths),
    Named("vert_lengths") = wrap(vert_lengths),
    Named("N")            = N
  );
}

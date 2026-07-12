// src/microstates.cpp — Recurrence microstates entropy (maximum-entropy
// threshold search)
//
// Direct C++ port of the core computation in Max_Entropy() from
// circadia-bio/maxEntropy (Entropy.jl), MIT licensed. See inst/COPYRIGHTS.
// Random pair sampling is done in R (recurrence_microstate_entropy(), in
// R/rqa.R) so results are reproducible via set.seed(); this file implements
// only the deterministic nested-loop core, which was validated against a
// Python transliteration of the original given identical index inputs
// (exact match on synthetic test data).
//
// Method: Corso G, Prado TL, dos Santos Lima GZ, Kurths J, Lopes SR.
// Quantifying entropy using recurrence matrix microstates. Chaos
// 2018;28(8):083108.

#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

//' @keywords internal
// [[Rcpp::export]]
List recurrence_microstate_entropy_cpp(NumericVector serie, IntegerVector x_idx,
                                        IntegerVector y_idx, int block,
                                        double eps_min, double eps_max,
                                        int frac, int frac2) {
  int n_samples = x_idx.size();
  int max_micro = 1 << (block * block);

  std::vector<double> stats(max_micro, 0.0);
  std::vector<double> stats_max(max_micro, 0.0);

  std::vector<long long> pow_vec(block * block);
  for (int i = 0; i < block * block; ++i) pow_vec[i] = 1LL << i;

  double S_max = 0.0, threshold_max = 0.0;
  double threshold = eps_min;
  double var_eps = (eps_max - eps_min) / frac;

  for (int i = 1; i <= frac2; ++i) {
    if (i > 1) {
      threshold = threshold_max - var_eps;
      var_eps = 2.0 * var_eps / frac;
    }
    for (int j = 1; j <= frac; ++j) {
      std::fill(stats.begin(), stats.end(), 0.0);
      for (int count = 0; count < n_samples; ++count) {
        long long add = 0;
        int xi = x_idx[count], yi = y_idx[count];
        for (int cy = 0; cy < block; ++cy) {
          for (int cx = 0; cx < block; ++cx) {
            double diff = std::fabs(serie[xi + cx] - serie[yi + cy]);
            int bit = (diff <= threshold) ? 1 : 0;
            add += (long long)bit * pow_vec[cx + cy * block];
          }
        }
        stats[add] += 1.0;
      }
      double S = 0.0;
      for (int k = 0; k < max_micro; ++k) {
        if (stats[k] > 0) {
          double p = stats[k] / n_samples;
          S += -(p * std::log(p));
        }
      }
      if (S > S_max) {
        S_max = S;
        threshold_max = threshold;
        for (int k = 0; k < max_micro; ++k) stats_max[k] = stats[k] / n_samples;
      }
      threshold += var_eps;
    }
  }

  NumericVector stats_max_out(stats_max.begin(), stats_max.end());
  return List::create(_["microstate_probs"] = stats_max_out,
                       _["entropy_max"] = S_max,
                       _["eps_max"] = threshold_max);
}

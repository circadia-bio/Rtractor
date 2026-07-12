# Recurrence quantification analysis (RQA)
#
# Family: rqa
# Wraps established reference implementations for recurrence plots and their
# derived quantification measures. Signal-agnostic (see entropy.R note).
#
# Planned functions (reference implementation to confirm before wrapping):
#   - recurrence_matrix()    Binary/distance recurrence matrix from a
#                            reconstructed trajectory
#   - rqa_measures()         Determinism, laminarity, recurrence rate,
#                            trapping time, mean/max diagonal line length,
#                            entropy of diagonal line lengths, etc.
#   - plot_recurrence()      Recurrence plot visualisation (Rtractor palette)
#
# Requires phase-space reconstruction first (embedding dimension + time
# delay) — see embed.R.
#
# Native counterpart(s) expected in: src/rqa.cpp

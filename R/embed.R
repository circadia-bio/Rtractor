# Phase-space embedding utilities
#
# Family: embed
# Shared reconstruction utilities consumed by the lyapunov and rqa families
# (and optionally fractal, for box-counting on reconstructed attractors).
# Kept as its own family rather than duplicated internal helpers, since
# multiple downstream families depend on the same reconstruction step.
#
# Planned functions (reference implementation to confirm before wrapping):
#   - embed_time_series()    Time-delay (Takens) embedding: builds an
#                            (n x m) matrix from dimension m and delay tau
#   - estimate_delay()       Time delay tau via first minimum of mutual
#                            information (or autocorrelation, as fallback)
#   - estimate_embed_dim()   Embedding dimension m via false nearest
#                            neighbours (Kennel et al. 1992)
#
# Native counterpart(s) expected in: src/embed.cpp

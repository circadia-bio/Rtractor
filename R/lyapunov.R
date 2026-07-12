# Lyapunov exponents
#
# Family: lyapunov
# Wraps established reference implementations for estimating sensitivity to
# initial conditions from a reconstructed state-space trajectory.
# Signal-agnostic (see entropy.R note).
#
# Planned functions (reference implementation to confirm before wrapping):
#   - lyapunov_rosenstein()  Largest Lyapunov exponent, Rosenstein et al. 1993
#   - lyapunov_wolf()        Largest Lyapunov exponent, Wolf et al. 1985
#   - lyapunov_spectrum()    Full Lyapunov spectrum (if reference code
#                            supports it; scope TBD)
#
# All functions in this family require phase-space reconstruction first
# (embedding dimension + time delay) — see embed.R.
#
# Native counterpart(s) expected in: src/lyapunov.cpp

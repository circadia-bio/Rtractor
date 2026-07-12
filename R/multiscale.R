# Multiscale metrics
#
# Family: multiscale
# Wraps established reference implementations for complexity measures
# computed across a range of temporal scales (coarse-graining). Signal-
# agnostic (see entropy.R note).
#
# Planned functions (reference implementation to confirm before wrapping):
#   - mse()                  Multiscale entropy, Costa et al. 2002
#   - rcmse()                Refined composite multiscale entropy,
#                            Wu et al. 2014
#   - mse_complexity_index() Area-under-curve summary of an MSE profile
#
# Built on the entropy family (entropy.R) applied at each coarse-grained
# scale, rather than a separate algorithmic base.
#
# Native counterpart(s) expected in: src/multiscale.cpp

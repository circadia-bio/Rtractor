# Entropy measures
#
# Family: entropy
# Wraps established reference implementations for entropy-based complexity
# measures on scalar time series. Signal-agnostic: accepts any numeric
# vector regardless of source (EEG, actigraphy, BOLD, HRV, etc.), consistent
# with the ecosystem's staging-agnosticism principle.
#
# Planned functions (reference implementation to confirm before wrapping):
#   - sample_entropy()       Sample entropy (SampEn), Richman & Moorman 2000
#   - approx_entropy()       Approximate entropy (ApEn), Pincus 1991
#   - perm_entropy()         Permutation entropy, Bandt & Pompe 2002
#   - shannon_entropy()      Shannon entropy of a discretised/binned signal
#   - renyi_entropy()        Renyi entropy, parametrised by order q
#
# Native counterpart(s) expected in: src/entropy.cpp

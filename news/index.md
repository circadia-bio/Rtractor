# Changelog

## Rtractor 0.0.0.9000

### 🧭 Ported from mrpheus

- Added
  [`perm_entropy()`](https://rtractor.circadia-lab.uk/reference/perm_entropy.md),
  [`petrosian_fd()`](https://rtractor.circadia-lab.uk/reference/petrosian_fd.md),
  [`hjorth_parameters()`](https://rtractor.circadia-lab.uk/reference/hjorth_parameters.md),
  and
  [`num_zerocross()`](https://rtractor.circadia-lab.uk/reference/num_zerocross.md)
  — four nonlinear/complexity features centralised from Lucas França’s
  own `mrpheus` package (its AASM staging feature pipeline, itself a
  validated R/C++ port of the YASA/antropy Python feature set).
  Re-validated directly against antropy 0.2.2 during this port: exact
  match for `perm_entropy` (normalized value), `petrosian_fd`, and
  `num_zerocross`. `hjorth_parameters` agrees closely but not exactly by
  design — it uses Bessel-corrected (ddof=1) variance, matching
  mrpheus’s/R’s convention, rather than antropy’s population-variance
  (ddof=0) convention; the difference shrinks as signal length grows.
  Rtractor’s existing
  [`higuchi_fd()`](https://rtractor.circadia-lab.uk/reference/higuchi_fd.md)
  was cross-checked against mrpheus/antropy’s version and found to agree
  exactly for realistic epoch lengths.

### 🌿 Fractal family

- Added [`dfa()`](https://rtractor.circadia-lab.uk/reference/dfa.md) —
  Detrended Fluctuation Analysis. Direct C++ port of PhysioNet’s
  reference `dfa.c` (Mietus, Peng & Moody), validated to reproduce the
  compiled reference binary’s output exactly on synthetic test data.
  Package license changed to `GPL (>= 2)` as a result (see
  `inst/COPYRIGHTS`).
- Added
  [`higuchi_fd()`](https://rtractor.circadia-lab.uk/reference/higuchi_fd.md)
  — Higuchi Fractal Dimension. Clean-room C++ reimplementation from
  Higuchi (1988), validated against a MATLAB reference implementation on
  synthetic test data (max abs. difference ~1e-11).
- Added [`mfdma()`](https://rtractor.circadia-lab.uk/reference/mfdma.md)
  — multifractal detrending moving average (Gu & Zhou 2010). Clean-room
  C++ reimplementation from the published algorithm, segment-fluctuation
  core validated against a Python transliteration of the reference
  MATLAB implementation on synthetic test data (exact match to displayed
  precision).
- Added
  [`chhabra_jensen()`](https://rtractor.circadia-lab.uk/reference/chhabra_jensen.md)
  — multifractal spectrum via the Chhabra-Jensen box-counting method
  (1989). Clean-room C++ reimplementation, moments core validated
  against a Python transliteration of the reference implementation
  (exact match to displayed precision).

### 🔁 RQA family

- Added
  [`recurrence_microstate_entropy()`](https://rtractor.circadia-lab.uk/reference/recurrence_microstate_entropy.md)
  — recurrence microstates maximum-entropy threshold search (Corso et
  al. 2018; Prado et al. 2018). Direct C++ port of the deterministic
  core of `Max_Entropy()` from `circadia-bio/maxEntropy` (Entropy.jl,
  MIT licensed), the script used in Ferre et al. (2024, PLOS ONE).
  Random pair sampling moved to R for set.seed()-reproducibility; core
  loop validated against a Python transliteration of the original (exact
  match on synthetic test data).

### 🚀 Initial scaffold

- Repository initialised: package skeleton (`DESCRIPTION`, `.Rprofile`,
  `R/`, `src/`, `tests/testthat/`) with no exported functions yet.
- Added Rtractor colour palette
  ([`rtractor_palette()`](https://rtractor.circadia-lab.uk/reference/rtractor_palette.md),
  [`rtractor_palettes()`](https://rtractor.circadia-lab.uk/reference/rtractor_palettes.md)),
  discrete/continuous ggplot2 scales
  ([`scale_colour_rtractor()`](https://rtractor.circadia-lab.uk/reference/scale_colour_rtractor.md),
  [`scale_fill_rtractor()`](https://rtractor.circadia-lab.uk/reference/scale_fill_rtractor.md),
  [`scale_colour_rtractor_c()`](https://rtractor.circadia-lab.uk/reference/scale_colour_rtractor_c.md),
  [`scale_fill_rtractor_c()`](https://rtractor.circadia-lab.uk/reference/scale_fill_rtractor_c.md)),
  and
  [`theme_rtractor()`](https://rtractor.circadia-lab.uk/reference/theme_rtractor.md).
- Planned function families stubbed as file headers only (no
  implementations yet, pending reference code inventory): `entropy`,
  `fractal`, `lyapunov`, `multiscale`, `rqa`, `embed`.

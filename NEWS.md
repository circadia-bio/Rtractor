## Rtractor 0.0.0.9000

### 🎲 Simulate family

* Added `pmodel()` — multifractal binomial cascade generator (Meneveau &
  Sreenivasan 1987), with optional Fourier-domain fractional integration
  (Davis et al. 1997). Clean-room R reimplementation of `pmodel.m`
  (Victor Venema). Validated against two RNG-independent structural
  properties: exact constant output at p=0.5, and exact mass
  conservation at every cascade level for any p.

### 🧭 Ported from mrpheus

* Added `perm_entropy()`, `petrosian_fd()`, `hjorth_parameters()`, and
  `num_zerocross()` — four nonlinear/complexity features centralised
  from Lucas França's own `mrpheus` package (its AASM staging feature
  pipeline, itself a validated R/C++ port of the YASA/antropy Python
  feature set). Re-validated directly against antropy 0.2.2 during this
  port: exact match for `perm_entropy` (normalized value), `petrosian_fd`,
  and `num_zerocross`. `hjorth_parameters` agrees closely but not exactly
  by design — it uses Bessel-corrected (ddof=1) variance, matching
  mrpheus's/R's convention, rather than antropy's population-variance
  (ddof=0) convention; the difference shrinks as signal length grows.
  Rtractor's existing `higuchi_fd()` was cross-checked against
  mrpheus/antropy's version and found to agree exactly for realistic
  epoch lengths.

### 🌿 Fractal family

* Added `dfa()` — Detrended Fluctuation Analysis. Direct C++ port of
  PhysioNet's reference `dfa.c` (Mietus, Peng & Moody), validated to
  reproduce the compiled reference binary's output exactly on synthetic
  test data. Package license changed to `GPL (>= 2)` as a result (see
  `inst/COPYRIGHTS`).
* Added `higuchi_fd()` — Higuchi Fractal Dimension. Clean-room C++
  reimplementation from Higuchi (1988), validated against a MATLAB
  reference implementation on synthetic test data (max abs. difference
  ~1e-11).
* Added `mfdma()` — multifractal detrending moving average (Gu & Zhou
  2010). Clean-room C++ reimplementation from the published algorithm,
  segment-fluctuation core validated against a Python transliteration of
  the reference MATLAB implementation on synthetic test data (exact
  match to displayed precision).
* Added `chhabra_jensen()` — multifractal spectrum via the Chhabra-Jensen
  box-counting method (1989). Clean-room C++ reimplementation, moments
  core validated against a Python transliteration of the reference
  implementation (exact match to displayed precision).

### 🔁 RQA family

* Added `recurrence_microstate_entropy()` — recurrence microstates
  maximum-entropy threshold search (Corso et al. 2018; Prado et al. 2018).
  Direct C++ port of the deterministic core of `Max_Entropy()` from
  `circadia-bio/maxEntropy` (Entropy.jl, MIT licensed), the script used
  in Ferre et al. (2024, PLOS ONE). Random pair sampling moved to R for
  set.seed()-reproducibility; core loop validated against a Python
  transliteration of the original (exact match on synthetic test data).

### 🚀 Initial scaffold

* Repository initialised: package skeleton (`DESCRIPTION`, `.Rprofile`,
  `R/`, `src/`, `tests/testthat/`) with no exported functions yet.
* Added Rtractor colour palette (`rtractor_palette()`, `rtractor_palettes()`),
  discrete/continuous ggplot2 scales (`scale_colour_rtractor()`,
  `scale_fill_rtractor()`, `scale_colour_rtractor_c()`,
  `scale_fill_rtractor_c()`), and `theme_rtractor()`.
* Planned function families stubbed as file headers only (no
  implementations yet, pending reference code inventory): `entropy`,
  `fractal`, `lyapunov`, `multiscale`, `rqa`, `embed`.

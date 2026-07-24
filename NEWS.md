## Rtractor (development version)

### 🐛 Bug fixes

* Corrected hex sticker logo proportions to match the ecosystem-standard
  aspect ratio.
* Regenerated favicons from the corrected logo.

## Rtractor 0.1.1  (2026-07)

### 🧭 Embedding family (new)

* Added `embed_time_series()` — time-delay (Takens 1981) phase-space
  reconstruction.
* Added `estimate_delay()` — delay `tau` via first local minimum of
  average mutual information (Fraser & Swinney 1986), with an
  autocorrelation-based fallback.
* Added `estimate_embed_dim()` — embedding dimension `m` via false
  nearest neighbours (Kennel, Brown & Abarbanel 1992).
  All three are clean-room implementations from the published methods
  -- no reference implementation was specified for this family. See
  `inst/COPYRIGHTS`.

### 🔁 RQA family

* Added `recurrence_matrix()` — binary recurrence matrix from a (raw or
  embedded) time series, supporting both fixed-radius and fixed-
  recurrence-rate thresholding (Eckmann, Kamphorst & Ruelle 1987).
* Added `rqa_measures()` — the standard RQA measures computed from a
  recurrence matrix: RR, DET, L_mean, L_max, ENTR, LAM, TT (Zbilut &
  Webber 1992; Marwan et al. 2007 conventions). Both are clean-room
  implementations -- no reference implementation was specified. See
  `inst/COPYRIGHTS`.
* `plot_recurrence()` remains planned, not yet implemented.

### 🌿 Fractal family

* Added `sda()` — Hurst exponent via windowed Root-Mean-Square deviation
  from the local mean, no detrending ("Standard Deviation Analysis";
  Russ 1994). This is the method used to compute the Hurst exponent of
  body-sway time series in da Silva Costa et al. (2017), a fibromyalgia
  gait/balance study co-authored by Lucas França. Ported from França's
  own legacy `mrug` C++ tool (used to produce that paper's figures), with
  three known bugs fixed rather than reproduced, since the exact code
  snapshot used for the 2017 results could not be identified:
    1. An uninitialized class member (`int k`) was shadowed by a
       same-named local variable, so the first sliding window at every
       scale used an undefined window size instead of the intended one.
    2. A `while (!file.eof())` read loop duplicated the final data row
       and inflated the point count by one.
    3. Windowing always started at index 1, not 0, silently dropping the
       series' first sample from every window.
  See `inst/COPYRIGHTS` for the full write-up. Because of this deviation,
  treat `sda()`'s output as a *corrected* re-implementation of the
  published method, not a bit-exact reproduction of the 2017 paper's
  numbers.

## Rtractor 0.1.0  (2026-07)

### 📉 Multiscale family

* Added `multiscale_entropy()` — multiscale entropy (MSE; Costa,
  Goldberger & Peng 2002), the first function in the previously-empty
  multiscale family. Direct C++ port of the counting core in PhysioNet's
  reference `mse.c` (Costa), validated to reproduce the compiled
  reference binary's output exactly (to its own displayed precision) on
  synthetic test data, both at the sample-entropy level and for the full
  multiscale sweep.
* Added `sample_entropy()` — sample entropy (SampEn; Richman & Moorman
  2000), filling a previously-planned gap in the entropy family and
  serving as the building block `multiscale_entropy()` applies at each
  coarse-grained scale.

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

### 🎨 Visual identity

* Hex sticker logo (`man/figures/logo.svg`), matching the Circadia Lab
  ecosystem's visual conventions: a dense woven strange-attractor motif
  (a nod to "Rtractor" = R + attractor) on a dotted-paper background,
  with the Rtractor colour palette carried through the border and
  wordmark. Wordmark baked as vector outlines rather than live text, so
  it renders correctly regardless of font availability in whatever tool
  processes the SVG (favicon generators, etc.).
* pkgdown site live at <https://rtractor.circadia-lab.uk>, with the full
  favicon set generated via `pkgdown::build_favicons()`.

### 📚 Documentation

* `vignette("getting-started")` — tour of every implemented family.
* `vignette("multifractal-methods")` — `mfdma()` vs `chhabra_jensen()`,
  validated against known `pmodel()` ground truth rather than just
  checking both run without error.
* `vignette("entropy-and-complexity")` — `perm_entropy()`/
  `sample_entropy()`/`multiscale_entropy()`, including the classic Costa
  et al. white-noise-vs-correlated-signal multiscale entropy
  demonstration.

## Rtractor 0.0.0.9000

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
* `mfdma()` and `chhabra_jensen()` (multifractal spectrum) reference
  MATLAB sources received but not yet ported — pending.

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

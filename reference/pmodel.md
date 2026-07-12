# Simulate a Multifractal Time Series Using the p-Model

Generates a multiplicative binomial cascade ("p-model") time series
(Meneveau & Sreenivasan 1987), optionally fractionally integrated in
Fourier space to a prescribed power-spectrum slope (Davis et al. 1997)
to produce a continuous, nonstationary series superposed on the
multifractal cascade.

## Usage

``` r
pmodel(n_values = 256L, p = 0.375, slope = NULL, seed = 42L)
```

## Arguments

- n_values:

  Integer. Desired length of the output series. Internally rounded up to
  the next power of two for the cascade construction, then truncated
  back to `n_values`. Default `256`.

- p:

  Numeric in `(0, 1)`. The p-model parameter: values near `0` or `1`
  produce a very peaked (strongly multifractal) series; values near
  `0.5` produce a calm series (`p = 0.5` gives an exactly constant unit
  vector – no multifractality at all). Default `0.375`.

- slope:

  Numeric or `NULL`. If supplied, the p-model series is fractionally
  integrated in Fourier space to have this power-spectrum slope,
  producing a continuous, nonstationary series (slopes between -1 and -3
  give stationary increments; steeper than -3 is differentiable) – see
  Davis et al. (1997) for the regimes. If `NULL` (default), no
  fractional integration is applied.

- seed:

  Integer or `NULL`. Seeds the random sign draws at *every* cascade
  level. The original MATLAB implementation reseeds to a fixed seed
  (`42`) at every level unconditionally, making the cascade fully
  deterministic by design (each level's random signs are literally a
  prefix of the next level's, since the RNG stream restarts from the
  same point every time). This port preserves that same
  reseed-every-level structure, with `seed` (default `42`, matching the
  original) controlling it. Pass `seed = NULL` to draw fresh randomness
  at every level instead (no reseeding) – unlike the original, which has
  no way to disable this.

## Value

If `slope` is `NULL`: a numeric vector of length `n_values`, the p-model
series itself. If `slope` is supplied: a list with

- x:

  The fractionally integrated series.

- y:

  The underlying p-model series before integration.

## Details

Clean-room R reimplementation of `pmodel.m`, written by the late Victor
Venema (no license header present in the source file; see
`inst/COPYRIGHTS`). The underlying random sign draws use R's own RNG
rather than attempting to replicate MATLAB's specific generator
bit-for-bit (not feasible across languages, and not meaningful for a
synthetic-data generator characterised by its cascade *rule* rather than
particular RNG bytes) – see the `seed` argument below for how the
original's reseed-every-level behaviour is preserved. Validated against
two RNG-independent structural properties that hold for *any* random
sign sequence: at `p = 0.5` the cascade is provably an exact constant
vector of `1`s, and total mass is conserved exactly at every cascade
level (`sum(y)` doubles per level) for any `p` – see
`tests/testthat/test-pmodel.R`.

## References

Meneveau C, Sreenivasan KR. Simple multifractal cascade model for fully
developed turbulence. Phys Rev Lett 1987;59:1424-1427.

Davis A, Marshak A, Cahalan R, Wiscombe W. The Landsat scale break in
stratocumulus as a three-dimensional radiative transfer effect:
implications for cloud remote sensing. J Atmos Sci 1997;54(2):241-260.

Venema V, Bachner S, Rust HW, Simmer C. Statistical characteristics of
surrogate data based on geophysical measurements. Nonlin Processes
Geophys 2006;13:449-466. (Requested citation for this code by its
original author, Victor Venema.)

## Examples

``` r
y <- pmodel(1024, p = 0.3)
range(y)
#> [1]  0.006046618 28.925465498

# p = 0.5 always gives a constant series, regardless of seed:
all(pmodel(64, p = 0.5) == 1)
#> [1] TRUE
```

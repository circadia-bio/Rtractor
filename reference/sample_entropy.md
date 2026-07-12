# Sample Entropy (SampEn)

Estimates the complexity of a time series using sample entropy (Richman
& Moorman 2000): the negative log ratio of template matches of length
`m + 1` to template matches of length `m`, using a fixed
Chebyshev-distance tolerance. Direct C++ port of the counting core in
PhysioNet's reference `mse.c` (Costa), validated to reproduce the
compiled reference binary's output exactly (to its own displayed
precision) on synthetic test data. See `inst/COPYRIGHTS`.

## Usage

``` r
sample_entropy(x, m = 2L, r = 0.15)
```

## Arguments

- x:

  Numeric vector. The time series to analyse.

- m:

  Integer \>= 1. Template length. Default `2` (the standard MSE
  convention, per Costa et al.).

- r:

  Numeric \> 0. Tolerance, as a fraction of `x`'s own standard deviation
  (i.e. the actual Chebyshev-distance tolerance used is `sd(x) * r`).
  Default `0.15` (the standard MSE convention).

## Value

A length-1 numeric: the sample entropy. If there are too few valid
template pairs to compare (`length(x) - m < 2`), returns `NA`. If there
are zero matches at either template length (undefined ratio), returns
the conventional fallback `-log(1 / (N*(N-1)))` used by the reference
implementation, where `N = length(x) - m`.

## Details

This is also the building block
[`multiscale_entropy()`](https://rtractor.circadia-lab.uk/reference/multiscale_entropy.md)
applies at each coarse-grained scale – see `R/multiscale.R`.

## References

Richman JS, Moorman JR. Physiological time-series analysis using
approximate entropy and sample entropy. Am J Physiol Heart Circ Physiol
2000;278(6):H2039-H2049.

## Examples

``` r
set.seed(1)
sample_entropy(rnorm(1000))
#> [1] 2.450377
```

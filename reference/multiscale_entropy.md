# Multiscale Entropy (MSE)

Estimates the complexity of a time series across a range of temporal
scales (Costa, Goldberger & Peng 2002): the series is coarse-grained
(non-overlapping block-averaged) at each scale factor, and
[`sample_entropy()`](https://rtractor.circadia-lab.uk/reference/sample_entropy.md)
is computed on each coarse-grained series, using a tolerance held
*fixed* relative to the original (not the coarse-grained) series'
standard deviation – this is the standard MSE convention and essential
for entropy values to be comparable across scales. Direct C++ port of
the coarse-graining + sample-entropy core in PhysioNet's reference
`mse.c` (Costa), validated to reproduce the compiled reference binary's
output exactly (to its own displayed precision) on synthetic test data.
See `inst/COPYRIGHTS`.

## Usage

``` r
multiscale_entropy(x, scale_max = 20L, m = 2L, r = 0.15)
```

## Arguments

- x:

  Numeric vector. The time series to analyse.

- scale_max:

  Integer \>= 1. Largest scale factor to evaluate; MSE is computed at
  every integer scale from `1` to `scale_max`. Default `20` (the
  standard MSE convention).

- m:

  Integer \>= 1. Template length, passed to
  [`sample_entropy()`](https://rtractor.circadia-lab.uk/reference/sample_entropy.md).
  Default `2` (the standard MSE convention).

- r:

  Numeric \> 0. Tolerance as a fraction of the *original* series'
  standard deviation, passed to
  [`sample_entropy()`](https://rtractor.circadia-lab.uk/reference/sample_entropy.md).
  Default `0.15` (the standard MSE convention).

## Value

A list with:

- scale:

  The scale factors evaluated, `1:scale_max`.

- mse:

  Sample entropy at each scale (may contain `NA` at large scales, where
  the coarse-grained series becomes too short to estimate reliably).

- m, r:

  The parameters used, echoed back for reference.

## References

Costa M, Goldberger AL, Peng CK. Multiscale entropy analysis of complex
physiologic time series. Phys Rev Lett 2002;89:068102.

## Examples

``` r
set.seed(1)
res <- multiscale_entropy(rnorm(2000), scale_max = 10)
plot(res$scale, res$mse, type = "b")

```

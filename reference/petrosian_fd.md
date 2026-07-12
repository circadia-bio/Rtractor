# Petrosian Fractal Dimension

Estimates fractal dimension from the rate of sign changes in a time
series' first difference (a fast proxy for signal irregularity). Ported
from Lucas Franca's own `mrpheus` package (AASM staging feature
pipeline), itself validated for exact parity against the `antropy`
Python library; re-validated here directly against `antropy` 0.2.2 on
synthetic test data (exact match). See `inst/COPYRIGHTS`.

## Usage

``` r
petrosian_fd(x)
```

## Arguments

- x:

  Numeric vector. The time series to analyse.

## Value

A length-1 numeric: the Petrosian fractal dimension.

## References

Petrosian A. Kolmogorov complexity of finite sequences and recognition
of different preictal EEG patterns. Proceedings of the Eighth IEEE
Symposium on Computer-Based Medical Systems 1995:212-217.

## Examples

``` r
set.seed(1)
petrosian_fd(rnorm(1000))
#> [1] 1.035499
```

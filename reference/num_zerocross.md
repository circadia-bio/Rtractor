# Number of Zero Crossings

Counts sign changes in a time series – a simple time-domain proxy for
dominant frequency / signal roughness, commonly reported alongside
Hjorth parameters and fractal dimension in EEG complexity work. Ported
from Lucas Franca's own `mrpheus` package (AASM staging feature
pipeline), itself validated for exact parity against the `antropy`
Python library; re-validated here directly against `antropy` 0.2.2 on
synthetic test data (exact match). See `inst/COPYRIGHTS`.

## Usage

``` r
num_zerocross(x)
```

## Arguments

- x:

  Numeric vector. The time series to analyse.

## Value

A length-1 integer: the number of zero crossings.

## Examples

``` r
set.seed(1)
num_zerocross(rnorm(1000))
#> [1] 516
```

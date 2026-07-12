# Hjorth Mobility and Complexity

Computes the Hjorth parameters (Hjorth 1970): mobility, a proxy for mean
frequency derived from the variance ratio of the first difference to the
signal itself; and complexity, a proxy for bandwidth/irregularity
derived from the second difference. Ported from Lucas Franca's own
`mrpheus` package (AASM staging feature pipeline). Uses Bessel-corrected
(`ddof = 1`) variance throughout, deliberately matching R's
[`var()`](https://rdrr.io/r/stats/cor.html) convention – this was
mrpheus's original design choice, and differs slightly from the
`antropy` Python library's population-variance (`ddof = 0`) convention.
The two converge as `length(x)` grows but differ meaningfully for short
series; see `inst/COPYRIGHTS`.

## Usage

``` r
hjorth_parameters(x)
```

## Arguments

- x:

  Numeric vector. The time series to analyse.

## Value

A named list with `mobility` and `complexity`.

## References

Hjorth B. EEG analysis based on time domain properties.
Electroencephalogr Clin Neurophysiol 1970;29(3):306-310.

## Examples

``` r
set.seed(1)
hjorth_parameters(rnorm(1000))
#> $mobility
#> [1] 1.441764
#> 
#> $complexity
#> [1] 1.203156
#> 
```

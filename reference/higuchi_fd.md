# Higuchi Fractal Dimension

Estimates the fractal dimension of a time series using Higuchi's (1988)
curve-length algorithm. A clean-room C++ reimplementation, validated
against a MATLAB reference implementation on synthetic test data — see
`inst/COPYRIGHTS`.

## Usage

``` r
higuchi_fd(x, k_max = 10L)
```

## Arguments

- x:

  Numeric vector. The time series to analyse.

- k_max:

  Integer (\>= 2). Maximum sub-series interval `k`. Choose based on
  where the log-log curve plateaus for your signal (see Doyle et al. for
  guidance in postural-sway / physiological contexts).

## Value

A list with:

- k:

  The `k` values evaluated, `1:k_max`.

- L:

  Average curve length at each `k`.

- hfd:

  The Higuchi Fractal Dimension: negative slope of `log(L)` on
  `log(1/k)`.

## References

Higuchi T. Approach to an irregular time series on the basis of the
fractal theory. Physica D 1988;31(2):277-283.

## Examples

``` r
set.seed(1)
higuchi_fd(rnorm(1000), k_max = 10)
#> $k
#>  [1]  1  2  3  4  5  6  7  8  9 10
#> 
#> $L
#>  [1] 1204.12194  297.33968  131.25692   73.58719   46.49776   34.08144
#>  [7]   23.10696   18.47379   14.18394   11.63339
#> 
#> $hfd
#> [1] 2.015343
#> 
```

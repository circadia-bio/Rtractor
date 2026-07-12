# Multifractal Spectrum via the Chhabra-Jensen Method

Estimates the multifractal spectrum f(alpha) and generalised dimension
spectrum D(q) of a strictly positive time series using the direct
box-counting method of Chhabra & Jensen (1989). Clean-room C++
reimplementation from the published algorithm (the reference MATLAB
implementation consulted, `ChhabraJensen_Yuj_w0.m`, co-authored by L.
Franca, had no license header; see `inst/COPYRIGHTS`). The moments core
was validated against a Python transliteration of that reference on
synthetic test data (exact match to displayed precision).

## Usage

``` r
chhabra_jensen(
  x,
  q_values = c(seq(-10, -0.1, by = 0.1), seq(1, 10, by = 0.1)),
  scales = NULL
)
```

## Arguments

- x:

  Numeric vector, strictly positive (treat it as a measure — e.g. apply
  a sigmoid transform first if your data can be negative). `length(x)`
  must be evenly divisible by `2^scale` for every value in `scales`,
  i.e. dyadic lengths (powers of two) work best.

- q_values:

  Numeric vector of multifractal orders. Per the original author's
  guidance, values strictly between 0 and 1 (exclusive) are numerically
  unstable for this method and best avoided – a warning is issued if any
  are supplied. Default skips that range:
  `c(seq(-10, -0.1, by = 0.1), seq(1, 10, by = 0.1))`.

- scales:

  Integer vector of box-counting scale exponents; the box size at each
  scale is `2^scale`. Default `1:floor(log2(length(x)) - 2)`, which
  keeps at least 4 points per box at the coarsest scale.

## Value

A list with:

- alpha:

  Singularity strength alpha(q).

- falpha:

  Multifractal spectrum f(alpha(q)).

- Dq:

  Generalised dimension spectrum D(q).

- r_squared_alpha, r_squared_falpha, r_squared_Dq:

  R-squared of the linear fit underlying each of the above, per `q` –
  inspect these before trusting a given q value's estimate.

- q:

  The `q_values` used.

- mu_scale, Ma, Mf, Md:

  The underlying regression inputs, included for completeness.

## References

Chhabra A, Jensen RV. Direct determination of the f(alpha) singularity
spectrum. Phys Rev Lett 1989;62:1327-1330.

## Examples

``` r
set.seed(1)
x <- abs(rnorm(1024)) + 0.01
res <- chhabra_jensen(x, scales = 1:6)
plot(res$alpha, res$falpha, type = "b")

```

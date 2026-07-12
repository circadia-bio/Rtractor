# Permutation Entropy

Estimates the complexity of a time series using permutation entropy
(Bandt & Pompe 2002): the Shannon entropy of the distribution of ordinal
patterns ("motifs") of length `order` found in the series. Ported from
Lucas Franca's own `mrpheus` package (part of its AASM sleep-staging
feature pipeline), itself validated for exact parity against the
`antropy` Python library; re-validated here directly against `antropy`
0.2.2 on synthetic test data (exact match to displayed precision). See
`inst/COPYRIGHTS`.

## Usage

``` r
perm_entropy(x, order = 3L, delay = 1L, normalize = TRUE)
```

## Arguments

- x:

  Numeric vector. The time series to analyse.

- order:

  Integer \>= 2. Length of the ordinal pattern (embedding dimension).
  Default `3`.

- delay:

  Integer \>= 1. Time delay between pattern elements. Default `1`.

- normalize:

  Logical. If `TRUE` (default), divide by `log(order!)` so the result
  falls in `[0, 1]`. If `FALSE`, return the raw Shannon entropy in nats
  (natural log). Note this differs from the `antropy` Python library,
  which computes the raw value in bits (log base 2); the normalized
  value is base-independent and matches `antropy` exactly, but the raw
  values are not directly comparable between the two. See
  `inst/COPYRIGHTS`.

## Value

A length-1 numeric: the permutation entropy.

## References

Bandt C, Pompe B. Permutation entropy: a natural complexity measure for
time series. Phys Rev Lett 2002;88:174102.

## Examples

``` r
set.seed(1)
perm_entropy(rnorm(1000))
#> [1] 0.9988401
```

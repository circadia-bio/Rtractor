# Detrended Fluctuation Analysis (DFA)

Estimates the scaling exponent alpha of a time series using Detrended
Fluctuation Analysis (Peng et al. 1994). This is a direct C++ port of
the reference `dfa.c` implementation distributed by PhysioNet (Mietus,
Peng & Moody), validated to reproduce the original compiled binary's
output exactly on synthetic test data — see `inst/COPYRIGHTS`.

## Usage

``` r
dfa(
  x,
  order = 1,
  min_box = NULL,
  max_box = NULL,
  integrate = TRUE,
  sliding_window = FALSE
)
```

## Arguments

- x:

  Numeric vector. The time series to analyse.

- order:

  Integer \>= 1. Order of the polynomial detrending fit (`1` = linear
  detrending, the DFA default; `2` = quadratic, etc.). Corresponds to
  `nfit - 1` in the original `dfa.c`.

- min_box, max_box:

  Integer or `NULL`. Smallest/largest box size to evaluate. Defaults (as
  in the original): `min_box = 2 * (order + 1)`,
  `max_box = length(x) / 4`.

- integrate:

  Logical. If `TRUE` (default), `x` is treated as an increment series
  and cumulatively summed before analysis (the usual DFA convention —
  e.g. pass RR-interval deviations, not a cumulative profile). Set
  `FALSE` if `x` is already an integrated/cumulative profile.

- sliding_window:

  Logical. Use overlapping (sliding) windows instead of non-overlapping
  boxes. Default `FALSE`.

## Value

A list with:

- n:

  Box sizes evaluated (integer vector).

- F:

  RMS fluctuation at each box size (numeric vector).

- alpha:

  The DFA scaling exponent: the slope of `log10(F)` on `log10(n)`.

## References

Peng C-K, Buldyrev SV, Havlin S, Simons M, Stanley HE, Goldberger AL.
Mosaic organization of DNA nucleotides. Phys Rev E 1994;49:1685-1689.

## Examples

``` r
set.seed(1)
dfa(rnorm(2000))
#> $n
#>  [1]   4   5   6   7   8   9  10  11  12  13  15  16  17  19  21  23  25  27  29
#> [20]  32  35  38  41  45  49  54  59  64  70  76  83  91  99 108 117 128 140 152
#> [39] 166 181 197 215 235 256 279 304 332 362 395 431 470
#> 
#> $F
#>  [1] 0.4554567 0.5446745 0.6155053 0.6958255 0.7302992 0.7841880 0.8242549
#>  [8] 0.8766239 0.9070949 0.9452574 1.0106687 1.0071432 1.0714474 1.1490677
#> [15] 1.2219190 1.2730067 1.3295481 1.3719023 1.4229023 1.4231487 1.5385946
#> [22] 1.5277759 1.6423866 1.7074346 1.8878224 2.0854952 2.1509697 2.1512357
#> [29] 2.3388256 2.4742135 2.5473760 2.6644721 2.6056816 2.7942607 2.8933307
#> [36] 3.1581955 3.5089603 3.3040542 3.3747209 3.9294164 4.1017177 4.0889159
#> [43] 4.4931074 4.1049745 4.5669029 4.2412576 4.7121078 4.7095386 4.8189230
#> [50] 4.6130934 5.4615918
#> 
#> $alpha
#> [1] 0.5088395
#> 
```

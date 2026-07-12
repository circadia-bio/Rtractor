# Recurrence Microstates Entropy (maximum-entropy threshold search)

Estimates the complexity of a (windowed) time series using the entropy
of recurrence microstates (Corso et al. 2018), searching over the
vicinity threshold `epsilon` to find the value that *maximises* the
Shannon entropy of the microstate distribution — a parameter-free
alternative to fixing `epsilon` by hand, related to the
optimal-threshold method of Prado et al. (2018).

## Usage

``` r
recurrence_microstate_entropy(
  x,
  block = 3L,
  n_samples = 18000L,
  eps_min = 1e-09,
  eps_max = 0.499999999,
  frac = 10L,
  frac2 = 10L,
  seed = NULL
)
```

## Arguments

- x:

  Numeric vector. The (already windowed) time series segment to analyse.
  In the original study, `x` was min-max normalised to `[0, 1]` per
  window before calling this function — do the same if you want
  comparable `eps_max` values across windows/subjects.

- block:

  Integer. Side length of the square recurrence microstate block; the
  microstate space has `2^(block^2)` possible patterns. Default `3`
  (matching the reference implementation).

- n_samples:

  Integer. Number of random point-pairs sampled to estimate the
  recurrence microstate distribution. Default `18000`.

- eps_min, eps_max:

  Numeric. Search range for the vicinity threshold. Defaults match the
  reference implementation.

- frac, frac2:

  Integer. Number of steps in the coarse search pass and number of
  refinement iterations, respectively. Defaults match the reference
  implementation.

- seed:

  Optional integer. Seeds the random pair sampling for reproducibility.
  The original Julia script reseeds from system entropy on every run
  (`Random.seed!()` with no argument) and is therefore not reproducible
  run-to-run by design; pass a `seed` here if you need reproducible
  results, or leave `NULL` to match the original's
  fresh-randomness-every-call behaviour.

## Value

A list with:

- microstate_probs:

  Probability of each of the `2^(block^2)` microstates at the
  entropy-maximising threshold.

- entropy_max:

  The maximum Shannon entropy found (`S_max`).

- eps_max:

  The vicinity threshold at which `entropy_max` was achieved.

## Details

Ported from `circadia-bio/maxEntropy` (`Entropy.jl`, MIT licensed), the
script used to produce the entropy estimates in Ferre et al. (2024, PLOS
ONE). The deterministic core loop was validated against a Python
transliteration of the original on synthetic test data (exact match on
`S_max` and `eps_max` given identical sampled index pairs). Only the
reusable `Max_Entropy()` routine was ported — the original script's
`main()` (file I/O over named subject directories) is study-specific
scaffolding, not part of the general-purpose method.

## References

Corso G, Prado TL, dos Santos Lima GZ, Kurths J, Lopes SR. Quantifying
entropy using recurrence matrix microstates. Chaos 2018;28(8):083108.

Prado TL, dos Santos Lima GZ, Lobao-Soares B, do Nascimento GC, Corso G,
Fontenele-Araujo J, Kurths J, Lopes SR. Optimizing the detection of
nonstationary signals by using recurrence analysis. Chaos
2018;28(8):085703.

Ferre IBS, Corso G, dos Santos Lima GZ, Lopes SR, Leocadio-Miguel MA,
Franca LGS, de Lima Prado T, Araujo JF. Cycling reduces the entropy of
neuronal activity in the human adult cortex. PLOS ONE
2024;19(10):e0298703.

## Examples

``` r
set.seed(1)
x <- (sin(seq(0, 20 * pi, length.out = 300)) + 1) / 2
recurrence_microstate_entropy(x, seed = 1)
#> $microstate_probs
#>   [1] 5.048889e-01 1.544444e-02 0.000000e+00 1.611111e-03 1.688889e-02
#>   [6] 5.555556e-05 2.055556e-03 4.555556e-03 0.000000e+00 2.000000e-03
#>  [11] 0.000000e+00 9.388889e-03 0.000000e+00 0.000000e+00 0.000000e+00
#>  [16] 4.000000e-03 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#>  [21] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#>  [26] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#>  [31] 0.000000e+00 2.666667e-03 0.000000e+00 0.000000e+00 0.000000e+00
#>  [36] 0.000000e+00 2.500000e-03 0.000000e+00 1.050000e-02 4.388889e-03
#>  [41] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#>  [46] 0.000000e+00 0.000000e+00 2.222222e-04 0.000000e+00 0.000000e+00
#>  [51] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#>  [56] 2.666667e-03 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#>  [61] 0.000000e+00 0.000000e+00 0.000000e+00 3.444444e-03 1.627778e-02
#>  [66] 2.222222e-04 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#>  [71] 0.000000e+00 0.000000e+00 2.388889e-03 4.166667e-03 0.000000e+00
#>  [76] 4.722222e-03 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#>  [81] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#>  [86] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#>  [91] 0.000000e+00 2.166667e-03 0.000000e+00 0.000000e+00 0.000000e+00
#>  [96] 1.177778e-02 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [101] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [106] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [111] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [116] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [121] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [126] 0.000000e+00 0.000000e+00 4.055556e-03 0.000000e+00 0.000000e+00
#> [131] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [136] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [141] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [146] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [151] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [156] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [161] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [166] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [171] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [176] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [181] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [186] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [191] 0.000000e+00 0.000000e+00 2.555556e-03 0.000000e+00 0.000000e+00
#> [196] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [201] 8.833333e-03 3.833333e-03 0.000000e+00 1.666667e-04 0.000000e+00
#> [206] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [211] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [216] 0.000000e+00 0.000000e+00 2.888889e-03 0.000000e+00 3.944444e-03
#> [221] 0.000000e+00 0.000000e+00 0.000000e+00 4.611111e-03 0.000000e+00
#> [226] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [231] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [236] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [241] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [246] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [251] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e-03
#> [256] 1.694444e-02 1.661111e-02 0.000000e+00 0.000000e+00 0.000000e+00
#> [261] 1.111111e-04 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [266] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [271] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [276] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [281] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [286] 0.000000e+00 0.000000e+00 0.000000e+00 2.944444e-03 0.000000e+00
#> [291] 0.000000e+00 0.000000e+00 4.944444e-03 0.000000e+00 4.388889e-03
#> [296] 1.666667e-04 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [301] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [306] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [311] 2.833333e-03 1.088889e-02 0.000000e+00 0.000000e+00 0.000000e+00
#> [316] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 4.166667e-03
#> [321] 5.555556e-05 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [326] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [331] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [336] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [341] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [346] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [351] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [356] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [361] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [366] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [371] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [376] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [381] 0.000000e+00 0.000000e+00 0.000000e+00 1.666667e-04 2.500000e-03
#> [386] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [391] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [396] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [401] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [406] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [411] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [416] 0.000000e+00 9.111111e-03 0.000000e+00 0.000000e+00 0.000000e+00
#> [421] 4.722222e-03 0.000000e+00 1.111111e-04 0.000000e+00 0.000000e+00
#> [426] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [431] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [436] 0.000000e+00 2.833333e-03 0.000000e+00 4.666667e-03 5.111111e-03
#> [441] 0.000000e+00 0.000000e+00 0.000000e+00 1.222222e-03 0.000000e+00
#> [446] 0.000000e+00 0.000000e+00 1.916667e-02 4.555556e-03 0.000000e+00
#> [451] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [456] 0.000000e+00 4.166667e-03 3.888889e-04 0.000000e+00 0.000000e+00
#> [461] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [466] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [471] 0.000000e+00 0.000000e+00 2.500000e-03 1.150000e-02 0.000000e+00
#> [476] 4.888889e-03 0.000000e+00 0.000000e+00 0.000000e+00 2.222222e-04
#> [481] 4.388889e-03 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [486] 0.000000e+00 0.000000e+00 0.000000e+00 2.222222e-04 0.000000e+00
#> [491] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> [496] 0.000000e+00 2.944444e-03 0.000000e+00 0.000000e+00 0.000000e+00
#> [501] 1.183333e-02 0.000000e+00 4.888889e-03 1.666667e-04 4.277778e-03
#> [506] 4.888889e-03 0.000000e+00 1.800000e-02 4.388889e-03 1.666667e-04
#> [511] 1.783333e-02 1.232222e-01
#> 
#> $entropy_max
#> [1] 2.392682
#> 
#> $eps_max
#> [1] 0.19624
#> 
```

# Retrieve an Rtractor palette

Returns a named character vector of hex colour codes for the requested
palette. Suitable for direct use in
[`ggplot2::scale_fill_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html)
or
[`ggplot2::scale_colour_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html).

## Usage

``` r
rtractor_palette(palette = "main", n = NULL, reverse = FALSE)
```

## Arguments

- palette:

  Name of the palette. One of `"main"`, `"core"`, `"diverging"`,
  `"seq_blue"`, `"seq_sage"`, `"seq_coral"`. Defaults to `"main"`.

- n:

  Number of colours to return. If `NULL` (default) all colours in the
  palette are returned. If `n` is smaller than the palette length the
  first `n` colours are returned; if larger an error is thrown.

- reverse:

  Logical. Reverse the order of the palette? Default `FALSE`.

## Value

A named (or unnamed for gradient palettes) character vector of hex
colour codes.

## Examples

``` r
rtractor_palette()
#>      coral      cream       sage steel_blue        ink 
#>  "#FFB6A6"  "#FFEBD3"  "#9BCEC1"  "#67A2C5"  "#23475C" 
rtractor_palette("core")
#>      coral      cream       sage steel_blue 
#>  "#FFB6A6"  "#FFEBD3"  "#9BCEC1"  "#67A2C5" 
rtractor_palette("seq_blue", n = 3)
#> [1] "#D9E8F0" "#9DC1D6" "#67A2C5"
rtractor_palette("diverging", reverse = TRUE)
#> [1] "#67A2C5" "#9BCEC1" "#CCE3D3" "#FFEBD3" "#FFD3BE" "#FFB6A6"
```

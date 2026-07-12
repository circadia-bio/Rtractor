# Rtractor continuous colour scale

Interpolates across an Rtractor palette for continuous `colour` data.
The `"diverging"` palette is recommended for signed quantities (e.g.
Lyapunov exponent sign, spectrum asymmetry); the `seq_*` palettes for
unipolar data (e.g. entropy magnitude).

## Usage

``` r
scale_colour_rtractor_c(palette = "seq_blue", reverse = FALSE, ...)

scale_color_rtractor_c(palette = "seq_blue", reverse = FALSE, ...)
```

## Arguments

- palette:

  Palette name. Default `"seq_blue"`.

- reverse:

  Logical. Reverse the palette? Default `FALSE`.

- ...:

  Additional arguments passed to
  [`ggplot2::scale_colour_gradientn()`](https://ggplot2.tidyverse.org/reference/scale_gradient.html).

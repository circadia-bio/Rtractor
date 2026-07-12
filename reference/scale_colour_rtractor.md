# Rtractor discrete colour scale

Applies an Rtractor qualitative palette to the `colour` aesthetic.
Defaults to `"core"` (excludes the dark `ink` colour, which is reserved
for text/axis annotation rather than data encoding).

## Usage

``` r
scale_colour_rtractor(palette = "core", reverse = FALSE, ...)

scale_color_rtractor(palette = "core", reverse = FALSE, ...)
```

## Arguments

- palette:

  Palette name. Default `"core"`.

- reverse:

  Logical. Reverse the palette? Default `FALSE`.

- ...:

  Additional arguments passed to
  [`ggplot2::discrete_scale()`](https://ggplot2.tidyverse.org/reference/discrete_scale.html).

## Examples

``` r
library(ggplot2)
ggplot(mtcars, aes(wt, mpg, colour = factor(cyl))) +
  geom_point(size = 3) +
  scale_colour_rtractor() +
  theme_rtractor()

```

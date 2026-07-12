# Rtractor discrete fill scale

Applies an Rtractor qualitative palette to the `fill` aesthetic.

## Usage

``` r
scale_fill_rtractor(palette = "core", reverse = FALSE, ...)
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
ggplot(mpg, aes(class, fill = drv)) +
  geom_bar() +
  scale_fill_rtractor() +
  theme_rtractor()

```

# Rtractor ggplot2 theme

A clean ggplot2 theme built on
[`ggplot2::theme_minimal()`](https://ggplot2.tidyverse.org/reference/ggtheme.html),
visually distinct from `theme_circadia()` (softer palette, same
structural conventions) so that Rtractor figures are recognisable at a
glance.

## Usage

``` r
theme_rtractor(
  base_size = 14,
  base_family = "",
  grid = "none",
  legend_position = "right"
)
```

## Arguments

- base_size:

  Base font size in points. Default `14`.

- base_family:

  Base font family. Default `""` (ggplot2 default).

- grid:

  Which grid lines to show. One of `"none"` (default), `"y"` (horizontal
  only), `"xy"` (both), `"x"` (vertical only).

- legend_position:

  Position of the legend passed to
  [`ggplot2::theme()`](https://ggplot2.tidyverse.org/reference/theme.html).
  Default `"right"`.

## Value

A
[`ggplot2::theme()`](https://ggplot2.tidyverse.org/reference/theme.html)
object.

## Examples

``` r
library(ggplot2)
ggplot(mtcars, aes(wt, mpg, colour = factor(cyl))) +
  geom_point(size = 2) +
  scale_colour_rtractor() +
  theme_rtractor()

```

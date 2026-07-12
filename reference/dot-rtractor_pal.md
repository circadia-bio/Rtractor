# Build a palette function from an Rtractor palette

Returns a function that interpolates `n` colours from the named palette.
Used internally by the scale constructors.

## Usage

``` r
.rtractor_pal(palette = "core", reverse = FALSE)
```

## Arguments

- palette:

  Palette name passed to
  [`rtractor_palette()`](https://rtractor.circadia-lab.uk/reference/rtractor_palette.md).

- reverse:

  Logical. Reverse the palette?

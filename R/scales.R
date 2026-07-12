# Rtractor ggplot2 scales
#
# ggplot2 is a Suggests-only dependency (see theme.R note). Every function
# here guards with requireNamespace() and calls ggplot2:: explicitly rather
# than using @importFrom, so no hard Imports dependency is created.

# ---- Colour interpolation helper -------------------------------------------

#' Build a palette function from an Rtractor palette
#'
#' Returns a function that interpolates `n` colours from the named palette.
#' Used internally by the scale constructors.
#'
#' @param palette Palette name passed to [rtractor_palette()].
#' @param reverse Logical. Reverse the palette?
#' @keywords internal
.rtractor_pal <- function(palette = "core", reverse = FALSE) {
  pal <- rtractor_palette(palette, reverse = reverse)
  function(n) {
    if (n <= length(pal)) {
      unname(pal[seq_len(n)])
    } else {
      grDevices::colorRampPalette(pal)(n)
    }
  }
}

# ---- Discrete scales --------------------------------------------------------

#' Rtractor discrete colour scale
#'
#' Applies an Rtractor qualitative palette to the `colour` aesthetic.
#' Defaults to `"core"` (excludes the dark `ink` colour, which is reserved
#' for text/axis annotation rather than data encoding).
#'
#' @param palette Palette name. Default `"core"`.
#' @param reverse Logical. Reverse the palette? Default `FALSE`.
#' @param ... Additional arguments passed to [ggplot2::discrete_scale()].
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg, colour = factor(cyl))) +
#'   geom_point(size = 3) +
#'   scale_colour_rtractor() +
#'   theme_rtractor()
#'
#' @export
scale_colour_rtractor <- function(palette = "core", reverse = FALSE, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for scale_colour_rtractor().", call. = FALSE)
  }
  ggplot2::discrete_scale(
    aesthetics = "colour",
    palette    = .rtractor_pal(palette, reverse),
    ...
  )
}

#' @rdname scale_colour_rtractor
#' @export
scale_color_rtractor <- scale_colour_rtractor

#' Rtractor discrete fill scale
#'
#' Applies an Rtractor qualitative palette to the `fill` aesthetic.
#'
#' @inheritParams scale_colour_rtractor
#'
#' @examples
#' library(ggplot2)
#' ggplot(mpg, aes(class, fill = drv)) +
#'   geom_bar() +
#'   scale_fill_rtractor() +
#'   theme_rtractor()
#'
#' @export
scale_fill_rtractor <- function(palette = "core", reverse = FALSE, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for scale_fill_rtractor().", call. = FALSE)
  }
  ggplot2::discrete_scale(
    aesthetics = "fill",
    palette    = .rtractor_pal(palette, reverse),
    ...
  )
}

# ---- Continuous scales -------------------------------------------------------

#' Rtractor continuous colour scale
#'
#' Interpolates across an Rtractor palette for continuous `colour` data.
#' The `"diverging"` palette is recommended for signed quantities (e.g.
#' Lyapunov exponent sign, spectrum asymmetry); the `seq_*` palettes for
#' unipolar data (e.g. entropy magnitude).
#'
#' @param palette Palette name. Default `"seq_blue"`.
#' @param reverse Logical. Reverse the palette? Default `FALSE`.
#' @param ... Additional arguments passed to [ggplot2::scale_colour_gradientn()].
#'
#' @export
scale_colour_rtractor_c <- function(palette = "seq_blue", reverse = FALSE, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for scale_colour_rtractor_c().", call. = FALSE)
  }
  colours <- unname(rtractor_palette(palette, reverse = reverse))
  ggplot2::scale_colour_gradientn(colours = colours, ...)
}

#' @rdname scale_colour_rtractor_c
#' @export
scale_color_rtractor_c <- scale_colour_rtractor_c

#' Rtractor continuous fill scale
#'
#' Interpolates across an Rtractor palette for continuous `fill` data.
#'
#' @inheritParams scale_colour_rtractor_c
#'
#' @export
scale_fill_rtractor_c <- function(palette = "seq_blue", reverse = FALSE, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for scale_fill_rtractor_c().", call. = FALSE)
  }
  colours <- unname(rtractor_palette(palette, reverse = reverse))
  ggplot2::scale_fill_gradientn(colours = colours, ...)
}

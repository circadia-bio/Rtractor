# Rtractor colour palette
#
# Palette definition and helper functions for retrieving colours.
# Deliberately independent from the `circadia` package (Rtractor should run
# standalone, per the ecosystem's isolation principle) even though the
# function names mirror it for consistency across Circadia Lab / CoDe-Neuro
# Lab packages.

#' Rtractor colour palette
#'
#' Named list of all palettes available in the Rtractor package.
#'
#' @format A named list of character vectors, each containing hex colour codes.
#' @keywords internal
.rtractor_palettes <- list(

  # ---- Qualitative ---------------------------------------------------
  # Four brand colours plus a dark "ink" colour for text/axes/annotation.
  # ink is a darkened derivative of steel_blue, kept in-family rather than
  # a generic near-black.
  main = c(
    coral      = "#FFB6A6",
    cream      = "#FFEBD3",
    sage       = "#9BCEC1",
    steel_blue = "#67A2C5",
    ink        = "#23475C"
  ),

  # Four-colour subset without ink (for fills/points where a dark anchor
  # isn't wanted)
  core = c(
    coral      = "#FFB6A6",
    cream      = "#FFEBD3",
    sage       = "#9BCEC1",
    steel_blue = "#67A2C5"
  ),

  # ---- Diverging -------------------------------------------------------
  # warm (coral) -> neutral (cream) -> cool (sage/blue).
  # Maps naturally onto signed quantities: Lyapunov exponent sign,
  # multifractal spectrum asymmetry, entropy deviation from a null model.
  diverging = c(
    "#FFB6A6", "#FFD3BE", "#FFEBD3", "#CCE3D3", "#9BCEC1", "#67A2C5"
  ),

  # ---- Sequential --------------------------------------------------------
  # Monochromatic ramps of each core colour, light -> saturated/dark.
  seq_blue = c(
    "#D9E8F0", "#9DC1D6", "#67A2C5", "#3D6E8C", "#23475C"
  ),

  seq_sage = c(
    "#E3F2EC", "#BEE0D5", "#9BCEC1", "#5FA592", "#3D6E60"
  ),

  seq_coral = c(
    "#FFEAE3", "#FFD0C1", "#FFB6A6", "#F58A72", "#D65F45"
  )
)

# ---- Public functions ------------------------------------------------------

#' Retrieve an Rtractor palette
#'
#' Returns a named character vector of hex colour codes for the requested
#' palette. Suitable for direct use in [ggplot2::scale_fill_manual()] or
#' [ggplot2::scale_colour_manual()].
#'
#' @param palette Name of the palette. One of `"main"`, `"core"`,
#'   `"diverging"`, `"seq_blue"`, `"seq_sage"`, `"seq_coral"`.
#'   Defaults to `"main"`.
#' @param n Number of colours to return. If `NULL` (default) all colours in
#'   the palette are returned. If `n` is smaller than the palette length the
#'   first `n` colours are returned; if larger an error is thrown.
#' @param reverse Logical. Reverse the order of the palette? Default `FALSE`.
#'
#' @return A named (or unnamed for gradient palettes) character vector of
#'   hex colour codes.
#'
#' @examples
#' rtractor_palette()
#' rtractor_palette("core")
#' rtractor_palette("seq_blue", n = 3)
#' rtractor_palette("diverging", reverse = TRUE)
#'
#' @export
rtractor_palette <- function(palette = "main", n = NULL, reverse = FALSE) {
  pal <- .rtractor_palettes[[palette]]
  if (is.null(pal)) {
    stop(
      sprintf(
        "Unknown palette '%s'. Choose from: %s.",
        palette,
        paste(names(.rtractor_palettes), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!is.null(n)) {
    if (n > length(pal)) {
      stop(
        sprintf(
          "Palette '%s' has only %d colours; %d requested.",
          palette, length(pal), n
        ),
        call. = FALSE
      )
    }
    pal <- pal[seq_len(n)]
  }

  if (reverse) pal <- rev(pal)
  pal
}

#' List all available Rtractor palettes
#'
#' Prints the names and sizes of all palettes defined in the package.
#'
#' @return A named integer vector of palette lengths, invisibly.
#'
#' @examples
#' rtractor_palettes()
#'
#' @export
rtractor_palettes <- function() {
  sizes <- lengths(.rtractor_palettes)
  cat("Available Rtractor palettes:\n")
  cat("  Qualitative:\n")
  for (nm in c("main", "core")) {
    cat(sprintf("    %-20s  %d colours\n", nm, sizes[[nm]]))
  }
  cat("  Diverging:\n")
  cat(sprintf("    %-20s  %d colours\n", "diverging", sizes[["diverging"]]))
  cat("  Sequential:\n")
  for (nm in c("seq_blue", "seq_sage", "seq_coral")) {
    cat(sprintf("    %-20s  %d colours\n", nm, sizes[[nm]]))
  }
  invisible(sizes)
}

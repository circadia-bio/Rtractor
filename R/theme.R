# Rtractor ggplot2 theme
#
# ggplot2 is a Suggests-only dependency (per the isolation principle —
# Rtractor's core metric functions don't need it). Every function here
# guards with requireNamespace() and calls ggplot2:: explicitly rather than
# using @importFrom, so no hard Imports dependency is created.

#' Rtractor ggplot2 theme
#'
#' A clean ggplot2 theme built on [ggplot2::theme_minimal()], visually
#' distinct from `theme_circadia()` (softer palette, same structural
#' conventions) so that Rtractor figures are recognisable at a glance.
#'
#' @param base_size Base font size in points. Default `14`.
#' @param base_family Base font family. Default `""` (ggplot2 default).
#' @param grid Which grid lines to show. One of `"none"` (default),
#'   `"y"` (horizontal only), `"xy"` (both), `"x"` (vertical only).
#' @param legend_position Position of the legend passed to
#'   [ggplot2::theme()]. Default `"right"`.
#'
#' @return A [ggplot2::theme()] object.
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg, colour = factor(cyl))) +
#'   geom_point(size = 2) +
#'   scale_colour_rtractor() +
#'   theme_rtractor()
#'
#' @export
theme_rtractor <- function(
    base_size       = 14,
    base_family     = "",
    grid            = "none",
    legend_position = "right"
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for theme_rtractor().", call. = FALSE)
  }

  grid <- match.arg(grid, c("y", "xy", "x", "none"))

  show_x <- grid %in% c("xy", "x")
  show_y <- grid %in% c("xy", "y")

  ink       <- "#23475C"
  grid_line <- ggplot2::element_line(colour = "#E0E0E0", linewidth = 0.4)
  no_line   <- ggplot2::element_blank()
  axis_line <- ggplot2::element_line(colour = ink, linewidth = 0.5)

  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      # --- Text ---
      plot.title    = ggplot2::element_text(
        colour = ink, face = "bold",
        size   = ggplot2::rel(1.2), margin = ggplot2::margin(b = 8)
      ),
      plot.subtitle = ggplot2::element_text(
        colour = "#555555",
        size   = ggplot2::rel(1.0), margin = ggplot2::margin(b = 10)
      ),
      plot.caption  = ggplot2::element_text(
        colour = "#888888", size = ggplot2::rel(0.8),
        margin = ggplot2::margin(t = 8), hjust = 1
      ),
      axis.title    = ggplot2::element_text(colour = ink, size = ggplot2::rel(1.0)),
      axis.text     = ggplot2::element_text(colour = "#444444", size = ggplot2::rel(0.95)),
      strip.text    = ggplot2::element_text(
        colour = ink, face = "bold", size = ggplot2::rel(1.0)
      ),
      legend.title  = ggplot2::element_text(colour = ink, size = ggplot2::rel(0.95)),
      legend.text   = ggplot2::element_text(colour = "#555555", size = ggplot2::rel(0.9)),

      # --- Panel ---
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
      panel.border     = ggplot2::element_blank(),

      # --- Axis lines (left and bottom only) ---
      axis.line.x.bottom = axis_line,
      axis.line.y.left   = axis_line,
      axis.line.x.top    = no_line,
      axis.line.y.right  = no_line,

      # --- Grid ---
      panel.grid.major.x = if (show_x) grid_line else no_line,
      panel.grid.major.y = if (show_y) grid_line else no_line,
      panel.grid.minor   = no_line,

      # --- Axes ---
      axis.ticks        = ggplot2::element_line(colour = "#AAAAAA", linewidth = 0.4),
      axis.ticks.length = ggplot2::unit(4, "pt"),

      # --- Legend ---
      legend.position   = legend_position,
      legend.key        = ggplot2::element_rect(fill = NA, colour = NA),
      legend.background = ggplot2::element_rect(fill = NA, colour = NA),

      # --- Facets ---
      strip.background = ggplot2::element_rect(fill = "#F0F4F0", colour = NA),

      # --- Spacing ---
      plot.margin = ggplot2::margin(12, 16, 10, 12)
    )
}

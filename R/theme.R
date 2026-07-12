# Rtractor ggplot2 theme

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
#' @importFrom ggplot2 theme_minimal theme element_text element_line
#'   element_blank element_rect margin rel unit
#' @export
theme_rtractor <- function(
    base_size       = 14,
    base_family     = "",
    grid            = "none",
    legend_position = "right"
) {
  grid <- match.arg(grid, c("y", "xy", "x", "none"))

  show_x <- grid %in% c("xy", "x")
  show_y <- grid %in% c("xy", "y")

  ink       <- "#23475C"
  grid_line <- element_line(colour = "#E0E0E0", linewidth = 0.4)
  no_line   <- element_blank()
  axis_line <- element_line(colour = ink, linewidth = 0.5)

  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      # --- Text ---
      plot.title    = element_text(
        colour = ink, face = "bold",
        size   = rel(1.2), margin = margin(b = 8)
      ),
      plot.subtitle = element_text(
        colour = "#555555",
        size   = rel(1.0), margin = margin(b = 10)
      ),
      plot.caption  = element_text(
        colour = "#888888", size = rel(0.8),
        margin = margin(t = 8), hjust = 1
      ),
      axis.title    = element_text(colour = ink, size = rel(1.0)),
      axis.text     = element_text(colour = "#444444", size = rel(0.95)),
      strip.text    = element_text(
        colour = ink, face = "bold", size = rel(1.0)
      ),
      legend.title  = element_text(colour = ink, size = rel(0.95)),
      legend.text   = element_text(colour = "#555555", size = rel(0.9)),

      # --- Panel ---
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.border     = element_blank(),

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
      axis.ticks        = element_line(colour = "#AAAAAA", linewidth = 0.4),
      axis.ticks.length = unit(4, "pt"),

      # --- Legend ---
      legend.position   = legend_position,
      legend.key        = element_rect(fill = NA, colour = NA),
      legend.background = element_rect(fill = NA, colour = NA),

      # --- Facets ---
      strip.background = element_rect(fill = "#F0F4F0", colour = NA),

      # --- Spacing ---
      plot.margin = margin(12, 16, 10, 12)
    )
}

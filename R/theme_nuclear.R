#' Nuclear Theme
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License version 3 as published by
#the Free Software Foundation.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#'
#' @inheritParams ggplot2::theme_bw
#'
#' @family themes
#' @export
#' @importFrom ggplot2 theme_bw
theme_nuclear  <- function () {
  theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(colour = "white", fill = "white"),
          plot.background = element_rect(colour = "white", fill = "white"),
          panel.background = element_rect(colour = "white", fill = "white"),
          legend.key = element_rect(colour = "white", fill = "white"),
          axis.title = element_text(size=22),
          axis.text  = element_text(size=18),
          axis.ticks = element_line(size = 0.45),
          axis.text.y = element_text(size = 20,
                         margin = unit(c(t = 0, r = 5, b = 0, l = 0), "mm")),
          axis.line = element_line(size = 0.5, linetype = "solid"))
}

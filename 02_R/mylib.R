theme_Publication <- function(base_size = 14, base_family = "helvetica") {
  library(ggthemes)
  library(ggplot2)
  ggthemes::theme_foundation(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      legend.background = element_rect(fill = "white", colour = NA),
      panel.borde = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(),
      axis.title = element_text(face = "bold", size = rel(1)),
      axis.title.y = element_text(angle = 90, vjust = 2),
      axis.title.x = element_text(vjust = -0.2),
      panel.grid.major = element_line(colour = "#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(fill = "white", colour = NA),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.size = unit(0.2, "cm"),
      legend.margin = margin(t = 0, unit = "cm"),
      legend.title = element_text(face = "italic"),
      strip.background = element_rect(fill = "white", colour = NA),
      strip.text = element_text(face = "bold"),
      plot.margin = unit(c(10, 5, 5, 5), "mm")
    )
}

export_plot <- function(myplot, file, width = 8, height = 6,
                        units = "in", type = "png", bg = "white") {
  require(Cairo)
  Cairo::Cairo(width = width, height = height, file = file,
               type = type, bg = bg, units = units, dpi = 300)
  plot(myplot)
  dev.off()
  message("Exported: ", file)
  invisible(myplot)
}
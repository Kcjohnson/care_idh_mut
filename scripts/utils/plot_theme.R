# Standard stripped down plotting theme.
plot_theme    <- theme(text = element_text(size = 8), 
                       axis.text = element_text(size = 8),
                       axis.title = element_text(size = 8),
                       strip.text = element_text(size = 8),
                       panel.background = element_rect(fill = "transparent"),
                       axis.line = element_blank(),
                       strip.background = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), 
                       panel.border = element_blank(),
                       axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                       axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))

# Variant theme
theme_pvsr <- function(theme = theme_pubr(), 
                       axis.title = element_text(size = 13.5), 
                       axis.text.x = element_text(size = 13.5), 
                       axis.text.y = element_text(size = 13.5), 
                       axis.ticks.x = element_line(), 
                       axis.ticks.y = element_line(), 
                       legend.title = element_text(size = 13.5), 
                       legend.text = element_text(size = 13.5), 
                       legend.position = "top", 
                       plot.title = element_text(size = 13.5), 
                       strip.text = element_text(size = 13.5), 
                       panel.spacing = unit(.0, "lines"),  
                       panel.border = element_rect(color = "black", fill = NA), 
                       strip.background = element_rect(color = "black", fill = NA)) {
  
  res <- theme +
    theme(axis.title = axis.title,
          axis.text.x = axis.text.x,
          axis.text.y = axis.text.y,
          axis.ticks.x = axis.ticks.x, axis.ticks.y = axis.ticks.y,
          legend.title = legend.title, legend.text = legend.text, legend.position = legend.position,
          plot.title = plot.title,
          strip.text = strip.text) +
    theme(panel.spacing = panel.spacing,
          panel.border = panel.border,
          strip.background = strip.background)
  
  return(res)
}
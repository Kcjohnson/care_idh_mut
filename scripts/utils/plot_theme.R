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

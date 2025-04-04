# """
# Project: track_genome
# April 2025 - XQ - LUMC
# visualization of genomic regions in R
# """


library(ggplot2)

# --------------------------------------------------
# Function to plot BigWig data with optional smoothing and customizable features
# span: controls how much of the data is used for smoothing
# level: interval for confidence interval
## Wide Interval: More uncertainty, broader shading (e.g., sparse data, high variability, higher level).
## Narrow Interval: Less uncertainty, tighter shading (e.g., dense data, low variability, lower level).


# --------------------------------------------------
plot_bigwig <- function(data, 
                        figure_dir, 
                        file_name, 
                        colour_set,
                        title = "chr:start-end",
                        feature1 = "day", 
                        feature2 = "week", 
                        smooth = TRUE, 
                        vline_values = NULL,    # New parameter for vertical lines
                        span = 0.1,             # Default span
                        level = 0.95,           # Default confidence level
                        se = FALSE) {           # Default for confidence interval display
  
  # Check if data is provided
  if (is.null(data)) stop("Providing data is required.")
  if (!all(c(feature1, feature2) %in% names(data))) {
  stop("Both 'feature1' and 'feature2' must be present in the data.")
  }
  
  # Message: Starting the plot creation process
  message("Starting the plot creation process.")
  
  # Start the ggplot object
  p <- ggplot(data, aes(x = start, y = score))
  
  p <- p + geom_point(aes(color = .data[[feature1]]), 
                      size = 1,
                      alpha = 0.15)

  # Add geom based on smoothing choice
if (smooth) {
    message("Adding smoothed lines with confidence intervals.")
    p <- p + geom_smooth(aes(color = .data[[feature1]], 
                            fill = .data[[feature1]]), 
                         method = "loess", span = span, 
                         size = 1.5, 
                         level = level, 
                         alpha = 0.75, 
                         se = se)
  }
# Message: Adding labels, scales, and faceting
  message("Adding labels, color scales, and facet wrapping.")
  
  # Add labels, scales, faceting, and theme
  p <- p +
    labs(x = "Genomic Position", y = "Coverage", 
         title = title) +
    scale_color_manual(values = setNames(data[[colour_set]], data[[feature1]])) +
    scale_fill_manual(values = setNames(data[[colour_set]], data[[feature1]])) +
    coord_cartesian(ylim = c(0, ceiling(max(data$score, na.rm = TRUE) / 500) * 500)) +
    facet_wrap(as.formula(paste("~", feature2)), ncol = 1) +
    geom_vline(xintercept = vline_values, linetype = "dashed", color = "red") +
    theme_minimal()
  # Customize legend based on se
  # Customize legend based on smooth and se, with alpha = 1
  if (!smooth) {
      p <- p + guides(fill = "none", 
                      color = guide_legend(override.aes = list(linetype = 1, 
                      size = 3, 
                      alpha = 1)))
    }

  # Message: Saving the plot
  message("Saving the plot to the specified directory.")
  
  # Save the plot
  ggsave(paste0(figure_dir, file_name), 
         plot = p, bg = 'white', height = 10, width = 16.2, dpi = 320)
  
  # Message: Plotting completed
  message("Plotting completed. Plot saved.")
  
  # Return the plot object
  return(p)
}

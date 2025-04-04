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

plot_bigwig <- function(data, 
                        figure_dir, 
                        file_name, 
                        colour_set,
                        title = "chr:start-end",
                        feature1 = "day", 
                        feature2 = "week", 
                        do_smooth = TRUE) {
  
  # Check if data is provided
  if (is.null(data)) stop("Providing data is required.")
  
  # Message: Starting the plot creation process
  message("Starting the plot creation process.")
  
  # Start the ggplot object
  p <- ggplot()
  
  # Add geom based on smoothing choice
  if (do_smooth) {
    message("Adding smoothed lines with confidence intervals.")
    p <- p + geom_smooth(data = data, 
                         aes(x = start, y = score, 
                             color = .data[[feature1]], 
                             fill = .data[[feature1]]),  # Dynamic color/fill
                         method = "loess", span = 0.1, 
                         size = 1.5,  # Thicker main line
                         level = 0.95,  # 95% confidence interval
                         alpha = 0.75,  # Light shading for visibility
                         se = FALSE)  # Include confidence interval
  } else {
    message("Adding points with smoothed lines (without confidence intervals).")
    p <- p + geom_point(data = data, 
                        aes(x = start, y = score, 
                            color = .data[[feature1]]),  # Dynamic color for points
                        size = 2,
                        alpha = 0.15) + # Point size
            geom_smooth(data = data, 
                         aes(x = start, y = score, 
                             color = .data[[feature1]], 
                             fill = .data[[feature1]]),  # Dynamic color/fill
                         method = "loess", span = 0.1, 
                         size = 1.5,  # Thicker main line
                         level = 0.95,  # 95% confidence interval
                         alpha = 0.75,  # Light shading for visibility
                         se = FALSE)  # Include confidence interval
  }
  
  # Message: Adding labels, scales, and faceting
  message("Adding labels, color scales, and facet wrapping.")
  
  # Add labels, scales, faceting, and theme
  p <- p +
    labs(x = "Genomic Position", y = "Coverage", 
         title = title) +
    scale_color_manual(values = setNames(data[[colour_set]], data[[feature1]])) +  # Dynamic color mapping
    scale_fill_manual(values = setNames(data[[colour_set]], data[[feature1]])) +   # Dynamic fill mapping
    coord_cartesian(ylim = c(0, max(data$score))) +  # Adjust y-axis
    facet_wrap(as.formula(paste("~", feature2)), ncol = 1) +  # Dynamic faceting
    theme_minimal()
  
  # Message: Saving the plot
  message("Saving the plot to the specified directory.")
  
  # Save the plot
  ggsave(paste0(figure_dir, file_name), 
         plot = p, bg = 'white', height = 10, width = 20, dpi = 320)
  
  # Message: Plotting completed
  message("Plotting completed. Plot saved.")
  
  # Return the plot object for previewing if desired
  return(p)
}

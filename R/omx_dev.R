library(ggplot2)


# dev: sp--------------------------------------------------
# Function to quickly save a ggplot object for checking
sp <- function(plot,
               figure_dir = "./figure/testing/",   
               filename = 'test', 
               width = 10, height = 5, dpi = 300, bg = "white") {
  # Check if the directory exists, if not, create it
  if (!dir.exists(figure_dir)) {
    dir.create(figure_dir, recursive = TRUE)
  }
  
  # Generate the full file path
  file_path <- file.path(figure_dir, paste0(filename, ".png"))
  
  # Check if the file exists and increment the filename if necessary
  if (file.exists(file_path)) {
    counter <- 1
    repeat {
      new_filename <- paste0(filename, "_", counter)
      file_path <- file.path(figure_dir, paste0(new_filename, ".png"))
      if (!file.exists(file_path)) {
        filename <- new_filename
        break
      }
      counter <- counter + 1
    }
  }
  
  # Save the plot
  ggsave(filename = file_path, 
         plot = plot, bg = bg, 
         width = width, height = height, dpi = dpi)
  
  # Print a message indicating where the plot was saved
  message(paste("Plot saved to", file_path))
}

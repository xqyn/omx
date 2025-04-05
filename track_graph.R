# """
# Project: track_genome
# April 2025 - XQ - LUMC
# visualization of genomic regions in R
# """


library(ggplot2)


# """
# bigWig file processing and plotting
# """

# --------------------------------------------------
# Function to plot BigWig data with optional smoothing and customizable features
# span: controls how much of the data is used for smoothing
# level: interval for confidence interval
## Wide Interval: More uncertainty, broader shading (e.g., sparse data, high variability, higher level).
## Narrow Interval: Less uncertainty, tighter shading (e.g., dense data, low variability, lower level).


# plot_bigwig --------------------------------------------------
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


# """
# BED file processing and isoform mapping
# """

# spacing_yline --------------------------------------------------
spacing_yline <- function(df) {
    # --------------------------------------------------
    #' Add spacing between genomic regions based on the number of rows
    # --------------------------------------------------
    n <- nrow(df)

    if (n >= 1 && n <= 100) {
    # determine the spacing tier (each 5 rows = 1 tier)
    tier <- ceiling(n / 5)

    # calculate spacing: starts at 0.5, increases by 0.25 per tier
    spacing_value <- 0.5 + 0.25 * (tier - 1)

    # cap the maximum spacing to 5
    spacing <- min(spacing_value, 5)
    } else {
    spacing <- c(NA, NA)
    message("Reduce number of genes in the region")
    }
    # assign y_line values based on the number of rows and spacing
    df$y_line <- 1:nrow(df) * spacing

    return(df)
}

# isoform_mapping --------------------------------------------------
isoform_mapping <- function(isoform, 
                            arrow_dist = 300, 
                            base_offset = 50) {
    #--------------------------------------------------
    #' Isoform block and arrow processing
    #'
    #' Given a single isoform row (from a BED-like dataframe),
    #' returns a list containing exon block annotations and
    #' arrow symbols representing exon connections.
    #'
    #' @param isoform single-row dataframe containing isoform (gene in BED file) info. 
    #'        Required columns: start, end, strand, blocks, y_line
    #'        y_line: coordinates line for plotting each isoform
    #'        blocks: list of blocks (exons) for the isoform
    #' @param arrow_dist bases distance between arrows (default: 300)
    #' @param base_offset extra bases spacing for arrow positioning for prevent arrow overlapping exon (default: 5) 
    #'
    #' @return A list with:
    #'   - block: dataframe with exon coordinates and arrow positions
    #'   - arrows: dataframe with arrow symbols and positions
    #--------------------------------------------------

    required_cols <- c("start", "end", "strand", "blocks", "y_line")
    missing_cols <- setdiff(required_cols, names(isoform))
    if (length(missing_cols) > 0) {
      stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
    }
    
    # --- Extract block information
    # block: one block is one exon
    block <- as.data.frame(isoform$blocks[[1]])
    # --- add strand and y_line to block
    block$strand <- isoform$strand
    block$y_line <- isoform$y_line
    
    # --- exon postions: |[=====]_i_th----[==]----[======]----[=====]|
    # exon start = isoform start + block start: | + [_i_th
    # exon end = isoform start + block end: | + ]_i_th
    block$exon_start <- isoform$start + block$start
    block$exon_end   <- isoform$start + block$end

    # --- arrow positions (if not exon -> arrow) to block
    block$arrow_start <- NA
    block$arrow_end   <- NA
    if (nrow(block) > 1) {
      # calculate arrow start and end positions for each block
      block$arrow_start[1:(nrow(block)-1)] <- block$exon_end[1:(nrow(block)-1)] + base_offset
      block$arrow_end[1:(nrow(block)-1)]   <- block$exon_start[2:nrow(block)] - base_offset
      
      # add final arrow position before last exon (so that it doesn't overlap)
      #block$arrow_start[nrow(block)]       <- block$exon_start[nrow(block)] - arrow_dist
      #block$arrow_end[nrow(block)]         <- block$exon_start[nrow(block)] - arrow_dist
      ## note: shall I added the last arrow position or keep it NA?
    } else {
      # if only one block (one exon), set arrow positions to NA
      block$arrow_start <- NA
      block$arrow_end   <- NA
    }

    # generate arrow positions
    arrow_position <- unlist(lapply(seq_len(nrow(block)), function(j) {
      # check if arrow_start and arrow_end are NA, then skip
      if (is.na(block$arrow_start[j]) || is.na(block$arrow_end[j])) return(NULL)
      
      # else, calculate the sequence of positions
      start <- block$arrow_start[j]
      end   <- block$arrow_end[j]
      by_val <- if (start < end) arrow_dist else -arrow_dist
      seq(from = start, to = end, by = by_val)
    }))

    # Create arrow data frame
    arrow_df <- if (length(arrow_position) > 0) {
      data.frame(
        arrow = rep(ifelse(isoform$strand == "+", ">", "<"), length(arrow_position)),
        arrow_pos = arrow_position,
        y_line = rep(isoform$y_line, length(arrow_position))
      )
    } else {
      NULL
    }
    
    return(list(block = block, arrows = arrow_df))
}

# plot_isoform_map --------------------------------------------------
plot_isoform_map <- function(bed, 
                            block_df, 
                            arrow_df, region_pad = 2000) {
    #--------------------------------------------------
    #' plot Isoform containing exons and arrows
    #'
    #' visualizes exon structure of isoforms using block (exon) and arrow
    #'
    #' @param bed dataframe containing isoform information including 
    #'        start, end, y_line, name, and strand
    #' @param block_df dataframe of exon block positions (from isoform_mapping function)
    #' @param arrow_df dataframe of arrow annotations (from isoform_mapping function)
    #' @param region_pad base padding added to left and right of plotting region (default = 2000)
    #'
    #' @return a ggplot object visualizing the isoform structure
    #--------------------------------------------------
    p <- ggplot() +
      # horizontal base line per gene
      geom_segment(data = bed,
                    aes(x = start, xend = end,
                        y = y_line, yend = y_line),
                    size = 0.15, alpha = 0.5,
                    color = 'black') +

      # block boxed for exons
      geom_segment(data = block_df, 
                    aes(x = exon_start, xend = exon_end, 
                        y = y_line, yend = y_line, 
                        color = strand),  
                    size = 4) +  
      scale_color_manual(values = c("+" = "#017075", "-" = "#F48D79")) +  # strand color

      # arrows between exons
      geom_text(data = arrow_df, 
                aes(x = arrow_pos, y = y_line, label = arrow),
                size = 1.5, color = "black", vjust = 0.3) + 

      # add gene names
      geom_text(data = bed, aes(x = start , y = y_line, label = name),
                hjust = ifelse(bed$strand == "+", 1.15, 1.15),
                vjust = 0.5, size = 2, color = "#0040B0") + 

      # axis coordinates
      xlim(min(bed$start) - region_pad, max(bed$end) + region_pad) + 
      ylim(-1, max(bed$y_line) + 2) + 

      # remove plot labels
      labs(title = NULL, x = NULL, y = NULL) +  
      theme_classic()

    return(p)
}

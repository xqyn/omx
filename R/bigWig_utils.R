# bigWig_utils.R
# project: omx
# april 2025 - XQ - LUMC
# modifying genomic regions format


# process_bigwig --------------------------------------------------
#' Process BigWig files into a merged data frame
#'
#' This function imports BigWig files as GRanges objects for a specified genomic region with
#' flanking regions, converts them to a data frame, and optionally merges them with a mapping
#' data frame using a specified column for merging.
#'
#' @param bw_files Character vector of paths to BigWig files to be processed.
#' @param region_str Character string specifying the genomic region in "chr:start-end" format
#'                   (default: "X:73851081-73851581").
#' @param flank_left Integer specifying the left flanking region size in base pairs (default: 2000).
#' @param flank_right Integer specifying the right flanking region size in base pairs
#'                    (default: \code{NULL}, which sets it equal to \code{flank_left}).
#' @param mapping Optional data frame containing mapping information to merge with the processed data.
#' @param by Character string specifying the column name to use for merging with the \code{mapping}
#'           data frame (default: "code").
#'
#' @return A data frame containing the combined data from all BigWig files, optionally merged
#'         with mapping information if provided.
#'
#' @details This function processes BigWig files by importing them into \code{GRanges} objects
#' within a specified genomic region (adjusted by flanking regions), converts them into a single
#' data frame with a source identifier, and optionally merges with a provided mapping data frame
#' using the specified \code{by} column. It dynamically sources the \code{chr_region} function
#' from a separate script if not already available.
#'
#' @examples
#' \dontrun{
#' bw_files <- c("path/to/file1.bw", "path/to/file2.bw")
#' mapping <- data.frame(sample_id = c("file1", "file2"), condition = c("A", "B"))
#' df <- process_bigwig(bw_files, "chr1:1000-2000", flank_left = 500, 
#'                      mapping = mapping, by = "sample_id")
#' }
#'
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom tools file_path_sans_ext
#' @export

process_bigwig <- function(bw_files,
                           region_str = "X:73851081-73851581",
                           flank_left = 0,
                           flank_right = NULL,
                           mapping = NULL,
                           by = "code") {
    
    # Create genomic region with flanks
    message("Creating genomic region with flanks...")
    if (is.null(region_str)) stop("Please provide a region_str in format 'chr:start-end'.")
    if (is.null(flank_right)) flank_right <- flank_left

    # split region string and parse components
    region_parts <- strsplit(region_str, "[:-]")[[1]]
    chr <- region_parts[1]
    start <- as.numeric(region_parts[2])
    end <- as.numeric(region_parts[3])

    # construct GRanges object
    region <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges = IRanges::IRanges(start = start - flank_left, end = end + flank_right)
    )
    
    # Import BigWig files as GRanges objects
    message("Importing BigWig files...")
    bw_region <- lapply(bw_files, function(bw_file) {
        rtracklayer::import(bw_file, which = region, as = "GRanges")
    })
    
    # Name the GRanges list using file names without extensions
    message("Naming GRanges list using file names...")
    names(bw_region) <- sub("\\..*$", "", basename(bw_files))
    
    # Convert GRanges to data frame with source identifier
    message("Converting GRanges to data frame...")
    bw_all <- do.call(rbind, lapply(names(bw_region), function(name) {
        bw_df <- as.data.frame(bw_region[[name]])
        bw_df[[by]] <- name
        return(bw_df)
    }))
    
    # Merge with mapping data frame if provided
    if (!is.null(mapping)) {
        message("Merging with mapping...")
        bw_all <- merge(bw_all, mapping, by = by, all.x = TRUE)
    }
    
    return(bw_all)
}



# --------------------------------------------------
#' Plot BigWig data with customizable features
#'
#' This function generates a ggplot2-based visualization of BigWig data, displaying coverage
#' across genomic positions with options for smoothing, faceting, vertical lines, and custom
#' colors. It saves the plot to a specified directory and returns the plot object.
#'
#' @param data Data frame containing BigWig data with columns for genomic position (\code{start}),
#'             coverage (\code{score}), and features for coloring and faceting.
#' @param colour_set Character string naming the column in \code{data} with color values for mapping.
#'                   If invalid or missing, colors are auto-generated based on unique \code{feature_1} levels.
#' @param title Character string for the plot title (default: "chr:start-end").
#' @param feature_1 Character string naming the column in \code{data} for color grouping (default: "day").
#' @param feature_2 Character string naming the column in \code{data} for faceting (default: "week").
#' @param smooth Logical indicating whether to add smoothed lines with \code{geom_smooth} (default: \code{TRUE}).
#' @param span Numeric value controlling the smoothness of the \code{loess} fit when \code{smooth = TRUE} (default: 0.1).
#' @param level Numeric value between 0 and 1 specifying the confidence level for smoothed intervals (default: 0.95).
#' @param se Logical indicating whether to display standard error confidence intervals with smoothing (default: \code{FALSE}).
#' @param vline_values Numeric vector of x-intercept values for vertical dashed lines (default: \code{NULL}, no lines).
#' @param file_name Character string specifying the base file name for the saved plot (default: "bigwig_plot").
#'                  If \code{smooth = TRUE}, appended with smoothing parameters (e.g., "_smooth_span_0.1_level_0.95").
#' @param figure_dir Character string specifying the directory to save the plot (default: "./").
#' @param alpha_point Numeric value for the transparency of points in \code{geom_point} (default: 0.15).
#' @param size_point Numeric value for the size of points in \code{geom_point} (default: 1).
#' @param alpha_smooth Numeric value for the transparency of smoothed lines in \code{geom_smooth} (default: 0.75).
#' @param size_smooth Numeric value for the size of smoothed lines in \code{geom_smooth} (default: 0.75).
#' @param alpha_guides Numeric value for the transparency of legend guides (default: 1).
#' @param size_guides Numeric value for the size of legend guides (default: 3).
#'
#' @return A \code{ggplot} object representing the generated plot, which is also saved to \code{figure_dir}.
#'
#' @details This function creates a scatter plot of coverage (\code{score}) over genomic positions (\code{start}),
#' colored by \code{feature_1}, with optional smoothing and faceting by \code{feature_2}. It validates \code{colour_set}
#' against \code{feature_1} levels, auto-generating colors if needed, and adjusts the file name based on smoothing parameters.
#'
#' @examples
#' \dontrun{
#' data <- data.frame(start = 1:100, score = rnorm(100), day = rep(c("Mon", "Tue"), 50), 
#'                    week = rep(c("W1", "W2"), each = 50), colors = c("red", "blue"))
#' p <- plot_bigwig(data, colour_set = "colors", title = "chr1:1-100", smooth = TRUE, 
#'                  vline_values = c(50), span = 0.2, alpha_smooth = 0.5)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth labs scale_color_manual scale_fill_manual
#' @importFrom ggplot2 coord_cartesian facet_wrap geom_vline theme_minimal theme guides guide_legend
#' @importFrom ggplot2 element_blank ggsave
#' @export
plot_bigwig <- function(data,
                        colour_set = 'color_shade',
                        title = "chr:start-end",
                        feature_1 = "day",
                        feature_2 = "week",
                        smooth = TRUE,
                        span = 0.1,
                        level = 0.95,
                        se = FALSE,
                        vline_values = NULL,
                        file_name = "bigwig_plot",
                        figure_dir = NULL,
                        alpha_point = 0.15,
                        size_point = 1,
                        alpha_smooth = 0.75,
                        size_smooth = 0.75,
                        alpha_guides = 1,
                        size_guides = 3) {
    
    # Check if data is provided and features exist
    if (is.null(data)) stop("Providing data is required.")
    if (!all(c(feature_1, feature_2) %in% names(data))) {
        stop("Both 'feature_1' and 'feature_2' must be present in the data.")
    }
    
    # Message: Starting the plot creation process
    message("Starting the plot creation process.")
   
    
    # Start the ggplot object
    p <- ggplot2::ggplot(data, ggplot2::aes(x = start, y = score)) +
        ggplot2::geom_point(ggplot2::aes(color = .data[[feature_1]]), 
                            size = size_point, 
                            alpha = alpha_point)
    
    # Add smoothed lines if requested
    if (smooth) {
        message("Adding smoothed lines with confidence intervals.")
        p <- p + ggplot2::geom_smooth(ggplot2::aes(color = .data[[feature_1]], 
                                                   fill = .data[[feature_1]]), 
                                      method = "loess", 
                                      span = span, 
                                      size = size_smooth, 
                                      level = level, 
                                      alpha = alpha_smooth, 
                                      se = se)
    }

    
    # Add labels, scales, faceting, vertical lines, and theme
    p <- p +
        ggplot2::labs(x = "Genomic Position", y = "Coverage", title = title) +
        ggplot2::coord_cartesian(ylim = c(0, ceiling(max(data$score, na.rm = TRUE) / 500) * 500)) +
        ggplot2::facet_wrap(stats::as.formula(paste("~", feature_2)), ncol = 1) +
        ggplot2::geom_vline(xintercept = vline_values, linetype = "dashed", color = "red") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_blank())
    
    # Add manual color scales only if colour_set is provided and valid
    if (!is.null(colour_set) && colour_set != "" && colour_set %in% names(data)) {
        p <- p +
            ggplot2::scale_color_manual(values = stats::setNames(data[[colour_set]], data[[feature_1]])) +
            ggplot2::scale_fill_manual(values = stats::setNames(data[[colour_set]], data[[feature_1]]))
    }

    # Customize legend based on smoothing
    if (!smooth) {
        p <- p + ggplot2::guides(fill = "none",
                                 color = ggplot2::guide_legend(override.aes = list(linetype = 1,
                                                                                   size = size_guides,
                                                                                   alpha = alpha_guides)))
    }
    
    # Message: Saving the plot
    message("Saving the plot to the specified directory.")
    
    # Construct file name based on smoothing
    full_file_name <- if (smooth) {
        paste(file_name, "smooth", paste0("span_", span), paste0("level_", level), sep = "_")
    } else {
        file_name
    }
    
    # Save the plot if figure_dir is provided
    if (!is.null(figure_dir)) {
        ggplot2::ggsave(paste0(figure_dir, full_file_name, ".png"), 
                        plot = p, bg = "white", height = 10, width = 16.2, dpi = 320)
    }
    
    # Message: Plotting completed
    message("Plotting completed.")
    return(p)
}
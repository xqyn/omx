# bed_utils.R
# project: omx
# april 2025 - XQ - LUMC
# plot bed graph

library(ggplot2)
library(rtracklayer)


# spacing_yline --------------------------------------------------
#' Calculate spacing lines for genomic isoform
#'
#' Adds spacing (y_line) between genomic isoform based on the number of rows,
#' visualization purposes.
#'
#' @param bed data.frame object or data frame with genomic regions
#'
#' @return The input object with added y_line column containing spacing values
#' @export
#'
#' @examples
#' \dontrun{
#'   gr <- GRanges(seqnames = c("chr1", "chr1"),
#'                 ranges = IRanges(start = c(100, 200), end = c(150, 250)))
#'   gr_with_spacing <- spacing_yline(as.data.frame(gr))
#' }

spacing_yline <- function(bed,
                          y_line = 'y_line') {
  n <- nrow(bed)
  
  if (n >= 1 && n <= 100) {
    # Determine spacing tier (each 5 rows = 1 tier)
    tier <- ceiling(n / 5)
    
    # Calculate spacing: starts at 0.5, increases by 0.25 per tier
    spacing_value <- 0.5 + 0.25 * (tier - 1)
    
    # Cap maximum spacing at 5
    spacing <- min(spacing_value, 5)
    
    # Assign y_line values
    bed[[y_line]] <- seq_len(n) * spacing
  } else {
    bed[[y_line]] <- rep(NA, n)
    message("Reduce number of genes in the region")
  }
  
  return(bed)
}


# process_bed_extract --------------------------------------------------
#' Extract a Genomic Region with Flanking and Intersect with BED
#'
#' Given a genomic region (`chr:start-end` string or `GRanges` object),
#' applies left and right flanking, filters entries from a BED file using bedtools.
#'
#' @param region_str A region in `"chr:start-end"` format or a `GRanges` object.
#' @param flank_left Number of bases to extend to the left (default 0).
#' @param flank_right Number of bases to extend to the right (default = `flank_left` if NULL).
#' @param bed_file Path to the BED file.
#' @param bedtools_path Path to the bedtools executable.
#' @param tmp_dir Directory for temporary files (default is "./tmp/").
#' @param add_spacing_yline Logical, whether to apply `spacing_yline()` to the output (default TRUE).
#'
#' @return A `GRanges` object containing the intersected BED regions.
#' @export

process_bed_extract <- function(region_str = NULL,
                                flank_left = 0,
                                flank_right = NULL,
                                bed_file = NULL,
                                bedtools_path = NULL,
                                tmp_dir = "./tmp/",
                                add_spacing_yline = TRUE) {
    
  if (is.null(region_str)) stop("Require 'region_str' in format 'chr:start-end' or a GRanges object.")
  if (is.null(bed_file)) stop("Required valid 'bed_file' path.")
  if (is.null(bedtools_path)) stop("Require path to the 'bedtools' executable.")
  if (is.null(flank_right)) flank_right <- flank_left

  # convert GRanges to string if needed
  if (inherits(region_str, "GRanges")) {
    chr <- as.character(GenomicRanges::seqnames(region_str))
    start <- GenomicRanges::start(region_str)
    end <- GenomicRanges::end(region_str)
    region_str <- paste0(chr, ":", start, "-", end)
  }

  # parse region string
  region_parts <- unlist(strsplit(region_str, "[:-]"))
  chr <- region_parts[1]
  start <- as.numeric(region_parts[2])
  end <- as.numeric(region_parts[3])

  # apply flanking
  flanked_start <- max(0, start - flank_left)  # Ensure non-negative
  flanked_end <- end + flank_right
    
  # create temporary files
  if (!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = TRUE)
  temp_file <- tempfile(tmpdir = tmp_dir, fileext = ".bed")
  temp_region_file <- file.path(tmp_dir, "temp_region.bed")

  # write region to temp file
  region_bed <- paste(chr, flanked_start, flanked_end, sep = "\t")
  writeLines(region_bed, temp_region_file)

  # run bedtools intersect
  cmd <- paste(shQuote(bedtools_path), "intersect -a", shQuote(bed_file), "-b", shQuote(temp_region_file), ">", shQuote(temp_file))
  system(cmd)

  # read the result
  bed_region <- rtracklayer::import(temp_file, format = "BED")

  # optionally apply spacing
  if (add_spacing_yline && exists("spacing_yline")) {
    bed_region <- spacing_yline(as.data.frame(bed_region))
  }

  # clean up
  file.remove(temp_region_file)
  file.remove(temp_file)

  return(bed_region)
}


# isoform_mapping --------------------------------------------------
#' Isoform Exon and Arrow Mapping
#'
#' Given a single isoform row (from a BED-like dataframe),
#' returns exon block annotations and connecting arrow symbols
#' representing intron connections.
#'
#' @param isoform A single-row data.frame containing isoform (gene) info. 
#'   Required columns: \code{start}, \code{end}, \code{strand}, \code{blocks}, \code{y_line}. 
#'   \code{blocks} should be a list of exon coordinates (start and end relative to isoform start).
#' @param arrow_dist Distance (in bases) between arrows along the intron (default: 300).
#' @param base_offset Base spacing to prevent arrows overlapping with exon ends (default: 50).
#'
#' @return A list with:
#' \describe{
#'   \item{block}{A data.frame of exon block coordinates and arrow boundary positions.}
#'   \item{arrows}{A data.frame with arrow symbols and their plotting coordinates.}
#' }
#' @export
isoform_mapping <- function(isoform, 
                            arrow_dist = 300, 
                            base_offset = 50) {
  
  # check for required columns
  required_cols <- c("start", "end", "strand", "blocks", "y_line")
  missing_cols <- setdiff(required_cols, names(isoform))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
  }

  # --- Handle MULTIPLE isoforms (recursive call)
  if (nrow(isoform) > 1) {
    iso_block <- lapply(seq_len(nrow(isoform)), function(i) {
      isoform_mapping(isoform[i, ],
                      arrow_dist = arrow_dist,
                      base_offset = base_offset)
    })

    block <- do.call(rbind, lapply(iso_block, `[[`, "block"))
    arrow <- do.call(rbind, lapply(iso_block, `[[`, "arrow"))

    return(list(block = block, arrow = arrow))
  }

  # --- BLOCK
  # Extract block information
  # block: one block is one exon
  block <- as.data.frame(isoform$blocks[[1]])
  #  add strand and y_line to block
  block$strand <- isoform$strand
  block$y_line <- isoform$y_line

  # exon postions: |[=====]_i_th----[==]----[======]----[=====]|
  # exon start = isoform start + block start: | + [_i_th
  # exon end = isoform start + block end: | + ]_i_th
  block$exon_start <- isoform$start + block$start
  block$exon_end   <- isoform$start + block$end

  # arrow positions (if not exon -> arrow) to block
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

  # --- ARROW
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
    arrow <- if (length(arrow_position) > 0) {
      data.frame(
        arrow_sign = rep(ifelse(isoform$strand == "+", ">", "<"), length(arrow_position)),
        arrow_pos = arrow_position,
        y_line = rep(isoform$y_line, length(arrow_position))
      )
    } else {
      NULL
    }
    
    return(list(block = block, arrow = arrow))
}


# plot_isoform_map --------------------------------------------------
#' Plot Isoform Map with Exons and Arrows
#'
#' Visualizes the exon structure of isoforms using block (exon) and arrow annotations.
#'
#' @param bed A dataframe containing isoform information, including `start`, `end`, `y_line`, `name`, and `strand`.
#' @param block_df A dataframe of exon block positions, typically generated from an `isoform_mapping()` function.
#' @param arrow_df A dataframe of arrow annotations, also from `isoform_mapping()`.
#' @param region_pad Integer. Base padding added to the left and right of the plotting region (default = 2000).
#'
#' @return A `ggplot` object visualizing the isoform structure with exon blocks and direction arrows.
#' 
#' @importFrom ggplot2 aes geom_segment geom_text ggplot scale_color_manual
#' @importFrom ggplot2 ylim labs scale_x_continuous theme_classic theme
#' @export

plot_isoform_map <- function(bed, 
                            block_df, 
                            arrow_df, 
                            region_pad = 0) {
  
  # Helper function to check required columns
  check_columns <- function(df, required_cols, df_name) {
    missing_cols <- setdiff(required_cols, names(df))
    if (length(missing_cols) > 0) {
      stop(sprintf("Missing required columns in %s: %s", 
                  df_name, 
                  paste(missing_cols, collapse = ", ")))
    }
  }
  
  # Validate input data frames
  check_columns(bed, c("start", "end", "y_line", "name", "strand"), "bed")
  check_columns(block_df, c("exon_start", "exon_end", "y_line", "strand"), "block_df")
  check_columns(arrow_df, c("arrow_pos", "y_line", "arrow_sign"),"arrow_df")
  
  # Create the plot
  p <- ggplot2::ggplot() +

    # Base horizontal lines per gene
    ggplot2::geom_segment(data = bed,
      ggplot2::aes(x = start, xend = end,
                   y = y_line, yend = y_line),
      size = 0.15, 
      alpha = 0.5,
      color = "black"
    ) +

    # block boxed for exons
    ggplot2::geom_segment(data = block_df,
      ggplot2::aes(x = exon_start, xend = exon_end,
                   y = y_line, yend = y_line,
                   color = strand),
      size = 4) +
    ggplot2::scale_color_manual(values = c("+" = "#017075", "-" = "#F48D79")) +

    # arrows between exons
    ggplot2::geom_text(data = arrow_df,
      ggplot2::aes(x = arrow_pos, y = y_line, label = arrow_sign),
      vjust = 0.3,
      size = 1.5,
      color = "black") +

    # Gene names
    ggplot2::geom_text(data = bed,
      ggplot2::aes(x = start, y = y_line, label = name),
      hjust = ifelse(bed$strand == "+", 1.15, 1.15),
      vjust = 0.5,
      size = 2, 
      color = "#0040B0") +
    
    # axis limits
    #ggplot2::xlim(min(bed$start) - region_pad, max(bed$end) + region_pad) +
    ggplot2::ylim(0.25, max(bed$y_line) + 0.75) +
    
    # labels and scaling
    ggplot2::labs(title = NULL, x = "Genomic Position (Kb)", y = NULL) +
    ggplot2::scale_x_continuous(
      limits = c(min(bed$start) - region_pad, max(bed$end) + region_pad),
      labels = function(x) {
        labels <- scales::comma(x * 0.001)
        labels[length(labels)] <- paste0(labels[length(labels)], " Kb")
        labels}
    ) +
    
    # theme adjustments
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank()
    )
  return(p)
}
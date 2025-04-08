# """
# Project: track_genome
# April 2025 - XQ - LUMC
# visualization of genomic regions in R
# """


library(ggplot2)
library(rtracklayer)

# """
# bigWig file processing and plotting
# """

# --------------------------------------------------
# Function to plot BigWig data with optional smoothing and customizable features
# span: controls how much of the data is used for smoothing
# level: interval for confidence interval
## Wide Interval: More uncertainty, broader shading (e.g., sparse data, high variability, higher level).
## Narrow Interval: Less uncertainty, tighter shading (e.g., dense data, low variability, lower level).



# """
# BED file processing and isoform mapping
# """



# bed_region --------------------------------------------------
bed_region <- function(region_string=NULL, 
                    flank_left=0, 
                    flank_right=NULL, 
                    bed_file=NULL,
                    bedtools_path=NULL,
                    tmp_dir = "./tmp/") {
    
    #--------------------------------------------------
    #' Extract genomic region and apply flanking
    #'
    #' Given a  region: chr:start:end string,
    #' applies left and right flanking, and filters BED file entries using bedtools intersect.
    #'
    #' @param region_string Genomic region in "chr:start-end" format
    #' @param flank_left Number of bases to extend on the left
    #' @param flank_right Number of bases to extend on the right
    #' @param bedtools_path Path to the bedtools executable
    #'
    #' @return A GRanges object containing the filtered BED regions
    #--------------------------------------------------
        
    # Check if required arguments are provided
    if (is.null(region_string)) {
        stop("Please provide a valid 'region_string' in the format 'chr:start-end'.")}
    
    if (is.null(bed_file)) {
        stop("Please provide a valid 'bed_file' path.")}
    
    if (is.null(bedtools_path)) {
        stop("Please provide the path to the 'bedtools' executable.")}
    
    # # Check if region_string is in the correct "chr:start-end" format
    # region_pattern <- "^([A-Za-z0-9]+):([0-9]+)-([0-9]+)$"
    # if (!grepl(region_pattern, region_string)) {
    #     stop("Invalid 'region_string' format. It should be in the 'chr:start-end' format.")}

    # Check if the region is a GRanges object
    if (inherits(region_string, "GRanges")) {
        region_string <- granges_to_string(region_string)    }

    # Extract chromosome, start, and end from the input string
    region_parts <- strsplit(region_string, ":|-")[[1]]
    chr <- region_parts[1]
    start <- as.numeric(region_parts[2])
    end <- as.numeric(region_parts[3])
    
    # If flank_right is NULL, set it to be equal to flank_left
    if (is.null(flank_right)) {
    flank_right <- flank_left}

    # Calculate the flanked coordinates
    flanked_start <- start - flank_left
    flanked_end <- end + flank_right
    
    # Create a temporary file in /tmp/ with the region of interest using bedtools intersect
    # Ensure ./tmp directory exists
    if (!dir.exists(tmp_dir)) {
    dir.create(tmp_dir)
    }

    # Create temporary file paths inside ./tmp
    temp_file <- tempfile(tmpdir = tmp_dir, fileext = ".bed")
    temp_region_file <- file.path(tmp_dir, "temp_region.bed")

    # Construct the BED format string and write to file
    region_bed <- paste(chr, flanked_start, flanked_end, sep = "\t")
    writeLines(region_bed, temp_region_file)

    # Use bedtools to intersect and create a new BED file with only the region of interest
    system(paste(bedtools_path, "intersect -a", bed_file, "-b", temp_region_file, ">", temp_file))

    # Import the filtered BED file
    filtered_data <- rtracklayer::import(temp_file, format = "BED")

    # Clean up temporary files
    #file.remove(temp_region_file)
    #file.remove(temp_file)
    
    # Convert to GRanges object
    region <- GRanges(
        seqnames = chr,
        ranges = IRanges(start = flanked_start, end = flanked_end)
    )
    
    return(filtered_data)
}


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
                            arrow_df, 
                            region_pad = 2000) {
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
      ylim(0.25, max(bed$y_line) + 0.75) + 

      # remove plot labels
      labs(title = NULL, x = "Genomic Position", y = NULL) +  
      theme_classic() +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank())

    return(p)
}

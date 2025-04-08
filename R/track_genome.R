# """
# Project: track_genome
# April 2025 - XQ - LUMC
# visualization of genomic regions in R
# """


# april 3, 2025 
# chr_region --------------------------------------------------
chr_region <- function(region_string, flank_left = 0, flank_right = NULL) {
    #--------------------------------------------------
    #' Extract genomic region and apply flanking
    #'
    #' Given a region: chr:start-end string,
    #' applies left and right flanking, and returns a GRanges object
    #' representing the modified region with flanks.
    #'
    #' @param region_string Genomic region in "chr:start-end" format
    #' @param flank_left Number of bases to extend on the left
    #' @param flank_right Number of bases to extend on the right
    #'
    #' @return A GRanges object containing the modified genomic region
    #'         with the applied flanking regions.
    #--------------------------------------------------
    # Check if region_string is provided
    if (is.null(region_string)) stop("Please provide a valid 'region_string' in the format 'chr:start-end'.")
    if (is.null(flank_right)) flank_right <- flank_left

    # Split region string and parse components
    region_parts <- strsplit(region_string, "[:-]")[[1]]
    chr <- region_parts[1]
    start <- as.numeric(region_parts[2])
    end <- as.numeric(region_parts[3])

    # Construct GRanges object
    region <- GRanges(
        seqnames = chr,
        ranges = IRanges(start = start - flank_left, end = end + flank_right)
    )

    cat(sprintf("chr: %s  start: %d  end: %d\n", chr, start - flank_left, end + flank_right))
    return(region)
}


# granges_to_string --------------------------------------------------
granges_to_string <- function(granges_obj, sep = ":") {
    #--------------------------------------------------
    #' Convert GRanges object to a formatted string.
    #'
    #' Given a GRanges object, this function extracts the chromosome, start, and end positions,
    #' and returns a string with a custom separator in the "chr:start-end" format (or custom separator).
    #'
    #' @param granges_obj A GRanges object that contains genomic ranges.
    #'                    The object should have seqnames (chromosome), start, and end values.
    #' @param sep A character string used as the separator between the chromosome, start, and end values.
    #'            Default is ":". You can pass other separators like "_", "-", etc.
    #'
    #' @return A character string representing the genomic region, formatted as "chr:start-end" or
    #'         using the specified separator.
    #--------------------------------------------------

    # Check if the input is a GRanges object
    if (!inherits(granges_obj, "GRanges")) {
        stop("The provided object is not a GRanges object.")
    }
    
    # Extract chromosome, start, and end
    chr <- as.character(seqnames(granges_obj))
    start <- start(granges_obj)
    end <- end(granges_obj)
    
    # Return the region as a chr:start-end string or using the custom separator
    region_string <- paste(chr, start, end, sep = sep)
    return(region_string)
}


# process_bigwig_data --------------------------------------------------
process_bigwig_data <- function(bw_files, 
                               region_str = "X:73851081-73851581",
                               flank_left = 2000, 
                               flank_right = NULL, 
                               mapping = NULL) {
  
    #' Process BigWig files into a merged data frame
    #'
    #' This function takes a set of BigWig files, imports them as GRanges objects for a specified
    #' genomic region with flanking regions, converts them to data frames, and optionally
    #' merges them with a mapping data frame.
    #'
    #' @param bw_files Character vector of paths to BigWig files to be processed
    #' @param region_str Character string specifying the genomic region in "chr:start-end" format
    #'                   (default: "X:73851081-73851581")
    #' @param flank_left Integer specifying the left flanking region size in base pairs (default: 2000)
    #' @param flank_right Integer specifying the right flanking region size in base pairs 
    #'                    (default: NULL, which sets it equal to flank_left)
    #' @param mapping Optional data frame containing mapping information to merge with the processed data
    #'
    #' @return A data frame containing the combined data from all BigWig files,
    #'         optionally merged with mapping information if provided
    #'
    #' @importFrom rtracklayer import
    #' @importFrom GenomicRanges GRanges
    #' @importFrom IRanges IRanges
    #' @importFrom tools file_path_sans_ext
    
    #--- create genomic region with flanks
    message("creating genomic region with flanks...")
    region <- chr_region(region_str, flank_left, flank_right)

    #--- import BigWig files as GRanges objects
    message("importing BigWig files...")
    rg_vars <- lapply(bw_files, function(bw_file) {
                        rtracklayer::import(bw_file, 
                                            which = region, 
                                            as = "GRanges")
    })

    #--- name the GRanges list using file names without extensions
    message("naming GRanges list using file names...")
    names(rg_vars) <- sub("\\..*", "", basename(bw_files))

    #--- convert GRanges to data frame with source identifier
    message("converting GRanges to data frame...")
    combined_df <- do.call(rbind, lapply(names(rg_vars), function(name) {
        df <- as.data.frame(rg_vars[[name]])
        df$code <- name
        return(df)
    }))

    #--- Merge with mapping data frame if provided
    message("merging with mapping...")
    if (!is.null(mapping)) {
    combined_df <- merge(combined_df, mapping, by = "code", all.x = TRUE)
    }

    return(combined_df)
}


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



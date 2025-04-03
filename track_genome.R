# """
# Project: track_genome
# XQ - LUMC
# visualization of genomic regions in R
# """


# april 3, 2025 
# --------------------------------------------------

chr_region <- function(region_string, 
                       flank_left = 0, 
                       flank_right = NULL) {

    #--------------------------------------------------
    #' Extract genomic region and apply flanking
    #'
    #' Given a region: chr:start:end string,
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

    # Extract chromosome, start, and end from the input string
    region_parts <- strsplit(region_string, ":|-")[[1]]

    chr <- region_parts[1]
    start <- as.numeric(region_parts[2])
    end <- as.numeric(region_parts[3])

    # If flank_right is NULL, set it to be equal to flank_left
    if (is.null(flank_right)) {
    flank_right <- flank_left
    }

    # Create the GRanges object with flanking regions
    region <- GRanges(
    seqnames = chr,
    ranges = IRanges(start = start - flank_left, 
                    end = end + flank_right))

    return(region)
}

# --------------------------------------------------

granges_to_string <- function(granges_obj) {
    #--------------------------------------------------
    #' Convert GRanges object to chr:start-end string
    #'
    #' Given a GRanges object, this function extracts the chromosome, start, and end positions,
    #' and returns a string in the "chr:start-end" format.
    #'
    #' @param granges_obj A GRanges object that contains genomic ranges.
    #'                    The object should have seqnames (chromosome), start, and end values.
    #'
    #' @return A character string in the format "chr:start-end" representing the genomic region.
    #--------------------------------------------------

    # Check if the input is a GRanges object
    if (!inherits(granges_obj, "GRanges")) {
        stop("The provided object is not a GRanges object.")
        }
    # Extract chromosome, start, and end
    chr <- as.character(seqnames(granges_obj))
    start <- start(granges_obj)
    end <- end(granges_obj)
    
    # Return the region as a chr:start-end string
    region_string <- paste(chr, start, end, sep=":")
    return(region_string)
}


# --------------------------------------------------
bed_region <- function(region_string=NULL, 
                    flank_left=0, 
                    flank_right=NULL, 
                    bed_file=NULL,
                    bedtools_path=NULL) {
    
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
    if (inherits(region, "GRanges")) {
        region_string <- granges_to_string(region)    }

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
    
    # Create a temporary file with the region of interest using bedtools intersect
    temp_file <- tempfile(fileext = ".bed")
    temp_region_file <- "temp_region.bed"
    
    # Construct the BED format string and write to file
    region_bed <- paste(chr, flanked_start, flanked_end, sep="\t")
    writeLines(region_bed, temp_region_file)
    
    # Use bedtools to intersect and create a new BED file with only the region of interest
    system(paste(bedtools_path, "intersect -a", bed_file, "-b", temp_region_file, ">", temp_file))
    
    # Import the filtered BED file
    filtered_data <- import(temp_file, format = "BED")
    
    # Clean up temporary files
    file.remove(temp_region_file)
    file.remove(temp_file)
    
    # Convert to GRanges object
    region <- GRanges(
        seqnames = chr,
        ranges = IRanges(start = flanked_start, end = flanked_end)
    )
    
    return(filtered_data)
}



# region_utils.R
# project: omx
# april 2025 - XQ - LUMC
# modifying genomic regions format

# chr_region --------------------------------------------------
#' Extract genomic region and apply flanking
#'
#' Given a region in "chr:start-end" format, applies left and right flanking, 
#' returns a GRanges object representing the modified region with flanks.
#'
#' @param region_string genomic region in "chr:start-end" format (e.g., "chr1:1000-2000").
#' @param flank_left number of bases to extend on the left (default: 0).
#' @param flank_right number of bases to extend on the right (defaults: \code{flank_left} if not specified).
#'
#' @return a \code{GRanges} object containing the modified genomic region with the applied flanking regions.
#'
#' @details This function parses a genomic region string, adjusts the start and end coordinates based on the
#' specified flanking values, and constructs a \code{GRanges} object. It assumes the input string is well-formed.
#'
#' @examples
#' chr_region("chr1:1000-2000", flank_left = 100, flank_right = 200)
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @export

chr_region <- function(region_string, flank_left = 0, flank_right = NULL) {

    if (is.null(region_string)) stop("Please provide a region_string in format 'chr:start-end'.")
    if (is.null(flank_right)) flank_right <- flank_left

    # split region string and parse components
    region_parts <- strsplit(region_string, "[:-]")[[1]]
    chr <- region_parts[1]
    start <- as.numeric(region_parts[2])
    end <- as.numeric(region_parts[3])

    # Construct GRanges object
    region <- GenomicRanges::GRanges(
        seqnames = chr,
        ranges = IRanges::IRanges(start = start - flank_left, end = end + flank_right)
    )

    cat(sprintf("chr: %s  start: %d  end: %d\n", chr, start - flank_left, end + flank_right))
    return(region)
}


# granges_to_string --------------------------------------------------
#' Convert GRanges object to a formatted string
#'
#' Given a GRanges object, this function extracts the chromosome, start, and end positions,
#' and returns a string with a custom separator in the "chr:start-end" format (or using a specified separator).
#'
#' @param granges_obj a \code{GRanges} object containing genomic ranges with seqnames (chromosome), start, and end values.
#' @param sep A character string used as the separator between chromosome, start, and end values (default: ":").
#'
#' @return A character string representing the genomic region, formatted as "chr:start-end" or using the specified separator.
#'
#' @details This function extracts the necessary components from a \code{GRanges} object and constructs a string.
#' It checks that the input is a valid \code{GRanges} object before proceeding.
#'
#' @examples
#' library(GenomicRanges)
#' gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1000, end = 2000))
#' granges_to_string(gr, sep = "-")  # Returns "chr1-1000-2000"
#'
#' @importFrom GenomicRanges seqnames start end
#' @export

granges_to_string <- function(granges_obj, sep = ":") {
    # Check if the input is a GRanges object
    if (!inherits(granges_obj, "GRanges")) {
        stop("The provided object is not a GRanges object.")
    }
    
    # Extract chromosome, start, and end
    chr <- as.character(GenomicRanges::seqnames(granges_obj))
    start <- GenomicRanges::start(granges_obj)
    end <- GenomicRanges::end(granges_obj)
    
    # Return the region as a chr:start-end string or using the custom separator
    region_string <- paste(chr, start, end, sep = sep)
    return(region_string)
}
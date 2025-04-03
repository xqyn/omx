"""
XQ - LUMC
visualization of genomic regions in R
"""


# --------------------------------------------------
# april 3

chr_region <- function(
    region_string, 
    flank_left, 
    flank_right, 
    bed_file,
    bedtools_path) {
  
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
  
  # Extract chromosome, start, and end from the input string
  region_parts <- strsplit(region_string, ":|-")[[1]]
  chr <- region_parts[1]
  start <- as.numeric(region_parts[2])
  end <- as.numeric(region_parts[3])
  
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
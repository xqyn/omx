# """
# Project: track_genome
# April 2025 - XQ - LUMC
# visualization of genomic regions in R
# """

# colour
# --------------------------------------------------
library(RColorBrewer)
# global varibles


colour.hm <- c(brewer.pal(9,'Blues')[9:3],brewer.pal(9,'Reds')[3:9])
colours <- c("#1D263B","#E03E3E","#53BB79","#EF8228","#937264",
            "#3043A2","#C25D7B","#D88C4C","dodgerblue","#FEB447",
            "#06CDE0","#1498BE","#0D5B11", '#4c00b0', '#ffaf00')
colour.day <- c("#A6E0FF","#4DB1FF","#0070F2","#0040B0")
colour.cell <- c("#1B90FF","#7858FF","#F31DED","#FA4F96")
colour.rep <- c('#F48D79', '#017075')


# --------------------------------------------------

generate_scheme <- function(base_color, 
                            n = 5, 
                            light_range = c(0.4, 0.1), 
                            dark_range = c(0.1, 0.4)) {
    #--------------------------------------------------
    #' Generate a Monochromatic Color Scheme
    #'
    #' Return a vector of lighter and darker shades based on a base color.
    #' 
    #' @param base_color    hex color string (e.g., "#7FC97F").
    #' @param n             mumber of colors scheme to generate  (default is 5).
    #' @param light_range   A numeric vector of length 2 specifying the range of lightening (default = c(0.4, 0.1)).
    #' @param dark_range    A numeric vector of length 2 specifying the range of darkening (default = c(0.1, 0.4)).
    #'
    #' @return A character vector of hex color strings representing the color scheme.
    #' @export
    #'
    #' @examples
    #' generate_scheme("#7FC97F")
    #' generate_scheme("#BEAED4", n = 7, light_range = c(0.6, 0.2), dark_range = c(0.2, 0.6))    

    n_light <- ceiling(n / 2)
    n_dark <- floor(n / 2)
    
    lighten_vals <- seq(light_range[1], light_range[2], length.out = n_light)
    darken_vals <- seq(dark_range[1], dark_range[2], length.out = n_dark)
    
    light_colors <- sapply(lighten_vals, function(a) colorspace::lighten(base_color, amount = a))
    dark_colors <- sapply(darken_vals, function(a) colorspace::darken(base_color, amount = a))
    
    scheme <- c(light_colors, dark_colors)
    return(scheme)
}

#colour <- c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#C8C8C8", "#06CDE0")
# Generate color schemes
#colour_schemes <- lapply(colour, generate_scheme)
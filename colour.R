# """
# Project: track_genome
# April 2025 - XQ - LUMC
# visualization of genomic regions in R
# """

# colour
# --------------------------------------------------
library(RColorBrewer)
# global varibles
colours <- c("#1D263B","#E03E3E","#53BB79","#EF8228","#937264",
            "#3043A2","#C25D7B","#D88C4C","dodgerblue","#FEB447",
            "#06CDE0","#1498BE","#0D5B11", '#4c00b0', '#ffaf00')
colour.day <- c("#A6E0FF","#4DB1FF","#0070F2","#0040B0")
colour.cell <- c("#1B90FF","#7858FF","#F31DED","#FA4F96")
colour.rep <- c('#F48D79', '#017075')
colour.hm <- c(brewer.pal(9,'Blues')[9:3],brewer.pal(9,'Reds')[3:9])

# --------------------------------------------------
generate_scheme <- function(base_color, 
                            n = 5, 
                            light_range = c(0.4, 0.1), 
                            dark_range = c(0.1, 0.4)) {
        # --------------------------------------------------
        #' Generate Monochromatic Color Schemes
        #'
        #' Return one or more vectors of lighter and darker shades based on one or more base colors.
        #' If multiple base colors are supplied, a list of color schemes will be returned.
        #'
        #' @param base_color  A hex color string or a vector of hex strings (e.g., "#7FC97F" or c("#7FC97F", "#BEAED4")).
        #' @param n           Number of colors in each scheme (default is 5).
        #' @param light_range A numeric vector of length 2 specifying the range of lightening (default = c(0.4, 0.1)).
        #' @param dark_range  A numeric vector of length 2 specifying the range of darkening (default = c(0.1, 0.4)).
        #'
        #' @return A character vector of hex colors if one base color is provided, or a named list of such vectors if multiple.
        #' @export
        #'
        #' @examples
        #' generate_scheme("#7FC97F")
        #' generate_scheme(c("#7FC97F", "#BEAED4"), n = 7, light_range = c(0.6, 0.2), dark_range = c(0.2, 0.6))
        # --------------------------------------------------
        # internal function for a single color
        generate_single <- function(clr) {
        n_light <- ceiling(n / 2)
        n_dark <- floor(n / 2)

        lighten_vals <- seq(light_range[1], light_range[2], length.out = n_light)
        darken_vals  <- seq(dark_range[1], dark_range[2], length.out = n_dark)

        light_colors <- sapply(lighten_vals, function(a) colorspace::lighten(clr, amount = a))
        dark_colors  <- sapply(darken_vals,  function(a) colorspace::darken(clr, amount = a))

        c(light_colors, dark_colors)
        }

        # check if base_color is a vector of length 1 or more
        if (length(base_color) == 1) {
        return(generate_single(base_color))
        } else {
        schemes <- lapply(base_color, generate_single)
        names(schemes) <- base_color
        return(schemes)
        }
}

#colour <- c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#C8C8C8", "#06CDE0")
# Generate color schemes
#colour_schemes <- lapply(colour, generate_scheme)
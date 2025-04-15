# fractalome.R
# project: omx
# april 15 2025 - XQ - LUMC
# circulizing the linear


# compute_radius --------------------------------------------------
#' Calculate angles, radii, and reference radii for data points based on position and feature score
#'
#' @param data A data frame containing chromosome, position, and feature columns.
#' @param arc Character string specifying the current arc (e.g., chromosome "1").
#' @param theta_coor Numeric value for the starting angle of chromosome arc (in radians).
#' @param angle_per_bp Numeric value for radians per base pair.
#' @param inner_radius Numeric value for the inner edge of the plotting band.
#' @param outer_radius Numeric value for the outer edge of the plotting band.
#' @param arc_col Character string for the name of the chromosome column (default: "chr").
#' @param unit_col Character string for the name of the position column (default: "pos").
#' @param feature Character string for the name of the feature column (default: "pc1").
#' @param y_line Numeric vector of values to compute reference radii for (default: 1).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{data_angle}: A data frame with columns 'angle', 'radius', 'x', 'y'.
#'   \item \code{y_line}: A list of tuples with reference values and their corresponding radii.
#' }
#'
#' @examples
#' data <- data.frame(chr = c("1", "1"), pos = c(100, 200), pc1 = c(0.5, 0.7))
#' result <- compute_radius(data, chrom = "1", theta_coor = 0, angle_per_bp = 0.01,
#'                         inner_radius = 0.5, outer_radius = 1.0)
#'
#' @importFrom dplyr filter select mutate
#' @export
compute_radius <- function(data, 
                        arc, 
                        theta_coor, 
                        angle_per_bp, 
                        inner_radius, 
                        outer_radius,
                        arc_col = "chr", 
                        unit_col = "pos", 
                        feature = "pc1", 
                        y_line = c(0)) {
  
  # Validate inputs
  if (!is.data.frame(data)) stop("data must be a data frame")
  if (!(arc_col %in% names(data))) stop("arc_col not found in data")
  if (!(unit_col %in% names(data))) stop("unit_col not found in data")
  if (!(feature %in% names(data))) stop("feature not found in data")
  
  # Filter data for this chromosome
  data_angle <- data %>%
    dplyr::filter(!!sym(arc_col) == arc) %>%
    dplyr::select(!!sym(unit_col), !!sym(feature))
  
  # Calculate angle: position along chromosome -> angle (counterclockwise)
  data_angle <- data_angle %>%
    dplyr::mutate(angle = theta_coor - !!sym(unit_col) * angle_per_bp)
  
  # Scale score to radius: map feature range to [inner_radius, outer_radius]
  score_min <- min(data[[feature]], na.rm = TRUE)
  score_max <- max(data[[feature]], na.rm = TRUE)
  
  if (score_max == score_min) {
    # Avoid division by zero
    data_angle <- data_angle %>%
      dplyr::mutate(radius = (inner_radius + outer_radius) / 2)
    y_line <- lapply(y_line, function(val) list(value = val, radius = (inner_radius + outer_radius) / 2))
  } else {
    # Scale scores to radius
    data_angle <- data_angle %>%
      dplyr::mutate(radius = inner_radius + 
                      ((!!sym(feature) - score_min) * (outer_radius - inner_radius) / 
                       (score_max - score_min)))
    
    # Calculate reference radii for specified values
    radius_range <- outer_radius - inner_radius
    score_range <- score_max - score_min
    y_line <- lapply(y_line, function(val) {
      list(value = val, 
           radius = inner_radius + (val - score_min) * radius_range / score_range)
    })
  }
  
  # Convert to x, y coordinates
  data_angle <- data_angle %>%
    dplyr::mutate(x = radius * cos(angle),
                  y = radius * sin(angle))
  
  # Return results
  return(list(data_angle = data_angle, y_line = y_line))
}


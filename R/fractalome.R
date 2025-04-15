# fractalome.R
# project: omx
# april 15 2025 - XQ - LUMC
# circulizing the linear


# compute_radius --------------------------------------------------
#' Calculate angles, inner and outer band
#'
#' @param data a data frame containing arc_col, proj_col, and feat_col columns.
#' @param arc character string specifying the current arc of arc_col (e.g., chromosome "1").
#' @param theta_coor numeric value for starting angle of chromosome arc (in radians).
#' @param angular_unit angular unit for each unit of arc (for e.g.: length (base) of the chromosome)
#' @param inner_radius numeric value for inner edge of the plotting band.
#' @param outer_radius numeric value for outer edge of the plotting band.
#' @param arc_col Character string for the name of the chromosome column (default: "chr").
#' @param proj_col projection column (default: position column  "pos").
#' @param feat_col feature column (default: principal component 1 "pc1").
#' @param y_line number of lines on the projected y-axis (default: 1).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{data_angle}:  data frame with columns 'angle', 'radius', 'x_theta', 'y_theta'.
#'   \item \code{y_line}: data frame with y_line value and radius
#' }
#'
#' @examples
#' data <- data.frame(chr = c("1", "1"), pos = c(100, 200), pc1 = c(0.5, 0.7))
#' result <- compute_radius(data, arc = "1", theta_coor = 0, angular_unit = 0.01,
#'                          inner_radius = 0.5, outer_radius = 1.0)
#'
#' @importFrom dplyr filter select mutate
#' @export
compute_radius <- function(data, 
                           arc, 
                           theta_coor, 
                           angular_unit, 
                           inner_radius, 
                           outer_radius,
                           arc_col = "chr", 
                           proj_col = "pos", 
                           feat_col = "pc1", 
                           y_line = c(0)) {
    
    # Validate inputs
    if (!is.data.frame(data)) stop("data must be a data frame")
    if (!(arc_col %in% names(data))) stop("arc_col not found in data")
    if (!(proj_col %in% names(data))) stop("proj_col not found in data")
    if (!(feat_col %in% names(data))) stop("feat_col not found in data")
    if (!is.numeric(theta_coor) || length(theta_coor) != 1) stop("theta_coor must be a single numeric value")
    if (!is.numeric(angular_unit) || length(angular_unit) != 1) stop("angular_unit must be a single numeric value")
    if (!is.numeric(inner_radius) || inner_radius <= 0) stop("inner_radius must be a positive numeric value")
    if (!is.numeric(outer_radius) || outer_radius <= inner_radius) stop("outer_radius must be greater than inner_radius")
    
    
    # subset dataset for given arc and keep two proj_col feat_cols
    data_angle <- data[data[[arc_col]] == arc, c(arc_col, proj_col, feat_col)]

    # calculate angle: position along chromosome -> angle (counterclockwise)
    data_angle[['angle']] <- theta_coor - data_angle[[proj_col]] * angular_unit
    
    # scale score to radius: map feat_col range to [inner_radius, outer_radius]
    score_min <- min(data[[feat_col]], na.rm = TRUE)
    score_max <- max(data[[feat_col]], na.rm = TRUE)
    
    if (score_max == score_min) {
      # Avoid division by zero
      data_angle[['radius']] <- (inner_radius + outer_radius) / 2
      y_line <- lapply(y_line, function(val) list(value = val, radius = (inner_radius + outer_radius) / 2))
    } else {
      # Scale scores to radius
      data_angle[['radius']] <- inner_radius + 
        ((data_angle[[feat_col]] - score_min) * (outer_radius - inner_radius) / 
        (score_max - score_min))

      # Calculate reference radii for specified values
      # y_line <- lapply(y_line, function(val) {
      #   list(value = val,
      #       radius = inner_radius + (val - score_min) * (outer_radius - inner_radius) / (score_max - score_min))
      # })
      y_line <- data.frame(
        y_line = unlist(y_line),
        radius = inner_radius + (unlist(y_line) - score_min) * (outer_radius - inner_radius) / (score_max - score_min))
    }

    # Convert to x, y coordinates
    data_angle[['x_theta']] <- data_angle[['radius']] * cos(data_angle[['angle']])
    data_angle[['y_theta']] <- data_angle[['radius']] * sin(data_angle[['angle']])
          
    # cover y ti data frame
    #y_line <- data.frame()

  # Return results
  return(list(data_angle = data_angle, y_line = y_line))
}


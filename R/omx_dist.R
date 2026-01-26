# project: omx
# XQ - LUMC
# jan 26 - 2026
# distance, similarity and plot corr

#' Calculate Pairwise Similarity Matrix
#'
#' Computes pairwise similarity between samples based on their top-ranked features
#' or features above a threshold. Supports multiple distance/similarity metrics.
#'
#' @param input_table Data frame or matrix with features as rows and samples as columns.
#'   All columns must be numeric.
#' @param ntop Integer specifying the number of top features to select per sample.
#'   Only used when `lowest_score = NULL` (default: 100).
#' @param lowest_score Numeric threshold for feature selection. If provided, features
#'   with values >= this threshold are selected instead of using `ntop` (default: NULL).
#' @param method Character string specifying the similarity/distance metric. Options:
#'   \itemize{
#'     \item "jaccard": Jaccard index (intersection over union)
#'     \item "overlap": Overlap coefficient (intersection over minimum set size)
#'     \item "dice": Dice coefficient (2 * intersection over sum of set sizes)
#'     \item "cosine": Cosine similarity (for numeric vectors)
#'     \item "pearson": Pearson correlation coefficient
#'     \item "spearman": Spearman rank correlation coefficient
#'     \item "euclidean": Euclidean distance (converted to similarity)
#'     \item "manhattan": Manhattan distance (converted to similarity)
#'   }
#'   (default: "jaccard").
#'
#' @return A symmetric matrix of pairwise similarity/correlation values between samples.
#'
#' @examples
#' # Using top N features
#' dist_similarity(expr_data, ntop = 100, method = "jaccard")
#' 
#' # Using threshold
#' dist_similarity(expr_data, lowest_score = 5, method = "dice")
#' 
#' # Using correlation
#' dist_similarity(expr_data, method = "pearson")
#'
#' @export
dist_similarity <- function(input_table, 
                            ntop = 100, 
                            lowest_score = NULL,
                            method = "jaccard") {
  
  # Validate method
  valid_methods <- c("jaccard", "overlap", "dice", "cosine", "pearson", "spearman", 
                     "euclidean", "manhattan")
  if (!method %in% valid_methods) {
    stop("Invalid method. Choose from: ", paste(valid_methods, collapse = ", "))
  }
  
  # Initialize results matrix
  results <- data.frame(matrix(ncol = ncol(input_table), 
                               nrow = ncol(input_table)))
  rownames(results) <- colnames(input_table)
  colnames(results) <- colnames(input_table)
  
  # For correlation methods, compute directly
  if (method %in% c("pearson", "spearman")) {
    message(paste0("Computing ", method, " correlation..."))
    results <- as.data.frame(cor(input_table, method = method))
    return(results)
  }
  
  # For distance methods (euclidean, manhattan), compute directly
  if (method %in% c("euclidean", "manhattan")) {
    message(paste0("Computing ", method, " distance..."))
    dist_matrix <- as.matrix(dist(t(input_table), method = method))
    # Convert distance to similarity (smaller distance = higher similarity)
    max_dist <- max(dist_matrix)
    results <- as.data.frame(1 - (dist_matrix / max_dist))
    return(results)
  }
  
  # For set-based and cosine methods
  message(paste0("Computing ", method, " similarity using ",
                ifelse(is.null(lowest_score), 
                      paste0("top ", ntop, " features"),
                      paste0("threshold >= ", lowest_score)), "..."))
  
  for (col1 in colnames(input_table)) {
    # Select features for column 1
    if (is.null(lowest_score)) {
      r1 <- rownames(input_table)[order(input_table[, col1], 
                                        decreasing = TRUE)[1:ntop]]
      v1 <- input_table[r1, col1]
    } else {
      r1 <- rownames(input_table)[input_table[, col1] >= lowest_score]
      v1 <- input_table[r1, col1]
    }
    
    for (col2 in colnames(input_table)) {
      # Select features for column 2
      if (is.null(lowest_score)) {
        r2 <- rownames(input_table)[order(input_table[, col2], 
                                          decreasing = TRUE)[1:ntop]]
        v2 <- input_table[r2, col2]
      } else {
        r2 <- rownames(input_table)[input_table[, col2] >= lowest_score]
        v2 <- input_table[r2, col2]
      }
      
      # Calculate similarity based on method
      value <- switch(method,
        "jaccard" = {
          length(intersect(r1, r2)) / length(unique(c(r1, r2)))
        },
        "overlap" = {
          length(intersect(r1, r2)) / min(length(r1), length(r2))
        },
        "dice" = {
          2 * length(intersect(r1, r2)) / (length(r1) + length(r2))
        },
        "cosine" = {
          common_features <- union(r1, r2)
          vec1 <- rep(0, length(common_features))
          vec2 <- rep(0, length(common_features))
          names(vec1) <- common_features
          names(vec2) <- common_features
          vec1[r1] <- v1
          vec2[r2] <- v2
          sum(vec1 * vec2) / (sqrt(sum(vec1^2)) * sqrt(sum(vec2^2)))
        }
      )
      
      results[col1, col2] <- value
    }
  }
  
  return(results)
}


#' Plot Correlation/Similarity Matrix as Heatmap
#'
#' Creates a heatmap visualization of a correlation or similarity matrix using ggplot2.
#' The plot is saved to a specified file path.
#'
#' @param corr_table Data frame or matrix containing correlation/similarity values.
#'   Row names should contain sample identifiers.
#' @param min Numeric minimum value for the color scale (default: 0).
#' @param max Numeric maximum value for the color scale (default: 1).
#' @param value Character string specifying the name for the legend and value column
#'   (default: "value").
#' @param figure_dir Character string specifying the output file path for saving the plot
#'   (default: "./test.png").
#'
#' @return A ggplot object representing the correlation heatmap.
#'
#' @examples
#' # Basic usage
#' plot_corr(similarity_matrix)
#' 
#' # Custom parameters
#' plot_corr(similarity_matrix, 
#'           min = 0, 
#'           max = 1, 
#'           value = "Jaccard Index",
#'           figure_dir = "./results/similarity_heatmap.png")
#'
#' @export
plot_corr <- function(corr_table,
                      min = 0,
                      max = 1,
                      value = "value",
                      figure_dir = "./test.png") {
  
  # Prepare data frame with row names as a column
  corr_df <- corr_table
  corr_df$Var1 <- rownames(corr_df)
  
  # Melt data frame to long format for ggplot2
  melted_corr <- reshape2::melt(corr_df,
                                id.vars = "Var1",
                                variable.name = "Var2",
                                value.name = value)
  
  # Create heatmap
  p <- ggplot(melted_corr, aes(x = Var2, y = Var1, fill = .data[[value]])) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(low = "white", 
                         high = "darkblue",
                         mid = "lightblue",
                         midpoint = (min + max) / 2,
                         limits = c(min, max),
                         name = paste0(value, " Index")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          legend.position = "right") +
    coord_fixed()
  
  # Save plot to file
  ggsave(figure_dir, plot = p, width = 10, height = 10)
  
  return(p)
}
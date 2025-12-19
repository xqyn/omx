##############################
# 2024 - XQ - LUMC
# DEG analysis
##############################

# check_load <- function(pkg) {
#   if (!requireNamespace(pkg, quietly = TRUE)) {
#     message(paste("Package", pkg, "is not installed. Installing now..."))
#     install.packages(pkg, dependencies = TRUE)
#   }
  
#   # Load the package silently
#   suppressPackageStartupMessages(library(pkg, character.only = TRUE))
#   message(paste("Package", pkg, "loaded successfully."))
# }

# # Custumize function for DEseq2 and PCA
# #check_load('docstring')
library('ggplot2')

#--------------------------------------------------
# return annotation UP and Down genes from DEseq2
deg <- function(
            df, 
            lfc = 1, 
            padj = 0.05){

    #--------------------------------------------------
    #' DEseq2 annotation
    #'
    #' Return filtering of UP and DOWN genes from DEseq2 in diffexpressed columns
    #'
    #' @param df dataframe from DEseq2 results
    #' @param lfc log2FoldChange threshold
    #' @param padj p-adjusted threshold
    #'
    #' @return dataframe with annotation of UP and DOWN genes
    #--------------------------------------------------

    df_de <- as.data.frame(df)                      # transfer to dataframe
    df_de <- df_de[!is.na(df_de$log2FoldChange),]   # filter lfc 
    df_de <- df_de[!is.na(df_de$padj),]             # filter p adjusted
    df_de$deg <- 'NO DE'                  # re-filter for no DEG
    df_de$deg[df_de$log2FoldChange > lfc & df_de$padj < padj] <- "UP"
    df_de$deg[df_de$log2FoldChange < -lfc & df_de$padj < padj] <- "DOWN"
    df_de <- df_de[order(df_de$log2FoldChange, decreasing=TRUE),] 
    return(df_de)
}


# Merge function --------------------------------------------------
merge_func <- function(
    df_1,
    df_2,
    merge_by = "row.names",
    sort_column = "log2FoldChange",
    sort_decreasing = TRUE
) {
  #'
  #' Merge two data frames with sorting
  #'
  #' Merges two data frames by a specified column and sorts by another column
  #'
  #' @param df_1 First input data frame
  #' @param df_2 Second input data frame
  #' @param merge_by Column name to merge by (default: "row.names")
  #' @param sort_column Column name to sort by (default: "log2FoldChange")
  #' @param sort_decreasing Boolean for descending sort (default: TRUE)
  #'
  #' @return Merged and sorted data frame
  #'
  
  df_1 <- as.data.frame(df_1)
  df_2 <- as.data.frame(df_2)
  
  # Merge the data frames
  df_merge <- merge(df_1, df_2, by = merge_by)
  
  # Sort the merged data frame
  if (!sort_column %in% colnames(df_merge)) {
    stop("The sort_column does not exist in the merged data frame.")
  }
  df_merge <- df_merge[order(df_merge[[sort_column]], 
                             decreasing = sort_decreasing), ]
  
  return(df_merge)
}


# PCA analysis --------------------------------------------------
pca_analysis <- function(
                  mtx,
                  pca_setting = NULL,
                  pca_title = NULL,
                  condition = NULL,
                  replicate = NULL,
                  cond_title = 'feature1',
                  rep_title = 'feature2',
                  ntop = 0.1,
                  scale.unit = FALSE,
                  ncp = 25,
                  fig_dir = NULL,
                  colour_set = NULL,
                  height = 6,
                  width = 7.5,
                  pdf = FALSE){

    #--------------------------------------------------
    #' PCA analysis
    #'
    #' Return PCA analysis output and figures
    #'
    #' @param mtx matrix array with colnames for samples
    #' @param pca_setting naming the pca run
    #' @param pca_title naming the pca title when plotting
    #' @param condition sample conditions
    #' @param replicate sample replicates
    #' @param ntop percentage of top variants
    #' @param scale.unit scale data when perform PCA
    #' @param ncp number of PCs for PCA analysis
    #' @param fig_dir saving directory for figures
    #' @param colours colour for conditions
    #' @param height height for PCA plot
    #' @param width width for PCA plot
    #' @param pdf to print pdf format or not
    #'
    #' @return list of attributes for PCA: pca dataframe, variance explain, and attributes
    #--------------------------------------------------
    
    if (is.null(pca_setting)) stop("Please provide PCA setting")
    if (is.null(pca_title)) {pca_title <- pca_setting}
    if (is.null(condition)) stop("Please provide condition")
    if (is.null(replicate)) stop("Please provide replicate")
    message('PCA analysis processing...')

    mtx <- as.matrix(mtx)
    sample <- colnames(mtx)
    replicate <- as.factor(replicate)
    condition <- as.factor(condition)
    # print out the samples
    message('Run samples:', toString(sample))
    message(cond_title, ': ', toString(condition))
    message(rep_title, ': ', toString(replicate))

    # Determine the number of top variant variables
    ntop <- round(nrow(mtx)*ntop,0)                     # choose ntop % most variants
    row_var <- matrixStats::rowVars(mtx)                # calculate the variance for each rows
    top_var <- order(row_var, decreasing=TRUE)[1:ntop]  # select the ntop variants by variance

    # perform a PCA FactoMineR 
    res.pca <- FactoMineR::PCA(t(mtx[top_var,]), 
                    scale.unit = scale.unit, 
                    ncp = ncp,
                    graph = FALSE)
    
    # PCA mtx and add infor from meta file to pca matrix
    pca <- as.data.frame(res.pca$ind$coord)
    pca[[cond_title]] <-  condition
    pca[[rep_title]] <-  replicate

    # variances explain from PCA mtx
    var <- as.data.frame(res.pca$eig)
    colnames(var) <- c('eigen', 'per_var', 'cumulative')
    var$per_var <- round(var$per_var,2)
    
    # features contribute to PC
    attri <- as.data.frame(res.pca$var$contrib)
    colnames(attri) <- sub('Dim.','pc', colnames(attri))
    attri <- attri[order(attri$pc1, decreasing=TRUE),] # sort on pc1
    attri_ori <- as.data.frame(res.pca$var$contrib)

    # loadings 
    loadings <- as.data.frame(res.pca$var$coord)
    colnames(loadings) <- sub('Dim.','pc', colnames(loadings))
    
    if (is.null(fig_dir)) {
        message("Not saving figures")
    } else {
        message("Ploting PCA figures...")
        # PCA plot
        plot_pca <- ggplot2::ggplot(pca, ggplot2::aes(x = Dim.1, y = Dim.2)) +
                        ggplot2::geom_point(ggplot2::aes(shape=.data[[rep_title]],
                                    colour = .data[[cond_title]]), 
                                    size = 5, 
                                    alpha=0.75, 
                                    show.legend = TRUE) +
                        ggplot2::labs(x = paste0("PC1: ", var$per_var[1], "%"),
                            y = paste0("PC2: ", var$per_var[2], "%"),
                            color = cond_title,
                            shape = rep_title
                            ) +
                        ggplot2::scale_color_manual(values = colour_set) + 
                        ggplot2::ggtitle(paste0('PCA: ', pca_setting)) +
                        ggplot2::theme_minimal()
        # save
        ggplot2::ggsave(paste0(fig_dir, pca_setting, '_pca.png'),
            plot=plot_pca, bg = 'white', width=width, height=height, units="in", dpi=600)
        if (pdf == TRUE){
        ggplot2::ggsave(paste0(fig_dir, pca_setting, '_pca.pdf'),
            plot=plot_pca, bg = "white", width=width, height=height, units="in", dpi= 1200, device = grDevices::cairo_pdf)
            }
            

        # Varaince plot
        plot_var <- ggplot2::ggplot(var, ggplot2::aes(x = factor(1:nrow(var)), y = per_var)) +
                        ggplot2::geom_bar(stat = "identity") +
                        ggplot2::labs(x = "PC", 
                            y = "Percent variance", 
                            title = paste0("Percent varience explained: ", pca_setting)) +
                        ggplot2::theme_minimal()
        ## save
        ggplot2::ggsave(paste0(fig_dir, pca_setting, '_pca_var.png'), 
            plot=plot_var, bg = 'white', width=6.5, height=6, units="in", dpi=600)
        if (pdf == TRUE) {
            ggplot2::ggsave(paste0(fig_dir, pca_setting, '_pca_var.pdf'), 
                plot=plot_var, bg = "white", width=6.5, height=6,units="in", dpi= 1200, device = grDevices::cairo_pdf)
                }
    }
    return(list(pca = pca, var = var, attri = attri, attri_ori = attri_ori, loadings = loadings))  
  }


# find_percentile --------------------------------------------------
#' Find Percentile Thresholds
#'
#' Calculates threshold values for top and bottom percentiles of a numeric vector.
#'
#' @param x A numeric vector to analyze.
#' @param threshold A numeric value between 0 and 1 specifying the percentile threshold (as a decimal) for both top and bottom percentiles. Defaults to \code{0.1} (10\%).
#'
#' @return A named numeric vector with two elements: \code{bottom} (value at the lower percentile) and \code{top} (value at the upper percentile).
#'
#' @examples
#' \dontrun{
#'   x <- rnorm(100)
#'   thresholds <- find_percentile(x, threshold = 0.1)
#'   print(thresholds)
#' }
#'
#' @export
find_percentile <- function(x, threshold = 0.1) {
  # Validate inputs
  if (!is.numeric(x)) stop("Input must be a numeric vector.")
  if (!is.numeric(threshold) || threshold <= 0 || threshold >= 1) {
    stop("Threshold must be a numeric value between 0 and 1.")
  }

  # Sort the values
  sorted_values <- sort(x)

  # Calculate thresholds
  n <- length(sorted_values)
  bottom_th <- sorted_values[ceiling(n * threshold)]
  top_th <- sorted_values[floor(n * (1 - threshold))]

  # Return thresholds as a named vector
  return(c(bottom = bottom_th, top = top_th))
}
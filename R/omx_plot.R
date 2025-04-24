library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(dendsort)
library(RColorBrewer)


# global variables --------------------------------------------------
colour.hm <- c(brewer.pal(9,'Blues')[9:3],brewer.pal(9,'Reds')[3:9])

# volcano --------------------------------------------------
# Plot a volcano plot
volcano <- function(
                df_de, 
                setting = 'add_titles',
                title = NULL,
                fig_dir = NULL,
                height=8,
                width=7){
    #--------------------------------------------------
    #' Plot volcano
    #'
    #' Ploting volcano for DEseq from deg function
    #'
    #' @param df_de dataframe from DEseq2 results
    #' @param setting file title of the plot
    #' @param title title of the plot
    #'
    #' @return dataframe with annotation of UP and DOWN genes
    #--------------------------------------------------    
    if (is.null(title)) {title <- setting}
    plot_vol <- ggplot(data=df_de, aes(x=log2FoldChange, y=-log10(padj), col=.data[['deg']])) + 
        geom_point() + 
        ggtitle(title)+
        theme_minimal()
    
    if (is.null(fig_dir)) {
            message("Not saving volcano figure...")
        } else {
            message("Saving volcano figure...")
            # Plot
            ggsave(paste0(fig_dir, setting, '_volcano.png'),
                plot=plot_vol, bg = "white", width=width, height=height, units="in", dpi= 600)
        }
  return(plot_vol)
}



#--------------------------------------------------
bar_plot <- function(
    data_frame,
    target_gene = 'PKP2',
    bar_fill = bar_color,
    plot_width = 8,
    plot_height = 6
) {
    #--------------------------------------------------
    #' Bar plot generation for protein abundance
    #'
    #' Creates a bar plot showing protein abundance across samples for a specific gene
    #'
    #' @param data_frame Input data frame containing protein expression data
    #' @param target_gene Gene symbol to plot (default: 'PKP2')
    #' @param bar_fill Color for the bars (default: royalblue)
    #' @param plot_width Width of the plot in inches (default: 8)
    #' @param plot_height Height of the plot in inches (default: 6)
    #'
    #' @return ggplot2 bar plot object
    #--------------------------------------------------
    
    # Extract data for the specified gene
    gene_data <- data_frame[which(data_frame$genes == target_gene), samples]
    gene_data <- as.data.frame(t(gene_data))
    colnames(gene_data) <- target_gene
    gene_data$samples <- factor(rownames(gene_data), levels = samples)
    
    # Generate the bar plot
    plot <- ggplot(gene_data, 
                  aes(x = samples, y = !!sym(target_gene))) +
        geom_bar(stat = 'identity', fill = bar_fill) +
        theme_minimal() +
        xlab('Samples') +
        ylab('Protein abundance') +
        ggtitle(target_gene)
    
    # Set plot dimensions
    ggsave(plot, width = plot_width, height = plot_height)
    
    return(plot)
}

#--------------------------------------------------
# Making heatmap
colour.hm <- c(brewer.pal(9,'Blues')[9:3],brewer.pal(9,'Reds')[3:9])
ht_opt$message = FALSE

heatmap_cm <- function(
                object,
                hm_setting = 'hm_setting',
                hm_title = NULL,
                top_anno_fig = NULL,
                left_anno_fig = NULL,
                ntop = 2000,
                colour_hm = colour.hm,
                fig_dir = NULL,
                width = 8,
                height = 18,
                cluster_columns = FALSE,
                cluster_rows = TRUE,
                show_column_names = FALSE,
                show_row_names = FALSE,
                kmeans = 1,
                #cluster = FALSE,
                use_raster = TRUE 
                ){
    
    #--------------------------------------------------
    #' heatmap plotting qq
    #'
    #' Return ComplexHeatmap analysis figures
    #'
    #' @param object matrix array with colnames for samples
    #' @param hm_setting naming the heatmap run
    #' @param hm_title naming  title when plotting
    #' @param top_anno_fig HeatmapAnnotation objects
    #' @param left_anno_fig HeatmapAnnotation objects
    #' @param ntop percentage of top variants
    #' @param colour_hm colour heatmap
    #' @param fig_dir saving directory for figures
    #' @param height height for heatmap plot
    #' @param width width for heatmap plot
    #' @param show_column_names columnnames
    #' @param show_row_names row names
    #' @param kmeans number of k-means clusters
    ##' @param cluster logical, whether to return row clusters
    #' @param use_raster suppress the message and create the heatmap
    #'
    #' @return heatmap object (and optionally row clusters if cluster = TRUE)
    #--------------------------------------------------

    #if (is.null(figure_dir)) stop("Please provide a colour heatmap.")
    if (is.null(hm_title)) {hm_title <- hm_setting}

    # check if nrop < nrow object
    if (nrow(object) < ntop) {
      ntop <- nrow(object)  # Keep the current value if it's less than 2000
      } else {
      ntop <- ntop  # Set ntop to 2000 if it's 2000 or more
      }
    # 
    if (!ht_opt$message) message('suppress message')


    row_var <- matrixStats::rowVars(object)   # calculate the variance for each gene
    top1_var <- order(row_var, decreasing=TRUE)[1:ntop]  # select the ntop genes by variance
    vsd_hm <- object[top1_var,]
    vsd_hm_scale <- t(scale(t(vsd_hm)))    # scale
    
    set.seed(920)
    # plot heatmap
    plot_vsd_hm <- Heatmap(vsd_hm_scale, 
                        name = paste(hm_setting, '_', nrow(vsd_hm_scale)),
                        col = colour_hm,
                        cluster_columns=cluster_columns,
                        cluster_rows = cluster_rows,
                        clustering_distance_rows = "spearman",
                        clustering_method_rows = "ward.D2",
                        #cluster_columns = dendsort(hclust(dist(t(vsd_hm_n)))),
                        show_column_names = show_column_names,
                        show_row_names = show_row_names,
                        row_names_side = "right",  # Show row names on the left
                        row_names_gp = gpar(fontsize = 5),  # Adjust the font size
                        show_row_dend = TRUE,
                        heatmap_legend_param = list(direction = "horizontal"),
                        km = kmeans,
                        #row_split = km,
                        #row_title = NULL,
                        top_annotation = top_anno_fig,
                        left_annotation = left_anno_fig,
                        heatmap_width = unit(width, "in"), 
                        heatmap_height = unit(height, "in")
                        )
                        
    plot_vsd_hm <- draw(plot_vsd_hm, merge_legend = TRUE, 
                        heatmap_legend_side = "right", 
                        annotation_legend_side = "bottom")
    # Save figures if fig_dir is provided
    if (!is.null(fig_dir)) {
        png_path <- paste0(fig_dir, hm_setting, '_heatmap_top_', as.character(ntop), '.png')
        pdf_path <- paste0(fig_dir, hm_setting, '_heatmap_top_', as.character(ntop), '.pdf')
        png(png_path, width = width + 4, height = height + 2, units = "in", res = 320)
        draw(plot_vsd_hm, merge_legend = TRUE, 
             heatmap_legend_side = "right", 
             annotation_legend_side = "bottom")
        dev.off()
        
        pdf(pdf_path, width + 4, height = height + 2)
        draw(plot_vsd_hm, merge_legend = TRUE, 
             heatmap_legend_side = "right", 
             annotation_legend_side = "bottom")
        dev.off()
    } else {
        message("Not saving figures")
    }
    
    # If cluster_rows = TRUE, compute and return row clusters as a table
    if (cluster_rows) {
        row_clusters <- row_order(plot_vsd_hm)
        row_names <- rownames(vsd_hm_scale)
        # Create a data frame: row names and their cluster assignments
        cluster_table <- data.frame(
            Row.names = unlist(lapply(row_clusters, function(x) row_names[x])),
            Cluster = rep(names(row_clusters), times = sapply(row_clusters, length)),
            stringsAsFactors = FALSE
        )
        return(list(heatmap = plot_vsd_hm, clusters = cluster_table))
    } else {
        return(plot_vsd_hm)
    }
}


# --------------------------------------------------
## Volcano --------------------------------------------------
# Predefined colors for differential expression
volcano_colors <- c("#FB6A4A", "grey", "#08519C")
# Define default volcano colors
volcano_colors <- c("#FB6A4A", "grey", "#08519C")  # Down-regulated, neutral, up-regulated

#' Volcano Plot for Differential Expression
#'
#' Creates a volcano plot of log2 fold change vs. -log10 adjusted p-value, 
#' with optional labels for the top percentage of genes by absolute log2 fold change,
#' connected by random-length sticks if stick_range is specified.
#'
#' @param df_de Data frame with differential expression results. Must contain 
#'   columns: `log2FoldChange`, `padj`, and the column specified by `col_deg`.
#' @param plot_title Character string for the plot title (default: "Differential Expression Volcano Plot").
#' @param col_deg Character string specifying the column name for DE status (default: "deg").
#' @param point_size Numeric, size of scatter points (default: 1.5).
#' @param plot_width Numeric, width of the plot in inches (default: 8).
#' @param plot_height Numeric, height of the plot in inches (default: 6).
#' @param alpha Numeric, transparency of points (default: 0.5).
#' @param colors Character vector of length 3 for colors (down-regulated, neutral, up-regulated; 
#'   default: `volcano_colors`).
#' @param label_top Logical, whether to label the top percentage of points by absolute log2 fold change 
#'   (default: TRUE).
#' @param label_col Character string, column name for labels (default: "SYMBOL" or "gene" if present).
#' @param stick_range Numeric vector of length 2, min and max stick length range (default: NULL, no sticks).
#' @param ntop Numeric, proportion of top points to label (default: 0.05, i.e., top 5%).
#' @param save_path Character string, file path to save the plot (e.g., "volcano.png"; default: NULL, no save).
#' @return A `ggplot2` object representing the volcano plot.
#' @examples
#' \dontrun{
#'   df <- data.frame(log2FoldChange = rnorm(100), padj = 10^-runif(100), 
#'                    deg = sample(c("Down", "Neutral", "Up"), 100, replace = TRUE))
#'   volcano_func(df, plot_title = "Basic Volcano", save_path = "volcano_basic.png")
#'   volcano_func(df, stick_range = c(0.5, 2), ntop = 0.1, save_path = "volcano_sticks.png")
#' }

volcano_plot <- function(
    df_de,
    plot_title = "Differential Expression Volcano Plot",
    col_deg = "deg",
    point_size = 1.5,
    plot_width = 8,
    plot_height = 6,
    alpha = 0.5,
    colors = volcano_colors,
    label_top = TRUE,
    label_col = NULL,
    stick_range = NULL,  # Optional, defaults to NULL
    ntop = 0.05,
    save_path = NULL,
    name = 'name'
) {
  
  # Load required library
  library(ggplot2)
  
  # Input validation
  if (!is.data.frame(df_de)) stop("Error: 'df_de' must be a data frame.")
  required_cols <- c("log2FoldChange", "padj", col_deg)
  missing_cols <- setdiff(required_cols, colnames(df_de))
  if (length(missing_cols) > 0) {
    stop("Error: Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  if (length(colors) != 3) stop("Error: 'colors' must be a vector of 3 colors.")
  if (!is.null(stick_range)) {
    if (length(stick_range) != 2 || !is.numeric(stick_range) || any(stick_range < 0)) {
      stop("Error: 'stick_range' must be a numeric vector of length 2 with non-negative values.")
    }
  }
  if (!is.numeric(ntop) || ntop <= 0 || ntop > 1) {
    stop("Error: 'ntop' must be a numeric value between 0 and 1.")
  }
  
  # Clean data for plotting
  df_de <- df_de[!is.na(df_de$log2FoldChange) & !is.infinite(df_de$log2FoldChange) &
                   !is.na(df_de$padj) & !is.infinite(df_de$padj), ]
  
  # Determine label column if not specified
  if (is.null(label_col)) {
    if ("SYMBOL" %in% colnames(df_de)) {
      label_col <- "SYMBOL"
    } else if ("gene" %in% colnames(df_de)) {
      label_col <- "gene"
    } else {
      warning("No 'SYMBOL' or 'gene' column found. No labels will be added.")
      label_top <- FALSE
    }
  }
  
  # Generate base volcano plot
  plot <- ggplot(data = df_de, 
                 aes(x = log2FoldChange, 
                     y = -log10(padj), 
                     col = .data[[col_deg]])) + 
    geom_point(size = point_size, alpha = alpha) + 
    ggtitle(plot_title) +
    labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
    theme_minimal() +
    scale_color_manual(values = colors, name = "diff") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      legend.position = "right"
    )
  
  # Add sticks and labels for top ntop% if requested and stick_range is provided
  if (label_top && !is.null(label_col)) {
    # Calculate absolute LFC for ranking
    df_de$abs_lfc <- abs(df_de$log2FoldChange)
    
    # Calculate top ntop% threshold
    n_top <- ceiling(nrow(df_de) * ntop)
    top_df <- df_de[order(df_de$abs_lfc, decreasing = TRUE), ][1:n_top, ]
    
    # Add sticks and labels only if stick_range is specified
    if (!is.null(stick_range)) {
      top_df$stick_length <- runif(n_top, min = stick_range[1], max = stick_range[2])
      top_df$y_label <- -log10(top_df$padj) + sign(top_df$log2FoldChange) * top_df$stick_length
      
      plot <- plot +
        geom_segment(data = top_df, 
                     aes(x = log2FoldChange, xend = log2FoldChange, 
                         y = -log10(padj), yend = y_label),
                     color = "black", size = 0.5) +  # Sticks with random lengths
        geom_text(data = top_df, 
                  aes(x = log2FoldChange, y = y_label, label = .data[[label_col]]),
                  size = 3, vjust = ifelse(top_df$log2FoldChange > 0, -0.5, 1.5),
                  check_overlap = TRUE,
                  show.legend = FALSE)
    } else {
      # If no stick_range, just add labels without sticks
      plot <- plot +
        geom_text(data = top_df, 
                  aes(label = .data[[label_col]]),
                  size = 3, vjust = -0.5, 
                  check_overlap = TRUE,
                  show.legend = FALSE)
    }
  }
  
  # Save plot if save_path is provided
  if (!is.null(save_path)) {
    figire_path <- paste0(save_path, name, '_volcano.png')
    ggsave(figire_path, plot, bg = 'white', width = plot_width, height = plot_height, dpi = 300)
    message("Plot saved to: ", save_path)
  }
  
  return(plot)
}

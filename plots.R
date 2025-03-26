library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(dendsort)
library(RColorBrewer)

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
                cluster = FALSE  # New parameter
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
    #' @param cluster logical, whether to return row clusters
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
    
    # If cluster = TRUE, compute and return row clusters as a table
    if (cluster) {
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
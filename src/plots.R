library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(dendsort)
library(RColorBrewer)

#--------------------------------------------------
bar_plot <- function(df, gene = 'PKP2'){
  df_gene <- df_protein[which(df_protein$genes == gene),samples]
  df_gene <- as.data.frame(t(df_gene))
  colnames(df_gene) <- gene
  df_gene$samples <- factor(rownames(df_gene), levels=samples)
  # plot 
  plot <- ggplot(df_gene, aes(x=samples, y=!!sym(gene))) +
    geom_bar(stat = 'identity', fill = 'royalblue')+
    theme_minimal() +
    xlab('Samples') +
    ylab('Protein abundance') +
    ggtitle(gene)
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
                ntop = 2000,
                colour_hm = colour.hm,
                fig_dir = NULL,
                width = 8,
                height = 18,
                cluster_columns = FALSE,
                cluster_rows = TRUE,
                show_column_names = FALSE,
                show_row_names = FALSE,
                kmeans = 1
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
    #' @param replicate sample replicates
    #' @param ntop percentage of top variants
    #' @param colour_hm colour heatmap
    #' @param fig_dir saving directory for figures
    #' @param colours colour for conditions
    #' @param height height for heatmap plot
    #' @param width width for heatmap plot
    #' @param show_column_names columnnames
    #' @param show_row_names row names
    #'
    #' @return heatmap objects
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
    

    # plot heatmap
    plot_vsd_hm <- Heatmap(vsd_hm_scale, name = paste("ATAC-seq - ", hm_setting),
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
                        heatmap_width = unit(width, "in"), 
                        heatmap_height = unit(height, "in")
                        )
    if (is.null(fig_dir)) {
        message("Not saving figures")
        } else {
        png_path <- paste0(figure_dir, hm_setting, '_heatmap_top_', as.character(ntop), '.png')
        pdf_path <- paste0(figure_dir, hm_setting, '_heatmap_top_', as.character(ntop), '.pdf')
        png(png_path, width=width+4, height=height+2, units="in", res=320)
        draw(plot_vsd_hm, merge_legend = TRUE, 
              heatmap_legend_side = "right", 
              annotation_legend_side = "bottom")
        dev.off()
        
        pdf(pdf_path, width+4, height=height+2)
        draw(plot_vsd_hm, merge_legend = TRUE, 
              heatmap_legend_side = "right", 
              annotation_legend_side = "bottom")
        dev.off()}
    
    return(plot_vsd_hm)
}
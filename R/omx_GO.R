##############################
# april - XQ - LUMC
# GO analysis
##############################


# go_analysis--------------------------------------------------
go_analysis <- function(setting, 
                           ont, 
                           DEup = NULL, 
                           DEdown = NULL, 
                           figure_dir = "figures/",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05,
                           go_id = TRUE) {

  if (is.null(DEup) && is.null(DEdown)) stop("Provide at least one UP or DOWN list")
  if (!is.null(DEup) && !is.character(DEup)) stop("Provide UP as character vector")
  if (!is.null(DEdown) && !is.character(DEdown)) stop("Provide DOWN as character vector")

  ontology <- c(MF = "Molecular Function",
                CC = "Cellular Component",
                BP = "Biological Process")
  
  message(paste0('Running ', ontology[[ont]], ' for: ', setting))
  
  # Initialize empty data frames
  go_data_up <- data.frame()
  go_data_down <- data.frame()
  
  # Perform GO enrichment for upregulated genes if provided
  if (!is.null(DEup)) {
    goDEMF_up <- enrichGO(DEup,
                       org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = ont,
                       pvalueCutoff = pvalueCutoff,
                       pAdjustMethod = "fdr",
                       qvalueCutoff = qvalueCutoff)
  if (is.null(goDEMF_up) || nrow(as.data.frame(goDEMF_up)) == 0) {
    message("No UP GO terms found.")
  } else {
    go_data_up <- as.data.frame(goDEMF_up)
    go_data_up$GeneRatio <- sapply(go_data_up$GeneRatio, function(x) eval(parse(text = x)))
    go_data_up$Regulation <- "Upregulated"
  }
  }
  
  # Perform GO enrichment for downregulated genes if provided
  if (!is.null(DEdown)) {
    goDEMF_down <- enrichGO(DEdown,
                            org.Hs.eg.db,
                            keyType = "SYMBOL",
                            ont = ont,
                            pvalueCutoff = pvalueCutoff,
                            pAdjustMethod = "fdr",
                            qvalueCutoff = qvalueCutoff)
    if (is.null(goDEMF_down || nrow(as.data.frame(goDEMF_down)) == 0)) {
    message("No DOWN GO terms found.")
  } else {
    go_data_down <- as.data.frame(goDEMF_down)
    go_data_down$GeneRatio <- sapply(go_data_down$GeneRatio, function(x) eval(parse(text = x)))
    go_data_down$Regulation <- "Downregulated"
    go_data_down$GeneRatio <- -go_data_down$GeneRatio  # Negative for leftward bars
  }
  }
    
    # Check if we have any data
  if (nrow(go_data_up) == 0 && nrow(go_data_down) == 0) {
    stop("No GO term found. Considering increase the pvalueCutoff and pvalueCutoff")
  }
  
  # Combine the top 10 from each dataset (if they exist)
  combined_data <- data.frame()
  if (nrow(go_data_up) > 0) {
    combined_data <- rbind(combined_data, go_data_up[1:min(10, nrow(go_data_up)), ])
  }
  if (nrow(go_data_down) > 0) {
    combined_data <- rbind(combined_data, go_data_down[1:min(10, nrow(go_data_down)), ])
  }
  combined_data <- combined_data %>% filter(Count != "NA")
  
  # add go id:
  if (go_id) {
    combined_data$Description <- paste0(combined_data$Description, ' (', combined_data$ID, ')')
  }

  # Function to format scientific notation
  format_scientific_expr <- function(x) {
  sci <- format(x, scientific = TRUE)
  coef <- as.numeric(sub("e.*", "", sci))
  exp <- as.numeric(sub(".*e", "", sci))
  parse(text = sprintf("%g %%*%% 10^%d", coef, exp))
  }

  # Create the plot
  plot_goterm <- ggplot(combined_data, aes(x = GeneRatio, y = reorder(Description, GeneRatio), fill = p.adjust)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Count, hjust = ifelse(Regulation == "Upregulated", -0.2, 1.2)), 
              size = 3.5, color = "black") +
    scale_fill_gradientn(
      colors = c('#08519C', '#9ECAE1'),
      name = "Adjusted P-value",
      labels = function(x) sapply(x, format_scientific_expr),
      guide = guide_colorbar(
        label.hjust = 0.5,
        label.theme = element_text(size = 10)
      )
    ) +

    labs(title = paste0(setting, ": GO-DAR ", ontology[[ont]]),
         x = "Gene Ratio",
         y = "GO Terms") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10, face = "bold"),
          axis.text.x = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
    coord_cartesian(xlim = c(min(combined_data$GeneRatio) * .75, max(combined_data$GeneRatio) * 1.25)) +
    geom_vline(xintercept = 0, color = "black", linetype = "longdash")
  
  # Add annotations only if both up and down data are present
  if (!is.null(DEup) && !is.null(DEdown)) {
    plot_goterm <- plot_goterm +
      annotate("rect", xmin = -0.06, xmax = -0.005, ymin = -0.75, ymax = -0.3, fill = "#F48D79") +
      annotate("rect", xmin = 0.06, xmax = 0.005, ymin = -0.75, ymax = -0.3, fill = "#017075") +
      annotate("text", x = -0.02, y = -0.1, label = "DOWN", size = 3.5) +
      annotate("text", x = 0.02, y = -0.1, label = "UP", size = 3.5)
  }
  
  # Save the plot with pv and qv in filename
  ggsave(paste0(figure_dir, setting, '_', ont, '_pv_', pvalueCutoff, '_qv_', qvalueCutoff, '_go_term.png'),
         plot = plot_goterm, 
         bg = 'white',
         width = 10,
         height = 7,
         units = "in",
         dpi = 320)
  
  # Return the plot object
  message(paste0('Finish ', ontology[[ont]], ' for: ', setting))
  return(plot_goterm)
}





# perform enrichment (old) --------------------------------------------------
# note: this function is redundence, consider removed
# make GO term
#' Perform GO Enrichment Analysis and Visualization
#'
#' @param deg_data DEG data frame containing differential expression results
#' @param layer Character string specifying the layer name
#' @param direction Character string specifying gene direction ('UP' or 'DOWN', default: 'UP')
#' @param lfc Log fold change cutoff (default: 1)
#' @param padj P-value adjustment cutoff (default: 0.05)
#' @param pvalueCutoff P-value cutoff for enrichment (default: 0.05)
#' @param qvalueCutoff Q-value cutoff for enrichment (default: 0.05)
#' @param figure_dir Directory path for saving the output plot
#' @return A grid plot object and saves PNG file
#' @export
perform_go_enrichment <- function(deg_data,
                                 layer,
                                 direction = "UP",
                                 lfc = 1,
                                 padj = 0.05,
                                 pvalueCutoff = 0.05,
                                 qvalueCutoff = 0.05,
                                 width = 20,
                                 height = 14,
                                 figure_dir = "./") {
  
  # Required packages
  require(clusterProfiler)
  require(org.Hs.eg.db)
  require(ggplot2)
  require(cowplot)
  
  # Validate direction parameter
  direction <- toupper(direction)
  if (!direction %in% c("UP", "DOWN")) {
    stop("direction must be either 'UP' or 'DOWN'")
  }
  
  message("Starting GO enrichment analysis...")
  
  # Calculate DEG
  message("Calculating differential expression...")
  deg_results <- deg(deg_data, lfc, padj)
  
  # Extract genes based on direction
  message(sprintf("Extracting %s-regulated genes...", direction))
  DEgenes <- deg_results[deg_results$deg == direction, 'SYMBOL']
  
  # Perform GO enrichment for all three ontologies
  message("Performing Molecular Function enrichment...")
  goDEMF <- enrichGO(DEgenes,
                     org.Hs.eg.db,
                     keyType = "SYMBOL",
                     ont = "MF",
                     pvalueCutoff = pvalueCutoff,
                     pAdjustMethod = "fdr",
                     qvalueCutoff = qvalueCutoff)
  
  message("Performing Cellular Component enrichment...")
  goDECC <- enrichGO(DEgenes,
                     org.Hs.eg.db,
                     keyType = "SYMBOL",
                     ont = "CC",
                     pvalueCutoff = pvalueCutoff,
                     pAdjustMethod = "fdr",
                     qvalueCutoff = qvalueCutoff)
  
  message("Performing Biological Process enrichment...")
  goDEBP <- enrichGO(DEgenes,
                     org.Hs.eg.db,
                     keyType = "SYMBOL",
                     ont = "BP",
                     pvalueCutoff = pvalueCutoff,
                     pAdjustMethod = "fdr",
                     qvalueCutoff = qvalueCutoff)
  
  # Create visualizations
  message("Generating plots...")
  goDEMFplot <- dotplot(goDEMF, showCategory=10) + 
    ggtitle(paste0(direction, " Genes in ", layer, "\n Molecular Function - FDR < ", pvalueCutoff))
  
  goDECCplot <- dotplot(goDECC, showCategory=10) + 
    ggtitle(paste0(direction, " Genes in ", layer, "\n Cellular Component - FDR < ", pvalueCutoff))
  
  goDEBPplot <- dotplot(goDEBP, showCategory=10) + 
    ggtitle(paste0(direction, " Genes in ", layer, "\n Biological Process - FDR < ", pvalueCutoff))
  
  # Create combined plot
  message("Combining plots...")
  title <- paste0("Over Representation Analysis of ", layer, " (", direction, " genes)")
  DElayerGO <- plot_grid(goDEMFplot,
                        goDECCplot,
                        goDEBPplot,
                        labels = title,
                        ncol = 3)
  
  # Save plot
  message("Saving plot to file...")
  output_file <- paste0(figure_dir, 'go_', layer, "_", tolower(direction), ".png")
  ggsave(output_file,
         plot = DElayerGO,
         width = width,
         height = height,
         units = "in",
         dpi = 320)
  
  message("Analysis complete!")
  
  # Return the plot object
  return(DElayerGO)
}

# Example usage:
# For upregulated genes:
# result_up <- perform_go_enrichment(
#   deg_data = dar_deg_w3_w1,
#   layer = "w3_e2",
#   direction = "UP",
#   figure_dir = "/path/to/figures/"
# )

# For downregulated genes:
# result_down <- perform_go_enrichment(
#   deg_data = dar_deg_w3_w1,
#   layer = "w3_e2",
#   direction = "DOWN",
#   figure_dir = "/path/to/figures/"
# )

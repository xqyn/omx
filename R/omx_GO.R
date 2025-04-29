# project omx: GO analysis 
# XQ - LUMC
# april 1 2025
# update:
#   - april 25 2025: 
#'     - add maths
#'     - split kind-func
#'     - add warning + control message for step
#'     - new colourset + colour control coordination
#'     - happy Koningsdag
#'  - april 29 2025:
#      - make_fill_scale


#' @title Global ontology vector
#' @description vector GO ontologies for mapping
ontology <- c(MF = "Molecular Function",
              CC = "Cellular Component",
              BP = "Biological Process")


# --- format_scientific_expr
#' Format Scientific Notation for Plot Labels
#'
#' Converts a numeric value into a formatted scientific notation expression suitable for plot labels.
#'
#' @param x Numeric value to be formatted.
#'
#' @return A parsed expression representing the value in scientific notation (e.g., 1.23 * 10^4).
#'
#' @examples
#' \dontrun{
#'   format_scientific_expr(0.000123)  # Returns expression for 1.23 * 10^-4
#' }

format_scientific_expr <- function(x) {
  sci <- format(x, scientific = TRUE)
  coef <- as.numeric(sub("e.*", "", sci))
  exp <- as.numeric(sub(".*e", "", sci))
  parse(text = sprintf("%g %%*%% 10^%d", coef, exp))
}


# --- go_enrichment
#' Process GO Enrichment Analysis
#'
#' Performs Gene Ontology (GO) enrichment analysis for a given set of genes 
#'  and processes the results.
#'
#' @param genes character vector of gene symbols.
#' @param regulation character string indicating regulation direction ("Upregulated" or "Downregulated").
#' @param ont character string specifying the GO ontology ("MF", "CC", or "BP").
#' @param pvalueCutoff numeric value for p-value cutoff for GO enrichment. (default: 0.05).
#' @param qvalueCutoff numeric value for the q-value cutoff for GO enrichment. (default: 0.05).
#'
#' @return A data frame with processed GO enrichment results or an empty data frame if no terms are found.
#'
#' @importFrom clusterProfiler enrichGO
#' @importFrom org.Hs.eg.db org.Hs.eg.db

go_enrichment <- function(genes = NULL,
                          regulation = "Upregulated", 
                          ont = "MF",
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.05) {

  if (is.null(genes)) return(data.frame())
  
  message(paste0("GO ", regulation, " for ", ontology[[ont]], "..."))

  go_result <- clusterProfiler::enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = ont,
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = "fdr",
    qvalueCutoff = qvalueCutoff
  )
  
  if (is.null(go_result) || nrow(as.data.frame(go_result)) == 0) {
    message(paste0("No ", regulation, " GO terms found for: ", ontology[[ont]], "!."))
    return(data.frame())
  }
  
  go_data <- as.data.frame(go_result)
  go_data$GeneRatio <- sapply(go_data$GeneRatio, function(x) eval(parse(text = x)))
  go_data$Regulation <- regulation
  if (regulation == "Downregulated") {
    go_data$GeneRatio <- -go_data$GeneRatio  # Negative for leftward bars
  }
  return(go_data)
}


# --- combine_go_data
#' Combine GO Enrichment Data
#'
#' Combines the top GO terms from up- and down-regulated gene datasets, filters out invalid counts,
#' and optionally appends GO IDs to term descriptions.
#'
#' @param go_up dataframe containing GO enrichment results for up-regulated genes. (from go_enrichment)
#' @param go_down dataframe containing GO enrichment results for down-regulated genes. (from go_enrichment)
#' @param ont character string specifying the GO ontology ("MF", "CC", or "BP").
#' @param top top GO-term to be plotted. (default: 10).
#' @param go_id nogical indicating whether to append GO IDs to term descriptions in the plot.
#'
#' @return A data frame combining the top GO terms from both datasets, 
#'         with filtered counts and optional GO IDs.

combine_go_data <- function(go_up = NULL,
                            go_down = NULL, 
                            ont = 'MF',
                            top = 10, 
                            go_id = TRUE) {
  
  # if input is NULL, replace as empty data frame
  if (is.null(go_up)) go_up <- data.frame()
  if (is.null(go_down)) go_down <- data.frame()

  # if both is empty, return empty
  if (nrow(go_up) == 0 && nrow(go_down) == 0) {
    warning(paste0("No GO term found in ", ontology[[ont]], " for both GO up and down."))
    message("Consider increasing the pvalueCutoff and qvalueCutoff")
    return(data.frame())
  }
  
  go_combine <- data.frame()
  if (nrow(go_up) > 0) {
    go_combine <- rbind(go_combine, go_up[1:min(top, nrow(go_up)), ])}
  if (nrow(go_down) > 0) {
    go_combine <- rbind(go_combine, go_down[1:min(top, nrow(go_down)), ])}
  
  go_combine <- go_combine[!is.na(go_combine$Count), ]
  
  if (go_id) {
    go_combine$Description <- paste0(go_combine$Description, ' (', go_combine$ID, ')')}
  
  return(go_combine)
}


# ---make_fill_scale
#' Create a ggplot2 fill scale with log10-transformed p-values and optional custom breaks
#'
#' This function generates a `scale_fill_gradientn()` object with a log10 color scale,
#' customized labels, limits, and colorbar appearance for visualizing adjusted p-values.
#' If a valid `breaks` vector is provided (length ≥ 2), it is added to the scale;
#' otherwise, breaks are omitted and ggplot2 defaults are used.
#'
#' @param breaks numeric vector of break points for the fill scale (e.g., log10 p-value thresholds).
#'               If NULL or fewer than two elements, no custom breaks are applied.
#' @param setting character string for plot title and file naming (default: "go_analysis").
#' @param ontology_name string for ontology (e.g., "BP", "MF", or "CC") shown in the legend title.
#' @param pvals numeric vector of p-values (typically `go_combine$p.adjust`) used to set scale limits.
#' @param pvalueCutoff numeric p-value cutoff for file naming (default: 0.05).
#' @param qvalueCutoff numeric q-value cutoff for file naming (default: 0.05).

#' @return A `ggplot2` scale object (`scale_fill_gradientn`) that can be added to a plot.
#'
#' @examples
#' # ased inside a ggplot2 plot
#' ggplot(data, aes(x, y, fill = p.adjust)) +
#'   geom_bar(stat = "identity") +
#'   make_fill_scale(breaks, "Treatment A", "BP", data$p.adjust)

make_fill_scale <- function(breaks, setting, ontology_name, pvals, 
                            pvalueCutoff = 0.05, qvalueCutoff = 0.05) {
  base_args <- list(
    colors = c('#05474A', '#8CD0D5'),
    name = paste0(
      setting, ": ", "\n",
      ontology_name, "\n",
      "(p_value < ", pvalueCutoff, "\n",
      " q_value < ", qvalueCutoff, ")", "\n",
      "\n",
      "\n",
      "adj P-value"),
    labels = function(x) sapply(x, format_scientific_expr),
    trans = "log10",
    limits = range(pvals, na.rm = TRUE),
    guide = ggplot2::guide_colorbar(
      label.hjust = 0.5,
      label.theme = ggplot2::element_text(size = 8),
      barheight = 8,
      nbin = 100)
  )

  # add breaks only if valid
  if (!is.null(breaks) && length(breaks) >= 2) {
    base_args$breaks <- breaks
  }

  do.call(ggplot2::scale_fill_gradientn, base_args)
}


# --- plot_go_terms
#' Plot GO Terms as a Bar Plot
#'
#' Creates a bar plot to visualize the top Gene Ontology (GO) terms from enrichment analysis.
#'
#' @param go_combine dataframe with GO enrichment results, including columns for GeneRatio, Description, p.adjust, Count, and Regulation.
#' @param setting character string for plot title and file naming (default: "go_analysis").
#' @param ont character string specifying GO ontology: "MF" (Molecular Function), "CC" (Cellular Component), or "BP" (Biological Process) (default: "MF").
#' @param pvalueCutoff numeric p-value cutoff for file naming (default: 0.05).
#' @param qvalueCutoff numeric q-value cutoff for file naming (default: 0.05).
#' @param GeneRatio character string specifying the column name for GeneRatio values (default: "GeneRatio").
#' @param figure_dir directory to save the plot; if NULL, the plot is not saved (default: NULL).
#' @param width plot width in inches (default: 10).
#' @param height plot height in inches (default: 7).
#' @param wrap_width numeric value specifying the maximum character width for wrapping y-axis labels (default: 50).
#' @param fig_format save the figure as which format
#'
#' @return A ggplot2 object representing the GO enrichment results as a bar plot.
#'
#' @importFrom ggplot2 ggplot aes geom_bar geom_text scale_fill_gradientn labs theme_minimal
#' @importFrom ggplot2 theme coord_cartesian geom_vline annotate guide_colorbar element_text
#' @importFrom ggplot2 ggsave
#'
#' @examples
#' \dontrun{
#' plot_go_terms(go_combine, setting = "go_analysis", ont = "BP", figure_dir = "./plots/")
#' }

plot_go_terms <- function(go_combine,
                          setting = 'go_analysis',
                          ont = 'MF',
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05,
                          GeneRatio = "GeneRatio",
                          figure_dir = NULL,
                          width = 10,
                          height = 7,
                          wrap_width = 50,
                          fig_format = c("png", "pdf")) {
  
  # check required columns in go_combine
  required_columns <- c("Description", "p.adjust", "Count", "Regulation")
  missing_columns <- setdiff(required_columns, colnames(go_combine))
  if (length(missing_columns) > 0) {
    stop("go_combine is missing required columns: ", paste(missing_columns, collapse = ", "))
  }
  
  # calculate max and min for GeneRatio x-axis
  max_gene_ratio <- max(go_combine[[GeneRatio]])
  min_gene_ratio <- min(go_combine[[GeneRatio]])
  max_gene_ratio <- ceiling(max_gene_ratio * 100) / 100 # Round up to nearest 0.01
  min_gene_ratio <- floor(min_gene_ratio * 100) / 100   # Round down to nearest 0.01
  
  # setting range for x-axis for GeneRatio + remove 0 for visually favor
  # adjust min_gene_ratio if both are positive or negative
  if (min_gene_ratio > 0 & max_gene_ratio > 0) {
    min_gene_ratio <- 0.01
  } else if (min_gene_ratio < 0 & max_gene_ratio < 0) {
    max_gene_ratio <- -0.01
  }
  # generate the x-axis range
  range_x_axis <- seq(min_gene_ratio, max_gene_ratio, by = 0.01)
  range_x_axis <- range_x_axis[range_x_axis != 0]

  # setting breaks for p-value legend
  min_p_adjust <- min(go_combine[["p.adjust"]])
  min_exponent <- floor(log10(min_p_adjust) / 2) * 2
  breaks_exponents <- seq(min_exponent, -2, by = 2)
  breaks <- 10^breaks_exponents
  
  # plot
  plot_goterm <- ggplot2::ggplot(
    go_combine, 
    ggplot2::aes(
      x = .data[[GeneRatio]], 
      y = reorder(stringr::str_wrap(Description, width = wrap_width), 
                  .data[[GeneRatio]]),
      fill = p.adjust)) +
    ggplot2::geom_bar(stat = "identity") +
    # text for count on top of bars
    ggplot2::geom_text(
      ggplot2::aes(label = Count),
      hjust = ifelse(go_combine$Regulation == "Upregulated", -0.2, 1.2),
      size = 3.5, 
      color = "black") +
    # make_fill function (ggplot objects)
    make_fill_scale(
      breaks = breaks,
      setting = setting,
      ontology_name = ontology[[ont]],
      pvalueCutoff = pvalueCutoff,
      qvalueCutoff = qvalueCutoff,
      pvals = go_combine$p.adjust
    ) +
    # setting scale for x-axis for GeneRatio
    ggplot2::coord_cartesian(
      xlim = c(if (min(go_combine[[GeneRatio]]) < 0) min(go_combine[[GeneRatio]]) * 1.25 else 0,
               if (max(go_combine[[GeneRatio]]) > 0) max(go_combine[[GeneRatio]]) * 1.25 else 0)) +
    ggplot2::scale_x_continuous(labels = abs, breaks = range_x_axis) +

    # coordinate information
    ggplot2::geom_vline(xintercept = 0, color = "black", linetype = "longdash") +
    ggplot2::labs(#title = paste0(setting, ": ", ontology[[ont]]),
                  x = "Gene Ratio",
                  y = "GO Terms") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 10, angle = 45, vjust = 0.5),
      axis.ticks.x = ggplot2::element_line(linewidth = 0.5),
      axis.ticks.length.x = grid::unit(0.2, "cm"),
      axis.text.y = ggplot2::element_text(size = 10, face = "bold"),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(10, 10, 10, 20, "mm"),
      legend.position = "right")
    
  # add annotation bars if GeneRatio contains both positive and negative values
  if (any(go_combine[[GeneRatio]] > 0) && any(go_combine[[GeneRatio]] < 0)) {
    plot_goterm <- plot_goterm +
      ggplot2::annotate("rect", xmin = min_gene_ratio - 0.01, xmax = -0.0025, 
                        ymin = -.85, ymax = -.15, fill = "#F48D79") +
      ggplot2::annotate("rect", xmin = max_gene_ratio + 0.01, xmax = 0.0025, 
                        ymin = -.85, ymax = -.15, fill = "#08519C") +
      ggplot2::annotate("text", label = "UP", x = 0.02, y = - .5, size = 3, color = "white", fontface = "bold") +
      ggplot2::annotate("text", label = "DOWN", x = - 0.02, y = -.5, size = 3, color = "white", fontface = "bold")
  }

  # save the plot with pv and qv in filename if figure_dir is provided
  # fig_format <- match.arg(fig_format)
  if (!is.null(figure_dir)) {
    for (ext in fig_format) {
      filepath <- paste0(figure_dir, setting, '_', ont, 
                        '_pv_', pvalueCutoff, 
                        '_qv_', qvalueCutoff,
                        '_go_term.', ext)
      ggplot2::ggsave(filepath,
                      plot = plot_goterm, 
                      bg = 'white',
                      width = width,
                      height = height,
                      units = "in",
                      dpi = 320)
      message(paste0('Saved figure to: ', filepath))
    }
  }
  return(plot_goterm)
}

# --- go_analysis
#' Gene Ontology Enrichment Analysis and Visualization
#'
#' Performs Gene Ontology (GO) enrichment analysis for up- and/or down-regulated genes
#' and generates a bar plot visualizing the top GO terms.
#'
#' @param setting character string for plot title and file naming (default: "go_analysis").
#' @param ont character string specifying GO ontology to analyze. Must be one of:
#'   "MF" (Molecular Function), "CC" (Cellular Component), or "BP" (Biological Process).
#' @param DEup character vector of gene symbols for up-regulated genes. (default: NULL).
#' @param DEdown character vector of gene symbols for down-regulated genes. (default: NULL).
#' @param pvalueCutoff numeric value for p-value cutoff for GO enrichment. (default: 0.05).
#' @param qvalueCutoff numeric value for the q-value cutoff for GO enrichment. (default: 0.05).
#' @param top top GO-term to be plotted. (default: 10).
#' @param go_id logical indicating whether to append GO IDs to term descriptions in the plot.
#' @param figure_dir directory to save the output plot. (default: "figures/").
#' @param width Plot width in inches (default: 10).
#' @param height Plot height in inches (default: 7).
#'
#' @return A ggplot2 plot object representing the GO enrichment results as a bar plot.
#'
#' @details
#' This function uses the `clusterProfiler` package to perform GO enrichment analysis
#' on provided lists of up- and/or down-regulated genes. It generates a bar plot
#' showing the top GO terms for each regulation direction, with GeneRatio on the
#' x-axis and GO terms on the y-axis. The plot is saved as a PNG file in the specified
#' directory, and the plot object is returned.
#'
#' @examples
#' \dontrun{
#'   DEup <- c("GENE1", "GENE2", "GENE3")
#'   DEdown <- c("GENE4", "GENE5")
#'   go_plot <- go_analysis(
#'     setting = "Condition1",
#'     ont = "BP",
#'     DEup = DEup,
#'     DEdown = DEdown,
#'     figure_dir = "output/",
#'     pvalueCutoff = 0.05,
#'     qvalueCutoff = 0.05,
#'     top = 10,
#'     go_id = TRUE
#'   )
#' }
#'
#' @importFrom clusterProfiler enrichGO
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom ggplot2 ggplot aes geom_bar geom_text scale_fill_gradientn labs theme_minimal
#' @importFrom ggplot2 theme coord_cartesian geom_vline annotate guide_colorbar element_text
#' @importFrom ggplot2 ggsave
#' @importFrom dplyr rbind filter
#' @export

go_analysis <- function(setting = 'go_analysis', 
                        ont = 'MF', 
                        DEup = NULL, 
                        DEdown = NULL, 
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        top = 10,
                        go_id = TRUE,
                        figure_dir = "figures/",
                        width = 10,
                        height = 7,
                        wrap_width = 50,
                        fig_format = c("png", "pdf")) {

  if (is.null(DEup) && is.null(DEdown)) stop("Provide at least one UP or DOWN list")
  if (!is.null(DEup) && !is.character(DEup)) stop("Provide UP as character vector")
  if (!is.null(DEdown) && !is.character(DEdown)) stop("Provide DOWN as character vector")
  
  message("")
  message(paste0('Running ', ontology[[ont]], ' for: ', setting))
  
  # initialize empty data frames
  go_up <- data.frame()
  go_down <- data.frame()
  
  # perform GO enrichment for upregulated genes:
  go_up <- go_enrichment(genes = DEup,
                        regulation = "Upregulated",
                        ont = ont, 
                        pvalueCutoff = pvalueCutoff, 
                        qvalueCutoff = qvalueCutoff)
  
  go_down <- go_enrichment(genes = DEdown,
                        regulation = "Downregulated",
                        ont = ont,
                        pvalueCutoff = pvalueCutoff,
                        qvalueCutoff = qvalueCutoff)

  # combine up and down go data 
  go_combine <- combine_go_data(
                        go_up = go_up, 
                        go_down = go_down,
                        ont = ont, 
                        top = top, 
                        go_id = go_id)
  
  # check if combined_data is empty
  if (nrow(go_combine) == 0) {
    message(paste0("NO valid GO terms were found ", ontology[[ont]], ". Skipping plot"))
    return(NULL)
  }

  # generate and save the plot
  plot_goterm <- plot_go_terms(
    go_combine = go_combine,
    setting = setting,
    ont = ont,
    pvalueCutoff = pvalueCutoff,
    qvalueCutoff = qvalueCutoff,
    figure_dir = figure_dir,
    width = width,
    height = height,
    fig_format = fig_format
  )
  
  # return the plot object
  message(paste0('Finish ', ontology[[ont]], ' for: ', setting))
  message("----------")
  return(plot_goterm)
}
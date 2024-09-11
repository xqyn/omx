# Custumize function for DEseq2 and PCA
library(ggplot2)


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

    return(df_de)
}

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

colours <- c("#1D263B","#E03E3E","#53BB79","#EF8228","#937264",
            "#3043A2","#C25D7B","#D88C4C","dodgerblue","#FEB447",
            "#06CDE0","#1498BE","#0D5B11", '#4c00b0', '#ffaf00')

# PCA analysis
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
                  colour_set = colours,
                  height = 6,
                  width = 7.5){

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
    
    ## fratures contribute to PC
    attri <- as.data.frame(res.pca$var$contrib)
    colnames(attri) <- sub('Dim.','pc', colnames(attri))
    attri <- attri[order(attri$pc1, decreasing=TRUE),] # sort on pc1
    
    if (is.null(fig_dir)) {
        message("Not saving figures")
    } else {
        message("Ploting PCA figures...")
        # Plot
        ## PCA plot
        plot_pca <- ggplot(pca, aes(x = Dim.1, y = Dim.2)) +
                        geom_point(aes(shape=.data[[rep_title]],
                                    colour = .data[[cond_title]]), 
                                    size = 5, 
                                    alpha=0.75, 
                                    show.legend = TRUE) +
                        labs(x = paste0("PC1: ", var$per_var[1], "%"),
                            y = paste0("PC2: ", var$per_var[2], "%"),
                            color = cond_title,
                            shape = rep_title
                            ) +
                        scale_color_manual(values = colour_set) + 
                        ggtitle(paste0('PCA: ', pca_setting)) +
                        theme_minimal()
        ## save
        ggsave(paste0(fig_dir, pca_setting, '_pca.png'),
            plot=plot_pca, bg = 'white', width=width, height=height, units="in", dpi=600)

        ggsave(paste0(fig_dir, pca_setting, '_pca.pdf'),
            plot=plot_pca, bg = "white", width=width, height=height, units="in", dpi= 1200, device = cairo_pdf)

        ## Varaince plot
        plot_var <- ggplot(var, aes(x = factor(1:nrow(var)), y = per_var)) +
                        geom_bar(stat = "identity") +
                        labs(x = "PC", 
                            y = "Percent variance", 
                            title = paste0("Percent varience explained: ", pca_setting)) +
                        theme_minimal()
        ## save
        ggsave(paste0(fig_dir, pca_setting, '_pca_var.png'), 
            plot=plot_var, bg = 'white', width=6.5, height=6, units="in", dpi=600)
        
        ggsave(paste0(fig_dir, pca_setting, '_pca_var.pdf'), 
            plot=plot_var, bg = "white", width=6.5, height=6,units="in", dpi= 1200, device = cairo_pdf)
    }
    
    return(list(pca = pca, var = var, attri = attri))  
  }
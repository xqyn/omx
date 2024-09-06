# Custumize function for DEseq2 and PCA
library(ggplot2)

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
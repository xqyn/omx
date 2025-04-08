

# GAMLSS function --------------------------------------------------
# library(ggplot2)
# library(dplyr)
# library(MASS)
# library(viridis)
# library(reshape2)
# library(mgcv)  # for gam and shash

plot_quantile_model <- function(df, title) {
  # Helper functions
  params_to_quantiles_shash <- function(quantiles, params, qshash){
    as.data.frame(sapply(quantiles, 
                         function(q){
                           qshash(p=q, mu=params)
                         }))
  }
  reshape_quantiles_to_long <- function(quantiles_df, x_var){
    quantiles_df$x <- x_var
    return(reshape2::melt(quantiles_df, id.vars = c('x')))
  }
  
  get_density <- function(x, y, n = 100) {
    dens <- MASS::kde2d(x = x, y = y, n = n)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
  }
  

  
  # Step 1: Fit the GAMLSS model
  colnames(df) <- c('x', 'y')
  df$x <- as.numeric(as.character(df$x))
  
  message("Fitting GAMLSS model...")
  model <- gam(list(y ~ s(x),
                    ~ s(x), 
                    ~ 1,
                    ~ 1),
               family = shash(),
               data = df)
  
  # Step 2: Generate predictions
  message("Generating model predictions...")
  predictions_params_model <- predict(model, newdata = df)
  quantiles <- pnorm(c(-2:2))
  qshash <- model$family$qf
  
  # Step 3: Compute quantiles
  message("Computing quantiles from predictions...")
  predictions_quantiles_model <- params_to_quantiles_shash(quantiles, 
                                                           predictions_params_model,
                                                           qshash)
  
  # Step 4: Reshape predictions to long format
  message("Reshaping quantile data to long format...")
  predictions_quantiles_model_long <- reshape_quantiles_to_long(
    predictions_quantiles_model, df$x)
  
  # Step 5: Calculate density
  message("Calculating density for contour plotting...")
  k <- 100  # Fixed value for density grid
  x <- df$x
  y <- df$y
  contour_cols <- viridis(k, alpha = 0.5)
  dens <- get_density(x, y, k)
  
  # Step 6: Assign quantile levels
  message("Assigning quantile levels...")
  quantile_full <- c('75th', '50th', 'mean', '50th', '75th')
  quantile_levels <- unique(quantile_full)
  predictions_quantiles_model_long$quantile <- quantile_full[predictions_quantiles_model_long$variable]
  predictions_quantiles_model_long$quantile <- factor(predictions_quantiles_model_long$quantile, 
                                                      levels = quantile_levels)
  
  # Step 7: Create the plot
  message("Generating the plot...")
  p_xq <- ggplot(df) +
    geom_jitter(aes(x = x, y = y), 
                width = 0.75, height = 0.75, size = 1, 
                col = contour_cols[findInterval(dens, seq(0, max(dens), length.out = k))], 
                pch = 16, 
                alpha = 0.25) +
    geom_line(data = predictions_quantiles_model_long, 
              aes(x = x, 
                  y = value, 
                  group = variable,
                  linetype = quantile,    # Map linetype to quantile
                  color = quantile,       # Map color to quantile
                  size = quantile)) +     # Map size to quantile
    scale_linetype_manual(values = c("dotted", "longdash", "dashed", "longdash", "dotted")) +
    scale_color_manual(values = c('#ff7b7b', '#ff5252', '#a70000', '#ff5252', '#ff7b7b')) +
    scale_size_manual(values = c(0.5, 0.35, 0.5, 0.35, 0.5)) +
    labs(title = title,
         subtitle = expression(paste('GAMLSS | y ~ D(',mu,',',sigma,',',nu,',',tau,')')),
         x = 'day',
         y = 'Log10(readcount)',
         linetype = "Quantile",   # Legend title for linetype
         color = "Quantile",      # Legend title for color
         size = "Quantile") +     # Legend title for size
    theme_minimal()
  
  message("Plot generation complete!")
  return(p_xq)
}

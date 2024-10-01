
plot_permanova <- function(ps_obj, formula, transformation = "compositional", method = "bray", show_plot = TRUE,level,det,prev,taxa_num,year,lab.size=lab.size,size=15,cols_meta,
                           palette) {
  
  # Convert sample data to a data frame
  metadf <- as.data.frame(as.matrix(ps_obj@sam_data))
  
  if (method == 'euclidean'){
    transformation <- 'clr'
    ps_obj <- microbiome::transform(ps_obj,transformation)
  }
  if (method == 'uunifrac'){
    transformation <- '-'
  }
  
  if ((transformation != "None" )& (method != 'euclidean') &(method != 'uunifrac')) {
    ps_obj <- microbiome::transform(ps_obj, transformation)
  }
  
  # Calculate distance matrix based on the specified method
  dist <- phyloseq::distance(ps_obj, method =method)
  
  # Create formula from the string
  formula1 <- as.formula(paste("dist ~", formula))
  
  # Perform adonis analysis with the specified formula
  permanova_result <- adonis2(formula1, data = metadf,by='term',parallel = 30)
  
  # Process adonis results for plotting
  res <- as.data.frame(permanova_result)
  res$meta <- row.names(res)
  res <- subset(res, res$meta != 'Total')
  res$meta <- gsub('Residuals', 'other', res$meta)
  res[which(res$`Pr(>F)` > 0.05), 'signif'] <- ' '
  res[which(res$`Pr(>F)` <= 0.05), 'signif'] <- '*'
  res[which(res$`Pr(>F)` <= 0.01), 'signif'] <- '**'
  res[which(res$`Pr(>F)` <= 0.001), 'signif'] <- '***'
  res$Feature <- res$meta
  res$Feature <- factor(res$Feature, levels = res$Feature[order(res$R2, decreasing = FALSE)])
  res <- rownames_to_column(res, var = "new_name")
  res <- res %>%
    left_join(cols_meta, by = 'new_name')
  res <- res %>%
    filter(new_name != "Residual")
  # Create the ggplot object for plotting
  main_plot <- ggplot(res, aes(x = Feature, y = R2*100, fill = category)) +
    geom_col() +  # Use geom_col() to plot actual values
    labs(title = "",
         x = "",
         y = "R2") +
    theme_classic() +  # Optional: use a minimal theme
    coord_flip()+
    geom_text(aes(label=signif),  # new
              size=3)+ 
    theme(
      panel_background=element_rect(fill="white"),
      legend.position = 'top',               # Position legend at top left corner          # Anchor legend's top-left corner
      legend.title = element_blank(), 
      legend.text = element_text(),
      legend.direction = 'horizontal',
      legend.justification = "left",
      axis_title_y=element_blank(),
      axis_line_x=element_line(color="black"),
      axis_line_y=element_blank(),
      axis.text.y=element_text(size = 14, hjust = 1, family = "sans"),
      axis.text.x=element_text(size = 14, hjust = 1, family = "sans"),
      axis_ticks_major_y=element_blank(),
      panel_grid=element_blank(),
      panel_border=element_blank(),
    )+ 
    facet_grid(category ~ ., scales = "free", space = "free")+
    theme(strip.text.x = element_blank(),  # Remove x strip labels
          strip.text.y = element_blank())+
    guides(fill = guide_legend(nrow = 1))+
    
    # Use pnw_colors for filling
    scale_fill_manual(values = palette)
  main_plot
  
  total_effect_size <- round(sum(res$R2)*100,3)
  unexplained_variance <- 100 - total_effect_size
  
  
  # Prepare data for circular barplot
  circular_data <- data.frame(
    category = c("Explained Variance", "Unexplained Variance"),
    value = c(total_effect_size, unexplained_variance)
  )
  
  
  circular_data$fraction <- circular_data$value / sum(circular_data$value)
  
  # Compute the cumulative percentages (top of each rectangle)
  circular_data$ymax <- cumsum(circular_data$fraction)
  
  # Compute the bottom of each rectangle
  circular_data$ymin <- c(0, head(circular_data$ymax, n=-1))
  
  # Compute label position
  circular_data$labelPosition <- (circular_data$ymax + circular_data$ymin) / 2
  
  # Compute a good label
  circular_data$label <- paste0(circular_data$category, "\n value: ", circular_data$value)
  
  # Make the plot
  circular_plot <- ggplot(circular_data, aes(ymax=ymax, ymin=ymin, xmax=3.5, xmin=3, fill=category)) +
    geom_rect() +
    geom_label( x=4, aes(y=labelPosition, label=label), size=5) +
    scale_fill_manual(values = pnw_palette("Shuksan2", n=2))+
    coord_polar(theta="y") +
    xlim(c(2, 4)) +
    theme_void() +
    theme(legend.position = "none")
  
  combined_plot <- main_plot + circular_plot + plot_layout(ncol=2)
  
  # Check if plot should be generated
  if (show_plot) {
    print(combined_plot)
  }
  # Return both adonis results and the ggplot object
  return(list(permanova_result, combined_plot,res))
}
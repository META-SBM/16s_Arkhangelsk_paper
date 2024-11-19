
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(gridExtra)

library(phyloseq)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(lemon)
library(ggsignif)
library(rstatix)

# Define the function to create the ordination plots
create_ordination_plots <- function(ps_obj, method = "PCoA", distance_method = "wunifrac", group, size = 10,palette) {
  
  # Calculate distance matrix
  dist <- phyloseq::distance(ps_obj, method = distance_method)
  
  # Perform ordination
  ordination <- ordinate(ps_obj, method = method, distance = dist)
  P <- cbind(as.data.frame(ordination$vectors), as.data.frame(as.matrix(sample_data(ps_obj))))
  
  # Check if the group column exists
  if (!group %in% colnames(P)) {
    stop(paste("Column", group, "does not exist in the sample data."))
  }
  
  # Calculate percentage of variance explained by the first two axes
  variance_explained <- 100 * ordination$values$Relative_eig[1:2]
  variance_explained <- 100 * ordination$values$Relative_eig[1:2]
  # Calculate the means for each unique value of the group
  means <- P %>%
    group_by(!!rlang::sym(group)) %>%
    summarise(mean_Axis1 = mean(Axis.1), mean_Axis2 = mean(Axis.2))
  
  # Merge means back into the original data
  P <- P %>%
    dplyr::left_join(means, by = group)
  anosim_res <- anosim(dist, P[[group]],parallel = 30)
  anosim_text <- paste("ANOSIM R:", round(anosim_res$statistic, 3), "p-value:", round(anosim_res$signif, 3))
  
  # Create scatter plot with ellipses and means
  pl <-  ggplot(P, aes_string(x = "Axis.1", y = "Axis.2", color = group)) +
    geom_point(size = 1.5, alpha = 0.8) +
    stat_ellipse() +
    geom_segment(aes(x = Axis.1, y = Axis.2, xend = mean_Axis1, yend = mean_Axis2), alpha = 0.3) +
    geom_label(data = means, aes(x = mean_Axis1, y = mean_Axis2, label = !!rlang::sym(group)), fill = "white", color = "black", fontface = "bold", size = 3) +
    geom_text(x = Inf, y = Inf, label = anosim_text, hjust = 1.1, vjust = 1.1, size = size * 0.3, color = "black") +
    scale_color_manual(values=palette) +
    xlab(label = paste0('Axis.1 [', round(variance_explained[1], 1), '%]')) +
    ylab(label = paste0('Axis.2 [', round(variance_explained[2], 1), '%]')) +
    theme(
      axis.text.y = element_text(color = "black", size = size),
      axis.text.x = element_text(angle = 90, hjust = 1, size = size, color = 'black'),
      legend.position = "none",
      axis.title.y = element_text(color = "black", size = size, angle = 90),
      axis.title.x = element_text(color = "black", size = size),
      text = element_text(size = size, color = 'black'),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.background = element_rect(fill = "white", color = "black"),
      panel.grid.major = element_line(size = 0.1, linetype = 'solid', color = "gray"),
      panel.grid.minor = element_line(size = 0.1, linetype = 'solid', color = "gray")
    )
  # Perform Wilcoxon test for Axis.1
  # Perform pairwise comparisons for Axis.1
  # Define pairwise comparisons
  unique_vals <- unique(P[[group]])
  my_comparisons <- lapply(combn(unique_vals, 2, simplify = FALSE), as.vector)
  anno_df <- ggpubr::compare_means(as.formula(paste("Axis.1 ~", group)), data = P, method = "wilcox.test",p.adjust.method = "BH")%>%
    add_significance("p.adj") %>%
    add_x_position()%>%
    add_y_position( data=P,formula = as.formula(paste("Axis.1 ~", group)),step.increase = 0.2)
  # Create violin plot for Axis.1
  x_dens <- ggplot(P, aes_string(x = group, y = "Axis.1",color=group)) +
    geom_violin(trim = FALSE, alpha = 0.1) +
    geom_boxplot(width = 0.5, alpha = 0.75, position = position_dodge(0.9)) +
    geom_jitter(size = 1.5, alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9)) +
    scale_color_manual(values=palette) +
    theme(
      axis.text.y = element_text(color = "black", size = size),
      axis.text.x = element_text(angle = 90, hjust = 1, size = size, color = 'black'),
      legend.position = "none",
      axis.title.y = element_text(color = "black", size = size, angle = 90),
      axis.title.x = element_text(color = "black", size = size),
      panel.background = element_rect(fill = "white", color = "black"),
      text = element_text(size = size, color = 'black'),
      panel.grid.major = element_line(size = 0.1, linetype = 'solid', color = "gray"),
      panel.grid.minor = element_line(size = 0.1, linetype = 'solid', color = "gray")
    )+
    stat_pvalue_manual(
      anno_df,  label = "p.adj.signif", tip.length = 0.02,
      step.increase = 0.05,coord.flip = TRUE
    )+
    coord_flip()
  
  # Perform Wilcoxon test for Axis.2
  test_result_axis2 <- ggpubr::compare_means(as.formula(paste("Axis.2 ~", group)), data = P, method = "wilcox.test")
  anno_df <- ggpubr::compare_means(as.formula(paste("Axis.2 ~", group)), data = P, method = "wilcox.test",p.adjust.method = "BH")%>%
    add_significance("p.adj") %>%
    add_x_position()%>%
    add_y_position( data=P,formula = as.formula(paste("Axis.2 ~", group)),step.increase = 0.2)
  # Create violin plot for Axis.2
  y_dens <- ggplot(P, aes_string(x = group, y = "Axis.2", color = group)) +
    geom_violin(trim = FALSE, alpha = 0.1) +
    geom_boxplot(width = 0.5, alpha = 0.75, position = position_dodge(0.9)) +
    geom_jitter(size = 1.5, alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9)) +
    scale_color_manual(values=palette) +
    theme(
      axis.text.y = element_text(color = "black", size = size),
      axis.text.x = element_text(angle = 90, hjust = 1, size = size, color = 'black'),
      legend.position = "none",
      axis.title.y = element_text(color = "black", size = size, angle = 90),
      axis.title.x = element_text(color = "black", size = size),
      plot.title = element_text(size = 25),
      panel.background = element_rect(fill = "white", color = "black"),
      text = element_text(size = size, color = 'black'),
      panel.grid.major = element_line(size = 0.1, linetype = 'solid', color = "gray"),
      panel.grid.minor = element_line(size = 0.1, linetype = 'solid', color = "gray")
    ) +
    stat_pvalue_manual(
      anno_df,  label = "p.adj.signif", tip.length = 0.02,
      step.increase = 0.05,coord.flip = FALSE
    )
  
  l <- ggplot(P, aes_string(x = "Axis.1", y = "Axis.2", color = group)) +
    geom_point(size = 2.5, alpha = 0.8) +
    stat_ellipse() +
    geom_segment(aes(xend = mean(Axis.1), yend = mean(Axis.2)), alpha = 0.3) +
    geom_label(aes(x = mean(Axis.1), y = mean(Axis.2), label = !!rlang::sym(group)), fill = "white", color = "black", fontface = "bold", size = 3) +
    scale_color_manual(values=palette) +
    xlab(label = paste0('Axis.1 [', round(variance_explained[1], 1), '%]')) +
    ylab(label = paste0('Axis.2 [', round(variance_explained[2], 1), '%]')) +
    theme(
      axis.text.y = element_text(color = "black", size = size),
      axis.text.x = element_text(angle = 90, hjust = 1, size = size, color = 'black'),
      legend.position = "right",
      axis.title.y = element_text(color = "black", size = size, angle = 90),
      axis.title.x = element_text(color = "black", size = size),
      text = element_text(size = size, color = 'black'),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.background = element_rect(fill = "white", color = "black"),
      panel.grid.major = element_line(size = 0.01, linetype = 'solid', color = "gray"),
      panel.grid.minor = element_line(size = 0.01, linetype = 'solid', color = "gray")
    ) 
  # Create a blank plot
  legend <- g_legend(l)
  
  # Make grobs
  gA <- ggplotGrob(x_dens)
  gB <- ggplotGrob(pl)
  gD <- ggplotGrob(y_dens)
  gL <-  legend
  # Get width
  xWidth = unit.pmax( gA$widths[2:4], gB$widths[2:4])
  yHeight = unit.pmax(gB$heights[4:5], gD$heights[4:5])
  
  # Set the widths
  gA$widths[2:3] <- xWidth
  gB$widths[2:3] <- xWidth
  # Set the heights
  gB$heights[4:5] <- yHeight
  gD$heights[3:5] <- yHeight
  
  p =grid.arrange( gD,gB,gL ,gA,ncol=2, nrow=2, widths=c(2, 5), heights=c(5, 2))
  
}

plot_alpha_div <- function(
    ps_obj,
    group ,
    color ,
    measure 
) {
  # Generate boxplot of alpha diversity without transformation
  p <- plot_richness(ps_obj, x = group, color = color, measures = measure)
  
  return(p)
}
create_alpha_plots <- function(ps_obj,col,measure = 'Shannon',method = 'wilcox.test',color,my_comparisons,size=10){
  
  p <- plot_alpha_div(ps_obj, group = col, color = col, measure = measure) +
    geom_violin(trim=F, alpha=0.1) +
    geom_boxplot(width=0.5, alpha=0.75, position=position_dodge(0.9)) +
    geom_jitter(size=1.5, alpha=0.5, position=position_jitterdodge(jitter.width=0.1, dodge.width=0.9)) + scale_color_manual(values=color)+
    theme_bw(base_size=20)+
    theme(
      plot.title = element_text(size = 25),
      axis.text.y = element_text(color = "black", size = size),
      axis.text.x = element_text(angle=90, hjust=1,size=15,color = 'black'),
      legend.position = "none",
      axis.title.y  = element_text(color = "black", size = size,angle=90),
      axis.title.x  = element_text(color = "black", size = size),
      #legend.key.size = unit(0.5, 'cm'),
      text = element_text(size = size,colour ='black'))+
    scale_color_manual(values=color)
  
  group <- col
  df <- as.data.frame(p$data)
  anno_df <- ggpubr::compare_means(as.formula(paste("value ~", group)), data = df, method = "wilcox.test",p.adjust.method = "BH")%>%
    add_significance("p.adj") %>%
    add_x_position()%>%
    add_y_position( data=df,formula = as.formula(paste("value ~", group)),step.increase = 0.2)
  
  p <- p +
    stat_pvalue_manual(
      anno_df,  label = "p.adj.signif", tip.length = 0.02,
      step.increase = 0.05,coord.flip = FALSE
    )
  return(p)
}

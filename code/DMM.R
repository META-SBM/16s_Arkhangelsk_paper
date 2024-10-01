library(ggplot2)
library(dplyr)
library(gridExtra)
library(paletteer)
library(tidyr)  # For pivot_longer
library(stringr)  # For str_wrap
library(cowplot)


# Define a function to create the plots
create_plot_rel <- function(cols, palette,ps_meta) {
  dataex <- ps_meta %>%
    dplyr::group_by_at(cols) %>%
    dplyr::summarise(N = n()) %>%
    dplyr::mutate(rel_count = N / sum(N)) 
  
  # Perform chi-squared test
  chisq_result <- chisq.test(ps_meta[[cols[1]]], ps_meta[[cols[2]]])
  chisq_pvalue <- format(chisq_result$p.value, digits = 3)
  
  ggplot(dataex, aes(x = !!sym(cols[[1]]), y = rel_count, fill = !!sym(cols[[2]]))) +
    geom_bar(position = "stack", stat = "identity") +
    geom_text(aes(label = scales::percent(rel_count, accuracy = 1)), position = position_stack(vjust = 0.5), size = 6)  +
    labs(title = paste("Chi-squared p-value:", chisq_pvalue ), fill = cols[[2]],tag = 'B') +
    scale_fill_manual(values = palette) +
    theme_linedraw()+
    theme(axis.text.x = element_text(hjust = 1),plot.tag=element_text(size=14,face="bold"),
          legend.text =  element_text(size = 14, family = "Fira Sans"),
          legend.title = element_text(size = 14,family = "Fira Sans",face = "bold"))
}

# Define a function to create the plots
create_plot_ab <- function(cols, palette,ps_meta) {
  dataex <- ps_meta %>%
    dplyr::group_by_at(cols) %>%
    dplyr::summarise(N = n()) %>%
    dplyr::mutate(rel_count = N / sum(N)) 
  
  # Perform chi-squared test
  chisq_result <- chisq.test(ps_meta[[cols[1]]], ps_meta[[cols[2]]])
  chisq_pvalue <- format(chisq_result$p.value, digits = 3)
  
  ggplot(dataex, aes(x = !!sym(cols[[1]]), y = N, fill = !!sym(cols[[2]]))) +
    geom_bar(position = "stack", stat = "identity") +
    geom_text(aes(label = N), position = position_stack(vjust = 0.5), size = 6)  +
    labs(title = paste("Chi-squared p-value:", chisq_pvalue ), fill = cols[[2]],tag = 'B') +
    scale_fill_manual(values = palette) +
    theme_linedraw()+
    theme(axis.text.x = element_text(hjust = 1),plot.tag=element_text(size=14,face="bold"),
          legend.text =  element_text(size = 14, family = "Fira Sans"),
          legend.title = element_text(size = 14,family = "Fira Sans",face = "bold"))
}

barplot_DMM <- function(fit, best, ps_obj, column_combinations,size=c(1,3), tag = 'None',top=10,annotation_table = FALSE,
                        plots,tt) {
  
  p0 <- fitted(fit[[1]], scale=TRUE)
  pk <- fitted(best, scale=TRUE)
  
  diff <- rowSums(abs(pk - as.vector(p0)))
  o <- order(diff, decreasing=TRUE)
  cdiff <- cumsum(diff[o]) / sum(diff)
  df <- head(cbind(Mean=p0[o], pk[o,], diff=diff[o], cdiff), top)
  df <- as.data.frame(df)
  
  if(tag == 'other') {
    # Compute values for 'other' row
    other_row <- 1 - colSums(df)
    df <- rbind(df, other = other_row)
  }
  
  df$plot <- tt[rownames(df), "plot"]
  df$plot[is.na(df$plot)] <- "other"
  
  df_selected <- df %>% select(-c(Mean,diff))
  colnames(df_selected) <- c(seq(1,ncol(pk)), 'cdiff', 'plot')
  
  df_long <- df_selected %>%
    pivot_longer(cols = -c(plot, cdiff), names_to = "Sample", values_to = "Count") %>%
    group_by(Sample) %>%
    mutate(Sample = reorder(Sample, -Count))
  df_long$plot <- str_wrap(df_long$plot, width = 10)  # Adjust width as needed
  
  # Create cdiff table
  cdiff_table <- df %>% select(cdiff)
  cdiff_table$cdiff <- round(cdiff_table$cdiff, 3)
  
  if(tag == 'other') {
    cdiff_table <- subset(cdiff_table, rownames(cdiff_table) != "other")
  }
  
  first <- cdiff_table[1,1]
  cdiff_table <- cdiff_table %>%
    mutate(cdiff = cdiff - lag(cdiff))
  
  cdiff_table[1,1] <- first
  colnames(cdiff_table) <- 'abs(cdiff)'
  
  # Create the main plot
  p <- ggplot(df_long, aes(x = Sample, y = Count, fill = plot)) +
    geom_bar(position = "stack", stat="identity") +
    scale_fill_manual("legend", values = palette, name = 'ASV') +
    geom_text(aes(label = round(Count,3)), position = position_stack(vjust = 0.5), size = 6) +
    labs(x = "Number of Dirichlet mixture components", y = "Relative abundance top 10 drivers",tag = 'C') +
    theme_linedraw() +
    theme(axis.text.x = element_text(hjust = 1),plot.tag=element_text(size=14,face="bold"),
          legend.text =  element_text(size = 14, family = "Fira Sans"),
          legend.title = element_text(size = 14,family = "Fira Sans",face = "bold"))
  
  # Add annotation table to the plot
  if( annotation_table == TRUE){
    p <- p + annotation_custom(
      grob = tableGrob(cdiff_table),
      xmin = ncol(pk) + 1.5,
      xmax = ncol(pk) + 1.5,
      ymin = sum(df_long$Count)/ncol(pk),
      ymax = sum(df_long$Count)/ncol(pk)
    )
  }
  

  # Create additional plots
  grob_list <- lapply(plots, ggplotGrob)
  gA <- grob_list[[1]]
  gB <- ggplotGrob(p)
  
  
  p <- plot_grid(gA,gB,align = "v",nrow = 2, rel_heights = c(1/3, 2/3))
  p
  

  #png('/mnt/disk1/PROJECTS/DAVID_WGS/different_figures/results_16s_for_publishing/det_5_prev_0.1_count_DMMclustering.png',height = 15, width = 15, units="in",res=500)
  #grid.arrange(plot_grid(p1,labels='A'), p, ncol=1, heights=size)
  #dev.off()
  print(df)
  return(p)
}

create_dmm_plot <- function(ps_obj,fit,best,column_combinations,relative = TRUE,
                            size=c(1,3),top=10,palette,tt,annotation_table = FALSE){

  asss <- mixture(best, assign = TRUE)
  sample_data(ps_obj)$best <- as.factor(asss)
  
  ps_meta <- as(sample_data(ps_obj), 'data.frame')
  
  if(relative == FALSE){
    plots <- lapply(column_combinations, create_plot_ab, palette = paletteer_c("ggthemes::Sunset-Sunrise Diverging",10),ps_meta)
    
  }else{
    plots <- lapply(column_combinations, create_plot_rel, palette = pnw_palette("Spring", n=5),ps_meta)
    
  }
  dmm_plot <- barplot_DMM(fit,best = best,ps_obj,column_combinations,plots = plots,tt = tt,annotation_table = annotation_table)

  return(dmm_plot)

}


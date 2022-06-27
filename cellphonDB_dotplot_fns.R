## downloaded from cellphoneDB github
## https://github.com/Teichlab/cellphonedb/blob/master/cellphonedb/src/plotters/R/plot_dot_by_column_name.R
## ------------------------------------------------------------------------------------ ##
library(ggplot2)
library(dplyr)
dot_plot = function(selected_rows = NULL, selected_columns = NULL, selected_interactions_save = as.logical(F), selected_interactions_fname = NULL, filename = 'plot.pdf', 
                    width = NULL, height = NULL, means_path = './means.txt', pvalues_path = './pvalues.txt', 
                    means_separator = '\t', pvalues_separator = '\t', output_extension = '.pdf', plotDir = getwd(), pvalue=0.05, min.mean = NULL, max.mean = NULL ){
  
  all_pval = read.table(pvalues_path, header=T, stringsAsFactors = F, sep=means_separator, comment.char = '', check.names=F)
  all_means = read.table(means_path, header=T, stringsAsFactors = F, sep=pvalues_separator, comment.char = '', check.names=F)
  if (!is.null(selected_rows)) {
    annotation = all_pval %>% dplyr::filter(interacting_pair %in% selected_rows ) %>% select(1:11)
    print(sprintf("Out of %s selected_rows, %s interaction pairs selected", length(selected_rows), dim(annotation)[1]))
    if (selected_interactions_save) write.table(annotation, selected_interactions_fname, quote = F, sep = '\t', row.names = F, col.names = T)
  }
  
  intr_pairs = all_pval$interacting_pair
  all_pval = all_pval[,-c(1:11)] ##remove first 11 column description, only results are left, the same for all_means
  all_means = all_means[,-c(1:11)]
  
  if(is.null(selected_rows)){
    selected_rows = intr_pairs
  }
  
  if(is.null(selected_columns)){
    selected_columns = colnames(all_pval)
  }
  
  sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
  sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]
  ## ----
  print(sprintf("%s interaction pairs will be plotted", dim(sel_pval)[1]))
  df_names = expand.grid(selected_rows, selected_columns)
  sel_pval[is.na(sel_pval)] = 1
  pval = unlist(sel_pval)
  # pval[pval==0] = 0.0009
  print(sprintf("%s pvalue originally are 0", sum(pval==0)))
  pval[pval==0] = 1e-5
  pval[pval>pvalue] = NA
  print("ETERTETEAWTW")
  print(sprintf("%s pvalue < %s", sum(!is.na(pval)), pvalue ))
  plot.data = cbind(df_names,pval)
  pr = unlist(as.data.frame(sel_means))
  pr[pr==0] = 1
  plot.data = cbind(plot.data,log2(pr))
  colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')
  ## mean value of p-vaue<0.05 change into NA
  plot.data$mean[is.na(plot.data$pvalue)] = NA
  my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)
  
  if (is.null(min.mean) & is.null(max.mean)) {
    print(sprintf("Default max.mean is %s, min.mean is %s ", max(plot.data$mean, na.rm = T), min(plot.data$mean, na.rm = T)))
    
    ggplot(plot.data,aes(x=clusters,y=pair)) +
      geom_point(aes(size=-log10(pvalue),color=mean)) +
      scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) + 
      # scale_size_area() + ##make zero value to have zero size
      theme_bw() +
      theme(panel.grid.minor = element_blank(),
            # panel.grid.major = element_blank(),
            axis.text=element_text(size=14, colour = "black"),
            axis.text.x = element_text(angle = 90, hjust = 1),
            # axis.text.y = element_text(size=12, colour = "black"),
            axis.text.y = element_text(size=12, colour = "red"),
            axis.title=element_blank(),
            panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
  } else {
    if (is.null(min.mean)) {
      min.mean = min(plot.data$mean, na.rm = T)
    } else {
      min.mean = min.mean
    }
    
    if (is.null(max.mean)) {
      max.mean = max(plot.data$mean, na.rm = T)
    } else {
      max.mean = max.mean
    }
    print(sprintf("specified max.mean is %s, min.mean is %s ", max.mean, min.mean))
    ggplot(plot.data,aes(x=clusters,y=pair)) +
      geom_point(aes(size=-log10(pvalue),color=mean)) +
      scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette, limits = c(min.mean, max.mean) ) + 
      # scale_size_area() + ##make zero value to have zero size
      theme_bw() +
      theme(panel.grid.minor = element_blank(),
            # panel.grid.major = element_blank(), 
            axis.text=element_text(size=14, colour = "black"),
            axis.text.x = element_text(angle = 90, hjust = 1),
            axis.text.y = element_text(size=12, colour = "black"),
            axis.title=element_blank(),
            panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
  }
  
  if (is.null(width)) width  = 3+dim(sel_pval)[2]
  if (is.null(height)) height = round(0.25*dim(sel_pval)[1], digits = 0)
  if (output_extension == '.pdf') {
    ggsave(filename = file.path(plotDir, filename), width = width, height = height, device = "pdf", limitsize=F)
  }
  else {
    ggsave(filename = file.path(plotDir, filename), width = width, height = height, limitsize=F)
  }
  return(list('selected_rows' = selected_rows, 'plot.data' = plot.data))
}
## ------------------------------------------------------------------------------------ ##

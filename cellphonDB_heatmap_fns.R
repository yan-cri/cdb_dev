## This script is cellphoneDB heatmap fn, disable the table writing, re-make the heatmap plot without log transformation
## nameUpdate: whether to update results of row/column names (the same orthogonal) from 'name2update' to 'updatedName' 
## maxNo: if provided, use the same max number for the plot.
## orderNames: if provided, re-order heatmap row/column names order.
## ------------------------------------------------------------------------------------ ##
library(pheatmap)
heatmaps_plot = function(meta_file, nameUpdate = F, name2update, updatedName, orderNames = NULL, maxNo = NA, pvalues_file, count_filename, log_filename, count_network_separator, interaction_count_separator, show_rownames = T, show_colnames = T,
                         scale="none", cluster_cols = T, border_color='white', cluster_rows = T, fontsize_row=11,
                         fontsize_col = 11, main = '',treeheight_row=0, family='Arial', treeheight_col = 0,
                         col1 = "dodgerblue4", col2 = 'peachpuff', col3 = 'deeppink4', meta_sep='\t', pvalues_sep='\t', pvalue=0.05, plot.width = 7, plot.height = 7){
  #######   Network
  
  meta = read.csv(meta_file, comment.char = '', sep=meta_sep)
  
  all_intr = read.table(pvalues_file, header=T, stringsAsFactors = F, sep=pvalues_sep, comment.char = '', check.names = F)
  intr_pairs = all_intr$interacting_pair
  all_intr = all_intr[,-c(1:11)]
  
  
  split_sep = '\\|'
  join_sep = '|'
  
  pairs1_all = unique(meta[,2])
  
  pairs1 = c()
  for (i in 1:length(pairs1_all)) # would apply functions be more efficient computationally
    for (j in 1:length(pairs1_all))
      pairs1 = c(pairs1,paste(pairs1_all[i],pairs1_all[j],sep=join_sep))
  
  all_count = matrix(ncol=3)
  colnames(all_count) = c('SOURCE','TARGET','count')
  count1 = c()
  for(i in 1:length(pairs1))
  {
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]
    
    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]
    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]
    
    if(p1!=p2)
      count1 = length(unique(n1))+length(unique(n2))
    else
      count1 = length(unique(n1))
    
    new_count = c(p1,p2,count1)
    names(new_count) = c('SOURCE','TARGET','count')
    all_count = rbind(all_count, new_count)
  }
  
  all_count = all_count[-1,]
  # write.table(all_count, count_network_filename, sep=count_network_separator, quote=F, row.names = F)
  
  #######   count interactions
  
  count1 = c()
  for(i in 1:length(pairs1))
  {
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]
    
    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]
    
    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]
    if(p1!=p2)
      count1 = c(count1,length(unique(n1))+length(unique(n2)))
    else
      count1 = c(count1,length(unique(n1)))
    
  }
  
  if (any(count1)>0)
  {
    count_matrix = matrix(count1, nrow=length(unique(meta[,2])), ncol=length(unique(meta[,2])))
    rownames(count_matrix)= unique(meta[,2])
    colnames(count_matrix)= unique(meta[,2])
    print(head(count_matrix))
    print(dim(count_matrix))
    print('-=-=-=-=-')
    
    if (nameUpdate) {
      rownames(count_matrix) = gsub(pattern = name2update, replacement = updatedName, x = rownames(count_matrix))
      colnames(count_matrix) = gsub(pattern = name2update, replacement = updatedName, x = colnames(count_matrix))
    }
    
    all_sum = rowSums(count_matrix)
    all_sum = cbind(names(all_sum), all_sum)
    # write.table(all_sum, file=interaction_count_filename, quote=F, sep=count_network_separator, row.names=F)
    print(head(count_matrix))
    if(is.na(maxNo)) {
      col.heatmap <- colorRampPalette(c(col1,col2,col3 ))( 1000 )
    } else {
      breaksList  <- seq(0, maxNo, by = 0.1)
      col.heatmap <- colorRampPalette(c(col1,col2,col3 ))(length(breaksList))
    }
    # col.heatmap <- colorRampPalette(c("navy", "white", "firebrick3"))(1000)
    print(sprintf("maximum number of count is %s", max(count_matrix)))
    if (!is.null(orderNames)) {
      count_matrix = count_matrix[match(orderNames, rownames(count_matrix)),match(orderNames, colnames(count_matrix))]
    }
    
    if(is.na(maxNo)) {
      pheatmap(count_matrix, show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, 
               cluster_cols = cluster_cols, cluster_rows = cluster_rows, 
               treeheight_row = treeheight_row, treeheight_col = treeheight_col, 
               border_color=border_color, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
               main = main, family = family,color = col.heatmap, filename = count_filename, 
               width = plot.width, height = plot.height)
    } else {
      print(sprintf("maximum number of legend is %s", maxNo))
      pheatmap(count_matrix, show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, 
               cluster_cols = cluster_cols, cluster_rows = cluster_rows, 
               treeheight_row = treeheight_row, treeheight_col = treeheight_col, 
               border_color=border_color, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
               main = main, family = family,color = col.heatmap, filename = count_filename,
               breaks = breaksList, 
               width = plot.width, height = plot.height)
    }
    
    
    # pheatmap(log(count_matrix+1), show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, cluster_cols = cluster_cols,
    #          border_color=border_color, cluster_rows = cluster_rows, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
    #          main = main, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col, filename = log_filename)
  } else {
    stop("There are no significant results using p-value of: ", pvalue, call.=FALSE)
  }
}
## ------------------------------------------------------------------------------------ ##
library(ggplot2)
dot_plot = function(selected_rows = NULL,
                    selected_columns = NULL,
                    filename = 'plot.png',
                    width = 8,
                    height = 10,
                    means_path = './means.txt',
                    pvalues_path = './pvalues.txt',
                    means_separator = '\t',
                    pvalues_separator = '\t',
                    output_extension = '.png'
){
  
  all_pval = read.table(pvalues_path, header=T, stringsAsFactors = F, sep=means_separator, comment.char = '', check.names=F)
  all_means = read.table(means_path, header=T, stringsAsFactors = F, sep=pvalues_separator, comment.char = '', check.names=F)
  
  intr_pairs = all_pval$interacting_pair
  all_pval = all_pval[,-c(1:11)]
  all_means = all_means[,-c(1:11)]
  
  if(is.null(selected_rows)){
    selected_rows = intr_pairs
  }
  
  if(is.null(selected_columns)){
    selected_columns = colnames(all_pval)
  }
  
  sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
  sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]
  
  df_names = expand.grid(selected_rows, selected_columns)
  pval = unlist(sel_pval)
  pval[pval==0] = 0.0009
  plot.data = cbind(df_names,pval)
  pr = unlist(as.data.frame(sel_means))
  pr[pr==0] = 1
  plot.data = cbind(plot.data,log2(pr))
  colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')
  
  my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)
  
  ggplot(plot.data,aes(x=clusters,y=pair)) +
    geom_point(aes(size=-log10(pvalue),color=mean)) +
    scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text=element_text(size=14, colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.title=element_blank(),
          panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
  
  if (output_extension == '.pdf') {
    ggsave(filename, width = width, height = height, device = cairo_pdf, limitsize=F)
  }
  else {
    ggsave(filename, width = width, height = height, limitsize=F)
  }
}
## ------------------------------------------------------------------------------------ ##

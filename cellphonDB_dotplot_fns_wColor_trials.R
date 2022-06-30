## downloaded from cellphoneDB github
## https://github.com/Teichlab/cellphonedb/blob/master/cellphonedb/src/plotters/R/plot_dot_by_column_name.R
## selected_rows.nameReorder: whether to re-order ligand-receptor names based on DB inputs
## selected.resSave: whether to save selected results including p-value, meta-value and annotation results into file name with 'selected.resFnamePrefix' prefix.
## selected_rows: select ligand-receptor pairs
## selected_columns: select cell interactions, consistent with column names
## selected_columns_rename, subnames: whether to rename column names with 'subnames'
## ------------------------------------------------------------------------------------ ##
library(ggplot2)
library(dplyr)
dot_plot = function(selected_rows = NULL, selected_rows.nameReorder = as.logical(F), selected.resSave = as.logical(F), selected.resFnamePrefix = NULL, 
                    selected_columns = NULL, selected_columns_rename = as.logical(F), subnames = NULL, filename = 'plot.pdf', 
                    nonSig.removal = as.logical(F), 
                    width = NULL, height = NULL, means_path = './means.txt', pvalues_path = './pvalues.txt', 
                    means_separator = '\t', pvalues_separator = '\t', output_extension = '.pdf', plotDir = getwd(), pvalue=0.05, min.mean = NULL, max.mean = NULL, yAxisCol = NULL, debug = as.logical(F)){
  ## read in 'pvalues_path' and 'means_path' as 'all_pval' and 'all_means'
  all_pval = read.table(pvalues_path, header=T, stringsAsFactors = F, sep=means_separator, comment.char = '', check.names=F)
  all_means = read.table(means_path, header=T, stringsAsFactors = F, sep=pvalues_separator, comment.char = '', check.names=F)
  
  # use the first 4 lines
  first4rows = all_pval[c(1:4),]
  # deliberately add some duplicates
  first4rows = rbind(first4rows, first3rows[c(1, 4, 2),])
  first4unique = first4rows[!duplicated(first4rows),]
  # it works!
  
  # so lets try this on the entire dataset
  unique_pval = all_pval[!duplicated(all_pval),]
  # Holy shamoley, I got the entire dataset back!
  
  # Let's alphabetize the IDs and check for duplicates
  sorted_pval = order(unique_pval$id_cp_interaction)
  # WHAAAAAT? All that was returned were the indexes of the columns with IDs in alphabetical order
  # Well then, let's count those numbers in order and look for dupes
  sort(sorted_pval)
  # Look ma, no dupes!
  # So, the data must be clean, huh?
  
  
  ## ------------------ ##
  ## if selected_rows!=NULL, select certain rows (ligand-receptor interactions) and output 'annotation' variable.
  if (!is.null(selected_rows)) {
    annotation = all_pval %>% dplyr::filter(interacting_pair %in% selected_rows ) %>% select(1:11)
  } else {
    annotation = all_pval %>% select(1:11)
  }
  ## ------------------ ##
  ## whether to reorder the selected ligand-receptor interactions based on column 'receptor_b'
  if (selected_rows.nameReorder) {
    annotation$interacting_pair_reName <- NA
    ## add column 'gene_a2' and 'gene_b2' to update gene_a/b with partner_a/b with complex: names
    annotation$gene_a2 <- annotation$gene_a
    annotation$gene_a2[which(annotation$gene_a2=="")] <- gsub('complex:', '', annotation$partner_a[which(annotation$gene_a2=="")])
    annotation$gene_b2 <- annotation$gene_b
    annotation$gene_b2[which(annotation$gene_b2=="")] <- gsub('complex:', '', annotation$partner_b[which(annotation$gene_b2=="")])
    ## update 'interacting_pair_reName' based on 'receptor_b'
    annotation$interacting_pair_reName[which(annotation$receptor_b=='True')]  <- annotation$interacting_pair[which(annotation$receptor_b=='True')]
    annotation$interacting_pair_reName[which(annotation$receptor_b=='False')] <- paste(annotation$gene_b2[which(annotation$receptor_b=='False')], annotation$gene_a2[which(annotation$receptor_b=='False')], sep = '_')
    if (any(is.na(annotation$interacting_pair_reName))) {
      print(sprintf("NOTE: %s pairs (%s) with no receptor marked from DB", sum(is.na(annotation$interacting_pair_reName)), paste(annotation$interacting_pair[which(is.na(annotation$interacting_pair_reName))], collapse = ', ')))
      annotation$interacting_pair_reName[which(is.na(annotation$interacting_pair_reName))] <- annotation$interacting_pair[which(is.na(annotation$interacting_pair_reName))]
    }
    if(debug) print(head(annotation))
  }
  ## ------------------ ##
  ## Printing out message based on different conditions
  if (!is.null(selected_rows)) {
    print(sprintf("Out of %s selected_rows, %s interaction pairs selected", length(selected_rows), dim(annotation)[1]))
  } else {
    if (!nonSig.removal) {
      print(sprintf("All %s interaction pairs will be plotted", dim(annotation)[1]))
    } else {
      print(sprintf("Some non-significant interactions will be removed from %s interaction pairs.", dim(annotation)[1]))
    }
  }
  ## ------------------ ##
  ## setup several variables, including 'intr_pairs' (all interacting_pair), 'all_pval.anno', 'all_means.anno', 'all_pval', 'all_means', 'selected_rows', 'selected_columns'
  intr_pairs = all_pval$interacting_pair
  all_pval.anno = all_pval[,1:11]
  all_means.anno =all_means[,1:11]
  all_pval = all_pval[,-c(1:11)] ##remove first 11 colum description, only results are left, the same for all_means
  all_means = all_means[,-c(1:11)]
  ## ------------------ ##
  ## remove duplicated 'intr_pairs' if 'selected_rows==NULL'
  if(is.null(selected_rows)){
    selected_rows = unique(intr_pairs)
    intr_pairs    = unique(intr_pairs)
  }
  ## if not specify 'selected_columns' all columns (cell-cell interactions) in 'all_pval' are selected
  if(is.null(selected_columns)){
    selected_columns = colnames(all_pval)
  }
  if (debug) {
    print('00000000000000000')
    print("below 2 variables should have the same length")
    print(sprintf("%s interactions rows", length(selected_rows)))
    print(sprintf("%s interactions in 'intr_pairs'.", length(intr_pairs)))
    print('00000000000000000')
  }
  ## ------------------ ##
  ## sort 'all_pval' and 'all_means' based on 'selected_rows' into variables 'sel_pval' and 'sel_means, meanwhile remove duplicated rows
  sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
  sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]
  if (debug) {
    print('00000000000000000')
    print("below 3 variables should have the same length.")
    print(sprintf("%s interactions in 'intr_pairs'.", length(intr_pairs)))
    print(sprintf("'all_pval' dimensions rows = %s", dim(sel_pval)[1]))
    print(sprintf("'all_means' dimensions rows = %s", dim(sel_means)[1]))
    print('00000000000000000')
  }
  ## ------------------ ##
  ## remove rows with non-significant interactions
  if (nonSig.removal) {
    ## if sel_pval > defined 'pvalue' (significant level), all values set to NA, then use NA numbers in each row to define whether rows are need to be removed or not.
    sel_pval2 <- sel_pval
    sel_pval2[sel_pval2 > pvalue] <- NA
    row2remove <- rowSums(is.na(sel_pval2))==dim(sel_pval2)[2]
    ## remove non-significant rows to update variables 'sel_pval', 'sel_means', 'all_pval.anno', 'all_means.anno', 'selected_rows', 'annotation', 'intr_pairs'
    sel_pval  = sel_pval2[!row2remove,]
    sel_means = sel_means[!row2remove,]
    all_pval.anno = all_pval.anno[!row2remove,]
    all_means.anno = all_means.anno[!row2remove,]
    selected_rows  = selected_rows[!row2remove]
    annotation = annotation[!row2remove,]
    intr_pairs  = intr_pairs[!row2remove]
    if (debug) {
      print('00000000000000000')
      print("After 'nonSig.removal', below 3 variables should have the same length.")
      print(sprintf("%s interactions rows", length(selected_rows)))
      print(sprintf("%s interactions in 'intr_pairs'.", length(intr_pairs)))
      print(sprintf("'all_pval' dimensions rows = %s", dim(sel_pval)[1]))
      print(sprintf("'all_means' dimensions rows = %s", dim(sel_means)[1]))
      print('00000000000000000')
    }
  }
  ## --------------------- ##
  ## whether to save results in 'annotation', 'sel_pval', 'sel_means'
  if (selected.resSave) {
    sel_means[sel_means==0] = 1
    sel_pval[sel_pval > pvalue] = NA
    sel_means2 <- log2(sel_means)
    pval.res <- cbind(all_pval.anno[match(selected_rows, intr_pairs), ], sel_pval)
    mean.res <- cbind(all_means.anno[match(selected_rows, intr_pairs),], sel_means2)
    if (selected_rows.nameReorder) {
      pval.res$interacting_pair_reName <- NA
      pval.res$interacting_pair_reName[which(pval.res$receptor_b=='True')]  <- pval.res$interacting_pair[which(pval.res$receptor_b=='True')]
      pval.res$interacting_pair_reName[which(pval.res$receptor_b=='False')] <- paste(pval.res$gene_b[which(pval.res$receptor_b=='False')], pval.res$gene_a[which(pval.res$receptor_b=='False')], sep = '_')
      ## -
      mean.res$interacting_pair_reName <- NA
      mean.res$interacting_pair_reName[which(mean.res$receptor_b=='True')]  <- mean.res$interacting_pair[which(mean.res$receptor_b=='True')]
      mean.res$interacting_pair_reName[which(mean.res$receptor_b=='False')] <- paste(mean.res$gene_b[which(mean.res$receptor_b=='False')], mean.res$gene_a[which(mean.res$receptor_b=='False')], sep = '_')
      ## -
    }
    # pval.res2 <- pval.res[order(nrow(pval.res):1),]
    # mean.res2 <- mean.res[order(nrow(mean.res):1),]
    write.table(x =annotation, file = sprintf('%s_selected_annotations.txt', selected.resFnamePrefix), quote = F, sep = '\t', row.names = F, col.names = T)
    write.table(x = pval.res, file = sprintf('%s_selected_pval.txt', selected.resFnamePrefix), quote = F, sep = '\t', row.names = F, col.names = T)
    write.table(x = mean.res, file = sprintf('%s_selected_means.txt', selected.resFnamePrefix), quote = F, sep = '\t', row.names = F, col.names = T)
  }
  ## if 'selected_columns_rename'=T, use 'subnames' (if provided) to change cell-cell interaction names.
  ## if 'subnames' provided, 'subnames' and '|" will be removed from cell-cell interaction names defined in 'selected_columns'.
  ## otherwise cell-cell interaction names are displayed as '(cell1)_(cell2)'
  if (selected_columns_rename) {
    if (!is.null(subnames)) {
      selected_columns = gsub(subnames, '', selected_columns)
      selected_columns = gsub('[|]', '_', selected_columns)
    } else {
      selected_columns = 
        sprintf('%s(%s)_%s(%s)', 
                sapply(strsplit(sapply(strsplit(selected_columns, split = '[|]'), '[[', 1), split = '[_]'), '[[', 2), 
                sapply(strsplit(sapply(strsplit(selected_columns, split = '[|]'), '[[', 1), split = '[_]'), '[[', 1), 
                sapply(strsplit(sapply(strsplit(selected_columns, split = '[|]'), '[[', 2), split = '[_]'), '[[', 2), 
                sapply(strsplit(sapply(strsplit(selected_columns, split = '[|]'), '[[', 2), split = '[_]'), '[[', 1) )
    }
  }
  ## print message for final plot interactions
  print(sprintf("Finally, %s interaction pairs will be plotted", dim(sel_pval)[1]))
  ## ---------------------------------- ##
  ## define 'df_names' used for plot
  df_names = expand.grid(selected_rows, selected_columns)
  sel_pval[is.na(sel_pval)] = 1 ##NA sel_pval set to 1 as non-signifcant
  pval = unlist(sel_pval) ## change 'sel_pval' df/matrix into a vector variable
  print(sprintf("%s pvalue originally are 0", sum(pval==0)))
  pval[pval==0] = 1e-5 ## change pval=0 into 1e-5 so that after log transformation, this will be a large dot
  pval[pval> pvalue ] = NA ## pval > defined 'pvalue(0.05)' is set to NA, so that dot-plot will not plot this point
  print(sprintf("%s pvalue < %s", sum(!is.na(pval)), pvalue))
  ## setup 'plot.data' variable for plot
  plot.data = cbind(df_names,pval)
  pr = unlist(as.data.frame(sel_means))
  pr[pr==0] = 1
  plot.data = cbind(plot.data,log2(pr))
  colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')
  ## mean value of p-vaue<0.05 change into NA
  plot.data$mean[is.na(plot.data$pvalue)] = NA
  my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)

  if (selected_rows.nameReorder) {
    # print('99999999999999')
    # print(levels(factor(plot.data$pair)))
    for (r in 1:length(annotation$interacting_pair)) {
      plot.data$pair =  gsub(pattern = sprintf('^%s$', annotation$interacting_pair[r]), replacement = annotation$interacting_pair_reName[r], x = plot.data$pair)
    }
    # print(levels(factor(plot.data$pair)))
    # print('66666666666666666666')
    selected_rows = annotation$interacting_pair_reName[match(selected_rows, annotation$interacting_pair)]
  }
  
  if (yAxisCol == 'red') {
    plot.data$pair <- gsub(pattern = '_', replacement = '_              ', plot.data$pair)
    selected_rows <- gsub(pattern = '_', replacement = '_              ', selected_rows)
  } else if (yAxisCol == 'black'){
    plot.data$pair <- gsub(pattern = '_', replacement = '             ', plot.data$pair)
    selected_rows <- gsub(pattern = '_', replacement = '             ', selected_rows)
  } else if (yAxisCol == 'org'){
    plot.data$pair <- plot.data$pair
    selected_rows <- selected_rows
  } 
  
  
  if (is.null(min.mean) & is.null(max.mean)) {
    print(sprintf("Default max.mean is %s, min.mean is %s ", max(plot.data$mean, na.rm = T), min(plot.data$mean, na.rm = T)))
    
    if (yAxisCol == 'red') {
      ggplot(plot.data,aes(x=clusters,y=factor(pair, levels = selected_rows))) +
        geom_point(aes(size=-log10(pvalue),color=mean)) +
        scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) + 
        # scale_size_area() + ##make zero value to have zero size
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              # panel.grid.major = element_blank(),
              axis.text=element_text(size=14, colour = "black", face="bold"),
              axis.text.x = element_text(angle = 90, hjust = 1, face="bold"),
              # axis.text.y = element_text(size=12, colour = "black"),
              axis.text.y = element_text(size=12, colour = "red", face="bold"),
              axis.title=element_blank(),
              panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black")) + 
        theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
              plot.background = element_rect(fill = "transparent", color = NA)) # bg of the plot) 
    } else {
      ggplot(plot.data,aes(x=clusters,y=factor(pair, levels = selected_rows))) +
        geom_point(aes(size=-log10(pvalue),color=mean)) +
        scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) + 
        # scale_size_area() + ##make zero value to have zero size
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              # panel.grid.major = element_blank(),
              axis.text=element_text(size=14, colour = "black", face="bold"),
              axis.text.x = element_text(angle = 90, hjust = 1, face="bold"),
              axis.text.y = element_text(size=12, colour = "black", face="bold"),
              axis.title=element_blank(),
              panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black")) 
    }

    # g2 <- ggplot(plot.data,aes(x=clusters,y=pair)) +
    #   geom_point(aes(size=-log10(pvalue),color=mean)) +
    #   scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) + 
    #   # scale_size_area() + ##make zero value to have zero size
    #   theme_bw() + 
    #   theme(panel.grid.minor = element_blank(),
    #         panel.grid.major = element_blank(),
    #         panel.background = element_rect(fill = "transparent"),
    #         axis.text=element_text(size=14, colour = "black"),
    #         axis.text.x = element_text(angle = 90, hjust = 1),
    #         # axis.text.y = element_text(size=12, colour = "black"),
    #         axis.text.y = element_text(size=12, colour = "red", ),
    #         axis.title=element_blank(),
    #         panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black")) + 
    #   theme(
    #     panel.background = element_rect(fill = "transparent"), # bg of the panel
    #     plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    #     panel.grid.major = element_blank(), # get rid of major grid
    #     panel.grid.minor = element_blank()
    #   )  + theme(axis.text.y = element_text(size=12, colour = "red", face="bold"))
    
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
    if (yAxisCol == 'red') {
      ggplot(plot.data,aes(x=clusters,y=factor(pair, levels = selected_rows))) +
        geom_point(aes(size=-log10(pvalue),color=mean)) +
        scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette, limits = c(min.mean, max.mean) ) + 
        # scale_size_area() + ##make zero value to have zero size
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              axis.text=element_text(size=14, colour = "black", face="bold"),
              axis.text.x = element_text(angle = 90, hjust = 1, ),
              axis.text.y = element_text(size=12, colour = "red", face="bold"),
              axis.title=element_blank(),
              panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))+ 
        theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
              plot.background = element_rect(fill = "transparent", color = NA)) # bg of the plot) 
    } else {
      ggplot(plot.data,aes(x=clusters,y=factor(pair, levels = selected_rows))) +
        geom_point(aes(size=-log10(pvalue),color=mean)) +
        scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette, limits = c(min.mean, max.mean) ) + 
        # scale_size_area() + ##make zero value to have zero size
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              # panel.grid.major = element_blank(),
              axis.text=element_text(size=14, colour = "black", face="bold"),
              axis.text.x = element_text(angle = 90, hjust = 1, face="bold"),
              axis.text.y = element_text(size=12, colour = "black", face="bold"),
              axis.title=element_blank(),
              panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
    }
    
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

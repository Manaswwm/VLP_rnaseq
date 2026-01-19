#### this script is the fourth (optional) step in the microarray analysis pipeline ####
## this script serves the purpose of making the volcano plots and heatmaps 
## to visualize the DGEs

## the function will take in two values - the tophits DGE table
## and the p-value and lfc cutoff (in that order) for heatmap
step4_makeVolcanoHeatplots = function(dge_info, expr_mtx, comp){
  
  
  #@print check
  print(paste("Comparing - ", comp, sep = ""))
  
  #@print check
  print("(Step 1) - plotting number of DE genes")
  
  #extracting the number of upregulated genes - pval cutoff - 0.05 and lfc cutoff - 0.1
  upreg_genes_num = dim(dge_info[dge_info$lfc > 0 & dge_info$p_val < 0.05,])[1]
  downreg_genes_num = dim(dge_info[dge_info$lfc < 0 & dge_info$p_val < 0.05,])[1]
  
  #making a dataframe of upreg and downreg
  de_genes = data.frame(genes = c(upreg_genes_num, downreg_genes_num), type = c("Upregulated", "Downregulated"))
  
  #making a barplot capturing the number of DE up and down regulated genes
  de_genes_plot = ggplot(de_genes, aes(y = type, x = as.numeric(genes), fill = type))+
    geom_bar(position = "dodge", stat = "identity", width = 0.4)+
    ylab("Type")+xlab("Number of genes")+geom_text(aes(label = genes), position = position_dodge(width = 0.9), hjust = -0.25, size = 6)+
    scale_fill_manual("Groups", values = c("#0A0DAE", "#F57A00"))+ ggtitle("Number of DEGs (pval < 0.05)")+
    theme_bw()+theme(plot.title = element_text(size = 16, hjust = 0.5), axis.title.y = element_blank(), axis.title.x = element_text(size = 16),
                     legend.title = element_text(size = 16), legend.text = element_text(size = 14),
                     axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_text(size = 14))
  
  #saving the pca plot
  ggsave(plot = de_genes_plot, filename = paste("dge_nums/de_genes_", comp, ".jpeg", sep = ""),
         dpi = 300, height = 6, width = 15, device = "jpeg", bg = "white")

    
  #### @@ Part 2 - clustering samples through volcano plot @ ####
  
  #@print check
  print("(Step 2) - performing volcano plot to visualize the DEGs")
  
  #filtering out genes that are not significant hits
  dge_info_significant = dge_info[dge_info$p_val < 0.05,]
  
  #extracting upreg and downreg info
  upreg_genes = dge_info_significant[dge_info_significant$lfc > 0,]
  downreg_genes = dge_info_significant[dge_info_significant$lfc < 0,]
  
  #sorting and keeping the top 5 entries only
  upreg_genes = upreg_genes[order(upreg_genes$lfc, decreasing = TRUE),][1:5,]
  downreg_genes = downreg_genes[order(downreg_genes$lfc, decreasing = FALSE),][1:5,]
  
  #adding labels
  upreg_genes$trend = rep("Up-regulated")
  downreg_genes$trend = rep("Down-regulated")
  
  #merging
  dge_highlight = rbind(upreg_genes, downreg_genes)
  
  #visualizing the dge using a volcano plot
  volcano_plot = ggplot(dge_info, aes(y = -log10(p_val), x = lfc, colour = lfc > 0))+
    geom_point(size = 1, alpha = 0.6)+ geom_hline(yintercept = 1.30, linetype = "dashed", color = "black", linewidth = 1)+
    scale_color_manual(values = setNames(c("#F57A00", "#0A0DAE"), c(TRUE, FALSE)))+ theme_cowplot()+
    ylab("-log10(p_value)") + xlab("log2(fold_change)")+
    theme(axis.title.x = element_text(size = 22), axis.title.y = element_text(size = 22), legend.position = "none", 
          plot.title = element_text(size = 24, hjust = 0.5), axis.text.x =element_text(size = 20), axis.text.y = element_text(size = 20))
  
  #adding labels of top 5 up and downregulated genes
  volcano_plot = volcano_plot + geom_label_repel(data = dge_highlight, 
                                  mapping = aes(x = lfc, y = -log10(p_val), label = gene_symbol), 
                                  size = 8, max.overlaps = 10)
  
  
  #saving the ggplot to a file - I have to add a white backgroun - strange
  ggsave(plot = volcano_plot, filename = paste("volcano_plot/volcano_plot_", comp,".jpeg", sep = ""),
         dpi = 300, height = 10, width = 12, device = "jpeg", bg = "white")
  
  
  #### @@ Part 4 - plotting DEGs through heatmap @ ####

  #@print check
  print("(Step 3) - performing heatmap plot to visualize the DEGs")

  ##upregulated genes - logfc is positive
  upreg_geneinfo = dge_info[dge_info$lfc > 0 & dge_info$p_val < 0.05,] #setting p-val cutoff

  #ordering on increase in logfc
  upreg_geneinfo = upreg_geneinfo[order(upreg_geneinfo$lfc, decreasing = TRUE),]

  #taking the top 20 genes
  upreg_geneinfo = upreg_geneinfo$gene_symbol[1:15]

  ##downregulated genes - logfc is positive
  downreg_geneinfo = dge_info[dge_info$lfc < 0 & dge_info$p_val < 0.05,]  #setting p-val cutoff

  #ordering on increase in logfc
  downreg_geneinfo = downreg_geneinfo[order(downreg_geneinfo$lfc, decreasing = FALSE),]

  #taking the top 20 genes
  downreg_geneinfo = downreg_geneinfo$gene_symbol[1:15]

  #subsetting the diff_exp_genes expression matrix to contain only the top 10 up and down regulated genes
  diff_exp_genes_subset = expr_mtx[rownames(expr_mtx) %in% c(upreg_geneinfo, downreg_geneinfo),]
  #diff_exp_genes_subset = diff_exp_genes[rownames(diff_exp_genes) %in% c(upreg_geneinfo, downreg_geneinfo),]

  #clustering columns
  #we calculate pairwise distances between genes (logically this is the strength of spearman)
  #calculating distance by 1-cor so that the resulting score is in the scale of 0-2
  clust_cols = hclust(as.dist(1-cor(diff_exp_genes_subset, method = "spearman")), method = "complete")

  #clustering rows
  clust_rows = hclust(as.dist(1-cor(t(diff_exp_genes_subset), method = "pearson")), method = "complete")

  #cutting the tree - stating how many samples we have in the study - clustering the rows accordingly
  #here k is the number of samples that we have (diseased and healthy)
  module_assign = cutree(clust_rows, k=2)

  #assigning colour to each module
  #here too - the colours are influenced by the number of samples (disease vs healthy)
  #getting the colours
  module_colour = rainbow(n = length(unique(module_assign)), start = 0.1, end = 0.9)

  #making copies of the colours on the basis of the group assignment
  module_colour = module_colour[as.vector(module_assign)]

  #defining a colour pallette
  col_pal = colorRampPalette(colors = c("red", "white", "green"))(10)

  #making row annotation
  my_gene_col = data.frame(cluster = ifelse(test = module_assign == 1, yes = "upregulated", no = "downregulated"))

  #plotting heatmap
  heatmap = pheatmap(mat = diff_exp_genes_subset, legend = TRUE, col = col_pal, scale = "row", cluster_cols = FALSE, cluster_rows = TRUE)

  #saving with ggsave
  ggsave(filename = paste("heatmap/heatmap_plot_", comp,".jpg", sep = ""), plot = heatmap,
         dpi = 300, height = 10, width = 12, device = "jpeg", bg = "white")

  #closing the dev
  dev.off()

  #@print check
  print("Done making the plots - saved in the backup folder")
  
  #returning the plots back
  #return(volcano_plot)
  

}

## code dump
#### @ chekcing probes and genes with highest log fold increase @####
# upreg_genes_info = dgelist_tophits[order(dgelist_tophits$logFC, decreasing = TRUE),]
# upreg_genes_info$logPval = -log10(upreg_genes_info$P.Value)
# upreg_genes_info = upreg_genes_info[upreg_genes_info$logPval>3 & upreg_genes_info$logFC > 0.0,]
# upreg_genes_info$gene_names = mapIds(clariomshumantranscriptcluster.db, keys = upreg_genes_info$probeID, column = "SYMBOL", keytype = "PROBEID")
# upreg_genes_info$gene_names = upreg_genes_names[!is.na(upreg_genes_names)]
# chris_ids = c("GBP2", "UBE2L6", "EPSTI1", "PR2RY14", "FCGR1BP", "GBP5", "GBP1", "FCGR1" , "TRIM22")
# 
# chris_ids[chris_ids %in% upreg_genes_names]

#### check done ####

#we use hierarchical cluster - unsupervised clustering algorithm

# #clustering columns
# #we calculate pairwise distances between genes (logically this is the strength of spearman)
# #calculating distance by 1-cor so that the resulting score is in the scale of 0-2
# clust_cols = hclust(as.dist(1-cor(diff_exp_genes, method = "spearman")), method = "complete")
# 
# #clustering rows
# clust_rows = hclust(as.dist(1-cor(t(diff_exp_genes), method = "pearson")), method = "complete")
# 
# #cutting the tree - stating how many samples we have in the study - clustering the rows accordingly
# #here k is the number of samples that we have (diseased and healthy)
# module_assign = cutree(clust_rows, k=2)
# 
# #assigning colour to each module
# #here too - the colours are influenced by the number of samples (disease vs healthy)
# #getting the colours
# module_colour = rainbow(n = length(unique(module_assign)), start = 0.1, end = 0.9)
# 
# #making copies of the colours on the basis of the group assignment
# module_colour = module_colour[as.vector(module_assign)]
# 
# #defining a colour pallette
# col_pal = colorRampPalette(colors = c("red", "white", "green"))(100)
# 
# #opening the heatmap save code
# png(paste(backup_file, "heatmap_plot_DEG.png"), width = 12, height = 10, units = "in", res = 300)
# 
# #plotting the heatmap
# #x indicates the exp values of DEGs, Rowv and Colv indicates the clustering of the rows and columns
# #RowSideColours indicates the colours to the row cluster, col indicates the colours for columns
# #by scaling for row I ensure that the comparison is normalised based and not magnitude based
# #for example difference of 100-10 is given the same importance as 10-1
# heatmap.2(x = diff_exp_genes, Rowv = as.dendrogram(clust_rows), Colv = as.dendrogram(clust_cols),
#           RowSideColors = module_colour, col = bluered(80), scale = "row",
#           labRow = NA, density.info = "none", trace = "none", cexRow = 1, cexCol = 1, keysize = 2, 
#           dendrogram = "column")

# #adding labels - T2 vs T0
# dge_tophits$genelabels = ifelse(dge_tophits$gene_symbol == "FCGR1A" |
#                                 dge_tophits$gene_symbol == "GBP5"|
#                                 dge_tophits$gene_symbol == "GBP1"|
#                                 dge_tophits$gene_symbol == "FCGR1BP"|
#                                 dge_tophits$gene_symbol == "GBP2"|
#                                 dge_tophits$gene_symbol == "UBE2L6"|
#                                 dge_tophits$gene_symbol == "EPSTI1"|
#                                 dge_tophits$gene_symbol == "STAT1"|
#                                 dge_tophits$gene_symbol == "TRIM22"|
#                                 dge_tophits$gene_symbol == "P2RY14", TRUE, FALSE)

#### extracting labels of top 5 up and downregulated genes ####

#checking if dev is off
#dev.off()

# #opening the heatmap save code
# png(paste(backup_file, "heatmap_plot_DEG.png", sep = ""), width = 12, height = 10, units = "in", res = 300)

#plotting the heatmap
#x indicates the exp values of DEGs, Rowv and Colv indicates the clustering of the rows and columns
#RowSideColours indicates the colours to the row cluster, col indicates the colours for columns
#by scaling for row I ensure that the comparison is normalised based and not magnitude based
#for example difference of 100-10 is given the same importance as 10-1
# pheatmap(mat = diff_exp_genes_subset, Rowv = as.dendrogram(clust_rows), Colv = as.dendrogram(clust_cols),
#           RowSideColors = module_colour, col = col_pal, scale = "row",
#           labRow = NA, density.info = "none", trace = "none", cexRow = 1, cexCol = 1, keysize = 2, 
#           dendrogram = "column", annotation_row = my_gene_col, annotation_col = my_sample_col, labels_col = rep("", dim(my_sample_col)[1]))

#clustering columns
# pheatmap(mat = diff_exp_genes_subset, legend = TRUE,
#          RowSideColors = module_colour, col = col_pal, scale = "row", cluster_cols = TRUE, cluster_rows = TRUE,
#          labRow = NA, density.info = "none", trace = "none", cexRow = 1, keysize = 2, 
#         annotation_row = my_gene_col, annotation_col = my_sample_col, labels_col = rep("", dim(my_sample_col)[1]), cutree_rows = 2)
# 


# #vital age distribution of the sample - extracting VAG from the metainfor dataframe by using expression matching
# my_sample_col$timepoint = unlist(lapply(sub("\\_.*", "", my_sample_col$sample_id), 
#                                      function(x) {metainfo$timepoint[grep(x = metainfo$microarray_id, pattern = x)]}))
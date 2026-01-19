#### this script will create venn diagrams to compare the number of up and down regulated genes ####
#### choice of conditions is - VLP and VLPCL ####

##importing libraries
library(VennDiagram)
library(ggvenn)

##VLPCL#ggvenn#VLPCL

#reading dge information
vlpcl_dge = read.delim("input_files/misc/NS_vs_VLPCL.all.xls", header = TRUE, sep = "\t")

#removing entries that are termed "NewGene" - do not have any annotations
newgene_entries = vlpcl_dge$Symbol[grep(pattern = "NewGene", vlpcl_dge$Symbol)]
vlpcl_dge = vlpcl_dge[!vlpcl_dge$Symbol %in% newgene_entries,]

#removing duplicated entries - example AKAP17A, AMSTL etc - they mostly have similar lfc and p-value
#there are 20 genes that are duplicated - only taking unique
vlpcl_dge = vlpcl_dge[!duplicated(vlpcl_dge$Symbol),]

### splitting this into two ###

## first - dgeinfo ##
vlpcl_dge = vlpcl_dge[,c("Symbol", "Pvalue", "log2FC")]

#changing column names
colnames(vlpcl_dge) = c("gene_symbol", "p_val", "lfc")

#shortlisting genes that are differentially expressed
vlpcl_dge = vlpcl_dge[vlpcl_dge$p_val < 0.05,]

#shortlisting genes that are up and downregulated
vlpcl_upreg = vlpcl_dge$gene_symbol[vlpcl_dge$lfc > 0]
vlpcl_downreg = vlpcl_dge$gene_symbol[vlpcl_dge$lfc < 0]

##VLP

#reading dge information
vlp_dge = read.delim("input_files/misc/NS_vs_VLP.all.xls", header = TRUE, sep = "\t")

#removing entries that are termed "NewGene" - do not have any annotations
newgene_entries = vlp_dge$Symbol[grep(pattern = "NewGene", vlp_dge$Symbol)]
vlp_dge = vlp_dge[!vlp_dge$Symbol %in% newgene_entries,]

#removing duplicated entries - example AKAP17A, AMSTL etc - they mostly have similar lfc and p-value
#there are 20 genes that are duplicated - only taking unique
vlp_dge = vlp_dge[!duplicated(vlp_dge$Symbol),]

### splitting this into two ###

## first - dgeinfo ##
vlp_dge = vlp_dge[,c("Symbol", "Pvalue", "log2FC")]

#changing column names
colnames(vlp_dge) = c("gene_symbol", "p_val", "lfc")

#shortlisting genes that are differentially expressed
vlp_dge = vlp_dge[vlp_dge$p_val < 0.05,]

#shortlisting genes that are up and downregulated
vlp_upreg = vlp_dge$gene_symbol[vlp_dge$lfc > 0]
vlp_downreg = vlp_dge$gene_symbol[vlp_dge$lfc < 0]


#### making venn diagrams ####

#making lists for venn_diagram
upreg_genes = list(VLP_Upreg = vlp_upreg, VLPCL_Upreg = vlpcl_upreg)
downreg_genes = list(VLP_Downreg = vlp_downreg, VLPCL_Downreg = vlpcl_downreg)
total_genes = list(VLP_Total = c(vlp_upreg, vlp_downreg), VLPCL_Total = c(vlpcl_upreg, vlpcl_downreg))

## upreg diagram
upreg_venn = ggvenn(upreg_genes)
ggsave(filename = "venn_plots/upreg_venn.jpg", plot = upreg_venn, width = 8, height = 6, dpi = 300, device = "jpeg")

## downreg diagram
downreg_venn = ggvenn(downreg_genes)
ggsave(filename = "venn_plots/downreg_venn.jpg", plot = downreg_venn, width = 8, height = 6, dpi = 300, device = "jpeg")

## upreg diagram
total_venn = ggvenn(total_genes)
ggsave(filename = "venn_plots/total_venn.jpg", plot = total_venn, width = 8, height = 6, dpi = 300, device = "jpeg")

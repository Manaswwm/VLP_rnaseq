#### this script will analyse the dge information passed to me by Shuran for her project ####
## the goal of this script is to give primary dge outputs (heatmap and volcano plot) and tmod analysis ##
## the input information is in the form equivalent to dge_tophits dataframe and expr_mtx matrix ##

##importing relevant packages
library(ggplot2)
library(tidyverse)
library(tmod)
library(gplots)
library(plotly)
library(pheatmap)
library(readxl)
library(stringr)
library(cowplot)
library(ggrepel)

##importing relevant functions
source("step4_makeVolcanoHeatplots.R")
source("step5_tmodDGEmods.R")

#declaring the file paths seperately for dge_tophits and expr_mtx
dge_tophits_path = list.files(path = "input_files/misc/", full.names = TRUE)

#declaring the comparisons
compare = c("CL_vs_VLPCL.all", "NS_vs_CL.all", "NS_vs_VLP.all", "VLP_vs_VLPCL.all", "NS_vs_VLPCL.all")

#declaring a function that will perform the tasks and create outputs
for(comp in compare){
  
  #extracting the dge_tophits and expr_mtx associated with the comparison to be made
  dge_tophits = dge_tophits_path[grep(pattern = comp, dge_tophits_path)]
  
  #loading the dge_tophits
  dge_tophits = read.delim(file = dge_tophits, header = TRUE, sep = "\t")
  
  #removing .all from the string names
  comp = str_remove(comp, pattern = ".all")
  
  #removing entries that are termed "NewGene" - do not have any annotations
  newgene_entries = dge_tophits$Symbol[grep(pattern = "NewGene", dge_tophits$Symbol)]
  dge_tophits = dge_tophits[!dge_tophits$Symbol %in% newgene_entries,]
  
  #removing duplicated entries - example AKAP17A, AMSTL etc - they mostly have similar lfc and p-value
  #there are 20 genes that are duplicated - only taking unique
  dge_tophits = dge_tophits[!duplicated(dge_tophits$Symbol),]
  
  ### splitting this into two ###
  
  ## first - dgeinfo ##
  dge_info = dge_tophits[,c("Symbol", "Pvalue", "log2FC")]
  
  #changing column names
  colnames(dge_info) = c("gene_symbol", "p_val", "lfc")
  
  ## second - expression matrix ##
  expr_mtx = dge_tophits[,seq(9,14,by=1)]
  rownames(expr_mtx) = dge_tophits$Symbol
  
  #removing FPKM from the column IDs
  colnames(expr_mtx) = str_remove(string = colnames(expr_mtx), pattern = "_FPKM")
  
  ### determining column order ###
  ## dirty trick
  cond1 = strsplit(comp, "_")[[1]][1]
  cond2 = strsplit(comp, "_")[[1]][3]
  
  #setting colnames order
  col_order = c(paste(cond1,"4", sep = ""), paste(cond1,"5", sep = ""), paste(cond1,"6", sep = ""), 
                paste(cond2,"4", sep = ""), paste(cond2,"5", sep = ""), paste(cond2,"6", sep = ""))
  
  #setting expr_mtx column names order
  expr_mtx = expr_mtx[,col_order]
  
  #### performing initial de genes, volcano plot and heatmap analyses ####
  step4_makeVolcanoHeatplots(dge_info = dge_info, expr_mtx = expr_mtx, comp = comp)
  
  #### performing tmod analysis ####
  step5_tmodDGEmods(dge_info = dge_info)

}
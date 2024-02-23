library('shiny')
library(scClustViz)
library(Seurat)

transplant_Integrated_scCLustViz_object <- "~/LiverTransplant/Objects/transplant_Integrated_scCLustViz_object.RData"
load(transplant_Integrated_scCLustViz_object)


transplant_data$cluster = as.character(sCVdata_list$RNA_snn_res.0.2@Clusters)
Idents(transplant_data) <- paste0('cluster_', as.character(transplant_data$cluster))
cluster_names <-  levels(transplant_data)
### finding the markers while removing the mito-genes
Cluster_markers <- sapply(1:length(cluster_names), 
                          function(i) FindMarkers(transplant_data, 
                                                  ident.1=cluster_names[i]
                                                  #logfc.threshold = 0,
                                                  #min.pct=0,
                                                  #min.cells.feature = 1,
                                                  #min.cells.group = 1
                          ), 
                          simplify = FALSE)
names(Cluster_markers) <- cluster_names


P_value_thr = 0.02
Cluster_markers_final <- sapply(1:length(Cluster_markers), function(i) {
  
  ## selecting the cluster of interest's marker dataframe (x)
  x = Cluster_markers[[i]]
  a_cluster_name <- names(Cluster_markers)[i]
  
  ## sort rows of x based on log-fold change
  x = x[order(x$avg_log2FC, decreasing = T),]
  
  ## sort based on adj p-value and sign of logFC  
  # x$ranking_score=-log10(x$p_val_adj+.Machine$double.xmin)*sign(x$avg_log2FC)
  # x = x[order(x$ranking_score, decreasing = T),]
  
  ## add the average expression of each gene as a column
  selected_cells = Idents(transplant_data) == a_cluster_name
  data <- GetAssayData(transplant_data)[rownames(x),]
  x$avg_exp <- rowSums(data[,selected_cells])/sum(selected_cells)
  
  ## filtering genes with adj-p-value higher than 0.05
  #x = x[x$p_val_adj<P_value_thr,]
  
  return(x)
}, simplify = F)

names(Cluster_markers_final) <- names(Cluster_markers)
#lapply(Cluster_markers_final, head)
#lapply(Cluster_markers_final, dim)

#Cluster_markers <- readRDS('Results/old_samples/Cluster_markers_mergedOldSamples.rds')
#saveRDS(Cluster_markers_final, 'Results/new_samples/Cluster_markers_final.rds')
saveRDS(Cluster_markers_final, '~/LiverTransplant/Results/merged_transplant_june26_Cluster_markers_res.0.2.rds')


#### saving markers #####
dir.create(paste0('~/LiverTransplant/Results/final_markers_merged_transplant_res.0.2/'))
for(i in 1:length(Cluster_markers_final)){
  #df <- data.frame(genes=rownames(Cluster_markers[[i]]),
  #                 score=Cluster_markers[[i]]$ranking_score)
  df <- data.frame(Cluster_markers_final[[i]])
  print(head(df, 25))
  file_name <- names(Cluster_markers_final)[i]
  write.csv(df, paste0('~/LiverTransplant/Results/final_markers_merged_transplant_res.0.2/final_markers_transplant_', file_name,'.txt'), 
            col.names = T, row.names = T, quote = F)
}


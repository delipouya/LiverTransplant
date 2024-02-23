library(plyr)
library(stats)
library(ggpubr)

library(RColorBrewer)
library(viridis)
library(scales)

source('~/RatLiver/Codes/Functions.R')
Initialize()

############ Importing the new data to be annotated #########
Syngeneic = readRDS('Objects/Integrated_Syngeneic_mt15_lib1500_harmony.rds')
Rejection = readRDS('Objects/Integrated_Rejection_mt15_lib1500_harmony.rds')
Tolerance = readRDS('Objects/Integrated_Tolerance_mt15_lib1500_harmony.rds')

annot.syn = readRDS('Objects/Integrated_Syngeneic_res0.5_set1_based_annotation.rds')
annot.rej = readRDS('Objects/Integrated_Rejection_res0.3_set2_based_annotation.rds')
annot.tol = readRDS('Objects/Integrated_Tolerance_res0.2_set2_based_annotation.rds')

transplant_data <- merge(Syngeneic, c(Rejection, Tolerance), # 
                  add.cell.ids = c('Syngeneic', 'Rejection', 'Tolerance'), 
                  project = 'All', 
                  merge.data = TRUE)


transplant_data$state = sapply(str_split(colnames(transplant_data), '_'), '[[', 1)
table(transplant_data$POD)
table(transplant_data$sample_name)
table(transplant_data$state)
table(transplant_data$ID)

transplant_data =  FindVariableFeatures(transplant_data)
transplant_data <- ScaleData(transplant_data)
transplant_data <- RunPCA(transplant_data, verbose = FALSE) #all_common_features for varimax
plot(100 * transplant_data@reductions$pca@stdev^2 / transplant_data@reductions$pca@misc$total.variance,
     pch=20,xlab="Principal Component",ylab="% variance explained",log="y")
PC_NUMBER = 30

transplant_data <- RunHarmony(transplant_data, group.by.vars = "state", assay.use="RNA")

transplant_data <- RunUMAP(transplant_data, reduction = "harmony", dims = 1:PC_NUMBER)

df_umap <- data.frame(UMAP_1=getEmb(transplant_data, 'umap')[,1], 
                      UMAP_2=getEmb(transplant_data, 'umap')[,2], 
                      clusters=transplant_data$cluster,
                      POD=transplant_data$POD,
                      state=transplant_data$state,
                      sample=transplant_data$sample_name,
                      umi=colnames(transplant_data))


ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=POD))+geom_point(alpha=0.3)+theme_bw()
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=state))+geom_point(alpha=0.3)+theme_bw()
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=sample))+geom_point(alpha=0.3)+theme_bw()
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=clusters))+geom_point(alpha=0.3)+theme_bw()


PC_NUMBER = 30
transplant_data = readRDS('~/LiverTransplant//Objects/Integrated_all.rds')

transplant_data@meta.data = transplant_data@meta.data[,-grep('RNA_snn_', colnames(transplant_data@meta.data))]
transplant_data@meta.data = transplant_data@meta.data[,-grep('cluster', colnames(transplant_data@meta.data))]

transplant_data <- FindNeighbors(transplant_data, 
                                  reduction = "harmony", 
                                  dims = 1:PC_NUMBER)

# sCVdata_list <- get_cluster_object(seur = transplant_data2, 
#                                    max_seurat_resolution = 2, ## change this to higher values
#                                    FDRthresh = 0.05, # FDR threshold for statistical tests
#                                    min_num_DE = 10,
#                                    seurat_resolution = 1, # Starting resolution is this plus the jump value below.
#                                    seurat_resolution_jump = 0.5,
#                                    DE_bw_clust = TRUE)

#saveRDS(sCVdata_list, 'Objects/sCVdata_list_Integrated_all.rds') ### find the results on run1

resolutions = c(0.2, 0.5, 0.7, 1)
for (i in 1:length(resolutions)){
  res = resolutions[i]
  transplant_data <- FindClusters(transplant_data, resolution = res)
}


## creating the meta.data dataframe
head(transplant_data@meta.data)
your_cluster_results <- data.frame(transplant_data@meta.data[,grep('RNA_snn_', colnames(transplant_data@meta.data))]) 
rownames(your_cluster_results)  = rownames(transplant_data@meta.data)
head(your_cluster_results)
your_cluster_results = readRDS('~/LiverTransplant/Objects/Integrated_all_cluster_df.rds')

### calculating the differentially expressed marker genes
sCVdata_list <- CalcAllSCV(
  inD=transplant_data,
  clusterDF=your_cluster_results,
  assayType='RNA', #specify assay slot of data
  DRforClust="harmony",#reduced dimensions for silhouette calc
  exponent=exp(1), #log base of normalized data
  pseudocount=1,
  DRthresh=0.3, #gene filter - minimum detection rate
  testAll=F, #stop testing clusterings when no DE between clusters
  FDRthresh=0.01,
  calcSil=F, #use cluster::silhouette to calc silhouette widths
  calcDEvsRest=T,
  calcDEcombn=F
)

transplantation_data_scCLustViz = 'Objects/Integrated_all_scClustViz.RData'
save(transplant_data2,sCVdata_list,file=transplantation_data_scCLustViz) ## new data scClustViz object


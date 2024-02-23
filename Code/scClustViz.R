library(plyr)
library(stats)
library(ggpubr)
library(RColorBrewer)
library(viridis)
library(scales)
source('~/RatLiver/Codes/Functions.R')
Initialize()

Syngeneic = readRDS('Objects/Integrated_Syngeneic_mt15_lib1500_harmony.rds')
Rejection = readRDS('Objects/Integrated_Rejection_mt15_lib1500_harmony.rds') ## checksum-3 killed
Tolerance = readRDS('Objects/Integrated_Tolerance_mt15_lib1500_harmony.rds')

annot.syn = readRDS('Objects/Integrated_Syngeneic_res0.5_set1_based_annotation.rds')
annot.rej = readRDS('Objects/Integrated_Rejection_res0.3_set2_based_annotation.rds')
annot.tol = readRDS('Objects/Integrated_Tolerance_res0.2_set2_based_annotation.rds')

sample_name = 'POD3'
transplant_data = readRDS('~/LiverTransplant/Objects/Integrated_POD3_mt3_lib1500_harmony.rds')



transplant_data = Rejection # Syngeneic Rejection
sample_name = 'Rejection' # 'Rejection' 'Tolerance' Syngeneic
sample_state = 'R' # S T

transplant_data <- FindNeighbors(transplant_data, reduction = "harmony", dims = 1:30)

resolution = 0.5
if(sample_state=='R') resolution = 0.3
if(sample_state=='T') resolution = 0.6

resolution = 0.8
transplant_data <- FindClusters(transplant_data, resolution = resolution)


#transplant_data = readRDS('~/LiverTransplant//Objects/Integrated_all.rds')
#your_cluster_results = readRDS('~/LiverTransplant/Objects/Integrated_all_cluster_df.rds')
#your_cluster_results = your_cluster_results[,c(1,2)]

#sample_indices = sample(1:ncol(transplant_data), size = 1000)
#transplant_data_sub = transplant_data[,sample_indices]
#your_cluster_results_sub = your_cluster_results[sample_indices,]

# saveRDS(list(transplant_data_sub, your_cluster_results_sub), '~/LiverTransplant/transplant_data_sub.rds')
#input_list = readRDS('~/LiverTransplant/transplant_data_sub.rds')
#transplant_data_sub = input_list[[1]]
#your_cluster_results_sub = input_list[[2]]

your_cluster_results = data.frame(res.0.3=transplant_data$RNA_snn_res.0.3,
                                    res.0.5=transplant_data$RNA_snn_res.0.5,
                                    res.0.8=transplant_data$RNA_snn_res.0.8)
### calculating the differentially expressed marker genes

#### The altered version used for the Rejection sample
sCVdata_list <- CalcAllSCV(
  inD=transplant_data,
  clusterDF=your_cluster_results,
  assayType='RNA', #specify assay slot of data
  DRforClust="harmony",#reduced dimensions for silhouette calc
  exponent=exp(1), #log base of normalized data
  pseudocount=1,
  DRthresh=0.5, #gene filter - minimum detection rate
  testAll=T, #stop testing clusterings when no DE between clusters
  FDRthresh=0.005,
  calcSil=F, #use cluster::silhouette to calc silhouette widths
  calcDEvsRest=T,
  calcDEcombn= F #
)


sCVdata_list <- CalcAllSCV(
  inD=transplant_data,
  clusterDF=your_cluster_results,
  assayType='RNA', #specify assay slot of data
  DRforClust="harmony",#reduced dimensions for silhouette calc
  exponent=exp(1), #log base of normalized data
  pseudocount=1,
  DRthresh=0.3, #gene filter - minimum detection rate
  testAll=T, #stop testing clusterings when no DE between clusters
  FDRthresh=0.01,
  calcSil=F, #use cluster::silhouette to calc silhouette widths
  calcDEvsRest=T,
  calcDEcombn=T #
)

.libPaths( c(.libPaths(), "/home/delaram/R/x86_64-pc-linux-gnu-library/4.2" ,
                             "/usr/local/lib/R/site-library",
                             "/usr/lib/R/site-library",
                             "/usr/lib/R/library") ) 
        
#saveRDS(sCVdata_list, '~/LiverTransplant/Objects/sCVdata_list_Integrated_all.rds') ### find the results on run1
#transplant_Integrated_scCLustViz_object <- "~/LiverTransplant/Objects/transplant_Integrated_scCLustViz_object.RData"
sample_name = 'POD3'
transplant_data_scCLustViz_object <- paste0("~/LiverTransplant/Objects/transplant_",sample_name,"_scCLustViz_object.RData")


save(transplant_data,sCVdata_list,
     file=transplant_data_scCLustViz_object) ## new data scClustViz object




load(transplant_data_scCLustViz_object)

runShiny(
  ## write the path to the file bellow:
  filePath= transplant_data_scCLustViz_object,
  outPath="./",
  # Save any further analysis performed in the app to the
  # working directory rather than library directory.
  annotationDB="org.Rn.eg.db",
  # This is an optional argument, but will add annotations.
  imageFileType="png"
  #Set the file format of any saved figures from the app.
)




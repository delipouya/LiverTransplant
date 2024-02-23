# In this script, we will read the 10X data files into session and
# we will filter the data based on the fraction of mitochondrial transcripts and 
# also library size, then we will scale and normalize the data and
# perform dimension reduction. Finally, we'll visualize the data using tsne and umap

# Tallulah: 50% mt threshold across all my samples plus >=250 detected genes


get_folder_name <- function(input_type, input_POD) {
  path_lists = list.files(paste0('Data/' , input_type))
  path_lists = path_lists[grepl(pattern = input_POD, path_lists)]
  return(path_lists)}

### loading the required libraries
source('../RatLiver/Codes/Functions.R')
Initialize()

## Define cut-off values
MIT_CUT_OFF = 5
LIB_SIZE_CUT_OFF = 1500
NUM_GENES_DETECTED = 250

####### merging samples based on POD
input_POD = 'POD7'#'POD3'
path_10x_folders_syn = list.files(paste0('Data/Syngeneic/' ),full.names = T)
path_10x_folders_rej = list.files(paste0('Data/Rejection/' ),full.names = T)
path_10x_folders_tol = list.files(paste0('Data/Tolerance/' ),full.names = T)
path_10x_folders <- c(path_10x_folders_rej, path_10x_folders_syn, path_10x_folders_tol)
path_10x_folders = path_10x_folders[grepl(input_POD, path_10x_folders)]

#folder_names = get_folder_name(input_type= 'Syngeneic', input_POD)

sample_names_POD3 = c( 'R003_POD3_RAT_CST_3pr_v3',  'R009_POD3_CST',
                 'XC116_POD3_RAT_CST_3pr_v3', 'XC130_POD3_CST',
                 'TOL01_POD3_CST_3pr_v3', 'TOL02_POD3_CST_3pr_v3')

sample_names_POD7= c('R001_POD7_RAT_CST_3pr_v3', 'R005_POD7_RAT_CST_3pr_v3',
                     'XC121_POD7_RAT_CST_3pr_v3', 'XC125_POD7_RAT_CST_3pr_v3',
                     'TOL03_POD7_CST_3pr_v3', 'TOL04_POD7_CST_3pr_v3')

outcome_type = c('Rejection', 'Rejection', 'Syngeneic', 'Syngeneic', 'Tolerance', 'Tolerance')
sample_names= ''
if (input_POD=='POD3') sample_names = sample_names_POD3
if (input_POD=='POD7') sample_names = sample_names_POD7
sample_names
#### naming the samples
sample_PODs = ifelse(grepl(pattern = 'POD3',path_10x_folders), 'POD3', 'POD7')

seur_raw_list = lapply(path_10x_folders, 
                       function(input_from_10x) CreateSeuratObject(counts=Read10X(input_from_10x, gene.column = 2),
                                                                   min.cells=0,min.features=1, project = input_POD))

seur_raw_list = sapply(1:length(seur_raw_list), 
                       function(i) {
                         seur_raw_list[[i]]$state = outcome_type[i]
                         seur_raw_list[[i]]$POD = sample_PODs[i]
                         seur_raw_list[[i]]$sample_name = sample_names[i]
                         seur_raw_list[[i]]$ID = paste0(sample_names, '_', sample_PODs)[i]
                         return(seur_raw_list[[i]])
                       },simplify = F)

names(seur_raw_list) = sample_names #paste0(sample_names, '_', sample_PODs)
head(seur_raw_list[[1]])
lapply(seur_raw_list, ncol)

seur_raw <- merge(seur_raw_list[[1]], c(seur_raw_list[[2]], seur_raw_list[[3]], 
                                        seur_raw_list[[4]], seur_raw_list[[5]], seur_raw_list[[6]] ), # 
                  add.cell.ids = names(seur_raw_list), 
                  project = names(seur_raw_list), 
                  merge.data = TRUE)
dim(seur_raw)
dim(seur_raw_list[[4]])
Idents(seur_raw) = seur_raw$ID

################ Choosing the QC metric thresholds #############################
i = 1
input_from_10x = path_10x_folders[i]
INPUT_NAME = input_POD

seur_genes_df <- read.delim(paste0(input_from_10x,'/features.tsv.gz'), header = F)
seur_raw[['RNA']] <- AddMetaData(seur_raw[['RNA']], seur_genes_df$V2, col.name = 'symbol')
seur_raw[['RNA']] <- AddMetaData(seur_raw[['RNA']], seur_genes_df$V1, col.name = 'ensembl')
libSize <- colSums(GetAssayData(seur_raw@assays$RNA))

##### Mt-percentage calculation ##### 
MIT_PATTERN = '^Mt-'
mito_genes_index = grep(pattern = MIT_PATTERN, rownames(seur_raw) )
rownames(seur_raw)[mito_genes_index]

seur_mt_perc <- PercentageFeatureSet(seur_raw, features = mito_genes_index)
summary(seur_mt_perc)

manual_mt_perc = colSums(seur_raw[mito_genes_index,])*100/colSums(seur_raw)
summary(manual_mt_perc)

seur_raw$mito_perc <- manual_mt_perc
####################################  

print(paste0('Total number of cells: ', ncol(seur_raw)))

to_drop_mito <- seur_raw$mito_perc > MIT_CUT_OFF

print(paste0('to_drop_mito: ',sum(to_drop_mito)))
print(paste0('to_drop_mito percentage: ', round(sum(to_drop_mito)*100/ncol(seur_raw),2) ))

LIB_SIZE_CUT_OFF_MAX = 60000
to_drop_lib_size <- seur_raw$nCount_RNA < LIB_SIZE_CUT_OFF | seur_raw$nCount_RNA > LIB_SIZE_CUT_OFF_MAX
print(paste0('to_drop_lib_size: ', sum(to_drop_lib_size)))
print(paste0('to_drop_lib_size percentage: ', round( sum(to_drop_lib_size)*100/ncol(seur_raw),2) ))

print(paste0('remained after both filters: ', sum(!to_drop_lib_size & !to_drop_mito)))
print(paste0('Percentage of remained after both filters: ', 
             round(  sum(!to_drop_lib_size & !to_drop_mito)*100/ncol(seur_raw),2) ))



to_drop_num_genes <- seur_raw$nFeature_RNA < NUM_GENES_DETECTED
print(paste0('to_drop_num_genes: ', sum(to_drop_num_genes)))
print(paste0('to_drop_num_genes percentage: ', round( sum(to_drop_num_genes)*100/ncol(seur_raw),2) ))

print(paste0('remained after both filters: ', sum(!to_drop_num_genes & !to_drop_mito)))
print(paste0('Percentage of remained after both filters: ', 
             round(  sum(!to_drop_num_genes & !to_drop_mito)*100/ncol(seur_raw),2) ))


df = data.frame(library_size= seur_raw$nCount_RNA, 
                mito_perc=seur_raw$mito_perc , 
                n_expressed=seur_raw$nFeature_RNA,
                POD=seur_raw$POD,
                ID=seur_raw$ID,
                sample_name=seur_raw$sample_name)


#dir.create(paste0('Plots/QC/',INPUT_NAME))
## Visualization of QC metrics
pdf(paste0('Plot/',INPUT_NAME,'/QC_',INPUT_NAME,'_',
           'mito_',MIT_CUT_OFF,'_lib_',LIB_SIZE_CUT_OFF,'.pdf'))

ggplot(data.frame(seur_raw$nCount_RNA), aes(seur_raw.nCount_RNA))+
  geom_histogram(bins = 60,color='black',fill='pink',alpha=0.5)+
  theme_bw()+ggtitle('library size for all cells (before filter)')+xlab('Library sizes')+
  ylab('Number of cells')+labs(caption = INPUT_NAME)

ggplot(data.frame(seur_raw$nFeature_RNA), aes(seur_raw.nFeature_RNA))+
  geom_histogram(bins = 60,color='black',fill='blue',alpha=0.3)+
  theme_bw()+ggtitle('# expressed genes for all cells(before filtering)')+xlab('Number of expressed genes')+
  ylab('Number of cells')+labs(caption = INPUT_NAME)

ggplot(data.frame(seur_raw$mito_perc), aes(seur_raw.mito_perc))+
  geom_histogram(bins = 60,color='black',fill='green',alpha=0.3)+
  theme_bw()+ggtitle('proportion of reads mapped to Mt genes(before filtering)')+xlab('Mitochondrial proportion (%)')+
  ylab('Number of cells')+labs(caption = INPUT_NAME)

ggplot(df, aes(x=library_size, y=mito_perc, color=n_expressed))+geom_point()+scale_color_viridis()+
  theme_bw()+xlab('Library Size')+ylab('Mitochondrial transcript percent')+
  geom_hline(yintercept= MIT_CUT_OFF, linetype="dashed", color = "red")+labs(caption = INPUT_NAME)+
  geom_vline(xintercept = LIB_SIZE_CUT_OFF, linetype="dashed", color = "red3", size=0.5)+
  ggtitle(paste0('mito threshold: ', MIT_CUT_OFF,' , library size threshold: ', LIB_SIZE_CUT_OFF, ' (before filter)'))

head(df)
ggplot(df, aes(x=library_size, y=mito_perc, color=POD))+geom_point(size=2,alpha=0.7)+
  theme_classic()+xlab('Library size')+ylab('Mitochondrial transcript percent')+
  geom_hline(yintercept= MIT_CUT_OFF, linetype="dashed", color = "red")+labs(caption = INPUT_NAME)+
  geom_vline(xintercept = LIB_SIZE_CUT_OFF, linetype="dashed", color = "red3", size=0.5)+
  ggtitle('Library size and mitochondrial transcript cutoffs')

ggplot(df, aes(x=library_size, y=mito_perc, color=ID))+geom_point(size=2,alpha=0.4)+
  theme_classic()+xlab('Library size')+ylab('Mitochondrial transcript percent')+
  geom_hline(yintercept= MIT_CUT_OFF, linetype="dashed", color = "red")+labs(caption = INPUT_NAME)+
  geom_vline(xintercept = LIB_SIZE_CUT_OFF, linetype="dashed", color = "red3", size=0.5)+
  ggtitle('Library size and mitochondrial transcript cutoffs')

ggplot(df, aes(x=n_expressed, y=mito_perc, color=library_size))+geom_point()+scale_color_viridis()+
  theme_bw()+xlab('num detected genes')+ylab('Mitochondrial transcript percent')+
  geom_hline(yintercept= MIT_CUT_OFF, linetype="dashed", color = "red")+labs(caption = INPUT_NAME)+
  geom_vline(xintercept = NUM_GENES_DETECTED, linetype="dashed", color = "red3", size=0.5)+
  ggtitle(paste0('mito threshold: ', MIT_CUT_OFF,' , num detected genes threshold: ', NUM_GENES_DETECTED, ' (before filter)'))


ggplot(df, aes(x=library_size, y=n_expressed, color=mito_perc))+geom_point()+labs(caption = INPUT_NAME)+
  theme_bw()+xlab('Library Size')+ylab('Number of expressed genes')+scale_color_viridis(option = 'magma')+ggtitle('before filter')

ggplot(df, aes(x=library_size, y=n_expressed))+geom_point(color='darkblue')+theme_bw()+xlab('Library Size')+
  ylab('Number of expressed genes')+geom_point(data=df[to_drop_mito,],pch=4,color="red")+labs(caption = INPUT_NAME)+
  scale_color_viridis(option='magma', direction = 1)+ggtitle('before filter, labels: high-mito cells')+labs(caption = INPUT_NAME)

to_drop_mito <- seur_raw$mito_perc > MIT_CUT_OFF
to_drop_lib_size <- seur_raw$nCount_RNA < LIB_SIZE_CUT_OFF | seur_raw$nCount_RNA > LIB_SIZE_CUT_OFF_MAX
seur <- seur_raw[,!to_drop_mito & !to_drop_lib_size & !seur_raw$nFeature_RNA<NUM_GENES_DETECTED]

df_filt = data.frame(library_size= seur$nCount_RNA, mito_perc=seur$mito_perc , n_expressed=seur$nFeature_RNA)
ggplot(df_filt, aes(x=library_size, y=mito_perc, color=library_size))+geom_point()+scale_color_viridis()+
  theme_bw()+xlab('Library Size')+ylab('Mitochondrial transcript percent')+labs(caption = INPUT_NAME)+
  ggtitle(paste0('mito threshold: ', MIT_CUT_OFF,' , library size threshold: ', LIB_SIZE_CUT_OFF, ' (after filter)'))

ggplot(df_filt, aes(x=n_expressed, y=mito_perc, color=library_size))+geom_point()+scale_color_viridis()+
  theme_bw()+xlab('Number of detected genes')+ylab('Mitochondrial transcript percent')+labs(caption = INPUT_NAME)+
  ggtitle(paste0('mito threshold: ', MIT_CUT_OFF,' , num detected genes threshold: ', LIB_SIZE_CUT_OFF, ' (after filter)'))

################################################################################
############## Apply the QC to each element of the list #################

qc_df = data.frame(sample=names(seur_raw_list),
                   filt_frac=rep(0, length(seur_raw_list)))

MIT_CUT_OFF = 3
LIB_SIZE_CUT_OFF = 1500
NUM_GENES_DETECTED = 250


#x = seur_raw_list[[1]]
seur_raw_l <- sapply(1:length(seur_raw_list), FUN = function(i) {
  x = seur_raw_list[[i]]
  ### Mt calculation for each dataset
  MIT_PATTERN = '^Mt-'
  mito_genes_index = grep(pattern = MIT_PATTERN, rownames(x) )
  x$mito_perc = colSums(x[mito_genes_index,])*100/colSums(x)
  
  to_drop_mito <- x$mito_perc > MIT_CUT_OFF
  to_drop_lib_size <- x$nCount_RNA < LIB_SIZE_CUT_OFF | x$nCount_RNA > LIB_SIZE_CUT_OFF_MAX
  num_cells_1 = ncol(x)
  x <- x[,!to_drop_mito & !to_drop_lib_size & !x$nFeature_RNA<NUM_GENES_DETECTED]
  print(paste0('Percentage of cells removed: ', (num_cells_1 - ncol(x))/num_cells_1))
  qc_df$filt_frac[i] <<- (num_cells_1 - ncol(x))/num_cells_1
  
  return(x)
}, simplify = F)

qc_df$num_cells_after = sapply(seur_raw_l, function(x) dim(x)[2])
dev.off()
gridExtra::grid.table(qc_df)
dev.off()

all_common_features = Reduce(intersect, lapply(seur_raw_l, rownames))
length(all_common_features)
################################################################################
# normalize and identify variable features for each dataset independently
seur_norm_list <- lapply(X = seur_raw_l, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = seur_norm_list)
length(features)

################################################################################
################## Merging the datasets without/with Harmony correction ########################

concat_samples <- merge(seur_norm_list[[1]], c(seur_norm_list[[2]], seur_norm_list[[3]], 
                                               seur_norm_list[[4]],seur_norm_list[[5]] ,seur_norm_list[[6]]  ), # 
                        add.cell.ids = names(seur_raw_list), 
                        project = names(seur_raw_list), 
                        merge.data = TRUE)

concat_samples <- ScaleData(concat_samples)
VariableFeatures(concat_samples) = features

concat_samples <- RunPCA(concat_samples, npcs = 30, verbose = FALSE, features = features) #all_common_features for varimax
#saveRDS(concat_samples,  paste0('Objects/Integrated_',INPUT_NAME,'_mt' , MIT_CUT_OFF,'_lib', LIB_SIZE_CUT_OFF,'_all_features.rds'))
#saveRDS(concat_samples,  paste0('Objects/Integrated_',INPUT_NAME,'_mt' , MIT_CUT_OFF,'_lib', LIB_SIZE_CUT_OFF,'.rds'))
'~/LiverTransplant/Objects/Integrated_Syngeneic_mt3_lib1500.rds'


####################################
MIT_CUT_OFF = 3
LIB_SIZE_CUT_OFF = 1500
NUM_GENES_DETECTED = 250
INPUT_NAME = 'POD7'#'POD7'

concat_samples <- readRDS(paste0('Objects/Integrated_',INPUT_NAME,'_mt' , MIT_CUT_OFF,'_lib', LIB_SIZE_CUT_OFF,'.rds'))
head(concat_samples)
concat_samples <- RunHarmony(concat_samples, group.by.vars = "ID", assay.use="RNA")

plot(100 * concat_samples@reductions$pca@stdev^2 / concat_samples@reductions$pca@misc$total.variance,
     pch=20,xlab="Principal Component",ylab="% variance explained",log="y")

#concat_samples <- RunUMAP(concat_samples, reduction = "pca", dims = 1:30)
concat_samples <- RunUMAP(concat_samples, reduction = "harmony", dims = 1:30)

#concat_samples <- FindNeighbors(concat_samples, reduction = "pca", dims = 1:30)
concat_samples <- FindNeighbors(concat_samples, reduction = "harmony", dims = 1:30)

resolution = 0.5
if(sample_state=='R') resolution = 0.3
if(sample_state=='T') resolution = 0.6
concat_samples <- FindClusters(concat_samples, resolution = resolution)


TITLE_umap = INPUT_NAME
df_umap <- data.frame(UMAP_1=getEmb(concat_samples, 'umap')[,1], 
                      UMAP_2=getEmb(concat_samples, 'umap')[,2], 
                      library_size= concat_samples$nCount_RNA, 
                      mito_perc=concat_samples$mito_perc, 
                      n_expressed=concat_samples$nFeature_RNA,
                      POD=concat_samples$POD,
                      ID= concat_samples$ID,
                      cluster=concat_samples$RNA_snn_res.0.5,
                      sample_name=concat_samples$sample_name)

markers = c('Alb', 'Apoa2', 'Apoc1', 'Sparc', 'Lyve1', 'Col3a1', 'Ecm1', 
            'Ptprc', 'Marco', 'Lyz2',	'Vsig4', 'Cd74', 'Cd79b', 
            'Irf8', 'Plac8', 'Gzmk', 'Gzma')

i = 8
gene_name = markers[i]
df_umap$gene = GetAssayData(concat_samples)[gene_name,] #, assay='RNA'
#df_umap$gene = GetAssayData(merged_samples)[gene_name,] 
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=gene))+geom_point(alpha=0.5)+theme_bw()+
  scale_color_viridis(direction = -1)+ggtitle(gene_name)+labs(caption = INPUT_NAME)

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=library_size))+geom_point(alpha=0.4)+theme_bw()+
  scale_color_viridis(direction = -1)+ggtitle(TITLE_umap)+labs(caption = INPUT_NAME)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=mito_perc))+geom_point(alpha=0.4)+theme_bw()+
  scale_color_viridis(direction = -1)+ggtitle(TITLE_umap)+labs(caption = INPUT_NAME)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=n_expressed))+geom_point(alpha=0.4)+theme_bw()+
  scale_color_viridis(direction = -1)+ggtitle(TITLE_umap)+labs(caption = INPUT_NAME)

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=POD))+geom_point(alpha=0.3)+theme_bw()+
  ggtitle(TITLE_umap)+labs(caption = INPUT_NAME)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=ID))+geom_point(alpha=0.3)+theme_bw()+
  ggtitle(TITLE_umap)+labs(caption = INPUT_NAME)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=cluster))+geom_point(alpha=0.4)+theme_bw()+
  ggtitle(paste0('resolution=', resolution))+labs(caption = INPUT_NAME)

merged_samples = concat_samples
rm(concat_samples)
rm(seur_raw_l)
rm(seur_raw_list)
gc()


merged_samples$cluster = df_umap$cluster

saveRDS(merged_samples, paste0('Objects/Integrated_',INPUT_NAME,'_mt', 
                               MIT_CUT_OFF,'_lib', LIB_SIZE_CUT_OFF,'_harmony.rds'))
head(colnames(merged_samples))




##################### Applying correction to the merged datasets  ###################

anchors <- FindIntegrationAnchors(object.list = seur_norm_list, anchor.features = features)
cell_names = sapply(1:length(seur_norm_list), function(i) paste0(colnames(seur_norm_list[[i]]), '_',names(seur_norm_list)[i]))
# this command creates an 'integrated' data assay
merged_samples <- IntegrateData(anchorset = anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(merged_samples) <- "integrated"

#### checking if the order of the names haven't changed
sum(unlist(lapply(seur_norm_list, colnames)) != unlist(lapply(str_split(colnames(merged_samples), '_'), '[[', 1)))
colnames(merged_samples) = cell_names

Syn_dict = c('1'='XC116_POD3', '2'='XC121_POD7', '3'='XC125_POD7', '4'='XC130_POD3')

# Run the standard workflow for visualization and clustering
merged_samples <- ScaleData(merged_samples, verbose = FALSE)

#saveRDS(merged_samples, paste0('Objects/Integrated_',INPUT_NAME,'_mt' ,MIT_CUT_OFF,'_lib', LIB_SIZE_CUT_OFF,'.rds'))


merged_samples <- RunPCA(merged_samples, npcs = 30, verbose = FALSE)
plot(100 * merged_samples@reductions$pca@stdev^2 / merged_samples@reductions$pca@misc$total.variance,
     pch=20,xlab="Principal Component",ylab="% variance explained",log="y")

merged_samples <- RunUMAP(merged_samples, reduction = "pca", dims = 1:30)
merged_samples <- FindNeighbors(merged_samples, reduction = "pca", dims = 1:30)
merged_samples <- FindClusters(merged_samples, resolution = 0.1)

merged_samples$mito_perc = seur_raw$mito_perc
TITLE_umap = INPUT_NAME
df_umap <- data.frame(UMAP_1=getEmb(merged_samples, 'umap')[,1], 
                      UMAP_2=getEmb(merged_samples, 'umap')[,2], 
                      library_size= merged_samples$nCount_RNA, 
                      mito_perc=seur_raw$mito_perc, 
                      n_expressed=merged_samples$nFeature_RNA,
                      POD=merged_samples$POD,
                      ID= merged_samples$ID,
                      #cluster=merged_samples$seurat_clusters,
                      sample_name=merged_samples$sample_name)


df_umap$Alb = GetAssayData(merged_samples, assay='RNA')['Alb',]
df_umap$Ptprc = GetAssayData(merged_samples)['Ptprc',]
df_umap$Cd68 = GetAssayData(merged_samples, assay='RNA')['Cd68',]
df_umap$Marco = GetAssayData(merged_samples)['Marco',]
df_umap$Lyz2 = GetAssayData(merged_samples)['Lyz2',]
df_umap$Cd74 = GetAssayData(merged_samples)['Cd74',]
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=Ptprc))+geom_point(alpha=0.6)+theme_bw()+
  scale_color_viridis(direction = -1)+ggtitle(TITLE_umap)+labs(caption = INPUT_NAME)

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=library_size))+geom_point()+theme_bw()+
  scale_color_viridis(direction = -1)+ggtitle(TITLE_umap)+labs(caption = INPUT_NAME)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=mito_perc))+geom_point()+theme_bw()+
  scale_color_viridis(direction = -1)+ggtitle(TITLE_umap)+labs(caption = INPUT_NAME)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=n_expressed))+geom_point()+theme_bw()+
  scale_color_viridis(direction = -1)+ggtitle(TITLE_umap)+labs(caption = INPUT_NAME)

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=cluster))+geom_point()+theme_bw()+
  ggtitle(TITLE_umap)+labs(caption = INPUT_NAME)

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=POD))+geom_point(alpha=0.6)+theme_bw()+
  ggtitle(TITLE_umap)+labs(caption = INPUT_NAME)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=ID))+geom_point(alpha=0.6)+theme_bw()+
  ggtitle(TITLE_umap)+labs(caption = INPUT_NAME)



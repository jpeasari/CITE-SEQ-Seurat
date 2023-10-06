#### Importing R Packages

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(writexl)
library(readxl)
library(tidyverse)
library(gridExtra)
library(data.table)


### Setting Working Directory

dir <- getwd()
setwd(dir)
dir

sample <- "output/Project"  
sample_name <- "N13M3"

sample <- as.character(sample)


## Creating a new directory
dir_1 <- paste(sample,"Output",sep="_")
dir_2 <- paste(sample,"Plots",sep="_")
dir_1

if (!dir.exists(dir_1)){
  dir.create(dir_1)
}else{
  print("dir exists")
}

if (!dir.exists(dir_2)){
  dir.create(dir_2)
}else{
  print("dir exists")
}







### Reading 10X genomics input data
### "s2880 = NR1 " "s3018 = NR2" "s3120 = R1" "s3184 = R2" and so on...........

dirs_group01 <- list.dirs(path = 'data/Group01/', recursive = F, full.names = F)
dirs_group01

dirs_group02 <- list.dirs(path = 'data/Group02/', recursive = F, full.names = F)
dirs_group02

flag_x <- 0
flag_y <- 0

name_x <- "R"
name_y <- "NR"

list_x <- list() # To store seurat objects
list_y <- list()

list_all_x <- list() # TO store sample names
list_all_y <- list()


######################

#### Function to create Seurat object for all the samples
makeObject <- function(dataRead,sampleName){
  
  al1 <- CreateSeuratObject(counts = dataRead$`Gene Expression`,project = sampleName)
  al1[["ADT"]] <- CreateAssayObject(counts = dataRead$`Antibody Capture`)
  al1 <- PercentageFeatureSet(al1, "^MT-", col.name = "percent_mito")
  al1 <- PercentageFeatureSet(al1, "^RP[SL]", col.name = "percent_ribo")
  al1 <- PercentageFeatureSet(al1, "^HB[^(P)]", col.name = "percent_hb")
  al1 <- PercentageFeatureSet(al1, "PECAM1|PF4", col.name = "percent_plat")
  selected_c <- WhichCells(al1, expression = nFeature_RNA > 200 & nCount_RNA < 10000 & percent_mito < 20 & percent_ribo > 5)
  al1 <- subset(al1, cells = selected_c)
  
}


##########



for (x in dirs_group01) {
  list_all_x[[flag_x + 1]] <- paste0(name_x,"_",as.character(flag_x+1))
  namex <- paste0(name_x,as.character(flag_x+1))
  namex <- Read10X(data.dir = paste0("data/Group01/", x))
  namexx <- paste0(name_x,"_",as.character(flag_x+1))
  namexx <- makeObject(namex,paste0(name_x,"_",as.character(flag_x+1)))
  list_x[[flag_x + 1]] <- namexx
  flag_x <- flag_x + 1
  
}
flag_x <- 0



for (y in dirs_group02) {
  list_all_y[[flag_y + 1]] <- paste0(name_y,"_",as.character(flag_y+1))
  namey <- paste0(name_y,as.character(flag_y+1))
  namey <- Read10X(data.dir = paste0("data/Group02/", y))
  nameyy <- paste0(name_y,"_",as.character(flag_y+1))
  nameyy <- makeObject(namey,paste0(name_y,"_",as.character(flag_y+1)))
  list_y[[flag_y + 1]] <- nameyy
  flag_y <- flag_y + 1
  
  
}
flag_y <- 0




all_samples_names <- unlist(c(list_all_x,list_all_y))
all_sample_objects <- unlist(c(list_x,list_y))
rest_samples <- all_sample_objects[2:length(all_sample_objects)]

group01_objects <- unlist(list_all_x)
group02_objects <- unlist(list_all_y)


merged_seurat <- merge(all_sample_objects[[1]], y = rest_samples,add.cell.ids = all_samples_names,project = 'DanaFarber')
View(merged_seurat@meta.data)
# create a sample empty column in merged_seurat
merged_seurat$sample <- rownames(merged_seurat@meta.data)

# split sample column into three more columns separated by _
###                                   sample  No  Bar code
###   NR_1_AAACCTGAGATGTGTA-1 becomes   NR      1   AAACCTGAGATGTGTA-1
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('sample', 'No', 'Barcode'), 
                                    sep = '_')
#############################################################################################
# calculate mitochondrial percentage
## We don't need to run below two lines again. We already ran  ran this  step at the begining in line 35
##merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern='^MT-')
##merged_seurat_filtered <- subset(merged_seurat, subset = nCount_RNA > 800 & nFeature_RNA > 500 & mitoPercent < 20)
#############################################################################################

merged_seurat_filtered <- merged_seurat

# perform standard workflow steps to figure out if we see any batch effects --------
merged_seurat_filtered <- NormalizeData(object = merged_seurat_filtered)
merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered)
merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered)
merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered)

ElbowPlot(merged_seurat_filtered)

merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:20)
merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered)
merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:20)

tiff(paste(dir_2,'dimplot_00.tiff',sep="/"),units="in", width=7, height=6, res=300)
p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'sample', label = TRUE, repel = TRUE,raster=FALSE)
p1
dev.off()

# perform integration to correct for batch effects ------
obj.list <- SplitObject(merged_seurat_filtered, split.by = 'sample')
for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
  
}

# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 1000)

# find integration anchors (CCA)
anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features,dims = 1:30,reduction = "rpca",reference = c(1, length(obj.list)),k.anchor = 5)

# integrate data
seurat.integrated <- IntegrateData(anchorset = anchors, dims = 1:20)
seurat.integrated_copy <- seurat.integrated
seurat.integrated <- seurat.integrated_copy

DefaultAssay(seurat.integrated) <- "integrated"


# Scale data, run PCA and UMAP and visualize integrated data
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated,npcs = 50)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:20)
seurat.integrated <- FindNeighbors(seurat.integrated, reduction = "pca", dims = 1:20)
seurat.integrated <- FindClusters(seurat.integrated, resolution = 1)

seurat.integrated_1 <- seurat.integrated

DefaultAssay(seurat.integrated_1) <- "RNA"
##DefaultAssay(seurat.integrated_1) <- "ADT" antibody-derived tags (ADT)


tiff(paste(dir_2,'dimplot_01.tiff',sep="/"),units="in", width=7, height=6, res=300)
p3 <- DimPlot(seurat.integrated_1, reduction = 'umap', group.by = 'sample',raster=FALSE, label = TRUE, repel = TRUE)
p3
dev.off()


tiff(paste(dir_2,'dimplot_02.tiff',sep="/"),units="in", width=7, height=6, res=300)
p4 <- DimPlot(seurat.integrated_1, reduction = 'umap', label = TRUE, repel = TRUE,raster=FALSE)
p4
dev.off()


tiff(paste(dir_2,'dimplot_03.tiff',sep="/"),units="in", width=7, height=6, res=300)
p5 <- DimPlot(seurat.integrated_1, reduction = 'umap', split.by = 'sample',raster=FALSE, repel = TRUE)
p5
dev.off()


tiff(paste(dir_2,'dimplot_04.tiff',sep="/"),units="in", width=7, height=6, res=300)
p6 <- DimPlot(seurat.integrated_1, reduction = "umap", split.by = "orig.ident",label = TRUE,repel = TRUE)
p6
dev.off()


tiff(paste(dir_2,'dimplot_05.tiff',sep="/"),units="in", width=7, height=6, res=300)
combined <- grid.arrange(p3,p4, ncol = 2, nrow = 1)
combined
dev.off()





seurat.integrated_1[['cell_labels']] <- seurat.integrated_1@active.ident
seurat.integrated_1 <- SetIdent(seurat.integrated_1, value = 'orig.ident')

tiff(paste(dir_2,'dimplot_06.tiff',sep="/"),units="in", width=7, height=6, res=300)
p7 <- FeaturePlot(subset(seurat.integrated_1, idents=group01_objects), features = c("GZMB"))
p7
dev.off()

tiff(paste(dir_2,'dimplot_07.tiff',sep="/"),units="in", width=7, height=6, res=300)
p8 <- FeaturePlot(subset(seurat.integrated_1, idents= group02_objects), features = c("GZMB"))
p8
dev.off()

tiff(paste(dir_2,'dimplot_08.tiff',sep="/"),units="in", width=7, height=6, res=300)
FeaturePlot(seurat.integrated_1, features = c("GZMB","GZMA"), split.by = "sample")
dev.off()


tiff(paste(dir_2,'dimplot_09.tiff',sep="/"),units="in", width=10, height=6, res=300)
VlnPlot(seurat.integrated_1, features = c("GZMB","GZMA"), split.by = "sample")
dev.off()



tiff(paste(dir_2,'dimplot_10.tiff',sep="/"),units="in", width=10, height=6, res=300)
DotPlot(seurat.integrated_1, features = c("GZMB","GZMA"), split.by = "sample") + RotatedAxis()
dev.off()

tiff(paste(dir_2,'dimplot_11.tiff',sep="/"),units="in", width=10, height=6, res=300)
VlnPlot(seurat.integrated_1, features = c("GZMB","GZMA"), pt.size = 0, group.by = "sample") +
  xlab("cluster_id") +
  NoLegend()
dev.off()




### Cell Proporting in each sample group

PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    
    results = lapply(list, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}

calc_helper <- function(object,genes){
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
  }else{return(NA)}
}




################

md <- seurat.integrated_1@meta.data %>% as.data.table
cell_in_clusters <- md[, .N, by = c("orig.ident", "seurat_clusters")]
df <- as.data.frame(cell_in_clusters)
df <- df[order(df$seurat_clusters, decreasing = FALSE),]
dir_3 <- paste(dir_1,"cells_in_clusters.csv",sep = "/")
dir_4 <- paste(dir_1,"cells_prop_each_sample.csv",sep = "/")

grp_by <- md[, .N, by = c("orig.ident", "seurat_clusters")] %>% dcast(., orig.ident ~ seurat_clusters, value.var = "N")
grp_by <- as.data.frame(grp_by)
write.csv(df,dir_3)
cells <- PrctCellExpringGene(seurat.integrated_1 ,genes =c("GZMB"), group.by = "sample")
cells_by <- as.data.frame(cells)
write.csv(cells_by,dir_4)

#################


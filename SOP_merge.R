library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)

rm(list = ls())
setwd('/mnt/data/home/tycloud/20240912_宫颈癌公共单细胞数据/')
merge <- readRDS('merge.rds')

### 单独处理每个数据
step1 <- function(){
  
  t1 <- Sys.time()
  
  ### 单个创建seurat对象
  create_seu_obj <- function(data_path, data_name, grade){
    sc_data <- Read10X(data_path)
    sc_data <- CreateSeuratObject(counts = sc_data, 
                                  project = data_name, 
                                  min.cells = 3, 
                                  min.features = 200)
    sc_data$'dataset' <- data_name
    sc_data$'grade' <- grade
    sc_data[["percent.mt"]] <- PercentageFeatureSet(sc_data, pattern = "^MT-")
    sc_data <- subset(sc_data, subset = nCount_RNA > 200 & percent.mt <= 10)
    return(sc_data)
  }
  ### dataset1
  d1s1 <- create_seu_obj('Dataset1/Normal', 'd1s1', 'NC')
  d1s2 <- create_seu_obj('Dataset1/Tumor', 'd1s2', 'CESC')
  ### dataset2
  d2s1 <- create_seu_obj('Dataset2/sample1', 'd2s1', 'CESC')
  d2s2 <- create_seu_obj('Dataset2/sample2', 'd2s2', 'CESC')
  d2s3 <- create_seu_obj('Dataset2/sample3', 'd2s3', 'CESC')
  d2s4 <- create_seu_obj('Dataset2/sample4', 'd2s4', 'NC')
  d2s5 <- create_seu_obj('Dataset2/sample5', 'd2s5', 'NC')
  d2s6 <- create_seu_obj('Dataset2/sample6', 'd2s6', 'NC')
  ### dataset3
  d3s1 <- create_seu_obj('Dataset3/sample1', 'd3s1', 'NC')
  d3s2 <- create_seu_obj('Dataset3/sample2', 'd3s2', 'NC')
  d3s3 <- create_seu_obj('Dataset3/sample3', 'd3s3', 'NC')
  d3s4 <- create_seu_obj('Dataset3/sample4', 'd3s4', 'HSIL')
  d3s5 <- create_seu_obj('Dataset3/sample5', 'd3s5', 'HSIL')
  d3s6 <- create_seu_obj('Dataset3/sample6', 'd3s6', 'CESC')
  d3s7 <- create_seu_obj('Dataset3/sample7', 'd3s7', 'CESC')
  d3s8 <- create_seu_obj('Dataset3/sample8', 'd3s8', 'CESC')
  d3s9 <- create_seu_obj('Dataset3/sample9', 'd3s9', 'CESC')
  d3s10 <- create_seu_obj('Dataset3/sample10', 'd3s10', 'MLN')
  ### dataset4
  d4s1 <- create_seu_obj('Dataset4/sample1', 'd4s1', 'NC')
  d4s2 <- create_seu_obj('Dataset4/sample2', 'd4s2', 'NC')
  d4s3 <- create_seu_obj('Dataset4/sample3', 'd4s3', 'NC')
  d4s4 <- create_seu_obj('Dataset4/sample4', 'd4s4', 'NC')
  d4s5 <- create_seu_obj('Dataset4/sample5', 'd4s5', 'HSIL')
  d4s6 <- create_seu_obj('Dataset4/sample6', 'd4s6', 'HSIL')
  d4s7 <- create_seu_obj('Dataset4/sample7', 'd4s7', 'CESC')
  d4s8 <- create_seu_obj('Dataset4/sample8', 'd4s8', 'CESC')
  d4s9 <- create_seu_obj('Dataset4/sample9', 'd4s9', 'CESC')
  ### dataset6
  d6s1 <- create_seu_obj('Dataset6/sample1', 'd6s1', 'CESC')
  d6s2 <- create_seu_obj('Dataset6/sample2', 'd6s2', 'CESC')
  d6s3 <- create_seu_obj('Dataset6/sample3', 'd6s3', 'CESC')
  d6s4 <- create_seu_obj('Dataset6/sample4', 'd6s4', 'CESC')
  d6s5 <- create_seu_obj('Dataset6/sample5', 'd6s5', 'CESC')
  d6s6 <- create_seu_obj('Dataset6/sample6', 'd6s6', 'CESC')
  d6s7 <- create_seu_obj('Dataset6/sample7', 'd6s7', 'CESC')
  ### dataset5
  create_seu_obj <- function(data_path, data_name, grade){
    mx <- read.csv(data_path, row.names = 1, check.names = F)
    sc_data <- CreateSeuratObject(counts = mx, 
                                  project = data_name, 
                                  min.cells = 3, 
                                  min.features = 200)
    sc_data$'dataset' <- data_name
    sc_data$'grade' <- grade
    sc_data[["percent.mt"]] <- PercentageFeatureSet(sc_data, pattern = "^MT-")
    sc_data <- subset(sc_data, subset = nCount_RNA > 200 & percent.mt <= 10)
    return(sc_data)
  }
  d5s1 <- create_seu_obj('Dataset5/sample1/AD1.csv', 'd5s1', 'CEAD')
  d5s2 <- create_seu_obj('Dataset5/sample2/AD2.csv', 'd5s2', 'CEAD')
  d5s3 <- create_seu_obj('Dataset5/sample3/AD3.csv', 'd5s3', 'CEAD')
  d5s4 <- create_seu_obj('Dataset5/sample4/SCC1.csv', 'd5s4', 'CESC')
  d5s5 <- create_seu_obj('Dataset5/sample5/SCC2.csv', 'd5s5', 'CESC')
  d5s6 <- create_seu_obj('Dataset5/sample6/SCC3.csv', 'd5s6', 'CESC')
  
  sc_data_list <- list(d1s1, d1s2,
                       d2s1, d2s2, d2s3, d2s4, d2s5, d2s6, 
                       d3s1, d3s2, d3s3, d3s4, d3s5, d3s6, d3s7, d3s8, d3s9, d3s10, 
                       d4s1, d4s2, d4s3, d4s4, d4s5, d4s6, d4s7, d4s8, d4s9, 
                       d5s1, d5s2, d5s3, d5s4, d5s5, d5s6, 
                       d6s1, d6s2, d6s3, d6s4, d6s5, d6s6, d6s7)
  t2 <- Sys.time()
  print(t2-t1)
  return(sc_data_list)
  
}

sc_data_list <- step1() # 24 min

### 整合
step2 <- function(){
  t1 <- Sys.time()
  sc_data_list <- lapply(X = sc_data_list,
                         FUN = SCTransform,
                         assay = "SCT",
                         return.only.var.genes = FALSE,
                         verbose = F)
  var.features <- SelectIntegrationFeatures(object.list = sc_data_list, nfeatures = 3000)
  merge <- merge(x = sc_data_list[[1]], y = sc_data_list[2:length(sc_data_list)], merge.data = T)
  VariableFeatures(merge) <- var.features
  merge <- RunPCA(merge, verbose = FALSE)
  merge <- RunHarmony(merge, assay.use = "SCT", group.by.vars = "dataset")
  merge <- RunUMAP(merge, reduction = "harmony", dims = 1:30)
  DimPlot(merge, group.by = 'dataset')
  t2 <- Sys.time()
  print(t2-t1)
  return(merge)
}

### 聚类
merge <- FindNeighbors(merge, dims = 1:30)
merge <- FindClusters(merge, resolution = 0.03)
DimPlot(merge, label = T, label.box = T)
merge <- PrepSCTFindMarkers(merge)  # 30 min
deg <- FindAllMarkers(merge, only.pos = TRUE, logfc.threshold = 1)
saveRDS(merge, file = 'merge.rds')

### 看marker在各亚群的表达
p <- DimPlot(merge, label = T, repel = T, raster = F)
ggsave('umap.png', p, width = 10, height = 8, device = 'png', bg = 'white', dpi = 300)
FeaturePlot(merge, features = 'CD3E')
deg %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC) -> top10

epithelial_pan_marker <- c('EPCAM', 'CDH1', 'CLDN3', 'CLDN7', 'OCLN', 'KRT19',
                           'KRT5', 'KRT14', 'TP63',
                           'KRT7', 'KRT8', 'KRT18', 'MUC1')
VlnPlot(merge, features = epithelial_pan_marker, pt.size = 0, flip = T, stack = T)
fibroblast_pan_marker <- c('VIM', 'COL1A1', 'PDGFRA', 'ACTA2', 'THY1', 'SPARC', 'MMP2', 'S100A4', 'TNC')
VlnPlot(merge, features = fibroblast_pan_marker, pt.size = 0, flip = T, stack = T)
immune_cells_pan_marker <- c('CD3E', 'CD4', 'CD8A', 'FOXP3', 'GZMB', 'PRF1',
                             'NCAM1', 'KLRD1', 'NKG7', 'FCGR3A',
                             'CD19', 'MS4A1', 'CD79A', 'CD79B', 'IGHM',
                             'ITGAX', 'HLA-DRA', 'CCR7', 'CLEC9A',
                             'CD68', 'CD14', 'CSF1R', 'MRC1', 'ITGAM')
VlnPlot(merge, features = immune_cells_pan_marker, pt.size = 0, flip = T, stack = T)
endothelial_cells_pan_marker <- c('PECAM1', 'VWF', 'VEGFR1', 'VEGFR2', 'ENG', 'ENG', 'CLDN5', 
                                  'TIE1', 'TIE2', 'ICAM1', 'ESAM', 'CDH5', 'PLVAP')
VlnPlot(merge, features = endothelial_cells_pan_marker, pt.size = 0, flip = T, stack = T)
smooth_muscle_cells_pan_marker <- c('ACTA2', 'MYH11', 'TAGLN', 'CNN1', 'DES', 'CALD1', 
                                    'MYLK', 'LMOD1', 'PDGFRB', 'SMTN')
VlnPlot(merge, features = smooth_muscle_cells_pan_marker, pt.size = 0, flip = T, stack = T)

### class注释
new.cluster.ids <- c('Epithelial_cells','Unknown','T_NK_cells','Fibroblasts_smooth_muscle_cells',
                     'Endothelial_cells', 'DC_macrophages', 'Neutrophils', 'Epithelial_cells',
                     'B_cells', 'Mast_cells', 'Epithelial_cells')
names(new.cluster.ids) <- levels(merge)
merge <- RenameIdents(merge, new.cluster.ids)
p <- DimPlot(merge, reduction = "umap", label = TRUE, pt.size = 0.5, repel = T,
             alpha = 0.3, raster = F) + NoLegend()
ggsave('umap.png', p, width = 10, height = 10, device = 'png', bg = 'white', dpi = 300)

### Unknown亚群单独重聚类注释
unknown <- subset(merge, idents = "Unknown")
unknown <- FindClusters(unknown, resolution = 0.03)
DimPlot(unknown, reduction = "umap", label = TRUE, pt.size = 0.5, 
        repel = T, raster = T, alpha = 1)
VlnPlot(unknown, features = epithelial_pan_marker, pt.size = 0, flip = T, stack = T) # 0，4
VlnPlot(unknown, features = fibroblast_pan_marker, pt.size = 0, flip = T, stack = T) # 2
VlnPlot(unknown, features = immune_cells_pan_marker, pt.size = 0, flip = T, stack = T) # 3 
VlnPlot(unknown, features = endothelial_cells_pan_marker, pt.size = 0, flip = T, stack = T)
VlnPlot(unknown, features = smooth_muscle_cells_pan_marker, pt.size = 0, flip = T, stack = T) # 2
new.cluster.ids_unknwon <- c('Epithelial_cells','T_NK_cells','Fibroblasts_smooth_muscle_cells',
                             'Epithelial_cells','Epithelial_cells')
names(new.cluster.ids_unknwon) <- levels(unknown)
unknown <- RenameIdents(unknown, new.cluster.ids_unknwon)
DimPlot(unknown, reduction = "umap", label = TRUE, pt.size = 0.5, 
        repel = T, raster = T, alpha = 1) + NoLegend()

### 替换Unknown亚群
merge$'class' <- Idents(merge)
epi_id <- Cells(subset(unknown, idents = "Epithelial_cells"))
fibro_id <- Cells(subset(unknown, idents = "Fibroblasts_smooth_muscle_cells"))
t_nk_id <- Cells(subset(unknown, idents = "T_NK_cells"))
merge@meta.data[epi_id,'class'] <- 'Epithelial_cells'
merge@meta.data[fibro_id,'class'] <- 'Fibroblasts_smooth_muscle_cells'
merge@meta.data[t_nk_id,'class'] <- 'T_NK_cells'
merge$class <- factor(merge$class, levels = c('B_cells','DC_macrophages','Epithelial_cells',
                                              'Endothelial_cells', 'Fibroblasts_smooth_muscle_cells', 
                                              'Mast_cells', 'Neutrophils', 'T_NK_cells'))
Idents(merge) <- merge$class
p <- DimPlot(merge, reduction = "umap", label = FALSE, pt.size = 0.5, repel = T,
             alpha = 0.3, raster = F) + NoLegend()
ggsave('umap.png', p, width = 10, height = 10, device = 'png', bg = 'white', dpi = 300)

### pan marker umap
FeaturePlot(merge, features = c('KRT19', 'KRT5', 'KRT8', 'KRT18', 'EPCAM'), raster = T, ncol = 2)
FeaturePlot(merge, features = c('COL1A1', 'ACTA2', 'THY1', 'MYH11', 'TAGLN','CALD1'), raster = T, ncol = 2)
FeaturePlot(merge, features = c('PECAM1', 'VWF', 'ENG', 'CLDN5', 'PLVAP'), raster = T, ncol = 2)
FeaturePlot(merge, features = c('CD3E', 'CD8A', 'GZMB', 'NKG7'), raster = T, ncol = 2)
FeaturePlot(merge, features = c('CD79A', 'MZB1', 'IGHG1'), raster = T, ncol = 2)
FeaturePlot(merge, features = c('CD68', 'CD14', 'CSF1R', 'MRC1'), raster = T, ncol = 2)
FeaturePlot(merge, features = c('S100A8', 'S100A9', 'SOD2'), raster = T, ncol = 2)
FeaturePlot(merge, features = c('KIT', 'TPSAB1'), raster = T, ncol = 2)
genes <- c('KRT19', 'KRT5', 'KRT8', 'KRT18', 'EPCAM',
           'COL1A1', 'ACTA2', 'THY1', 'MYH11', 'TAGLN','CALD1',
           'PECAM1', 'VWF', 'ENG', 'CLDN5', 'PLVAP',
           'CD3E', 'CD8A', 'GZMB', 'NKG7',
           'CD79A', 'MZB1', 'IGHG1',
           'CD68', 'CD14', 'CSF1R', 'MRC1',
           'S100A8', 'S100A9', 'SOD2',
           'KIT', 'TPSAB1')
for (gene in genes){
  p <- FeaturePlot(merge, features = gene, raster = T, cols = c('floralwhite', 'firebrick1'), alpha = 1)
  ggsave(paste0(gene, '.png'), p, width = 6, height = 6, device = 'png', bg = 'white', dpi = 300)
}

### 把八个class的细胞单独保存为rds
epi <- subset(merge, idents = "Epithelial_cells")
saveRDS(epi, file = 'epi.rds')

fibro_smooth <- subset(merge, idents = "Fibroblasts_smooth_muscle_cells")
saveRDS(fibro_smooth, file = 'fibro_smooth.rds')

endo <- subset(merge, idents = "Endothelial_cells")
saveRDS(endo, file = 'endo.rds')

t_nk <- subset(merge, idents = "T_NK_cells")
saveRDS(t_nk, file = 't_nk.rds')

b_cell <- subset(merge, idents = "B_cells")
saveRDS(b_cell, file = 'b_cell.rds')

mast <- subset(merge, idents = "Mast_cells")
saveRDS(mast, file = 'mast.rds')

neutro <- subset(merge, idents = "Neutrophils")
saveRDS(neutro, file = 'neutro.rds')

dc_macro <- subset(merge, idents = "DC_macrophages")
saveRDS(dc_macro, file = 'dc_macro.rds')







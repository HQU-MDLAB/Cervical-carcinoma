library(Seurat)
library(ggplot2)
library(dplyr)
library(dittoSeq)
library(monocle)
library(distdimscr)
library(tidyr)
library(ggalluvial)
library(pheatmap)
library(lsa) # cosine similarity
library(scales)
library(corrplot)
library(RColorBrewer)

rm(list = setdiff(ls(), c('epi', 'nc_hsil')))
rm(list = ls())
setwd('/mnt/data/home/tycloud/20250528_CJY_HPV_MS/')
merged <- readRDS('/mnt/data/home/tycloud/20240912_宫颈癌公共单细胞数据/merge.rds')
epi <- readRDS('/mnt/data/home/tycloud/20240912_宫颈癌公共单细胞数据/epi.rds')
nc_hsil <- readRDS('/mnt/data/home/tycloud/20240912_宫颈癌公共单细胞数据/上皮细胞/nc_hsil.rds')


# Figure 1 ----------------------------------------------------------------

### 初步注释
new.cluster.ids <- c('Epithelial_cells','Unknown','T_NK_cells','Fibroblasts_smooth_muscle_cells',
                     'Endothelial_cells', 'DC_macrophages', 'Neutrophils', 'Epithelial_cells',
                     'B_cells', 'Mast_cells', 'Epithelial_cells')
names(new.cluster.ids) <- levels(merged)
merged <- RenameIdents(merged, new.cluster.ids)

### Unknown亚群单独重聚类注释
unknown <- subset(merged, idents = "Unknown")
unknown <- FindClusters(unknown, resolution = 0.03)
new.cluster.ids_unknwon <- c('Epithelial_cells','T_NK_cells','Fibroblasts_smooth_muscle_cells',
                             'Epithelial_cells','Epithelial_cells')
names(new.cluster.ids_unknwon) <- levels(unknown)
unknown <- RenameIdents(unknown, new.cluster.ids_unknwon)

### 替换Unknown亚群
merged$'class' <- Idents(merged)
epi_id <- Cells(subset(unknown, idents = "Epithelial_cells"))
fibro_id <- Cells(subset(unknown, idents = "Fibroblasts_smooth_muscle_cells"))
t_nk_id <- Cells(subset(unknown, idents = "T_NK_cells"))
merged@meta.data[epi_id,'class'] <- 'Epithelial_cells'
merged@meta.data[fibro_id,'class'] <- 'Fibroblasts_smooth_muscle_cells'
merged@meta.data[t_nk_id,'class'] <- 'T_NK_cells'
merged$class <- factor(merged$class, levels = c('B_cells','DC_macrophages','Epithelial_cells',
                                                'Endothelial_cells', 'Fibroblasts_smooth_muscle_cells', 
                                                'Mast_cells', 'Neutrophils', 'T_NK_cells'))
Idents(merged) <- merged$class

### Figure. 1B
p <- DimPlot(merged, reduction = "umap", label = FALSE, pt.size = 0.5, repel = T,
             alpha = 0.3, raster = F) + NoLegend() +
  theme(
    axis.line = element_blank(),         
    axis.text = element_blank(),         
    axis.ticks = element_blank(),        
    axis.title = element_blank()         
  )
ggsave('umap.pdf', p, width = 8, height = 8, device = 'pdf', bg = 'white', dpi = 300)
ggsave('umap.tif', p, width = 8, height = 8, device = 'tif', bg = 'white', dpi = 300)
ggsave('umap.jpeg', p, width = 8, height = 8, device = 'jpeg', bg = 'white', dpi = 300)

### Figure. 1C
genes <- c('KRT19', 'KRT5', 'KRT8', 'KRT18', 'EPCAM',
           'COL1A1', 'ACTA2', 'THY1', 'MYH11', 'TAGLN','CALD1',
           'PECAM1', 'VWF', 'ENG', 'CLDN5', 'PLVAP',
           'CD3E', 'CD8A', 'GZMB', 'NKG7',
           'CD79A', 'MZB1', 'IGHG1',
           'CD68', 'CD14', 'CSF1R', 'MRC1',
           'S100A8', 'S100A9', 'SOD2',
           'KIT', 'TPSAB1')
for (gene in genes){
  p <- FeaturePlot(merged, features = gene, raster = F,
                   cols = c('gray20', 'red'), alpha = 1) + NoLegend() +
    labs(title = '') +
    theme(
      axis.line = element_blank(),         
      axis.text = element_blank(),         
      axis.ticks = element_blank(),        
      axis.title = element_blank(),
    )
  ggsave(paste0(gene, '.jpeg'), p, width = 8, height = 8, device = 'jpeg', bg = 'white', dpi = 300)
}

# Figure 2 -----------------------------------------------------------------

### Figure. 2A
epi <- FindNeighbors(epi, dims = 1:30)
epi <- FindClusters(epi, resolution = 0.1)
epi <- RunUMAP(epi, dims = 1:30)
p <- DimPlot(epi, label = F, repel = T, raster = F, pt.size = 1,alpha = 0.3,
             cols = c('#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', 
                      '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5', 
                      '#8c564b', 'gold')) + NoLegend() +
  theme(
    axis.line = element_blank(),         
    axis.text = element_blank(),         
    axis.ticks = element_blank(),        
    axis.title = element_blank(),
  )
ggsave('umap_epi_unlabeled.jpeg', p, width = 8, height = 8, device = 'jpeg', 
       bg = 'white', dpi = 300)

### Figure. 2B
epi$grade <- factor(epi$grade, levels = c('NC', 'HSIL', 'CESC', 'CEAD', 'MLN'))
p <- dittoBarPlot(epi, main = '', x.labels.rotate = F, var = "ident", 
                  group.by = "grade",retain.factor.levels = T,
                  color.panel = c('#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', 
                                  '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5', 
                                  '#8c564b', 'gold'))+
  theme(
    axis.text.x = element_text(size = 0),
    axis.text.y = element_text(size = 0),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  ) + NoLegend() +
  xlab('')+
  ylab('')
ggsave('cluster_bar.pdf', p, width = 8, height = 6, device = 'pdf', bg = 'white', dpi = 300)

### Figure. 2C
p <- dittoBarPlot(epi, main = '', x.labels.rotate = T, var = "ident", 
                  group.by = "dataset",retain.factor.levels = T,
                  color.panel = c('#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', 
                                  '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5', 
                                  '#8c564b', 'gold'))+
  theme(
    axis.text.x = element_text(size = 0),
    axis.text.y = element_text(size = 0),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  ) + NoLegend() +
  xlab('')+
  ylab('')
ggsave('cluster_dataset_bar.pdf', p, width = 8, height = 6, device = 'pdf', bg = 'white', dpi = 300)

### Figure. 2D
deg <- FindAllMarkers(epi, only.pos = TRUE, logfc.threshold = 1, recorrect_umi = FALSE)
write.csv(deg, 'DEG_epi.csv', row.names = F)
deg %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
levels(Idents(epi)) <- paste0('epi', levels(Idents(epi)))
mean_exp <- AggregateExpression(epi, assays = 'SCT', features = top5$gene, 
                                group.by = 'ident', return.seurat = T)
p <- dittoHeatmap(mean_exp,
                  assay = "SCT", 
                  cluster_cols = F, cluster_rows = F,
                  scale = "row", show_colnames = F, show_rownames = F,
                  angle = 0, gaps_row = seq(5,60,5), border_color = 'gray20')
ggsave('deg_heatmap.pdf', p, width = 8.5, height = 13, device = 'pdf', bg = 'white', dpi = 300)

### Figure. 2E

### 拟时分析
# nc_hsil <- subset(epi, grade %in% c('NC', 'HSIL'))
# nc_hsil <- subset(nc_hsil, idents = c('epi1', 'epi3', 'epi4', 'epi8', 'epi9', 'epi11'))
# nc_hsil_cds <- as.CellDataSet(nc_hsil, assay = 'SCT')
# nc_hsil_cds <- estimateSizeFactors(nc_hsil_cds)
# nc_hsil_cds <- estimateDispersions(nc_hsil_cds)
# nc_hsil_cds <- detectGenes(nc_hsil_cds, min_expr = 0.1)
# expressed_genes <- row.names(subset(fData(nc_hsil_cds),
#                                     num_cells_expressed >= 10))
# diff_test_res <- differentialGeneTest(nc_hsil_cds[expressed_genes,],
#                                       fullModelFormulaStr = "~grade") # 1.5 h
# ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
# nc_hsil_cds <- setOrderingFilter(nc_hsil_cds, ordering_genes)
# plot_ordering_genes(nc_hsil_cds)
# nc_hsil_cds <- reduceDimension(nc_hsil_cds, max_components = 2, method = 'DDRTree') # 6 h
# nc_hsil_cds <- orderCells(nc_hsil_cds)

nc_hsil_cds <- readRDS('/mnt/data/home/tycloud/20240912_宫颈癌公共单细胞数据/上皮细胞/nc_hsil_cds.rds')
p <- plot_cell_trajectory(nc_hsil_cds, color_by = "State", cell_size = 2,
                          show_branch_points = F) +
  theme_void() + NoLegend()
ggsave('monocle_state.jpeg', p, width = 6, height = 6, device = 'jpeg', bg = 'white', dpi = 300)

### Figure. 2F

# table(pData(nc_hsil_cds)$State, pData(nc_hsil_cds)$seurat_clusters)
# table(pData(nc_hsil_cds)$State, pData(nc_hsil_cds)$grade)
# nc_hsil_cds <- orderCells(nc_hsil_cds, root_state = 3)

p <- plot_cell_trajectory(nc_hsil_cds, color_by = "Pseudotime", cell_size = 2,
                          show_branch_points = F)+
  scale_color_viridis_c(option = "magma") +
  theme_void() + NoLegend()
ggsave('monocle_time.jpeg', p, width = 6, height = 6, device = 'jpeg', bg = 'white', dpi = 300)

# saveRDS(nc_hsil_cds, 'nc_hsil_cds.rds')

### Figure. 2G
library(scales)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(DO.db)
R.utils::setOption("clusterProfiler.download.method",'auto') #解决clusterProfiler网络连接问题

### state差异分析
Idents(nc_hsil) <- 'state'
colnames(nc_hsil@meta.data)
deg <- FindAllMarkers(nc_hsil, only.pos = TRUE, logfc.threshold = 1, recorrect_umi = FALSE)
deg %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
write.csv(deg, file = 'DEG_nc_hsil_state.csv', row.names = F)

### clusterProfiler
clusterProfiler_GOsum <- function(state){
  
  df <- read.csv('DEG_nc_hsil_state.csv')
  df <- df %>% filter(cluster == state)
  df <- df$gene
  df <- bitr(df,
             fromType = "SYMBOL",
             toType = "ENTREZID",
             OrgDb = org.Hs.eg.db)
  colnames(df) <- c('gene','entrez')
  
  ### MF
  GO_MF <- enrichGO(gene = df$gene,
                    keyType = "SYMBOL",
                    OrgDb = org.Hs.eg.db,
                    ont = "MF",
                    pAdjustMethod = "fdr",
                    pvalueCutoff = 0.1,
                    qvalueCutoff = 0.1,
                    readable = FALSE)
  result1 <- GO_MF@result %>% slice_min(n = 5,order_by = pvalue)
  
  ### BP
  GO_BP <- enrichGO(gene = df$gene,
                    keyType = "SYMBOL",
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    pAdjustMethod = "fdr",
                    pvalueCutoff = 0.1,
                    qvalueCutoff = 0.1,
                    readable = FALSE)
  result2 <- GO_BP@result %>% slice_min(n = 5,order_by = pvalue)
  
  ### CC
  GO_CC <- enrichGO(gene = df$gene,
                    keyType = "SYMBOL",
                    OrgDb = org.Hs.eg.db,
                    ont = "CC",
                    pAdjustMethod = "fdr",
                    pvalueCutoff = 0.1,
                    qvalueCutoff = 0.1,
                    readable = FALSE)
  result3 <- GO_CC@result %>% slice_min(n = 5,order_by = pvalue)
  
  ### KEGG
  kegg <- enrichKEGG(gene = df$entrez,
                     organism = 'hsa',
                     pvalueCutoff = 0.1)
  result4 <- kegg@result %>% slice_min(n = 5,order_by = pvalue)
  result4 <- result4[,-c(1:2)]
  
  ### Finalization
  result <- rbind(result1, result2, result3,result4)
  result$'classification' <- c(rep("MF",nrow(result1)), rep("BP",nrow(result2)), 
                               rep("CC",nrow(result3)), rep("KEGG",nrow(result4)))
  
  ### 将GeneRatio转换为小数
  result$GeneRatio <- as.numeric(str_extract(result$GeneRatio, "^[^/]+")) / as.numeric(str_extract(result$GeneRatio, "(?<=/).*"))
  write.csv(result, file = paste0('function_state', as.character(state), '.csv'), row.names = F)
  
  return(result)
}
result1 <- clusterProfiler_GOsum(state = 1)
result2 <- clusterProfiler_GOsum(state = 2)
result3 <- clusterProfiler_GOsum(state = 3)
result4 <- clusterProfiler_GOsum(state = 4)
result5 <- clusterProfiler_GOsum(state = 5)

result1$'state' <- '1'
result2$'state' <- '2'
result3$'state' <- '3'
result4$'state' <- '4'
result5$'state' <- '5'

result <- do.call(rbind, list(result1, result2, result3, result4, result5))

featured_pathway <- c(
  
  ### state 3
  'epidermal cell differentiation',
  'structural constituent of skin epidermis',
  'ribosomal subunit',
  ### state 4
  'regulation of innate immune response',
  'NOD-like receptor signaling pathway',
  'autophagosome',
  'macroautophagy',
  'ubiquitin-protein transferase activity',
  'ubiquitin-like protein ligase binding',
  'viral life cycle',
  ### state 2
  'Viral protein interaction with cytokine and cytokine receptor',
  'chemokine activity',
  'CXCR chemokine receptor binding',
  'antimicrobial humoral response',
  'IL-17 signaling pathway',
  ### state 5
  'Antigen processing and presentation',
  'MHC protein complex assembly',
  'MHC protein complex binding',
  'amide binding',
  ### state 1
  'ATP-dependent chromatin remodeling',
  'Polycomb repressive complex',
  'double-strand break repair',
  'histone modifying activity',
  'catalytic activity, acting on RNA',
  'nuclear speck'
)
result <- filter(result, result$Description %in% featured_pathway)

p <- ggplot(result, aes(x = factor(Description, levels = featured_pathway), y = -log10(p.adjust))) +
  geom_point(aes(color = state), size = 10) +  
  geom_bar(aes(fill = state), stat = "identity", width = 0.1)+
  geom_text(aes(label = Count, size = 8))+
  coord_flip()+
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15)
  )+
  ylab('-log10 adjusted P')+
  xlab('')+
  scale_x_discrete(labels = label_wrap(40))+
  scale_y_continuous(expand = c(0.08, 0))+NoLegend()
ggsave('function_state.pdf', p, width = 9, height = 8, device = 'pdf', bg = 'white', dpi = 300)

# Figure 3 ----------------------------------------------------------------

### Figure. 3A
featured_genes <- c(
  'GOLGA4', 'DDX17', 'SON', # 1
  'CXCL2', 'CXCL3', 'CXCL8', 'RGCC', 'CSF3', 'RARRES1', 'HP', # 2
  'KRT14', 'KRT6B', 'LGALS7B', 'DAPL1', 'AKR1B10', 'CRYAB', 'LGALS7', 'VSNL1', 'SBSN', # 3
  'MUC21', 'PRSS27', 'ECM1', 'TMPRSS11B', 'KRT78', 'B3GALT5-AS1', # 4
  'FCGRT', 'IFITM2' # 5
)
p <- VlnPlot(nc_hsil, features = featured_genes, group.by = 'state', 
             flip = T, stack = T, combine = T, fill.by = 'ident') + 
  NoLegend()
ggsave('state_marker_vln.pdf', p, width = 6, height = 8, device = 'pdf', bg = 'white', dpi = 300)

### Figure. 3B
visium_score <- function(data_path, sample_name){
  
  WK1 <- Load10X_Spatial(data.dir = data_path,
                         filename = "filtered_feature_bc_matrix.h5", 
                         assay = 'Spatial')
  WK1 <- SCTransform(WK1, return.only.var.genes = F, assay = 'Spatial')
  
  p <- SpatialFeaturePlot(WK1, features = c('KRT5'), crop = T, image.alpha = 0, pt.size.factor = 2)
  ggsave(paste0(sample_name, '_visium_KRT5.pdf'), p, width = 6, height = 6, device = 'pdf', bg = 'white', dpi = 300)
  
  WK1 <- AddModuleScore(
    object = WK1,
    features = list(c('GOLGA4', 'DDX17', 'SON')),
    ctrl = 100,
    name = 'state1_'
  )
  WK1 <- AddModuleScore(
    object = WK1,
    features = list(c('CXCL2', 'CXCL3', 'CXCL8', 'RGCC', 'CSF3', 'RARRES1', 'HP')),
    ctrl = 100,
    name = 'state2_'
  )
  WK1 <- AddModuleScore(
    object = WK1,
    features = list(c('KRT14', 'KRT6B', 'LGALS7B', 'DAPL1', 'AKR1B10', 'CRYAB', 'LGALS7', 'VSNL1', 'SBSN')),
    ctrl = 100,
    name = 'state3_'
  )
  WK1 <- AddModuleScore(
    object = WK1,
    features = list(c('MUC21', 'PRSS27', 'ECM1', 'TMPRSS11B', 'KRT78', 'B3GALT5-AS1')),
    ctrl = 100,
    name = 'state4_'
  )
  WK1 <- AddModuleScore(
    object = WK1,
    features = list(c('FCGRT', 'IFITM2')),
    ctrl = 100,
    name = 'state5_'
  )
  
  p <- SpatialFeaturePlot(WK1, features = c('state1_1'), crop = T, image.alpha = 0, pt.size.factor = 2)
  ggsave(paste0(sample_name, '_visium_state1.pdf'), p, width = 6, height = 6, device = 'pdf', bg = 'white', dpi = 300)
  
  p <- SpatialFeaturePlot(WK1, features = c('state2_1'), crop = T, image.alpha = 0, pt.size.factor = 2)
  ggsave(paste0(sample_name, '_visium_state2.pdf'), p, width = 6, height = 6, device = 'pdf', bg = 'white', dpi = 300)
  
  p <- SpatialFeaturePlot(WK1, features = c('state3_1'), crop = T, image.alpha = 0, pt.size.factor = 2)
  ggsave(paste0(sample_name, '_visium_state3.pdf'), p, width = 6, height = 6, device = 'pdf', bg = 'white', dpi = 300)
  
  p <- SpatialFeaturePlot(WK1, features = c('state4_1'), crop = T, image.alpha = 0, pt.size.factor = 2)
  ggsave(paste0(sample_name, '_visium_state4.pdf'), p, width = 6, height = 6, device = 'pdf', bg = 'white', dpi = 300)
  
  p <- SpatialFeaturePlot(WK1, features = c('state5_1'), crop = T, image.alpha = 0, pt.size.factor = 2)
  ggsave(paste0(sample_name, '_visium_state5.pdf'), p, width = 6, height = 6, device = 'pdf', bg = 'white', dpi = 300)
  
}
visium_score(data_path = '/mnt/data/home/tycloud/20240912_宫颈癌公共单细胞数据/空转数据/WK1/outs/', sample_name = 'WK1')
visium_score(data_path = '/mnt/data/home/tycloud/20240912_宫颈癌公共单细胞数据/空转数据/WK2/outs/', sample_name = 'WK2')
visium_score(data_path = '/mnt/data/home/tycloud/20240912_宫颈癌公共单细胞数据/空转数据/WK3/outs/', sample_name = 'WK3')
visium_score(data_path = '/mnt/data/home/tycloud/20240912_宫颈癌公共单细胞数据/空转数据/WK4/outs/', sample_name = 'WK4')

### Figure. 3C

### 两个类别之间的相关性热图
cluster_exp <- AggregateExpression(nc_hsil, assays = 'SCT', 
                                   group.by = 'seurat_clusters', 
                                   return.seurat = F)
cluster_exp <- t(cluster_exp$SCT)

state_exp <- AggregateExpression(nc_hsil, assays = 'SCT', 
                                 group.by = 'state', 
                                 return.seurat = F)
state_exp <- t(state_exp$SCT)

similarity_matrix <- matrix(0, nrow = nrow(cluster_exp), ncol = nrow(state_exp))
for (i in 1:nrow(cluster_exp)) {
  for (j in 1:nrow(state_exp)) {
    similarity_matrix[i, j] <- cosine(cluster_exp[i, ], state_exp[j, ])
  }
}
rownames(similarity_matrix) <- paste0('epi', c(1,3,4,8,9,11))
colnames(similarity_matrix) <- paste0('state', c(1,2,3,4,5))

p <- pheatmap(similarity_matrix,display_numbers = T,border_color = 'gray55',
              show_rownames = F, show_colnames = F, 
              angle_col = 0, fontsize_number = 25,legend = F)
ggsave('matrix_cosine_similarity.pdf', p, width = 7, height = 6, device = 'pdf', bg = 'white', dpi = 300)

### Figure. 3E

### deg火山图
epi3_8 <- subset(nc_hsil, seurat_clusters %in% c(3,8)) # 9238 cells
Idents(epi3_8) <- 'SCT_snn_res.0.1'
deg_38 <- FindAllMarkers(epi3_8, only.pos = TRUE, logfc.threshold = 1, recorrect_umi = FALSE)
deg_38[which(deg_38$cluster == 8),'avg_log2FC'] <- -deg_38[which(deg_38$cluster == 8),'avg_log2FC']

### 挑选兴趣基因
featured_panel <- c(
  ### epi3
  ### 跨膜蛋白酶，可能参与病毒感染的防御
  'TMPRSS11E', 'TMPRSS2', 'TMPRSS11B',
  ### 氧化应激，激活免疫反应
  'DUOX2', 'DUOXA2',
  ### 抗菌肽，参与固有免疫
  'RNASE7',
  ### 正向调控自噬体的形成
  'TP53INP2',
  ### 
  'CDKN2A',
  
  ### epi8
  ### chemokine
  'CCL2', 'CCL14', 'CCL26', 
  ### HSP heat shock proteins
  'HSPA1A', 'HSPA1B', 'HSPA6', 'DNAJB1', 'CRYAB'
)

### 火山图
p <- EnhancedVolcano(deg_38,
                     lab = rownames(deg_38),
                     x = 'avg_log2FC',
                     y = 'p_val',
                     FCcutoff = 2,
                     pCutoff = 10e-50,
                     pointSize = 5,
                     labSize = 7,
                     subtitle = 'epi8 (left) epi3 (right)',
                     selectLab = featured_panel,
                     boxedLabels = TRUE,
                     drawConnectors = T)
ggsave('volcano_epi3_epi8.pdf', p, width = 10, height = 10, device = 'pdf', bg = 'white', dpi = 300)

### Figure. 3F

### 全部亚群的top100 deg (只有这样才能做出有差异的生存分析)
Idents(nc_hsil) <- 'seurat_clusters'
deg <- FindAllMarkers(nc_hsil, only.pos = TRUE, logfc.threshold = 1, recorrect_umi = FALSE)
deg %>% 
  group_by(cluster) %>% 
  top_n(100, wt = avg_log2FC) -> deg
write.csv(deg, 'deg_nc_hsil_cluster.csv')

### 生存分析
library(TCGAbiolinks)
library(SummarizedExperiment)

### 下载数据
### RNA表达数据
query <- GDCquery(
  project = "TCGA-CESC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query)
data <- GDCprepare(query) # 2 h
saveRDS(data, 'TCGA_CESC_SE.rds')
### 更改基因名
data <- readRDS('/mnt/data/home/tycloud/20240912_宫颈癌公共单细胞数据/上皮细胞/TCGA_CESC_SE.rds')
expr_matrix <- assay(data)
gene_meta <- rowData(data)$gene_name
names(gene_meta) <- rownames(rowData(data))
rownames(expr_matrix) <- gene_meta[rownames(expr_matrix)]

### 根据基因/甲基化表达获得分组信息
### epi3和epi8的生存分析
epi3_marker <- filter(deg, cluster == 'epi3') %>% dplyr::select(gene) %>% unlist()
epi3_marker <- base::intersect(epi3_marker, rownames(expr_matrix))
epi3_expression <- expr_matrix[epi3_marker, ]
epi3_expression <- apply(epi3_expression, 2, mean)

epi8_marker <- filter(deg, cluster == 'epi8') %>% dplyr::select(gene) %>% unlist()
epi8_marker <- base::intersect(epi8_marker, rownames(expr_matrix))
epi8_expression <- expr_matrix[epi8_marker, ]
epi8_expression <- apply(epi8_expression, 2, mean)

patient_meta <- colData(data)

library(survival)
library(survminer)

surv_data <- data.frame(
  barcode = patient_meta$barcode,
  time = ifelse(is.na(patient_meta$days_to_death), 
                patient_meta$days_to_last_follow_up / 30, 
                patient_meta$days_to_death / 30),
  status = ifelse(patient_meta$vital_status == "Dead", 1, 0),
  age = patient_meta$age_at_diagnosis,
  stage = patient_meta$figo_stage,
  gender = patient_meta$gender,
  epi3 = epi3_expression, 
  epi8 = epi8_expression
)
surv_data <- dplyr::filter(surv_data, time <= 150)

cutpoint <- surv_cutpoint(
  surv_data,
  time = "time",
  event = "status",
  variables = c('epi3', 'epi8'),
  minprop = 0.1,
  progressbar = TRUE
)

surv_data$'epi3_group' <- ifelse(surv_data$epi3 > cutpoint$epi3$estimate, "High", "Low")
surv_data$'epi8_group' <- ifelse(surv_data$epi8 > cutpoint$epi8$estimate, "High", "Low")

surv_object <- Surv(surv_data$time, surv_data$status)
fit <- survfit(surv_object ~ epi8_group, data = surv_data)
p <- ggsurvplot(
  fit, 
  data = surv_data,
  conf.int = TRUE,          # 显示置信区间
  conf.int.style = "ribbon", # 使用阴影带（默认）
  conf.int.alpha = 0.2,     # 设置透明度（0-1）
  pval = TRUE,              # 添加p值
  risk.table = FALSE,        # 添加风险表
  palette = c("#E41A1C", "#377EB8"), # 自定义颜色
  xlab = "Time (Months)", 
  ylab = "Survival Probability"
)
ggsave('survival_epi8.pdf', p$plot, width = 6, height = 6, device = 'pdf', bg = 'white', dpi = 300)

surv_object <- Surv(surv_data$time, surv_data$status)
fit <- survfit(surv_object ~ epi3_group, data = surv_data)
p <- ggsurvplot(
  fit, 
  data = surv_data,
  conf.int = TRUE,          # 显示置信区间
  conf.int.style = "ribbon", # 使用阴影带（默认）
  conf.int.alpha = 0.2,     # 设置透明度（0-1）
  pval = TRUE,              # 添加p值
  risk.table = FALSE,        # 添加风险表
  palette = c("#E41A1C", "#377EB8"), # 自定义颜色
  xlab = "Time (Months)", 
  ylab = "Survival Probability"
)
ggsave('survival_epi3.pdf', p$plot, width = 6, height = 6, device = 'pdf', bg = 'white', dpi = 300)

### Figure. 3G

### inferCNV
library(infercnv)

### 矩阵文件
### reference
get_ref_mx <- function(ref_rds){
  reference <- readRDS(ref_rds)
  reference <- subset(reference, grade == 'NC')
  reference_counts <- GetAssayData(reference, layer = 'counts')
  print(dim(reference_counts))
  return(reference_counts)
}
b_cells <- get_ref_mx('/mnt/data/home/tycloud/20240912_宫颈癌公共单细胞数据/b_cell.rds')
mast_cells <- get_ref_mx('/mnt/data/home/tycloud/20240912_宫颈癌公共单细胞数据/mast.rds')
neutrophils <- get_ref_mx('/mnt/data/home/tycloud/20240912_宫颈癌公共单细胞数据/neutro.rds')
t_nk <- get_ref_mx('/mnt/data/home/tycloud/20240912_宫颈癌公共单细胞数据/t_nk.rds')
### epi3_8
epi3_8_counts <- GetAssayData(epi3_8, layer = 'counts')
### 矩阵拼接
counts <- do.call(cbind, list(b_cells, mast_cells, neutrophils, t_nk,
                              epi3_8_counts))
dim(counts)

### 细胞类型注释文件
epi3_8$'new' <- paste0('epi', epi3_8$seurat_clusters)
annot <- data.frame(
  'barcode' = colnames(counts),
  'celltype' = c(rep('b_cells', ncol(b_cells)), 
                 rep('mast_cells', ncol(mast_cells)),
                 rep('neutrophils', ncol(neutrophils)),
                 rep('t_nk_cells', ncol(t_nk)),
                 epi3_8$new)
)
colnames(annot) <- NULL
write.table(annot, file = '../inferCNV/annot.txt', sep = '\t', row.names = F)

### 创建inferCNV对象
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = counts,
                                     annotations_file = '../inferCNV/annot.txt',
                                     delim = "\t",
                                     gene_order_file = '../inferCNV/hg38_gencode_v27.txt',
                                     min_max_counts_per_cell = c(100, +Inf),
                                     ref_group_names = c("b_cells", 'mast_cells', 
                                                         'neutrophils', 't_nk_cells'))

### 运行inferCNV
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff = 0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir = "../inferCNV/上皮细胞1", 
                             cluster_by_groups = T,
                             denoise = T,
                             HMM = F,
                             num_threads = 60,
                             leiden_resolution = 0.00001
)

### 计算CNV score
cnvScore <- function(data){
  data <- data %>% 
    as.matrix() %>%
    t() %>% 
    scale() %>% 
    rescale(to = c(-1, 1)) %>% 
    t()
  
  cnv_score <- as.data.frame(colSums(data * data))
  return(cnv_score)
}

infercnv_obj <-  readRDS('/mnt/data/home/tycloud/20240912_宫颈癌公共单细胞数据/inferCNV/上皮细胞1/run.final.infercnv_obj')
cnv_score <- cnvScore(infercnv_obj@expr.data)
cnv_score[colnames(counts[,infercnv_obj@reference_grouped_cell_indices$b_cells]),'group'] <- 'b_cells'
cnv_score[colnames(counts[,infercnv_obj@reference_grouped_cell_indices$mast_cells]),'group'] <- 'mast_cells'
cnv_score[colnames(counts[,infercnv_obj@reference_grouped_cell_indices$neutrophils]),'group'] <- 'neutrophils'
cnv_score[colnames(counts[,infercnv_obj@reference_grouped_cell_indices$t_nk_cells]),'group'] <- 't_nk_cells'
cnv_score[colnames(counts[,infercnv_obj@observation_grouped_cell_indices$epi3]),'group'] <- 'epi3'
cnv_score[colnames(counts[,infercnv_obj@observation_grouped_cell_indices$epi8]),'group'] <- 'epi8'
colnames(cnv_score) <- c('value', 'group')
cnv_score$group <- factor(cnv_score$group, levels = c("epi3", "epi8", "b_cells", 
                                                      "mast_cells", "neutrophils", "t_nk_cells"))

cnv_score %>% 
  group_by(group) %>% 
  summarise(new = median(log1p(value)))

p <- ggplot(cnv_score,aes(x = group, y = log1p(value)))+
  geom_boxplot(aes(fill = group))+
  geom_signif(
    comparisons = list(c("epi3", "epi8"), 
                       c("epi3", "b_cells"),
                       c("epi3", "mast_cells"),
                       c("epi3", "neutrophils"),
                       c("epi3", "t_nk_cells")),
    map_signif_level = TRUE, textsize = 10, step_increase = 0.15
  )+
  geom_hline(yintercept = 3.24, linetype="solid", color = 'red3', linewidth = 2)+
  ylim(NA, 7)+
  theme_classic()+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )+
  coord_flip()+
  ylab('')+
  xlab('')+NoLegend()+
  scale_fill_manual(values = c('red3', 'red3', 'gray', 'gray', 'gray', 'gray'))
ggsave('infercnv_CNVscore.pdf', p, width = 6, height = 6, device = 'pdf', bg = 'white', dpi = 300)

# Figure 4 ----------------------------------------------------------------

### Figure. 4A

### hdWGCNA

library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
theme_set(theme_cowplot())
set.seed(2024)
enableWGCNAThreads(nThreads = 60)

rm(list = setdiff(ls(), c('nc_hsil')))
Idents(nc_hsil) <- factor(paste0('epi', Idents(nc_hsil)),
                          levels = c('epi1', 'epi3', 'epi4', 'epi8', 'epi9', 'epi11'))
nc_hsil$new <- Idents(nc_hsil)

### 新建seurat对象，否则ModuleEigengenes报错
seurat_obj <- CreateSeuratObject(nc_hsil@assays$SCT@counts)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)
seurat_obj$'new' <- nc_hsil$new
seurat_obj$'dataset' <- nc_hsil$dataset
seurat_obj@reductions$'umap' <- nc_hsil@reductions$umap
Idents(seurat_obj) <- 'new'

### seurat object setup
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "cc" # the name of the hdWGCNA experiment
)

### construct metacells in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = "new", # specify the columns in seurat_obj@meta.data to group by
  reduction = 'umap', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'new' # set the Idents of the metacell seurat object
)
seurat_obj <- NormalizeMetacells(seurat_obj)

### expression matrix setup
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = c('epi1', 'epi3', 'epi4', 'epi8', 'epi9', 'epi11'), # the name of the group of interest in the group.by column
  group.by = 'new', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using SCT assay
  slot = 'data' # using normalized data
)

### select soft-power threshold
# Test different soft powers
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)
# plot the results
plot_list <- PlotSoftPowers(seurat_obj)
# assemble with patchwork
wrap_plots(plot_list, ncol = 2)

### construct co-expression network
t1 <- Sys.time()
seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name = 'new', # name of the topoligical overlap matrix written to disk
)
t2 <- Sys.time()
t2-t1 # 7 min

png('hdWGCNA_dendrogram.png', width = 8, height = 8, res = 300,units = 'in')
PlotDendrogram(seurat_obj, main = 'hdWGCNA Dendrogram')
dev.off()

### inspect the topoligcal overlap matrix (TOM)
TOM <- GetTOM(seurat_obj)
table(colnames(TOM) %in% VariableFeatures(seurat_obj))

### compute harmonized module eigengenes
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars = "dataset"
)
# harmonized module eigengenes
hMEs <- GetMEs(seurat_obj)
# module eigengenes
MEs <- GetMEs(seurat_obj, harmonized = FALSE)

### compute module connectivity
# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'new', group_name = c('epi1', 'epi3', 'epi4', 'epi8', 'epi9', 'epi11')
)
# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "Module"
)
# plot genes ranked by kME for each module
pdf('hdWGCNA_kME.pdf', width = 12, height = 8)
PlotKMEs(seurat_obj, ncol = 4, text_size = 4.5)
dev.off()

### getting the module assignment table
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')
head(modules[,1:6])
# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)
head(hub_df)

### ModuleRadarPlot
pdf('hdWGCNA_radar.pdf', width = 8, height = 8)
ModuleRadarPlot(
  seurat_obj,
  group.by = 'new',
  axis.label.size=4,
  grid.label.size=4,
  ncol = 4
)
dev.off()

### compute hub gene signature scores
library(UCell)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 10,
  method = 'Seurat'
)

### GO
clusterProfiler_GOsum <- function(gene_list){
  
  df <- gene_list
  df <- bitr(df,
             fromType = "SYMBOL",
             toType = "ENTREZID",
             OrgDb = org.Hs.eg.db)
  colnames(df) <- c('gene','entrez')
  
  ### MF
  GO_MF <- enrichGO(gene = df$gene,
                    keyType = "SYMBOL",
                    OrgDb = org.Hs.eg.db,
                    ont = "MF",
                    pAdjustMethod = "fdr",
                    pvalueCutoff = 0.1,
                    qvalueCutoff = 0.1,
                    readable = FALSE)
  result1 <- GO_MF@result %>% slice_min(n = 5,order_by = pvalue)
  
  ### BP
  GO_BP <- enrichGO(gene = df$gene,
                    keyType = "SYMBOL",
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    pAdjustMethod = "fdr",
                    pvalueCutoff = 0.1,
                    qvalueCutoff = 0.1,
                    readable = FALSE)
  result2 <- GO_BP@result %>% slice_min(n = 5,order_by = pvalue)
  
  ### CC
  GO_CC <- enrichGO(gene = df$gene,
                    keyType = "SYMBOL",
                    OrgDb = org.Hs.eg.db,
                    ont = "CC",
                    pAdjustMethod = "fdr",
                    pvalueCutoff = 0.1,
                    qvalueCutoff = 0.1,
                    readable = FALSE)
  result3 <- GO_CC@result %>% slice_min(n = 5,order_by = pvalue)
  
  ### KEGG
  kegg <- enrichKEGG(gene = df$entrez,
                     organism = 'hsa',
                     pvalueCutoff = 0.1)
  result4 <- kegg@result %>% slice_min(n = 5,order_by = pvalue)
  result4 <- result4[,-c(1:2)]
  
  ### Finalization
  result <- rbind(result1, result2, result3,result4)
  result$'classification' <- c(rep("MF",nrow(result1)), rep("BP",nrow(result2)), 
                               rep("CC",nrow(result3)), rep("KEGG",nrow(result4)))
  
  ### 将GeneRatio转换为小数
  result$GeneRatio <- as.numeric(str_extract(result$GeneRatio, "^[^/]+")) / as.numeric(str_extract(result$GeneRatio, "(?<=/).*"))
  # write.csv(result, file = paste0('function_state', as.character(state), '.csv'), row.names = F)
  
  return(result)
}

module1_hub <- filter(hub_df, module == 'Module1') %>% dplyr::select(gene_name) %>% unlist()
module3_hub <- filter(hub_df, module == 'Module3') %>% dplyr::select(gene_name) %>% unlist()
module5_hub <- filter(hub_df, module == 'Module5') %>% dplyr::select(gene_name) %>% unlist()

module1_hub <- clusterProfiler_GOsum(module1_hub)
module3_hub <- clusterProfiler_GOsum(module3_hub)
module5_hub <- clusterProfiler_GOsum(module5_hub)

p <- ggplot(module1_hub, aes(x = factor(Description, levels = Description), 
                             y = -log10(p.adjust))) +
  geom_point(aes(color = classification), size = 10) +  
  geom_bar(aes(fill = classification), stat = "identity", width = 0.1)+
  geom_text(aes(label = Count, size = 8))+
  coord_flip()+
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15)
  )+
  ylab('-log10 adjusted P')+
  xlab('')+
  scale_x_discrete(labels = label_wrap(40))+
  scale_y_continuous(expand = c(0.08, 0))+NoLegend()
ggsave('hdWGCNA_function_module1.pdf', p, width = 12, height = 8, 
       device = 'pdf', bg = 'white', dpi = 300)

### Figure. 4B

### 找TCGA里哪些基因的表达和FAM19A4和MIR124-2的甲基化呈正相关
data <- readRDS('/mnt/data/home/tycloud/20240912_宫颈癌公共单细胞数据/上皮细胞/TCGA_CESC_SE.rds')
expr_matrix <- assay(data)
gene_meta <- rowData(data)$gene_name
names(gene_meta) <- rownames(rowData(data))
rownames(expr_matrix) <- gene_meta[rownames(expr_matrix)]

data <- readRDS('/mnt/data/home/tycloud/20240912_宫颈癌公共单细胞数据/上皮细胞/TCGA_CESC_Methylation.rds')
met_matrix <- assay(data)
fam19a4 <- read.table('/mnt/data/home/tycloud/20240912_宫颈癌公共单细胞数据/上皮细胞/GDCdata/FAM19A4.tsv', sep = '\t', header = T,row.names = 1)
mir1242 <- read.table('/mnt/data/home/tycloud/20240912_宫颈癌公共单细胞数据/上皮细胞/GDCdata/MIR124-2.tsv', sep = '\t', header = T,row.names = 1)
met_matrix <- met_matrix[c(colnames(fam19a4), colnames(mir1242)),]
df <- data.frame(
  row.names = colnames(met_matrix),
  'FAM19A4' = apply(met_matrix[colnames(fam19a4),],2,mean),
  'MIR124_2' = apply(met_matrix[colnames(mir1242),],2,mean)
)
df[is.na(df)] <- 0
df <- t(df)

colnames(expr_matrix) <- substr(colnames(expr_matrix), 1, 12)
colnames(df) <- substr(colnames(df), 1, 12)
sample <- intersect(colnames(expr_matrix),colnames(df))
df <- df[,sample]

marker <- c(
  ### epi3
  'CDKN2A', 
  'TMPRSS11E', 'TMPRSS2', 'TMPRSS11B',
  'DUOX2', 'DUOXA2',
  'RNASE7',
  'TP53INP2',
  'CDC42EP5', 'TMSB4X', 'TMSB10', 'SAT1', 'FTH1',
  'ADGRF1', 'PSCA', 'PRSS27', 'MUC21', 'SCEL', 'ERO1A', 'ECM1', 'TTC9',
  
  ### epi8
  'CCL2', 'CCL14','CCL26',
  'HSPA1A', 'HSPA1B', 'HSPA6', 'DNAJB1', 'CRYAB',
  'RPLP1', 'RPL13', 'RPS8', 'RPL32', 'RPS3', 'RPS14', 'RPS24', 'RPS15', 'RPS7', 'RPS15A',
  'H1-2', 'H1-5', 'H2BC18', 'H1-4', 'HMGN2',
  'H4C3', 'H3C4', 'SLC25A5', 'H4C1', 'DEK',
  'H3C7', 'H4C4', 'H2AC11', 'TUBB', 'HMGB2',
  'H2BC12', 'H2AC14', 'H2AZ1', 'H2BC9', 'H2AC7'
)

df2 <- expr_matrix[marker, sample]
normalize_expr_matrix <- function(expr_matrix) {
  min_val <- min(expr_matrix, na.rm = TRUE)
  max_val <- max(expr_matrix, na.rm = TRUE)
  normalized_matrix <- (expr_matrix - min_val) / (max_val - min_val)
  return(normalized_matrix)
}
df2 <- normalize_expr_matrix(df2)
df2 <- rbind(df2, df)

df3 <- cor(t(df2))
testRes <-  cor.mtest(t(df2), conf.level = 0.95)

### 不带星号的相关性热图
pdf('corrplot_epi3_8_final.pdf', width = 12, height = 12)
COL2 <- colorRampPalette(rev(brewer.pal(100, "RdBu")))
corrplot(df3, order = 'hclust', 
         method = "color", 
         tl.col = "black", 
         addrect = 5,
         col = COL2(100),
         tl.cex = 1,
         cl.cex = 1)
dev.off()

### Figure. 4C
markers <- c('RPLP1', 'RPL13', 
             'NME2', 'GAS5', 'SNHG5', 'SNHG8', 'SNHG29', 
             'HSPA1B', 'CRYAB')
Idents(epi3_8) <- factor(paste0('epi', Idents(epi3_8)),
                         levels = c('epi8', 'epi3'))
p <- VlnPlot(epi3_8, markers, stack = T,fill.by = 'ident', flip = T,
             cols = c('#9467bd', '#ffbb78'))+ NoLegend()
ggsave('epi8_markers.pdf', p, width = 6, height = 10, device = 'pdf',
       bg = 'white', dpi = 300)


# Figure 6 ----------------------------------------------------------------

### Figure. 6
data <- readRDS('/mnt/data/home/tycloud/20240912_宫颈癌公共单细胞数据/上皮细胞/TCGA_CESC_SE.rds')
expr_matrix <- assay(data)
gene_meta <- rowData(data)$gene_name
names(gene_meta) <- rownames(rowData(data))
rownames(expr_matrix) <- gene_meta[rownames(expr_matrix)]

### 进展组和逆转组
expr_progression <- apply(expr_matrix[c('CRYAB', 'GAS5'), ], 2, mean)
expr_regression <- apply(expr_matrix[c('RPLP1', 'SNHG5'), ], 2, mean)

patient_meta <- colData(data)
# colnames(patient_meta)

surv_data <- data.frame(
  barcode = patient_meta$barcode,
  time = ifelse(is.na(patient_meta$days_to_death), 
                patient_meta$days_to_last_follow_up / 30, 
                patient_meta$days_to_death / 30),
  status = ifelse(patient_meta$vital_status == "Dead", 1, 0),
  age = patient_meta$age_at_diagnosis,
  stage = patient_meta$figo_stage,
  gender = patient_meta$gender,
  progression = expr_progression,
  regression = expr_regression
)
surv_data <- dplyr::filter(surv_data, time <= 150) 

### 阈值确定
cutpoint <- surv_cutpoint(
  surv_data,
  time = "time",
  event = "status",
  variables = c('progression', 'regression'),
  minprop = 0.1,
  progressbar = TRUE
)

surv_data$'progression_group' <- ifelse(surv_data$progression > cutpoint$progression$estimate, "High", "Low")
surv_data$'regression_group' <- ifelse(surv_data$regression > cutpoint$regression$estimate, "High", "Low")

### 生存分析
surv_object <- Surv(surv_data$time, surv_data$status)
fit <- survfit(surv_object ~ progression_group, data = surv_data)
p <- ggsurvplot(
  fit, 
  data = surv_data,
  conf.int = TRUE,          # 显示置信区间
  conf.int.style = "ribbon", # 使用阴影带（默认）
  conf.int.alpha = 0.2,     # 设置透明度（0-1）
  pval = TRUE,              # 添加p值
  risk.table = FALSE,        # 添加风险表
  palette = c("#E41A1C", "#377EB8"), # 自定义颜色
  xlab = "Time (Months)", 
  ylab = "Survival Probability"
)
ggsave('survival_progression.pdf', p$plot, width = 6, height = 6, device = 'pdf', bg = 'white', dpi = 300)

fit <- survfit(surv_object ~ regression_group, data = surv_data)
p <- ggsurvplot(
  fit, 
  data = surv_data,
  conf.int = TRUE,          # 显示置信区间
  conf.int.style = "ribbon", # 使用阴影带（默认）
  conf.int.alpha = 0.2,     # 设置透明度（0-1）
  pval = TRUE,              # 添加p值
  risk.table = FALSE,        # 添加风险表
  palette = c("#E41A1C", "#377EB8"), # 自定义颜色
  xlab = "Time (Months)", 
  ylab = "Survival Probability"
)
ggsave('survival_regression.pdf', p$plot, width = 6, height = 6, device = 'pdf', bg = 'white', dpi = 300)

### 拟合Cox模型
surv_data$progression_group <- factor(surv_data$progression_group, levels = c("Low", "High"))
surv_data$regression_group <- factor(surv_data$regression_group, levels = c("Low", "High"))

cox_model <- coxph(
  Surv(time, status) ~ progression_group + regression_group,
  data = surv_data
)
cox_summary <- summary(cox_model)

### 可视化
cox_summary$coefficients[, c("exp(coef)", "Pr(>|z|)")] # 输出HR和p值
p <- ggforest(cox_model, data = surv_data, main = "Hazard Ratios from Cox Model",
              fontsize = 0.5)
ggsave('HR.pdf', p, width = 8, height = 8, device = 'pdf', bg = 'white', dpi = 300)



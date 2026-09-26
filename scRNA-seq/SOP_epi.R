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

rm(list = setdiff(ls(), c('nc_hsil')))
rm(list = ls())
# nc_hsil <- readRDS('nc_hsil.rds')
setwd('E:/20250527_HPV文章/单细胞数据/')
epi <- readRDS('epi.rds')

# 聚类 ----------------------------------------------------------------------

### 聚类
epi <- FindNeighbors(epi, dims = 1:30)
epi <- FindClusters(epi, resolution = 0.1)
epi <- RunUMAP(epi, dims = 1:30)
p <- DimPlot(epi, label = F, repel = T, raster = F, pt.size = 1,alpha = 0.5,
             cols = c('#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', 
                      '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5', 
                      '#8c564b', 'gold')) + NoLegend()
ggsave('umap_epi_unlabeled.png', p, width = 8, height = 8, device = 'png', bg = 'white', dpi = 300)
deg <- FindAllMarkers(epi, only.pos = TRUE, logfc.threshold = 1, recorrect_umi = FALSE)
VlnPlot(epi, features = epithelial_pan_marker, pt.size = 0, flip = T, stack = T)
deg %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5

### deg热图
levels(Idents(epi)) <- paste0('epi', levels(Idents(epi)))
mean_exp <- AggregateExpression(epi, assays = 'SCT', features = top5$gene, 
                                group.by = 'ident', return.seurat = T)
p <- dittoHeatmap(mean_exp,
                  assay = "SCT", 
                  cluster_cols = F, cluster_rows = F,
                  scale = "row", show_colnames = T, angle = 0, fontsize = 10,
                  gaps_row = seq(5,60,5), border_color = 'gray67')
ggsave('deg_heatmap.png', p, width = 8, height = 10, device = 'png', bg = 'white', dpi = 300)

### 亚群数量分布图
epi$grade <- factor(epi$grade, levels = c('NC', 'HSIL', 'CESC', 'CEAD', 'MLN'))
p <- dittoBarPlot(epi, main = '', x.labels.rotate = F, var = "ident", 
                  group.by = "grade",retain.factor.levels = T,
                  color.panel = c('#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', 
                                  '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5', 
                                  '#8c564b', 'gold'))+
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15)
  )
ggsave('cluster_bar.png', p, width = 8, height = 6, device = 'png', bg = 'white', dpi = 300)
p <- dittoBarPlot(epi, main = '', x.labels.rotate = T, var = "ident", 
                  group.by = "dataset",retain.factor.levels = T,
                  color.panel = c('#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', 
                                  '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5', 
                                  '#8c564b', 'gold'))+
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15)
  )
ggsave('cluster_dataset_bar.png', p, width = 8, height = 6, device = 'png', bg = 'white', dpi = 300)

# 研究NC到HSIL的中间态 -----------------------------------------------------------

nc_hsil <- subset(epi, grade %in% c('NC', 'HSIL'))

### 各亚群细胞数量变化趋势图
new <- as.data.frame(table(nc_hsil$grade, Idents(nc_hsil)))
colnames(new) <- c('grade', 'cluster', 'value')
new <- new[new$value>0,]
new <- new %>% 
  group_by(grade) %>%                              
  mutate(total_value = sum(value)) %>%              
  group_by(grade, cluster) %>%                      
  summarise(cluster_value = sum(value),             
            percent = (cluster_value / total_value) * 100) %>% # 计算百分比
  ungroup()

p <- ggplot(new, aes(x = grade, y = percent, fill = cluster,
                     stratum = cluster, alluvium = cluster)) + 
  geom_col(width = 0.5, color = 'gray22')+
  geom_flow(width = 0.5, alpha = 0.4, knot.pos = 0)+
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15)
  )+
  labs(x = 'Grade',y = 'Percent')+
  scale_fill_manual(values = c('#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', 
                               '#98df8a', '#ff9896', '#9467bd', '#c5b0d5', 
                               '#8c564b', 'gold'))
ggsave('nc_hsil_sankey.png', p, width = 6, height = 6, device = 'png', bg = 'white', dpi = 300)

### Bhatt距离 (HSIL vs. NC)
nc_hsil <- subset(nc_hsil, idents = c('epi1', 'epi2', 'epi3', 'epi4'))
cellType <- levels(Idents(nc_hsil))
bhatt.dist <- bhatt.dist.rand <- as.data.frame(matrix(NA, ncol = length(cellType), nrow = 100))
names(bhatt.dist) <- cellType
names(bhatt.dist.rand) <- cellType

for (CT in cellType){
  
  for (j in 1:10){
    
    set.seed(j)
    
    # n = length(which(nc_hsil$grade == "NC" & Idents(nc_hsil) == CT))
    n = 100
    cells.HSIL <- sample(colnames(nc_hsil)[which(nc_hsil$grade == "HSIL" & Idents(nc_hsil) == CT)],n)
    # cells.NC <- colnames(nc_hsil)[Idents(nc_hsil) == CT & nc_hsil$grade == "NC"]
    cells.NC <- sample(colnames(nc_hsil)[which(nc_hsil$grade == "NC" & Idents(nc_hsil) == CT)],n)
    
    tmp <- nc_hsil@reductions$harmony@cell.embeddings
    
    cells.NC.pca <- tmp[cells.NC,]
    cells.HSIL.pca <- tmp[cells.HSIL,]
    
    for (i in 1:10) {
      
      d <- (j-1)*10 + i
      
      bhatt.dist[d,CT] <- dim_dist(embed_mat_x = cells.NC.pca, 
                                   embed_mat_y = cells.HSIL.pca,
                                   dims_use = 1:30,
                                   num_cells_sample = 50,
                                   distance_metric = "bhatt_dist",
                                   random_sample = FALSE)
      
      bhatt.dist.rand[d,CT] <- dim_dist(embed_mat_x = cells.NC.pca, 
                                        embed_mat_y = cells.HSIL.pca,
                                        dims_use = 1:30,
                                        num_cells_sample = 50,
                                        distance_metric = "bhatt_dist",
                                        random_sample = TRUE)
      
    }
  }
}

### 可视化
bhatt.dist <- pivot_longer(bhatt.dist, cols = 1:4)
p <- ggplot(bhatt.dist)+
  geom_violin(aes(x = name, y = value, fill = name))+
  geom_boxplot(aes(x = name, y = value, fill = name))+
  scale_fill_manual(values = c('#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c'))+
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15)
  )+
  xlab('Cluster')+
  ylab('Bhattacharyya distance')+
  NoLegend()
ggsave('bhatt_violn.png', p, width = 6, height = 6, device = 'png', bg = 'white', dpi = 300)

### 拟时分析
nc_hsil <- subset(epi, grade %in% c('NC', 'HSIL'))
nc_hsil <- subset(nc_hsil, idents = c('epi1', 'epi3', 'epi4', 'epi8', 'epi9', 'epi11'))
nc_hsil_cds <- as.CellDataSet(nc_hsil, assay = 'SCT')
nc_hsil_cds <- estimateSizeFactors(nc_hsil_cds)
nc_hsil_cds <- estimateDispersions(nc_hsil_cds)
nc_hsil_cds <- detectGenes(nc_hsil_cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(nc_hsil_cds),
                                    num_cells_expressed >= 10))
diff_test_res <- differentialGeneTest(nc_hsil_cds[expressed_genes,],
                                      fullModelFormulaStr = "~grade") # 1.5 h
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
nc_hsil_cds <- setOrderingFilter(nc_hsil_cds, ordering_genes)
plot_ordering_genes(nc_hsil_cds)
nc_hsil_cds <- reduceDimension(nc_hsil_cds, max_components = 2, method = 'DDRTree') # 6 h
nc_hsil_cds <- orderCells(nc_hsil_cds)
p <- plot_cell_trajectory(nc_hsil_cds, color_by = "State", cell_size = 0.6)+
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15)
  )+
  guides(color = guide_legend(override.aes = list(size = 6)))
ggsave('monocle_state.png', p, width = 6, height = 6, device = 'png', bg = 'white', dpi = 300)
table(pData(nc_hsil_cds)$State, pData(nc_hsil_cds)$seurat_clusters)
table(pData(nc_hsil_cds)$State, pData(nc_hsil_cds)$grade)
nc_hsil_cds <- orderCells(nc_hsil_cds, root_state = 3)
p <- plot_cell_trajectory(nc_hsil_cds, color_by = "Pseudotime", cell_size = 0.6)+
  scale_color_viridis_c(option = "magma")+
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15)
  )
ggsave('monocle_time.png', p, width = 6, height = 6, device = 'png', bg = 'white', dpi = 300)
saveRDS(nc_hsil_cds, 'nc_hsil_cds.rds')

### state桑基图
new <- as.data.frame(table(pData(nc_hsil_cds)$'State', 
                           pData(nc_hsil_cds)$'grade',
                           pData(nc_hsil_cds)$'seurat_clusters'))
colnames(new) <- c('state', 'grade', 'cluster', 'value')
new$cluster <- paste0('epi', new$cluster)
new$cluster <- factor(new$cluster, levels = c('epi1', 'epi3', 'epi4', 'epi8', 'epi9', 'epi11'))
new <- new[new$value>0,]
p <- ggplot(data = new,
            aes(axis1 = state, axis2 = grade, axis3 = cluster, y = value)) +
  xlab("Demographic") +
  geom_alluvium(aes(fill = cluster), alpha = 0.8) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 5)+
  theme_void()+
  scale_fill_manual(values = c('#aec7e8', '#ffbb78', '#2ca02c',
                               '#9467bd', '#c5b0d5','gold'))+NoLegend()
ggsave('state_sankey.png', p, width = 8, height = 6, device = 'png', bg = 'white', dpi = 300)

### 加metadata
nc_hsil <- subset(epi, grade %in% c('NC', 'HSIL'))
nc_hsil <- subset(nc_hsil, idents = c('epi1', 'epi3', 'epi4', 'epi8', 'epi9', 'epi11'))
monocle_meta <- pData(nc_hsil_cds)
table(rownames(monocle_meta) == rownames(nc_hsil@meta.data))
nc_hsil$'state' <- monocle_meta$State
nc_hsil$'pseudotime' <- monocle_meta$Pseudotime

### state差异分析
Idents(nc_hsil) <- 'state'
colnames(nc_hsil@meta.data)
deg <- FindAllMarkers(nc_hsil, only.pos = TRUE, logfc.threshold = 1, recorrect_umi = FALSE)
deg %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
write.csv(deg, file = 'deg_nc_hsil_state.csv', row.names = F)

### 功能富集
library(scales)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(DO.db)
R.utils::setOption("clusterProfiler.download.method",'auto') #解决clusterProfiler网络连接问题

### clusterProfiler
clusterProfiler_GOsum <- function(state){
  
  df <- read.csv('deg_nc_hsil_state.csv')
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
  geom_point(aes(color = state), size = 8) +  
  geom_bar(aes(fill = state), stat = "identity", width = 0.1)+
  geom_text(aes(label = Count))+
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
  scale_y_continuous(expand = c(0.03, 0))+NoLegend()
ggsave('function_state.png', p, width = 10, height = 8, device = 'png', bg = 'white', dpi = 300)

### seurat对象转为h5ad
MuDataSeurat::WriteH5AD(nc_hsil, "nc_hsil.h5ad", assay = "RNA")

### 为每个state挑选deg
deg <- read.csv('deg_nc_hsil_state.csv')
deg %>% 
  filter(cluster == 1) %>% 
  top_n(30, wt = pct.1) -> top10
featured_genes <- c(
  'GOLGA4', 'DDX17', 'SON', # 1
  'CXCL2', 'CXCL3', 'CXCL8', 'RGCC', 'CSF3', 'RARRES1', 'HP', # 2
  'KRT14', 'KRT6B', 'LGALS7B', 'DAPL1', 'AKR1B10', 'CRYAB', 'LGALS7', 'VSNL1', 'SBSN', # 3
  'MUC21', 'PRSS27', 'ECM1', 'TMPRSS11B', 'KRT78', 'B3GALT5-AS1', # 4
  'FCGRT', 'IFITM2' # 5
)
p <- VlnPlot(nc_hsil, features = featured_genes, group.by = 'state', 
             flip = T, stack = T, combine = T, fill.by = 'ident') + NoLegend()
ggsave('state_marker_vln.png', p, width = 6, height = 8, device = 'png', bg = 'white', dpi = 300)

### AddModuleScore映射到空转
visium_score <- function(data_path, sample_name){
  WK1 <- Load10X_Spatial(data.dir = data_path,
                         filename = "filtered_feature_bc_matrix.h5", 
                         assay = 'Spatial')
  WK1 <- SCTransform(WK1, return.only.var.genes = F, assay = 'Spatial')
  
  p <- SpatialFeaturePlot(WK1, features = c('KRT5'), crop = T, image.alpha = 0, pt.size.factor = 2)
  ggsave(paste0(sample_name, '_visium_KRT5.png'), p, width = 6, height = 6, device = 'png', bg = 'white', dpi = 300)
  
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
  ggsave(paste0(sample_name, '_visium_state1.png'), p, width = 6, height = 6, device = 'png', bg = 'white', dpi = 300)
  
  p <- SpatialFeaturePlot(WK1, features = c('state2_1'), crop = T, image.alpha = 0, pt.size.factor = 2)
  ggsave(paste0(sample_name, '_visium_state2.png'), p, width = 6, height = 6, device = 'png', bg = 'white', dpi = 300)
  
  p <- SpatialFeaturePlot(WK1, features = c('state3_1'), crop = T, image.alpha = 0, pt.size.factor = 2)
  ggsave(paste0(sample_name, '_visium_state3.png'), p, width = 6, height = 6, device = 'png', bg = 'white', dpi = 300)
  
  p <- SpatialFeaturePlot(WK1, features = c('state4_1'), crop = T, image.alpha = 0, pt.size.factor = 2)
  ggsave(paste0(sample_name, '_visium_state4.png'), p, width = 6, height = 6, device = 'png', bg = 'white', dpi = 300)
  
  p <- SpatialFeaturePlot(WK1, features = c('state5_1'), crop = T, image.alpha = 0, pt.size.factor = 2)
  ggsave(paste0(sample_name, '_visium_state5.png'), p, width = 6, height = 6, device = 'png', bg = 'white', dpi = 300)
  
}
visium_score(data_path = '../空转数据/WK1/outs/', sample_name = 'WK1')
visium_score(data_path = '../空转数据/WK2/outs/', sample_name = 'WK2')
visium_score(data_path = '../空转数据/WK3/outs/', sample_name = 'WK3')
visium_score(data_path = '../空转数据/WK4/outs/', sample_name = 'WK4')

# 挖掘epi3和epi8 -------------------------------------------------------------

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
data <- readRDS('TCGA_CESC_SE.rds')
expr_matrix <- assay(data)
gene_meta <- rowData(data)$gene_name
names(gene_meta) <- rownames(rowData(data))
rownames(expr_matrix) <- gene_meta[rownames(expr_matrix)]

### DNA甲基化数据
query <- GDCquery(
  project = "TCGA-CESC",
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value"
)
GDCdownload(query)
data <- GDCprepare(query) # 10 min
saveRDS(data, 'TCGA_CESC_Methylation.rds')
### 从UCSC得到FAM19A4和MIR124-2的甲基化位点序号
data <- readRDS('TCGA_CESC_Methylation.rds')
expr_matrix <- assay(data)
fam19a4 <- read.table('GDCdata/FAM19A4.tsv', sep = '\t', header = T,row.names = 1)
mir1242 <- read.table('GDCdata/MIR124-2.tsv', sep = '\t', header = T,row.names = 1)
expr_matrix <- expr_matrix[c(colnames(fam19a4), colnames(mir1242)),]

### 根据基因/甲基化表达获得分组信息
### epi3和epi8的生存分析
epi3_marker <- filter(deg, cluster == 3) %>% select(gene) %>% unlist()
epi3_marker <- base::intersect(epi3_marker, rownames(expr_matrix))
epi3_expression <- expr_matrix[epi3_marker, ]
epi3_expression <- apply(epi3_expression, 2, mean)

epi8_marker <- filter(deg, cluster == 8) %>% select(gene) %>% unlist()
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
p <- ggsurvplot(fit, data = surv_data, pval = TRUE)
ggsave('survival_epi3.png', p$plot, width = 6, height = 6, device = 'png', bg = 'white', dpi = 300)

### 挑出epi3_8对象进行差异分析

### deg火山图
epi3_8 <- subset(nc_hsil, seurat_clusters %in% c(3,8)) # 9238 cells
deg_38 <- FindAllMarkers(epi3_8, only.pos = TRUE, logfc.threshold = 1, recorrect_umi = FALSE)
deg_38[which(deg_38$cluster == 8),'avg_log2FC'] <- -deg_38[which(deg_38$cluster == 8),'avg_log2FC']

library(EnhancedVolcano)
p <- EnhancedVolcano(deg_38,
                     lab = rownames(deg_38),
                     x = 'avg_log2FC',
                     y = 'p_val',
                     FCcutoff = 2,
                     pCutoff = 10e-50,
                     pointSize = 3,
                     labSize = 4.5,
                     subtitle = 'epi8 (left) epi3 (right)')

### 找到火山图里FC_P组和deg的交集，进一步确定和预后相关的基因集
deg_38 <- p$data
epi3_intersect <- intersect(
  deg %>% filter(cluster == '3') %>% select(gene) %>% unlist,
  deg_38 %>% filter(cluster == '3', Sig == 'FC_P') %>% select(gene) %>% unlist
)
epi8_intersect <- intersect(
  deg %>% filter(cluster == '8') %>% select(gene) %>% unlist,
  deg_38 %>% filter(cluster == '8', Sig == 'FC_P') %>% select(gene) %>% unlist
)
VlnPlot(epi3_8, epi8_intersect, stack = T, flip = T)+NoLegend()

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
VlnPlot(epi3_8, featured_panel, stack = T, flip = T)+NoLegend()
DotPlot(epi3_8, featured_panel, dot.scale = 7.5, cols = c('gray', 'steelblue1'))+
  coord_flip()

### 重画火山图
p <- EnhancedVolcano(deg_38,
                     lab = rownames(deg_38),
                     x = 'avg_log2FC',
                     y = 'p_val',
                     FCcutoff = 2,
                     pCutoff = 10e-50,
                     pointSize = 3,
                     labSize = 5,
                     subtitle = 'epi8 (left) epi3 (right)',
                     selectLab = featured_panel,
                     boxedLabels = TRUE,
                     drawConnectors = T)
ggsave('volcano_epi3_epi8.png', p, width = 10, height = 10, device = 'png', bg = 'white', dpi = 300)

### GSEApy结果可视化
hallmark <- read.csv('../hallmark.csv')
c2 <- read.csv('../C2.csv')
c5 <- read.csv('../C5.csv')
c6 <- read.csv('../C6.csv')
c7 <- read.csv('../C7.csv')
c8 <- read.csv('../C8.csv')

featured_functions <- c(
  ### epi3
  'HALLMARK_INTERFERON_GAMMA_RESPONSE',
  'HALLMARK_INTERFERON_ALPHA_RESPONSE',
  'WP_TYPE_I_INTERFERON_INDUCTION_AND_SIGNALING_DURING_SARSCOV2_INFECTION',
  'BROWNE_INTERFERON_RESPONSIVE_GENES',
  'GOBP_INTERFERON_MEDIATED_SIGNALING_PATHWAY',
  'AKT_UP.V1_UP',
  'ERBB2_UP.V1_UP',
  'MYC_UP.V1_DN',
  'EGFR_UP.V1_DN',
  
  ### epi8
  'HALLMARK_TNFA_SIGNALING_VIA_NFKB',
  'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
  'FOROUTAN_INTEGRATED_TGFB_EMT_UP',
  'NAGASHIMA_EGF_SIGNALING_UP',
  'GOBP_MESENCHYMAL_CELL_PROLIFERATION',
  'GOBP_POSITIVE_REGULATION_OF_MESENCHYMAL_CELL_PROLIFERATION',
  'BCAT_GDS748_UP',
  'RELA_DN.V1_UP',
  'TGFB_UP.V1_UP',
  'LEF1_UP.V1_UP'
)
GSEA <- do.call(rbind, list(hallmark, c2, c5, c6, c7, c8))
GSEA <- GSEA %>% 
  filter(Term %in% featured_functions)
GSEA$Term <- factor(GSEA$Term, levels = featured_functions)
GSEA$group <- ifelse(GSEA$NES > 0, '1', '2')
p <- ggplot(GSEA)+
  geom_col(aes(x = Term, y = NES, fill = group))+
  coord_flip()+
  theme_classic()+
  scale_fill_manual(values = c('royalblue', 'indianred4'))+
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )+ NoLegend()
ggsave('gsea_epi3_8.png', p, width = 10, height = 5, device = 'png', bg = 'white', dpi = 300)

### SpaGene
### 读取空转数据
WK1 <- Load10X_Spatial(data.dir = '../空转数据/WK1/outs/',
                       filename = "filtered_feature_bc_matrix.h5", 
                       assay = 'Spatial')
WK1 <- SCTransform(WK1, return.only.var.genes = F, assay = 'Spatial')
wk1_niche <- read.csv('../cell2loc_WK1_niche.csv', row.names = 1)
wk1_niche_exp <- read.csv('../cell2loc_WK1.csv', row.names = 1)
meta <- merge(wk1_niche, wk1_niche_exp, by = 0)
rownames(meta) <- meta$Row.names
WK1@meta.data <- meta
Idents(WK1) <- 'region_cluster'

WK2 <- Load10X_Spatial(data.dir = '../空转数据/WK2/outs/',
                       filename = "filtered_feature_bc_matrix.h5", 
                       assay = 'Spatial')
WK2 <- SCTransform(WK2, return.only.var.genes = F, assay = 'Spatial')
wk2_niche <- read.csv('../cell2loc_WK2_niche.csv', row.names = 1)
wk2_niche_exp <- read.csv('../cell2loc_WK2.csv', row.names = 1)
meta <- merge(wk2_niche, wk2_niche_exp, by = 0)
rownames(meta) <- meta$Row.names
WK2@meta.data <- meta
Idents(WK2) <- 'region_cluster'

# new <- subset(WK1, region_cluster %in% c(0,2,12))
new <- subset(WK2, region_cluster %in% c(1,3,4))

### SpaGene
mx <- GetAssayData(new, layer = "counts", assay = "SCT")
coord <- data.frame(-new$array_row, new$array_col)

### main function
spagene_svgs <- SpaGene(mx, coord[,c(1,2)], knn = 24, normalize = T)

### sort by adjp
spagene_svgs_reorder <- spagene_svgs$spagene_res[order(spagene_svgs$spagene_res$adjp),]

### find patterns
pattern <- FindPattern(spagene_svgs, nPattern = 3)

### plot patterns
p1 <- PlotPattern(pattern, coord[, c(1:2)], pt.size = 1.5, max.cutoff = 1)

### plot heatmap
top5 <- apply(pattern$genepattern,2, function(x){names(x)[order(x,decreasing=T)][1:20]})
p2 <- pheatmap(pattern$genepattern[rownames(pattern$genepattern)%in%top5,],
               angle_col = 0,fontsize = 10)

### save results
# write.csv(spagene_svgs$spagene_res, file = 'SpaGene_result.csv')
ggsave(p1, filename = 'SpaGene_pattern.png', device = "png", height = 4,width = 15,dpi = 300, bg = 'white')
ggsave(p2, filename = 'SpaGene_heatmap.png', device = "png", height = 8,width = 8,dpi = 300, bg = 'white')

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

# ### HIST基因要重命名否则map不到功能
# top5[,'Pattern1'] <- c('H1-2', 'H1-5', 'H2BC18', 'H1-4', 'HMGN2',
#                        'H4C3', 'H3C4', 'SLC25A5', 'H4C1', 'DEK',
#                        'H3C7', 'H4C4', 'H2AC11', 'TUBB', 'HMGB2',
#                        'H2BC12', 'H2AC14', 'H2AZ1', 'H2BC9', 'H2AC7')

pattern1 <- clusterProfiler_GOsum(top5[,'Pattern1'])
pattern2 <- clusterProfiler_GOsum(top5[,'Pattern2'])
pattern3 <- clusterProfiler_GOsum(top5[,'Pattern3'])

p <- ggplot(pattern1, aes(x = factor(Description, levels = pattern1$Description), 
                          y = -log10(p.adjust))) +
  geom_point(aes(color = classification), size = 8) +  
  geom_bar(aes(fill = classification), stat = "identity", width = 0.1)+
  geom_text(aes(label = Count))+
  coord_flip()+
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15)
  )+
  ylab('-log10 adjusted P')+
  xlab('')+
  # scale_x_discrete(labels = label_wrap(50))+
  scale_y_continuous(expand = c(0.15, 0))+NoLegend()
ggsave('SpaGene_function_pattern1.png', p, width = 10, height = 8, device = 'png', bg = 'white', dpi = 300)

### 单细胞亚群和空转pattern基因集打分的相关性热图
### 先AddModuleScore
new <- AddModuleScore(
  object = new,
  features = list(top5[,'Pattern1']),
  ctrl = 100,
  name = 'Pattern1_'
)
new <- AddModuleScore(
  object = new,
  features = list(top5[,'Pattern2']),
  ctrl = 100,
  name = 'Pattern2_'
)
new <- AddModuleScore(
  object = new,
  features = list(top5[,'Pattern3']),
  ctrl = 100,
  name = 'Pattern3_'
)
### 再画corrplot
epi_pattern_cor <- function(patho_file, patho_name){
  
  df <- read.csv(patho_file)
  meta <- new@meta.data
  meta$'pathology' <- 0
  meta[df[df$pathology == 'HSIL','Barcode'],'pathology'] <- 'HSIL'
  meta[df[df$pathology == 'normal_LSIL','Barcode'],'pathology'] <- 'normal_LSIL'
  meta <- meta %>% filter(pathology != 0)
  intersect_id <- intersect(meta$Row.names, colnames(new))
  new <- new[,intersect_id]
  meta <- meta[intersect_id,]
  new@meta.data <- meta
  Idents(new) <- 'pathology'
  
  df2 <- meta
  df2 <- filter(meta, pathology == patho_name)
  df2 <- df2 %>% dplyr::select(c(q05cell_abundance_w_sf_3, q05cell_abundance_w_sf_8,
                                 Pattern1_1, Pattern2_1, Pattern3_1))
  df2 <- cor(df2)
  df2 <- df2[-c(3:5),-c(1:2)]
  rownames(df2) <- c('Abundance_epi3', 'Abundance_epi8')
  colnames(df2) <- c('Pattern1_score', 'Pattern2_score', 'Pattern3_score')
  library(corrplot)
  library(RColorBrewer)
  COL2 <- colorRampPalette(rev(brewer.pal(100, "RdBu")))
  png(paste0('SpaGene_', patho_name,'_epi_pattern_cor.png'), width = 8, 
      height = 6, res = 300,units = 'in')
  corrplot(df2, method = "circle", tl.col = "black",col = COL2(100),
           addCoef.col = "black", number.cex = 1,tl.srt = 0,tl.cex = 1.5)
  dev.off()
}
epi_pattern_cor(patho_file = 'WK1_pathology.csv', patho_name = 'HSIL')
epi_pattern_cor(patho_file = 'WK1_pathology.csv', patho_name = 'normal_LSIL')
epi_pattern_cor(patho_file = 'WK2_pathology.csv', patho_name = 'HSIL')
epi_pattern_cor(patho_file = 'WK2_pathology.csv', patho_name = 'normal_LSIL')

### CellChat
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

### 数据准备
data.input <-  GetAssayData(new, layer = "data", assay = "SCT") # normalized data matrix
meta <- data.frame(labels = paste0('Niche', Idents(new)), 
                   samples = "WK1", 
                   row.names = names(Idents(new)))
meta$samples <- factor(meta$samples)
spatial.locs <-  GetTissueCoordinates(new, scale = NULL, cols = c("imagerow", "imagecol")) 
spatial.locs <- spatial.locs[,-3]
scalefactors <-  jsonlite::fromJSON(txt = '../空转数据/WK1/outs/spatial/scalefactors_json.json')
spot.size <-  65
conversion.factor <-  spot.size/scalefactors$spot_diameter_fullres
spatial.factors <-  data.frame(ratio = conversion.factor, tol = spot.size/2)
d.spatial <- computeCellDistance(coordinates = spatial.locs, 
                                 ratio = spatial.factors$ratio, 
                                 tol = spatial.factors$tol)
min(d.spatial[d.spatial!=0]) # this value should approximately equal 100um for 10X Visium data
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs, 
                           spatial.factors = spatial.factors)

### 选择数据库
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact", 
                           key = "annotation")
cellchat@DB <- CellChatDB.use

### 数据前处理
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = F)

### 主程序
t1 <- Sys.time()
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                              distance.use = TRUE, interaction.range = 250, 
                              scale.distance = 0.01,contact.dependent = TRUE, 
                              contact.range = 100)
t2 <- Sys.time()
(t2-t1)
cellchat <- filterCommunication(cellchat, min.cells = 10)

### 结果表格
df.net <- subsetCommunication(cellchat)
df.net <- subsetCommunication(cellchat, 
                              sources.use = 'Niche2',
                              targets.use = 'Niche0')
cellchat <- computeCommunProbPathway(cellchat)
cellchat@net

### 可视化
netVisual_bubble(cellchat, remove.isolate = F, 
                 sort.by.source = T, 
                 sort.by.source.priority = T,
                 sources.use = c('Niche12', 'Niche0'),
                 targets.use = c('Niche0', 'Niche2'),
                 signaling = c('DESMOSOME', 'NOTCH'))

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
par(mfrow=c(10,10))
netAnalysis_signalingRole_network(cellchat, signaling = 'CDH', 
                                  width = 8, height = 8, font.size = 10)
netAnalysis_signalingRole_scatter(cellchat)

###
featured_lr <- c(
  
)


spatialFeaturePlot(cellchat, features = c("GJA1"), point.size = 2, 
                   color.heatmap = "Reds", direction = 1)
spatialFeaturePlot(cellchat, pairLR.use = "DSC2_DSG3", point.size = 1, 
                   do.binary = TRUE, cutoff = 0.05, enriched.only = F, 
                   color.heatmap = "Reds", direction = 1)
netAnalysis_contribution(cellchat, signaling = c('GAP', 'APP', 'DESMOSOME'))

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
b_cells <- get_ref_mx('../b_cell.rds')
mast_cells <- get_ref_mx('../mast.rds')
neutrophils <- get_ref_mx('../neutro.rds')
t_nk <- get_ref_mx('../t_nk.rds')
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

# infercnv_obj <-  readRDS('../inferCNV/上皮细胞/run.final.infercnv_obj')
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
    map_signif_level = TRUE, textsize = 4, step_increase = 0.1
  )+
  geom_hline(yintercept = 3.24, linetype="longdash", color = 'red3')+
  ylim(NA, 7)+
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15)
  )+
  coord_flip()+
  ylab('')+
  xlab('')+NoLegend()+
  scale_fill_manual(values = c('red3', 'red3', 'gray', 'gray', 'gray', 'gray'))
ggsave('infercnv_CNVscore.png', p, width = 6, height = 6, device = 'png', bg = 'white', dpi = 300)

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
png('hdWGCNA_kME.png', width = 12, height = 8, res = 300,units = 'in')
PlotKMEs(seurat_obj, ncol = 4, text_size = 4)
dev.off()

### getting the module assignment table
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')
head(modules[,1:6])
# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)
head(hub_df)

### ModuleRadarPlot
png('hdWGCNA_radar.png', width = 8, height = 8, res = 300,units = 'in')
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

p <- ggplot(module5_hub, aes(x = factor(Description, levels = Description), 
                             y = -log10(p.adjust))) +
  geom_point(aes(color = classification), size = 8) +  
  geom_bar(aes(fill = classification), stat = "identity", width = 0.1)+
  geom_text(aes(label = Count))+
  coord_flip()+
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15)
  )+
  ylab('-log10 adjusted P')+
  xlab('')+
  # scale_x_discrete(labels = label_wrap(50))+
  scale_y_continuous(expand = c(0.05, 0))+NoLegend()
ggsave('hdWGCNA_function_module5.png', p, width = 12, height = 8, device = 'png', bg = 'white', dpi = 300)

# ### 差异ME分析 (1 vs 1)
# group1 <- seurat_obj@meta.data %>% subset(new == 'epi3') %>% rownames
# group2 <- seurat_obj@meta.data %>% subset(new == 'epi8') %>% rownames
# DMEs <- FindDMEs(
#   seurat_obj,
#   barcodes1 = group1,
#   barcodes2 = group2,
#   test.use = 'wilcox',
#   wgcna_name = 'cc'
# )
# 
# head(DMEs)
# PlotDMEsLollipop(
#   seurat_obj, 
#   DMEs, 
#   wgcna_name = 'cc', 
#   pvalue = "p_val_adj"
# )
# 
# ### 差异ME分析 (1 vs all)
# DMEs_all <- FindAllDMEs(
#   seurat_obj,
#   group.by = 'new',
#   wgcna_name = 'cc'
# )
# p + PlotDMEsVolcano(
#   seurat_obj,
#   DMEs_all,
#   wgcna_name = 'cc',
#   plot_labels=FALSE,
#   show_cutoff=FALSE,mod_point_size = 4
# )
# p + facet_wrap(~group, ncol=3)

# ### 箱线图
# df <- data.frame(
#   'cell' = rownames(seurat_obj@meta.data),
#   'cluster' = seurat_obj$new,
#   'module1' = seurat_obj@misc$cc$module_scores$Module3
# )
# ggplot(df)+
#   geom_boxplot(aes(x = cluster, y = module1, fill = cluster))
# 
# ### dotplot
# # get hMEs from seurat object
# hMEs <- GetMEs(seurat_obj, harmonized = TRUE)
# modules <- GetModules(seurat_obj)
# mods <- levels(modules$module)
# mods <- mods[mods != 'grey']
# # add hMEs to Seurat meta-data
# seurat_obj@meta.data <- cbind(seurat_obj@meta.data, hMEs)
# DotPlot(seurat_obj, features = mods, group.by = 'new')
# 
# ### heatmap
# annot_row <- data.frame(
#   row.names = rownames(seurat_obj@meta.data),
#   'group' = seurat_obj$new
# )
# pheatmap(seurat_obj@misc$cc$module_scores, show_rownames = F, 
#          annotation_row = annot_row, cluster_rows = F)

### 找TCGA里哪些基因的表达和FAM19A4和MIR124-2的甲基化呈正相关
data <- readRDS('TCGA_CESC_Methylation.rds')
met_matrix <- assay(data)
fam19a4 <- read.table('GDCdata/FAM19A4.tsv', sep = '\t', header = T,row.names = 1)
mir1242 <- read.table('GDCdata/MIR124-2.tsv', sep = '\t', header = T,row.names = 1)
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

### 带星号的相关性热图
corrplot(df3,order = 'hclust',method = "circle", tl.col = "black", addrect = 6,
         p.mat = testRes$p, sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9, 
         insig = 'label_sig', pch.col = 'grey20')

### 不带星号的相关性热图
png('corrplot_epi3_8_final.png', width = 12, height = 12, res = 300,units = 'in')
COL2 <- colorRampPalette(rev(brewer.pal(100, "RdBu")))
corrplot(df3, order = 'hclust', method = "color", tl.col = "black", addrect = 5,col = COL2(100))
dev.off()

### epi3_8挑marker——生存分析——与tCIN2和pCIN2的联系
surv_function <- function(epi3_marker, epi8_marker){
  
  if(all(epi3_marker %in% rownames(expr_matrix)) & all(epi8_marker %in% rownames(expr_matrix))){
    
    epi3_marker <- base::intersect(epi3_marker, rownames(expr_matrix))
    epi3_expression <- expr_matrix[epi3_marker, ]
    epi3_expression <- apply(epi3_expression, 2, mean)
    
    epi8_marker <- base::intersect(epi8_marker, rownames(expr_matrix))
    epi8_expression <- expr_matrix[epi8_marker, ]
    epi8_expression <- apply(epi8_expression, 2, mean)
    
    patient_meta <- colData(data)
    
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
    
    # fit <- survfit(surv_object ~ epi3_group, data = surv_data)
    # p <- ggsurvplot(fit, data = surv_data, pval = TRUE)
    # ggsave('survival_epi3.png', p$plot, width = 6, height = 6, device = 'png', bg = 'white', dpi = 300)
    
    fit <- survfit(surv_object ~ epi8_group, data = surv_data)
    ggsurvplot(fit, data = surv_data, pval = TRUE)
    # ggsave('survival_epi8.png', p$plot, width = 6, height = 6, device = 'png', bg = 'white', dpi = 300)
    
  }
  
  else{print('Gene names error!')}
}
surv_function(
  ### epi3
  c('CDKN2A', 'RNASE7'), 
  
  ### epi8
  c('CDKN2A', 'MKI67')
)

# 其他 ----------------------------------------------------------------------

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
              fontsize = 15, angle_col = 0, fontsize_number = 15)
ggsave('matrix_cosine_similarity.png', p, width = 7, height = 6, device = 'png', bg = 'white', dpi = 300)







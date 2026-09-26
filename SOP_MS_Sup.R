
# Figure. S1 --------------------------------------------------------------

### Fig. S1A
p <- dittoBarPlot(merge, main = '', x.labels.rotate = F, var = "ident", 
                  group.by = "grade", retain.factor.levels = T, 
                  color.panel = hue_pal()(8))+coord_flip()+
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15)
  ) +
  xlab('Group')+
  ylab('Percent')
ggsave('Sup_cluster_bar.pdf', p, width = 8, height = 10, device = 'pdf', bg = 'white', dpi = 300)

result <- prop.table(table(merge$class, merge$grade), margin = 2) * 100
write.csv(result, 'Sup_cluster_table.csv')

### Fig. S1B
p <- dittoBarPlot(merge, main = '', x.labels.rotate = F, var = "ident", 
                  group.by = "dataset", retain.factor.levels = T, 
                  color.panel = hue_pal()(8))+ coord_flip()+
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15)
  ) +
  xlab('Dataset')+
  ylab('Percent')
ggsave('Sup_dataset_bar.pdf', p, width = 8, height = 10, device = 'pdf', bg = 'white', dpi = 300)

result <- prop.table(table(merge$class, merge$dataset), margin = 2) * 100
write.csv(result, 'Sup_dataset_table.csv')

# Figure. S2 ---------------------------------------------------------------

### Figure. S2A

genes <- c('KRT5',
           'RPLP1', 'RPS19',
           'HSPB1','DNAJB1',
           'MCM2','ALDH1A1', 'SOX2', 'CXCR4')
p <- VlnPlot(epi, genes, stack = T, flip = T,
             fill.by = 'ident', raster = F, pt.size = 0, 
             cols = c('#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', 
                      '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5', 
                      '#8c564b', 'gold')) + NoLegend() +
  theme(
    axis.title = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20, angle = 0),
    axis.title.y = element_text(size = 20)
  ) 
ggsave('Sup_epi8_genes.pdf', p, width = 8, height = 10, device = 'pdf', bg = 'white', dpi = 300)

### Figure. S2B

new <- as.data.frame(table(epi$grade, Idents(epi)))
colnames(new) <- c('grade', 'cluster', 'value')
new <- new %>% filter(grade %in% c('NC', 'HSIL'))
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
                               '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5', 
                               '#8c564b', 'gold'))
ggsave('Sup_nc_hsil_sankey.pdf', p, width = 6, height = 6, device = 'pdf', bg = 'white', dpi = 300)

### Figure. S2C

new <- as.data.frame(table(nc_hsil$state, 
                           nc_hsil$grade,
                           nc_hsil$seurat_clusters))
colnames(new) <- c('state', 'grade', 'cluster', 'value')
new$cluster <- paste0('epi', new$cluster)
new$state <- paste0('State', new$state)
new$cluster <- factor(new$cluster, levels = c('epi1', 'epi3', 'epi4', 'epi8', 'epi9', 'epi11'))
new$grade <- factor(new$grade, levels = c('NC', 'HSIL'))
new <- new[new$value>0,]

p <- ggplot(data = new,
            aes(axis1 = state, axis2 = grade, axis3 = cluster, y = value)) +
  xlab("Demographic") +
  geom_alluvium(aes(fill = cluster), alpha = 1) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 7)+
  theme_void()+
  scale_fill_manual(values = c('#aec7e8', '#ffbb78', '#2ca02c',
                               '#9467bd', '#c5b0d5', 'gold'))
ggsave('Sup_state_sankey.pdf', p, width = 12, height = 8, device = 'pdf', 
       bg = 'white', dpi = 300)










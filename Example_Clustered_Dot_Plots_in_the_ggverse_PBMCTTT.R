# https://davemcg.github.io/post/lets-plot-scrna-dotplots/

library(tidyverse)
library(ggdendro)
library(cowplot)
library(ggtree)
library(patchwork) 

#### Load data #####
  ## Example
  # GO_cluster <- read_tsv('https://github.com/davemcg/davemcg.github.io/raw/master/content/post/scRNA_dotplot_data.tsv.gz')
  # GO_cluster %>% sample_n(5)
  # GO_cluster_Ori <- GO_cluster
  GO_cluster <- GSEA_Fin.df2
  
  GO_cluster$PhenoType <-str_replace_all(GO_cluster$PhenoType,"\\+","")

#### Ori Dotplot ####
markers <- GO_cluster$GO %>% unique()

GO_cluster %>% filter(GO %in% markers) %>% 
  ggplot(aes(x=PhenoType, y = GO, color = NES, size = -log10(padj))) + 
  geom_point() 


#### Remove dots where there is zero (or near zero expression) ####
GO_cluster %>% filter(GO %in% markers)  %>% 
  filter(abs(NES) > 0, -log10(padj) > 1) %>% 
  ggplot(aes(x=PhenoType, y = GO, color = NES, size = -log10(padj))) + 
  geom_point() 


# #### Better color, better theme, rotate x axis labels ####
# GO_cluster %>% filter(GO %in% markers) %>% 
#   mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100) %>% 
#   filter(count > 0, `% Expressing` > 1) %>% 
#   ggplot(aes(x=cluster, y = GO, color = count, size = `% Expressing`)) + 
#   geom_point() + 
#   scale_color_viridis_c(name = 'log2 (count + 1)') + 
#   cowplot::theme_cowplot() + 
#   theme(axis.line  = element_blank()) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   ylab('') +
#   theme(axis.ticks = element_blank()) 
# 
# #### Tweak color scaling ####
# GO_cluster %>% filter(GO %in% markers) %>% 
#   mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100) %>% 
#   filter(count > 0, `% Expressing` > 1) %>% 
#   ggplot(aes(x=cluster, y = GO, color = count, size = `% Expressing`)) + 
#   geom_point() + 
#   cowplot::theme_cowplot() + 
#   theme(axis.line  = element_blank()) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   ylab('') +
#   theme(axis.ticks = element_blank()) +
#   scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)')

#### Hey look: ggtree ####
# make data square to calculate euclidean distance
mat <- GO_cluster %>% 
  filter(GO %in% markers) %>% 
  select(PhenoType, NES, GO) %>%  # drop unused columns to faciliate widening
  pivot_wider(names_from = PhenoType, values_from = NES) %>% 
  data.frame() # make df as tibbles -> matrix annoying
row.names(mat) <- mat$GO  # put GO in `row`
mat <- mat[,-1] #drop GO column as now in rows
mat[is.na(mat)] <- 0
clust <- hclust(dist(mat %>% as.matrix())) # hclust with distance matrix


ddgram <- as.dendrogram(clust) # create dendrogram
ggtree_plot <- ggtree::ggtree(ddgram)
ggtree_plot

#### Let’s glue them together with cowplot ####
dotplot <- GO_cluster %>% filter(GO %in% markers)  %>% 
  filter(abs(NES) > 0, -log10(padj) > 1) %>% 
  select(PhenoType, NES, GO,padj) %>%  # drop unused columns to faciliate widening
  
  ggplot(aes(x=PhenoType, y = GO, color = NES, size = -log10(padj))) + 
  geom_point() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'NES')



plot_grid(ggtree_plot, dotplot, nrow = 1, rel_widths = c(0.5,2), align = 'h')

#### How do we do better? ####
dotplot <- GO_cluster %>% filter(GO %in% markers) %>% 
  mutate( GO = factor(GO, levels = clust$labels[clust$order])) %>% 
  #filter(count > 0, `% Expressing` > 1) %>% 
  ggplot(aes(x=PhenoType, y = GO, color = NES, size = -log10(padj))) + 
  geom_point() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'NES')

plot_grid(ggtree_plot, NULL, dotplot, nrow = 1, rel_widths = c(0.5,-0.05, 2), align = 'h')


## Two more tweak options if you are having trouble: 
ggtree_plot_yset <- ggtree_plot + ylim2(dotplot)
plot_grid(ggtree_plot_yset, NULL, dotplot, nrow = 1, rel_widths = c(0.5,-0.05, 2), align = 'h')

## One more adjust
dotplot <- GO_cluster %>% filter(GO %in% markers) %>% 
  mutate( GO = factor(GO, levels = clust$labels[clust$order])) %>% 
  #filter(count > 0, `% Expressing` > 1) %>% 
  ggplot(aes(x=PhenoType, y = GO, color = NES, size = -log10(padj))) + 
  geom_point() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'NES') +
  scale_y_discrete(position = "right")
#################################################

plot_grid(ggtree_plot, NULL, dotplot, nrow = 1, rel_widths = c(0.5,-0.05, 2), align = 'h')

#### Moonshot ####
# make data square to calculate euclidean distance
mat <- GO_cluster %>% 
  filter(GO %in% markers) %>% 
  select(PhenoType, NES, GO) %>%  # drop unused columns to faciliate widening
  pivot_wider(names_from = PhenoType, values_from = NES) %>% 
  data.frame() # make df as tibbles -> matrix annoying
row.names(mat) <- mat$GO  # put GO in `row`
mat <- mat[,-1] #drop GO column as now in rows
mat[is.na(mat)] <- 0
v_clust <- hclust(dist(mat %>% as.matrix() %>% t())) # hclust with distance matrix
############ NOTICE THE t() above)

ddgram_col <- as.dendrogram(v_clust)
ggtree_plot_col <- ggtree(ddgram_col) + layout_dendrogram()

dotplot <- GO_cluster %>% filter(GO %in% markers) %>% 
  mutate(GO = factor(GO, levels = clust$labels[clust$order]),
         PhenoType = factor(PhenoType, levels = v_clust$labels[v_clust$order])) %>% 
  filter(abs(NES) > 0, -log10(padj) > 1) %>% 
  ggplot(aes(x=PhenoType, y = GO, color = NES, size = -log10(padj))) + 
  geom_point() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  #scale_color_gradientn(colours = viridis::viridis(20), limits = c(-4,4), oob = scales::squish, name = 'NES') +
  scale_colour_gradient2(low = "#2d76e3", mid = "white", high = "#e32dc2",guide = "colourbar",midpoint = 0)+
  scale_y_discrete(position = "right")
#################################################
ggtree_plot_col <- ggtree_plot_col + xlim2(dotplot)

ggtree_plot <- ggtree_plot + ylim2(dotplot)

library(RColorBrewer)
colourCount = length(unique(GO_cluster$PhenoType))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
labels <- ggplot(GO_cluster %>% 
                   mutate(Pathway = Representative,
                          PhenoType = factor(PhenoType, levels = v_clust$labels[v_clust$order])), 
                 aes(x = PhenoType, y = 1, fill = PhenoType)) + 
  geom_tile() + 
  scale_fill_manual(values = getPalette(colourCount)) + 
  theme_nothing() +
  xlim2(dotplot)

labels
legend <- plot_grid(get_legend(labels + theme(legend.position="bottom")))


# mat2 <- data.frame(GO=row.names(mat),mat)
# mat2 <- left_join(mat2,GSEA_Fin.df2,by="GO")

colourCountY = length(unique(GO_cluster$Representative))
getPalette = colorRampPalette(brewer.pal(9, "Set3"))
labelsY <- ggplot(GO_cluster %>% 
                 mutate(GO = factor(GO, levels = clust$labels[clust$order]),
                        SumPathway = factor(Representative)), 
                  
                 aes(x = 1, y = GO, fill = SumPathway)) + 
  geom_tile() + 
  scale_fill_manual(values = getPalette(colourCountY)) + 
  theme_nothing() +
  ylim2(dotplot)

# Try
# labelsY <- ggplot(mat2 %>% 
#                     mutate(GO = factor(GO, levels = clust$labels[clust$order]),
#                            SumPathway = factor(Representative)), 
#                   
#                   aes(x = 1, y = GO, fill = SumPathway)) + 
#   geom_tile() + 
#   scale_fill_manual(values = getPalette(colourCountY)) + 
#   theme_nothing() +
#   ylim2(dotplot)

labelsY

legendY <- plot_grid(get_legend(labelsY + theme(legend.position="right"))) #bottom

plot_spacer() + plot_spacer() + ggtree_plot_col +plot_spacer() +
 plot_spacer() + plot_spacer() + labels +  plot_spacer() +
  plot_spacer() + plot_spacer() + plot_spacer() +plot_spacer() +
  ggtree_plot + labelsY + dotplot + plot_spacer() + 
  plot_spacer() + plot_spacer() + legend + legendY+
  plot_layout(ncol = 4, widths = c(0.7, -0.1, 4,1), heights = c(0.9, 0.1, -0.1, 4, 1))

legendY <- plot_grid(get_legend(labelsY + theme(legend.position="right")))+plot_layout(nrow =1)
legendY


pdf(file = paste0(Save.Path,"/PBMC_GSEA_Bubble_Annotation.pdf"),width = 17, height = 20 )
# GO term
plot_spacer() + plot_spacer() + ggtree_plot_col +plot_spacer() +
  plot_spacer() + plot_spacer() + labels +  plot_spacer() +
  plot_spacer() + plot_spacer() + plot_spacer() +plot_spacer() +
  ggtree_plot + labelsY + dotplot+
  theme(axis.title.y = element_blank(),
        #axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + plot_spacer() +
  plot_layout(ncol = 4, widths = c(0.7, -0.1, 4,1), heights = c(0.9, 0.1, -0.1, 30, 1))

# Base
plot_spacer() + plot_spacer() + ggtree_plot_col +plot_spacer() +
  plot_spacer() + plot_spacer() + plot_spacer() +  plot_spacer() +
  plot_spacer() + plot_spacer() + plot_spacer() +plot_spacer() +
  ggtree_plot + labelsY + dotplot+
  theme(axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + plot_spacer() +
  plot_layout(ncol = 4, widths = c(0.7, -0.1, 4,1), heights = c(0.9, 0.1, -0.1, 30, 1))

legendY
dev.off()

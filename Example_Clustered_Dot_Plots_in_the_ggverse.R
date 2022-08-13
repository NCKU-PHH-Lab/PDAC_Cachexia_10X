# https://davemcg.github.io/post/lets-plot-scrna-dotplots/

library(tidyverse)
library(ggdendro)
library(cowplot)
library(ggtree)
library(patchwork) 

#### Load data #####
gene_cluster <- read_tsv('https://github.com/davemcg/davemcg.github.io/raw/master/content/post/scRNA_dotplot_data.tsv.gz')
gene_cluster %>% sample_n(5)


#### Ori Dotplot ####
markers <- gene_cluster$Gene %>% unique()

gene_cluster %>% filter(Gene %in% markers) %>% 
  mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100) %>% 
  ggplot(aes(x=cluster, y = Gene, color = count, size = `% Expressing`)) + 
  geom_point() 


#### Remove dots where there is zero (or near zero expression) ####
gene_cluster %>% filter(Gene %in% markers) %>% 
  mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100) %>% 
  filter(count > 0, `% Expressing` > 1) %>% 
  ggplot(aes(x=cluster, y = Gene, color = count, size = `% Expressing`)) + 
  geom_point() 


#### Better color, better theme, rotate x axis labels ####
gene_cluster %>% filter(Gene %in% markers) %>% 
  mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100) %>% 
  filter(count > 0, `% Expressing` > 1) %>% 
  ggplot(aes(x=cluster, y = Gene, color = count, size = `% Expressing`)) + 
  geom_point() + 
  scale_color_viridis_c(name = 'log2 (count + 1)') + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) 

#### Tweak color scaling ####
gene_cluster %>% filter(Gene %in% markers) %>% 
  mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100) %>% 
  filter(count > 0, `% Expressing` > 1) %>% 
  ggplot(aes(x=cluster, y = Gene, color = count, size = `% Expressing`)) + 
  geom_point() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)')

#### Hey look: ggtree ####
# make data square to calculate euclidean distance
mat <- gene_cluster %>% 
  filter(Gene %in% markers) %>% 
  select(-cell_ct, -cell_exp_ct, -Group) %>%  # drop unused columns to faciliate widening
  pivot_wider(names_from = cluster, values_from = count) %>% 
  data.frame() # make df as tibbles -> matrix annoying
row.names(mat) <- mat$Gene  # put gene in `row`
mat <- mat[,-1] #drop gene column as now in rows
clust <- hclust(dist(mat %>% as.matrix())) # hclust with distance matrix


ddgram <- as.dendrogram(clust) # create dendrogram
ggtree_plot <- ggtree::ggtree(ddgram)
ggtree_plot

#### Let’s glue them together with cowplot ####
dotplot <- gene_cluster %>% filter(Gene %in% markers) %>% 
  mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100) %>% 
  filter(count > 0, `% Expressing` > 1) %>% 
  ggplot(aes(x=cluster, y = Gene, color = count, size = `% Expressing`)) + 
  geom_point() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)')



plot_grid(ggtree_plot, dotplot, nrow = 1, rel_widths = c(0.5,2), align = 'h')

#### How do we do better? ####
dotplot <- gene_cluster %>% filter(Gene %in% markers) %>% 
  mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100,
         Gene = factor(Gene, levels = clust$labels[clust$order])) %>% 
  #filter(count > 0, `% Expressing` > 1) %>% 
  ggplot(aes(x=cluster, y = Gene, color = count, size = `% Expressing`)) + 
  geom_point() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)')

plot_grid(ggtree_plot, NULL, dotplot, nrow = 1, rel_widths = c(0.5,-0.05, 2), align = 'h')


## Two more tweak options if you are having trouble: 
ggtree_plot_yset <- ggtree_plot + ylim2(dotplot)
plot_grid(ggtree_plot_yset, NULL, dotplot, nrow = 1, rel_widths = c(0.5,-0.05, 2), align = 'h')

## One more adjust
dotplot <- gene_cluster %>% filter(Gene %in% markers) %>% 
  mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100,
         Gene = factor(Gene, levels = clust$labels[clust$order])) %>% 
  filter(count > 0, `% Expressing` > 1) %>% 
  ggplot(aes(x=cluster, y = Gene, color = count, size = `% Expressing`)) + 
  geom_point() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)') +
  scale_y_discrete(position = "right")
#################################################

plot_grid(ggtree_plot, NULL, dotplot, nrow = 1, rel_widths = c(0.5,-0.05, 2), align = 'h')

#### Moonshot ####
# make data square to calculate euclidean distance
mat <- gene_cluster %>% 
  filter(Gene %in% markers) %>% 
  select(-cell_ct, -cell_exp_ct, -Group) %>%  # drop unused columns to faciliate widening
  pivot_wider(names_from = cluster, values_from = count) %>% 
  data.frame() # make df as tibbles -> matrix annoying
row.names(mat) <- mat$Gene  # put gene in `row`
mat <- mat[,-1] #drop gene column as now in rows
v_clust <- hclust(dist(mat %>% as.matrix() %>% t())) # hclust with distance matrix
############ NOTICE THE t() above)

ddgram_col <- as.dendrogram(v_clust)
ggtree_plot_col <- ggtree(ddgram_col) + layout_dendrogram()

dotplot <- gene_cluster %>% filter(Gene %in% markers) %>% 
  mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100,
         Gene = factor(Gene, levels = clust$labels[clust$order]),
         cluster = factor(cluster, levels = v_clust$labels[v_clust$order])) %>% 
  filter(count > 0, `% Expressing` > 1) %>% 
  ggplot(aes(x=cluster, y = Gene, color = count, size = `% Expressing`)) + 
  geom_point() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)') +
  scale_y_discrete(position = "right")
#################################################
ggtree_plot_col <- ggtree_plot_col + xlim2(dotplot)

ggtree_plot <- ggtree_plot + ylim2(dotplot)

labels <- ggplot(gene_cluster %>% 
                   mutate(`Cell Type` = Group,
                          cluster = factor(cluster, levels = v_clust$labels[v_clust$order])), 
                 aes(x = cluster, y = 1, fill = `Cell Type`)) + 
  geom_tile() + 
  scale_fill_brewer(palette = 'Set1') + 
  theme_nothing() +
  xlim2(dotplot)

legend <- plot_grid(get_legend(labels + theme(legend.position="bottom")))

plot_spacer() + plot_spacer() + ggtree_plot_col +
  plot_spacer() + plot_spacer() + labels + 
  plot_spacer() + plot_spacer() + plot_spacer() +
  ggtree_plot + plot_spacer() + dotplot + 
  plot_spacer() + plot_spacer() + legend + 
  plot_layout(ncol = 3, widths = c(0.7, -0.1, 4), heights = c(0.9, 0.1, -0.1, 4, 1))

rm(list = ls(all.names = TRUE))

library(ggplot2)
library(cowplot)

# Spearman correlations of molecular evolution metrics vs % intact singletons
pangenome <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/pangenome_and_related_metrics_filt.tsv.gz",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

dnds_vs_si <- ggplot(data = pangenome,
                            aes(x = median_MG94_dnds,
                                y = mean_percent_singletons_per9)) +
  geom_point() +
  theme_bw() +
  ylab(expression(paste('Mean percent intact singletons (', s[i], ')', sep = ''))) +
  xlab('Median dN/dS across core genes') +
  ggpubr::stat_cor(method = 'spearman', label.y.npc = 'bottom') +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10')

ds_vs_si <- ggplot(data = pangenome,
                     aes(x = median_pairwise_ds,
                         y = mean_percent_singletons_per9)) +
  geom_point() +
  theme_bw() +
  ylab(expression(paste('Mean percent intact singletons (', s[i], ')', sep = ''))) +
  xlab('Median pairwise dS across core genes') +
  ggpubr::stat_cor(method = 'spearman') +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10')

ds_vs_dnds <- ggplot(data = pangenome,
                     aes(x = median_pairwise_ds,
                         y = median_MG94_dnds)) +
  geom_point() +
  theme_bw() +
  ylab('Median dN/dS across core genes') +
  xlab('Median pairwise dS across core genes') +
  ggpubr::stat_cor(method = 'spearman', label.y.npc = 'bottom') +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10')

combined <- plot_grid(dnds_vs_si,
                      ds_vs_si,
                      ds_vs_dnds,
                      labels = c('a', 'b', 'c'))

ggsave(filename = '/home/gdouglas/scripts/pangenome_pseudogene_null/display_items/extended_figure_si_vs_evolultion_metrics.png',
       plot = combined,
       device = 'png',
       dpi = 400,
       width = 10,
       height = 8)
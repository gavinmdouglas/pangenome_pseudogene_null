rm(list = ls(all.names = TRUE))

library(ggplot2)
library(cowplot)

# Spearman correlation of dnds vs syn pi
pangenome <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/pangenome_and_related_metrics_filt.tsv.gz",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

ds_vs_dnds <- ggplot(data = pangenome,
                     aes(x = median_pairwise_ds,
                         y = median_MG94_dnds)) +
  geom_point() +
  theme_bw() +
  ylab('Median dN/dS across core genes') +
  xlab(expression(paste('Median ', pi[s], ' across core genes', sep = ''))) +
  ggpubr::stat_cor(aes(label = paste(r.label, gsub("p", "P", ..p.label..), sep = "~`,`~")),
                   method = 'spearman', label.y.npc = 'bottom', cor.coef.name = 'rho') +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10')

ggsave(filename = '/home/gdouglas/scripts/pangenome_pseudogene_null/display_items/extended_figure_dNdS_vs_syn_pi.pdf',
       plot = ds_vs_dnds,
       device = 'pdf',
       dpi = 400,
       width = 6,
       height = 6)
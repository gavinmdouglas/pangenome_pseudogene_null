rm(list = ls(all.names = TRUE))

library(ggplot2)
library(cowplot)

# Spearman correlations of molecular evolution metrics vs si and si/sp
pangenome <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/pangenome_and_related_metrics_filt.tsv.gz",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

dnds_vs_si <- ggplot(data = pangenome,
                            aes(x = median_MG94_dnds,
                                y = mean_percent_singletons_per9)) +
  geom_point() +
  theme_bw() +
  ylab(expression(paste('Mean percent intact singletons (', s[i], ')', sep = ''))) +
  xlab('Median dN/dS across core genes') +
  ggpubr::stat_cor(aes(label = paste(r.label, gsub("p", "P", ..p.label..), sep = "~`,`~")),
                   method = 'spearman', label.y.npc = 'bottom', cor.coef.name = 'rho') +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10')

ds_vs_si <- ggplot(data = pangenome,
                     aes(x = median_pairwise_ds,
                         y = mean_percent_singletons_per9)) +
  geom_point() +
  theme_bw() +
  ylab(expression(paste('Mean percent intact singletons (', s[i], ')', sep = ''))) +
  xlab(expression(paste('Median ', pi[s], ' across core genes', sep = ''))) +
  ggpubr::stat_cor(aes(label = paste(r.label, gsub("p", "P", ..p.label..), sep = "~`,`~")),
                   method = 'spearman', label.y.npc = 'bottom', cor.coef.name = 'rho') +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10')

dnds_vs_si_sp <- ggplot(data = pangenome,
                     aes(x = median_MG94_dnds,
                         y = si_sp)) +
  geom_point() +
  theme_bw() +
  ylab(expression(paste('Mean percent singleton ratio (', s[i], '/', s[p], ')', sep = ''))) +
  xlab('Median dN/dS across core genes') +
  ggpubr::stat_cor(aes(label = paste(r.label, gsub("p", "P", ..p.label..), sep = "~`,`~")),
                   method = 'spearman', label.y.npc = 'bottom', cor.coef.name = 'rho') +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10')

ds_vs_si_sp <- ggplot(data = pangenome,
                   aes(x = median_pairwise_ds,
                       y = si_sp)) +
  geom_point() +
  theme_bw() +
  ylab(expression(paste('Mean percent singleton ratio (', s[i], '/', s[p], ')', sep = ''))) +
  xlab(expression(paste('Median ', pi[s], ' across core genes', sep = ''))) +
  ggpubr::stat_cor(aes(label = paste(r.label, gsub("p", "P", ..p.label..), sep = "~`,`~")),
                   method = 'spearman', label.y.npc = 'bottom', cor.coef.name = 'rho') +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10')


combined <- plot_grid(dnds_vs_si,
                      ds_vs_si,
                      dnds_vs_si_sp,
                      ds_vs_si_sp,
                      labels = c('a', 'b', 'c', 'd'))

ggsave(filename = '/home/gdouglas/scripts/pangenome_pseudogene_null/display_items/main_si_and_sisp_vs_evolution_metrics.pdf',
       plot = combined,
       device = 'pdf',
       dpi = 600,
       width = 10,
       height = 8)

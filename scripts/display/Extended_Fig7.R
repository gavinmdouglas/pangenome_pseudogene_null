rm(list = ls(all.names = TRUE))

library(ggplot2)
library(cowplot)

additional <- read.table(file = "/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/genomic_fluidity_and_related_additional_subsamples.tsv",
                         header = TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

# Spearman correlations of molecular evolution metrics vs % intact singletons
pangenome <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/pangenome_and_related_metrics.tsv.gz",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

pangenome <- pangenome[rownames(additional), ]

all.equal(pangenome[, c("mean_num_genes", "mean_num_pseudo", "num_panaroo_genomes", "num_pseudofinder_genomes")],
          additional[, c("mean_num_genes", "mean_num_pseudo", "num_panaroo_genomes", "num_pseudofinder_genomes")])

additional <- additional[, -which(colnames(additional) %in% c("mean_num_genes", "mean_num_pseudo", "num_panaroo_genomes", "num_pseudofinder_genomes"))]

pangenome <- cbind(pangenome, additional)

pangenome$mean_percent_singletons_per3 <- 100 * (pangenome$mean_num_singletons_per3 / pangenome$mean_num_genes)
pangenome$mean_percent_singletons_pseudo_per3 <- 100 * (pangenome$mean_num_singletons_pseudo_per3 / pangenome$mean_num_pseudo)
pangenome$si_sp_per3 <- pangenome$mean_percent_singletons_per3 / pangenome$mean_percent_singletons_pseudo_per3

pangenome$mean_percent_singletons_per20 <- 100 * (pangenome$mean_num_singletons_per20 / pangenome$mean_num_genes)
pangenome$mean_percent_singletons_pseudo_per20 <- 100 * (pangenome$mean_num_singletons_pseudo_per20 / pangenome$mean_num_pseudo)
pangenome$si_sp_per20 <- pangenome$mean_percent_singletons_per20 / pangenome$mean_percent_singletons_pseudo_per20

dnds_vs_si_per3 <- ggplot(data = pangenome,
                            aes(x = dnds,
                                y = mean_percent_singletons_per3)) +
  geom_point() +
  theme_bw() +
  ylab(expression(paste(s[i], ' (three genomes)', sep = ''))) +
  xlab(expression(paste('dN/dS (', log[10], ')', sep = ''))) +
  ggpubr::stat_cor(method = 'spearman', label.x.npc = 'centre', label.y.npc = 'top', cor.coef.name = 'rho') +
  scale_x_continuous(trans='log10')

dnds_vs_si_per20 <- ggplot(data = pangenome,
                          aes(x = dnds,
                              y = mean_percent_singletons_per20)) +
  geom_point() +
  theme_bw() +
  ylab(expression(paste(s[i], ' (20 genomes)', sep = ''))) +
  xlab(expression(paste('dN/dS (', log[10], ')', sep = ''))) +
  ggpubr::stat_cor(method = 'spearman', label.x.npc = 'centre', label.y.npc = 'top', cor.coef.name = 'rho') +
  scale_x_continuous(trans='log10')


dnds_vs_sisp_per3 <- ggplot(data = pangenome,
                        aes(x = dnds,
                            y = si_sp_per3)) +
  geom_point() +
  theme_bw() +
  ylab(expression(paste(s[i], '/', s[p], ' (three genomes)', sep = ''))) +
  xlab(expression(paste('dN/dS (', log[10], ')', sep = ''))) +
  ggpubr::stat_cor(method = 'spearman', cor.coef.name = 'rho', label.x.npc = 'centre', label.y.npc = 'top') +
  scale_x_continuous(trans='log10')

dnds_vs_sisp_per20 <- ggplot(data = pangenome,
                         aes(x = dnds,
                             y = si_sp_per20)) +
  geom_point() +
  theme_bw() +
  ylab(expression(paste(s[i], '/', s[p], ' (20 genomes)', sep = ''))) +
  xlab(expression(paste('dN/dS (', log[10], ')', sep = ''))) +
  ggpubr::stat_cor(method = 'spearman', cor.coef.name = 'rho', label.x.npc = 'centre', label.y.npc = 'top') +
  scale_x_continuous(trans='log10')


combined <- plot_grid(dnds_vs_si_per3, dnds_vs_si_per20,
                      dnds_vs_sisp_per3, dnds_vs_sisp_per20,
                      labels = c('a', 'b', 'c', 'd'),
                      nrow = 2, ncol = 2)

ggsave(filename = '/home/gdouglas/scripts/pangenome_pseudogene_null/display_items/Douglas_ED_Fig7.pdf',
       plot = combined,
       device = 'pdf',
       dpi = 400,
       width = 10,
       height = 8)

# Write out source data:
source_out <- pangenome[, c('dnds',
                            'mean_percent_singletons_per3', 'mean_percent_singletons_pseudo_per3', 'si_sp_per3',
                            'mean_percent_singletons_per20', 'mean_percent_singletons_pseudo_per20', 'si_sp_per20')]
orig_col <- colnames(source_out)
source_out$species <- rownames(source_out)
source_out <- source_out[, c('species', orig_col)]
write.table(x = source_out,
            file = "/home/gdouglas/scripts/pangenome_pseudogene_null/display_source_data/ED_Fig7.tsv",
            col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)

rm(list = ls(all.names = TRUE))

library(ggplot2)
library(cowplot)
library(ppcor)

options(scipen = 10)

# Correlations of dN/dS vs si, si/sp, genomic fluidity, and the mean # of genes.
pangenome <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/pangenome_and_related_metrics_filt.tsv.gz",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

# Source functions to get Spearman correlation utility script.
source('/home/gdouglas/scripts/pangenome_pseudogene_null/scripts/functions.R')


# Mean number of genes
dnds_vs_mean_num_genes_cor_string <- spearman_cor_string(x = pangenome$dnds,
                                                          y = pangenome$mean_num_genes)

dnds_vs_mean_num_genes <- ggplot(data = pangenome,
                                 aes(x = dnds,
                                     y = mean_num_genes)) +
                          geom_point() +
                          theme_bw() +
                          ylab('Mean number of genes') +
                          xlab(expression(paste('dN/dS (', log[10], ')', sep = ''))) +
                          annotate(
                            geom = "text",
                            x = 1,
                            y = 7500,
                            label = dnds_vs_mean_num_genes_cor_string,
                            parse = TRUE) +
                            scale_x_continuous(trans='log10')
  
# Genomic fluidity
dnds_vs_fluidity_cor_string <- spearman_cor_string(x = pangenome$dnds,
                                                   y = pangenome$genomic_fluidity)

dnds_vs_fluidity <- ggplot(data = pangenome,
                                 aes(x = dnds,
                                     y = genomic_fluidity)) +
                            geom_point() +
                            theme_bw() +
                            ylab('Genomic fluidity') +
                            xlab(expression(paste('dN/dS (', log[10], ')', sep = ''))) +
                            annotate(
                              geom = "text",
                              x = 1,
                              y = 0.2,
                              label = dnds_vs_fluidity_cor_string,
                              parse = TRUE) +
                            scale_x_continuous(trans='log10')

### Si
dnds_vs_si_cor_string <- spearman_cor_string(x = pangenome$dnds,
                                             y = pangenome$mean_percent_singletons_per9)

dnds_vs_si <- ggplot(data = pangenome,
                           aes(x = dnds,
                               y = mean_percent_singletons_per9)) +
              geom_point() +
              theme_bw() +
  ylab(expression(paste('Mean percent singleton intact genes (', s[i], ')', sep = ''))) +
              xlab(expression(paste('dN/dS (', log[10], ')', sep = ''))) +
              annotate(
                geom = "text",
                x = 1,
                y = 9,
                label = dnds_vs_si_cor_string,
                parse = TRUE) +
              scale_x_continuous(trans='log10')


# si/sp
dnds_vs_sisp_cor_string <- spearman_cor_string(x = pangenome$dnds,
                                               y = pangenome$si_sp)

dnds_vs_sisp <- ggplot(data = pangenome,
                     aes(x = dnds,
                         y = si_sp)) +
                    geom_point() +
                    theme_bw() +
                    ylab(expression(paste('Mean percent singleton ratio (', s[i], '/', s[p], ')', sep = ''))) +
                    xlab(expression(paste('dN/dS (', log[10], ')', sep = ''))) +
                    annotate(
                      geom = "text",
                      x = 1,
                      y = 0.4,
                      label = dnds_vs_sisp_cor_string,
                      parse = TRUE) +
                    scale_x_continuous(trans='log10')

combined <- plot_grid(dnds_vs_mean_num_genes,
                      dnds_vs_fluidity,
                      dnds_vs_si,
                      dnds_vs_sisp,
                      labels = c('a', 'b', 'c', 'd'))

ggsave(filename = '/home/gdouglas/scripts/pangenome_pseudogene_null/display_items/Douglas_Fig4.pdf',
       plot = combined,
       device = 'pdf',
       dpi = 600,
       width = 10,
       height = 8)

options(scipen = 0)

# Get exact P-values to report in figure legend:
cor.test(pangenome$dnds, pangenome$genomic_fluidity, method = 'spearman')
cor.test(pangenome$dnds, pangenome$mean_percent_singletons_per9, method = 'spearman')
cor.test(pangenome$dnds, pangenome$si_sp, method = 'spearman')

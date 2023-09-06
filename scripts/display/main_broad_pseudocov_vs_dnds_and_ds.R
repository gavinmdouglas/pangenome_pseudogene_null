rm(list = ls(all.names = TRUE))

library(cowplot)
library(ggplot2)

# Source functions to get Spearman correlation utility script.
source('/home/gdouglas/scripts/pangenome_pseudogene_null/scripts/functions.R')

pangenome <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/pangenome_and_related_metrics_filt.tsv.gz",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

# Compare % pseudogene coverage vs dN/dS (and fill in si/sp)
pseudogene_percent_mean_by_species <- read.table('/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/element_info/element_percent_coverage/pseudogene_mean_percent_coverage_by_species.tsv.gz',
                                                 header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

pangenome$mean_percent_coverage_pseudogenes <- pseudogene_percent_mean_by_species[rownames(pangenome), 'mean_percent']

ppcor_dnds_vs_coverage_string <- spearman_ppcor_string(x = pangenome$mean_percent_coverage_pseudogenes,
                                                       y = pangenome$dnds,
                                                       z = pangenome$ds)

cor_dnds_vs_coverage_string <- spearman_cor_string(x = pangenome$mean_percent_coverage_pseudogenes,
                                                   y = pangenome$dnds)

dnds_vs_mean_percent_coverage_pseudogenes <- ggplot(data = pangenome,
                                                    aes(y = mean_percent_coverage_pseudogenes,
                                                        x = dnds,
                                                        colour = log10(si_sp))) +
                                                        labs(colour = expression(paste(log[10], '(', s[i], '/', s[p], ')', sep = ''))) +
                                                        geom_point() +
                                                        theme_bw() +
                                                        ylab('Mean\npercent\nof genome\ncovered by\npseudogenes') +
                                                        xlab(expression(paste('dN/dS (', log[10], ')', sep = ''))) +
                                                        theme(axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5)) +
                                                        scale_x_continuous(trans='log10') +
                                                        annotate(
                                                          geom = "text",
                                                          y = 7,
                                                          x = 1,
                                                          label = ppcor_dnds_vs_coverage_string,
                                                          parse = TRUE) +
                                                          annotate(
                                                            geom = "text",
                                                            y = 6.5,
                                                            x = 1,
                                                            label = cor_dnds_vs_coverage_string,
                                                            parse = TRUE)

# Similar, but for dS:
cor_ds_vs_coverage_string <- spearman_cor_string(x = pangenome$mean_percent_coverage_pseudogenes,
                                                 y = pangenome$ds)

ds_vs_mean_percent_coverage_pseudogenes <- ggplot(data = pangenome,
                                                    aes(y = mean_percent_coverage_pseudogenes,
                                                        x = ds,
                                                        colour = log10(si_sp))) +
                                                        labs(colour = expression(paste(log[10], '(', s[i], '/', s[p], ')', sep = ''))) +
                                                        geom_point() +
                                                        theme_bw() +
                                                        ylab('Mean\npercent\nof genome\ncovered by\npseudogenes') +
                                                        xlab('dS') +
                                                        theme(axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5)) +
                                                        annotate(
                                                          geom = "text",
                                                          y = 7,
                                                          x = 0.15,
                                                          label = cor_ds_vs_coverage_string,
                                                          parse = TRUE)

combined_plot <- cowplot::plot_grid(dnds_vs_mean_percent_coverage_pseudogenes,
                                    ds_vs_mean_percent_coverage_pseudogenes,
                                    labels = c('a', 'b'),
                                    nrow = 2)

ggsave(filename = '/home/gdouglas/scripts/pangenome_pseudogene_null/display_items/main_figure_broad_ds_and_dnds_vs_pseudocov.pdf',
       plot = combined_plot,
       device = 'pdf',
       dpi = 400,
       width = 8,
       height = 10)

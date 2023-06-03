rm(list = ls(all.names = TRUE))

library(ggplot2)

pangenome <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/pangenome_and_related_metrics_filt.tsv.gz",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

dnds_vs_ds <- ggplot(data = pangenome,
                      aes(y = ds,
                          x = dnds,
                          colour = log10(si_sp))) +
                      geom_point() +
                      theme_bw() +
                      ylab('        dS  ') +
                      xlab('dN/dS') +
                      labs(colour = expression(paste(log[10], '(', s[i], '/', s[p], ')', sep = ''))) +
                      theme(axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5)) +
                      scale_x_continuous(trans='log10') +
                      scale_y_continuous(trans='log10')


# Run Partial Spearman's rho of % pseudogene coverage vs dN/dS (while controlling for dS).
pseudogene_percent_mean_by_species <- read.table('/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/element_info/element_percent_coverage/pseudogene_mean_percent_coverage_by_species.tsv.gz',
                                                 header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

pangenome$mean_percent_coverage_pseudogenes <- pseudogene_percent_mean_by_species[rownames(pangenome), 'mean_percent']

ppcor_dnds_vs_coverage <- ppcor::pcor.test(x = pangenome$mean_percent_coverage_pseudogenes,
                                           y = pangenome$dnds,
                                           z = pangenome$ds,
                                           method = 'spearman')

cor_dnds_vs_coverage <- cor.test(x = pangenome$mean_percent_coverage_pseudogenes,
                                 y = pangenome$dnds,
                                 method = 'spearman')

ppcor_dnds_vs_coverage_estimate <- format(round(ppcor_dnds_vs_coverage$estimate, digits=4), nsmall = 4)
ppcor_dnds_vs_coverage_p <- format(round(ppcor_dnds_vs_coverage$p.value, digits=4), nsmall = 4)

ppcor_dnds_vs_coverage_string <- paste0('Partial ~ rho ~', "'= '*", ppcor_dnds_vs_coverage_estimate,"*','", '~ italic(P) ~', "'= '*", ppcor_dnds_vs_coverage_p)

dnds_vs_mean_percent_coverage_pseudogenes <- ggplot(data = pangenome,
                                                    aes(y = mean_percent_coverage_pseudogenes,
                                                        x = dnds,
                                                        colour = log10(si_sp))) +
                                                        geom_point() +
                                                        theme_bw() +
                                                        ylab('Mean\npercent\ngenome\ncovered by\npseudogenes') +
                                                        xlab('dN/dS') +
                                                        labs(colour = expression(paste(log[10], '(', s[i], '/', s[p], ')', sep = ''))) +
                                                        theme(axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5)) +
                                                        scale_x_continuous(trans='log10') +
                                                        scale_y_continuous(trans='log10') +
                                                        annotate(
                                                          geom = "text",
                                                          y = 5,
                                                          x = 0.8,
                                                          label = ppcor_dnds_vs_coverage_string,
                                                          parse = TRUE)

combined_plot <- cowplot::plot_grid(dnds_vs_ds,
                                    dnds_vs_mean_percent_coverage_pseudogenes,
                                    labels = c('a', 'b'),
                                    nrow = 2)

ggsave(filename = '/home/gdouglas/scripts/pangenome_pseudogene_null/display_items/extended_figure_broad_dnds_vs_ds_and_pseudocov.pdf',
       plot = combined_plot,
       device = 'pdf',
       dpi = 400,
       width = 6,
       height = 8)

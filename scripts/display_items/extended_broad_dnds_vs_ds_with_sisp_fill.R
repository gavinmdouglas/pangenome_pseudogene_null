rm(list = ls(all.names = TRUE))

library(ggplot2)

pangenome <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/pangenome_and_related_metrics_filt.tsv.gz",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

dnds_vs_ds <- ggplot(data = pangenome,
                      aes(x = ds,
                          y = dnds,
                          colour = log10(si_sp))) +
                      geom_point() +
                      theme_bw() +
                      xlab('dS') +
                      ylab('dN/dS') +
                      labs(colour = expression(paste(log[10], '(', s[i], '/', s[p], ')', sep = ''))) +
                      theme(axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5))

ggsave(filename = '/home/gdouglas/scripts/pangenome_pseudogene_null/display_items/extended_figure_dnds_vs_ds_with_sisp_fill.pdf',
       plot = dnds_vs_ds,
       device = 'pdf',
       dpi = 400,
       width = 6,
       height = 4)

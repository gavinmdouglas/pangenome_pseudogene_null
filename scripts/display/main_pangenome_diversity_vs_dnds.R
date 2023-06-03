rm(list = ls(all.names = TRUE))

library(ggplot2)
library(cowplot)
library(ppcor)

options(scipen = 10)

# *Partial* Spearman correlations of dN/dS vs si, si/sp, genomic fluidity, and the mean # of genes.

pangenome <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/pangenome_and_related_metrics_filt.tsv.gz",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

### Mean num genes
dnds_vs_mean_num_genes_pcor <- pcor.test(x = pangenome$dnds,
                                         y = pangenome$mean_num_genes,
                                         z = pangenome$ds,
                                         method = 'spearman')

dnds_vs_mean_num_genes_pcor_estimate <- format(round(dnds_vs_mean_num_genes_pcor$estimate, digits=4), nsmall = 4)
dnds_vs_mean_num_genes_pcor_p <- format(round(dnds_vs_mean_num_genes_pcor$p.value, digits=4), nsmall = 4)

dnds_vs_mean_num_genes_pcor_string <- paste0('Partial ~ rho ~', "'='", dnds_vs_mean_num_genes_pcor_estimate,"*','", '~ italic(P) ~', "'= '*", dnds_vs_mean_num_genes_pcor_p)

dnds_vs_mean_num_genes <- ggplot(data = pangenome,
                                 aes(x = dnds,
                                     y = mean_num_genes)) +
                          geom_point() +
                          theme_bw() +
                          ylab('Mean number of genes') +
                          xlab('dN/dS') +
                          annotate(
                            geom = "text",
                            x = 0.2,
                            y = 700,
                            label = dnds_vs_mean_num_genes_pcor_string,
                            parse = TRUE) +
                          scale_x_continuous(trans='log10') +
                          scale_y_continuous(trans='log10')
  
### Genomic fluidity
dnds_vs_fluidity_pcor <- pcor.test(x = pangenome$dnds,
                                   y = pangenome$genomic_fluidity,
                                   z = pangenome$ds,
                                   method = 'spearman')

dnds_vs_fluidity_pcor_estimate <- format(round(dnds_vs_fluidity_pcor$estimate, digits=4), nsmall = 4)
dnds_vs_fluidity_pcor_p <- format(round(dnds_vs_fluidity_pcor$p.value, digits=4), nsmall = 4)
# NOTE -- I MANUALLY ADDED A TRAILING ZERO HERE, because I couldn't figure out how to make 'parse' include it.
dnds_vs_fluidity_pcor_string <- paste0('Partial ~ rho ~', "'='", dnds_vs_fluidity_pcor_estimate,"*'0,'", '~ italic(P) ~', "'= '*", dnds_vs_fluidity_pcor_p)

dnds_vs_fluidity <- ggplot(data = pangenome,
                                 aes(x = dnds,
                                     y = genomic_fluidity)) +
                    geom_point() +
                    theme_bw() +
                    ylab('Genomic fluidity') +
                    xlab('dN/dS') +
                    annotate(
                      geom = "text",
                      x = 0.2,
                      y = 0.00025,
                      label = dnds_vs_fluidity_pcor_string,
                      parse = TRUE) +
                    scale_x_continuous(trans='log10') +
                    scale_y_continuous(trans='log10')

### Si
dnds_vs_si_pcor <- pcor.test(x = pangenome$dnds,
                                   y = pangenome$mean_percent_singletons_per9,
                                   z = pangenome$ds,
                                   method = 'spearman')

dnds_vs_si_pcor_estimate <- format(round(dnds_vs_si_pcor$estimate, digits=4), nsmall = 4)
dnds_vs_si_pcor_p <- format(round(dnds_vs_si_pcor$p.value, digits=4), nsmall = 4)
if (dnds_vs_si_pcor_p >= 0.0001) {
  dnds_vs_si_pcor_string <- paste0('Partial ~ rho ~', "'='", dnds_vs_si_pcor_estimate,"*','", '~ italic(P) ~', "'= '*", dnds_vs_si_pcor_p)
} else {
  dnds_vs_si_pcor_string <- paste0('Partial ~ rho ~', "'='", dnds_vs_si_pcor_estimate,"*','", '~ italic(P) ~', "'< '* 0.0001")
}



dnds_vs_si <- ggplot(data = pangenome,
                           aes(x = dnds,
                               y = mean_percent_singletons_per9)) +
              geom_point() +
              theme_bw() +
              ylab(expression(paste('Mean percent intact singletons (', s[i], ')', sep = ''))) +
              xlab('dN/dS') +
              annotate(
                geom = "text",
                x = 0.2,
                y = 0.007,
                label = dnds_vs_si_pcor_string,
                parse = TRUE) +
              scale_x_continuous(trans='log10') +
              scale_y_continuous(trans='log10')



dnds_vs_sisp_pcor <- pcor.test(x = pangenome$dnds,
                                   y = pangenome$si_sp,
                                   z = pangenome$ds,
                                   method = 'spearman')
dnds_vs_sisp_pcor_estimate <- format(round(dnds_vs_sisp_pcor$estimate, digits=4), nsmall = 4)
dnds_vs_sisp_pcor_p <- format(round(dnds_vs_sisp_pcor$p.value, digits=4), nsmall = 4)

if (dnds_vs_sisp_pcor_p >= 0.0001) {
  dnds_vs_sisp_pcor_string <- paste0('Partial ~ rho ~', "'='", dnds_vs_sisp_pcor_estimate,"*','", '~ italic(P) ~', "'= '*", dnds_vs_sisp_pcor_p)
} else {
  dnds_vs_sisp_pcor_string <- paste0('Partial ~ rho ~', "'='", dnds_vs_sisp_pcor_estimate,"*','", '~ italic(P) ~', "'< '* 0.0001")
}

dnds_vs_sisp <- ggplot(data = pangenome,
                           aes(x = dnds,
                               y = si_sp)) +
                geom_point() +
                theme_bw() +
                ylab(expression(paste('Mean percent singleton ratio (', s[i], '/', s[p], ')', sep = ''))) +
                xlab('dN/dS') +
                annotate(
                  geom = "text",
                  x = 0.2,
                  y = 0.001,
                  label = dnds_vs_sisp_pcor_string,
                  parse = TRUE) +
                scale_x_continuous(trans='log10') +
                scale_y_continuous(trans='log10')

combined <- plot_grid(dnds_vs_mean_num_genes,
                      dnds_vs_fluidity,
                      dnds_vs_si,
                      dnds_vs_sisp,
                      labels = c('a', 'b', 'c', 'd'))

ggsave(filename = '/home/gdouglas/scripts/pangenome_pseudogene_null/display_items/main_si_and_sisp_vs_evolution_metrics.pdf',
       plot = combined,
       device = 'pdf',
       dpi = 600,
       width = 10,
       height = 8)

options(scipen = 0)

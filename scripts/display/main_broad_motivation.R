rm(list = ls(all.names = TRUE))

# Provide visual motivation and context for broad-scale analyses:
  # si vs sp
  # dN/dS vs dS
  # si vs dS
  # sp vs dS

library(cowplot)
library(ggplot2)
library(ggrepel)

# Source functions to get Spearman correlation utility script.
source('/home/gdouglas/scripts/pangenome_pseudogene_null/scripts/functions.R')

pangenome <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/pangenome_and_related_metrics_filt.tsv.gz",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

upper_annotation <- "More positive selection driving either\ngain in rare accessory genes or loss\nof pseudogenes?"
lower_annotation <- "Less positive selection driving either\ngain in rare accessory genes or loss\nof pseudogenes? Insufficient time?"

# Highlight a few taxa mention in main text.
pangenome$highlighted <- ''
pangenome['Escherichia_coli', 'highlighted'] <- 'Escherichia coli'
pangenome['Chlamydia_trachomatis', 'highlighted'] <- 'Chlamydia trachomatis'
pangenome['Rickettsia_prowazekii', 'highlighted'] <- 'Rickettsia prowazekii'

si_vs_sp <- ggplot(data = pangenome,
                   aes(x = mean_percent_singletons_pseudo_per9,
                       y = mean_percent_singletons_per9)) +
            geom_point() +
            theme_bw() +
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank()) +
            xlab(expression(paste('Mean ', 'percent ', 'pseudogene ', 'singletons (', s[p], ')', sep = ''))) +
            ylab(expression(paste('Mean ', 'percent ', 'intact ', 'singletons (', s[i], ')', sep = ''))) +
            ggrepel::geom_text_repel(aes(label = highlighted),
                                     col = 'grey25',
                                     hjust = -0.1,
                                     fontface = 'italic',
                                     size = 3,
                                     max.overlaps = 100,
                                     min.segment.length = 0) +
            geom_abline(slope = 3/20, intercept = 0, colour = "black", lty = 2, linewidth = 2, alpha = 0.3) +
            geom_segment(aes(x = 40, y = 5, xend = 56, yend = 2),
                         arrow = arrow(length = unit(0.5, "cm"), type = 'closed'),
                         inherit.aes = FALSE,
                         colour = 'red',
                         linejoin = 'mitre',
                         alpha = 0.005, 
                         linewidth = 5) +
            geom_segment(aes(x = 34, y = 6, xend = 18, yend = 9),
                         arrow = arrow(length = unit(0.5, "cm"), type = 'closed'),
                         inherit.aes = FALSE,
                         colour = 'blue',
                         linejoin = 'mitre',
                         alpha = 0.005, 
                         linewidth = 5) +
            annotate(geom = "text",
                     x = 0,
                     y = 10.5,
                     color = "blue",
                     label = upper_annotation,
                     hjust = 0,
                     parse = FALSE,
                     size = 3) +
            annotate(geom = "text",
                     x = 40,
                     y = 0.4,
                     color = "red",
                     label = lower_annotation,
                     hjust = 0,
                     parse = FALSE,
                     size = 3)


# dN/dS vs dS (coloured by si/sp)

dnds_vs_ds_cor_string <- spearman_cor_string(x = pangenome$dnds, y = pangenome$ds)

dnds_vs_ds <- ggplot(data = pangenome,
                     aes(x = dnds,
                         y = ds,
                         colour = log10(si_sp))) +
                     geom_point() +
                     theme_bw() +
                     ylab('dS') +
                     xlab(expression(paste('dN/dS (', log[10], ')', sep = ''))) +
                     labs(colour = expression(paste(log[10], '(', s[i], '/', s[p], ')', sep = ''))) +
                     theme(axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5)) +
                     scale_x_continuous(trans='log10') +
                      annotate(geom = "text",
                               y = 0.085,
                               x = 1,
                               label = dnds_vs_ds_cor_string,
                               parse = TRUE) +
                     theme(legend.position = c(0.85, 0.75),
                           legend.box.background = element_rect(colour = "black"))

# Compare % singleton intact genes and dS (while display dN/dS too).
cor_ds_vs_si_string <- spearman_cor_string(x = pangenome$mean_percent_singletons_per9, y = pangenome$ds)

ds_vs_si_pseudogenes <- ggplot(data = pangenome,
                               aes(y = mean_percent_singletons_per9,
                                   x = ds,
                                   colour = log10(dnds))) +
                              scale_color_gradient(low='grey70', high='grey30') +
                              geom_point() +
                              theme_bw() +
                              ylab(expression(s[i])) +
                              xlab('dS') +
                              labs(colour = 'dN/dS') +
                              theme(axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5)) +
                              annotate(geom = "text",
                                       y = 1,
                                       x = 0.2,
                                       label = cor_ds_vs_si_string,
                                       parse = TRUE) +
                              xlim(0, 0.25) +
                              theme(legend.position = c(0.85, 0.75),
                                    legend.box.background = element_rect(colour = "black"))


# Compare % singleton pseudogenes and dS (while displaying dN/dS).

cor_ds_vs_sp_string <- spearman_cor_string(x = pangenome$mean_percent_singletons_pseudo_per9, y = pangenome$ds)

ds_vs_sp_pseudogenes <- ggplot(data = pangenome,
                               aes(y = mean_percent_singletons_pseudo_per9,
                                   x = ds,
                                   colour = log10(dnds))) +
                              scale_color_gradient(low='grey70', high='grey30') +
                              geom_point() +
                              theme_bw() +
                              ylab(expression(s[p])) +
                              xlab('dS') +
                              labs(colour = 'dN/dS') +
                              theme(axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5)) +
                              annotate(geom = "text",
                                       y = 6.6666,
                                       x = 0.2,
                                       label = cor_ds_vs_sp_string,
                                       parse = TRUE) +
                              xlim(0, 0.25) +
                              theme(legend.position = c(0.85, 0.75),
                                    legend.box.background = element_rect(colour = "black"))

combined <- plot_grid(si_vs_sp,
                     dnds_vs_ds,
                     ds_vs_si_pseudogenes,
                     ds_vs_sp_pseudogenes,
                     labels = c('a', 'b', 'c', 'd'),
                     nrow = 2)

ggsave(filename = '/home/gdouglas/scripts/pangenome_pseudogene_null/display_items/main_broad_motivation.pdf',
       plot = combined,
       device = 'pdf',
       dpi = 600,
       width = 11,
       height = 9)

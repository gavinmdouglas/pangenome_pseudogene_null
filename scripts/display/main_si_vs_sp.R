rm(list = ls(all.names = TRUE))

library(ggplot2)

pangenome <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/pangenome_and_related_metrics_filt.tsv.gz",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

upper_annotation <- "More positive selection driving either\ngain in rare accessory genes or loss\nof pseudogenes?"
lower_annotation <- "Less positive selection driving either\ngain in rare accessory genes or loss\nof pseudogenes? Insufficient time?"

si_vs_sp <- ggplot(data = pangenome,
                   aes(x = mean_percent_singletons_pseudo_per9,
                       y = mean_percent_singletons_per9)) +
                     # colour = log10(genome_size))) + ### Could add colour per dot indicating genome size, but decided it was too distracting.
            geom_point() +
            scale_colour_gradient(limits = c(5.75, 7.00)) +
            #labs(colour = expression(paste(log[10], '(Genome size)', sep = ''))) + ### Decided against adding in genome size.
            theme_bw() +
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank()) +
            xlab(expression(paste('Mean ', 'percent ', 'pseudogene ', 'singletons (', s[p], ')', sep = ''))) +
            ylab(expression(paste('Mean ', 'percent ', 'intact ', 'singletons (', s[i], ')', sep = ''))) +
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
                     parse = FALSE) +
            annotate(geom = "text",
                     x = 40,
                     y = 0.4,
                     color = "red",
                     label = lower_annotation,
                     hjust = 0,
                     parse = FALSE)

ggsave(filename = '/home/gdouglas/scripts/pangenome_pseudogene_null/display_items/main_si_vs_sp.pdf',
       plot = si_vs_sp,
       device = 'pdf',
       dpi = 600,
       width = 6,
       height = 5)

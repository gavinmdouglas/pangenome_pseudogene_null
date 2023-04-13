rm(list = ls(all.names = TRUE))

library(ggplot2)

pangenome <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/pangenome_and_related_metrics_filt.tsv.gz",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

upper_annotation <- "Higher ratios consistent with\nmore positive selection to\ngain rare accessory genes"
lower_annotation <- "Lower ratios consistent with\nless positive selection to\ngain rare accessory genes"
  
si_vs_sp <- ggplot(data = pangenome,
                   aes(x = mean_percent_singletons_pseudo_per9,
                       y = mean_percent_singletons_per9)) +
  geom_point() +
  theme_bw() +
  xlab(expression(paste('Mean ', 'percent ', 'pseudogene ', 'singletons (', s[p], ')', sep = ''))) +
  ylab(expression(paste('Mean ', 'percent ', 'intact ', 'singletons (', s[i], ')', sep = ''))) +
  geom_abline(slope = 3/20, intercept = 0, colour = "black", lty = 2, size = 2, alpha = 0.3) +
  geom_segment(aes(x = 40, y = 5, xend = 60, yend = 2),
               arrow = arrow(length = unit(0.5, "cm"), type = 'closed'),
               inherit.aes = FALSE,
               colour = 'red',
               linejoin = 'mitre',
               alpha = 0.005, 
               linewidth = 5) +
  geom_segment(aes(x = 34, y = 6, xend = 14, yend = 9),
               arrow = arrow(length = unit(0.5, "cm"), type = 'closed'),
               inherit.aes = FALSE,
               colour = 'blue',
               linejoin = 'mitre',
               alpha = 0.005, 
               linewidth = 5) +
  annotate(geom = "text",
           x = 0,
           y = 10.25,
           color = "blue",
           label = upper_annotation,
           hjust = 0,
           parse = FALSE) +
    annotate(geom = "text",
             x = 47,
             y = 0.75,
             color = "red",
             label = lower_annotation,
             hjust = 0,
             parse = FALSE)


ggsave(filename = '/home/gdouglas/scripts/pangenome_pseudogene_null/display_items/main_si_vs_sp.pdf',
       plot = si_vs_sp,
       device = 'pdf',
       dpi = 600,
       width = 6,
       height = 6)

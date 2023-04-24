rm(list = ls(all.names = TRUE))

library(ggplot2)

pangenome <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/pangenome_and_related_metrics_filt.tsv.gz",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

# First compute mean genome size per species.
broad_genome_accession <- read.table(file = "/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/genome_info/accessions.tsv.gz",
                                     header = TRUE, sep = "\t", stringsAsFactors = FALSE)
broad_genome_accession <- broad_genome_accession[which(broad_genome_accession$species %in% rownames(pangenome)), ]
broad_genome_accession <- broad_genome_accession[which(broad_genome_accession$could_download), ]
rownames(broad_genome_accession) <- broad_genome_accession$accession

genome_sizes <- read.table('/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/genome_info/genome_sizes.tsv.gz',
                           header = FALSE, sep = '\t', stringsAsFactors = FALSE)
colnames(genome_sizes) <- c('accession', 'length')
genome_sizes <- genome_sizes[which(genome_sizes$accession %in% broad_genome_accession$accession), ]
genome_sizes$species <- broad_genome_accession[genome_sizes$accession, 'species']
genome_size_means <- aggregate(length ~ species, data = genome_sizes, FUN = mean)
rownames(genome_size_means) <- genome_size_means$species
pangenome$genome_size <- genome_size_means$length

upper_annotation <- "Higher ratios consistent with\nmore positive selection to\ngain rare accessory genes"
lower_annotation <- "Lower ratios consistent with\nless positive selection to\ngain rare accessory genes"

si_vs_sp <- ggplot(data = pangenome,
                   aes(x = mean_percent_singletons_pseudo_per9,
                       y = mean_percent_singletons_per9,
                       colour = log10(genome_size))) +
  geom_point() +
  scale_colour_gradient(limits = c(5.75, 7.00)) +
  labs(colour = expression(paste(log[10], '(Genome size)', sep = ''))) +
  theme_bw() +
  xlab(expression(paste('Mean ', 'percent ', 'pseudogene ', 'singletons (', s[p], ')', sep = ''))) +
  ylab(expression(paste('Mean ', 'percent ', 'intact ', 'singletons (', s[i], ')', sep = ''))) +
  geom_abline(slope = 3/20, intercept = 0, colour = "black", lty = 2, linewidth = 2, alpha = 0.3) +
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
       width = 8,
       height = 6)

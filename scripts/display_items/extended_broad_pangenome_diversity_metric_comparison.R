rm(list = ls(all.names = TRUE))

library(ggplot2)
library(cowplot)

# Compare mean number and % singletons vs genome size.
# Also compare with genomic fluidity and show funky species with population substructure.
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

size_v_num_single <- ggplot(data = pangenome,
                            aes(x = genome_size,
                                y = mean_num_singletons_per9)) +
  geom_point() +
  theme_bw() +
  ylab('Mean number of singletons') +
  xlab('Genome size') +
  ggpubr::stat_cor()

size_v_percent_single <- ggplot(data = pangenome,
                                aes(x = genome_size,
                                    y = mean_percent_singletons_per9)) +
  geom_point() +
  theme_bw() +
  ylab('Mean percent of singletons') +
  xlab('Genome size') +
  ggpubr::stat_cor()

percent_single_vs_fluidity <- ggplot(data = pangenome,
                                     aes(x = genomic_fluidity,
                                         y = mean_percent_singletons_per9)) +
  geom_point() +
  theme_bw() +
  ylab('Mean percent of singletons') +
  xlab('Genomic fluidity') +
  ggpubr::stat_cor()


Mycoplasmopsis_bovis_panaroo_filename <- paste('/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/non.focal_genomes/panaroo_non.focal/',
                                               'Mycoplasmopsis_bovis', '/gene_presence_absence.csv.gz', sep = '')

Mycoplasmopsis_bovis_panaroo <- read.table(Mycoplasmopsis_bovis_panaroo_filename,
                                           header = TRUE, sep = ",",  stringsAsFactors = FALSE, quote = "", comment.char = "")

num_isolates_per_gene <- data.frame(rowSums(Mycoplasmopsis_bovis_panaroo[, 4:ncol(Mycoplasmopsis_bovis_panaroo), drop = FALSE] != ''))
colnames(num_isolates_per_gene) <- 'count'
Mycoplasmopsis_bovis_genefreq_dist <- ggplot(data = num_isolates_per_gene,
                                             aes(x = count)) +
                                             geom_histogram(binwidth = 1) +
                                             theme_bw() +
                                             ylab('Number of genes') +
                                             xlab('Number of genomes') +
                                             theme(plot.title = element_text(hjust = 0.5)) +
                                             ggtitle(expression(italic('Mycoplasmopsis bovis')*' gene distribution'))

combined_plot <- plot_grid(size_v_num_single,
                           size_v_percent_single,
                           percent_single_vs_fluidity,
                           Mycoplasmopsis_bovis_genefreq_dist,
                           labels = c('a', 'b', 'c', 'd'))

ggsave(filename = '/home/gdouglas/scripts/pangenome_pseudogene_null/display_items/extended_figure_pangenome_metrics_compare.png',
       plot = combined_plot,
       device = 'png',
       dpi = 400,
       width = 10,
       height = 8)

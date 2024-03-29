rm(list = ls(all.names = TRUE))

library(ggplot2)
library(cowplot)

# Compare mean number and % singletons vs mean number of genes.
# Also compare with genomic fluidity and show funky species with population substructure.
pangenome <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/pangenome_and_related_metrics_filt.tsv.gz",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

num_any_v_num_single <- ggplot(data = pangenome,
                               aes(x = mean_num_genes,
                               y = mean_num_singletons_per9)) +
                        geom_point() +
                        theme_bw() +
                        ylab('Mean number of intact singletons') +
                        xlab('Mean number of genes') +
                        ggpubr::stat_cor(aes(label = paste(r.label, gsub("p", "P", ..p.label..), sep = "~`,`~")),
                                         method = 'spearman', label.y.npc = 'top', cor.coef.name = 'rho',
                                         digits = 3)

num_any_v_percent_single <- ggplot(data = pangenome,
                                aes(x = mean_num_genes,
                                    y = mean_percent_singletons_per9)) +
                              geom_point() +
                              theme_bw() +
                              ylab(expression(paste('Mean ', 'percent ', 'intact ', 'singletons (', s[i], ')', sep = ''))) +
                              xlab('Mean number of genes') +
                              ggpubr::stat_cor(aes(label = paste(r.label, gsub("p", "P", ..p.label..), sep = "~`,`~")),
                                               method = 'spearman', label.y.npc = 'top', cor.coef.name = 'rho',
                                               digits = 3)

percent_single_vs_fluidity <- ggplot(data = pangenome,
                                     aes(x = genomic_fluidity,
                                         y = mean_percent_singletons_per9)) +
                              geom_point() +
                              theme_bw() +
                              ylab(expression(paste('Mean ', 'percent ', 'intact ', 'singletons (', s[i], ')', sep = ''))) +
                              xlab('Genomic fluidity') +
                              ggpubr::stat_cor(aes(label = paste(r.label, gsub("p", "P", ..p.label..), sep = "~`,`~")),
                                               method = 'spearman', label.y.npc = 'top', cor.coef.name = 'rho',
                                               digits = 3)

Mycoplasmopsis_bovis_panaroo <- read.table('/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/example_Mycoplasmopsis_bovis_panaroo_output.csv.gz',
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

combined_plot <- plot_grid(num_any_v_num_single,
                           num_any_v_percent_single,
                           percent_single_vs_fluidity,
                           Mycoplasmopsis_bovis_genefreq_dist,
                           labels = c('a', 'b', 'c', 'd'))

ggsave(filename = '/home/gdouglas/scripts/pangenome_pseudogene_null/display_items/Douglas_ED_Fig6.pdf',
       plot = combined_plot,
       device = 'pdf',
       dpi = 400,
       width = 10,
       height = 8)

# Write out source data:
source_out_a_b_c <- pangenome[, c('mean_num_genes', 'mean_num_singletons_per9', 'mean_percent_singletons_per9', 'genomic_fluidity')]
source_out_a_b_c$species <- rownames(source_out_a_b_c)
source_out_a_b_c <- source_out_a_b_c[, c('species', 'mean_num_genes', 'mean_num_singletons_per9', 'mean_percent_singletons_per9', 'genomic_fluidity')]
write.table(x = source_out_a_b_c,
            file = "/home/gdouglas/scripts/pangenome_pseudogene_null/display_source_data/ED_Fig6_abc.tsv",
            col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)


source_out_d <- data.frame(mycoplasmopsis_bovis_gene = Mycoplasmopsis_bovis_panaroo$Gene,
                           mycoplasmopsis_bovis_gene_annotation = Mycoplasmopsis_bovis_panaroo$Annotation,
                           num_isolates_w_gene = num_isolates_per_gene$count)
write.table(x = source_out_d,
            file = "/home/gdouglas/scripts/pangenome_pseudogene_null/display_source_data/ED_Fig6_d.tsv",
            col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)

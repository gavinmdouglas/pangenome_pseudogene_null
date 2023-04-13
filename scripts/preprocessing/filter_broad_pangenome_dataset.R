rm(list = ls(all.names = TRUE))

# Combine broad dataset information and filter

# Remove species with < 9 genomes and also those with median pairwise dS == 0.
pangenome <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/genomic_fluidity_and_related.tsv",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

Ne_related <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/dnds_summary.tsv",
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

pseudo <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/pseudogene_summary.tsv",
                     header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

if ((! identical(rownames(pseudo), rownames(Ne_related))) || (! identical(rownames(pseudo), rownames(pangenome)))) {
  print("WARNING - different rownames (or order of rownames)!") 
}

pangenome <- cbind(pangenome, Ne_related)
pangenome <- cbind(pangenome, pseudo)

pangenome$mean_percent_singletons_per9 <- (pangenome$mean_num_singletons_per9 / pangenome$mean_num_genes) * 100
pangenome$mean_percent_singletons_pseudo_per9 <- (pangenome$mean_num_singletons_pseudo_per9 / pangenome$mean_num_pseudo) * 100
pangenome$si_sp <- pangenome$mean_percent_singletons_per9 / pangenome$mean_percent_singletons_pseudo_per9

sp_taxonomy <- read.table('/data1/gdouglas/projects/accessory_vs_pseudogene/mapfiles/focal_and_non.focal_taxonomy.tsv',
                          row.names = 1, stringsAsFactors = FALSE, sep = '\t', header = TRUE)
pangenome$class <- sp_taxonomy[rownames(pangenome), "class"]


# Write out this table for sharing on FigShare:
write.table(x = pangenome,
           file = "/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/pangenome_and_related_metrics.tsv",
           col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t")

# And taxonomy table:
write.table(x = sp_taxonomy[rownames(pangenome), ],
           file = "/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/taxonomy.tsv",
           col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t")

# Decided to filter out two species with < 9 genomes (Micromonospora_arenicola and Micromonospora_oceanensis).
# These genomes must have been removed/changed for some reason in between the GTDB release and when I downloaded them.
pangenome <- pangenome[which(pangenome$num_panaroo_genomes >= 9), ]

# I decided to remove species with median dS == 0, as the dN/dS inferences definitely cannot be trusted with such little synonymous-site diversity.
pangenome <- pangenome[which(pangenome$median_pairwise_ds > 0), ]

# Save filtered table:
write.table(x = pangenome,
            file = "/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/pangenome_and_related_metrics_filt.tsv",
            col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t")

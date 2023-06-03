rm(list = ls(all.names = TRUE))

# Combine broad dataset information and filter

# Remove species with < 9 genomes and also those with median pairwise dS == 0.
pangenome <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/genomic_fluidity_and_related.tsv",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

Ne_related <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/combined_core_dnds_summary.tsv",
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

sp_taxonomy <- read.table('/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/taxonomy.tsv.gz',
                          row.names = 1, stringsAsFactors = FALSE, sep = '\t', header = TRUE)
pangenome$class <- sp_taxonomy[rownames(pangenome), "class"]

# Also compute mean genome size per species.
broad_genome_accession <- read.table(file = "/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/genome_info/accessions.tsv.gz",
                                     header = TRUE, sep = "\t", stringsAsFactors = FALSE)
broad_genome_accession <- broad_genome_accession[which(broad_genome_accession$species %in% rownames(pangenome)), ]
broad_genome_accession <- broad_genome_accession[which(broad_genome_accession$could_download), ]
rownames(broad_genome_accession) <- broad_genome_accession$accession
genome_sizes <- read.table('/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/genome_info/genome_sizes.tsv.gz',
                           header = FALSE, sep = '\t', stringsAsFactors = FALSE)
colnames(genome_sizes) <- c('accession', 'length')
genome_sizes <- genome_sizes[which(genome_sizes$accession %in% broad_genome_accession$accession), ]
genome_sizes$species <- broad_genome_accession[genome_sizes$accession, 'species']
genome_size_means <- aggregate(length ~ species, data = genome_sizes, FUN = mean)
rownames(genome_size_means) <- genome_size_means$species
pangenome$genome_size <- genome_size_means$length

# Write out this table for sharing on FigShare:
write.table(x = pangenome,
           file = "/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/pangenome_and_related_metrics.tsv",
           col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t")

# Decided to filter out two species with < 9 genomes (Micromonospora_arenicola and Micromonospora_oceanensis).
# These genomes must have been removed/changed for some reason in between the GTDB release and when I downloaded them.
pangenome <- pangenome[which(pangenome$num_panaroo_genomes >= 9), ]

# Save filtered table:
write.table(x = pangenome,
            file = "/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/pangenome_and_related_metrics_filt.tsv",
            col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t")

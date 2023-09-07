rm(list = ls(all.names = TRUE))

# Quick sanity check that a random species has the correct number of clusters in each frequency range reported in the heatmaps
# (Based on COG-annotated and not).

# Read in cluster membership table.
cluster_breakdown <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/indepth_10_species_analysis/cluster_member_breakdown.tsv.gz",
                                header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cluster_filt_lengths <- read.table(file = "/data1/gdouglas/projects/accessory_vs_pseudogene/summary/cluster_filt_lengths_and_additional.tsv.gz",
                                   header = TRUE, sep = "\t", row.names = 1)

# Restrict to clusters that remained after length filters.
cluster_breakdown <- cluster_breakdown[which(cluster_breakdown$cluster %in% rownames(cluster_filt_lengths)), ]

cluster_breakdown_intact <- cluster_breakdown[grep("pseudo", cluster_breakdown$gene, invert = TRUE), ]
cluster_breakdown_pseudo <- cluster_breakdown[grep("pseudo", cluster_breakdown$gene, invert = FALSE), ]

# Identify clusters in "both" types
# (note that this is based on a comparison of clusters across all 10 species, not just within one species).
clusters_both <- intersect(cluster_breakdown_intact$cluster, cluster_breakdown_pseudo$cluster)

# Remove clusters of "both" types from tables.
cluster_breakdown_intact <- cluster_breakdown_intact[which(! cluster_breakdown_intact$cluster %in% clusters_both), ]
cluster_breakdown_pseudo <- cluster_breakdown_pseudo[which(! cluster_breakdown_pseudo$cluster %in% clusters_both), ]

cluster_breakdown_mixed <- cluster_breakdown[which(cluster_breakdown$cluster %in% clusters_both), ]

# Subample to specific species.
sp_cluster_breakdown_intact <- cluster_breakdown_intact[which(cluster_breakdown_intact$species == 'Staphylococcus_epidermidis'), ]
sp_cluster_breakdown_pseudo <- cluster_breakdown_pseudo[which(cluster_breakdown_pseudo$species == 'Staphylococcus_epidermidis'), ]
sp_cluster_breakdown_mixed <- cluster_breakdown_mixed[which(cluster_breakdown_mixed$species == 'Staphylococcus_epidermidis'), ]

total_genomes <- length(unique(cluster_breakdown[which(cluster_breakdown$species == 'Staphylococcus_epidermidis'), 'accession']))

COG_category_annot <- readRDS("/data1/gdouglas/projects/accessory_vs_pseudogene/summary/all_filt_cluster_COG_annot.rds")
# # Save this file in the Zenodo repository, as I realized it wasn't included:
# COG_category_annot$cluster <- rownames(COG_category_annot)
# COG_category_annot <- COG_category_annot[, c('cluster', 'COG', 'COG_category')]
# write.table(x = COG_category_annot,
#             file = '/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/indepth_10_species_analysis/cluster_COG_annot.tsv',
#             col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')

sp_cluster_breakdown_mixed_annot <- sp_cluster_breakdown_mixed[which(sp_cluster_breakdown_mixed$cluster %in% rownames(COG_category_annot)), ]
sp_cluster_breakdown_intact_annot <- sp_cluster_breakdown_intact[which(sp_cluster_breakdown_intact$cluster %in% rownames(COG_category_annot)), ]
sp_cluster_breakdown_pseudo_annot <- sp_cluster_breakdown_pseudo[which(sp_cluster_breakdown_pseudo$cluster %in% rownames(COG_category_annot)), ]

quick_pangenome_partition_breakdown <- function(breakdown_df, total_genomes) {
  breakdown_df <- breakdown_df[, c('cluster', 'accession')]
  breakdown_df <- breakdown_df[which(! duplicated(breakdown_df)), ]
  cluster_tallies <- table(breakdown_df$cluster)
  ultra_rare <- length(which(cluster_tallies <= 2))
  other_rare <- length(which(cluster_tallies <= (0.15 * total_genomes))) - ultra_rare
  shell <- length(which(cluster_tallies > (0.15 * total_genomes) & cluster_tallies < (0.95 * total_genomes)))
  soft.core <- length(which(cluster_tallies >= 0.95 * total_genomes))
  out <- c(ultra_rare, other_rare, shell, soft.core)
  names(out) <- c('ultra_rare', 'other_rare', 'shell', 'soft.core')
  return(out)
}

# COG-annotated breakdown.
quick_pangenome_partition_breakdown(sp_cluster_breakdown_intact_annot, total_genomes)
quick_pangenome_partition_breakdown(sp_cluster_breakdown_mixed_annot, total_genomes)
quick_pangenome_partition_breakdown(sp_cluster_breakdown_pseudo_annot, total_genomes)

# All clusters.
quick_pangenome_partition_breakdown(sp_cluster_breakdown_intact, total_genomes)
quick_pangenome_partition_breakdown(sp_cluster_breakdown_mixed, total_genomes)
quick_pangenome_partition_breakdown(sp_cluster_breakdown_pseudo, total_genomes)


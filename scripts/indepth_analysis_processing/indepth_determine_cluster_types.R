rm(list = ls(all.names = TRUE))

# Determine whether each cluster is:
  # 1) Solely pseudogenes.
  # 2) Solely intact genes.
  # 3) A mix of pseudogene and intact gene sequences.

# Read in cluster membership table.
cluster_breakdown <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/summary/cluster_member_breakdown.tsv.gz",
                                header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cluster_breakdown_intact <- cluster_breakdown[grep("pseudo", cluster_breakdown$gene, invert = TRUE), ]
cluster_breakdown_pseudo <- cluster_breakdown[grep("pseudo", cluster_breakdown$gene, invert = FALSE), ]

# Identify clusters in "both" types.
clusters_both <- intersect(cluster_breakdown_intact$cluster, cluster_breakdown_pseudo$cluster)

# Remove clusters of "both" types from tables.
cluster_breakdown_intact <- cluster_breakdown_intact[which(! cluster_breakdown_intact$cluster %in% clusters_both), ]
cluster_breakdown_pseudo <- cluster_breakdown_pseudo[which(! cluster_breakdown_pseudo$cluster %in% clusters_both), ]


cluster_types <- list("intact" = unique(cluster_breakdown_intact$cluster),
                      "intergenic.pseudogene" = unique(cluster_breakdown_pseudo$cluster),
                      "both" = unique(clusters_both))

saveRDS(object = cluster_types,
        file = "/data1/gdouglas/projects/accessory_vs_pseudogene/summary/cluster_types.rds")

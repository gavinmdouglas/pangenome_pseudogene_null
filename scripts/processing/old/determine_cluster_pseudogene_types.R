# For clusters that include pseudogenes, return whether they are intergenic, ORF small/large, or fragments.

rm(list = ls(all.names = TRUE))

cluster_member_breakdown <- read.table("/data1/gdouglas/projects/hgt_fragments/pseudogenes/summary/cluster_member_breakdown.tsv.gz",
                                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cluster_pangenome_categories <- readRDS(file = "/data1/gdouglas/projects/hgt_fragments/pseudogenes/summary/cluster_pangenome_categories.rds")

clusters2keep <- as.character()

for (sp in names(cluster_pangenome_categories)) {
  for (pangenome_category in c("cloud", "soft.core", "shell")) {
    clusters2keep <- c(clusters2keep, cluster_pangenome_categories[[sp]][["pseudo"]][[pangenome_category]])
    clusters2keep <- c(clusters2keep, cluster_pangenome_categories[[sp]][["pseudo"]][[pangenome_category]])
  }
}

clusters2keep <- unique(clusters2keep)

cluster_member_breakdown <- cluster_member_breakdown[which(cluster_member_breakdown$cluster %in% clusters2keep), ]
rownames(cluster_member_breakdown) <- cluster_member_breakdown$gene


# Figure out type of pseudogene for all individual genes.
intergenic_ids <- as.character()
outlier.length_ids <- as.character()
fragments_ids <- as.character()

all_id_files <- list.files("/data1/gdouglas/projects/hgt_fragments/pseudogenes/processed/pseudo_filtered_ids", pattern = ".txt$", full.names = TRUE)

for (id_file in all_id_files) {

  ids_in <- read.table(id_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  intergenic_rows <- which(ids_in$V2 == "intergenic")
  if (length(intergenic_rows) > 0) {
    intergenic_ids <- c(intergenic_ids, ids_in[intergenic_rows, "V1"])
  }
  
  outlier.length_rows <- which(ids_in$V2 == "ORF_small" | ids_in$V2 == "ORF_large")
  if (length(outlier.length_rows) > 0) {
    outlier.length_ids <- c(outlier.length_ids, ids_in[outlier.length_rows, "V1"])
  }
  
  fragments_rows <- which(ids_in$V2 == "fragments")
  if (length(fragments_rows) > 0) {
    fragments_ids <- c(fragments_ids, ids_in[fragments_rows, "V1"])
  }
}

intergenic_clusters <- unique(cluster_member_breakdown[intergenic_ids, "cluster"])
outlier.length_clusters <- unique(cluster_member_breakdown[outlier.length_ids, "cluster"])
fragments_clusters <- unique(cluster_member_breakdown[fragments_ids, "cluster"])

# Save RDS containing cross-species cluster information.
saveRDS(object = cross_species_cluster_summary,
        file = "/data1/gdouglas/projects/hgt_fragments/pseudogenes/summary/cross_species_cluster_summary.rds")


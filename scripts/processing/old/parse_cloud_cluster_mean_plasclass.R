rm(list = ls(all.names = TRUE))

# Get mean plasclass of all underlying genes per cluster per species. Focused on *intact* clusters that in the cloud genome per species.

cluster_pangenome_categories <- readRDS("/data1/gdouglas/projects/hgt_fragments/pseudogenes/summary/cluster_pangenome_categories.rds")

all_intact_clusters <- as.character()
for (sp in names(cluster_pangenome_categories)) {
  all_intact_clusters <- c(all_intact_clusters, cluster_pangenome_categories[[sp]]$intact$cloud)
}
all_intact_clusters <- unique(all_intact_clusters)

cluster_info <- read.table("/data1/gdouglas/projects/hgt_fragments/pseudogenes/summary/cluster_member_breakdown_w_contig.tsv.gz",
                           header = TRUE, sep = "\t", stringsAsFactors = FALSE)

contig_id_map <- read.table("/data1/gdouglas/projects/hgt_fragments/mapfiles/contig_id_mapping.tsv.gz",
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(contig_id_map) <- contig_id_map$orig

contig_by_plasclass <- read.table("/data1/gdouglas/projects/hgt_fragments/plasclass/all_plasclass_out.tsv",
                                  header = FALSE, sep = "\t", stringsAsFactors = FALSE)
contig_by_plasclass <- contig_by_plasclass[which(contig_by_plasclass$V1 %in% rownames(contig_id_map)), ]

rownames(contig_by_plasclass) <- gsub("^.*Prokka\\|", "", contig_id_map[contig_by_plasclass$V1, "prokka_renamed"])

length(which(! cluster_info$contig %in% rownames(contig_by_plasclass)))

cluster_info <- cluster_info[which(cluster_info$cluster %in% all_intact_clusters), ]

cluster_info$plasclass <- contig_by_plasclass[cluster_info$contig, "V2"]

cloud_intact_plasclass <- list()
for (sp in names(cluster_pangenome_categories)) {
  cloud_intact_plasclass[[sp]] <- cluster_info[which(cluster_info$cluster %in% cluster_pangenome_categories[[sp]]$intact$cloud & cluster_info$species == sp),
                                               c("cluster", "gene", "contig", "accession", "plasclass")]
}

saveRDS(object = cloud_intact_plasclass,
        file = "/data1/gdouglas/projects/hgt_fragments/pseudogenes/summary/cloud_intact_plasclass.rds")


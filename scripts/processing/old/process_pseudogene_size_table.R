rm(list = ls(all.names = TRUE))

# Save sizes table as RDS after filtering to only intergenic pseudogenes present in filtered files and adding in species, cluster, and accession details.
# Do this so that it's much faster to regenerate the related figures!

library(cowplot)
library(ggplot2)
library(ggbeeswarm)
library(ggtext)

sizes <- read.table("/data1/gdouglas/projects/hgt_fragments/pseudogenes/summary/intergenic_pseudogene_size_and_vs_db.tsv.gz",
                    header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

region_info <- read.table("/data1/gdouglas/projects/hgt_fragments/pseudogenes/summary/cluster_member_breakdown_w_contig.tsv.gz",
                          header = TRUE, sep = "\t", row.names = 2, stringsAsFactors = FALSE)

# Keep only intergenic pseudogenes present in filtered files.
intersecting_rows <- rownames(sizes)[which(rownames(sizes) %in% rownames(region_info))]

# Subset to these intersecting rows.
sizes <- sizes[intersecting_rows, ]
region_info <- region_info[intersecting_rows, ]

sizes$species <- region_info$species
sizes$cluster <- region_info$cluster
sizes$accession <- region_info$accession

saveRDS(object = sizes, file = "/data1/gdouglas/projects/hgt_fragments/pseudogenes/summary/intergenic_pseudogene_size_and_vs_db_filt_w_additional.rds")
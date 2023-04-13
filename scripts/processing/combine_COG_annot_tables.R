rm(list = ls(all.names = TRUE))

## Read in COG annot information
intact_annot <- readRDS("/data1/gdouglas/projects/accessory_vs_pseudogene/summary/intact_and_both_cluster_annot.rds")
intact_annot <- intact_annot[-which(intact_annot$COG_category == "-"), ]
intact_annot <- intact_annot[-which(intact_annot$COG_category == ""), ]

## For intact clusters, filter out clusters outside size range (at least based on what clusters they're in).
cluster_filt_lengths <- read.table(file = "/data1/gdouglas/projects/accessory_vs_pseudogene/summary/cluster_filt_lengths_and_additional.tsv.gz",
                                   header = TRUE, sep = "\t", row.names = 1)

intact_annot <- intact_annot[which(rownames(intact_annot) %in% rownames(cluster_filt_lengths)), ]

pseudo_annot <- readRDS("/data1/gdouglas/projects/accessory_vs_pseudogene/summary/intergenic_pseudo_cluster_COG_annot.rds")
pseudo_annot <- pseudo_annot$pseudo
pseudo_annot <- pseudo_annot[-which(pseudo_annot$majority_COG_category == ""), ]

rownames(pseudo_annot) <- pseudo_annot$cluster
pseudo_annot <- pseudo_annot[, c(2, 3)]
colnames(pseudo_annot) <- c("COG", "COG_category")

intact_annot <- intact_annot[, c("all_COG", "COG_category")]
colnames(intact_annot) <- c("COG", "COG_category")

COG_category_annot <- rbind(intact_annot, pseudo_annot)

saveRDS(object = COG_category_annot,
        file = "/data1/gdouglas/projects/accessory_vs_pseudogene/summary/all_filt_cluster_COG_annot.rds")

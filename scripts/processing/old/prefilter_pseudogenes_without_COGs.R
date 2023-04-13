# Filter down to intergenic pseudogenes only and separate those that can be annotated with a COG id or not.

rm(list = ls(all.names = TRUE))

COG_category_annot <- readRDS("/data1/gdouglas/projects/hgt_fragments/pseudogenes/processed/cdhit_workflow/eggnog_mapper_annot/cluster_annot.rds")
colnames(COG_category_annot)[which(colnames(COG_category_annot) == "COG_category_all.COGs")] <- "COG_category"
COG_category_annot <- COG_category_annot[which(COG_category_annot$COG_category != "-"), ]

cluster_pseudo_types <- read.table("/data1/gdouglas/projects/hgt_fragments/pseudogenes/summary/cluster_pseudo_types.tsv.gz",
                                   header = TRUE, sep = "\t", row.names = 1)

intergenic_pseudo_clusters <- rownames(cluster_pseudo_types)[which(cluster_pseudo_types$intergenic > 0)]

COG_intergenic_pseudo_clusters <- intergenic_pseudo_clusters[which(intergenic_pseudo_clusters %in%  rownames(COG_category_annot))]
nonCOG_intergenic_pseudo_clusters <- intergenic_pseudo_clusters[which(! intergenic_pseudo_clusters %in%  rownames(COG_category_annot))]

saveRDS(object = COG_intergenic_pseudo_clusters,
        file = "/data1/gdouglas/projects/hgt_fragments/pseudogenes/summary/intergenic_pseudo_clusters_COG_annot.rds")

saveRDS(object = nonCOG_intergenic_pseudo_clusters,
        file = "/data1/gdouglas/projects/hgt_fragments/pseudogenes/summary/intergenic_pseudo_clusters_nonCOG_annot.rds")

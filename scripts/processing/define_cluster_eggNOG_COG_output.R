rm(list = ls(all.names = TRUE))

### Parse eggNOG annotation and define COG categories based on the presence of *any* COG amongst list of OGs that were hit.

COG2catgory_uniq <- readRDS("/data1/gdouglas/db/COG_definitions/cog-20.to_category_collapse.rds")

# Read in cluster membership table + definitions of cluster types.
cluster_breakdown <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/summary/cluster_member_breakdown.tsv.gz",
                                header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 2)

cluster_types <- readRDS(file = "/data1/gdouglas/projects/accessory_vs_pseudogene/summary/cluster_types.rds")


# Read in COG annotations for intergenic elements based on UniRef hit annotations.
intergenic_cluster_annot <- readRDS("/data1/gdouglas/projects/accessory_vs_pseudogene/summary/intergenic_pseudo_cluster_COG_annot.rds")

# Also read in full eggNOG annotations based on cluster-level as well.
intact_cluster_annot <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/eggnog_all.clusters/intact_seq_annot.emapper.annotations.tsv.gz",
                                   header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, quote = "", comment.char = "")

# Get COG info for clusters based on COG definition of any orthogroups in EggNOG hierarchy definition.
intact_cluster_annot$all_COG <- sapply(intact_cluster_annot$eggNOG_OGs, function(x) { paste(grep("^COG", unique(gsub("@.*$", "", strsplit(x, ",")[[1]])), value = TRUE), collapse = ",") })

intact_cluster_annot$cluster <- cluster_breakdown[rownames(intact_cluster_annot), "cluster"]

intact_cluster_annot$COG_category <- "-"

intact_i_with_COG <- which(intact_cluster_annot$all_COG != "")

intact_cluster_annot$COG_category[intact_i_with_COG] <- sapply(intact_i_with_COG,
                                                               function(i) {
                                                                 
                                                                 all_COGs <- intact_cluster_annot[i, "all_COG"]
                                                                 categories <- as.character()
                                                                 COGs <- strsplit(all_COGs, ",")[[1]]
                                                                 for (COG in COGs) {
                                                                   categories <- c(categories, COG2catgory_uniq[COG, "category"])
                                                                 }
                                                                 
                                                                 return(paste(sort(unique(categories)), collapse = ","))
                                                               })

# Compare COG category annotations for clusters in "both".
intergenic_based_both_annot <- intergenic_cluster_annot$both_intact_and_pseudo_uniref_approach
intact_based_both_annot <- intact_cluster_annot[which(intact_cluster_annot$cluster %in% intergenic_based_both_annot$cluster), ]

rownames(intergenic_based_both_annot) <- intergenic_based_both_annot$cluster
intergenic_based_both_annot_subset <- intergenic_based_both_annot[intact_based_both_annot$cluster, ]

intact_based_both_annot_called <- intact_based_both_annot[which(intact_based_both_annot$COG_category != "-"), ]
intergenic_based_both_annot_subset_called <- intergenic_based_both_annot_subset[which(intergenic_based_both_annot_subset$majority_COG_category != ""), ]
length(intersect(intact_based_both_annot_called$cluster, intergenic_based_both_annot_subset_called$cluster))
# Long-story short: the good news is that the annotations are actually really similar across both methods!
# The vast majority of the COG category annotations are exactly the same (of the clusters that can be annotated by either).


intact_cluster_annot <- intact_cluster_annot[, c("cluster", "Description", "all_COG", "COG_category")]

rownames(intact_cluster_annot) <- cluster_breakdown[rownames(intact_cluster_annot), "cluster"]

saveRDS(object = intact_cluster_annot,
        file = "/data1/gdouglas/projects/accessory_vs_pseudogene/summary/intact_and_both_cluster_annot.rds")

rm(list = ls(all.names = TRUE))

# Report mean and sd % non-COG annotated
# Also % in cloud partitions vs core
# Also compute Fisher's exact tests and report summary

all.clusters_percent_by_type_fill <- read.table(file = "/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/indepth_10_species_analysis/cluster_breakdown_tables/all.clusters_percent_by_type_fill.tsv.gz",
                                                      header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

cog.category.annot_percent_by_type_fill <- read.table(file = "/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/indepth_10_species_analysis/cluster_breakdown_tables/cog.category.annot_percent_by_type_fill.tsv.gz",
                                                       header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

# Get counts only:

get_counts_only <- function(in_df) {
  in_df <- apply(in_df, 2, FUN = function(x) { gsub('^.*\\(', '', x) } )
  in_df <- apply(in_df, 2, FUN = function(x) { gsub('\\)', '', x) } )
  in_df <- apply(in_df, 2, FUN = function(x) { as.integer(x) })
  return(data.frame(in_df))
}

all.clusters_counts <- get_counts_only(all.clusters_percent_by_type_fill)
cog.category.annot_counts <- get_counts_only(cog.category.annot_percent_by_type_fill)

cog_annot_intact = cog.category.annot_counts$intact_ultra.cloud + cog.category.annot_counts$intact_other.cloud + cog.category.annot_counts$intact_shell + cog.category.annot_counts$intact_soft.core
cog_annot_pseudo = cog.category.annot_counts$pseudo_ultra.cloud + cog.category.annot_counts$pseudo_other.cloud + cog.category.annot_counts$pseudo_shell + cog.category.annot_counts$pseudo_soft.core


total_intact = all.clusters_counts$intact_ultra.cloud + all.clusters_counts$intact_other.cloud + all.clusters_counts$intact_shell + all.clusters_counts$intact_soft.core
total_pseudo = all.clusters_counts$pseudo_ultra.cloud + all.clusters_counts$pseudo_other.cloud + all.clusters_counts$pseudo_shell + all.clusters_counts$pseudo_soft.core

# Mean and sd percent annot by element type:
round(mean((cog_annot_intact / total_intact) * 100), 2)
round(sd((cog_annot_intact / total_intact) * 100), 2)

round(mean((cog_annot_pseudo / total_pseudo) * 100), 2)
round(sd((cog_annot_pseudo / total_pseudo) * 100), 2)

ORs <- as.numeric()
Ps <- as.numeric()
for (i in 1:10) {
  contigency_table <- matrix(c(cog_annot_intact[i], total_intact[i], cog_annot_pseudo[i], total_pseudo[i]), nrow = 2, ncol = 2)
  fisher_out <- fisher.test(contigency_table)
  ORs <- c(ORs, fisher_out$estimate)
  Ps <- c(Ps, fisher_out$p.value)
}

round(mean(ORs), 2)
round(sd(ORs), 2)
length(which(Ps < 0.05))

rownames(all.clusters_percent_by_type_fill)[which((Ps < 0.05) & (ORs > 1))]
rownames(all.clusters_percent_by_type_fill)[which((Ps < 0.05) & (ORs < 1))]

# Percents by pangenome partition
all.clusters_percent_by_type <- read.table(file = "/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/indepth_10_species_analysis/cluster_breakdown_tables/all.clusters_percent_by_type.tsv.gz",
                                                header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

round(mean(all.clusters_percent_by_type$intact_ultra.cloud + all.clusters_percent_by_type$intact_other.cloud), 2)
round(sd(all.clusters_percent_by_type$intact_ultra.cloud + all.clusters_percent_by_type$intact_other.cloud), 2)

round(mean(all.clusters_percent_by_type$pseudo_ultra.cloud + all.clusters_percent_by_type$pseudo_other.cloud), 2)
round(sd(all.clusters_percent_by_type$pseudo_ultra.cloud + all.clusters_percent_by_type$pseudo_other.cloud), 2)

round(mean(all.clusters_percent_by_type$pseudo_soft.core), 2)
round(sd(all.clusters_percent_by_type$pseudo_soft.core), 2)

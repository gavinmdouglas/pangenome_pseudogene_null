rm(list = ls(all.names = TRUE))

# Get ranges 

pangenome <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/pangenome_and_related_metrics_filt.tsv.gz",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

# Intact genes
round(range(pangenome$genomic_fluidity), 3)
round(range(pangenome$mean_num_genes), 1)
round(range(pangenome$mean_num_singletons_per9), 2)
round(range(pangenome$mean_percent_singletons_per9), 2)

# Pseudogenes
round(range(pangenome$pseudogene_genomic_fluidity), 3)
round(range(pangenome$mean_num_pseudo), 1)
round(range(pangenome$mean_num_singletons_pseudo_per9), 2)
round(range(pangenome$mean_percent_singletons_pseudo_per9), 2)

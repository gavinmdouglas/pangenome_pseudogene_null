rm(list = ls(all.names = TRUE))

# Check data distributions for different taxonomic classes across main dataset.

pangenome <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/pangenome_and_related_metrics_filt.tsv.gz",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

pangenome_Bacilli <- pangenome[which(pangenome$class == 'Bacilli'), ]

boxplot(pangenome$ds)
boxplot(log10(pangenome$dnds))
boxplot(pangenome$mean_num_genes)
boxplot(pangenome$si_sp)

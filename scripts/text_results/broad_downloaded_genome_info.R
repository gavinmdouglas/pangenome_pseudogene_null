rm(list = ls(all.names = TRUE))

# Info on genomes for broad analysis.
pangenome <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/pangenome_and_related_metrics_filt.tsv.gz",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

# How many genomes total:
sum(pangenome$num_panaroo_genomes)


# Taxonomic class breakdown:
sp_taxonomy <- read.table('//data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/taxonomy.tsv.gz',
                          row.names = 1, stringsAsFactors = FALSE, sep = '\t', header = TRUE)
sp_taxonomy <- sp_taxonomy[rownames(pangenome), ]

table(sp_taxonomy$class)

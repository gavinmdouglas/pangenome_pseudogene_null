rm(list = ls(all.names = TRUE))

pangenome <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/pangenome_and_related_metrics_filt.tsv.gz",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

# Check three highlighted taxa values of si and sp.
pangenome['Escherichia_coli', c('mean_percent_singletons_per9', 'mean_percent_singletons_pseudo_per9')]
pangenome['Chlamydia_trachomatis', c('mean_percent_singletons_per9', 'mean_percent_singletons_pseudo_per9')]
pangenome['Rickettsia_prowazekii', c('mean_percent_singletons_per9', 'mean_percent_singletons_pseudo_per9')]

cor.test(pangenome$mean_percent_singletons_per9, pangenome$mean_percent_singletons_pseudo_per9, method = 'spearman')

# Shown in figure panels
cor.test(pangenome$dnds, pangenome$ds, method = 'spearman')
cor.test(pangenome$mean_percent_singletons_per9, pangenome$ds, method = 'spearman')
cor.test(pangenome$mean_percent_singletons_pseudo_per9, pangenome$ds, method = 'spearman')
cor.test(pangenome$mean_num_genes, pangenome$dnds, method = 'spearman')
cor.test(pangenome$genomic_fluidity, pangenome$dnds, method = 'spearman')
cor.test(pangenome$mean_percent_singletons_per9, pangenome$dnds, method = 'spearman')
cor.test(pangenome$si_sp, pangenome$dnds, method = 'spearman')

# dS vs si/sp
cor.test(pangenome$si_sp, pangenome$ds, method = 'spearman')

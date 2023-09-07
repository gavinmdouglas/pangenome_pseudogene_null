rm(list = ls(all.names = TRUE))

pangenome <- read.table('/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/pangenome_and_related_metrics_filt.tsv.gz',
                        header = TRUE, sep = '\t', row.names = 1)

metrics <- c('ds',
             'dnds',
             'mean_num_genes',
             'genomic_fluidity',
             'mean_percent_singletons_per9',
             'mean_percent_singletons_pseudo_per9',
             'si_sp')

pangenome <- pangenome[, metrics]

cor(pangenome, method = 'spearman')

pangenome_filt <- pangenome[which(pangenome$ds > 0.01 & pangenome$dnds < 0.5), ]
cor(pangenome_filt, method = 'spearman')

pangenome_metrics <- c('mean_num_genes',
                       'genomic_fluidity',
                       'mean_percent_singletons_per9',
                       'mean_percent_singletons_pseudo_per9',
                       'si_sp')
sapply(pangenome_metrics, function(x) { ppcor::pcor.test(x = pangenome$dnds, y = pangenome[, x], z = pangenome$ds, method = 'spearman')$estimate})

sapply(pangenome_metrics, function(x) { ppcor::pcor.test(x = pangenome_filt$dnds, y = pangenome_filt[, x], z = pangenome_filt$ds, method = 'spearman')$estimate})

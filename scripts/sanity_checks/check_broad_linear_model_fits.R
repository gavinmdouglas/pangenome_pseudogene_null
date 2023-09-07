rm(list = ls(all.names = TRUE))

out_models <- readRDS('/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/model_output/pangenome_linear_models.rds')

summary(out_models$mean_num_genes)

summary(out_models$genomic_fluidity)

summary(out_models$mean_percent_singletons_per9)

summary(out_models$si_sp)


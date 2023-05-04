rm(list = ls(all.names = TRUE))

out_models <- readRDS(file = '/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/model_output/pangenome_rank_based_linear_models.rds')

si_model_summary <- summary(out_models$si)
round(si_model_summary$R2, 3)

si_sp_model_summary <- summary(out_models$si_sp)
round(si_sp_model_summary$R2, 3)

round(si_model_summary$coefficients['dnds', 'Estimate'], 3)
round(si_model_summary$coefficients['dnds', 'p.value'], 4)

round(si_sp_model_summary$coefficients['dnds', 'Estimate'], 3)
round(si_sp_model_summary$coefficients['dnds', 'p.value'], 4)

# Number of significant classes:
length(which(si_model_summary$coefficients[grep("class", rownames(si_model_summary$coefficients)), 'p.value'] < 0.05))
length(which(si_sp_model_summary$coefficients[grep("class", rownames(si_sp_model_summary$coefficients)), 'p.value'] < 0.05))

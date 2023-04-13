rm(list = ls(all.names = TRUE))

out_models <- readRDS(file = '/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/model_output/pangenome_linear_models.rds')

si_model_summary <- summary(out_models$si$class_dS_dNdS)
round(si_model_summary$adj.r.squared, 2)

si_sp_model_summary <- summary(out_models$si_sp$class_dS_dNdS)
round(si_sp_model_summary$adj.r.squared, 2)


round(si_model_summary$coefficients['log_c_dS', 'Estimate'], 2)
round(si_model_summary$coefficients['log_c_dS', 'Pr(>|t|)'], 2)

round(si_model_summary$coefficients['log_c_dNdS', 'Estimate'], 2)
round(si_model_summary$coefficients['log_c_dNdS', 'Pr(>|t|)'], 2)


round(si_sp_model_summary$coefficients['log_c_dS', 'Estimate'], 2)
round(si_sp_model_summary$coefficients['log_c_dS', 'Pr(>|t|)'], 2)

round(si_sp_model_summary$coefficients['log_c_dNdS', 'Estimate'], 2)
round(si_sp_model_summary$coefficients['log_c_dNdS', 'Pr(>|t|)'], 2)

# Number of significant classes:
length(which(si_model_summary$coefficients[grep("class", rownames(si_model_summary$coefficients)), 'Pr(>|t|)'] < 0.05))
length(which(si_sp_model_summary$coefficients[grep("class", rownames(si_sp_model_summary$coefficients)), 'Pr(>|t|)'] < 0.05))

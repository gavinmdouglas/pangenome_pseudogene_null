rm(list = ls(all.names = TRUE))

ultra.cloud_model_w_additional <- readRDS("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/indepth_10_species_analysis/glmm_output/ultra.cloud_model_outputs_additional.rds")
other.cloud_model_w_additional <- readRDS("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/indepth_10_species_analysis/glmm_output/other.cloud_model_outputs_additional.rds")
shell_model_w_additional <- readRDS("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/indepth_10_species_analysis/glmm_output/shell_model_outputs_additional.rds")

ultra.cloud_model <- readRDS("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/indepth_10_species_analysis/glmm_output/ultra.cloud_model_outputs.rds")
other.cloud_model <- readRDS("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/indepth_10_species_analysis/glmm_output/other.cloud_model_outputs.rds")
shell_model <- readRDS("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/indepth_10_species_analysis/glmm_output/shell_model_outputs.rds")

# Number of elements per model
dim(ultra.cloud_model_w_additional$COG.only_species.interaction_redundant.both.interaction$frame)
dim(other.cloud_model_w_additional$COG.only_species.interaction_redundant.both.interaction$frame)
dim(shell_model_w_additional$COG.only_species.interaction_redundant.both.interaction$frame)


# AIC breakdown per input.
sapply(ultra.cloud_model, AIC)
sapply(ultra.cloud_model_w_additional, AIC)

sapply(other.cloud_model_w_additional, AIC)
sapply(other.cloud_model, AIC)

sapply(shell_model_w_additional, AIC)
sapply(shell_model, AIC)


# Species random effect coefficients per input.
glmmTMB::ranef(ultra.cloud_model$COG.only_species.interaction_redundant.both.interaction)$cond$species
glmmTMB::ranef(other.cloud_model$COG.only_species.interaction_redundant.both.interaction)$cond$species
glmmTMB::ranef(shell_model$COG.only_species.interaction_redundant.both.interaction)$cond$species


# Check model coefficients.
ultra.cloud_model_coef <- summary(ultra.cloud_model$COG.only_species.interaction_redundant.both.interaction)
other.cloud_model_coef <- summary(other.cloud_model$COG.only_species.interaction_redundant.both.interaction)
shell_model_model_coef <- summary(shell_model$COG.only_species.interaction_redundant.both.interaction)


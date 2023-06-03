rm(list = ls(all.names = TRUE))

# Parse out species random effects.
ultra.cloud_model <- readRDS("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/indepth_10_species_analysis/glmm_output/ultra.cloud_model_outputs_additional.rds")
random_effects_ultra.cloud_species <- ranef(ultra.cloud_model$COG.only_species.interaction_redundant.both.interaction)$cond$species
random_effects_ultra.cloud_species$species <- rownames(random_effects_ultra.cloud_species)
random_effects_ultra.cloud_species$partition <- "Ultra-cloud"

other.cloud_model <- readRDS("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/indepth_10_species_analysis/glmm_output/other.cloud_model_outputs_additional.rds")
random_effects_other.cloud_species <- ranef(other.cloud_model$COG.only_species.interaction_redundant.both.interaction)$cond$species
random_effects_other.cloud_species$species <- rownames(random_effects_other.cloud_species)
random_effects_other.cloud_species$partition <- "Other-cloud"

shell_model <- readRDS("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/indepth_10_species_analysis/glmm_output/shell_model_outputs_additional.rds")
random_effects_shell_species <- ranef(shell_model$COG.only_species.interaction_redundant.both.interaction)$cond$species
random_effects_shell_species$species <- rownames(random_effects_shell_species)
random_effects_shell_species$partition <- "Shell"

# Make combined table and write out.
random_effects_species <- rbind(random_effects_ultra.cloud_species, random_effects_other.cloud_species)
random_effects_species <- rbind(random_effects_species, random_effects_shell_species)

colnames(random_effects_species)[which(colnames(random_effects_species) == "(Intercept)")] <- "Intercept"

write.table(x = random_effects_species,
            file = "/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/indepth_10_species_analysis/glmm_output/species_random_effects.tsv",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

rm(list = ls(all.names = TRUE))

library(ggplot2)

random_effects_species <- read.table(file = "/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/indepth_10_species_analysis/glmm_output/species_random_effects.tsv.gz",
                                     header = TRUE, sep = "\t", stringsAsFactors = FALSE)

random_effects_species[which(random_effects_species$partition == "Ultra-cloud"), "partition"] <- "Ultra-rare"
random_effects_species[which(random_effects_species$partition == "Other-cloud"), "partition"] <- "Other-rare"

random_effects_species$partition <- factor(random_effects_species$partition, levels = c("Ultra-rare", "Other-rare", "Shell"))

colnames(random_effects_species)[which(colnames(random_effects_species) == "(Intercept)")] <- "Intercept"

random_effects_species$species <- gsub('_', ' ', random_effects_species$species)
random_effects_species$species <- factor(random_effects_species$species, levels = rev(sort(unique(random_effects_species$species))))

species_random_effects <- ggplot(data = random_effects_species, aes(x = Intercept, y = species)) +
                                  geom_bar(stat="identity") +
                                  theme_bw() +
                                  facet_wrap(. ~ partition) +
                                  ylab("Random effect") +
                                  xlab("Intercept") +
                                  geom_vline(xintercept = 0, linetype="dotted", 
                                             color = "black") +
                                  theme(axis.text.y = element_text(face = "italic"))

ggsave(filename = '/home/gdouglas/scripts/pangenome_pseudogene_null/display_items/extended_indepth_glmm_species_random_effects.pdf',
       plot = species_random_effects,
       device = 'pdf',
       dpi = 400,
       width = 8,
       height = 6)

### Explore frequency and species distributions of pseudogenes vs intact genes
### Also breakdown things like how often they are in the same cluster.
### Save files with key info for later.

rm(list = ls(all.names = TRUE))

library(ggplot2)
library(ggVennDiagram)
library(kableExtra)
library(knitr)
library(cowplot)

# Somewhat arbitrarily use < 15% as cut-off for "cloud" so that frequencies can be compared more sensibly.
cloud_cutoff <- list()
soft.core_cutoff <- list()

pseudo_filt_genomes <- gsub(".txt$", "", list.files("/data1/gdouglas/projects/hgt_fragments/pseudogenes/processed/pseudo_filtered_ids/"))

for (ACC_FILE in list.files(path = "/data1/gdouglas/projects/hgt_fragments/genome_accessions/", full.names = FALSE)) {
  
  sp <- gsub("_accessions.txt", "", ACC_FILE)
  
  sp_genomes <- read.table(paste("/data1/gdouglas/projects/hgt_fragments/genome_accessions/", ACC_FILE, sep = ""),
                           stringsAsFactors = FALSE)$V1
  
  num_filt_genomes <- length(which(sp_genomes %in% pseudo_filt_genomes))
  
  cloud_cutoff[[sp]] <- num_filt_genomes * 0.15
  
  soft.core_cutoff[[sp]] <- num_filt_genomes * 0.95
}

cluster_member_breakdown <- read.table("/data1/gdouglas/projects/hgt_fragments/pseudogenes/summary/cluster_member_breakdown.tsv.gz",
                                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)

pseudogene_clusters <- unique(cluster_member_breakdown$cluster[grep("_pseudo", cluster_member_breakdown$gene)])
intact_clusters <- unique(cluster_member_breakdown$cluster[grep("_pseudo", cluster_member_breakdown$gene, invert = TRUE)])

intersecting_clusters <- pseudogene_clusters[which(pseudogene_clusters %in% intact_clusters)]

# Breakdown of total tallies.
ggVennDiagram(list(Pseudogenes=pseudogene_clusters, Intact=intact_clusters)) +
  ggtitle("All species") +
  theme(plot.title = element_text(hjust = 0.5))



all_species <- sort(unique(cluster_member_breakdown$species))

sp_venndiagrams <- list()

for (sp in all_species) {

  sp_cluster_member_breakdown <- cluster_member_breakdown[which(cluster_member_breakdown$species == sp), ]

  sp_pseudogene_clusters <- unique(sp_cluster_member_breakdown$cluster[grep("_pseudo", sp_cluster_member_breakdown$gene)])
  sp_intact_clusters <- unique(sp_cluster_member_breakdown$cluster[grep("_pseudo", sp_cluster_member_breakdown$gene, invert = TRUE)])
  
  sp_venndiagrams[[sp]] <- ggVennDiagram(list(Pseudogenes=sp_pseudogene_clusters, Intact=sp_intact_clusters)) +
                                        ggtitle(sp) +
                                        theme(plot.title = element_text(hjust = 0.5))
}

plot_grid(sp_venndiagrams$Agrobacterium_tumefaciens,
          sp_venndiagrams$Enterococcus_faecalis,
          sp_venndiagrams$Escherichia_coli,
          sp_venndiagrams$Lactococcus_lactis,
          sp_venndiagrams$Pseudomonas_aeruginosa,
          sp_venndiagrams$Sinorhizobium_meliloti,
          sp_venndiagrams$Staphylococcus_epidermidis,
          sp_venndiagrams$Streptococcus_pneumoniae,
          sp_venndiagrams$Wolbachia_pipientis,
          sp_venndiagrams$Xanthomonas_oryzae)



pangenome_category_breakdown <- data.frame(matrix(NA, nrow = length(all_species), ncol = 9))
rownames(pangenome_category_breakdown) <- all_species
colnames(pangenome_category_breakdown) <- c("pseudo_cloud", "pseudo_shell", "pseudo_soft.core",
                                            "intact_cloud", "intact_shell", "intact_soft.core",
                                            "both_cloud", "both_shell", "both_soft.core")

cluster_pangenome_categories <- list()

for (sp in all_species) {
  
  cluster_pangenome_categories[[sp]] <- list("pseudo" = list(),
                                             "intact" = list(),
                                             "both" = list())
  
  # Restrict to hits in specific species only and get table with just cluster and accession ids.
  sp_cluster_member_breakdown <- cluster_member_breakdown[which(cluster_member_breakdown$species == sp), ]
  sp_cluster_member_breakdown_simple <- sp_cluster_member_breakdown[, c("cluster", "accession")]
  
  # Remove lines representing multiple instances in same genome.
  sp_cluster_member_breakdown_simple <- sp_cluster_member_breakdown_simple[-which(duplicated(sp_cluster_member_breakdown_simple)), ]
  
  cluster_frequency <- table(sp_cluster_member_breakdown_simple$cluster)
  
  cluster_frequency_both <- cluster_frequency[which(names(cluster_frequency) %in% intersecting_clusters)]
  
  cluster_frequency_nointersect <- cluster_frequency[which(! names(cluster_frequency) %in% intersecting_clusters)]
  cluster_frequency_nointersect_pseudo <- cluster_frequency_nointersect[which(names(cluster_frequency_nointersect) %in% pseudogene_clusters)]
  cluster_frequency_nointersect_intact <- cluster_frequency_nointersect[which(names(cluster_frequency_nointersect) %in% intact_clusters)]
  
  cluster_pangenome_categories[[sp]][["pseudo"]][["cloud"]] <- names(cluster_frequency_nointersect_pseudo)[which(cluster_frequency_nointersect_pseudo <= cloud_cutoff[[sp]])]
  cluster_pangenome_categories[[sp]][["pseudo"]][["soft.core"]] <- names(cluster_frequency_nointersect_pseudo)[which(cluster_frequency_nointersect_pseudo >= soft.core_cutoff[[sp]])]
  cluster_pangenome_categories[[sp]][["pseudo"]][["shell"]] <- names(cluster_frequency_nointersect_pseudo)[which(! names(cluster_frequency_nointersect_pseudo) %in% c(cluster_pangenome_categories[[sp]][["pseudo"]][["cloud"]], cluster_pangenome_categories[[sp]][["pseudo"]][["soft.core"]]))]
  
  cluster_pangenome_categories[[sp]][["intact"]][["cloud"]] <- names(cluster_frequency_nointersect_intact)[which(cluster_frequency_nointersect_intact <= cloud_cutoff[[sp]])]
  cluster_pangenome_categories[[sp]][["intact"]][["soft.core"]] <- names(cluster_frequency_nointersect_intact)[which(cluster_frequency_nointersect_intact >= soft.core_cutoff[[sp]])]
  cluster_pangenome_categories[[sp]][["intact"]][["shell"]] <- names(cluster_frequency_nointersect_intact)[which(! names(cluster_frequency_nointersect_intact) %in% c(cluster_pangenome_categories[[sp]][["intact"]][["cloud"]], cluster_pangenome_categories[[sp]][["intact"]][["soft.core"]]))]
  
  cluster_pangenome_categories[[sp]][["both"]][["cloud"]] <- names(cluster_frequency_both)[which(cluster_frequency_both <= cloud_cutoff[[sp]])]
  cluster_pangenome_categories[[sp]][["both"]][["soft.core"]] <- names(cluster_frequency_both)[which(cluster_frequency_both >= soft.core_cutoff[[sp]])]
  cluster_pangenome_categories[[sp]][["both"]][["shell"]] <- names(cluster_frequency_both)[which(! names(cluster_frequency_both) %in% c(cluster_pangenome_categories[[sp]][["both"]][["cloud"]], cluster_pangenome_categories[[sp]][["both"]][["soft.core"]]))]
  
  pangenome_category_breakdown[sp, "pseudo_cloud"] <- (length(cluster_pangenome_categories[[sp]][["pseudo"]][["cloud"]]) / length(cluster_frequency_nointersect_pseudo)) * 100
  pangenome_category_breakdown[sp, "pseudo_soft.core"] <- (length(cluster_pangenome_categories[[sp]][["pseudo"]][["soft.core"]]) / length(cluster_frequency_nointersect_pseudo)) * 100
  pangenome_category_breakdown[sp, "pseudo_shell"] <- 100 - pangenome_category_breakdown[sp, "pseudo_cloud"] - pangenome_category_breakdown[sp, "pseudo_soft.core"]
  
  pangenome_category_breakdown[sp, "intact_cloud"] <- (length(cluster_pangenome_categories[[sp]][["intact"]][["cloud"]]) / length(cluster_frequency_nointersect_intact)) * 100
  pangenome_category_breakdown[sp, "intact_soft.core"] <- (length(cluster_pangenome_categories[[sp]][["intact"]][["soft.core"]]) / length(cluster_frequency_nointersect_intact)) * 100
  pangenome_category_breakdown[sp, "intact_shell"] <- 100 - pangenome_category_breakdown[sp, "intact_cloud"] - pangenome_category_breakdown[sp, "intact_soft.core"]
  
  pangenome_category_breakdown[sp, "both_cloud"] <- (length(cluster_pangenome_categories[[sp]][["both"]][["cloud"]]) / length(cluster_frequency_both)) * 100
  pangenome_category_breakdown[sp, "both_soft.core"] <- (length(cluster_pangenome_categories[[sp]][["both"]][["soft.core"]]) / length(cluster_frequency_both)) * 100
  pangenome_category_breakdown[sp, "both_shell"] <- 100 - pangenome_category_breakdown[sp, "both_cloud"] - pangenome_category_breakdown[sp, "both_soft.core"]
  
}

write.table(x = pangenome_category_breakdown,
            file = "/data1/gdouglas/projects/hgt_fragments/pseudogenes/summary/category_pangenome_percents.tsv",
            col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)

pangenome_category_breakdown_rounded <- round(pangenome_category_breakdown, 2)

# Summary table of pangenome category %
kbl(pangenome_category_breakdown_rounded) %>%
  kable_styling()


# Save RDS containing breakdowns of which genes are in which categories.
saveRDS(object = cluster_pangenome_categories,
        file = "/data1/gdouglas/projects/hgt_fragments/pseudogenes/summary/cluster_pangenome_categories.rds")


# Identify cross-species clusters and summarize what species they are in and what category subtype they are (e.g., "both_cloud").
cluster_member_breakdown_by_species <- cluster_member_breakdown[, c("species", "cluster")]

cluster_member_breakdown_single_by_species <- cluster_member_breakdown[-which(duplicated(cluster_member_breakdown_by_species)), ]

cluster_member_breakdown_single_by_species_tally <- table(cluster_member_breakdown_single_by_species$cluster)

cross_species_clusters <- names(cluster_member_breakdown_single_by_species_tally)[which(cluster_member_breakdown_single_by_species_tally > 1)]

cluster_member_breakdown_single_by_species_cross.species <- cluster_member_breakdown_single_by_species[which(cluster_member_breakdown_single_by_species$cluster %in% cross_species_clusters), ]


cross_species_cluster_summary <- list()

for (cross_species_cluster in cross_species_clusters) {
  
  cross_species_cluster_subset <- cluster_member_breakdown_single_by_species_cross.species[which(cluster_member_breakdown_single_by_species_cross.species$cluster == cross_species_cluster), ]
  
  category_vec <- as.character()
  
  for (i in 1:nrow(cross_species_cluster_subset)) {
   
    shared_sp <- cross_species_cluster_subset[i, "species"]

    gene_type <- c("intact", "pseudo", "both")
    pangenome_type <- c("cloud", "soft.core", "shell")
    
    for (g in gene_type) {
      
      for (p in pangenome_type) {
     
        if (cross_species_cluster %in% cluster_pangenome_categories[[shared_sp]][[g]][[p]]) {
         
          category_vec <- c(category_vec, paste(g, p, sep = "|"))
          break
  
        }
      }
    }
  }

  names(category_vec) <- cross_species_cluster_subset$species
  
  cross_species_cluster_summary[[cross_species_cluster]] <- category_vec

}

# Save RDS containing cross-species cluster information.
saveRDS(object = cross_species_cluster_summary,
        file = "/data1/gdouglas/projects/hgt_fragments/pseudogenes/summary/cross_species_cluster_summary.rds")


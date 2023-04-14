rm(list = ls(all.names = TRUE))

source("~/scripts/accessory_vs_pseudogenes/R_scripts/functions.R")

# Summarize pseudogene clusters per species.

focal_species <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/mapfiles/focal_species.txt",
                            stringsAsFactors = FALSE, header = FALSE)$V1

nonfocal_species <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/mapfiles/non.focal_species.txt",
                               stringsAsFactors = FALSE, header = FALSE)$V1

focal_and_nonfocal_species <- sort(c(focal_species, nonfocal_species))

pseudo_summary <- data.frame(matrix(NA, nrow = length(focal_and_nonfocal_species), ncol = 6))
rownames(pseudo_summary) <- focal_and_nonfocal_species
colnames(pseudo_summary) <- c("mean_num_pseudo", "mean_num_singleton_pseudo",
                              "mean_num_singletons_pseudo_per9","sd_num_singletons_pseudo_per9",
                              "num_pseudofinder_genomes", "pseudogene_genomic_fluidity")

for (sp in focal_and_nonfocal_species) {
    
    infile <- NULL
    if (sp %in% focal_species) {
      infile <- paste('/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/focal_genomes_subsampled/pseudogenes/pseudogene_cluster_presence_tables/',
                      sp, '.tsv', sep = '')
    } else if (sp %in% nonfocal_species) {
      infile <- paste('/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/non.focal_genomes/pseudogenes/pseudogene_cluster_presence_tables/',
                      sp, '.tsv', sep = '')
    }
  
   if (! file.exists(infile)) {
    stop('No infile for this species: ', sp)
   }
   
   intab <- read.table(infile, header = TRUE, sep = "\t",  stringsAsFactors = FALSE, row.names = 1)
   
   num_isolates_per_pseudo <- rowSums(intab != '')
   intab_singletons <- intab[which(num_isolates_per_pseudo == 1), ]
   mean_num_singletons <- mean(colSums(intab_singletons != ""))
   
   mean_num_pseudo <- mean(colSums(intab != ''))
   
   mean_num_singletons_per9 <- compute_mean_num_singletons_per_combo(gene_presence = intab, k = 9)
   
   pseudogene_genomic_fluidity <- genomic_fluidity_from_df(intab)
   
   pseudo_summary[sp, ] <- c(mean_num_pseudo, mean_num_singletons,
                             mean_num_singletons_per9$mean, mean_num_singletons_per9$sd,
                             ncol(intab), pseudogene_genomic_fluidity)
  
}

write.table(x = pseudo_summary,
           file = "/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/pseudogene_summary.tsv",
           col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)

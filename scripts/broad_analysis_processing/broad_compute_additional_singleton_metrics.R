rm(list = ls(all.names = TRUE))

# Compute si and sp for 3 and 20 genomes as well.

source("/home/gdouglas/scripts/pangenome_pseudogene_null/scripts/functions.R")

focal_species <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/indepth_10_species_analysis/focal_species.txt",
                            stringsAsFactors = FALSE, header = FALSE)$V1

focal_and_nonfocal_species <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/focal_and_non.focal_full_to_short.tsv.gz",
                                          stringsAsFactors = FALSE, header = FALSE, sep = '\t')$V1

pangenome_summary <- data.frame(matrix(NA, nrow = length(focal_and_nonfocal_species), ncol = 12))
rownames(pangenome_summary) <- focal_and_nonfocal_species
colnames(pangenome_summary) <- c("mean_num_genes", "mean_num_pseudo",
                                 "mean_num_singletons_per3","sd_num_singletons_per3",
                                 "mean_num_singletons_per20","sd_num_singletons_per20",
                                 "mean_num_singletons_pseudo_per3","sd_num_singletons_pseudo_per3",
                                 "mean_num_singletons_pseudo_per20","sd_num_singletons_pseudo_per20",
                                 "num_panaroo_genomes", "num_pseudofinder_genomes")

set.seed(131)
for (sp in focal_and_nonfocal_species) {
    
    infile <- NULL
    pseudo_infile <- NULL
    if (sp %in% focal_species) {
      infile <- paste('/data1/gdouglas/projects/accessory_vs_pseudogene/focal_panaroo/',
                      sp, '/gene_presence_absence_ran.subset.csv', sep = '')
      pseudo_infile <- paste('/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/focal_genomes_subsampled/pseudogenes/pseudogene_cluster_presence_tables/',
                             sp, '.tsv', sep = '')
    } else if (! sp %in% focal_species) {
      infile <- paste('/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/non.focal_genomes/panaroo_non.focal/',
                      sp, '/gene_presence_absence.csv', sep = '')
      pseudo_infile <- paste('/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/non.focal_genomes/pseudogenes/pseudogene_cluster_presence_tables/',
                             sp, '.tsv', sep = '')
    }
  
   if (! file.exists(infile)) {
    stop('No infile for this species: ', sp)
   }
   
  if (! file.exists(pseudo_infile)) {
      stop('No pseudo_infile for this species: ', sp)
    }
    
   intab <- read.table(infile, header = TRUE, sep = ",",  stringsAsFactors = FALSE, quote = "", comment.char = "")
   
   accessions <- colnames(intab[, 4:ncol(intab), drop = FALSE])
   
   mean_num_genes <- mean(colSums(intab[, accessions, drop = FALSE] != ''))
   
   mean_num_singletons_per3 <- compute_mean_num_singletons_per_combo(gene_presence = intab[, accessions, drop = FALSE], k = 3)
   mean_num_singletons_per20 <- compute_mean_num_singletons_per_combo(gene_presence = intab[, accessions, drop = FALSE], k = 20)

   
   pseudo_intab <- read.table(pseudo_infile, header = TRUE, sep = "\t",  stringsAsFactors = FALSE, row.names = 1)

   mean_num_pseudo <- mean(colSums(pseudo_intab != ''))
   
   mean_num_pseudo_singletons_per3 <- compute_mean_num_singletons_per_combo(gene_presence = pseudo_intab, k = 3)
   mean_num_pseudo_singletons_per20 <- compute_mean_num_singletons_per_combo(gene_presence = pseudo_intab, k = 20)
   
   pangenome_summary[sp, ] <- c(mean_num_genes, mean_num_pseudo,
                                mean_num_singletons_per3$mean, mean_num_singletons_per3$sd,
                                mean_num_singletons_per20$mean, mean_num_singletons_per20$sd,
                                mean_num_pseudo_singletons_per3$mean, mean_num_pseudo_singletons_per3$sd,
                                mean_num_pseudo_singletons_per20$mean, mean_num_pseudo_singletons_per20$sd,
                                length(accessions), ncol(pseudo_intab))
}

write.table(x = pangenome_summary,
           file = "/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/genomic_fluidity_and_related_additional_subsamples.tsv",
           col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)

rm(list = ls(all.names = TRUE))

focal_species <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/mapfiles/focal_species.txt",
                            stringsAsFactors = FALSE, header = FALSE)$V1

nonfocal_species <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/mapfiles/non.focal_species.txt",
                               stringsAsFactors = FALSE, header = FALSE)$V1

focal_and_nonfocal_species <- sort(c(focal_species, nonfocal_species))

output <- data.frame(matrix(NA, nrow = length(focal_and_nonfocal_species), ncol = 4))
rownames(output) <- focal_and_nonfocal_species
colnames(output) <- c("median_pairwise_dn", "median_pairwise_ds", "median_pairwise_dnds", "median_MG94_dnds")

for (sp in focal_and_nonfocal_species) {
    
    dnds_pairwise_infile <- NULL
    if (sp %in% focal_species) {
      dnds_pairwise_infile <- paste('/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/focal_genomes_subsampled/core_gene_cds/pairwise_dnds_stats/',
                      sp, '.tsv', sep = '')
    } else if (sp %in% nonfocal_species) {
      dnds_pairwise_infile <- paste('/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/non.focal_genomes/panaroo_non.focal_cds/pairwise_dnds_stats/',
                      sp, '.tsv.gz', sep = '')
    }
   if (! file.exists(dnds_pairwise_infile)) {
    stop('No dnds_pairwise_infile for this species: ', sp)
   }
   pairwise <- read.table(dnds_pairwise_infile, header = FALSE, sep = "\t",  stringsAsFactors = FALSE, row.names = 1)
   colnames(pairwise) <- c('mean_n_subs', 'mean_n_sites', 'mean_s_subs', 'mean_s_sites', 'mean_dn', 'mean_ds', 'mean_dnds')
   
   dnds_MG94_infile <- NULL
   if (sp %in% focal_species) {
     dnds_MG94_infile <- paste('/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/focal_genomes_subsampled/core_gene_cds/tree_dnds/',
                                   sp, '.txt', sep = '')
   } else if (sp %in% nonfocal_species) {
     dnds_MG94_infile <- paste('/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/non.focal_genomes/panaroo_non.focal_cds/tree_dnds/',
                                   sp, '.txt', sep = '')
   }
   if (! file.exists(dnds_MG94_infile)) {
     stop('No dnds_MG94_infile for this species: ', sp)
   }
   MG94 <- read.table(dnds_MG94_infile, header = FALSE, sep = " ",  stringsAsFactors = FALSE, row.names = 2)
  
  
   output[sp, "median_pairwise_dn"] <- median(pairwise$mean_dn, na.rm = TRUE)
   output[sp, "median_pairwise_ds"] <- median(pairwise$mean_ds, na.rm = TRUE)
   output[sp, "median_pairwise_dnds"] <- median(pairwise$mean_dnds, na.rm = TRUE)
   output[sp, "median_MG94_dnds"] <- median(MG94$V3, na.rm = TRUE)

}

write.table(x = output,
           file = "/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/dnds_summary.tsv",
           col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)

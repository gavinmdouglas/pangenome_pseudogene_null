rm(list = ls(all.names = TRUE))

source("~/scripts/accessory_vs_pseudogenes/R_scripts/functions.R")

focal_species <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/mapfiles/focal_species.txt",
                            stringsAsFactors = FALSE, header = FALSE)$V1

nonfocal_species <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/mapfiles/non.focal_species.txt",
                               stringsAsFactors = FALSE, header = FALSE)$V1

focal_and_nonfocal_species <- sort(c(focal_species, nonfocal_species))

pangenome_summary <- data.frame(matrix(NA, nrow = length(focal_and_nonfocal_species), ncol = 6))
rownames(pangenome_summary) <- focal_and_nonfocal_species
colnames(pangenome_summary) <- c("mean_num_genes", "mean_num_singletons",
                                 "mean_num_singletons_per9","sd_num_singletons_per9",
                                 "genomic_fluidity", "num_panaroo_genomes")

set.seed(131)
for (sp in focal_and_nonfocal_species) {
    
    infile <- NULL
    if (sp %in% focal_species) {
      infile <- paste('/data1/gdouglas/projects/accessory_vs_pseudogene/focal_panaroo/',
                      sp, '/gene_presence_absence_ran.subset.csv', sep = '')
    } else if (sp %in% nonfocal_species) {
      infile <- paste('/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/non.focal_genomes/panaroo_non.focal/',
                      sp, '/gene_presence_absence.csv.gz', sep = '')
    }
  
   if (! file.exists(infile)) {
    stop('No infile for this species: ', sp)
   }
   
   intab <- read.table(infile, header = TRUE, sep = ",",  stringsAsFactors = FALSE, quote = "", comment.char = "")
   
   accessions <- colnames(intab[, 4:ncol(intab), drop = FALSE])
   
   num_isolates_per_gene <- rowSums(intab[, accessions, drop = FALSE] != '')
   
   intab_singletons <- intab[which(num_isolates_per_gene == 1), accessions, drop = FALSE]
   
   mean_num_genes <- mean(colSums(intab[, accessions, drop = FALSE] != ''))
   
   mean_num_singletons <- mean(colSums(intab_singletons != ""))
   
   mean_num_singletons_per9 <- compute_mean_num_singletons_per_combo(gene_presence = intab[, accessions, drop = FALSE], k = 9)
   
   genomic_fluidity <- genomic_fluidity_from_df(intab[, accessions, drop = FALSE])
   
   pangenome_summary[sp, ] <- c(mean_num_genes, mean_num_singletons,
                                mean_num_singletons_per9$mean, mean_num_singletons_per9$sd, 
                                genomic_fluidity, length(accessions))
}

write.table(x = pangenome_summary,
           file = "/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/genomic_fluidity_and_related.tsv",
           col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)

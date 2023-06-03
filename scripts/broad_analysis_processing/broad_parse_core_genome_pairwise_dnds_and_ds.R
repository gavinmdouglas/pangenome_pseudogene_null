rm(list = ls(all.names = TRUE))

focal_species <- read.table('/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/indepth_10_species_analysis/focal_species.txt',
                            stringsAsFactors = FALSE, header = FALSE)$V1

focal_and_nonfocal_species <- rownames(read.table('/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/taxonomy.tsv.gz',
                                                  stringsAsFactors = FALSE, header = TRUE, sep = '\t', row.names = 1))

focal_and_nonfocal_species <- sort(focal_and_nonfocal_species)

full_core_dnds <- data.frame(matrix(NA, nrow = length(focal_and_nonfocal_species), ncol = 3))
rownames(full_core_dnds) <- focal_and_nonfocal_species
colnames(full_core_dnds) <- c('dn', 'ds', 'dnds')

for (sp in focal_and_nonfocal_species) {
    
    dnds_pairwise_infile <- NULL
    if (sp %in% focal_species) {
      dnds_pairwise_infile <- paste('/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/focal_genomes_subsampled/core_gene_cds/core_genome_msa_dnds/',
                                    sp, '.fna.pairwise.tsv', sep = '')
    } else {
      dnds_pairwise_infile <- paste('/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/non.focal_genomes/panaroo_non.focal_cds/core_genome_msa_dnds/',
                                    sp, '.fna.pairwise.tsv', sep = '')
    }
   if (! file.exists(dnds_pairwise_infile)) {
    stop('No dnds_pairwise_infile for this species: ', sp)
   }
    
  raw_tab <- read.table(dnds_pairwise_infile, header = TRUE, sep = "\t",  stringsAsFactors = FALSE)
  raw_tab <- raw_tab[, c('dn', 'ds', 'dnds')]
  raw_tab <- raw_tab[which(rowSums(is.na(raw_tab)) == 0), ]
  full_core_dnds[sp, c('dn', 'ds', 'dnds')] <- as.numeric(colMeans(raw_tab))
   
}

write.table(x = full_core_dnds,
            file = "/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/combined_core_dnds_summary.tsv",
            col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)

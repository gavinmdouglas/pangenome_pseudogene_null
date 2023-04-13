rm(list = ls(all.names = TRUE))

# Summaries of numbers, sizes, etc. of called intact genes vs pseudogenes.

long_to_short <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/focal_and_non.focal_full_to_short.tsv.gz",
                            header = FALSE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)
colnames(long_to_short) <- 'short'

# First read in information based on the focal species subsampled to 20 genomes each.
focal_species <- read.table('/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/indepth_10_species_analysis/focal_species.txt',
                            header = FALSE, sep = "", stringsAsFactors = FALSE)$V1
focal_pseudogene_lengths_raw <- list()
focal_pseudogene_counts_raw <- list()
for (focal_sp in focal_species) {
  
  # Pseudogene numbers per accession.
  focal_pseudo_link_filename <- paste('/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/focal_genomes_subsampled/pseudogenes/filt_intergenic_pseudogene_ids/',
                                      focal_sp, '.txt', sep = '')
  pseudogene_to_accession_map <- read.table(file = focal_pseudo_link_filename, header = FALSE, sep = '\t', stringsAsFactors = FALSE)
  rownames(pseudogene_to_accession_map) <- paste(long_to_short[focal_sp, 'short'], pseudogene_to_accession_map$V1, sep = '_')
  pseudogenes_per_accession <- table(pseudogene_to_accession_map$V2)
  
  focal_pseudogene_counts_raw[[focal_sp]] <- data.frame(accession = names(pseudogenes_per_accession), pseudogenes = as.integer(pseudogenes_per_accession))
  
  # Pseudogene lengths.  
  focal_pseudo_length_filename <- paste('/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/focal_genomes_subsampled/pseudogenes/filt_intergenic_pseudogene_fastas_lengths/',
                                        focal_sp, '.fasta.tsv', sep = '')
  focal_pseudogene_lengths_raw[[focal_sp]] <- read.table(file = focal_pseudo_length_filename, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  
  focal_pseudogene_lengths_raw[[focal_sp]]$species <- focal_sp
  focal_pseudogene_lengths_raw[[focal_sp]]$accession <- pseudogene_to_accession_map[focal_pseudogene_lengths_raw[[focal_sp]]$sequence_id, 'V2']
}

focal_gene_sizes_raw <- list()
for (focal_acc_gene_sizes in list.files('/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/focal_genomes_subsampled/focal_genomes_subsampled_gene_sizes',
                                        full.names = TRUE,
                                        pattern = '.genesizes.txt')) {
   focal_acc <- gsub('.genesizes.txt$', '', basename(focal_acc_gene_sizes))
   
   focal_gene_sizes_raw[[focal_acc]] <- read.table(file = focal_acc_gene_sizes, header = FALSE, sep = ' ', stringsAsFactors = FALSE)
   
   focal_gene_sizes_raw[[focal_acc]]$accession <- focal_acc
}

focal_gene_sizes <- do.call(rbind, focal_gene_sizes_raw)
focal_pseudogene_counts <- do.call(rbind, focal_pseudogene_counts_raw)
focal_pseudogene_sizes <- do.call(rbind, focal_pseudogene_lengths_raw)

focal_gene_numbers <- read.table(file = '/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/focal_genomes_subsampled/gene_numbers.txt',
                                 sep = ' ', stringsAsFactors = FALSE, header = FALSE)
# Rename to match non-focal table.
focal_gene_numbers <- focal_gene_numbers[, -1]
colnames(focal_gene_numbers) <- c('V1', 'V2')

# Then do the same for nonfocal species.
nonfocal_species <- gsub('.pseudo.fasta.tsv', '', list.files('/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/non.focal_genomes/pseudogenes/non.focal_filt_pseudogenes_lengths'))
nonfocal_pseudogene_lengths_raw <- list()
nonfocal_pseudogene_counts_raw <- list()
for (nonfocal_sp in nonfocal_species) {
  
  # Pseudogene numbers per accession.
  nonfocal_pseudo_link_filename <- paste('/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/non.focal_genomes/pseudogenes/filt_intergenic_pseudo_ids/',
                                         nonfocal_sp, '.txt', sep = '')
  pseudogene_to_accession_map <- read.table(file = nonfocal_pseudo_link_filename, header = FALSE, sep = '\t', stringsAsFactors = FALSE)
  rownames(pseudogene_to_accession_map) <- paste(long_to_short[nonfocal_sp, 'short'], pseudogene_to_accession_map$V1, sep = '_')
  pseudogenes_per_accession <- table(pseudogene_to_accession_map$V2)
  
  nonfocal_pseudogene_counts_raw[[nonfocal_sp]] <- data.frame(accession = names(pseudogenes_per_accession), pseudogenes = as.integer(pseudogenes_per_accession))
  
  # Pseudogene lengths.  
  nonfocal_pseudo_length_filename <- paste('/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/non.focal_genomes/pseudogenes/non.focal_filt_pseudogenes_lengths/',
                                           nonfocal_sp, '.pseudo.fasta.tsv', sep = '')
  nonfocal_pseudogene_lengths_raw[[nonfocal_sp]] <- read.table(file = nonfocal_pseudo_length_filename, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  
  nonfocal_pseudogene_lengths_raw[[nonfocal_sp]]$species <- nonfocal_sp
  nonfocal_pseudogene_lengths_raw[[nonfocal_sp]]$accession <- pseudogene_to_accession_map[nonfocal_pseudogene_lengths_raw[[nonfocal_sp]]$sequence_id, 'V2']

}

nonfocal_gene_sizes_raw <- list()
for (nonfocal_acc_gene_sizes in list.files('/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/non.focal_genomes/non.focal_gff_gene_sizes',
                                        full.names = TRUE,
                                        pattern = '.genesizes.txt')) {
  nonfocal_acc <- gsub('.genesizes.txt$', '', basename(nonfocal_acc_gene_sizes))
  
  nonfocal_gene_sizes_raw[[nonfocal_acc]] <- read.table(file = nonfocal_acc_gene_sizes, header = FALSE, sep = ' ', stringsAsFactors = FALSE)
  
  nonfocal_gene_sizes_raw[[nonfocal_acc]]$accession <- nonfocal_acc
  
}

nonfocal_gene_sizes <- do.call(rbind, nonfocal_gene_sizes_raw)
nonfocal_pseudogene_counts <- do.call(rbind, nonfocal_pseudogene_counts_raw)
nonfocal_pseudogene_sizes <- do.call(rbind, nonfocal_pseudogene_lengths_raw)

nonfocal_gene_numbers <- read.table(file = '/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/non.focal_genomes/gene_numbers.txt',
                                 sep = ' ', stringsAsFactors = FALSE, header = FALSE)

accession_gene_numbers <- rbind(nonfocal_gene_numbers, focal_gene_numbers)
colnames(accession_gene_numbers) <- c("accession", "number_of_genes")

accession_pseudogene_numbers <- rbind(nonfocal_pseudogene_counts, focal_pseudogene_counts)
colnames(accession_pseudogene_numbers) <- c("accession", "number_of_pseudogenes")
rownames(accession_pseudogene_numbers) <- NULL

accession_numbers <- accession_pseudogene_numbers
rownames(accession_numbers) <- accession_numbers$accession

rownames(accession_gene_numbers) <- accession_gene_numbers$accession
accession_numbers$number_of_genes <- accession_gene_numbers[rownames(accession_numbers), "number_of_genes"]
rownames(accession_numbers) <- NULL

write.table(x = accession_numbers,
            file = '/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/element_info/element_counts.tsv',
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)


pseudogene_sizes <- rbind(nonfocal_pseudogene_sizes, focal_pseudogene_sizes)
rownames(pseudogene_sizes) <- NULL
pseudogene_sizes <- pseudogene_sizes[, c('species', 'accession', 'sequence_id', 'length')]

gene_sizes <- rbind(nonfocal_gene_sizes, focal_gene_sizes)
rownames(gene_sizes) <- NULL
gene_sizes <- gene_sizes[, c('accession', 'V1', 'V2')]
colnames(gene_sizes) <- c('accession', 'gene', 'length')

write.table(x = gene_sizes,
            file = '/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/element_info/gene_sizes.tsv',
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(x = pseudogene_sizes,
            file = '/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/element_info/pseudogene_sizes.tsv',
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)

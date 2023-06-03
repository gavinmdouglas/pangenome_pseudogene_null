rm(list = ls(all.names = TRUE))

# Build tables with gene and pseudogene % coverage per accession and averaged across accessions per species.

# First, need to get the genome accession ids of interest (i.e., of filtered species).
pangenome <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/pangenome_and_related_metrics_filt.tsv.gz",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

sp_taxonomy <- read.table('/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/taxonomy.tsv.gz',
                          row.names = 1, stringsAsFactors = FALSE, sep = '\t', header = TRUE)
sp_taxonomy <- sp_taxonomy[rownames(pangenome), ]

broad_genome_accession <- read.table(file = "/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/genome_info/accessions.tsv.gz",
                                     header = TRUE, sep = "\t", stringsAsFactors = FALSE)
broad_genome_accession <- broad_genome_accession[which(broad_genome_accession$species %in% rownames(pangenome)), ]
broad_genome_accession <- broad_genome_accession[which(broad_genome_accession$could_download), ]

# Get % of pseudogene elements across genomes.
element_counts <- read.table('/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/element_info/element_counts.tsv.gz',
                             header = TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)
element_counts <- element_counts[broad_genome_accession$accession, ]
element_counts$percent_pseudo <- (element_counts$number_of_pseudogenes / (element_counts$number_of_pseudogenes + element_counts$number_of_genes)) * 100

gene_sizes <- read.table('/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/element_info/gene_sizes.tsv.gz',
                         header = TRUE, sep = '\t', stringsAsFactors = FALSE)

pseudogene_sizes <- read.table('/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/element_info/pseudogene_sizes.tsv.gz',
                               header = TRUE, sep = '\t', stringsAsFactors = FALSE)

genome_sizes <- read.table('/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/genome_info/genome_sizes.tsv.gz',
                           header = FALSE, sep = '\t', stringsAsFactors = FALSE)


gene_sizes <- gene_sizes[which(gene_sizes$accession %in% broad_genome_accession$accession), ]
pseudogene_sizes <- pseudogene_sizes[which(pseudogene_sizes$accession %in% broad_genome_accession$accession), ]
genome_sizes <- genome_sizes[which(genome_sizes$V1 %in% broad_genome_accession$accession), ]
rownames(genome_sizes) <- genome_sizes$V1

gene_sizes$percent <- (gene_sizes$length / genome_sizes[gene_sizes$accession, 'V2']) * 100
pseudogene_sizes$percent <- (pseudogene_sizes$length / genome_sizes[pseudogene_sizes$accession, 'V2']) * 100

gene_percent_summed <- aggregate(percent ~ accession, data = gene_sizes, FUN = sum)
pseudogene_percent_summed <- aggregate(percent ~ accession, data = pseudogene_sizes, FUN = sum)

write.table(x = gene_percent_summed,
            file = '/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/element_info/element_percent_coverage/gene_percent_coverage_by_accession.tsv',
            col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)

write.table(x = pseudogene_percent_summed,
            file = '/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/element_info/element_percent_coverage/pseudogene_percent_coverage_by_accession.tsv',
            col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)

pseudogene_percent_mean_by_species <- pseudogene_percent_summed
rownames(broad_genome_accession) <- broad_genome_accession$accession
pseudogene_percent_mean_by_species$species <- broad_genome_accession[pseudogene_percent_mean_by_species$accession, 'species']
pseudogene_percent_sd_by_species <- aggregate(percent ~ species, data = pseudogene_percent_mean_by_species, FUN = sd)
pseudogene_percent_mean_by_species <- aggregate(percent ~ species, data = pseudogene_percent_mean_by_species, FUN = mean)

identical(pseudogene_percent_sd_by_species$species, pseudogene_percent_mean_by_species$species)
pseudogene_percent_mean_by_species$sd <- pseudogene_percent_sd_by_species$percent
colnames(pseudogene_percent_mean_by_species) <- c('species', 'mean_percent', 'sd_percent')

write.table(x = pseudogene_percent_mean_by_species,
            file = '/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/element_info/element_percent_coverage/pseudogene_mean_percent_coverage_by_species.tsv',
            col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)

gene_percent_mean_by_species <- gene_percent_summed
rownames(broad_genome_accession) <- broad_genome_accession$accession
gene_percent_mean_by_species$species <- broad_genome_accession[gene_percent_mean_by_species$accession, 'species']
gene_percent_sd_by_species <- aggregate(percent ~ species, data = gene_percent_mean_by_species, FUN = sd)
gene_percent_mean_by_species <- aggregate(percent ~ species, data = gene_percent_mean_by_species, FUN = mean)

identical(gene_percent_sd_by_species$species, gene_percent_mean_by_species$species)
gene_percent_mean_by_species$sd <- gene_percent_sd_by_species$percent
colnames(gene_percent_mean_by_species) <- c('species', 'mean_percent', 'sd_percent')

write.table(x = gene_percent_mean_by_species,
            file = '/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/element_info/element_percent_coverage/gene_mean_percent_coverage_by_species.tsv',
            col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)


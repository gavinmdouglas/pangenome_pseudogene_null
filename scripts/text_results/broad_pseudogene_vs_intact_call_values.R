rm(list = ls(all.names = TRUE))

# Summaries of numbers, sizes, etc. of called intact genes vs pseudogenes.

# First, need to get the genome accession ids of interest (i.e., of filtered species).
pangenome <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/pangenome_and_related_metrics_filt.tsv.gz",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

sp_taxonomy <- read.table('/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/taxonomy.tsv.gz',
                          row.names = 1, stringsAsFactors = FALSE, sep = '\t', header = TRUE)
sp_taxonomy <- sp_taxonomy[rownames(pangenome), ]

broad_genome_accession <- read.table(file = "/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/genome_info/accessions.tsv.gz",
                                     header = TRUE, sep = "\t", stringsAsFactors = FALSE)
broad_genome_accession <- broad_genome_accession[which(broad_genome_accession$species %in% rownames(pangenome)), ]
broad_genome_accession <- broad_genome_accession[which(broad_genome_accession$could_download), ]

# Get % of pseudogene elements across genomes.
element_counts <- read.table('/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/element_info/element_counts.tsv.gz',
                             header = TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)
element_counts <- element_counts[broad_genome_accession$accession, ]
element_counts$percent_pseudo <- (element_counts$number_of_pseudogenes / (element_counts$number_of_pseudogenes + element_counts$number_of_genes)) * 100
round(mean(element_counts$percent_pseudo), 2)
round(sd(element_counts$percent_pseudo), 2)
round(min(element_counts$percent_pseudo), 2)
round(max(element_counts$percent_pseudo), 2)

gene_sizes <- read.table('/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/element_info/gene_sizes.tsv.gz',
                         header = TRUE, sep = '\t', stringsAsFactors = FALSE)

pseudogene_sizes <- read.table('/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/element_info/pseudogene_sizes.tsv.gz',
                               header = TRUE, sep = '\t', stringsAsFactors = FALSE)

genome_sizes <- read.table('/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/genome_info/genome_sizes.tsv.gz',
                           header = FALSE, sep = '\t', stringsAsFactors = FALSE)


gene_sizes <- gene_sizes[which(gene_sizes$accession %in% broad_genome_accession$accession), ]
pseudogene_sizes <- pseudogene_sizes[which(pseudogene_sizes$accession %in% broad_genome_accession$accession), ]
genome_sizes <- genome_sizes[which(genome_sizes$V1 %in% broad_genome_accession$accession), ]
rownames(genome_sizes) <- genome_sizes$V1

gene_sizes$percent <- (gene_sizes$length / genome_sizes[gene_sizes$accession, 'V2']) * 100
pseudogene_sizes$percent <- (pseudogene_sizes$length / genome_sizes[pseudogene_sizes$accession, 'V2']) * 100

gene_percent_summed <- aggregate(percent ~ accession, data = gene_sizes, FUN = sum)
pseudogene_percent_summed <- aggregate(percent ~ accession, data = pseudogene_sizes, FUN = sum)

round(mean(gene_percent_summed$percent), 2)
round(sd(gene_percent_summed$percent), 2)

round(mean(pseudogene_percent_summed$percent), 2)
round(sd(pseudogene_percent_summed$percent), 2)

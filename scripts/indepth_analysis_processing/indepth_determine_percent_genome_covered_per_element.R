rm(list = ls(all.names = TRUE))

### Get these values computed per cluster:
###  - Mean % genome covered by intact clusters
###  - Mean % genome covered by pseudogene clusters  

# Read in cluster membership table + definitions of cluster types.
cluster_breakdown <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/summary/cluster_member_breakdown.tsv.gz",
                                header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 2)

element_sizes <- read.table(file = "/data1/gdouglas/projects/accessory_vs_pseudogene/cdhit_workflow_intact_and_intergenic/filt_intergenic_pseudogenes_and_intact.lengths.tsv.gz",
                            header = TRUE, sep = "\t")

element_sizes$cluster <- cluster_breakdown[element_sizes$sequence_id, "cluster"]

# Filter out elements outside size range (at least based on what clusters they're in).
cluster_filt_lengths <- read.table(file = "/data1/gdouglas/projects/accessory_vs_pseudogene/summary/cluster_filt_lengths_and_additional.tsv.gz",
                                   header = TRUE, sep = "\t", row.names = 1)

element_sizes_filt <- element_sizes[which(element_sizes$cluster %in% rownames(cluster_filt_lengths)), ]

# Subset cluster_breakdown to this same set, and add in length column.
cluster_breakdown <- cluster_breakdown[element_sizes_filt$sequence_id, ]
cluster_breakdown$length <- element_sizes_filt$length


# Read in genome sizes.
species <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/mapfiles/species.txt",
                      header = FALSE, stringsAsFactors = FALSE)$V1
raw_sizes <- list()
for (sp in species) {
  size_file <- paste("/data1/gdouglas/projects/accessory_vs_pseudogene/genome_sizes/", sp, "_sizes.txt", sep = "")
  raw_sizes[[sp]] <- read.table(file = size_file, header = FALSE, sep = " ", stringsAsFactors = FALSE)
}

genome_sizes <- do.call(rbind, raw_sizes)
rownames(genome_sizes) <- genome_sizes$V1
genome_sizes <- genome_sizes[, -1, drop = FALSE]
colnames(genome_sizes) <- "size"

cluster_breakdown$genome_percent <- (cluster_breakdown$length / genome_sizes[cluster_breakdown$accession, "size"]) * 100


cluster_types <- readRDS("/data1/gdouglas/projects/accessory_vs_pseudogene/summary/cluster_types.rds")

cluster_breakdown_intact <- cluster_breakdown[which(cluster_breakdown$cluster %in% cluster_types$intact), ]
cluster_breakdown_pseudo <- cluster_breakdown[which(cluster_breakdown$cluster %in% cluster_types$intergenic.pseudogene), ]

per_genome_coverages <- list()
per_genome_coverages[["intact"]] <- list()
per_genome_coverages[["intergenic.pseudogene"]] <- list()

for (sp in species) {
  intact_tmp <- cluster_breakdown_intact[which(cluster_breakdown_intact$species == sp), c("accession", "genome_percent")]
  pseudo_tmp <- cluster_breakdown_pseudo[which(cluster_breakdown_pseudo$species == sp), c("accession", "genome_percent")]
  
  per_genome_coverages[["intact"]][[sp]] <- aggregate(genome_percent ~ accession, data = intact_tmp, FUN = sum)
  per_genome_coverages[["intergenic.pseudogene"]][[sp]] <- aggregate(genome_percent ~ accession, data = pseudo_tmp, FUN = sum)
}

saveRDS(object = per_genome_coverages,
        file = "/data1/gdouglas/projects/accessory_vs_pseudogene/summary/per_genome_element.type_percent_coverages.rds")

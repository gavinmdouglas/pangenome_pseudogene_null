rm(list = ls(all.names = TRUE))

### Parse CARD annotation and get meaningful gene categories for downstream enrichment tests.
CARD_annot <- read.table("/data1/gdouglas/projects/hgt_fragments/pseudogenes/processed/card_annot/output/filt_pseudogenes_and_intact.txt",
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")

rep2cluster <- read.table("/data1/gdouglas/projects/hgt_fragments/pseudogenes/processed/cdhit_workflow/cluster2rep.tsv.gz",
                          header = FALSE, sep = "\t", row.names = 2, stringsAsFactors = FALSE)

# Realized after the fact that I had run CARD on all sequences rather than just the representative sequences.
# To keep it consistent with how eggNOG was run, restrict parsing to just the representative sequences, and ignore all others.
CARD_annot <- CARD_annot[which(CARD_annot$ORF_ID %in% rownames(rep2cluster)), ]

CARD_clusters <- rep2cluster[CARD_annot$ORF_ID, 1]


raw_unique_drug_classes <- unique(CARD_annot$Drug.Class) 
unique_drug_classes <- as.character()
for (raw_unique_drug_class in raw_unique_drug_classes) {
  unique_drug_classes <- c(unique_drug_classes, strsplit(raw_unique_drug_class, "; ")[[1]])
}
unique_drug_classes <- sort(unique(unique_drug_classes), decreasing = FALSE)

drug_class_clusters <- list()
for (drug_class in unique_drug_classes) {
  drug_class_genes <- CARD_annot[grep(drug_class, CARD_annot$Drug.Class), "ORF_ID"]
  drug_class_clusters[[drug_class]] <- rep2cluster[drug_class_genes, 1]
}


raw_unique_resistance_mechanisms <- unique(CARD_annot$Resistance.Mechanism) 
unique_resistance_mechanisms <- as.character()
for (raw_unique_resistance_mechanism in raw_unique_resistance_mechanisms) {
  unique_resistance_mechanisms <- c(unique_resistance_mechanisms, strsplit(raw_unique_resistance_mechanism, "; ")[[1]])
}
unique_resistance_mechanisms <- sort(unique(unique_resistance_mechanisms), decreasing = FALSE)

resistance_mechanism_clusters <- list()
for (resistance_mechanism in unique_resistance_mechanisms) {
    resistance_mechanism_genes <- CARD_annot[grep(resistance_mechanism, CARD_annot$Resistance.Mechanism), "ORF_ID"]
  resistance_mechanism_clusters[[resistance_mechanism]] <- rep2cluster[resistance_mechanism_genes, 1]
}

CARD_cluster_sets <- list(CARD = CARD_clusters,
                          resistance_mechanism = resistance_mechanism_clusters,
                          drug_class = drug_class_clusters)

saveRDS(object = CARD_cluster_sets,
        file = "/data1/gdouglas/projects/hgt_fragments/pseudogenes/summary/CARD_cluster_sets.rds")

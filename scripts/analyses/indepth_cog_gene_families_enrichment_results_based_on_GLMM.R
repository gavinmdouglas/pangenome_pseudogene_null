# Run enrichments as before, but on individual COGs rather than overall categories.
# Restrict this to cases where the COG category that includes
# each COG was significant in the GLMM.
# Note that this analysis focuses on individual elements rather than clusters.

rm(list = ls(all.names = TRUE))

source("~/scripts/accessory_vs_pseudogenes/R_scripts/functions.R")

identify_enriched_COG_ids <- function(COGs2test,
                                      genes,
                                      background,
                                      min_category_count = 10) {
  
  enrichments_out <- data.frame(matrix(NA, nrow = length(COGs2test), ncol = 8))
  rownames(enrichments_out) <- COGs2test
  colnames(enrichments_out) <- c("category", "genes_num_category", "genes_num_other",
                                 "background_num_category", "background_num_other", "OR", "p", "fdr")
  
  enrichments_out[COGs2test, "category"] <- COGs2test
  
  for (category in rownames(enrichments_out)) {
    
    genes_num_category <- length(which(genes == category))
    genes_num_other <- length(genes) - genes_num_category
    
    background_num_category <- length(which(background == category))
    background_num_other <- length(background) - background_num_category
    
    count_table <- matrix(c(genes_num_category, genes_num_other, background_num_category, background_num_other), nrow = 2, ncol = 2)
    
    enrichments_out[category, c("genes_num_category",
                                "genes_num_other",
                                "background_num_category",
                                "background_num_other")] <- c(genes_num_category,
                                                              genes_num_other,
                                                              background_num_category,
                                                              background_num_other)
    
    if (min(c(genes_num_category + background_num_category, genes_num_other + background_num_other)) >= min_category_count) {
      
      fisher_out <- fisher.test(count_table)
      
      enrichments_out[category, "p"] <- fisher_out$p.value
      
      if (genes_num_other > 0) {
        ratio_numer <- genes_num_category / genes_num_other
      } else {
        ratio_numer <- genes_num_category / 1 
      }
      
      if (background_num_other == 0) {
        ratio_denom <- 1
      } else if(background_num_category == 0) {
        ratio_denom <- 1 / background_num_other
      } else {
        ratio_denom <- background_num_category / background_num_other
      }
      
      
      enrichments_out[category, "OR"] <- ratio_numer / ratio_denom
      
    }
    
  }
  
  enrichments_out$fdr <- p.adjust(enrichments_out$p, "fdr")
  
  rownames(enrichments_out) <- NULL
  
  return(enrichments_out)
  
}

cluster_member_breakdown <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/summary/cluster_member_breakdown.tsv.gz",
                                       header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 2)

element_info <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/summary/element_glm_input.tsv.gz",
                           header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

element_info$element <- gsub("\\.\\d$", "", rownames(element_info))

element_info$cluster <- cluster_member_breakdown[element_info$element, "cluster"]

# Remove rows corresponding to really rare COG categories
to_ignore <- c("A", "B", "Y", "Z")
to_ignore_category_row_i <- which(element_info$COG_category %in% to_ignore)
#length(to_ignore_category_row_i)
element_info <- element_info[-to_ignore_category_row_i, ]

element_info$pseudogenized <- FALSE
element_info[grep("pseudo", rownames(element_info)), "pseudogenized"] <- TRUE

# Set K ("Transcription") as reference COG category.
COG_category_level_order <- sort(unique(element_info$COG_category))
COG_category_level_order <- c("K", COG_category_level_order[-which(COG_category_level_order == "K")])
element_info$COG_category <- factor(element_info$COG_category, levels = COG_category_level_order)

# Reverse redundant_intact_COG, as the majority of elements have a redundant COG, so this makes more sense as the reference.
element_info$NO_redundant_intact_COG <- ! element_info$redundant_intact_COG

COG_sig_dir <- list()
COG_sig_dir[["ultra.cloud"]] <- list()

COG_sig_dir[["ultra.cloud"]]$pos.enriched_COG_false.non.redundant <- c("C", "F", "J", "S", "X")
COG_sig_dir[["ultra.cloud"]]$neg.enriched_COG_false.non.redundant <- c("D")

COG_sig_dir[["ultra.cloud"]]$pos.enriched_COG_true.non.redundant <- as.character()
COG_sig_dir[["ultra.cloud"]]$neg.enriched_COG_true.non.redundant <- c("E", "F", "H", "J", "L", "N", "P", "Q" ,"R", "S", "U", "V", "W", "X")


COG_sig_dir[["other.cloud"]] <- list()

COG_sig_dir[["other.cloud"]]$pos.enriched_COG_false.non.redundant <- c("X")
COG_sig_dir[["other.cloud"]]$neg.enriched_COG_false.non.redundant <- as.character()

COG_sig_dir[["other.cloud"]]$pos.enriched_COG_true.non.redundant <- "C"
COG_sig_dir[["other.cloud"]]$neg.enriched_COG_true.non.redundant <- c("D", "E", "F", "G", "H", "I", "J", "L", "M", "N", "O", "P", "Q" ,"R", "S", "T", "U", "V", "W", "X")


COG_sig_dir[["shell"]] <- list()

COG_sig_dir[["shell"]]$pos.enriched_COG_false.non.redundant <- as.character()
COG_sig_dir[["shell"]]$neg.enriched_COG_false.non.redundant <- c("I", "N", "W")

COG_sig_dir[["shell"]]$pos.enriched_COG_true.non.redundant <- c("S", "V", "W")
COG_sig_dir[["shell"]]$neg.enriched_COG_true.non.redundant <- c("C", "E", "F", "G", "H", "I", "J", "L", "M", "N", "O", "P", "Q", "R", "T", "U", "X")


# Read in COG annot information
COG_annot <- readRDS("/data1/gdouglas/projects/accessory_vs_pseudogene/summary/all_filt_cluster_COG_annot.rds")


COG_gene_map <- read.table("/data1/gdouglas/db/COG_definitions/cog-20.to_category.tsv", header = FALSE, sep = "\t")

COG_category_to_COG <- list()
for (category in unique(COG_gene_map$V2)) {
  COG_category_to_COG[[category]] <- COG_gene_map[which(COG_gene_map$V2 == category), "V1"]
}

species <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/mapfiles/species.txt", sep = "\t")$V1

for (partition in c("ultra.cloud", "other.cloud", "shell")) {
  
    raw_COG_enrichment_out <- list()
    dummy_i <- 1
    
    element_info_subset <- element_info[which(element_info$partition == partition), ]
    
    for (sp in species) {
      
      element_info_subset_sp <- element_info_subset[which(element_info_subset$species == sp), , drop = FALSE]
      
      intact_element_info_subset_sp <- element_info_subset_sp[grep("pseudo", rownames(element_info_subset_sp), invert = TRUE), ]
      pseudo_element_info_subset_sp <- element_info_subset_sp[grep("pseudo", rownames(element_info_subset_sp)), ]
      
      intact_element_info_subset_sp_redundant <- intact_element_info_subset_sp[which(! intact_element_info_subset_sp$NO_redundant_intact_COG), ]
      intact_element_info_subset_sp_non.redundant <- intact_element_info_subset_sp[which(intact_element_info_subset_sp$NO_redundant_intact_COG), ]

      pseudo_element_info_subset_sp_redundant <- pseudo_element_info_subset_sp[which(! pseudo_element_info_subset_sp$NO_redundant_intact_COG), ]
      pseudo_element_info_subset_sp_non.redundant <- pseudo_element_info_subset_sp[which(pseudo_element_info_subset_sp$NO_redundant_intact_COG), ]
        
      sp_intact_COGs_non.redundant <- separate_all_categories(COG_annot[intact_element_info_subset_sp_non.redundant$cluster, "COG"], ",")      
      sp_pseudo_COGs_non.redundant <- separate_all_categories(COG_annot[pseudo_element_info_subset_sp_non.redundant$cluster, "COG"], ",")
      
      sp_intact_COGs_redundant <- separate_all_categories(COG_annot[intact_element_info_subset_sp_redundant$cluster, "COG"], ",")      
      sp_pseudo_COGs_redundant <- separate_all_categories(COG_annot[pseudo_element_info_subset_sp_redundant$cluster, "COG"], ",")
      

      for (cog_category in COG_sig_dir[[partition]]$pos.enriched_COG_false.non.redundant) {
        
        sp_intact_category_COGs <- sp_intact_COGs_redundant[which(sp_intact_COGs_redundant %in% COG_category_to_COG[[cog_category]])]
        
        sp_pseudo_category_COGs <- sp_pseudo_COGs_redundant[which(sp_pseudo_COGs_redundant %in% COG_category_to_COG[[cog_category]])]
        
        
        category_COGs_to_test <- unique(c(sp_pseudo_category_COGs, sp_intact_category_COGs))
        
        enrich_out <- identify_enriched_COG_ids(COGs2test = category_COGs_to_test,
                                                genes = sp_pseudo_COGs_redundant,
                                                background = sp_intact_COGs_redundant,
                                                min_category_count = 10)
        
        enrich_out$sp <- sp
        enrich_out$COG_category <- cog_category
        enrich_out$partition <- partition
        enrich_out$type <- "pos.enriched_COG_false.non.redundant"
        
        raw_COG_enrichment_out[[dummy_i]] <- enrich_out
        
        dummy_i <- dummy_i + 1
        
      }
      
      for (cog_category in COG_sig_dir[[partition]]$neg.enriched_COG_false.non.redundant) {
        
        sp_intact_category_COGs <- sp_intact_COGs_redundant[which(sp_intact_COGs_redundant %in% COG_category_to_COG[[cog_category]])]
        
        sp_pseudo_category_COGs <- sp_pseudo_COGs_redundant[which(sp_pseudo_COGs_redundant %in% COG_category_to_COG[[cog_category]])]
        
        
        category_COGs_to_test <- unique(c(sp_pseudo_category_COGs, sp_intact_category_COGs))
        
        enrich_out <- identify_enriched_COG_ids(COGs2test = category_COGs_to_test,
                                                genes = sp_pseudo_COGs_redundant,
                                                background = sp_intact_COGs_redundant,
                                                min_category_count = 10)
        
        enrich_out$sp <- sp
        enrich_out$COG_category <- cog_category
        enrich_out$partition <- partition
        enrich_out$type <- "neg.enriched_COG_false.non.redundant"
        
        raw_COG_enrichment_out[[dummy_i]] <- enrich_out
        
        dummy_i <- dummy_i + 1
        
      }
      
      for (cog_category in COG_sig_dir[[partition]]$pos.enriched_COG_true.non.redundant) {
        
        sp_intact_category_COGs <- sp_intact_COGs_non.redundant[which(sp_intact_COGs_non.redundant %in% COG_category_to_COG[[cog_category]])]
        
        sp_pseudo_category_COGs <- sp_pseudo_COGs_non.redundant[which(sp_pseudo_COGs_non.redundant %in% COG_category_to_COG[[cog_category]])]
        
        
        category_COGs_to_test <- unique(c(sp_pseudo_category_COGs, sp_intact_category_COGs))
        
        enrich_out <- identify_enriched_COG_ids(COGs2test = category_COGs_to_test,
                                                genes = sp_pseudo_COGs_non.redundant,
                                                background = sp_intact_COGs_non.redundant,
                                                min_category_count = 10)
        
        enrich_out$sp <- sp
        enrich_out$COG_category <- cog_category
        enrich_out$partition <- partition
        enrich_out$type <- "pos.enriched_COG_true.non.redundant"
        
        raw_COG_enrichment_out[[dummy_i]] <- enrich_out
        
        dummy_i <- dummy_i + 1
        
      }
      
      for (cog_category in COG_sig_dir[[partition]]$neg.enriched_COG_true.non.redundant) {
        
        sp_intact_category_COGs <- sp_intact_COGs_non.redundant[which(sp_intact_COGs_non.redundant %in% COG_category_to_COG[[cog_category]])]
        
        sp_pseudo_category_COGs <- sp_pseudo_COGs_non.redundant[which(sp_pseudo_COGs_non.redundant %in% COG_category_to_COG[[cog_category]])]
        
        
        category_COGs_to_test <- unique(c(sp_pseudo_category_COGs, sp_intact_category_COGs))
        
        enrich_out <- identify_enriched_COG_ids(COGs2test = category_COGs_to_test,
                                                genes = sp_pseudo_COGs_non.redundant,
                                                background = sp_intact_COGs_non.redundant,
                                                min_category_count = 10)
        
        enrich_out$sp <- sp
        enrich_out$COG_category <- cog_category
        enrich_out$partition <- partition
        enrich_out$type <- "neg.enriched_COG_true.non.redundant"
        
        raw_COG_enrichment_out[[dummy_i]] <- enrich_out
        
        dummy_i <- dummy_i + 1
        
      }
      
    }

  all_enrich_out <- do.call(rbind, raw_COG_enrichment_out)
  
  all_enrich_out <- all_enrich_out[, c("sp", "partition", "type", "COG_category", "category", "genes_num_category", "genes_num_other",
                                       "background_num_category", "background_num_other", "OR", "p", "fdr")]
  
  enrichment_outfile <- paste("/data1/gdouglas/projects/accessory_vs_pseudogene/summary/glmm_output/COG_enrichment_results/",
                              partition,
                              "-COG-gene-enrichments.tsv",
                              sep = "")
  
  write.table(x = all_enrich_out, file = enrichment_outfile,
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

}

spearman_cor_string <- function(x, y) {
  cor_out <- cor.test(x = x, y = y, method = 'spearman')
  estimate <- format(round(cor_out$estimate, digits=4), nsmall = 4)
  p <- format(round(cor_out$p.value, digits=4), nsmall = 4)

  if (p == '0.0000') {
    return(paste0('rho ~', "'= '*", estimate,"*','", '~ italic(P) ~', "'< 0.0001'"))
  } else {
    return(paste0('rho ~', "'= '*", estimate,"*','", '~ italic(P) ~', "'= '*", p))
  }
}


spearman_ppcor_string <- function(x, y, z) {
  cor_out <-  ppcor::pcor.test(x = x, y = y, z = z, method = 'spearman')
  estimate <- format(round(cor_out$estimate, digits=4), nsmall = 4)
  p <- format(round(cor_out$p.value, digits=4), nsmall = 4)
  
  if (p == '0.0000') {
    return(paste0('Partial ~ rho ~', "'= '*", estimate,"*','", '~ italic(P) ~', "'< 0.0001'"))
  } else {
    return(paste0('Partial ~ rho ~', "'= '*", estimate,"*','", '~ italic(P) ~', "'= '*", p))
  }
}

combine_by_dup_colval <- function(in_df, focal_col, id_col, delimiter = ",", num_cores = 1) {
  rownames(in_df) <- NULL
  dup_ids <- unique(in_df[which(duplicated(in_df[, id_col])), id_col])
  in_df_nondup <- in_df[which(! in_df[, id_col] %in% dup_ids), ]
  in_df_dup <- in_df[which(in_df[, id_col] %in% dup_ids), ]
  
  in_df_dup_merged_raw <- parallel::mclapply(dup_ids,
                                             function(dup_id) {
                                               dup_subset <- in_df_dup[which(in_df_dup[, id_col] == dup_id), ]
                                               all_focal <- as.character()
                                               for (i in 1:nrow(dup_subset)) {
                                                 all_focal <- c(all_focal, strsplit(dup_subset[i, focal_col], delimiter)[[1]])
                                               }
                                               
                                               if (length(all_focal) > 0) {
                                                 return(data.frame(focal_col = paste(sort(unique(all_focal)), collapse = delimiter), id_col = dup_id))
                                               } else {
                                                 return(data.frame(focal_col = "", id_col = dup_id))
                                               }
                                             },
                                             mc.cores = num_cores)
  
  in_df_dup_merged <- do.call(rbind, in_df_dup_merged_raw)
  
  colnames(in_df_dup_merged) <- c(focal_col, id_col)
  
  return(rbind(in_df_nondup, in_df_dup_merged))
}


split_multi_category_rows <- function(in_df, category_col, delimiter = ",", num_cores = 1) {
  
  multi_category_row_i <- grep(delimiter, in_df[, category_col])
  
  in_df_nonmulti <- in_df[-multi_category_row_i, , drop = FALSE]
  in_df_multi <- in_df[multi_category_row_i, , drop = FALSE]
  
  in_df_multi_merged_raw <- parallel::mclapply(1:nrow(in_df_multi),
                                               function(i) {
                                                 row_subset <- in_df_multi[i, , drop = FALSE]
                                                 category_split <- base::strsplit(x = row_subset[, category_col], split = ",")[[1]]
                                                 split_row <- row_subset[rep(1, length(category_split)), , drop = FALSE]
                                                 split_row[, category_col] <- category_split
                                                 return(split_row)
                                               },
                                               mc.cores = num_cores)
  
  in_df_multi_merged <- do.call(rbind, in_df_multi_merged_raw)
  
  return(rbind(in_df_nonmulti, in_df_multi_merged))
}

parse_unique_categories <- function(invec, separator) {
  unique_values <- as.character()
  
  for (raw_value in invec) {
    unique_values <- c(unique_values, strsplit(raw_value, separator)[[1]])
  }
  return(sort(unique(unique_values)))
}


separate_all_categories <- function(invec, separator) {
  all_values <- as.character()
  
  for (raw_value in invec) {
    all_values <- c(all_values, strsplit(raw_value, separator)[[1]])
  }
  return(sort(all_values))
}


genomic_fluidity_from_df <- function(gene_presence) {
  uniqueness_sum <- 0
  
  for (s1 in colnames(gene_presence[, -ncol(gene_presence), drop = FALSE])) {
    
    remaining_samples <- colnames(gene_presence)[(which(colnames(gene_presence) == s1):ncol(gene_presence))]
    
    for (s2 in remaining_samples) {
      
      s1_present <- which(gene_presence[, s1] != "")
      s2_present <- which(gene_presence[, s2] != "")
      
      uniqueness_sum <- uniqueness_sum + length(setdiff(s1_present, s2_present)) / (length(s1_present) + length(s2_present))
      
    }
  }
  
  num_genomes <- ncol(gene_presence)
  
  return((2 / (num_genomes * (num_genomes - 1))) * uniqueness_sum)
  
}


compute_mean_num_singletons_per_combo <- function(gene_presence, k) {
  
  if (k > ncol(gene_presence)) {
    return(list(mean = NA,
                sd = NA))
  }
  
  if (choose(ncol(gene_presence), k) > 100) {
    
    past_combos <- list()
    past_combos[['Start']] <- ''
    
    num_singletons <- as.integer()
    
    for (i in 1:100) {
      
      sample_string <- 'Start'
      
      while (sample_string %in% names(past_combos)) {
        sample_combo <- sample(colnames(gene_presence), k)
        sample_string <- paste(sort(sample_combo), collapse = ',')
      }
      
      past_combos[[sample_string]] <- ''
      
      num_singletons <- c(num_singletons, mean(colSums(gene_presence[which(rowSums(gene_presence[, sample_combo, drop = FALSE] != '') == 1), sample_combo, drop = FALSE] != '')))
      
    }
    
  } else {
    
    combos <- combn(colnames(gene_presence), k)
    num_singletons <- apply(combos, 2, function(x) { mean(colSums(gene_presence[which(rowSums(gene_presence[, x, drop = FALSE] != '') == 1), x, drop = FALSE] != '')) })
    
  }
  
  return(list(mean = mean(num_singletons),
              sd = sd(num_singletons)))
}

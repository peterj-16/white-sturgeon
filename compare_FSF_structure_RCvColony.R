# This function reconstructs full-sibling family structure from relationship coefficients and compares it with full-sib family structures from multiple Colony BestConfig files. Relationship coefficients should be in sample matrix (all X all) and a full-sibship threshold value must be input.

library(readxl)
library(tidyverse)


compare_Col_RC_FSFs <- function(RC_matrix, BestConfig_files, FS_RC_threshold) {
  # Convert the dataframe to matrix and remove row names
  RC_matrix <- RC_matrix %>% remove_rownames %>% column_to_rownames(var = "1")
  RCs_matrix <- as.matrix(RC_matrix)
  
  # Identify full-sib pairs
  full_sib_pairs <- which(RCs_matrix >= FS_RC_threshold, arr.ind = TRUE)
  
  # Convert pairs to a list of sets using individual IDs
  full_sib_pairs <- lapply(1:nrow(full_sib_pairs), function(i) {
    pair <- full_sib_pairs[i, ]
    return(c(rownames(RCs_matrix)[pair[1]], colnames(RCs_matrix)[pair[2]]))
  })
  
  # Function to combine overlapping groups
  combine_groups <- function(groups) {
    merged <- FALSE
    result <- list()
    for (group in groups) {
      added <- FALSE
      for (i in seq_along(result)) {
        if (length(intersect(result[[i]], group)) > 0) {
          result[[i]] <- unique(c(result[[i]], group))
          added <- TRUE
          merged <- TRUE
          break
        }
      }
      if (!added) {
        result <- append(result, list(group))
      }
    }
    return(list(result = result, merged = merged))
  }
  
  # Iteratively combine overlapping full-sib pairs to form full-sib families
  combined <- full_sib_pairs
  repeat {
    combined_result <- combine_groups(combined)
    combined <- combined_result$result
    if (!combined_result$merged) break
  }
  
  # Convert to a more readable format
  full_sibships_combined <- lapply(combined, function(group) {
    sort(unique(group))
  })
  
  # Transform list of full-sib families to df with sample ID and sibship index
  RC_FSFs <- full_sibships_combined
  RC_FSFs_trnsfrmd <- data.frame(FullSibshipIndex = rep(seq_along(RC_FSFs), sapply(RC_FSFs, length)),
                                    SampleID = unlist(RC_FSFs))
  RC_FSFs_trnsfrmd <- RC_FSFs_trnsfrmd %>%
    arrange(SampleID)
  
  # Extract the unique sample IDs
  sample_ids <- unique(RC_FSFs_trnsfrmd$SampleID)
  
  # Create an empty matrix to store the results
  RC_FSF_mtrx <- matrix("Different", nrow = length(sample_ids), ncol = length(sample_ids))
  rownames(RC_FSF_mtrx) <- sample_ids
  colnames(RC_FSF_mtrx) <- sample_ids
  
  # Create a named vector for quick lookup of FullSibshipIndex
  sibship_indices <- setNames(RC_FSFs_trnsfrmd$FullSibshipIndex, RC_FSFs_trnsfrmd$SampleID)
  
  # Iterate through each combination of samples
  for (i in sample_ids) {
    for (j in sample_ids) {
      if (sibship_indices[i] == sibship_indices[j]) {
        RC_FSF_mtrx[i, j] <- "Same"
      }
    }
  }
  
  RC_FSF_mtrx_upper_triangle <- RC_FSF_mtrx[upper.tri(RC_FSF_mtrx, diag = FALSE)]
  
  # Initialize a vector to store similarity proportions
  similarity_proportions <- numeric(length(BestConfig_files))
  
  # Iterate over each BestConfig file
  for (k in seq_along(BestConfig_files)) {
    BestConfig <- read.delim(BestConfig_files[k], sep = "")
    
    # Identify full-sibling groups
    full_sib_groups <- BestConfig %>%
      group_by(MotherID, FatherID) %>%
      summarise(OffspringIDs = list(OffspringID), .groups = 'drop') %>%
      pull(OffspringIDs)
    
    # Transform list of full-sib groups
    col_FSFs_trnsfrmd <- data.frame(
      FullSibshipIndex = rep(seq_along(full_sib_groups), sapply(full_sib_groups, length)),
      SampleID = unlist(full_sib_groups)
    )
    col_FSFs_trnsfrmd <- col_FSFs_trnsfrmd %>%
      arrange(SampleID)
    
    # Extract the unique sample IDs
    sample_ids <- unique(col_FSFs_trnsfrmd$SampleID)
    
    # Create an empty matrix to store the results
    col_FSF_mtrx <- matrix("Different", nrow = length(sample_ids), ncol = length(sample_ids))
    rownames(col_FSF_mtrx) <- sample_ids
    colnames(col_FSF_mtrx) <- sample_ids
    
    # Create a named vector for quick lookup of FullSibshipIndex
    sibship_indices <- setNames(col_FSFs_trnsfrmd$FullSibshipIndex, col_FSFs_trnsfrmd$SampleID)
    
    # Iterate through each combination of samples
    for (i in sample_ids) {
      for (j in sample_ids) {
        if (sibship_indices[i] == sibship_indices[j]) {
          col_FSF_mtrx[i, j] <- "Same"
        }
      }
    }
    
    col_FSF_mtrx_upper_triangle <- col_FSF_mtrx[upper.tri(col_FSF_mtrx, diag = FALSE)]
    
    # Calculate the number and proportion of matching dyads
    num_matching_dyads <- sum(RC_FSF_mtrx_upper_triangle == col_FSF_mtrx_upper_triangle)
    prop_matching_dyads <- num_matching_dyads / length(col_FSF_mtrx_upper_triangle)
    
    # Store the similarity proportion
    similarity_proportions[k] <- prop_matching_dyads
  }
  
  # Calculate the average similarity proportion
  avg_similarity_proportion <- mean(similarity_proportions)
  
  return(list(similarity_proportions = similarity_proportions, avg_similarity_proportion = avg_similarity_proportion))
}

# Input relationship coefficient threshold for full-siblings
FS_RC_threshold <- 0.42

# Input BestConfigs and RCs
BestConfigs <- list.files(path = "~/path/to/your/BestConfigs", full.names = TRUE)
RCs <- read_excel("~/path/to/your/RCs.xlsx", sheet = "RCs")

# Execute script. Gives similarity proportion between RC structure and each Colony structure, and the average among them
RC_Col_similarity <- compare_Col_RC_FSFs(RCs, BestConfigs, FS_RC_threshold)

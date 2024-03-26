# Written by Peter Johnson. 2023. 

# For tetraploid SNP genotypes, this script identifies samples with duplicate genotypes and removes those with lower proportions of SNPs genotyped. 

library(readxl)
library(tidyverse)
library(data.table)
library(tidyr)
library(dplyr)
library(stringr)

# load your genotype file that includes columns 'Sample', 'PercentGenotyped', and 4n SNP genotypes
genotypes <- read_excel("your_genotypes.xlsx")
genos <- genotypes

# set 0's from numeric to character
for (col_name in names(genos)) {
  if (all(genos[[col_name]] == 0)) {
    genos[[col_name]] <- as.character(genos[[col_name]])
  }
}

#split genotypes by allele
genos[genos == "0"] <- "0000"
setDT(genos)
cols <- grep("Atr", names(genos), value = TRUE)
for (i in cols) {
  temp <- tstrsplit(genos[[i]], "")
  set(genos, j = sprintf("%s_%d", i, seq_along(temp)), value = temp)
  set(genos, j = i, value = NULL)
}
genos_split <- genos
genos_split_df <- as.data.frame(genos_split)

data <- genos_split_df


# Function to calculate proportion difference btwn sample pairs, using only genotyped SNPs
calc_prop_diff <- function(genotype1, genotype2) {
  # Find markers that are non-"0" in both genotypes
  total_genotyped_markers <- which(genotype1 != "0" & genotype2 != "0")
  
  # Count markers w/ different genotypes, excluding "0" values
  different_genotyped_markers <- sum(genotype1[total_genotyped_markers] != genotype2[total_genotyped_markers])
  
  prop_diff <- different_genotyped_markers / length(total_genotyped_markers)
  return(prop_diff)
}

# Function to create a list of duplicate pairs based on similarity and PercentGenotyped
get_duplicate_pairs <- function(data, threshold) {
  genotype_cols <- colnames(data)[grepl("Atr", colnames(data))]
  genotype_matrix <- as.matrix(data[, genotype_cols])
  num_samples <- nrow(genotype_matrix)
  
  duplicate_pairs <- list()
  
  for (i in 1:(num_samples - 1)) {
    for (j in (i + 1):num_samples) {
      prop_diff <- calc_prop_diff(genotype_matrix[i, ], genotype_matrix[j, ])
      if (prop_diff <= threshold) {
        pair <- c(data$Sample[i], data$Sample[j])
        pair_with_prop_diff <- c(pair, prop_diff)
        
        # Find differing markers between the two genotypes
        differing_markers <- which(genotype_matrix[i, ] != "0" & genotype_matrix[j, ] != "0" & genotype_matrix[i, ] != genotype_matrix[j, ])
        differing_marker_names <- genotype_cols[differing_markers]
        
        pair_with_differences <- c(pair_with_prop_diff, differing_marker_names)
        duplicate_pairs[[length(duplicate_pairs) + 1]] <- pair_with_differences
      }
    }
  }
  
  return(list(duplicate_pairs = duplicate_pairs))
}

# Call function with threshold
duplicate_pairs <- get_duplicate_pairs(data, threshold = 0.1)

# Function to remove duplicates based on output from get_duplicate_pairs
remove_duplicates <- function(data, duplicate_pairs) {
  duplicates_to_remove <- character(0)
  
  for (pair_info in duplicate_pairs$duplicate_pairs) {
    sample1 <- pair_info[1]
    sample2 <- pair_info[2]
    prop_diff <- pair_info[3]
    
    if (data$PercentGenotyped[data$Sample == sample1] > data$PercentGenotyped[data$Sample == sample2]) {
      duplicates_to_remove <- c(duplicates_to_remove, sample2)
    } else {
      duplicates_to_remove <- c(duplicates_to_remove, sample1)
    }
  }
  
  filtered_data <- data %>% filter(!(Sample %in% duplicates_to_remove))
  return(filtered_data)
}

filtered_data <- remove_duplicates(genos_split_df, duplicate_pairs)


# Remove rows from original 'genotypes' file that are not present in 'filtered_data' based on the 'Sample' column, to get 4n genotypes of 'filtered_data'
filtered_data <- genotypes %>%
  filter(Sample %in% filtered_data$Sample)

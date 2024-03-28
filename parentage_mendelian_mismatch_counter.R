# Written by Peter Johnson. 2024. 

library(readxl)
library(dplyr)
library(tidyverse)

# COUNTS MENDELIAN MISMATCHES OF POSSIBLE PARENT-OFFSPRING PAIRS AND TRIOS (2 PARENTS - 1 OFFSPRING)
# FOR TETRAPLOID SNP GENOTYPES
# REFERENCES RELATIONSHIP COEFFICIENTS


# read in genotypes
offspring_genos <- read_excel("your_offspring_genotypes")
adult_genos <- read_excel("your_adult_genotypes")

# read in adult x offspring relationship coefficients (adult columns X juvenile rows)
adult_juv_RCs <- read.delim("your_relationship_coefficients")

# Initialize an empty dataframe to store pair mismatches
adult_juv_pair_mismatches <- data.frame(offspring = character(), adult = character(), mismatches = numeric(), RC = numeric(), stringsAsFactors = FALSE)

# Function to check if a genotype is missing 
is_missing <- function(x) any(x == "0")

# initialize list for pair mismatch SNPs
pair_mismatch_SNPs <- list()

# find number of mismatches for adult-juvenile pairs; fill 'adult_juv_pair_mismatches' dataframe
for (i in 1:(nrow(offspring_genos))) {
  offspring <- offspring_genos$Sample[i]
  offspring_genotype <- offspring_genos[i, -c(1)]
  
  for (j in 1:(nrow(adult_genos))) {
    adult <- adult_genos$Sample[j]
    adult_genotype <- adult_genos[j, -c(1)]
    
    # Check the relationship coefficient condition (minimum RC seen in known parent-offspring relationships)
    RC <- adult_juv_RCs[i, j]
    if (RC > 0.34) {
      
      # initialize mismatch vector per pair
      mismatch_snps <- character(0)
      
      # Initialize mismatch count for each pair
      mismatch_count <- 0
      
      # Iterate through each SNP
      for (snp in names(adult_genos)[grep("^Atr", names(adult_genos))]) {
        # Skip SNPs with '0' in one or both samples
        if (any(is_missing(offspring_genos[i, snp])) | any(is_missing(adult_genos[j, snp]))) {
          next
        }
        
        # Count Mendelian mismatches for each SNP
        shared_alleles <- sum(strsplit(as.character(offspring_genos[i, snp]), "")[[1]] %in% strsplit(as.character(adult_genos[j, snp]), "")[[1]])
        if (shared_alleles < 2) {
          mismatch_count <- mismatch_count + 1
          mismatch_snps <- c(mismatch_snps, snp)
        }
      }
      
      # Add the result to the adult_offspring_mismatches dataframe
      adult_juv_pair_mismatches <- rbind(adult_juv_pair_mismatches, data.frame(offspring = offspring, adult = adult, mismatches = mismatch_count, RC = RC))
      # Add mismatched SNPs per pair to final list
      pair_mismatch_SNPs[[paste(offspring, adult, sep = ".")]] <- mismatch_snps
    }
  }
}
pairs <- adult_juv_pair_mismatches

# filter adult-juv pairs conservatively 
pairs_filtered <- pairs[pairs$mismatches <= 3 & pairs$RC >= 0.4, ]
row.names(pairs_filtered) <- NULL

# function to make each possible pair of diploid genotypes from a tetraploid genotype
split_genos <- function(genotype) {
  n <- length(genotype)
  poss_splits <- character(n * (n - 1) / 2)
  k <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      split_1 <- paste0(genotype[i], genotype[j])
      split_2 <- paste0(genotype[setdiff(1:n, c(i, j))], collapse = "")
      poss_splits[k] <- paste(split_1, split_2, sep = "_")
      k <- k + 1
    }
  }
  return(poss_splits[1:(k-1)])
}

# make dataframe for trios (possible set of 2 parents and 1 offspring)
pairs_trios <- data.frame(
  offspring = character(),
  adult_1 = character(),
  adult_1_mismatches = numeric(),
  adult_1_RC = numeric(),
  adult_2 = character(),
  adult_2_mismatches = numeric(),
  adult_2_RC = numeric(),
  trio_mismatches = numeric(),
  stringsAsFactors = FALSE
)

# populate trios df w/ all possible pairs from pair info (any 2 pairs w/ same offspring)
for (i in 1:(nrow(pairs_filtered) - 1)) {
  for (j in (i + 1):nrow(pairs_filtered)) {
    
    # Check if the offspring strings are the same
    if (pairs_filtered$offspring[i] == pairs_filtered$offspring[j]) {
      
      # Create a new row for pairs_trios_noMM
      new_row <- data.frame(
        offspring = pairs_filtered$offspring[i],
        adult_1 = pairs_filtered$adult[i],
        adult_1_mismatches = pairs_filtered$mismatches[i],
        adult_1_RC = pairs_filtered$RC[i],
        adult_2 = pairs_filtered$adult[j],
        adult_2_mismatches = pairs_filtered$mismatches[j],
        adult_2_RC = pairs_filtered$RC[j],
        trio_mismatches = NA,
        stringsAsFactors = FALSE
      )
      
      # Bind the new row to pairs_trios_noMM
      pairs_trios <- rbind(pairs_trios, new_row)
    }
  }
}

# make list of genotype dataframes for each trio
trio_genotypes_list <- list()
for (i in 1:nrow(pairs_trios)) {
  offspring <- pairs_trios$offspring[i]
  
  # Extract genotypes for the offspring
  offspring_genotype <- offspring_genos[offspring_genos$Sample == offspring, ]
  
  # Extract genotypes for the first candidate
  first_candidate <- pairs_trios$adult_1[i]
  first_candidate_genotype <- adult_genos[adult_genos$Sample == first_candidate, ]
  
  # Extract genotypes for the second candidate
  second_candidate <- pairs_trios$adult_2[i]
  second_candidate_genotype <- adult_genos[adult_genos$Sample == second_candidate, ]
  
  # Create a dataframe for the trio genotypes with three rows
  trio_genotype_df <- data.frame(rbind(offspring_genotype, first_candidate_genotype, second_candidate_genotype))
  
  # Add trio genotypes dataframe to the list using a unique name
  list_entry_name <- paste(offspring, "trio", i, sep = "_")
  trio_genotypes_list[[list_entry_name]] <- trio_genotype_df
}

# initialize list for SNPs with mismatches at trio level
trio_mismatch_SNPs <- list()

# count trio mismatches for each trio and add to pairs_trios df
for (trio_data in trio_genotypes_list) {
  mismatch_snps <- character(0)
  offspring_genotype <- trio_data[1, ]
  first_candidate_genotype <- trio_data[2, ]
  second_candidate_genotype <- trio_data[3, ]
  
  # Initialize mismatch count and SNPs list for each trio
  mismatch_count <- 0
    
  # Iterate through each SNP
  for (snp in names(trio_data)[grep("^Atr", names(trio_data))]) {
    # Skip SNP if any missing genotypes are present
    if (is_missing(offspring_genotype[snp]) || is_missing(first_candidate_genotype[snp]) || is_missing(second_candidate_genotype[snp])) {
      next
    }
    
    # Get the offspring genotype for the current SNP
    split_offspring_genotype <- unlist(strsplit(as.character(offspring_genotype[snp]), ""))
    
    # Generate all possible splits of the offspring genotype
    poss_splits <- split_genos(split_offspring_genotype)
    
    # Initialize mismatch flag
    mismatch_found <- TRUE
    
    # Extract adult genotypes for the current SNP
    adult_1_genotype <- trio_data[2, -c(1)][[snp]]
    adult_2_genotype <- trio_data[3, -c(1)][[snp]]
    
    # Check each split for Mendelian mismatches
    for (split in poss_splits) {
      # Check if alleles on one side of the underscore are present in adult_1
      # and alleles on the other side are in adult_2
      alleles <- unlist(strsplit(split, "_"))
      if ((grepl(alleles[1], adult_1_genotype) && grepl(alleles[2], adult_2_genotype)) ||
          (grepl(alleles[1], adult_2_genotype) && grepl(alleles[2], adult_1_genotype))) {
        mismatch_found <- FALSE
        break  # No need to check further if a valid split is found
      }
    }
    
    # If none of the splits satisfy the condition, it's a mismatch
    if (mismatch_found) {
      mismatch_count <- mismatch_count + 1
      mismatch_snps <- c(mismatch_snps, snp)
    }
  }
  
  # Find corresponding row in pairs_trios dataframe and update trio_mismatches
  trio_index <- which(pairs_trios$offspring == offspring_genotype$Sample &
                        pairs_trios$adult_1 == first_candidate_genotype$Sample &
                        pairs_trios$adult_2 == second_candidate_genotype$Sample)
  
  if (length(trio_index) > 0) {
    # Update trio_mismatches for the corresponding trio in pairs_trios
    pairs_trios$trio_mismatches[trio_index] <- mismatch_count
  }
  
  trio_mismatch_SNPs[[paste(offspring_genotype$Sample, first_candidate_genotype$Sample, second_candidate_genotype$Sample, sep = ".")]] <- mismatch_snps
}

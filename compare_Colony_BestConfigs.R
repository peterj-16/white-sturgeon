# This script compares BestConfig files from replicate Colony runs for sibship inference, quantifying similarity as proportion of matching pairwise relationships. This is done by making a relationship table for each BestConfig file, which records each pairwise relationship as full-siblings (FS), half-siblings (HS), or unrelated (UR), and comparing each pair of relationship tables.
# Peter Johnson 2024

library(tidyverse)

# Function to create a relationship table for a single Colony BestConfig file
make_relationship_table <- function(file_path) {
  
  # Read in BestConfig 
  BestConfig <- read.delim(file_path, sep = "")
  
  # Create an empty dataframe to become relationship table
  relationship_table <- data.frame(OffspringID_1 = character(),
                                   OffspringID_2 = character(),
                                   Relationship = character(),
                                   stringsAsFactors = FALSE)
  
  # Nested loops to iterate through each offspring pair
  for (i in 1:(nrow(BestConfig) - 1)) {
    for (j in (i + 1):nrow(BestConfig)) {
      
      # Extract offspring and parent information
      offspring1 <- BestConfig$OffspringID[i]
      mother1 <- BestConfig$MotherID[i]
      father1 <- BestConfig$FatherID[i]
      
      offspring2 <- BestConfig$OffspringID[j]
      mother2 <- BestConfig$MotherID[j]
      father2 <- BestConfig$FatherID[j]
      
      # Determine relationship
      if (mother1 == mother2 && father1 == father2) {
        relationship <- "FS"  # Full-siblings
      } else if (mother1 == mother2 || father1 == father2) {
        relationship <- "HS"  # Half-siblings
      } else {
        relationship <- "UR"  # Unrelated
      }
      
      # Add the pair and their relationship to the relationship table
      relationship_table <- rbind(relationship_table, data.frame(OffspringID_1 = offspring1,
                                                                 OffspringID_2 = offspring2,
                                                                 Relationship = relationship,
                                                                 stringsAsFactors = FALSE))
    }
  }
  
  return(relationship_table)
}

# Function that applies previous function to all BestConfig files in a directory and create a relationship table for each 
process_all_BestConfigs <- function(directory_path) {
  
  # List all BestConfig files in directory
  BestConfigs <- list.files(directory_path, pattern = "\\.BestConfig.txt$", full.names = TRUE)
  
  # Loop through each BestConfig file and create a relationship table with appropriate name
  for (file in BestConfigs) {
    
    # Get the base file name without the directory and extension
    file_name <- basename(file)
    file_name <- str_replace(file_name, "\\.BestConfig\\.txt$", "")  
    
    # Generate relationship table name by appending 'relationship_table_' to BestConfig file name
    rt_name <- paste0("relationship_table_", file_name)
    
    # Apply make_relationship_table function to BestConfig file and name resulting relationship table, then put it environment 
    assign(rt_name, make_relationship_table(file), envir = .GlobalEnv)
  }
}


######
# Apply process_all_BestConfigs (which uses make_relationship_table) to directory of BestConfig files
directory_path <- "~/path/to/your/BestConfigs"
process_all_BestConfigs(directory_path)
######

# Function to calculate similarity (proportion of pairs with matching relationships) between two relationship tables 
calc_rt_similarity <- function(table1, table2) {
  # Join the two tables on OffspringID_1 and OffspringID_2
  comparison <- table1 %>%
    inner_join(table2, by = c("OffspringID_1", "OffspringID_2"), suffix = c("_1", "_2"))
  
  # Calculate the proportion of matching relationships
  matching_pairs <- sum(comparison$Relationship_1 == comparison$Relationship_2)
  total_pairs <- nrow(comparison)
  
  prop_matching <- matching_pairs / total_pairs
  return(prop_matching)
}


# Get all the relationship tables from the environment
relationship_tables <- mget(ls(pattern = "^relationship_table_"))

# Create an empty matrix to store similarity proportions between relationship tables
num_tables <- length(relationship_tables)
similarity_matrix <- matrix(NA, nrow = num_tables, ncol = num_tables, 
                            dimnames = list(names(relationship_tables), names(relationship_tables)))

# Loop through each pair of relationship tables and calculate the proportion of offspring pairs with the same relationship
for (i in 1:(num_tables - 1)) {
  for (j in (i + 1):num_tables) {
    table1 <- relationship_tables[[i]]
    table2 <- relationship_tables[[j]]
    
    # Calculate the proportion of matching relationships
    prop_similarity <- calc_rt_similarity(table1, table2)
    
    # Store the result in the matrix
    similarity_matrix[i, j] <- prop_similarity
    similarity_matrix[j, i] <- prop_similarity  # Mirror the result for the symmetric pair
  }
}

# View similarity matrix
similarity_matrix

# Pull similarity proportions for all pairwise comparisons of relationship tables (upper triangle of similarity matrix)
similarity_props_all <- similarity_matrix[upper.tri(similarity_matrix, diag = FALSE)]

# Get their basic stats
similarity_props_all_mean <- mean(similarity_props_all)
similarity_props_all_sd <- sd(similarity_props_all)

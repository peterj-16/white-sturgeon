# Function to find subset of individuals that are unrelated 
# Takes relationship coefficient matrix as input
# Iteratively removes individual that is related to the most other individuals until there are no pairs left with RC above selected threshold value
# Does not return the largest possible unrelated subset; that can be done with igraph - 'largest.independent.vertex.sets', but is much longer process

extract_unrelated_subset <- function(RC_matrix, threshold) {   
  # initialize remaining_samples as all and iteration as 0
  remaining_samples <- colnames(RC_matrix)   
  iteration <- 0    
  
  # loop through iterations until there are no more pairs above relatedness threshold
  while (TRUE) {     

    # set related_pairs as all pairs of individuals within remaining_samples that have RC above threshold
    # store as matrix with 2 cols, 1 for row index and 1 for col index of high RCs in submatrix of remaining samples
    related_pairs <- which(RC_matrix[remaining_samples, remaining_samples] > threshold, arr.ind = TRUE)

    # keep only the upper triangle (row index < col index), removing self-pairs and pair reciprocals
    related_pairs <- related_pairs[related_pairs[, 1] < related_pairs[, 2], , drop = FALSE]
    
    # stop if no more related pairs exist
    if (nrow(related_pairs) == 0) {
      cat("No more related pairs above threshold. Stopping at iteration", iteration, "\n")
      break
    }
    
    # extract sample names from related pairs; store as char vector; a sample is listed x times if present in x related pairs
    related_samples <- c(remaining_samples[related_pairs[, 1]], 
                         remaining_samples[related_pairs[, 2]])     
   
    # count number of related pairs each individual is in
    related_pair_counts <- table(related_samples)

    # identify individual with most relatives for removal
    most_related_sample <- names(which.max(related_pair_counts))  
    
    # remove that individual
    remaining_samples <- remaining_samples[remaining_samples != most_related_sample]
    
    # enumerate iteration
    iteration <- iteration + 1
  }   
  
  return(remaining_samples) 
}

unrelated_subset <- extract_unrelated_subset(your_RC_matrix, threshold = your_value)

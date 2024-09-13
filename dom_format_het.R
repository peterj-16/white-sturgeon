# Peter Johnson 2024
# function to calculate heterozygosity of polyploid SNPs in dominant marker format (alleles treated as pseudo-loci marked present (1), absent (2), or missing data (0))
# input is genotype dataframe with column per marker_allele combo (e.g. SNP_1_A, SNP_1_T, SNP_2_C, SNP_2_T)

dom_format_het <- function(x) {
  
  # get marker names (everything before 2nd underscore in this case)
  markers <- sapply(strsplit(names(df), "_"), function(x) paste(x[1:2], collapse = "_"))
  
  # make blank df w/ correct dimensions to fill with combined dom genotypes per marker
  combined_df <- data.frame(matrix(ncol = length(unique(markers)), nrow = nrow(x)))
  colnames(combined_df) <- unique(markers)
  
  # loop through marker names to fill in empty df 
  for (marker in unique(markers)) {
    cols_to_combine <- names(x)[markers == marker]
    combined_column <- apply(x[cols_to_combine], 1, function(row) paste(row, collapse = ""))
    combined_df[[marker]] <- combined_column
  }
  
  # sub-function to calculate heterozygosity for each marker
  calc_het <- function(column) {
    non_missing_genos <- sum(column != "00")
    het_genos <- sum(column == "11")
    het <- het_genos / non_missing_genos
    return(het)
  }
  
  # apply sub-function and get basic stats
  het <- sapply(combined_df, calc_het)
  het_mean <- mean(het)
  het_sd <- sd(het)
  
  # make output list
  dom_het_stats <- list(
    mean = het_mean,
    SD = het_sd,
    all = het
  )
  return(dom_het_stats)
}

genos_dom <- read.table("~/path/to/your/input.txt", header = TRUE)
results <- dom_format_het(genos_dom)

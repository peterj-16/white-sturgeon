# Written by Peter Johnson and Jeffrey Yen. 2021. 
# Transforms tetrasomic SNP genotypes to dominant marker format (alleles present or absent). 

library(tidyverse)

setwd("/your/working/directory/")

# read in allele header, one row with each SNP_allele (e.g. SNP1_A, SNP1_C, SNP2_C, SNP2_G,...)
allele_header <- read.delim("~/path/to/your/allele_header.txt", header=FALSE)

# read in tetrasomic SNP genotypes, with samples as row names and column names removed
genos <- read_excel("~/path/to/your/genotypes.xlsx", col_names = F)
genos <- genos[-1,]
genos <- data.frame(genos, row.names = 1)

# duplicate every column of genos df and add allele header; makes a column for each SNP_allele to check whether allele is in each genotype
genos_dup <- data.frame(do.call(cbind,rep(genos, each=2)), check.names = F)
allele_header_list <- allele_header[1,]
colnames(genos_dup) <- c(allele_header_list)

# add back row names (sample IDs)
samples_list <- rownames(genos)
rownames(genos_dup) <- c(samples_list)


# function that checks if alleles are present in genotypes and makes new df in dominant format (1 if allele is present, 2 if absent, 0 for missing genotype)
transform_genotypes <- function(data) {
  genos_dom <- c()
  
  for(col in 1:ncol(data)) {
    # extract allele from column name
    allele <- unlist(strsplit(colnames(data)[col], "_"))[3]
    
    # vector to store each transformed marker genotype
    temp_vec <- c()
    
    for(row in 1:nrow(data)) {
      # separate each genotype into respective 4 alleles
      alleles_vector <- unlist(strsplit(data[row, col], ""))  
      
      if(data[row, col] == "0") {  
        temp_vec <- c(temp_vec, "0")
      } else if(allele %in% alleles_vector) {  
        temp_vec <- c(temp_vec, "1")
      } else {  
        temp_vec <- c(temp_vec, "2")
      }
    }
    genos_dom <- cbind(genos_dom, temp_vec)
  }
  
  # Add back row and column names
  colnames(genos_dom) <- colnames(data)
  rownames(genos_dom) <- rownames(data)
  
  return(genos_dom)
}

# apply function to data
genos_dom <- transform_genotypes(genos_dup)

write.table(genos_dom, file = "genos_dom.txt", append = FALSE, quote = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)

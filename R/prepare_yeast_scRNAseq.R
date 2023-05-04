# Filename: prepare_yeast_scRNAseq.R
# Authors: Chris Allsman, Jaeweon Shin, Andrey Zaznaev
# Last edit date: Apr 2023
# Purpose: prepares the counts file to be used by generate_metacells.R from 
#   yeast single-cell RNA-seq experiment by Teixeira et al. 
#   (doi: 10.1093/nar/gkac1041)

# Read in the file
data <- read.table("GSE144820_GSE125162.tsv", 
                   header = TRUE, row.names = 1, sep = "\t")

# Transpose the matrix to have genes*cells dimensions
transposed_data <- t(data)

# Export the counts matrix
write.table(transposed_data, file = "yeast_counts.txt", 
            sep = "\t", quote = FALSE)
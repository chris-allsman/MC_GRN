# Filename: prepare_mouse_scRNAseq.R
# Authors: Chris Allsman, Jaeweon Shin, Andrey Zaznaev
# Last edit date: Apr 2023
# Purpose: prepares the counts file to be used by generate_metacells.R from 
#   mouse single-cell RNA-seq experiment by Jasso et al. 
#   (doi: 10.1371/journal.pbio.3001532)

# Load necessary libraries
library(Seurat)
library(Matrix)

# Import data
barcodes <- read.table("processed_barcodes.tsv",
                       header = FALSE, stringsAsFactors = FALSE)
features <- read.table("processed_features.tsv",
                       header = FALSE, stringsAsFactors = FALSE)
matrix <- readMM("processed_matrix.mtx")
metadata <- read.table("./SCP1711/metadata/metaData_modified.txt", 
                       sep="\t", quote = "", 
                       header = TRUE, stringsAsFactors = FALSE)

# Set row names to feature names and column names to barcode IDs
rownames(matrix) <- features$V2
colnames(matrix) <- barcodes$V1

# Create Seurat object
sc <- CreateSeuratObject(counts = matrix, project = "my_project")

# Add cell type annotations to metadata
metadata <- metadata[metadata$NAME %in% colnames(matrix), ]
names <- metadata$NAME
types <- metadata$cell_type

# Extract condition from barcode IDs
conditions <- sub("M[[:digit:]]+_(.*?)_(.*?)_.+", "\\2_\\1", colnames(matrix))

# Combine condition and annotation to create cell names
cell_names <- paste0(types, "_", conditions)
# Rename column names of matrix with new cell names
colnames(matrix) <- cell_names

# Select only H2O condition for Distal tissue
selected_columns <- grepl("_H2O_Distal$", colnames(matrix))
matrix_subset <- matrix[, selected_columns]
# Remove "_H2O_Distal" from column names, leaving only cell types
new_colnames <- gsub("_H2O_Distal$", "", colnames(matrix_subset))
# Assign the new column names to the sub_matrix
colnames(matrix_subset) <- new_colnames

# Export count matrix with cell names
write.table(matrix_subset, file = "counts_matrix_H2O_Distal_all_genes.txt", 
            quote = FALSE, sep = "\t", col.names = TRUE, row.names=TRUE)

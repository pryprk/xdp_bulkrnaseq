#HeatMap for common genes
# Load necessary libraries
library(pheatmap)
library(RColorBrewer)

# Read the CSV file
log2fc_data <- read.csv("/Users/prakap02/Documents/Postdoc/projects/xdp_paper_1/bulk_rna_seq/DESeq2_Results/STR/heatmap_common_genes.csv", header = TRUE)

# Clean the data: Remove rows with empty gene names and ensure uniqueness
log2fc_data <- log2fc_data[log2fc_data[, 1] != "", ]  # Remove rows with empty gene names
log2fc_data <- log2fc_data[!duplicated(log2fc_data[, 1]), ]  # Remove duplicate gene names

# Ensure unique row names by appending a unique identifier to duplicates
unique_gene_names <- make.unique(log2fc_data[, 1])
rownames(log2fc_data) <- unique_gene_names
log2fc_data <- log2fc_data[, -1]  # Remove the first column (gene names)

# Transpose the matrix to make the heatmap horizontal
log2fc_data_t <- t(log2fc_data)

# Define color palette with 50 shades
cell_colors <- colorRampPalette(c("#1c5894", "#91BFDB", "#F7F7F7", "#e15f58"))(50)

# Create the heatmap
pheatmap(log2fc_data_t,
         color = cell_colors,
         border_color = NA,
         scale = "none",  # No scaling
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Heatmap of Log2 Fold Changes",
         fontsize_row = 10,
         fontsize_col = 10,
         show_rownames = TRUE,
         show_colnames = TRUE,
         cellwidth=10,
         cellheight=20)  # Ensure column names are shown

# Load required libraries
library(ggvenn)
library(ggplot2)
library(pheatmap)
# Load expression matrices from edgeR and limma
# logCPM contains the expression matrix from edgeR (log2-CPM)
# logCPM_limma contains the expression matrix from limma (voom log2-CPM)
# Group is the same for both since we used the same samples
result.e = read.csv("edgeR_results.csv", row.names = 1)
result.l = read.csv("limma_results.csv", row.names = 1)
logc = read.csv("logcpm.csv", row.names = 1, header = T, check.names = F)
group = read.csv("group.csv")$group
group = factor(group)
# Identify top genes from edgeR results
top_edgeR <- rownames(results_edgeR[order(results_edgeR$FDR), ])[1:100]
# Identify top genes from limma results
top_limma <- rownames(results_limma[order(results_limma$adj.P.Val), ])[1:100]
#identify the common results
common_top <- intersect(top_edgeR, top_limma)
# Subset both matrices using the common genes
logCPM_common <- logCPM[common_top, ]
# Prepare annotation for heatmaps
annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(logCPM_common)
# Plot heatmap for common top genes
pheatmap(logCPM_common,
         annotation_col = annotation_col,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = FALSE,
         fontsize_col = 10,
         main = "Heatmap of Common Top Genes (edgeR + limma)")
#Perform PCA on common genes
pca_res <- prcomp(t(logCPM_common), scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  Group = group
)
## Plot PCA for common genes
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA - Common Top Genes",
       x = "PC1",
       y = "PC2") +
  theme(legend.position = "bottom")
identical(rownames(logCPM_sig), rownames(logCPM_limma))
genes_edgeR <- readLines("edgeR.significant.clean.txt")
genes_limma <- readLines("limma.significant.clean.txt")
genes_sig_edgeR <- rownames(sig_edgeR)
genes_sig_limma <- rownames(sig_limma) 
gene_lists <- list(
  edgeR = genes_sig_edgeR,
  limma = genes_sig_limma
)
ggvenn(gene_lists, 
       fill_color = c("skyblue", "palegreen"), 
       stroke_size = 0.5, 
       set_name_size = 4)
common_genes <- intersect(rownames(logCPM), rownames(logcpm_lima))
length(common_genes)  

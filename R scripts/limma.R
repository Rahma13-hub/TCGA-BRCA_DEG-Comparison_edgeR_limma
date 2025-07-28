# Load required libraries
library(limma)
library(edgeR)
library(ggplot2)
library(pheatmap)
library(recount3)
#Load project info for TCGA-BRCA
all_proj <- available_projects()
brca_proj <- subset(all_proj, project == "BRCA" & project_home == "data_sources/tcga")
#Download and load the RangedSummarizedExperiment object
rse_brca <- create_rse(brca_proj)
#Extract counts matrix and sample metadata
counts <- assay(rse_brca, "raw_counts")
metadata <- colData(rse_brca)
#Define sample groups Primary Tumor / Solid Tissue Normal
group <- as.character(colData(rse_brca)$tcga.cgc_sample_sample_type)
# 5. Keep only Tumor and Normal samples
valid_samples <- which(group %in% c("Primary Tumor", "Solid Tissue Normal"))
counts <- counts[, valid_samples]
metadata <- metadata[valid_samples, ]
group <- group[valid_samples]
table(group)
#Select 50 samples from each group
tumor_idx <- which(group == "Primary Tumor")[1:50]
normal_idx <- which(group == "Solid Tissue Normal")[1:50]
selected_idx <- c(tumor_idx, normal_idx)
#Filter to keep only selected samples
counts <- counts[, selected_idx]
metadata <- metadata[selected_idx, ]
group <- group[selected_idx]
#Rename group levels for clarity
group[group == "Primary Tumor"] <- "Tumor"
group[group == "Solid Tissue Normal"] <- "Normal"
group <- factor(group)
#Create DGEList object
dge = DGEList(counts = counts)
#Filter lowly expressed genes
keep <- filterByExpr(dge, group)
dge <- dge[keep, , keep.lib.sizes = FALSE]
#Normalize the counts using TMM
dge <- calcNormFactors(dge)
#Create design matrix without intercept
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
#Apply voom transformation
v <- voom(dge, design, plot = TRUE)
# Extract log2-CPM expression matrix from voom result for PCA, heatmap
logcpm_lima = v$E
#Fit the linear model
fit <- lmFit(v, design)
#Define contrast Tumor vs Normal
contrast <- makeContrasts(Tumor_vs_Normal = Tumor - Normal, levels = design)
#Apply contrast and empirical Bayes moderation
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
#Extract all results
results_limma <- topTable(fit2, number = Inf, sort.by = "P")
#Filter significant genes 
sig_limma <- subset(results_limma, adj.P.Val < 0.05 & abs(logFC) > 1)
#Save results
write.csv(results_limma, file = "limma_results.csv")
write.table(rownames(sig_limma), file = "limma_significant_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
#select top 100 genes
top_genes_limma <- rownames(head(sig_limma[order(-abs(sig_limma$logFC)), ], 100))
logCPM_sig_limma <- logCPM[top_genes_limma, ]
#Create sample annotation for heatmap
annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(logCPM_sig_limma)
#visualize the top 100 genes using heatmap
pheatmap(logCPM_sig_limma,
         annotation_col = annotation_col,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = F,
         fontsize_col = 10)
#use top 100 edgeR genes on limma expression matrix
top_genes_edgeR <- rownames(head(sig_edgeR[order(-abs(sig_edgeR$logFC)), ], 100))
logCPM_limma_edgeRGenes <- logCPM[top_genes_edgeR, ]
annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(logCPM_limma_edgeRGenes)
pheatmap(logCPM_limma_edgeRGenes,
         annotation_col = annotation_col,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = FALSE,
         fontsize_col = 10,
         main = "Heatmap (limma, using edgeR's top genes)")
#visualize DEGs genes using volcano plot
results_limma$significant <- with(results_limma, adj.P.Val < 0.05 & abs(logFC) > 1)
ggplot(results_limma, aes(x = logFC, y = -log10(adj.P.Val), color = significant)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot - limma",
       x = "log2 Fold Change",
       y = "-log10 Adjusted P-value") +
  theme(legend.position = "bottom")
write.csv(logCPM_limma_edgeRGenes, file = "logcpm.csv", row.names = T)
write.csv(data.frame(group), file = "group.csv", row.names = F)
#to make sure the genes are identical in both 
identical(rownames(logCPM), rownames(logcpm_lima))
#PCA
pca_limma <- prcomp(t(logcpm_lima), scale. = TRUE)
plot(pca_limma$x[,1:2], col = group, pch = 19,
     main = "PCA - limma", xlab = "PC1", ylab = "PC2")
legend("topright", legend = levels(group), col = 1:2, pch = 19)
genes = readLines("limma_significant_genes.txt")
genes_clean = sub("\\..*", "", genes)
writeLines(genes_clean, "limma.significant.clean.txt")

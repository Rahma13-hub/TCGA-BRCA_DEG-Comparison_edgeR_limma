#Differential expression analysis using edgeR on TCGA-BRCA data
#load required libraries
library(recount3)
library(edgeR)
library(ggplot2)
library(pheatmap)
#load project info for TCGA-BRCA, it  gets the files from recount3
all_proj <- available_projects()
brca_proj <- subset(all_proj, project == "BRCA" & project_home == "data_sources/tcga")
rse_brca <- create_rse(brca_proj)
#extract counts matrix and sample metadata
counts <- assay(rse_brca, "raw_counts")
metadata <- colData(rse_brca)
#define sample groups (tumor and normal)
group <- as.character(colData(rse_brca)$tcga.cgc_sample_sample_type)
valid_samples <- which(group %in% c("Primary Tumor", "Solid Tissue Normal"))
#filter the data set to include only tumor and normal tissues samples
counts <- counts[, valid_samples]
metadata <- metadata[valid_samples, ]
group <- group[valid_samples]
table(group)
valid_idx <- !is.na(group)
counts <- counts[, valid_idx]
group <- group[valid_idx]
#select the first 50 tumor and 50 normal samples
tumor_idx <- which(group == "Primary Tumor")[1:50]
normal_idx <- which(group == "Solid Tissue Normal")[1:50]
selected_idx <- c(tumor_idx, normal_idx)
#filter to keep only selected samples
counts = counts[, selected_idx]
metadata = metadata[selected_idx, ]
group = group[selected_idx]
#rename group levels
group[group == "Primary Tumor"] = "Tumor"
group[group == "Solid Tissue Normal"] = "Normal"
group = factor(group)
#create DGEList object for edgeR
dge <- DGEList(counts = counts, group = group)
#filter out lowly expressed genes 
keep <- filterByExpr(dge, group)
dge <- dge[keep, , keep.lib.sizes = FALSE]
#normalize the counts using TMM method
dge <- calcNormFactors(dge)
#create the design matrix for linear model
design = model.matrix(~ group)
#estimate dispersion values
dge <- estimateDisp(dge, design)
#fit the model to the data
fit <- glmFit(dge, design)
#perform likelihood ratio test
lrt <- glmLRT(fit, coef = 2)
#extract all differential expression results
results_edgeR <- topTags(lrt, n = Inf)$table
#subset significant genes
sig_edgeR <- subset(results_edgeR, FDR < 0.05 & abs(logFC) > 1)
#save results
write.csv(results_edgeR, file = "edgeR_results.csv")
write.table(rownames(sig_edgeR), file = "edgeR_significant_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
#calculate logcpm for heatmap and PCA
logCPM <- cpm(dge, log = TRUE, prior.count = 1)
#select top 100 genes
top_genes <- rownames(head(sig_edgeR[order(-abs(sig_edgeR$logFC)), ], 100))
logCPM_sig <- logCPM[top_genes, ] 
#create sample annotation for the heatmap
annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(logCPM_sig)
#visualize the top 100 genes using heatmap
pheatmap(logCPM_sig,
         annotation_col = annotation_col,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = FALSE,
         fontsize_col = 10)
#add significant classification
results_edgeR$threshold <- as.factor(
  ifelse(results_edgeR$FDR < 0.05 & abs(results_edgeR$logFC) > 1,
         ifelse(results_edgeR$logFC > 1, "Up", "Down"),
         "Not Significant")
)
#visualize DEGs genes using volcano plot
ggplot(results_edgeR, aes(x = logFC, y = -log10(FDR), color = threshold)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(title = "Volcano Plot (edgeR)",
       x = "log2 Fold Change",
       y = "-log10(FDR)") +
  theme_minimal()
#PCA plot to assess sample clustering 
pca_edger <- prcomp(t(logCPM), scale. = TRUE)
plot(pca_edger$x[,1:2], col = group, pch = 19,
     main = "PCA - edgeR", xlab = "PC1", ylab = "PC2")
legend("topright", legend = levels(group), col = 1:2, pch = 19)
genes = readLines("edgeR_significant_genes.txt")
genes_clean = sub("\\..*", "", genes)
writeLines(genes_clean, "edgeR.significant.clean.txt")

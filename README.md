# Transcriptomic Analysis Summary – TCGA-BRCA Dataset
This report summarizes the RNA-Seq analysis of breast cancer (TCGA-BRCA) samples, which compares tumor and normal tissues. Two differential expression tools, edgeR and limma, were used for the analyses. PCA, volcano plots, heatmaps, and pathway enrichment using g were then used.The profiler.
## Differential Gene Expression
Two statistical methods were used:
- edgeR: More conservative, focused on structural and chromatin-related genes.
- limma: Detected a broader range of DEGs, including more immune/metabolic pathways.
## Principal Component Analysis (PCA)
It was performed to evaluate separation between tumor and normal samples using:
- edgeR DEGs: Clearly distinct separation along PC1.
- limma DEGs: exhibited a strong separation.
- Shared DEGs: Cleanest and most powerful separation.
In conclusion, shared DEGs offer a stronger and more consistent signal for differentiating between tumor and normal profiles.
## Volcano Plots & Heatmaps

- Volcano plots: Displayed significantly up/down-regulated genes.
- Heatmaps:
  - edgeR: More neutral, moderate fold-changes.
  - limma: More intense expression patterns.
  - Shared genes: Produced sharp and clear sample clustering.
  ## g:Profiler Pathway Enrichment Analysis
  Overview:
- 30 pathways shared between edgeR & limma.
- 22 pathways unique to edgeR.
- 25 pathways unique to limma.
- 3854 uniqe gene found with edgeR (27.2%)
- 3968 uniqe gene found with limma (28%)
- 6363 shared genes between edgeR and limma (44.9%)
-  ## Top Shared Pathways & Their Relationship to Breast Cancer
| Pathway | Biological Role in BRCA |
|--------|--------------------------|
| Calcium ion binding | Regulates growth, apoptosis, and motility. Dysregulated calcium signaling is common in breast tumors. |
| Nucleosome assembly | Drives chromatin remodeling, often leading to oncogene activation and silencing of tumor suppressors. |
| Actin cytoskeleton | Promotes cell migration, and structural plasticity essential for metastasis. |
| Protein binding | Indicates altered protein interactions, essential for cancer signaling and progression. |
| Cell recognition | Cancer cells modify surface molecules to escape immune detection and enhance invasive potential. |
| Hemoglobin complex | May relate to hypoxia-driven adaptation in tumors, supporting angiogenesis. |
| Golgi lumen | Involved in processing and secretion of proteins that regulate cell signaling and extracellular communication. |
| Vesicle transport | Tumor cells use vesicles to modulate their environment and promote metastasis. |
| Myelin sheath / Axoneme | Though neural in origin, these structures reflect transport and polarity genes that are hijacked in tumors. |
| Amino acid transport | Supports metabolic rewiring in cancer cells by enhancing nutrient uptake for rapid proliferation. |
# Unique Pathways from edgeR & Their Roles
| Pathway | The Role |
| Structural constituent of muscle | Reflects cytoskeletal elements that may relate to cancer cell contractility and invasion. |
| Transmembrane receptor protein kinase activity | Represents signaling receptors that alterations may enhance oncogenic signaling. |
| Sensory system development | Reflects developmental reprogramming of cancer cells, which adopt stem-like properties. |
| Myelin sheath | Suggests altered intracellular trafficking, potentially enhancing migration. |
| Oxidation-reduction processes | Indicate oxidative stress and redox imbalance common in tumors. |
- edgeR revealed additional metabolic, neural, and structural processes, which may be a sign of modifications in the tumor microenvironment and cell architecture.
# Unique Pathways from limma & Their Roles
| Pathway | The Role |
| Immune response | Indicates active immune signaling, which is often modulated in breast cancer. |
| Hormone binding | Hormone receptors (e.g. estrogen) are central to BRCA classification and therapy. |
| Retinoid binding | Retinoic acid pathways can suppress or promote tumor progression depending on context. |
| Adrenergic receptor activity | May influence tumor stress responses and microenvironment signaling. |
| FGF receptor binding | Related to proliferation and angiogenesis, key target in breast cancer therapeutics. |
- limma has more immune and hormonal signaling, including sensitivity to regulatory and receptor-level changes.
## Biologically significant DEGs are successfully detected by edgeR and limma,
Their unique results offer different perspectives, while their shared findings capture the core concepts of tumor biology:
  - edgeR: Epigenetic and structural alterations.
 - Limma: Immune, hormonal, and metabolic processes.
# Combining the two tools (edgeR/limma) increases assurance and provides a more comprehensive view of tumor biology.
## Supplementary Files
- PCA_plots/: All PCA visuals (edgeR, limma, shared).
- Heatmaps/: Heatmaps of DEGs for (edgeR, limma, shared).
- Volcano_plots/: DEGs distribution visuals.
- gProfiler.pathways.edgeR_limma/: Full pathways.
- venn_diagram_pathways.png: Visual overlap of enriched pathways and genes.
- RScripts for edgeR and limma


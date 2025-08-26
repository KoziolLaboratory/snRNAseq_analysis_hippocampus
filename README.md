# Human snRNA-seq Analysis (GSE175814)

This repository contains an end-to-end R/Seurat workflow to analyze **ALKBH3** and **PINK1** expression in human single-nucleus RNA-seq data (four samples; AD vs Control).  
It loads matrices, merges samples, performs standard preprocessing (Normalize → HVGs → Scale → PCA → UMAP → Neighbors/Clusters), and generates publication-ready figures (UMAP overlays, feature maps, and violin/box plots), plus highlight overlays for **high-ALKBH3/low-PINK1** and **high-PINK1/low-ALKBH3** cell groups.



# Contents
- `snRNAseq-analysis.R` – main analysis script  
- `Figures2/` – auto-created folder containing all exported figures  



# Requirements
- **R ≥ 4.3**  
- R packages:  
  - `Seurat`  
  - `tidyverse`  
  - `ggpubr`  
  - `patchwork`  
  - `dplyr`  
  - `tidyr`  


# Data Layout
**Place the 10x-style matrices under your project root, e.g.:**

D:/Koziol_lab/AD_Project/GSE175814/
  
  GSM5348374_A1_barcodes.tsv.gz
  
  GSM5348374_A1_features.tsv.gz
 
  GSM5348374_A1_matrix.mtx.gz
 
  GSM5348375_A2_barcodes.tsv.gz
  
  GSM5348375_A2_features.tsv.gz
 
  GSM5348375_A2_matrix.mtx.gz
  
  GSM5348376_A3_barcodes.tsv.gz
  
  GSM5348376_A3_features.tsv.gz
  
  GSM5348376_A3_matrix.mtx.gz
  
  GSM5348377_A4_barcodes.tsv.gz
  
  GSM5348377_A4_features.tsv.gz
  
  GSM5348377_A4_matrix.mtx.gz
# Quick Start
Edit the working directory at the top of the script:
setwd("D:/Koziol_lab/AD_Project")
Run the script (analysis.R).
This will:

- `Load A1–A4, tag condition (AD/Control)`

- `Preprocess (Normalize/HVF/Scale/PCA/UMAP/Neighbors/Clusters)`

- `Save figures into Figures2/`

# Thresholds (tuning):

**The script uses simple, user-tunable thresholds:**

**Threshold high/low**

hi_thr  <- 1.0

low_thr <- 1.0

- `Increase/decrease to be more/less stringent.`

- `For robust cut-offs, consider quantile-based thresholds (e.g., q80, q20) or mixture models.`
# Cell-Type Annotation (optional extension)

**You can annotate clusters using your marker panel:**

marker_panels <- list(
  "Excitatory Neuron" = c("SLC17A7","CAMK2A","SATB2"),
  "Inhibitory Neuron" = c("GAD1","GAD2","SLC6A1"),
  "Astrocyte"         = c("AQP4","ALDH1L1","GFAP","SLC1A2"),
  "Microglia"         = c("P2RY12","CX3CR1","TYROBP"),
  "Oligodendrocyte"   = c("MBP","MOG","PLP1","MOBP"),
  "OPC"               = c("PDGFRA","VCAN","CSPG4"),
  "Endothelial"       = c("CLDN5","VWF","FLT1"),
  "Pericyte"          = c("PDGFRB","RGS5","MCAM")
)

Map cluster IDs to labels using canonical marker expression ( via AverageExpression, DotPlot, or manual inspection) and then repeat the ALKBH3/PINK1 summaries per cell type

# Reference Citation:

Soreq, L., Bird, H., Mohamed, W. and Hardy, J., 2023. Single-cell RNA sequencing analysis of human Alzheimer’s disease brain samples reveals neuronal and glial specific cells differential expression. PLoS One, 18(2), p.e0277630.


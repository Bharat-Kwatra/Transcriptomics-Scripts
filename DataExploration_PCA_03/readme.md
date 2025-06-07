# 📊 Step 3: Exploratory Data Analysis and PCA

This module demonstrates data exploration techniques using log2-normalized gene expression values (e.g., log2 CPM). It builds on the output of Step 2.

---

## 🔍 Goals
- Perform **Principal Component Analysis (PCA)** on transcriptomic data
- Generate **interactive plots** using Plotly and DT
- Create **publication-ready figures** using ggplot2 and gt
- Export gene loading values and coordinate matrices
- Conduct **hierarchical clustering** and **correlation analysis**

---

## 📂 Files

| File | Description |
|------|-------------|
| `03_explore_pca.R` | Main R script for PCA and exploration |
| `example_log2cpm_normalized.tsv` | Optional example input |
| `example_targets.tsv` | Sample design metadata |
| `pca_clusters.pdf` | PCA scatter plot with ellipse |
| `pca_loadings_plot.pdf` | Small multiple PCA bar plots |

---

## 📦 Key Libraries
- `tidyverse`, `plotly`, `DT`, `gt`
- `pheatmap`, `corrplot`

---

## 🚀 To Run

1. Ensure the object `log2.cpm.filtered.norm` from Step 2 is loaded
2. Run `03_explore_pca.R` script
3. Outputs:
   - Interactive PCA plots
   - Ranked gene tables
   - PCA loadings as bar charts
   - Heatmap and correlation matrix

---

## 📌 Output Highlights

- `pca_coordinates.tsv` — coordinates for each sample in PC space
- `pca_gene_loadings.tsv` — gene weights (rotation)
- `pca_clusters.pdf` — 2D PCA plot
- `pca_loadings_plot.pdf` — PC1–4 bar plots (small multiples)
- Heatmap of top 50 most variable genes
- Correlation matrix of normalized expression

---

## 🧠 Recommended Next Step

Move to **Step 4: Differential Expression Analysis (DESeq2 or edgeR)**.

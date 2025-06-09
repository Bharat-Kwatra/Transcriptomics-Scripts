# 🧬 Transcriptomics Analysis Pipeline in R

This repository contains a modular, 8-step pipeline for the analysis of RNA-sequencing data. The workflow is designed to guide users from raw transcript quantification through quality control, differential expression analysis, functional enrichment, and visualization using a series of commented R scripts.

## 📊 Workflow Overview

The pipeline is organized into a logical sequence of analysis steps:

1. **Data Import & Annotation:** Load transcript-level data and map it to genes.
2. **Quality Control & Normalization:** Clean, filter, and normalize expression data across samples.
3. **Exploratory Data Analysis:** Visualize data structure and sample relationships using PCA and other multivariate techniques.
4. **Public Data Integration:** Query and extract datasets from public repositories like ARCHS4.
5. **Differential Expression Analysis:** Identify statistically significant differentially expressed genes between experimental conditions.
6. **Gene Module & Heatmap Visualization:** Visualize gene expression patterns with heatmaps and identify co-expressed gene modules.
7. **Functional Enrichment Analysis:** Uncover the biological meaning behind gene lists using GO and Gene Set Enrichment Analysis (GSEA).
8. **Essentials Summary:** A condensed script that combines the core commands from all previous steps for efficient execution.

## 🛠️ Installation

Before starting, install the required R packages from CRAN and Bioconductor:

```r
# Install packages from CRAN
install.packages(c("tidyverse", "rhdf5", "tximport", "beepr", "matrixStats", "cowplot", "DT", "plotly", "gt", "RColorBrewer", "gplots", "devtools"))

# Install packages from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("ensembldb", "EnsDb.Hsapiens.v86", "biomaRt", "edgeR", "limma", "GSEABase", "Biobase", "GSVA", "gprofiler2", "clusterProfiler", "msigdbr", "enrichplot", "IsoformSwitchAnalyzeR"))

# Install packages from GitHub
# devtools::install_github("aljrico/gameofthrones") # Optional fun color palettes
# devtools::install_github("talgalili/d3heatmap") # Optional interactive heatmaps
```

## 🚀 Pipeline Steps

### Step 1: Data Import & Annotation

**Goal:** Import transcript quantification data and summarize it at the gene level with appropriate annotations.

**Key Concepts:**
- `tximport`: Summarize transcript data at the gene level.
- Annotation via Ensembl: Using `ensembldb` or `biomaRt`.

**Key Libraries:** `tximport`, `rhdf5`, `ensembldb`, `biomaRt`, `tidyverse`.

**Input:** `studydesign.txt`, `abundance.tsv` for each sample.  
**Output:** `Txi_gene` object.  
**Run:** `source("Step1_TxImport.R")`

### Step 2: Quality Control & Normalization

**Goal:** Filter low-expression genes and normalize across samples.

**Key Concepts:** CPM filtering, TMM normalization.

**Key Libraries:** `edgeR`, `matrixStats`, `cowplot`, `tidyverse`.

**Input:** `Txi_gene`  
**Output:** `myDGEList.filtered.norm`, `log2.cpm.filtered.norm`  
**Run:** `source("Step2_dataWrangling.R")`

### Step 3: Exploratory Data Analysis

**Goal:** PCA and clustering of sample relationships.

**Key Libraries:** `DT`, `plotly`, `gt`, `tidyverse`

**Input:** `log2.cpm.filtered.norm`  
**Output:** Interactive plots, dendrograms, summary tables  
**Run:** `source("Step3_multivariate.R")`

### Step 4: Public Data Integration

**Goal:** Query ARCHS4 and retrieve external samples.

**Key Libraries:** `rhdf5`, `edgeR`, `tidyverse`

**Input:** GEO/SRA sample IDs, ARCHS4 `.h5` file  
**Output:** Public expression matrix and metadata  
**Run:** `source("Step4_publicData.R")`

### Step 5: Differential Expression Analysis

**Goal:** Identify differentially expressed genes using limma-voom.

**Key Libraries:** `limma`, `edgeR`, `plotly`, `DT`, `gt`

**Input:** `myDGEList.filtered.norm`, contrast matrix  
**Output:** `myTopHits`, `diffGenes`, volcano plot  
**Run:** `source("Step5_diffGenes.R")`

### Step 6: Gene Modules & Heatmap Visualization

**Goal:** Cluster DEGs and visualize gene modules.

**Key Libraries:** `gplots`, `RColorBrewer`, `d3heatmap`

**Input:** `diffGenes`  
**Output:** Heatmaps, gene module lists  
**Run:** `source("Step6_modules.R")`

### Step 7: Functional Enrichment Analysis

**Goal:** GO analysis, GSEA, GSVA for pathway enrichment.

**Key Libraries:** `gprofiler2`, `clusterProfiler`, `GSVA`, `msigdbr`, `limma`, `enrichplot`

**Input:** DEG list or ranked gene list  
**Output:** GO plots, enrichment maps  
**Run:** `source("Step7_functionalEnrichment.R")`

### Step 8: Essentials Summary

**Goal:** Minimal, efficient pipeline summary script.

**Output:** Full analysis in one script  
**Run:** `source("Step8_essentials.R")`

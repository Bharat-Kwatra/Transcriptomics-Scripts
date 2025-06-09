🧬 Transcriptomics Analysis Pipeline in R
This repository contains a modular, 8-step pipeline for the analysis of RNA-sequencing data. The workflow is designed to guide users from raw transcript quantification through quality control, differential expression analysis, functional enrichment, and visualization using a series of commented R scripts.

📊 Workflow Overview
The pipeline is organized into a logical sequence of analysis steps:

Data Import & Annotation: Load transcript-level data and map it to genes.

Quality Control & Normalization: Clean, filter, and normalize expression data across samples.

Exploratory Data Analysis: Visualize data structure and sample relationships using PCA and other multivariate techniques.

Public Data Integration: Query and extract datasets from public repositories like ARCHS4.

Differential Expression Analysis: Identify statistically significant differentially expressed genes between experimental conditions.

Gene Module & Heatmap Visualization: Visualize gene expression patterns with heatmaps and identify co-expressed gene modules.

Functional Enrichment Analysis: Uncover the biological meaning behind gene lists using GO and Gene Set Enrichment Analysis (GSEA).

Essentials Summary: A condensed script that combines the core commands from all previous steps for efficient execution.

🛠️ Installation
Before starting, install the required R packages from CRAN and Bioconductor.

Run the following commands in your R console to install all dependencies:

# Install packages from CRAN
install.packages(c("tidyverse", "rhdf5", "tximport", "beepr", "matrixStats", "cowplot", "DT", "plotly", "gt", "RColorBrewer", "gplots", "devtools"))

# Install packages from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("ensembldb", "EnsDb.Hsapiens.v86", "biomaRt", "edgeR", "limma", "GSEABase", "Biobase", "GSVA", "gprofiler2", "clusterProfiler", "msigdbr", "enrichplot", "IsoformSwitchAnalyzeR"))

# Install packages from GitHub
# devtools::install_github("aljrico/gameofthrones") # Optional fun color palettes
# devtools::install_github("talgalili/d3heatmap") # Optional interactive heatmaps

🚀 Pipeline Steps
Step 1: Data Import & Annotation

Goal: Import transcript quantification data and summarize it at the gene level with appropriate annotations.

Key Concepts:

tximport: A function to read transcript abundance files and correctly summarize them to the gene level, accounting for transcript length.

Annotation: Mapping transcript IDs to gene names using databases like Ensembl (via ensembldb or biomaRt).

Key Libraries: tximport, rhdf5, ensembldb, biomaRt, tidyverse.

Input:

A studydesign.txt file detailing sample IDs and experimental groups.

Transcript abundance files (e.g., abundance.tsv) for each sample.

Output: A list object named Txi_gene containing matrices for gene-level abundance, counts, and length.

How to Run: source("Step1_TxImport.R")

Step 2: Quality Control & Normalization

Goal: Filter out low-expression genes, normalize data to make samples comparable, and visualize the data distributions.

Key Concepts:

Filtering: Removing genes with low counts across samples to improve statistical power. A common method is to keep genes with a Counts Per Million (CPM) > 1 in a minimum number of samples.

Normalization: The Trimmed Mean of M-values (TMM) method from edgeR is used to correct for differences in library size and RNA composition.

Visualization: Violin plots are used to compare sample expression distributions before and after normalization to assess its effectiveness.

Key Libraries: edgeR, matrixStats, cowplot, tidyverse.

Input: The Txi_gene object from Step 1.

Output:

A normalized DGEList object: myDGEList.filtered.norm.

A matrix of filtered, normalized log2 CPM values: log2.cpm.filtered.norm.

How to Run: source("Step2_dataWrangling.R")

Step 3: Exploratory Data Analysis

Goal: Perform unsupervised analysis using PCA and hierarchical clustering to explore sample-to-sample relationships and identify major sources of variation.

Key Concepts:

PCA: A dimensionality-reduction technique that visualizes sample clustering based on expression profiles.

Hierarchical Clustering: An alternative method to group samples in a dendrogram based on similarity.

Interactive Graphics: Using plotly and DT to create interactive plots and tables for deeper data exploration.

Key Libraries: DT, plotly, gt, tidyverse.

Input: The log2.cpm.filtered.norm matrix from Step 2.

Output: Interactive and static PCA plots, dendrograms, and tables of summary statistics (e.g., Log Fold Change between average group expression).

How to Run: source("Step3_multivariate.R")

Step 4: Public Data Integration

Goal: Query and extract public RNA-seq datasets from the ARCHS4 project, which contains uniformly processed data from GEO/SRA.

Key Concepts:

ARCHS4: A resource providing gene-level counts for hundreds of thousands of samples in a single HDF5 file, avoiding manual download and processing.

HDF5: A file format for storing large numerical datasets, which can be queried efficiently in R using the rhdf5 package.

Key Libraries: rhdf5, edgeR, tidyverse.

Input: A list of desired GEO/SRA sample identifiers and the ARCHS4 HDF5 file.

Output: An expression matrix and metadata file for the selected public samples, ready for downstream analysis.

How to Run: source("Step4_publicData.R")

Step 5: Differential Expression Analysis

Goal: Identify genes that are significantly up- or down-regulated between experimental conditions using a linear modeling framework.

Key Concepts:

limma-voom Framework: A robust method that transforms count data for linear modeling, borrowing information across genes to improve statistical power.

Contrast Matrix: A matrix specifying the comparisons of interest (e.g., Treatment - Control).

Volcano Plot: A standard visualization plotting statistical significance against fold change to visualize differentially expressed genes (DEGs).

Key Libraries: limma, edgeR, gt, DT, plotly.

Input: The myDGEList.filtered.norm object and the study design.

Output:

A ranked table of DEGs: myTopHits.

An interactive volcano plot to explore results.

A matrix of expression values for the significant DEGs: diffGenes.

How to Run: source("Step5_diffGenes.R")

Step 6: Gene Modules & Heatmap Visualization

Goal: Visualize expression patterns of DEGs with heatmaps and identify modules of co-expressed genes.

Key Concepts:

Heatmap: A grid where color represents the scaled expression level (Z-score) of genes (rows) across samples (columns).

Gene Modules: Genes with similar expression patterns, often sharing common functions, are identified via hierarchical clustering and can be visualized with colored sidebars.

Key Libraries: gplots, RColorBrewer, d3heatmap.

Input: The diffGenes matrix of DEG expression values from Step 5.

Output: Static and interactive heatmaps, and exported gene lists for each identified module.

How to Run: source("Step6_modules.R")

Step 7: Functional Enrichment Analysis

Goal: Interpret the biological significance of gene lists by testing for the enrichment of known biological pathways and functions.

Key Concepts:

Gene Ontology (GO) Analysis: An over-representation analysis to find GO terms that are significantly enriched in a gene list.

Gene Set Enrichment Analysis (GSEA): A method using a ranked list of all genes to determine if predefined gene sets are enriched at the top or bottom of the list.

GSVA: A single-sample GSEA method that transforms the data into a pathway-by-sample enrichment matrix, allowing for the analysis of pathway activity on a per-sample basis.

Key Libraries: gprofiler2, clusterProfiler, msigdbr, GSVA, limma, enrichplot.

Input: A list of DEGs or a ranked list of all genes.

Output: GO enrichment plots, GSEA enrichment plots, and heatmaps of pathway activity.

How to Run: source("Step7_functionalEnrichment.R")

Step 8: Essentials Summary

Goal: This script serves as a condensed summary, bringing together the minimal, essential commands from all previous steps into a single, runnable file.

Purpose: It acts as a quick-reference "cheatsheet" for the entire pipeline, allowing an experienced user to execute the full workflow efficiently.

How to Run: source("Step8_essentials.R")

# 📘 Step 1: Importing Quantification Data into R

This folder contains a template R script to demonstrate how transcript-level quantification data (e.g., from Kallisto) is imported and annotated in R using `tximport`.

## 🔍 Objective
- Load transcript abundance data from pseudoalignment tools
- Annotate transcripts using Ensembl or BiomaRt
- Format the output for downstream analysis with tools like DESeq2

## 🧪 Sample Files
- `example_studydesign.txt`: Simulated sample sheet with sample names
- `example_abundance.tsv`: Placeholder transcript quantification output

## 🛠️ Key Tools Used
- `tximport`, `rhdf5`, `ensembldb`, `EnsDb.Hsapiens.v86`, `biomaRt`

## 🚀 How to Run

Open the script file in RStudio and run all sections sequentially:

```r
source("import_tx_example.R")

# Introduction to the Step 1 script----
# Step 1 Learning Objectives:
# 1 - Learn the proper 'anatomy' for any R script.
# 2 - Install packages and load libraries into R.
# 3 - Understand file types (e.g., abundance data) and import into R.
# 4 - Learn basic annotation techniques.

# Load packages----
library(rhdf5)
library(tidyverse)
library(tximport)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(beepr)

# Read in study design ----
targets <- read_tsv("example_studydesign.txt")
path <- file.path(targets$sample, "example_abundance.tsv")
all(file.exists(path)) 

# Get annotations using organism-specific package ----
Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))
Tx <- as_tibble(Tx)
Tx <- dplyr::rename(Tx, target_id = tx_id)
Tx <- dplyr::select(Tx, "target_id", "gene_name")

# OPTIONAL: get annotations using BiomaRt----
library(biomaRt)
myMart <- useMart(biomart="ENSEMBL_MART_ENSEMBL")
available.datasets <- listDatasets(myMart)
dog.anno <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset = "clfamiliaris_gene_ensembl")
dog.attributes <- listAttributes(dog.anno)
Tx.dog <- getBM(attributes=c('ensembl_transcript_id_version',
                             'external_gene_name'),
                mart = dog.anno)
Tx.dog <- as_tibble(Tx.dog)
Tx.dog <- dplyr::rename(Tx.dog, target_id = ensembl_transcript_id_version, 
                        gene_name = external_gene_name)

# Import transcript counts using tximport ----
Txi_gene <- tximport(path, 
                     type = "kallisto", 
                     tx2gene = Tx, 
                     txOut = FALSE,
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE)
beep(sound = 6)
class(Txi_gene)
names(Txi_gene)
print("Step 1 complete!")

# Essentials (minimal runnable block) ----
library(tidyverse)
library(tximport)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
targets <- read_tsv("example_studydesign.txt")
path <- file.path(targets$sample, "example_abundance.tsv")
Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))
Tx <- as_tibble(Tx)
Tx <- dplyr::rename(Tx, target_id = tx_id)
Tx <- dplyr::select(Tx, "target_id", "gene_name")
Txi_gene <- tximport(path, 
                     type = "kallisto", 
                     tx2gene = Tx, 
                     txOut = FALSE,
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE)

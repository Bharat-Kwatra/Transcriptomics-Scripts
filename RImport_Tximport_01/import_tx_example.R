# Step 1: Importing and Annotating Transcriptomics Data in R ----

# 🧠 Learning Objectives:
# 1. Learn the structure and logic of an R-based RNA-seq import script.
# 2. Install and load required libraries for transcriptomics analysis.
# 3. Import transcript-level quantification data (from pseudoalignment tools like Kallisto).
# 4. Perform transcript-to-gene mapping using EnsDb or biomaRt.
# 5. Use tximport to summarize counts for downstream DE analysis.

# 📦 Load required packages ----
library(rhdf5) # Required to read HDF5-formatted files (used in bootstrap replicates from tools like Kallisto)
library(tidyverse) # Collection of R packages for data wrangling, visualization, and I/O
library(tximport) # Tool to import and summarize transcript quantification data
library(ensembldb) # Interface to Ensembl-based annotations
library(EnsDb.Hsapiens.v86) # Ensembl database for Homo sapiens (v86) — replace with organism-specific package
library(beepr) # Optional: provides audio notifications

# 📋 Step 1: Load study design ----
# Read a sample sheet that defines sample IDs and paths to abundance files
# This sheet should have a column named 'sample' with sample folder or ID

targets <- read_tsv("example_studydesign.txt")

# 🗂️ Generate file paths to abundance.tsv for each sample ----
path <- file.path(targets$sample, "example_abundance.tsv")

# 🔎 Validate that all files exist
all(file.exists(path)) # Returns TRUE if all paths exist; useful to debug path issues

# 🧬 Step 2: Annotation using EnsDb ----
# Extract transcript-to-gene mapping using EnsDb (preferred for human/mouse/rat)
Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))
Tx <- as_tibble(Tx) # Convert to tibble for tidy manipulation

# 🧾 Rename and rearrange columns as required by tximport
Tx <- dplyr::rename(Tx, target_id = tx_id) # tximport expects column to be named 'target_id'
Tx <- dplyr::select(Tx, "target_id", "gene_name") # Select only the needed columns

# 🔄 OPTIONAL: Annotation using biomaRt for non-model organisms ----
library(biomaRt) # Used to fetch Ensembl data programmatically from online

# Define biomart object (defaults to Ensembl main site)
myMart <- useMart(biomart="ENSEMBL_MART_ENSEMBL")

# View all available datasets (useful to find organism-specific dataset)
available.datasets <- listDatasets(myMart)

# Example: Fetch dog annotations
# Replace 'clfamiliaris_gene_ensembl' with appropriate dataset for your organism
dog.anno <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset = "clfamiliaris_gene_ensembl")
dog.attributes <- listAttributes(dog.anno) # View available attributes

# Fetch transcript-to-gene mappings from biomaRt for dog
Tx.dog <- getBM(attributes=c('ensembl_transcript_id_version', 'external_gene_name'),
                mart = dog.anno)
Tx.dog <- as_tibble(Tx.dog)
Tx.dog <- dplyr::rename(Tx.dog, target_id = ensembl_transcript_id_version,
                        gene_name = external_gene_name)

# 📥 Step 3: Import abundance files using tximport ----
# This function summarizes transcript-level quantifications to gene-level
Txi_gene <- tximport(path, 
                     type = "kallisto", # Set type as "kallisto"
                     tx2gene = Tx,      # Provide transcript-to-gene mapping
                     txOut = FALSE,     # Set to TRUE for transcript-level; FALSE for gene-level
                     countsFromAbundance = "lengthScaledTPM", # Use scaled TPM for counts
                     ignoreTxVersion = TRUE) # Remove version suffix in transcript IDs for consistency

# 🔔 Audio alert when done
beep(sound = 6)

# 🧾 Inspect structure of imported object
class(Txi_gene)  # Should be "list"
names(Txi_gene)  # Will include 'abundance', 'counts', 'length'

print("Step 1 complete!")

# 💾 OPTIONAL: Export counts or TPM ----
# Txi_trans <- as_tibble(Txi_trans$counts, rownames = "target_id")
# Txi_trans <- left_join(Txi_trans, Tx) # Attach gene names for reporting

# 🔁 The Essentials Block (Minimal runnable version of the above) ----
library(tidyverse)
library(tximport)
library(ensembldb)
library(EnsDb.Hsapiens.v86)

# Read study design
targets <- read_tsv("example_studydesign.txt")
path <- file.path(targets$sample, "example_abundance.tsv")

# Get annotation
Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))
Tx <- as_tibble(Tx)
Tx <- dplyr::rename(Tx, target_id = tx_id)
Tx <- dplyr::select(Tx, "target_id", "gene_name")

# Import quantification to R
Txi_gene <- tximport(path, 
                     type = "kallisto", 
                     tx2gene = Tx, 
                     txOut = FALSE, 
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE)

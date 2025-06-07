# Step 2: Quality Control, Filtering & Normalization

## 🎯 Goal
To explore, filter, and normalize transcriptomics data before downstream analysis.

## 🧠 Concepts
- **Gene expression is noisy**: Many genes have zero or very low expression across samples.
- **Filtering**: We remove genes that are unlikely to contribute to biological signal.
- **Normalization**: Accounts for differences in sequencing depth and composition.
- **Visualization**: We use violin plots to compare distributions before and after filtering/normalization.

## 🛠 Tools Used
- `edgeR::DGEList()`: A structured object to hold count data and metadata.
- `edgeR::cpm()`: Computes counts per million (raw/log2).
- `edgeR::calcNormFactors()`: Normalizes data using TMM method.
- `ggplot2`: Plots the distribution of gene expression values.
- `cowplot`: Combines multiple plots into a single figure.

## 🗂 Output
- **Violin plots**: Raw, filtered, and TMM-normalized gene expression profiles.
- **Filtered DGEList**: Ready for differential expression analysis.

## 🧪 Key Output Objects
- `myDGEList`: Raw DGEList object
- `myDGEList.filtered`: Filtered DGEList
- `myDGEList.filtered.norm`: Normalized DGEList
- `log2.cpm.filtered.norm`: Matrix of normalized log2 CPM values

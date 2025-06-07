# Step 2: Filtering, Normalization, and Visualization ----

# 🧠 Learning Objectives:
# 1. Filter and normalize your expression matrix
# 2. Visualize before/after using violin plots and summary statistics
# 3. Tidy the gene expression data to long format for visualization

# 📦 Load Required Libraries ----
library(tidyverse)      # For plotting and data wrangling
library(edgeR)          # For DGEList, CPM, normalization
library(matrixStats)    # Efficient row-/col-wise stats
library(cowplot)        # For combining multiple ggplots into one figure

# 📊 Extract counts and TPM ----
myTPM <- Txi_gene$abundance
myCounts <- Txi_gene$counts
colSums(myTPM)
colSums(myCounts)

# 📋 Use study design to capture sample labels ----
sampleLabels <- targets$sample

# 📈 Compute row-wise statistics (SD, mean, median) for TPM matrix ----
myTPM.stats <- transform(myTPM, 
                         SD=rowSds(myTPM), 
                         AVG=rowMeans(myTPM),
                         MED=rowMedians(myTPM))
head(myTPM.stats)

# 🎨 Create TPM scatter plot with ggplot2 ----
ggplot(myTPM.stats) + 
  aes(x = SD, y = MED) +
  geom_point(shape=25, size=3)

# 🔍 Add more layers to plot for better insight ----
ggplot(myTPM.stats) + 
  aes(x = SD, y = MED) +
  geom_point(shape=16, size=2) +
  geom_smooth(method=lm) +
  geom_hex(show.legend = FALSE) +
  labs(y="Median", x = "Standard deviation",
       title="Transcripts per million (TPM)",
       subtitle="unfiltered, non-normalized data",
       caption="DIYtranscriptomics - Spring 2020") +
  theme_classic() +
  theme_dark() + 
  theme_bw()

# 🧬 Construct a DGEList object from raw counts ----
myDGEList <- DGEList(myCounts)
myDGEList
save(myDGEList, file = "example_myDGEList.RData")
load(file = "example_myDGEList.RData")

# 📉 Convert counts to log2 CPM ----
cpm <- edgeR::cpm(myDGEList)
log2.cpm <- edgeR::cpm(myDGEList, log=TRUE)

log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
colnames(log2.cpm.df) <- c("geneID", sampleLabels)

# 🧹 Reshape CPM to long format for plotting ----
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, 
                                  cols = -1, 
                                  names_to = "samples", 
                                  values_to = "expression")

# 🎻 Visualize log2 CPM distribution (unfiltered) ----
ggplot(log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# 🔎 Check how many genes have 0 counts in all samples ----
table(rowSums(myDGEList$counts == 0) == ncol(myDGEList$counts))

# 🧽 Filter low-expression genes ----
keepers <- rowSums(cpm > 1) >= 5
myDGEList.filtered <- myDGEList[keepers,]
dim(myDGEList.filtered)

# 🪄 Recalculate log2 CPM after filtering ----
log2.cpm.filtered <- edgeR::cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)

log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df, 
                                           cols = -1, 
                                           names_to = "samples", 
                                           values_to = "expression")

# 🎻 Violin plot: filtered, unnormalized ----
p1 <- ggplot(log2.cpm.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# ⚖️ Normalize using TMM ----
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- edgeR::cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)

log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df, 
                                                cols = -1, 
                                                names_to = "samples", 
                                                values_to = "expression")

# 🎻 Violin plot: filtered and normalized ----
p2 <- ggplot(log2.cpm.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# 🖼️ Combine plots together ----
p3 <- ggplot(log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

plot_grid(p3, p1, p2, labels = c('A', 'B', 'C'), label_size = 12)

print("Step 2 complete!")

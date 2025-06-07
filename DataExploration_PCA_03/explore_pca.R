# Step 3: Exploratory Analysis and PCA Visualization ----

# 🧠 Learning Objectives:
# 1. Perform PCA and hierarchical clustering on normalized expression data.
# 2. Use dplyr to derive logFC and summary metrics.
# 3. Generate high-quality visual and interactive tables and plots.
# 4. Understand data structure visually and statistically before DE analysis.

# 📦 Load Required Libraries ----
library(tidyverse)  # Data wrangling and ggplot2
library(DT)         # Interactive HTML tables
library(plotly)     # Interactive ggplot-to-web plots
library(gt)         # Publication-quality tables
library(pheatmap)   # Heatmaps of variable genes
library(corrplot)   # Correlation matrix plots

# 📋 Define Group Variable from Study Design ----
group <- factor(targets$group)  # e.g., "Healthy" vs "Disease"

# ✅ Verify Input is Available
# log2.cpm.filtered.norm.df should be present from Step 2

# 🧭 Hierarchical Clustering (Distance & Agglomeration Method Sensitive) ----
distance <- dist(t(log2.cpm.filtered.norm), method = "maximum")  # Try 'euclidean' or 'manhattan'
clusters <- hclust(distance, method = "average")
plot(clusters, labels = sampleLabels, main = "Hierarchical Clustering of Samples")

# 🧠 PCA to Capture Global Expression Patterns ----
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale. = FALSE, retx = TRUE)
pc.var <- pca.res$sdev^2
pc.per <- round(pc.var / sum(pc.var) * 100, 1)  # % variance per PC
screeplot(pca.res, main = "Screeplot of Principal Components")

# 📈 Create PCA Plot ----
pca.res.df <- as_tibble(pca.res$x)
pca.plot <- ggplot(pca.res.df) +
  aes(x = PC1, y = PC2, label = sampleLabels, color = group) +
  geom_point(size = 4) +
  stat_ellipse() +
  xlab(paste0("PC1 (", pc.per[1], "% )")) +
  ylab(paste0("PC2 (", pc.per[2], "% )")) +
  labs(title = "PCA of Sample Clustering", caption = paste0("Generated on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

# 🌐 Interactive PCA
ggplotly(pca.plot)

# 🔍 PCA Loadings in Small Multiples ----
pca.small <- pca.res$x[, 1:4] %>%
  as_tibble() %>%
  add_column(sample = sampleLabels, group = group) %>%
  pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "loadings")

ggplot(pca.small) +
  aes(x = sample, y = loadings, fill = group) +
  geom_bar(stat = "identity") +
  facet_wrap(~PC) +
  coord_flip() +
  labs(title = "Sample Contributions to Top 4 PCs") +
  theme_bw()

# 🔬 Summary Stats: Healthy vs Disease Avg, LogFC ----
mydata.df <- log2.cpm.filtered.norm.df %>%
  mutate(healthy.AVG = rowMeans(select(., starts_with("HS"))),
         disease.AVG = rowMeans(select(., starts_with("CL"))),
         LogFC = disease.AVG - healthy.AVG) %>%
  mutate_if(is.numeric, round, 2)

# 🧪 Top Differential Genes (by LogFC) ----
mydata.sort <- mydata.df %>% arrange(desc(LogFC)) %>% select(geneID, LogFC)

# 🎯 Select Specific Biomarkers of Interest ----
mydata.filter <- mydata.df %>%
  filter(geneID %in% c("MMP1", "GZMB", "IL1B", "GNLY", "IFNG", "CCL4", "PRF1", "APOBEC3A", "UNC13A")) %>%
  select(geneID, healthy.AVG, disease.AVG, LogFC) %>%
  arrange(desc(LogFC))

# 🔍 Search Genes Using Pattern (e.g., Interferons or Chemokines)
mydata.grep <- mydata.df %>%
  filter(str_detect(geneID, "CXCL|IFI")) %>%
  select(geneID, healthy.AVG, disease.AVG, LogFC) %>%
  arrange(desc(geneID))

# 🏷️ Publish-Ready Summary Table Using gt ----
mydata.filter %>%
  gt() %>%
  fmt_number(columns = 2:4, decimals = 1) %>%
  tab_header(title = md("**Key Immunogenes**"), subtitle = md("*filtered from expression matrix*")) %>%
  tab_footnote("Potential druggable targets", locations = cells_body(columns = geneID, rows = 2:5)) %>%
  tab_source_note(md("Data from Example Pipeline v1"))

# 💡 Interactive Table for Exploratory Filtering ----
datatable(mydata.df[, c(1, 12:14)],
          extensions = c('KeyTable', "FixedHeader"),
          filter = 'top',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100")))

# 📊 Interactive Scatter Plot ----
myplot <- ggplot(mydata.df) +
  aes(x = healthy.AVG, y = disease.AVG, text = paste("Symbol:", geneID)) +
  geom_point(shape = 16, size = 1) +
  labs(title = "Healthy vs Disease Expression", x = "Healthy AVG", y = "Disease AVG") +
  theme_bw()

ggplotly(myplot)

# 🧬 Heatmap of Top 50 Most Variable Genes ----
topVarGenes <- order(rowSds(as.matrix(log2.cpm.filtered.norm)), decreasing = TRUE)[1:50]
heatmap_data <- log2.cpm.filtered.norm[topVarGenes, ]
heatmap_data_scaled <- t(scale(t(heatmap_data)))
pheatmap(heatmap_data_scaled,
         show_colnames = TRUE,
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         annotation_col = targets["group"],
         main = "Top 50 Variable Genes Heatmap")

# 📈 Correlation Matrix of Samples ----
corr_matrix <- cor(log2.cpm.filtered.norm)
corrplot(corr_matrix, method = "color", type = "upper", tl.cex = 0.7, title = "Sample Correlation Matrix")

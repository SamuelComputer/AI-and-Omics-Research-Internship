#===============================================================================
#           AI AND OMICS RESEARCH INTERNSHIP
#===============================================================================
# Module II – Introduction to Genomics Data Analysis
#-------------------------------------------------------------------------------
# Class 3C: Differential Gene Expression Analysis in R

#               GEO accession ID: GSE42568
#===============================================================================

# ASSIGNMENT 5

# Check working directory
getwd()
[1] "C:/Users/User/Desktop/AI_Omics_Internship_2025"

# Organizing Project Class 3C
# File menu > New Project > New Directory > New Project > Class 3C

# Set working directory
setwd("C:/Users/User/Desktop/AI_Omics_Internship_2025/Class 3C")

# Create project folders
dir.create("Raw_Data")
dir.create("Plots")
dir.create("Results")
dir.create("Scripts")

# Install and load packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "GEOquery",
  "limma",
  "AnnotationDbi",
  "hgu133plus2.db",
  "pheatmap",
  "ggplot2"
))

library(GEOquery)
library(limma)
library(AnnotationDbi)
library(hgu133plus2.db)
library(pheatmap)
library(ggplot2)

# Download GSE42568 from GEO
gse <- getGEO("GSE42568", GSEMatrix = TRUE)
length(gse)

# Extract the Expression Matrix
gse <- gse[[1]]   
exprs_matrix <- exprs(gse)

dim(exprs_matrix)

# Extract sample Metadata
pheno <- pData(gse)
head(pheno[, 1:5])

colnames(pheno)

# Define Sample Groups
table(pheno$`tissue:ch1`)

# Create the group factor (Cancer vs Normal)
group <- factor(ifelse(
  grepl("normal", pheno$`tissue:ch1`, ignore.case = TRUE),
  "Normal",
  "Cancer"
))

table(group)

# Create Design Matrix
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
design

# Map Probe IDs → Gene Symbols (AnnotationDbi)
probe_ids <- rownames(exprs_matrix)

gene_symbols <- mapIds(
  hgu133plus2.db,
  keys = probe_ids,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

# Remove probes without gene symbols
exprs_annot <- exprs_matrix[!is.na(gene_symbols), ]
gene_symbols <- gene_symbols[!is.na(gene_symbols)]

exprs_annot <- cbind(GeneSymbol = gene_symbols, exprs_annot)

# Create numeric indices
idx <- which(!is.na(gene_symbols))

# Check
length(idx)

# Subset BOTH objects using the SAME indices
exprs_filtered <- exprs_matrix[idx, ]
gene_symbols_filtered <- gene_symbols[idx]

# Verify
nrow(exprs_filtered)
length(gene_symbols_filtered)

# Create exprs_annot
exprs_annot <- data.frame(
  GeneSymbol = gene_symbols_filtered,
  exprs_filtered,
  check.names = FALSE
)

class(exprs_annot)
dim(exprs_annot)

# Check duplicate probes per gene
table(table(exprs_annot$GeneSymbol))[1:10]

# Collapse duplicate probes
exprs_unique <- avereps(
  as.matrix(exprs_annot[, -1]),
  ID = exprs_annot$GeneSymbol
)

dim(exprs_unique)

## Differential Expression Analysis with Limma
# Fit the linear model
fit <- lmFit(exprs_unique, design)

# Define the contrast (Cancer vs Normal)
contrast.matrix <- makeContrasts(
  Cancer_vs_Normal = Cancer - Normal,
  levels = design
)

# Apply contrast and empirical Bayes moderation
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract full DEG table
deg <- topTable(
  fit2,
  number = Inf,
  adjust.method = "BH"
)

head(deg)

# Save DEG Results (CSV Files)
write.csv(deg, "Results/DEG_complete.csv")

# Define thresholds
upregulated <- deg[deg$logFC > 1 & deg$adj.P.Val < 0.05, ]
downregulated <- deg[deg$logFC < -1 & deg$adj.P.Val < 0.05, ]

# Save results
write.csv(upregulated, "Results/DEG_upregulated.csv")
write.csv(downregulated, "Results/DEG_downregulated.csv")

# Count DEGs
nrow(upregulated)
nrow(downregulated)

# Volcano Plot
deg$Significance <- "Not Significant"
deg$Significance[deg$logFC > 1 & deg$adj.P.Val < 0.05] <- "Upregulated"
deg$Significance[deg$logFC < -1 & deg$adj.P.Val < 0.05] <- "Downregulated"

volcano <- ggplot(deg, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "Not Significant" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Cancer vs Normal")

ggsave("Plots/Volcano_plot.png", volcano, width = 7, height = 5)

# Heatmap of Top 25 DEGs
top25_genes <- rownames(deg)[1:25]
heatmap_data <- exprs_unique[top25_genes, ]

png("Plots/Heatmap_top25_DEGs.png", width = 1000, height = 800)
pheatmap(
  heatmap_data,
  scale = "row",
  show_rownames = TRUE,
  main = "Top 25 Differentially Expressed Genes"
)
dev.off()

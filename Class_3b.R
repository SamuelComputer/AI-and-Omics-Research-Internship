#-------------------------------------------------------------
# AI AND OMICS RESEARCH INTERNSHIP
# Module II – Class 3B
# Preprocessing and Normalization of Microarray Data in R
# Dataset: GSE3292 (Affymetrix CEL files)
#-------------------------------------------------------------

# ASSIGNMENT 4

# Check working directory
getwd()
[1] "C:/Users/User/Desktop/AI_Omics_Internship_2025"

# Set working directory
setwd("C:/Users/User/Desktop/AI_Omics_Internship_2025/Class 3B")

# Organizing Project Class 3B
#File menu > New Project > New Directory > New Project > Class 3B

# Create project folders
dir.create("Raw_Data")
dir.create("Plots_Class3B_Microarray_QC")
dir.create("Results")
dir.create("Scripts")

# Install & load packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("affy", "GEOquery", "limma"), ask = FALSE, update = FALSE)
install.packages(c("matrixStats", "dplyr"))

library(affy)
library(GEOquery)
library(limma)
library(matrixStats)
library(dplyr)

# Define CEL file path
cel_path <- "Raw_Data/CEL_Files"

# List CEL files
cel_files <- list.files(cel_path, pattern = "\\.CEL(\\.gz)?$", full.names = TRUE)
length(cel_files)

# Remove corrupted CEL files (VERY IMPORTANT)
good_files <- cel_files[
  sapply(cel_files, function(f) {
    tryCatch({
      affyio::read.celfile(f)
      TRUE
    }, error = function(e) FALSE)
  })
]

length(good_files)

# Read only valid CEL files
raw_data <- ReadAffy(filenames = good_files)

# 9. Quality Control – Raw data boxplot
png("Plots/raw_boxplot.png", width = 1200, height = 800)
boxplot(exprs(raw_data),
        outline = FALSE,
        las = 2,
        main = "Raw CEL Data – Boxplot")
dev.off()

# Normalize using RMA
norm_data <- rma(raw_data)

# Extract normalized expression matrix
expression_data <- exprs(norm_data)

# QC – Density plot after normalization
png("Plots/normalized_density.png", width = 1200, height = 800)
plotDensities(expression_data,
              main = "RMA Normalized Expression Density")
dev.off()

# PCA after normalization
pca <- prcomp(t(expression_data), scale. = TRUE)

pdf("Plots/pca_normalized.pdf", width = 10, height = 8)
plot(pca$x[,1], pca$x[,2],
     pch = 19,
     xlab = "PC1",
     ylab = "PC2",
     main = "PCA of RMA-Normalized Data")
grid()
dev.off()

# Filter low-expression probes
row_medians <- rowMedians(expression_data)
threshold <- 3.5

filtered_data <- expression_data[row_medians > threshold, ]
dim(filtered_data)

# QC – Filtered data boxplot
png("Plots/filtered_boxplot.png", width = 1200, height = 800)
boxplot(filtered_data,
        outline = FALSE,
        las = 2,
        main = "Filtered Expression Data")
dev.off()

# Save filtered expression matrix
write.csv(filtered_data,
          "Results/filtered_expression_data.csv")

# Histogram of row medians
png("Plots/row_median_hist.png", width = 1200, height = 800)
hist(row_medians,
     breaks = 100,
     freq = FALSE,
     main = "Row Median Distribution",
     xlab = "Median Expression")
abline(v = threshold, col = "red", lwd = 2)
dev.off()

# Count total, disease, normal samples
# Sample count summary
total_samples <- ncol(expression_data)

disease_samples <- total_samples
normal_samples <- 0

cat("Total samples:", total_samples, "\n")
cat("Disease samples:", disease_samples, "\n")
cat("Normal samples:", normal_samples, "\n")

# Find number of probes before filtering

library(affy)
library(matrixStats)

# Path to your CEL files
cel_path <- "C:/Users/User/Desktop/AI_Omics_Internship_2025/Class 3B/Raw_Data/CEL_Files"

# List all CEL files (skip any corrupted files if needed)
cel_files <- list.files(cel_path, pattern = "*.CEL.gz", full.names = TRUE)
cel_files <- cel_files[!grepl("GSM73683.CEL.gz", cel_files)]  # example

# Read CEL files
raw_data <- ReadAffy(filenames = cel_files)

# Normalize with RMA
normalized_data <- rma(raw_data)
processed_data <- exprs(normalized_data)

# Now you can get number of probes
n_probes_before <- nrow(processed_data)
cat("Number of probes before filtering:", n_probes_before, "\n")

# Boxplot Method
boxplot(processed_data, outline = TRUE, las = 2,
        main = "Boxplot of Normalized Data")

# PCA Method
pca <- prcomp(t(processed_data), scale. = TRUE)

plot(pca$x[,1], pca$x[,2],
     xlab = "PC1", ylab = "PC2",
     main = "PCA of Normalized Expression Data",
     pch = 19, col = "blue")
grid()

# Calculate z value
sample_medians <- apply(processed_data, 2, median)
z_scores <- scale(sample_medians)

which(abs(z_scores) > 3)  # samples with median >3 SD from mean

# list the outlier sample names
colnames(processed_data)[which(abs(z_scores) > 3)]

# Create and save boxplot of normalized data
png("Plots_Class3B_Microarray_QC/normalized_boxplot.png", width = 1200, height = 800)
boxplot(processed_data, outline = FALSE, las = 2,
        main = "Boxplot of Normalized Expression Data",
        col = "lightblue")
dev.off()

# Perform PCA on normalized data
pca <- prcomp(t(processed_data), scale. = TRUE)

# Save PCA plot as PDF
pdf("Plots_Class3B_Microarray_QC/normalized_pca_plot.pdf", width = 12, height = 8)
plot(pca$x[,1], pca$x[,2],
     xlab = "PC1", ylab = "PC2",
     main = "PCA of Normalized Expression Data",
     pch = 19, col = "blue")
grid()
dev.off()

# Filter low-intensity probes
row_median <- rowMedians(processed_data)
threshold <- 3.5
filtered_data <- processed_data[row_median > threshold, ]

# Number of transcripts after filtering
n_transcripts_after <- nrow(filtered_data)
cat("Number of transcripts after filtering:", n_transcripts_after, "\n")

# Example: all CEL files are tumor
groups <- factor(rep("disease", ncol(filtered_data)))

# Check relabeled groups
table(groups)
    

#-------------------------------------------------------------------------------
# Molecular Bridging of IBS and MDD Through Gut Dysbiosis
# STEP 3: Differential Expression Analysis
#-------------------------------------------------------------------------------

# Checking working directory
getwd()
[1] "C:/Users/User/Desktop/AI_Omics_Internship_2025"

# Set working directory
setwd("C:/Users/User/Desktop/AI_Omics_Internship_2025/MDD")

#Organizing Project
# File menu > New Project > New Directory > New Project > MDD - Step 3

# Load GEOquery
library(GEOquery)

gse <- getGEO(filename = "GSE201332_series_matrix.txt.gz")

class(gse)

expr <- exprs(gse)
pdata <- pData(gse)

head(expr[, 1:5])
head(pdata)

# Check if the metadata column contains "MDD" and "Control"
unique(pdata$source_name_ch1)

# Create the Group Factor (Control vs MDD)
pdata$source_name_ch1

# Create a clean group vector
group <- ifelse(pdata$source_name_ch1 == "human blood of MDD patients",
                "MDD", "Control")
group <- factor(group)
table(group)


# Prepare Data for Limma (Differential Expression Analysis)
library(limma)

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
design

# Fit Linear Model and Compute DEGs
fit <- lmFit(expr, design)

contrast.matrix <- makeContrasts(MDDvsControl = MDD - Control, levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

deg_results <- topTable(fit2, adjust.method = "fdr", number = Inf)
head(deg_results)

# Save the DEGs
write.csv(deg_results, "DEG_results_MDD_vs_Control.csv")
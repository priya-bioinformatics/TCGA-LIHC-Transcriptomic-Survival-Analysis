# Load packages
library(limma)

# Set working directory
setwd("D:/Bioinformatics project/TCGA-LIHC")

# Load the STAR gene-level expression matrix
counts <- read.table(
  "TCGA-LIHC.star_counts.tsv",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)
# Identify Tumor vs Normal samples from TCGA barcodes
sample_info <- data.frame(
  sample_id = colnames(counts),
  sample_type = ifelse(
    substr(colnames(counts), 14, 15) == "01",
    "Tumor",
    "Normal"
  )
)
# Set sample IDs as row names
rownames(sample_info) <- sample_info$sample_id
# Verify sample distribution
table(sample_info$sample_type)

# Create group factor

group <- factor(sample_info$sample_type,
                levels = c("Normal", "Tumor"))

# Build the design matrix
design <- model.matrix(~ group)
colnames(design)
# Fit linear model to every gene
fit <- lmFit(counts, design)
# Apply empirical Bayes moderation
fit <- eBayes(fit)
# Extract differential expression results
deg_results <- topTable(
  fit,
  coef = "groupTumor",
  number = Inf,
  adjust.method = "BH"
)

head(deg_results)

# Identify significant DEGs (FDR < 0.05)
sig_degs <- deg_results[deg_results$adj.P.Val < 0.05, ]
nrow(sig_degs)
# Apply biological effect size filter
sig_degs_strict <- sig_degs[abs(sig_degs$logFC) > 1, ]
nrow(sig_degs_strict)


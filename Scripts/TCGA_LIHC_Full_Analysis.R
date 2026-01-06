# Load the packages
library(limma)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(ggrepel)
library(data.table)
library(pheatmap)
library(RColorBrewer)
library(edgeR)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(ReactomePA)
library(dplyr)
library(survival)
library(survminer)

# Set Working directory
setwd("D:/Bioinformatics project/TCGA-LIHC")

# Load STAR gene-level count data
expr <- read.table("TCGA-LIHC.star_counts.tsv", 
                   header = TRUE, sep = "\t", 
                   row.names = 1, check.names = FALSE)
head(expr[, 1:5])

# ============================================================
# Differential Expression Analysis (TCGA-LIHC)
# - Identify Tumor vs Normal samples
# - Fit limma linear model with empirical Bayes moderation
# - Identify significant DEGs (FDR < 0.05)
# - Annotate genes (ENSEMBL → SYMBOL / ENTREZ)
# - Visualize results using a volcano plot
# ============================================================

# Identify Tumor vs Normal samples
sample_type <- substr(colnames(expr), 14, 15)

# Assign biological condition: 01 -> Tumor, 11 -> Normal
condition <- ifelse(sample_type == "01", "Tumor",
                    ifelse(sample_type == "11", "Normal", NA))

# Remove any samples that are not Tumor/Normal (Unwanted Samples)
expr <- expr[, !is.na(condition)]
condition <- condition[!is.na(condition)]

# Convert condition to factor
condition <- factor(condition, levels = c("Normal", "Tumor"))
table(condition)

# Create design matrix
design <- model.matrix(~condition)
colnames(design)

# Fit linear model
fit <- lmFit(expr, design)

# Apply empirical Bayes smoothing
fit <- eBayes(fit)
# Get top DE genes (adjusted p-value < 0.05)
top_genes <- topTable(fit, coef="conditionTumor", number=Inf, adjust.method="BH")
head(top_genes)
sig_genes <- top_genes[top_genes$adj.P.Val < 0.05, ]

# Separate up- and down-regulated genes
up_genes <- sig_genes[sig_genes$logFC > 0, ]
down_genes <- sig_genes[sig_genes$logFC < 0, ]
table(significant = sig_genes$logFC > 0)

# Remove ENSEMBL version numbers
ensembl_ids <- sub("\\..*", "", rownames(sig_genes))

# Map ENSEMBL → ENTREZ ID
gene_entrez <- mapIds(org.Hs.eg.db, keys = ensembl_ids,
                      column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")

# Map ENSEMBL → Gene Symbol
gene_symbol <- mapIds(org.Hs.eg.db, keys = ensembl_ids,
                      column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Add annotation columns to DE results
sig_genes$ENTREZID <- gene_entrez
sig_genes$SYMBOL <- gene_symbol

# Remove unmapped genes
sig_genes <- sig_genes[!is.na(sig_genes$ENTREZID), ]

# Define volcano plot thresholds
logFC_cutoff <- 1       # fold change threshold
adjP_cutoff <- 0.05     # adjusted p-value threshold

# Classify each gene (for coloring)
sig_genes$Significance <- "Not Significant"

sig_genes$Significance[
  sig_genes$logFC > logFC_cutoff & sig_genes$adj.P.Val < adjP_cutoff
] <- "Upregulated"

sig_genes$Significance[
  sig_genes$logFC < -logFC_cutoff & sig_genes$adj.P.Val < adjP_cutoff
] <- "Downregulated"

# Select top genes for labeling (most significant)
top_labels <- sig_genes[order(sig_genes$adj.P.Val), ]
top_labels <- rbind(
  head(top_labels[top_labels$Significance=="Upregulated", ], 10),
  head(top_labels[top_labels$Significance=="Downregulated", ], 10)
)


# Volcano plot
ggplot(sig_genes, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Downregulated"="blue", "Upregulated"="red", "Not Significant"="grey")) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype="dashed", color="black") +
  geom_hline(yintercept = -log10(adjP_cutoff), linetype="dashed", color="black") +
  geom_text_repel(data = top_labels, aes(label = SYMBOL), 
                  size = 3, max.overlaps = 15) +
  labs(title = "Volcano Plot: Tumor vs Normal",
       x = "log2 Fold Change",
       y = "-log10(adjusted P-value)",
       color = "Gene Regulation") +
  theme_minimal(base_size = 12)

# ============================================================
# Heatmap Visualization of Top Differentially Expressed Genes
# - Aggregate duplicate ENSEMBL gene IDs
# - Clean expression matrix (remove version numbers)
# - Select top 50 DEGs based on adjusted p-values
# - Annotate samples (Tumor vs Normal) and genes (Up/Down)
# - Generate clustered heatmap with row-wise Z-score scaling
# ============================================================

# Convert expression matrix to data.table
dt <- as.data.table(expr)
head(dt)

# Add ENSEMBL as a column
dt[, ENSEMBL := rownames(expr)]

# Aggregate duplicate ENSEMBL IDs
expr_agg_dt <- dt[, lapply(.SD, mean, na.rm = TRUE), by = ENSEMBL]

# Convert back to matrix with ENSEMBL as rownames
expr_agg_mat <- as.matrix(expr_agg_dt[, -1, with = FALSE])
rownames(expr_agg_mat) <- expr_agg_dt$ENSEMBL

# Quality Check
dim(expr_agg_mat)
head(expr_agg_mat)

# Remove version numbers from ENSEMBL IDs in the matrix
rownames(expr_agg_mat) <- sub("\\..*", "", rownames(expr_agg_mat))

# Take top 20 significant genes
top20 <- sig_genes[order(sig_genes$adj.P.Val), ][1:20, ]

# Remove version numbers from rownames to match matrix
top20_ids <- sub("\\..*", "", rownames(top20))
expr_subset <- expr_agg_mat[top20_ids, , drop = FALSE]

# Replace ENSEMBL IDs with gene symbols
rownames(expr_subset) <- top20$SYMBOL

# Create sample annotation
annotation_col <- data.frame(Condition = condition)
rownames(annotation_col) <- colnames(expr_subset)


# Select top 50 DEGs
top50 <- sig_genes[order(sig_genes$adj.P.Val), ][1:50, ]

# Remove ENSEMBL version numbers
top50_ids <- sub("\\..*", "", rownames(top50))

# Subset expression matrix
expr_subset <- expr_agg_mat[top50_ids, , drop = FALSE]  # aggregated expression matrix

# Replace row names with gene symbols for readability
rownames(expr_subset) <- top50$SYMBOL  

# Column annotation (sample condition)
annotation_col <- data.frame(Condition = condition)
rownames(annotation_col) <- colnames(expr_subset)

# Row annotation (Up/Down)
row_annotation <- data.frame(
  Regulation = ifelse(top50$logFC > 0, "Up", "Down")
)
rownames(row_annotation) <- top50$SYMBOL

# Heatmap colors
heat_colors <- colorRampPalette(c("blue", "white", "red"))(100)  # blue=low, red=high

# Plot heatmap
pheatmap(
  expr_subset,
  scale = "row",               # Z-score per gene
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = annotation_col,
  annotation_row = row_annotation,
  color = heat_colors,
  fontsize_row = 9,
  fontsize_col = 10,
  main = "Top 50 DEGs Heatmap",
  angle_col = 45
)

# ============================================================
# Principal Component Analysis (PCA)
# - Use top 50 differentially expressed genes
# - Log2-CPM normalize expression values
# - Perform PCA to assess Tumor vs Normal separation
# ============================================================
# PCA of Top 50 DEGs

# Prepare Top 50 gene IDs
top50_ids <- sub("\\..*", "", rownames(top50))  # remove Ensembl versions

# Log2 CPM transformation
log_expr_top50 <- cpm(expr_agg_mat[top50_ids, ], log=TRUE, prior.count=1)

# Perform PCA (samples in rows)
pca_res <- prcomp(t(log_expr_top50), scale.=TRUE)

# Calculate % variance explained
pca_var <- pca_res$sdev^2 / sum(pca_res$sdev^2)
percentVar <- round(100 * pca_var[1:2], 1)  # PC1 and PC2

# Build PCA plotting dataframe
pca_df <- data.frame(
  Sample = colnames(log_expr_top50),
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2],
  Condition = condition
)

# PCA plot
ggplot(pca_df, aes(x=PC1, y=PC2, color=Condition)) +
  geom_point(size=4) +
  stat_ellipse(aes(fill=Condition), geom="polygon", alpha=0.2, color=NA) +
  labs(title="PCA of Top 50 DEGs",
       x=paste0("PC1: ", percentVar[1], "% variance"),
       y=paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal(base_size=12) +
  scale_color_manual(values=c("Normal"="blue","Tumor"="red")) +
  scale_fill_manual(values=c("Normal"="blue","Tumor"="red"))

# ============================================================
# Functional Enrichment Analysis
# - Prepare significant DEG gene lists (ENTREZ IDs)
# - Perform GO enrichment (BP, CC, MF)
# - Perform KEGG pathway enrichment
# - Visualize enriched terms using dot plots
# ============================================================
# Prepare all significant gene lists
all_genes <- sig_genes$ENTREZID

# Up- and down-regulated
up_genes <- sig_genes$ENTREZID[sig_genes$logFC > logFC_cutoff & sig_genes$adj.P.Val < adjP_cutoff]
down_genes <- sig_genes$ENTREZID[sig_genes$logFC < -logFC_cutoff & sig_genes$adj.P.Val < adjP_cutoff]

# GO enrichment
ego_all <- enrichGO(
  gene         = all_genes,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENTREZID",
  ont          = "ALL",       # options: BP, CC, MF, or ALL
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable     = TRUE         # converts Entrez ID to gene symbols
)

# View GO enrichment result
head(ego_all)
# Plot top 20 GO terms
dotplot(ego_all, showCategory=20) + ggtitle("GO Enrichment: All DEGs")
# KEGG enrichment
ekegg_all <- enrichKEGG(
  gene         = all_genes,
  organism     = "hsa",      # human
  pvalueCutoff = 0.05
)

# Convert Entrez IDs → Gene symbols
ekegg_all <- setReadable(ekegg_all, OrgDb=org.Hs.eg.db, keyType="ENTREZID")

# View top KEGG pathways
head(ekegg_all)

# Plot top 20 KEGG pathways
dotplot(ekegg_all, showCategory=20) + ggtitle("KEGG Enrichment: All DEGs")

# ============================================================
# GO Term Simplification and Enrichment Map Visualization
# - Remove redundant GO terms based on semantic similarity
# - Retain most significant representative terms
# - Visualize relationships among GO terms using emapplot
# ============================================================
# Remove redundant GO terms (simplify())
ego_simp <- simplify(
  ego_all,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

# Compute term similarity matrix
ego_simp <- pairwise_termsim(ego_simp)

# Enrichment map plot
emapplot(ego_simp, showCategory = 20)

# ============================================================
# Reactome Pathway Enrichment Analysis
# - Prepare ENTREZ gene lists (all, upregulated, downregulated)
# - Perform Reactome enrichment analysis
# - Visualize enriched pathways using barplot and dotplot
# ============================================================
# Upregulated genes
up_entrez <- sig_genes$ENTREZID[sig_genes$logFC > 0]

# Downregulated genes
down_entrez <- sig_genes$ENTREZID[sig_genes$logFC < 0]

# All significant DEGs
all_entrez <- sig_genes$ENTREZID

# Reactome enrichment (ALL DEGs)
reactome_all <- enrichPathway(gene = all_entrez,
                              organism = "human",   # 'human' for TCGA
                              pvalueCutoff = 0.05,
                              readable = TRUE)      # Converts Entrez IDs to gene symbols

# Reactome enrichment (Up & Down separately)
reactome_up <- enrichPathway(gene = up_entrez, organism = "human", pvalueCutoff = 0.05, readable = TRUE)
reactome_down <- enrichPathway(gene = down_entrez, organism = "human", pvalueCutoff = 0.05, readable = TRUE)

# Inspect Reactome results
head(reactome_all)
as.data.frame(reactome_all)

# Bar plot of top enriched pathways
barplot(reactome_all, showCategory=20, title="Reactome Pathway Enrichment")

# Optional: Dot plot of top enriched pathways
dotplot(reactome_all, showCategory=20, title="Reactome Pathway Enrichment")
# ============================================================
# Reactome Gene Set Enrichment Analysis (GSEA)
# - Create ranked gene list using log2 fold change
# - Remove unmapped and duplicated ENTREZ IDs
# - Perform Reactome GSEA
# - Visualize enriched pathways and enrichment curves
# ============================================================
# Create a ranked gene list
gene_list <- sig_genes$logFC
names(gene_list) <- sig_genes$ENTREZID

# Remove genes without Entrez IDs
gene_list <- gene_list[!is.na(names(gene_list))]

# If duplicated ENTREZ IDs exist, keep the max logFC
gene_list <- tapply(gene_list, names(gene_list), max)

# Sort decreasing (Rank genes)
gene_list <- sort(gene_list, decreasing = TRUE)

# Check ranked gene list
head(gene_list)



# (Safety step) Ensure numeric + sorted
gene_list <- setNames(as.numeric(gene_list), names(gene_list))
gene_list <- sort(gene_list, decreasing = TRUE)

# Run Reactome GSEA
gsea_reactome <- gsePathway(
  geneList     = gene_list,
  organism     = "human",
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

# Check results
head(gsea_reactome)

# Dot plot
library(enrichplot)
dotplot(gsea_reactome, showCategory = 20) +
  ggtitle("Reactome GSEA – TCGA LIHC")

# GSEA score plot
gseaplot2(
  gsea_reactome,
  geneSetID = 1,
  title = gsea_reactome$Description[1]
)

# Save GSEA results to CSV
write.csv(
  as.data.frame(gsea_reactome),
  file = "Reactome_GSEA_TCGA_LIHC.csv",
  row.names = FALSE
)

# ============================================================
# Survival Analysis of Candidate Genes (TCGA-LIHC)
# - Match expression samples with patient survival data
# - Prioritize strong DEGs supported by GO/Reactome enrichment
# - Perform Kaplan–Meier survival analysis (High vs Low expression)
# - Fit Cox proportional hazards models
# - Identify and validate prognostic biomarkers
# ============================================================
#Survival analysis
# Load survival data
surv <- read.table(
  "TCGA-LIHC.survival.tsv",   # filename may vary slightly
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)
head(surv)
# Match expression samples with survival patients
expr_patients <- substr(colnames(expr), 1, 12)
common_patients <- intersect(expr_patients, surv$X_PATIENT)
length(common_patients)

# Subset expression matrix
expr_surv <- expr[, expr_patients %in% common_patients]

# Subset survival data
surv_sub <- surv[surv$X_PATIENT %in% common_patients, ]

# Reorder survival rows to match expression columns
surv_sub <- surv_sub[match(substr(colnames(expr_surv), 1, 12), surv_sub$X_PATIENT), ]

# Sanity check
all(substr(colnames(expr_surv), 1, 12) == surv_sub$X_PATIENT)

# Select strong DEGs for candidate genes
# Select one gene which is meeting all 3 criterias
strong_deg <- sig_genes[
  abs(sig_genes$logFC) > 2 &
    sig_genes$adj.P.Val < 1e-10,
]
# Count how many strong DEGs
nrow(strong_deg)

# View top candidate genes by effect size
head(
  strong_deg[order(abs(strong_deg$logFC), decreasing = TRUE),
             c("SYMBOL","logFC","adj.P.Val","ENTREZID")],
  20
)

# Extract genes from GO enrichment results
go_genes <- unique(unlist(
  strsplit(ego_all@result$geneID, "/")
))
head(go_genes)

# Extract genes from Reactome GSEA results
reactome_genes <- unique(unlist(
  strsplit(gsea_reactome@result$core_enrichment, "/")
))
head(reactome_genes)
# Find overlap
# Keep strong DEGs present in GO or Reactome
final_candidates <- strong_deg[
  strong_deg$ENTREZID %in% c(go_genes, reactome_genes),
]
# Sort by strongest effect
final_candidates <- final_candidates[
  order(abs(final_candidates$logFC), decreasing = TRUE),
]
# View the top candidate genes
head(final_candidates[, c("SYMBOL","logFC","adj.P.Val","ENTREZID")], 10)
# Select the top-ranked gene
final_gene <- final_candidates[1, ]
final_gene
# Get final gene identifiers
gene_symbol <- final_gene$SYMBOL
gene_symbol        # Selected CYP1A2 as your final candidate

# Get Ensembl ID for CYP1A2 gene
gene_ensembl <- rownames(sig_genes)[sig_genes$SYMBOL == gene_symbol]
gene_ensembl <- sub("\\..*", "", gene_ensembl)[1]
gene_ensembl

# Convert expression matrix to data.table
dt <- as.data.table(expr, keep.rownames = "ENSEMBL")

# Remove Ensembl version numbers
dt[, ENSEMBL := sub("\\..*", "", ENSEMBL)]

# Aggregate duplicated Ensembl IDs by mean
expr_agg_dt <- dt[, lapply(.SD, mean, na.rm = TRUE), by = ENSEMBL]
# Convert back to matrix
expr_agg_mat <- as.matrix(expr_agg_dt[, -1, with = FALSE])
rownames(expr_agg_mat) <- expr_agg_dt$ENSEMBL

# Check aggregated expression matrix
dim(expr_agg_mat)
head(expr_agg_mat)

# Extract expression of “CYP1A2”
gene_expr <- expr_agg_mat[gene_ensembl, , drop = TRUE]

# Check
head(gene_expr)

# Extract patient IDs from expression data
expr_patients <- substr(colnames(expr_agg_mat), 1, 12)
common_patients <- intersect(expr_patients, surv$X_PATIENT)

# Subset gene expression to matched patients
gene_expr <- gene_expr[substr(names(gene_expr), 1, 12) %in% common_patients]
# Subset survival data
surv_sub <- surv[surv$X_PATIENT %in% common_patients, ]

# Reorder survival table to match expression
surv_sub <- surv_sub[match(substr(names(gene_expr), 1, 12), surv_sub$X_PATIENT), ]
# Verify alignment (CRITICAL check)
all(substr(names(gene_expr), 1, 12) == surv_sub$X_PATIENT)  # should be TRUE
# Median split of gene expression
expr_group <- ifelse(gene_expr > median(gene_expr, na.rm = TRUE), "High", "Low")
# Check group sizes
table(expr_group)
# Create survival object
surv_obj <- Surv(time = surv_sub$OS.time, event = surv_sub$OS)

# Fit Kaplan-Meier model
fit <- survfit(surv_obj ~ expr_group)

# Plot
plot_data <- data.frame(
  expr_group = expr_group,
  OS.time = surv_sub$OS.time,
  OS = surv_sub$OS
)
ggsurvplot(
  fit,
  data = plot_data,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = TRUE,
  palette = c("blue", "red"),
  legend.title = gene_symbol,
  legend.labs = c("High", "Low"),
  title = paste("Overall Survival by", gene_symbol, "Expression")
)
# Cox model
cox_fit <- coxph(surv_obj ~ expr_group)
# View Cox model results
summary(cox_fit)
# Get Cox summary
cox_sum <- summary(cox_fit)

# Convert Cox model results into a data frame
cox_df <- data.frame(
  Variable = rownames(cox_sum$coefficients),
  coef = cox_sum$coefficients[, "coef"],
  exp_coef = cox_sum$coefficients[, "exp(coef)"],
  se_coef = cox_sum$coefficients[, "se(coef)"],
  z = cox_sum$coefficients[, "z"],
  p_value = cox_sum$coefficients[, "Pr(>|z|)"],
  exp_neg_coef = cox_sum$conf.int[, "exp(-coef)"],
  lower95 = cox_sum$conf.int[, "lower .95"],
  upper95 = cox_sum$conf.int[, "upper .95"],
  concordance = cox_sum$concordance[1],
  concordance_se = cox_sum$concordance[2]
)
# Save to CSV
write.csv(cox_df, "Cox_results_gene.csv", row.names = FALSE)

# Looping through the top 10 candidate genes and for each gene to do KM + Cox
top_genes <- head(final_candidates, 10)  # test top 10 strong DEGs
results <- data.frame()
for(i in 1:nrow(top_genes)) {
  gene_sym <- top_genes$SYMBOL[i]
  gene_ensembl <- sub("\\..*", "", rownames(sig_genes)[sig_genes$SYMBOL == gene_sym][1])
  
  # Make sure gene exists in aggregated expression matrix
  if(!gene_ensembl %in% rownames(expr_agg_mat)) next
  # Extract gene expression
  gene_expr <- expr_agg_mat[gene_ensembl, , drop = TRUE]
  
  # Match patients with survival data
  expr_patients <- substr(names(gene_expr), 1, 12)
  common_patients <- intersect(expr_patients, surv$X_PATIENT)
  gene_expr <- gene_expr[substr(names(gene_expr), 1, 12) %in% common_patients]
  surv_sub <- surv[surv$X_PATIENT %in% common_patients, ]
  surv_sub <- surv_sub[match(substr(names(gene_expr), 1, 12), surv_sub$X_PATIENT), ]
  
  # Median split
  expr_group <- ifelse(gene_expr > median(gene_expr, na.rm = TRUE), "High", "Low")
  
  # Cox model
  cox_fit <- coxph(Surv(OS.time, OS) ~ expr_group, data = surv_sub)
  cox_sum <- summary(cox_fit)
  
  # Record results
  results <- rbind(results, data.frame(
    SYMBOL = gene_sym,
    coef = cox_sum$coefficients[1, "coef"],
    HR = cox_sum$coefficients[1, "exp(coef)"],
    p_value = cox_sum$coefficients[1, "Pr(>|z|)"],
    lower95 = cox_sum$conf.int[1, "lower .95"],
    upper95 = cox_sum$conf.int[1, "upper .95"]
  ))
}

# Sort by lowest p-value
results <- results[order(results$p_value), ]
results

# STAB2- chosen for Survival analysis
# p = 0.0457 (significant at 0.05)
# HR = 0.732 → patients with High STAB2 expression have lower hazard, so potentially protective
# 95% CI = 0.539–0.994 → barely excludes 1, but still acceptable
# Extract STAB2 expression from the aggregated expression matrix
gene_expr <- expr_agg_mat[sub("\\..*", "", rownames(sig_genes)[sig_genes$SYMBOL=="STAB2"]), ]
# Subset only patients with survival data
expr_patients <- substr(names(gene_expr), 1, 12)
common_patients <- intersect(expr_patients, surv_sub$X_PATIENT)
gene_expr <- gene_expr[expr_patients %in% common_patients]
# Reorder survival table to match expression
surv_sub_matched <- surv_sub[match(substr(names(gene_expr), 1, 12), surv_sub$X_PATIENT), ]
# Sanity check to ensure patient IDs align
all(substr(names(gene_expr), 1, 12) == surv_sub_matched$X_PATIENT)  # should be TRUE
# Median split of expression into High vs Low groups for Kaplan-Meier
expr_group <- ifelse(gene_expr > median(gene_expr, na.rm = TRUE), "High", "Low")
table(expr_group)  # Check how many patients fall into each group
# Create survival object
surv_obj <- Surv(time = surv_sub_matched$OS.time, event = surv_sub_matched$OS)
# Kaplan-Meier fit
fit <- survfit(surv_obj ~ expr_group)
#  Plot KM curve with risk table and confidence interval
km_plot <- ggsurvplot(
  fit,
  data = data.frame(expr_group = expr_group),
  pval = TRUE,
  risk.table = TRUE,
  conf.int = TRUE,
  palette = c("blue", "red"),
  title = "Survival by STAB2 expression",
  risk.table.height = 0.25
)
print(km_plot)

# Cox Model
cox_fit <- coxph(surv_obj ~ expr_group)
summary(cox_fit)
# Convert Cox summary to a data frame for saving/reporting
cox_sum <- summary(cox_fit)
cox_df <- data.frame(
  Variable = rownames(cox_sum$coefficients),
  coef = cox_sum$coefficients[, "coef"],
  HR = cox_sum$coefficients[, "exp(coef)"],
  se_coef = cox_sum$coefficients[, "se(coef)"],
  z = cox_sum$coefficients[, "z"],
  p_value = cox_sum$coefficients[, "Pr(>|z|)"],
  lower95 = cox_sum$conf.int[, "lower .95"],
  upper95 = cox_sum$conf.int[, "upper .95"],
  concordance = cox_sum$concordance[1],
  concordance_se = cox_sum$concordance[2]
)
# Save Cox results to CSV
write.csv(cox_df, "Cox_results_STAB2.csv", row.names = FALSE)


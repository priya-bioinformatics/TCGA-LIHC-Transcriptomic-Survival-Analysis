# TCGA-LIHC-Transcriptomic-Survival-Analysis
**Project Summary**
This project analyzes TCGA-LIHC STAR-aligned RNA-seq data to identify differentially expressed genes, interpret dysregulated biological pathways, and evaluate prognostic biomarkers using survival analysis. The focus is on combining statistical rigor, biological interpretation, and clinical relevance.

**Dataset**
424 samples (371 tumor, 53 normal)
STAR-aligned, log₂-normalized gene-level expression

**Tools**
R, limma, edgeR, clusterProfiler, ReactomePA, enrichplot, survminer

**Analysis Overview**
Differential Expression
Tumor vs Normal comparison using limma
DEGs defined as FDR < 0.05 and |log₂FC| > 1
Clear tumor-normal separation observed in volcano plots, heatmaps, and PCA

Functional Enrichment
GO, KEGG, and Reactome pathway analysis
GSEA performed on ranked gene lists
Enrichment highlights cell-cycle activation, genomic instability, metabolic reprogramming, and immune dysregulation in LIHC

Candidate Gene Prioritization
Strong DEGs filtered using stringent thresholds (|log₂FC| > 2, FDR < 1e-10)
Further refined using GO and Reactome pathway support
Biologically relevant candidates identified (e.g., CYP1A2)

Survival Analysis
Expression data integrated with TCGA clinical survival information
Kaplan–Meier and Cox proportional hazards models applied
Systematic survival screening identified STAB2 as a protective prognostic marker
(HR = 0.73, p = 0.046; high expression associated with improved survival)

**Key Take-Home Message**
This analysis integrates differential expression, pathway enrichment, and survival modeling to identify biologically meaningful and clinically relevant molecular features of hepatocellular carcinoma.

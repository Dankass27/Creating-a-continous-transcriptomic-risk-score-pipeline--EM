################################################################################
# This is my code from September 2025 - May 2026
# MEDULLOBLASTOMA SURVIVAL ANALYSIS
# Dataset: GSE85217 (Taylor Lab, 763 pediatric MB samples)
#
# Student: Ewan McGibbon
# Date:   7/4/2026
#
# Description:
#   Analysis of gene expression and clinical data to identify prognostic
#   factors in pediatric medulloblastoma. Includes dimensionality reduction
#   (PCA/UMAP), Kaplan-Meier survival analysis, genome-wide Cox regression,
#   and penalized Cox models for risk score development.
#
################################################################################

# =============================================================================
# 0. SETUP
# =============================================================================

# -----------------------------------------------------------------------------
# 0.1 Install packages 
# -----------------------------------------------------------------------------

# install.packages("data.table")
# install.packages("ggplot2")
# install.packages("uwot")
# install.packages("patchwork")
# install.packages("survival")
# install.packages("survminer")
# install.packages("glmnet")
# BiocManager::install("biomaRt")

# -----------------------------------------------------------------------------
# 0.2 Load libraries
# -----------------------------------------------------------------------------

library(data.table)
library(ggplot2)
library(uwot)
library(biomaRt)
library(survival)
library(survminer)
library(patchwork)
library(glmnet)
library(org.Hs.eg.db)
library(pbapply)
library(ComplexHeatmap)
library(circlize)
library(ggridges)
library(ggbeeswarm)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(flextable)
library(officer)
library(gridExtra)


########################################################################################################################
########################################################################################################################
#
#                                   1. BUILD THE ANALYSIS MATRIX (Steps 1-3)
#
########################################################################################################################
########################################################################################################################
#
# OVERVIEW:
# The raw GEO expression file contains ~21,000 genes across 763 MB samples.
# Before running PCA/UMAP, we need to clean and filter this data:
#
#   - Remove sex-chromosome genes (X/Y) to prevent sex-driven clustering
#   - Select the most variable genes to focus on biologically meaningful signal
#     and reduce noise from lowly-expressed or invariant genes
#
# STEPS:
#   Step 1: Load expression data → expr_mat (genes x samples)
#   Step 2: Remove X/Y genes     → expr_noXY
#   Step 3: Select top 2000 variable genes → expr_top
#
# OUTPUT: expr_top (2000 genes x 763 samples)
#
########################################################################################################################

########################################
# STEP 1: Load expression data         #
########################################

# Aim: Load the raw GEO expression table and convert it into a clean numeric
# matrix that we can use for PCA/UMAP.
#
# Input:  Text file where rows = genes, columns = annotation + sample values
# Output: expr_mat (genes x samples), numeric, with Ensembl IDs as row names

# 1.1 Read in the expression data file
expr <- fread("GSE85217_M_exp_763_MB_SubtypeStudy_TaylorLab.txt")

# 1.2 Check dimensions and structure
dim(expr)
head(expr)
names(expr)[1:10]

# 1.3 Convert from data.table to data.frame
expr_df <- as.data.frame(expr)

# 1.4 Extract the expression matrix (first 5 columns are annotations)
expr_mat <- as.matrix(expr_df[, 6:ncol(expr_df)])

# 1.5 Set row names to Ensembl gene IDs
rownames(expr_mat) <- expr_df$EnsemblGeneID_from_ensemblv77

# 1.6 Verify the matrix looks right
dim(expr_mat)
expr_mat[1:5, 1:5]
##############################################
# STEP 2: Remove sex-chromosome genes (X/Y)  #
##############################################

# 2.1 Connect to Ensembl with mirror fallback
mirrors <- c("useast", "uswest", "asia")
ensembl <- NULL

for (mirror in mirrors) {
  ensembl <- tryCatch({
    useEnsembl(
      biomart = "genes",
      dataset = "hsapiens_gene_ensembl",
      mirror  = mirror
    )
  }, error = function(e) {
    message("Mirror '", mirror, "' failed: ", conditionMessage(e))
    NULL
  })
  if (!is.null(ensembl)) {
    message("Connected via mirror: ", mirror)
    break
  }
}

if (is.null(ensembl)) stop("All Ensembl mirrors failed. Try again later.")

# 2.2 Query chromosome locations for all genes in expr_mat
annot <- getBM(
  attributes = c("ensembl_gene_id", "chromosome_name"),
  filters    = "ensembl_gene_id",
  values     = rownames(expr_mat),
  mart       = ensembl
)

# 2.3 Check
message(nrow(annot), " genes annotated")

# 2.4 Keep only autosomal genes (chromosomes 1–22)
annot_noXY <- annot[annot$chromosome_name %in% as.character(1:22), ]
ens_keep   <- unique(annot_noXY$ensembl_gene_id)
message(length(ens_keep), " autosomal genes retained")

# 2.5 Filter expression matrix
expr_noXY <- expr_mat[rownames(expr_mat) %in% ens_keep, ]
dim(expr_noXY)

#################################################
# STEP 3: Select top 2000 most variable genes   #
#################################################

# Aim: Keep genes with highest variation across samples. These carry the
# strongest signal for separating samples in PCA/UMAP and reduce noise.

# 3.1 Calculate standard deviation for each gene
gene_sd <- apply(expr_noXY, 1, sd, na.rm = TRUE)

# 3.2 Rank genes from most to least variable
order_sd <- order(gene_sd, decreasing = TRUE)

# 3.3 Pick the top 2000 most variable genes
top2000_idx <- order_sd[1:2000]

# 3.4 Subset the expression matrix
expr_top <- expr_noXY[top2000_idx, ]
dim(expr_top)


########################################################################################################################
########################################################################################################################
#
#                                   2. PCA + UMAP WORKFLOW (Steps 4-8)
#
########################################################################################################################
########################################################################################################################
#
# OVERVIEW:
# Dimensionality reduction compresses the 2000-gene expression space into 2D
# so we can visualise how samples relate to each other. We use two methods:
#
#   - PCA: Linear method that captures global variance structure
#   - UMAP: Non-linear method that preserves local neighbourhood structure
#
# We first run both methods WITHOUT labels to check for natural clustering,
# then add molecular subgroup labels to see if the clusters match biology.
#
# STEPS:
#   Step 4: Unsupervised PCA and UMAP (no labels yet)
#   Step 5: Load clinical metadata
#   Step 6: Align metadata rows to expression columns
#   Step 7: PCA coloured by subgroup
#   Step 8: UMAP coloured by subgroup
#   → Figure 1: Two-panel PCA + UMAP with consensus colours
#
########################################################################################################################

###########################################
# STEP 4: Unsupervised PCA and UMAP       #
###########################################

# Aim: Compress the 2000-gene space down to 2D to visualise sample structure
# BEFORE adding any subgroup labels. This lets us check for natural clustering.

set.seed(123)

# --- 4.1 PCA ---

# PCA expects rows = samples, columns = genes, so we transpose expr_top
pca_res <- prcomp(t(expr_top), scale. = TRUE)

# Extract PC1 and PC2 for plotting
pca_df <- data.frame(
  Sample = rownames(pca_res$x),
  PC1    = pca_res$x[, 1],
  PC2    = pca_res$x[, 2]
)

# Plot (no colours yet)
ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(size = 2) +
  theme_classic() +
  labs(title = "PCA - unlabelled (top 2000 autosomal genes)")


# --- 4.2 UMAP ---

# UMAP is stochastic, so the seed ensures reproducible layouts
umap_res <- umap(
  t(expr_top),
  n_neighbors = 15,
  min_dist    = 0.3,
  metric      = "euclidean"
)

dim(umap_res)

# Build plotting data frame
umap_df <- data.frame(
  Sample = rownames(pca_res$x),
  UMAP1  = umap_res[, 1],
  UMAP2  = umap_res[, 2]
)

# Plot (no colours yet)
ggplot(umap_df, aes(UMAP1, UMAP2)) +
  geom_point(size = 2) +
  theme_classic() +
  labs(title = "UMAP - unlabelled (top 2000 autosomal genes)")


##########################################
# STEP 5: Load metadata                  #
##########################################

# Aim: Load sample-level metadata (subgroup labels, survival info, etc.)
# so we can colour plots and do downstream survival analysis.

meta <- read.csv("Taylor_Pheno_forDan.csv", stringsAsFactors = FALSE)

# Check what we're working with
dim(meta)
head(meta)
names(meta)

# 5.1 Identify which column contains sample IDs
# (Find the column with most overlap to expression sample names)
overlaps <- sapply(meta, function(col) length(intersect(col, colnames(expr_top))))
overlaps

sample_col_name <- names(overlaps)[which.max(overlaps)]
sample_col_name


###############################################
# STEP 6: Align expression data with metadata #
###############################################

# Aim: Make sure expression columns and metadata rows refer to the same
# samples in the same order. This is critical for all downstream analysis.

# 6.1 Find samples present in BOTH datasets
common_samples <- intersect(colnames(expr_top), meta[[sample_col_name]])
length(common_samples)

# 6.2 Subset expression matrix to matched samples
expr_top_sub <- expr_top[, common_samples]

# 6.3 Subset and reorder metadata to match expression column order
meta_sub <- meta[match(common_samples, meta[[sample_col_name]]), ]

# 6.4 Verify alignment (this MUST return TRUE)
all(colnames(expr_top_sub) == meta_sub[[sample_col_name]])


#############################################
# STEP 7: PCA coloured by subgroup          #
#############################################

# Aim: Rerun PCA on the aligned data and colour by molecular subgroup
# to see whether subgroup structure is captured in the main PCs.

set.seed(123)

pca_res2 <- prcomp(t(expr_top_sub), scale. = TRUE)

pca_df2 <- data.frame(
  Sample   = rownames(pca_res2$x),
  PC1      = pca_res2$x[, 1],
  PC2      = pca_res2$x[, 2],
  Subgroup = meta_sub$Subgroup
)

head(pca_df2)

# Define consensus subgroup colours (used throughout the analysis)
subgroup_cols <- c(
  "WNT"    = "#0000FF",
  "SHH"    = "#FF0000",
  "Group3" = "#FFFF00",
  "Group4" = "#008000"
)

# Plot PCA with subgroup colours
ggplot(pca_df2, aes(PC1, PC2, colour = Subgroup)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_colour_manual(values = subgroup_cols) +
  theme_classic() +
  labs(title = "PCA (top 2000 autosomal genes)")


#############################################
# STEP 8: UMAP coloured by subgroup         #
#############################################

# Aim: Run UMAP on the aligned data and colour by subgroup to check
# whether molecular subgroups form distinct clusters.

set.seed(123)

umap_res2 <- umap(
  t(expr_top_sub),
  n_neighbors = 15,
  min_dist    = 0.3,
  metric      = "euclidean"
)

umap_df2 <- data.frame(
  Sample   = colnames(expr_top_sub),
  UMAP1    = umap_res2[, 1],
  UMAP2    = umap_res2[, 2],
  Subgroup = meta_sub$Subgroup
)

head(umap_df2)

# Plot UMAP with subgroup colours
ggplot(umap_df2, aes(UMAP1, UMAP2, colour = Subgroup)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_colour_manual(values = subgroup_cols) +
  theme_classic() +
  labs(title = "UMAP (top 2000 autosomal genes)")


########################################################################################
# Figure 1: Two-panel PCA + UMAP                                                       #
########################################################################################

# Aim: Combine PCA and UMAP into a publication-ready two-panel figure
# with consistent styling and a shared legend.

# Panel A: PCA
p_pca <- ggplot(pca_df2, aes(PC1, PC2, colour = Subgroup)) +
  geom_point(size = 1.8, alpha = 0.85) +
  scale_colour_manual(values = subgroup_cols) +
  labs(
    title = "PCA\n(top 2000 autosomal genes)",
    x = "PC1",
    y = "PC2",
    tag = "A"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.margin = margin(6, 18, 6, 6)
  )

# Panel B: UMAP
p_umap <- ggplot(umap_df2, aes(UMAP1, UMAP2, colour = Subgroup)) +
  geom_point(size = 1.8, alpha = 0.85) +
  scale_colour_manual(values = subgroup_cols) +
  labs(
    title = "UMAP\n(top 2000 autosomal genes)",
    x = "UMAP1",
    y = "UMAP2",
    tag = "B"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.margin = margin(6, 6, 6, 18)
  )

# Combine panels with shared legend
# fig1 <- (p_pca | p_umap) +
#   plot_layout(guides = "collect") &
#   theme(legend.position = "right")

# print(fig1)

# ggsave("Figure1_PCA_UMAP.pdf", fig1, width = 10, height = 5, dpi = 300)

########################################################################################################################
########################################################################################################################
#
#                                   3. KAPLAN-MEIER SURVIVAL ANALYSIS
#
########################################################################################################################
########################################################################################################################
#
# OVERVIEW:
# Kaplan-Meier curves visualise survival probability over time for different
# patient groups. We compare overall survival across four clinical variables:
#
#   - Molecular subgroup (WNT, SHH, Group3, Group4)
#   - Metastasis status (M0 vs Metastatic)
#   - Histology (Classic, Desmoplastic, LCA, MBEN)
#   - Age group (0-3, 4-10, 10+ years)
#
# Each plot includes a log-rank p-value testing whether survival differs
# between groups, plus a risk table showing patients at risk over time.
#
# STEPS:
#   3.1: Create survival object
#   3.2: KM by molecular subgroup
#   3.3: KM by metastasis status
#   3.4: KM by histology
#   3.5: KM by age group
#   3.6: Four-panel Figure 2
#
# OUTPUT: Individual KM plots + Figure 2 (four-panel summary)
#
########################################################################################################################

###############################################################################
# 3.1 Create survival object
###############################################################################

# Define which columns contain survival time and event status
time_var   <- "OS..years."
status_var <- "Dead"

# Create survival object (used for all KM fits)
# Event: 1 = dead, 0 = alive/censored
surv_obj <- Surv(
  time  = meta_sub[[time_var]],
  event = meta_sub[[status_var]]
)


###############################################################################
# 3.2 KM plot: Molecular subgroup
###############################################################################

# Aim: Compare survival across WNT, SHH, Group3, Group4

table(meta_sub$Subgroup)

fit_subgroup <- survfit(surv_obj ~ Subgroup, data = meta_sub)

ggsurvplot(
  fit_subgroup,
  data         = meta_sub,
  risk.table   = TRUE,
  pval         = TRUE,
  conf.int     = FALSE,
  palette      = subgroup_cols,
  legend.title = "Subgroup",
  xlab         = "Time (years)",
  ylab         = "Overall survival probability",
  title        = "Kaplan-Meier curves by molecular subgroup",
  xlim         = c(0, 10)
)


###############################################################################
# 3.3 KM plot: Metastasis status
###############################################################################

# Aim: Compare survival between M0 (non-metastatic) and metastatic patients

table(meta_sub$Met.status..1.Met..0.M0., useNA = "ifany")

# Recode to readable factor
meta_sub$Metastasis_status <- factor(
  meta_sub$Met.status..1.Met..0.M0.,
  levels = c("0", "1"),
  labels = c("M0 / non-metastatic", "Metastatic")
)

table(meta_sub$Metastasis_status, useNA = "ifany")

fit_met <- survfit(surv_obj ~ Metastasis_status, data = meta_sub)

ggsurvplot(
  fit_met,
  data         = meta_sub,
  risk.table   = TRUE,
  pval         = TRUE,
  conf.int     = FALSE,
  legend.title = "Metastasis",
  legend.labs  = levels(meta_sub$Metastasis_status),
  xlab         = "Time (years)",
  ylab         = "Overall survival probability",
  title        = "Kaplan-Meier curves by metastasis status",
  xlim         = c(0, 10)
)


###############################################################################
# 3.4 KM plot: Histology
###############################################################################

# Aim: Compare survival across histological subtypes
# Factor levels: MBEN (reference), Classic, Desmoplastic, LCA
# Palette and legend.labs must match factor level order

table(meta_sub$histology)

meta_sub$histology <- relevel(factor(meta_sub$histology), ref = "MBEN")

fit_hist <- survfit(surv_obj ~ histology, data = meta_sub)

ggsurvplot(
  fit_hist,
  data         = meta_sub,
  risk.table   = TRUE,
  pval         = TRUE,
  conf.int     = FALSE,
  legend.title = "Histology",
  legend.labs  = c("MBEN","Classic","Desmoplastic","LCA"),
  palette      = c("#c51b7d","#878787","#762a83","#1b7837"),
  xlab         = "Time (years)",
  ylab         = "Overall survival probability",
  title        = "Kaplan-Meier curves by histological subtype",
  xlim         = c(0, 10)
)


###############################################################################
# 3.5 KM plot: Age group
###############################################################################

# Aim: Compare survival across age bands (0-3, 4-10, 10+ years)
# AgeGroup built from raw Age variable — previous coding missed 10+ patients

meta_sub$AgeGroup <- case_when(
  as.numeric(meta_sub$Age) < 3                                   ~ "0-3",
  as.numeric(meta_sub$Age) >= 3 & as.numeric(meta_sub$Age) < 10 ~ "4-10",
  as.numeric(meta_sub$Age) >= 10                                 ~ "10+",
  TRUE                                                            ~ NA_character_
)
meta_sub$AgeGroup <- factor(meta_sub$AgeGroup, levels = c("0-3","4-10","10+"))

table(meta_sub$AgeGroup, useNA = "always")

fit_age <- survfit(surv_obj ~ AgeGroup, data = meta_sub)

ggsurvplot(
  fit_age,
  data         = meta_sub,
  risk.table   = TRUE,
  pval         = TRUE,
  conf.int     = FALSE,
  legend.title = "Age group (years)",
  legend.labs  = c("0-3 years","4-10 years","10+ years"),
  xlab         = "Time (years)",
  ylab         = "Overall survival probability",
  title        = "Kaplan-Meier curves by age group",
  xlim         = c(0, 10)
)


###############################################################################
# 3.6 Figure 2: Four-panel KM summary
###############################################################################

library(gridExtra)

# Consensus subgroup colours — must match factor levels exactly
subgroup_cols <- c(
  "WNT"    = "#0000FF",
  "SHH"    = "#FF0000",
  "Group3" = "#FFFF00",
  "Group4" = "#008000"
)

# Set factor levels
meta_sub$Subgroup <- factor(
  meta_sub$Subgroup,
  levels = c("WNT", "SHH", "Group3", "Group4")
)

raw_met <- as.character(meta_sub$Met.status..1.Met..0.M0.)
raw_met[raw_met %in% c("0","0.0","0 ","0.00")] <- "M0"
raw_met[raw_met %in% c("1","1.0","1 ","1.00")] <- "Met"
meta_sub$Met_status <- relevel(factor(raw_met), ref = "M0")

# Refit models
fit_subgroup <- survfit(Surv(OS..years., Dead) ~ Subgroup,          data = meta_sub)
fit_met      <- survfit(Surv(OS..years., Dead) ~ Metastasis_status, data = meta_sub)
fit_hist     <- survfit(Surv(OS..years., Dead) ~ histology,         data = meta_sub)
fit_age      <- survfit(Surv(OS..years., Dead) ~ AgeGroup,          data = meta_sub)

# Panel A: Subgroup
km_subgroup <- ggsurvplot(
  fit_subgroup,
  data          = meta_sub,
  risk.table    = FALSE,
  pval          = TRUE,
  conf.int      = FALSE,
  legend.title  = "",
  legend.labs   = c("WNT", "SHH", "Group3", "Group4"),
  legend        = "top",
  palette       = subgroup_cols,
  xlab          = "Time (years)",
  ylab          = "Survival probability",
  title         = "Subgroup",
  xlim          = c(0, 10),
  break.time.by = 2,
  ggtheme       = theme_classic(base_size = 11)
)

# Panel B: Metastasis
km_met <- ggsurvplot(
  fit_met,
  data          = meta_sub,
  risk.table    = FALSE,
  pval          = TRUE,
  conf.int      = FALSE,
  legend.title  = "",
  legend.labs   = c("Non-metastatic", "Metastatic"),
  legend        = "top",
  xlab          = "Time (years)",
  ylab          = "Survival probability",
  title         = "Metastasis",
  xlim          = c(0, 10),
  break.time.by = 2,
  ggtheme       = theme_classic(base_size = 11)
)

# Panel C: Histology
# legend.labs and palette match factor levels: MBEN, Classic, Desmoplastic, LCA
km_hist <- ggsurvplot(
  fit_hist,
  data          = meta_sub,
  risk.table    = FALSE,
  pval          = TRUE,
  conf.int      = FALSE,
  legend.title  = "",
  legend.labs   = c("MBEN", "Classic", "Desmoplastic", "LCA"),
  legend        = "top",
  palette       = c("#c51b7d", "#878787", "#762a83", "#1b7837"),
  xlab          = "Time (years)",
  ylab          = "Survival probability",
  title         = "Histology",
  xlim          = c(0, 10),
  break.time.by = 2,
  ggtheme       = theme_classic(base_size = 11)
)

# Panel D: Age group
km_age <- ggsurvplot(
  fit_age,
  data          = meta_sub,
  risk.table    = FALSE,
  pval          = TRUE,
  conf.int      = FALSE,
  legend.title  = "",
  legend.labs   = c("0-3 years", "4-10 years", "10+ years"),
  legend        = "top",
  xlab          = "Time (years)",
  ylab          = "Survival probability",
  title         = "Age group",
  xlim          = c(0, 10),
  break.time.by = 2,
  ggtheme       = theme_classic(base_size = 11)
)

# Combine into 2x2 panel — Figure 2
fig2 <- arrangeGrob(
  km_subgroup$plot, km_met$plot,
  km_hist$plot,     km_age$plot,
  ncol = 2, nrow = 2
)
grid::grid.draw(fig2)
########################################################################################################################
########################################################################################################################
#
#                                   4. COX REGRESSION: CLINICAL PREDICTORS
#
########################################################################################################################
########################################################################################################################
#
# OVERVIEW:
# Quantify survival differences using Cox proportional hazards regression.
# Fit univariable Cox models for key clinical predictors across the full cohort.
# These models feed directly into Tables 3 and 4 in MB_tables.R.
#
# Reference categories:
# - Histology: MBEN (best prognosis)
# - Subgroup: Group3 (intermediate prognosis)
# - Metastasis: M0 (non-metastatic)
# - Age group: 4-10 years (better prognosis)
#
# OUTPUT:
#   - cox_hist, cox_subgroup, cox_met, cox_age (full cohort)
#   - cox_hist_g34, cox_age_g34, cox_met_g34  (G3/4 subgroup)
#   - cox_hist_shh, cox_age_shh, cox_met_shh  (SHH subgroup)
#
########################################################################################################################

#########################################
# (A) Prepare clinical covariates       #
#########################################

# Histology: set MBEN as reference
meta_sub$histology <- factor(meta_sub$histology)
meta_sub$histology <- relevel(meta_sub$histology, ref = "MBEN")

# Subgroup: set Group3 as reference
meta_sub$Subgroup <- factor(meta_sub$Subgroup)
meta_sub$Subgroup <- relevel(meta_sub$Subgroup, ref = "Group3")

# Metastasis: clean and set M0 as reference
raw_met <- as.character(meta_sub$Met.status..1.Met..0.M0.)
raw_met[raw_met %in% c("0", "0.0", "0 ", "0.00")] <- "M0"
raw_met[raw_met %in% c("1", "1.0", "1 ", "1.00")] <- "Met"
raw_met[grepl("M0",  raw_met, ignore.case = TRUE)] <- "M0"
raw_met[grepl("Met", raw_met, ignore.case = TRUE)] <- "Met"
meta_sub$Met_status <- factor(raw_met)
meta_sub$Met_status <- relevel(meta_sub$Met_status, ref = "M0")

# Age group: set 4-10 years as reference
meta_sub$AgeGroup <- factor(as.character(meta_sub$AgeGroup))
meta_sub$AgeGroup <- relevel(meta_sub$AgeGroup, ref = "4-10")

#########################################
# (B) Full cohort Cox models            #
#########################################

cox_hist     <- coxph(Surv(OS..years., Dead) ~ histology,  data = meta_sub)
cox_subgroup <- coxph(Surv(OS..years., Dead) ~ Subgroup,   data = meta_sub)
cox_met      <- coxph(Surv(OS..years., Dead) ~ Met_status, data = meta_sub)
cox_age      <- coxph(Surv(OS..years., Dead) ~ AgeGroup,   data = meta_sub)

summary(cox_hist)
summary(cox_subgroup)
summary(cox_met)
summary(cox_age)

#########################################
# (C) Subgroup subsetting               #
#########################################

# Split cohort into G3/4 and SHH for subgroup-specific analyses
meta_g34 <- meta_sub[meta_sub$Subgroup %in% c("Group3", "Group4"), ]
meta_shh <- meta_sub[meta_sub$Subgroup == "SHH", ]

cat("G3/4 samples:", nrow(meta_g34), "\n")
cat("SHH  samples:", nrow(meta_shh), "\n")

# Relevel age and metastasis within subgroups
meta_g34$AgeGroup  <- relevel(factor(as.character(meta_g34$AgeGroup)),  ref = "4-10")
meta_g34$Met_status <- relevel(factor(as.character(meta_g34$Met_status)), ref = "M0")
meta_g34$histology  <- relevel(factor(as.character(meta_g34$histology)),  ref = "MBEN")

meta_shh$AgeGroup  <- relevel(factor(as.character(meta_shh$AgeGroup)),  ref = "4-10")
meta_shh$Met_status <- relevel(factor(as.character(meta_shh$Met_status)), ref = "M0")
meta_shh$histology  <- relevel(factor(as.character(meta_shh$histology)),  ref = "MBEN")

#########################################
# (D) G3/4 subgroup Cox models          #
#########################################

cox_age_g34  <- coxph(Surv(OS..years., Dead) ~ AgeGroup,   data = meta_g34)
cox_hist_g34 <- coxph(Surv(OS..years., Dead) ~ histology,  data = meta_g34)
cox_met_g34  <- coxph(Surv(OS..years., Dead) ~ Met_status, data = meta_g34)

summary(cox_age_g34)
summary(cox_hist_g34)
summary(cox_met_g34)

#########################################
# (E) SHH subgroup Cox models           #
#########################################

cox_age_shh  <- coxph(Surv(OS..years., Dead) ~ AgeGroup,   data = meta_shh)
cox_hist_shh <- coxph(Surv(OS..years., Dead) ~ histology,  data = meta_shh)
cox_met_shh  <- coxph(Surv(OS..years., Dead) ~ Met_status, data = meta_shh)

summary(cox_age_shh)
summary(cox_hist_shh)
summary(cox_met_shh)


########################################################################################################################
########################################################################################################################
#                                   5. MYC/MYCN AMPLIFICATION SURVIVAL ANALYSIS
########################################################################################################################
#
# OVERVIEW:
# Tests whether MYC or MYCN amplification is associated with overall survival
# in the full cohort and within Group 3/4 and SHH subgroups.
#
# OUTPUT:
#   - km_myc_all         : full cohort 3-group KM (None / MYC / MYCN)
#   - km_myc_g34         : G3/4 MYC KM
#   - km_mycn_g34        : G3/4 MYCN KM
#   - km_mycn_shh        : SHH MYCN KM
#   - myc_amp_ids        : MYC-amplified sample IDs (used downstream)
#   - mycn_amp_ids       : MYCN-amplified sample IDs (used downstream)
########################################################################################################################

# -----------------------------------------------------------------------------
# 5.1 Load amplification IDs
# -----------------------------------------------------------------------------

myc_ids  <- read.csv("mycTayDems.csv",  stringsAsFactors = FALSE)
mycn_ids <- read.csv("mycnTayDems.csv", stringsAsFactors = FALSE)

vec_myc_ids  <- myc_ids$Study_ID
vec_mycn_ids <- mycn_ids$Study_ID

myc_amp_ids  <- vec_myc_ids
mycn_amp_ids <- vec_mycn_ids

# -----------------------------------------------------------------------------
# 5.2 Assign amplification status to subgroup metadata
# -----------------------------------------------------------------------------

# Full cohort: 3-level factor
meta_sub$MYC_group <- "None"
meta_sub$MYC_group[meta_sub$Study_ID %in% vec_myc_ids]  <- "MYC-amplified"
meta_sub$MYC_group[meta_sub$Study_ID %in% vec_mycn_ids] <- "MYCN-amplified"
meta_sub$MYC_group <- factor(
  meta_sub$MYC_group,
  levels = c("None", "MYC-amplified", "MYCN-amplified")
)

# Group 3/4: binary MYC and MYCN
meta_g34$MYC <- factor(
  ifelse(meta_g34$Study_ID %in% vec_myc_ids, "Amplified", "Not amplified"),
  levels = c("Not amplified", "Amplified")
)
meta_g34$MYCN <- factor(
  ifelse(meta_g34$Study_ID %in% vec_mycn_ids, "Amplified", "Not amplified"),
  levels = c("Not amplified", "Amplified")
)

# SHH: binary MYCN only (MYC amplification absent/rare in SHH)
meta_shh$MYCN <- factor(
  ifelse(meta_shh$Study_ID %in% vec_mycn_ids, "Amplified", "Not amplified"),
  levels = c("Not amplified", "Amplified")
)

# Check counts
message("Full cohort MYC_group:"); print(table(meta_sub$MYC_group))
message("G3/4 MYC:");              print(table(meta_g34$MYC))
message("G3/4 MYCN:");             print(table(meta_g34$MYCN))
message("SHH MYCN:");              print(table(meta_shh$MYCN))

# -----------------------------------------------------------------------------
# 5.3 Full cohort: MYC/MYCN KM
# -----------------------------------------------------------------------------

fit_myc_all <- survfit(Surv(OS..years., Dead) ~ MYC_group, data = meta_sub)

km_myc_all <- ggsurvplot(
  fit_myc_all, data = meta_sub,
  title             = "Overall survival by MYC/MYCN amplification status",
  legend.title      = "Amplification",
  legend.labs       = c("None", "MYC-amplified", "MYCN-amplified"),
  risk.table        = TRUE,
  risk.table.height = 0.28,
  xlim              = c(0, 10),
  xlab              = "Overall survival (years)",
  ylab              = "Survival probability",
  pval              = TRUE,
  pval.coord        = c(0.3, 0.15),
  palette           = c("grey60", "#d73027", "#fd8d3c"),
  ggtheme           = theme_classic(base_size = 11),
  tables.theme      = theme_cleantable()
)

print(km_myc_all)

# -----------------------------------------------------------------------------
# 5.4 Group 3/4: MYC and MYCN side-by-side
# -----------------------------------------------------------------------------

fit_myc_g34 <- survfit(Surv(OS..years., Dead) ~ MYC, data = meta_g34)
km_myc_g34  <- ggsurvplot(
  fit_myc_g34, data = meta_g34,
  title             = "Group 3/4: MYC amplification",
  legend.title      = "MYC",
  legend.labs       = c("Not amplified", "Amplified"),
  risk.table        = TRUE,
  risk.table.height = 0.25,
  xlim              = c(0, 10),
  xlab              = "Overall survival (years)",
  ylab              = "Survival probability",
  pval              = TRUE,
  pval.coord        = c(0.3, 0.15),
  palette           = c("grey60", "#d73027"),
  ggtheme           = theme_classic(base_size = 11),
  tables.theme      = theme_cleantable()
)

fit_mycn_g34 <- survfit(Surv(OS..years., Dead) ~ MYCN, data = meta_g34)
km_mycn_g34  <- ggsurvplot(
  fit_mycn_g34, data = meta_g34,
  title             = "Group 3/4: MYCN amplification",
  legend.title      = "MYCN",
  legend.labs       = c("Not amplified", "Amplified"),
  risk.table        = TRUE,
  risk.table.height = 0.25,
  xlim              = c(0, 10),
  xlab              = "Overall survival (years)",
  ylab              = "Survival probability",
  pval              = TRUE,
  pval.coord        = c(0.3, 0.15),
  palette           = c("grey60", "#fd8d3c"),
  ggtheme           = theme_classic(base_size = 11),
  tables.theme      = theme_cleantable()
)

arrange_ggsurvplots(
  list(km_myc_g34, km_mycn_g34),
  ncol = 2, nrow = 1
)

# -----------------------------------------------------------------------------
# 5.5 SHH: MYCN only
# -----------------------------------------------------------------------------

fit_mycn_shh <- survfit(Surv(OS..years., Dead) ~ MYCN, data = meta_shh)
km_mycn_shh  <- ggsurvplot(
  fit_mycn_shh, data = meta_shh,
  title             = "SHH: MYCN amplification",
  legend.title      = "MYCN",
  legend.labs       = c("Not amplified", "Amplified"),
  risk.table        = TRUE,
  risk.table.height = 0.25,
  xlim              = c(0, 10),
  xlab              = "Overall survival (years)",
  ylab              = "Survival probability",
  pval              = TRUE,
  pval.coord        = c(0.3, 0.15),
  palette           = c("grey60", "#fd8d3c"),
  ggtheme           = theme_classic(base_size = 11),
  tables.theme      = theme_cleantable()
)

print(km_mycn_shh)


########################################################################################################################
#                                   5.6 APPENDIX: Clinical variable KM panels
#
# PURPOSE:
#   Visualise the clinical variables from the univariable Cox analysis (Tables 3 & 4)
#   as Kaplan-Meier survival curves. These supplement the main results and are
#   referenced in section 2.4 of the dissertation.
#
# VARIABLES SHOWN:
#   - Metastatic stage (G3/4 and SHH)
#   - Age group (G3/4 and SHH) — three groups: 0-3, 4-10, 10+ years
#   - Intra-subgroup: Group 3 vs Group 4
#   - SHH histology: Classic vs Desmoplastic vs LCA vs MBEN
#
# OUTPUT:
#   - fig_supp: 3x2 panel cowplot figure (Supplementary Figure X)
#
# NOTE:
#   age_group is created here as a three-level factor (4-10, 0-3, 10+)
#   with 4-10 as reference, consistent with cox_age_g34 and cox_age_shh.
#   The 10+ group was absent from the original two-level coding and is
#   added here to fully represent the age distribution.
########################################################################################################################

# -----------------------------------------------------------------------------
# Age group factor — three levels with 4-10 as reference
# -----------------------------------------------------------------------------

meta_g34$age_group <- factor(meta_g34$AgeGroup, levels = c("4-10","0-3","10+"))
meta_shh$age_group <- factor(meta_shh$AgeGroup, levels = c("4-10","0-3","10+"))

# -----------------------------------------------------------------------------
# Metastatic stage KM — G3/4 and SHH
# -----------------------------------------------------------------------------

fit_meta_g34 <- survfit(Surv(OS..years., Dead) ~ Metastasis_status, data = meta_g34)
km_meta_g34  <- ggsurvplot(fit_meta_g34, data = meta_g34,
                           title = "Group 3/4: metastasis", legend.title = "Metastasis",
                           legend.labs = c("M0","Metastatic"), risk.table = TRUE, risk.table.height = 0.25,
                           xlim = c(0,10), xlab = "Overall survival (years)", ylab = "Survival probability",
                           pval = TRUE, pval.coord = c(0.3,0.15), palette = c("#F08080","#20B2AA"),
                           ggtheme = theme_classic(base_size = 11), tables.theme = theme_cleantable())

fit_meta_shh <- survfit(Surv(OS..years., Dead) ~ Metastasis_status, data = meta_shh)
km_meta_shh  <- ggsurvplot(fit_meta_shh, data = meta_shh,
                           title = "SHH: metastasis", legend.title = "Metastasis",
                           legend.labs = c("M0","Metastatic"), risk.table = TRUE, risk.table.height = 0.25,
                           xlim = c(0,10), xlab = "Overall survival (years)", ylab = "Survival probability",
                           pval = TRUE, pval.coord = c(0.3,0.15), palette = c("#F08080","#20B2AA"),
                           ggtheme = theme_classic(base_size = 11), tables.theme = theme_cleantable())

# -----------------------------------------------------------------------------
# Age group KM — G3/4 and SHH
# -----------------------------------------------------------------------------

fit_age_g34 <- survfit(Surv(OS..years., Dead) ~ age_group, data = meta_g34)
km_age_g34  <- ggsurvplot(fit_age_g34, data = meta_g34,
                          title = "Group 3/4: age group", legend.title = "Age group (years)",
                          legend.labs = c("4-10","0-3","10+"), risk.table = TRUE, risk.table.height = 0.25,
                          xlim = c(0,10), xlab = "Overall survival (years)", ylab = "Survival probability",
                          pval = TRUE, pval.coord = c(0.3,0.15), palette = c("#F08080","#20B2AA","#9370DB"),
                          ggtheme = theme_classic(base_size = 11), tables.theme = theme_cleantable())

fit_age_shh <- survfit(Surv(OS..years., Dead) ~ age_group, data = meta_shh)
km_age_shh  <- ggsurvplot(fit_age_shh, data = meta_shh,
                          title = "SHH: age group", legend.title = "Age group (years)",
                          legend.labs = c("4-10","0-3","10+"), risk.table = TRUE, risk.table.height = 0.25,
                          xlim = c(0,10), xlab = "Overall survival (years)", ylab = "Survival probability",
                          pval = TRUE, pval.coord = c(0.3,0.15), palette = c("#F08080","#20B2AA","#9370DB"),
                          ggtheme = theme_classic(base_size = 11), tables.theme = theme_cleantable())

# -----------------------------------------------------------------------------
# Intra-subgroup KM — Group 3 vs Group 4
# Significant in univariable Cox (HR=0.55, p=0.001)
# -----------------------------------------------------------------------------

fit_subgroup_g34 <- survfit(Surv(OS..years., Dead) ~ Subgroup, data = meta_g34)
km_subgroup_g34  <- ggsurvplot(fit_subgroup_g34, data = meta_g34,
                               title = "Group 3/4: intra-subgroup", legend.title = "Subgroup",
                               legend.labs = c("Group 3","Group 4"), risk.table = TRUE, risk.table.height = 0.25,
                               xlim = c(0,10), xlab = "Overall survival (years)", ylab = "Survival probability",
                               pval = TRUE, pval.coord = c(0.3,0.15), palette = c("#FFFF00","#008000"),
                               ggtheme = theme_classic(base_size = 11), tables.theme = theme_cleantable())

# -----------------------------------------------------------------------------
# SHH histology KM — Classic as reference (MBEN excluded from Cox due to 0 events)
# Desmoplastic significantly protective vs Classic (HR=0.22, p=0.006)
# -----------------------------------------------------------------------------

meta_shh$histology <- relevel(factor(as.character(meta_shh$histology)), ref = "Classic")
fit_hist_shh <- survfit(Surv(OS..years., Dead) ~ histology, data = meta_shh)
km_hist_shh  <- ggsurvplot(fit_hist_shh, data = meta_shh,
                           title = "SHH: histology", legend.title = "Histology",
                           legend.labs = c("Classic","Desmoplastic","LCA","MBEN"), risk.table = TRUE, risk.table.height = 0.25,
                           xlim = c(0,10), xlab = "Overall survival (years)", ylab = "Survival probability",
                           pval = TRUE, pval.coord = c(0.3,0.15), palette = c("#878787","#fc8d59","black","red4"),
                           ggtheme = theme_classic(base_size = 11), tables.theme = theme_cleantable())

# -----------------------------------------------------------------------------
# Combine into 3x2 supplementary figure using cowplot
# Row 1: Metastasis (G3/4 | SHH)
# Row 2: Age group  (G3/4 | SHH)
# Row 3: Intra-subgroup (G3/4) | SHH histology
# -----------------------------------------------------------------------------

library(cowplot)

fig_supp <- plot_grid(
  km_meta_g34$plot,     km_meta_shh$plot,
  km_age_g34$plot,      km_age_shh$plot,
  km_subgroup_g34$plot, km_hist_shh$plot,
  ncol   = 2, nrow = 3,
  align  = "hv",
  labels = c("A","B","C","D","E","F")
)

print(fig_supp)


########################################################################################################################
########################################################################################################################
#
#          6. PENALIZED COX REGRESSION (Elastic Net Risk Score) — GROUP 3/4
#
########################################################################################################################
########################################################################################################################
#
# OVERVIEW:
# Build a multi-gene prognostic risk score using penalized Cox regression.
# Unlike the univariable gene-by-gene approach in Section 4, penalized regression
# fits all genes simultaneously and uses regularization to:
#   - Handle high-dimensional data (more genes than samples)
#   - Perform automatic feature selection (lasso/elastic net)
#   - Prevent overfitting
#
# Approach:
#   - Test three penalty types: ridge (α=0), elastic net (α=0.5), lasso (α=1)
#   - Use 10-fold cross-validation to select optimal lambda
#   - Extract selected genes and compute continuous risk scores
#   - Validate risk score as a survival predictor
#
# OUTPUT:
#   - cv_fits: cross-validation results for each alpha
#   - selected_genes: genes with non-zero coefficients
#   - Risk score for each patient
#
########################################################################################################################

# 1) building a clinical dataset to model survival outcome which,
# keeps only patients with survival time + event status so its in the format 0/1
dat <- meta
dat <- dat[!is.na(dat$OS..years.) & !is.na(dat$Dead), ]
dat$Dead <- as.integer(dat$Dead)
stopifnot(all(dat$Dead %in% c(0, 1)))

#2) Building the expression matrix - this pulls out only the expression columns ,
#converts them to numeric values and makes a matrix that glmnet wants 
#- rows = samples, cols = genes

expr_df <- as.data.frame(expr)
sample_cols <- grep("^MB_SubtypeStudy_", names(expr_df), value = TRUE) 
gene_col <- "EnsemblGeneID_from_ensemblv77"
gene_ids <- expr_df[[gene_col]]

X_gxs <- as.matrix(sapply(expr_df[, sample_cols, drop = FALSE], as.numeric))
rownames(X_gxs) <- gene_ids
x_all <- t(X_gxs)
rownames(x_all) <- sample_cols

# Subset to G3/4 only
dat_g34 <- subset(dat, Subgroup %in% c("Group3", "Group4"))


#3) Align the expression to clinical data with the patients, same order
# As an outcome - this step ensures x_glm corresponds to the correct patient
# in dat_glm + removes genes that can be used by glmnet
id_col <- "Study_ID"  
common_ids <- intersect(dat_g34[[id_col]], rownames(x_all))
stopifnot(length(common_ids) > 10)

dat_glm <- dat_g34[match(common_ids, dat_g34[[id_col]]), ]
x_glm   <- x_all[common_ids, , drop = FALSE]

#Remove genes with any NA (glmnet cannot handle NA)
keep_genes <- colSums(is.na(x_glm)) == 0
x_glm <- x_glm[, keep_genes, drop = FALSE]

#4) create the survival outcome object (y) - which label the model is learning to predict 
y <- with(dat_glm, Surv(OS..years., Dead))

#5) setting up a 10-fold cross validation (folds)
set.seed(1) 
K <- 10
foldid <- sample(rep(1:K, length.out = nrow(x_glm)))

################################################################################
# Defining ALPHA values - this is where alpha is set 0 , 0.5 , 1 
# Controls the kind of penalty the apply to the coefficients
# 0 = ridge
# 0.5 = elastic net
# 1 = lasso
################################################################################

alphas <- c(0, 0.5, 1)

#1) Fit cv.glmnet for each alpha - within cv.glmnet it tries lots of Lambda values automatically
# For each lambda it fits a penalised cox model and evaluates prediction error via 10-fold CV

cat("Starting cv.glmnet alpha loop at:", as.character(Sys.time()), "\n")

cv_fits <- pblapply(alphas, function(a) {
  cat("  Running alpha =", a, " at:", as.character(Sys.time()), "\n")
  cv.glmnet(
    x_glm, y,
    family = "cox",
    alpha  = a,
    nfolds = K,
    foldid = foldid,
    type.measure = "deviance",
    trace.it = 1
  )
})

cat("Finished cv.glmnet alpha loop at:", as.character(Sys.time()), "\n")

names(cv_fits) <- paste0("alpha_", alphas)


################################################################################
# Lambda results and comparison of Alphas G3/4
################################################################################

cv_summary <- data.frame(
  alpha      = alphas,
  lambda_min = sapply(cv_fits, function(f) f$lambda.min),
  lambda_1se = sapply(cv_fits, function(f) f$lambda.1se),
  cvm_min    = sapply(cv_fits, function(f) min(f$cvm))
)
print(cv_summary)

best_alpha <- cv_summary$alpha[which.min(cv_summary$cvm_min)]
cat("G3/4 chosen alpha =", best_alpha, "\n")

# Based on these lambda results and comparison of alpha i will be using alpha = 0.5

################################################################################
# CV comparison curves of alpha 0 , 0.5 , 1 G3/4
################################################################################

# --- CV comparison curves across alpha (same y-scale) ---

fits_to_plot <- list(
  ridge = cv_fits[["alpha_0"]],
  enet  = cv_fits[["alpha_0.5"]],
  lasso = cv_fits[["alpha_1"]]
)

# y-limits across all three fits
yl <- range(unlist(lapply(fits_to_plot, function(f) f$cvm)), na.rm = TRUE)

par(mfrow = c(1,3), mar = c(4,4,2,1), oma = c(0,0,3,0))

plot(fits_to_plot$ridge, main = "", ylim = yl)
mtext("Alpha = 0 (ridge)", side = 3, line = -1, cex = 0.9)

plot(fits_to_plot$enet,  main = "", ylim = yl)
mtext("Alpha = 0.5 (elastic net)", side = 3, line = -1, cex = 0.9)

plot(fits_to_plot$lasso, main = "", ylim = yl)
mtext("Alpha = 1 (lasso)", side = 3, line = -1, cex = 0.9)

mtext("G3/4: 10-fold CV curves across alpha (ridge / elastic net / lasso)",
      outer = TRUE, cex = 1.1, line = 1)

par(mfrow = c(1,1))


################################################################################
# Alpha was chosen based on elastic net (alpha = 0.5): extract selected genes at lambda.min
# Re-fit an unpenalised Cox model on the selected genes
# 1) Obtaining conventional Cox coefficients (no penalty)
# 2) Compute a continuous linear predictor (risk score)
# 3) Validate that manual x-beta matches predict
################################################################################

# Alpha = 0.5 (elastic net) was chosen because CV deviance was comparable across alphas,
# while providing an interpretable model.
# - For this dataset, lambda.1se yields 0 non-zero coefficients for alpha=0.5,
# so lambda.min is used to avoid a null model.

best_fit_g34 <- cv_fits[["alpha_0.5"]]

best_fit_g34$lambda.min
best_fit_g34$lambda.1se

# 1. Extract coefficients at optimal lambda
coef_matrix <- coef(best_fit_g34, s = "lambda.min")

# 2. Identify selected genes (non-zero coefficients)
selected_coefs <- coef_matrix[coef_matrix[,1] != 0, , drop = FALSE]
selected_genes <- rownames(selected_coefs)

# Print out of the selected genes
print(paste("Selected", length(selected_genes), "genes"))
print(selected_genes)

#3 Build a modelling dataframe for a standard Cox model (unpenalised)
# This uses the same samples and the selected gene columns only.

y2 <- as.character(y)
# event status
OS_Status <- rep(1, length(y2))
OS_Status[grep("\\+", y2)] <- 0
# time: remove "+" censor marker and convert to numeric
OS_Time <- as.numeric(gsub("\\+","", y2))
# Building the modelling data frame with: time/event (survival outcome) - expression matrix on selected genes

cox_data <- data.frame(
  time  = OS_Time,
  event = OS_Status,
  x_glm[, selected_genes]  # subset to selected genes only
)

##############################################################
# 4) Build Cox formula using selected genes
##############################################################

formula_string <- paste("Surv(time, event) ~", 
                        paste(selected_genes, collapse = " + "))
cox_formula <- as.formula(formula_string)

##############################################################
# 5. Fit the unpenalised Cox model (refit step)
##############################################################

cox_model <- coxph(cox_formula, data = cox_data)

# 6. View results containing HR, SE, p-values on selected genes
summary(cox_model)

# 7. Extract coefficients for validation
final_coefs <- coef(cox_model)
print(final_coefs)

############################################################################################
# 8) Compute continuous risk scores (linear predictor) and validate equivalence
############################################################################################
# Risk score = weighted sum of gene expression
risk_score_discovery <- as.matrix(x_glm[, names(final_coefs)]) %*% final_coefs

# This is equivalent to:
risk_score_discovery2 <- predict(cox_model, type = "lp")  # linear predictor

#write.csv(risk_score_discovery, file="G3G4_predictor.csv")

################################################################################
# 9) Compare refit Cox risk score to glmnet risk score (optional diagnostic)
################################################################################

# Validation check: risk_score_discovery and risk_score_discovery2 should match closely
# Density plot of continuous risk score (with rug + N)
plot(density(risk_score_discovery2),
     main = "G3/4 risk score distribution",
     xlab = "Risk score (linear predictor)",
     ylab = "Density")
rug(risk_score_discovery2)
title(sub = paste0("N = ", length(risk_score_discovery2)))

##################################################
## Density plot + rug plot with labelled subgroup
##################################################
plot(density(risk_score_discovery2))

summary(coxph(Surv(OS_Time, OS_Status) ~ risk_score_discovery2))


######################################################################
# 
# Further analysis of G3/4
# 
# Forest Plot of Top 10 Prognostic Genes from Elastic Net Model
# Dataset: Group 3/4 tumours | Model: glmnet with alpha = 0.5
######################################################################

######################################################################
## Step 1: Get the best model and extract non-zero coefficients
######################################################################

best_fit <- cv_fits[["alpha_0.5"]]
stopifnot(!is.null(best_fit))

lambda_idx <- which.min(abs(best_fit$lambda - best_fit$lambda.min))

beta_matrix <- best_fit$glmnet.fit$beta
betas <- as.matrix(beta_matrix[, lambda_idx, drop = FALSE])
betas <- betas[betas[,1] != 0, , drop = FALSE]

if (nrow(betas) == 0) stop("No genes selected at this lambda - try a smaller value?")

n_top <- min(10, nrow(betas))
top_genes <- head(rownames(betas)[order(abs(betas[,1]), decreasing = TRUE)], n_top)

cat("Found", nrow(betas), "genes with non-zero coefficients\n")
cat("Selecting top", n_top, "by effect size\n")

######################################################################
## Step 2: Convert Ensembl IDs to gene symbols
######################################################################

get_gene_symbols <- function(ids) {
  if (grepl("^ENSG", ids[1])) {
    clean_ids <- gsub("\\..*", "", ids)
    symbols <- mapIds(org.Hs.eg.db,
                      keys = clean_ids,
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")
  } else if (grepl("^[0-9]+$", ids[1])) {
    symbols <- mapIds(org.Hs.eg.db,
                      keys = ids,
                      column = "SYMBOL",
                      keytype = "ENTREZID",
                      multiVals = "first")
  } else {
    symbols <- setNames(ids, ids)
  }
  
  symbols <- as.character(symbols)
  symbols[is.na(symbols)] <- ids[is.na(symbols)]
  return(symbols)
}

gene_names <- get_gene_symbols(top_genes)

######################################################################
## Step 3: Univariable Cox regression for each gene
######################################################################

fit_cox_model <- function(gene_id, gene_name) {
  df <- data.frame(
    time       = dat_glm$OS..years.,
    status     = dat_glm$Dead,
    expression = x_glm[, gene_id]
  )
  
  cox_fit <- coxph(Surv(time, status) ~ expression, data = df)
  cox_summary <- summary(cox_fit)
  
  data.frame(
    symbol      = gene_name,
    glmnet_beta = as.numeric(betas[gene_id, 1]),
    HR          = cox_summary$coef[1, "exp(coef)"],
    lower_95    = cox_summary$conf.int[1, "lower .95"],
    upper_95    = cox_summary$conf.int[1, "upper .95"],
    pval        = cox_summary$coef[1, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
}

results <- do.call(rbind, Map(fit_cox_model, top_genes, gene_names))
results$padj <- p.adjust(results$pval, method = "BH")
results <- results[order(results$HR), ]

cat("\n--- Cox Regression Results ---\n")
print(results, row.names = FALSE)

######################################################################
## Step 4: Forest plot — gene names in italics
######################################################################

png("G34_forest_plot.png", width = 10, height = 8, units = "in", res = 150)

par(mar = c(5, 10, 4, 6))
y_positions <- seq_len(nrow(results))

plot(results$HR, y_positions,
     xlim = range(c(results$lower_95, results$upper_95, 1), na.rm = TRUE),
     yaxt = "n", 
     xlab = "Hazard Ratio (95% CI)", 
     ylab = "",
     log  = "x", 
     pch  = 16, 
     cex  = 1.2,
     main = paste0("G3/4: Top ", n_top, " Prognostic Genes (Elastic Net \u03b1=0.5)"))

arrows(results$lower_95, y_positions, 
       results$upper_95, y_positions,
       angle = 90, code = 3, length = 0.04)

axis(2, at = y_positions, labels = results$symbol,
     las = 1, cex.axis = 0.9, font = 3)   # font = 3 = italic

abline(v = 1, lty = 2, col = "grey40")

text(results$upper_95, y_positions,
     labels = sprintf("HR=%.2f (p=%.1e)", results$HR, results$pval),
     pos = 4, cex = 0.7, xpd = NA)

dev.off()
cat("\nForest plot saved to: G34_forest_plot.png\n")

######################################################################
## Step 5: Save results
######################################################################

write.csv(results, "G34_top10_forest_df.csv", row.names = FALSE)
cat("Results saved to: G34_top10_forest_df.csv\n")


##############################################################################
# Stratify G3/4 patients into risk tertiles (Low / Mid / High)
# Tertiles chosen to give equal-sized groups without assuming a threshold.
##############################################################################

rs    <- as.numeric(risk_score_discovery2)
cuts  <- quantile(rs, probs = c(1/3, 2/3), na.rm = TRUE, type = 7)

dat_glm$risk_score <- rs
dat_glm$risk_tertile <- cut(
  rs,
  breaks         = c(-Inf, cuts[1], cuts[2], Inf),
  labels         = c("Low", "Mid", "High"),
  include.lowest = TRUE,
  right          = TRUE
)

# check: score vector must align with the dataframe row count
stopifnot("Risk score length doesn't match dat_glm rows" = length(rs) == nrow(dat_glm))

# Quick QC — check group sizes are roughly balanced, ranges don't overlap,
# and that events (Dead) are distributed across all three tertiles
message("G3/4 tertile sizes:")
print(table(dat_glm$risk_tertile))

message("Cutpoints (33rd / 67th percentile):")
print(cuts)

message("Score range per tertile:")
print(tapply(dat_glm$risk_score, dat_glm$risk_tertile, range))

message("Event counts by tertile:")
print(with(dat_glm, table(risk_tertile, Dead)))


###########################################################################
# G3/4 KM curves, one per risk tertile
###########################################################################

km_fit <- survfit(
  Surv(OS..years., Dead) ~ risk_tertile,
  data = dat_glm
)

# ── Log-rank p-value across all three groups 
lr_test <- survdiff(Surv(OS..years., Dead) ~ risk_tertile, data = dat_glm)
lr_p    <- 1 - pchisq(lr_test$chisq, df = length(lr_test$n) - 1)

# Plot
ggsurvplot(
  km_fit,
  data              = dat_glm,
  palette           = c("#2166ac", "#f4a582", "#d6604d"),  # Low=blue, Mid=orange, High=red
  conf.int          = TRUE,
  conf.int.alpha    = 0.12,
  risk.table        = TRUE,
  risk.table.height = 0.25,
  xlab              = "Time (years)",
  ylab              = "Overall survival probability",
  title             = "Kaplan-Meier survival by risk tertile \u2014 G3/4 cohort",
  legend.title      = "Risk group",
  legend.labs       = c("Low", "Mid", "High"),
  pval              = FALSE,
  ggtheme           = theme_classic(base_size = 13),
  xlim              = c(0, 10)
) |>
  (\(p) {
    p$plot <- p$plot +
      annotate(
        "text",
        x     = max(dat_glm$OS..years., na.rm = TRUE) * 0.65,
        y     = 0.92,
        label = paste0("Log-rank p = ", format.pval(lr_p, digits = 2, eps = 0.001)),
        size  = 4,
        hjust = 0
      )
    p
  })()


###########################################################################
# G3/4 cohort — Clinical annotation figure
# Page 1: Oncoprint ordered by ascending Cox risk score
# Page 2: Ridge plots of risk score by clinical variable
#
# MYC and MYCN were not originally in dat_glm — they are merged in
# from myc_amp_ids and mycn_amp_ids (character vectors of Study_IDs
# with confirmed amplification) at the start of this script.
#
# Adapted from my G3/4 oncoprint
# Input: dat_glm with risk_score and risk_tertile already computed.
###########################################################################

# ── 1. Merge MYC and MYCN amplification into dat_glm
dat_glm$MYC  <- ifelse(dat_glm$Study_ID %in% myc_amp_ids,  1, 0)
dat_glm$MYCN <- ifelse(dat_glm$Study_ID %in% mycn_amp_ids, 1, 0)

cat("MYC amplification in G3/4:\n");  print(table(dat_glm$MYC,  useNA = "ifany"))
cat("MYCN amplification in G3/4:\n"); print(table(dat_glm$MYCN, useNA = "ifany"))

# ── 2. Prepare data for oncoprint
dat_plot      <- dat_glm[order(dat_glm$risk_score), ]
dat_plot$Dead <- as.numeric(dat_plot$Dead)

# ── 3. Clean NAs ─────────────────────────────────────────────────────────────
dat_plot$Gender    <- ifelse(is.na(dat_plot$Gender),    "Unknown", as.character(dat_plot$Gender))
dat_plot$histology <- ifelse(is.na(dat_plot$histology) | dat_plot$histology == "", "Unknown", dat_plot$histology)
dat_plot$Subgroup  <- ifelse(is.na(dat_plot$Subgroup),  "Unknown", as.character(dat_plot$Subgroup))
dat_plot$Met.status..1.Met..0.M0. <- ifelse(
  is.na(dat_plot$Met.status..1.Met..0.M0.), "Unknown",
  as.character(dat_plot$Met.status..1.Met..0.M0.)
)
dat_plot$MYC  <- as.character(dat_plot$MYC)
dat_plot$MYCN <- as.character(dat_plot$MYCN)

# ── 4. Colour palettes ────────────────────────────────────────────────────────
tertile_col_fun  <- c("Low" = "#2166ac", "Mid" = "#f4a582", "High" = "#d6604d")
gender_col_fun   <- c("M" = "white",     "F" = "black",     "Unknown" = "lightgrey")
subgroup_col_fun <- c("Group3" = "#FFFF00", "Group4" = "#008000", "Unknown" = "lightgrey")
hist_col_fun     <- c("Classic" = "white", "Desmoplastic" = "#fc8d59",
                      "LCA" = "black",     "MBEN" = "red4", "Unknown" = "lightgrey")
mstage_col_fun   <- c("0" = "white", "1" = "black", "Unknown" = "lightgrey")
myc_col_fun      <- c("0" = "white", "1" = "#e31a1c")
mycn_col_fun     <- c("0" = "white", "1" = "#ff7f00")

# ── 5. Age dot annotation (capped at 20) ──────────────────────────────────────
Age_capped <- pmin(dat_plot$Age, 20)

age_dot_anno <- HeatmapAnnotation(
  Age = anno_points(
    Age_capped,
    ylim = c(0, 20),
    pch = ifelse(is.na(dat_plot$Age), 4, 20),
    gp = gpar(col = ifelse(is.na(dat_plot$Age), "red", "grey")),
    axis_param = list(at = c(0, 5, 10, 15, 20)),
    baseline = 3,
    baseline_gp = gpar(col = "grey", lty = 2, lwd = 1)
  ),
  annotation_name_side = "left",
  annotation_height = unit(2, "cm")
)

# ── 6. OS dot annotation (capped at 10 years) ────────────────────────────────
OS <- pmin(dat_plot$OS..years., 10)
pchSetup <- ifelse(dat_plot$Dead == 1 & dat_plot$OS..years. <  10,  4,
                   ifelse(dat_plot$Dead == 0 & dat_plot$OS..years. <  10, 20,
                          ifelse(dat_plot$Dead == 0 & dat_plot$OS..years. >= 10,  2, 2)))

OS_dot_anno <- HeatmapAnnotation(
  OS = anno_points(
    OS,
    pch = pchSetup,
    gp = gpar(col = ifelse(dat_plot$Dead == 1, "firebrick3", "lightgrey")),
    ylim = c(0, 10),
    axis_param = list(at = c(0, 5, 10)),
    baseline_gp = gpar(col = "lightgrey", lty = 2, lwd = 1)
  ),
  annotation_name_side = "left",
  annotation_height = unit(2, "cm")
)

# ── 7. Categorical bar annotations ───────────────────────────────────────────
cat_anno <- HeatmapAnnotation(
  `Risk tertile` = anno_simple(as.character(dat_plot$risk_tertile), col = tertile_col_fun,  border = TRUE),
  Gender         = anno_simple(dat_plot$Gender,                     col = gender_col_fun,   border = TRUE),
  Subgroup       = anno_simple(dat_plot$Subgroup,                   col = subgroup_col_fun, border = TRUE),
  Histology      = anno_simple(dat_plot$histology,                  col = hist_col_fun,     border = TRUE),
  `M-stage`      = anno_simple(dat_plot$Met.status..1.Met..0.M0.,  col = mstage_col_fun,   border = TRUE),
  MYC            = anno_simple(dat_plot$MYC,                        col = myc_col_fun,      border = TRUE),
  MYCN           = anno_simple(dat_plot$MYCN,                       col = mycn_col_fun,     border = TRUE),
  annotation_name_side = "left"
)

# ── 8. Base heatmap ───────────────────────────────────────────────────────────
rs_matrix <- matrix(dat_plot$risk_score, nrow = 1)

ht_g34 <- Heatmap(
  rs_matrix,
  name  = "Risk score",
  col   = colorRamp2(
    c(min(dat_plot$risk_score), 0, max(dat_plot$risk_score)),
    c("#2166ac", "white", "#d6604d")
  ),
  border              = TRUE,
  show_column_names   = FALSE,
  cluster_columns     = FALSE,
  cluster_rows        = FALSE,
  row_names_side      = "left",
  row_names_gp        = gpar(fontsize = 10),
  height              = unit(0.6, "cm"),
  show_heatmap_legend = TRUE,
  top_annotation      = c(cat_anno, age_dot_anno, OS_dot_anno)
)

# ── 9. Manual legends ────────────────────────────────────────────────────────
lgd_tertile  <- Legend(labels = c("Low", "Mid", "High"),
                       title = "Risk tertile",
                       legend_gp = gpar(fill = c("#2166ac", "#f4a582", "#d6604d")), border = TRUE)
lgd_subgroup <- Legend(labels = c("Group3", "Group4", "Unknown"),
                       title = "Subgroup",
                       legend_gp = gpar(fill = c("#FFFF00", "#008000", "lightgrey")), border = TRUE)
lgd_gender   <- Legend(labels = c("Male", "Female", "Unknown"),
                       title = "Gender",
                       legend_gp = gpar(fill = c("white", "black", "lightgrey")), border = TRUE)
lgd_hist     <- Legend(labels = c("Classic", "Desmoplastic", "LCA", "MBEN", "Unknown"),
                       title = "Histology",
                       legend_gp = gpar(fill = c("white", "#fc8d59", "black", "red4", "lightgrey")), border = TRUE)
lgd_mstage   <- Legend(labels = c("M0", "M+", "Unknown"),
                       title = "M-stage",
                       legend_gp = gpar(fill = c("white", "black", "lightgrey")), border = TRUE)
lgd_myc      <- Legend(labels = c("Not amplified", "Amplified"),
                       title = "MYC",
                       legend_gp = gpar(fill = c("white", "#e31a1c")), border = TRUE)
lgd_mycn     <- Legend(labels = c("Not amplified", "Amplified"),
                       title = "MYCN",
                       legend_gp = gpar(fill = c("white", "#ff7f00")), border = TRUE)

# ── 10. Ridge plots ───────────────────────────────────────────────────────────
fmt_p <- function(p) {
  if (p < 2.2e-16)    "p < 2.2e-16"
  else if (p < 0.001) paste0("p = ", formatC(p, format = "e", digits = 2))
  else                paste0("p = ", round(p, 3))
}

rs_range    <- range(dat_glm$risk_score, na.rm = TRUE)
ridge_theme <- theme_bw() + theme(legend.position = "none", panel.grid = element_blank())

# A — Subgroup
toPlot_A <- data.frame(risk_score = dat_glm$risk_score, group = dat_glm$Subgroup) %>%
  filter(!is.na(group), group != "Unknown")
p_A <- fmt_p(t.test(risk_score ~ group, data = toPlot_A)$p.value)

A <- ggplot(toPlot_A, aes(risk_score, group)) +
  geom_density_ridges(aes(fill = group), scale = 1.5, rel_min_height = 0.001, alpha = 0.45) +
  scale_fill_manual(values = c("Group3" = "#FFFF00", "Group4" = "#008000")) +
  xlim(rs_range) + ridge_theme +
  labs(x = "Risk score", y = "Subgroup") +
  annotate("text", x = rs_range[2] - 0.1, y = 3.3, label = p_A, hjust = 1, size = 3.5)

# B — Histology
toPlot_B <- data.frame(risk_score = dat_glm$risk_score, group = dat_glm$histology) %>%
  filter(!is.na(group), !group %in% c("Unknown", ""))
p_B <- fmt_p(anova(aov(risk_score ~ group, data = toPlot_B))$`Pr(>F)`[1])
hist_cols <- c("Classic" = "grey80", "Desmoplastic" = "#fc8d59", "LCA" = "black", "MBEN" = "red4")

B <- ggplot(toPlot_B, aes(risk_score, group)) +
  geom_density_ridges(aes(fill = group), scale = 1.5, rel_min_height = 0.001, alpha = 0.45) +
  scale_fill_manual(values = hist_cols) +
  xlim(rs_range) + ridge_theme +
  labs(x = "Risk score", y = "Histology") +
  annotate("text", x = rs_range[2] - 0.1, y = 4.3, label = p_B, hjust = 1, size = 3.5)

# C — M-stage
toPlot_C <- data.frame(risk_score = dat_glm$risk_score,
                       group = as.character(dat_glm$Met.status..1.Met..0.M0.)) %>%
  filter(!is.na(group), group != "Unknown")
p_C <- fmt_p(t.test(risk_score ~ group, data = toPlot_C)$p.value)

C <- ggplot(toPlot_C, aes(risk_score, group)) +
  geom_density_ridges(aes(fill = group), scale = 1.5, rel_min_height = 0.001, alpha = 0.45) +
  scale_fill_manual(values = c("0" = "grey80", "1" = "black")) +
  scale_y_discrete(labels = c("0" = "M0", "1" = "M+")) +
  xlim(rs_range) + ridge_theme +
  labs(x = "Risk score", y = "M-stage") +
  annotate("text", x = rs_range[2] - 0.1, y = 3.3, label = p_C, hjust = 1, size = 3.5)

# D — Gender
toPlot_D <- data.frame(risk_score = dat_glm$risk_score, group = dat_glm$Gender) %>%
  filter(!is.na(group), group != "Unknown")
p_D <- fmt_p(t.test(risk_score ~ group, data = toPlot_D)$p.value)

D <- ggplot(toPlot_D, aes(risk_score, group)) +
  geom_density_ridges(aes(fill = group), scale = 1.5, rel_min_height = 0.001, alpha = 0.45) +
  scale_fill_manual(values = c("F" = "#d73027", "M" = "#4575b4")) +
  xlim(rs_range) + ridge_theme +
  labs(x = "Risk score", y = "Gender") +
  annotate("text", x = rs_range[2] - 0.1, y = 3.3, label = p_D, hjust = 1, size = 3.5)

# E — MYC amplification
toPlot_E <- data.frame(risk_score = dat_glm$risk_score,
                       group = as.character(dat_glm$MYC)) %>%
  filter(!is.na(group))
p_E <- fmt_p(t.test(risk_score ~ group, data = toPlot_E)$p.value)

E <- ggplot(toPlot_E, aes(risk_score, group)) +
  geom_density_ridges(aes(fill = group), scale = 1.5, rel_min_height = 0.001, alpha = 0.45) +
  scale_fill_manual(values = c("0" = "grey80", "1" = "#e31a1c")) +
  scale_y_discrete(labels = c("0" = "Not amplified", "1" = "Amplified")) +
  xlim(rs_range) + ridge_theme +
  labs(x = "Risk score", y = "MYC") +
  annotate("text", x = rs_range[2] - 0.1, y = 3.3, label = p_E, hjust = 1, size = 3.5)

# F — MYCN amplification (named F_plot to avoid conflict with R's built-in F = FALSE)
toPlot_F <- data.frame(risk_score = dat_glm$risk_score,
                       group = as.character(dat_glm$MYCN)) %>%
  filter(!is.na(group))
p_F <- fmt_p(t.test(risk_score ~ group, data = toPlot_F)$p.value)

F_plot <- ggplot(toPlot_F, aes(risk_score, group)) +
  geom_density_ridges(aes(fill = group), scale = 1.5, rel_min_height = 0.001, alpha = 0.45) +
  scale_fill_manual(values = c("0" = "grey80", "1" = "#ff7f00")) +
  scale_y_discrete(labels = c("0" = "Not amplified", "1" = "Amplified")) +
  xlim(rs_range) + ridge_theme +
  labs(x = "Risk score", y = "MYCN") +
  annotate("text", x = rs_range[2] - 0.1, y = 3.3, label = p_F, hjust = 1, size = 3.5)

# ── 11. Combined PDF ──────────────────────────────────────────────────────────
pdf("G34_clinical_annotation.pdf", width = 14, height = 8)

draw(
  ht_g34,
  heatmap_legend_side    = "right",
  annotation_legend_side = "right",
  annotation_legend_list = list(lgd_tertile, lgd_subgroup, lgd_gender,
                                lgd_hist, lgd_mstage, lgd_myc, lgd_mycn),
  padding = unit(c(2, 30, 2, 2), "mm")
)

decorate_annotation("Age", {
  grid.lines(x = unit(c(0, 1), "npc"), y = unit(c(3, 3), "native"),
             gp = gpar(col = "black", lty = "dotted", lwd = 1))
})

decorate_annotation("OS", {
  grid.lines(x = unit(c(0, 1), "npc"), y = unit(c(5, 5), "native"),
             gp = gpar(col = "black", lty = "dotted", lwd = 1))
})

grid.newpage()

print(ggarrange(A, B, C, D, E, F_plot,
                ncol = 2, nrow = 3,
                labels = c("A", "B", "C", "D", "E", "F")))

dev.off()


########################################################################################################################
########################################################################################################################
#
#          7. PENALIZED COX REGRESSION (Elastic Net Risk Score) — SHH
#
########################################################################################################################
########################################################################################################################
#
# Build a multi-gene prognostic risk score using penalized Cox regression.
# Tests ridge (α=0), elastic net (α=0.5), and lasso (α=1) with 10-fold CV.
# Output: cv_fits_shh, selected genes, continuous risk scores (+ hazard scores)
#
########################################################################################################################

###########################
#1) Prepare clinical data
###########################
dat_all <- meta
dat_all <- dat_all[!is.na(dat_all$OS..years.) & !is.na(dat_all$Dead), ]
dat_all$Dead <- as.integer(dat_all$Dead)
stopifnot(all(dat_all$Dead %in% c(0, 1)))

###############################
#2) Build expression matrix
###############################
expr_df <- as.data.frame(expr)
sample_cols <- grep("^MB_SubtypeStudy_", names(expr_df), value = TRUE)
gene_ids <- expr_df[["EnsemblGeneID_from_ensemblv77"]]

X_gxs <- as.matrix(sapply(expr_df[, sample_cols, drop = FALSE], as.numeric))
rownames(X_gxs) <- gene_ids
x_all <- t(X_gxs)
rownames(x_all) <- sample_cols

######################
#3) Subset to SHH
######################
dat_shh <- subset(dat_all, Subgroup == "SHH")

############################################
#4) Align expression matrix to clinical data
############################################
id_col <- "Study_ID"
common_ids_shh <- intersect(dat_shh[[id_col]], rownames(x_all))
stopifnot(length(common_ids_shh) > 10)

dat_glm_shh <- dat_shh[match(common_ids_shh, dat_shh[[id_col]]), ]
x_glm_shh <- x_all[common_ids_shh, , drop = FALSE]

# sanity: alignment check
stopifnot(identical(dat_glm_shh[[id_col]], rownames(x_glm_shh)))

######################################################
#5) Remove genes with NA values (glmnet requirement)
######################################################
keep_genes_shh <- colSums(is.na(x_glm_shh)) == 0
x_glm_shh <- x_glm_shh[, keep_genes_shh, drop = FALSE]

############################
#6) Create survival outcome
############################
y_shh <- with(dat_glm_shh, Surv(OS..years., Dead))

######################
#7) Set up 10-fold CV
######################
set.seed(1)
K <- 10
foldid_shh <- sample(rep(1:K, length.out = nrow(x_glm_shh)))

#######################################################
#8) Fit penalized Cox models across alpha values
######################################################
alphas <- c(0, 0.5, 1)
names_alphas <- paste0("alpha_", alphas)

pboptions(type = "timer")
cat("SHH: Starting cv.glmnet |", format(Sys.time()), "\n")

cv_fits_shh <- pblapply(seq_along(alphas), function(i) {
  a <- alphas[i]
  cat("  α =", a, "started |", format(Sys.time()), "\n")
  
  fit <- cv.glmnet(
    x_glm_shh, y_shh,
    family = "cox",
    alpha = a,
    nfolds = K,
    foldid = foldid_shh,
    type.measure = "deviance",
    trace.it = 1
  )
  
  cat("  α =", a, "complete |", format(Sys.time()), "\n")
  fit
})
names(cv_fits_shh) <- names_alphas

cat("SHH: Finished cv.glmnet |", format(Sys.time()), "\n\n")

####################################
#9) Compare CV results across alpha
####################################
cv_summary_shh <- data.frame(
  alpha      = alphas,
  lambda_min = sapply(cv_fits_shh, function(f) f$lambda.min),
  lambda_1se = sapply(cv_fits_shh, function(f) f$lambda.1se),
  cvm_min    = sapply(cv_fits_shh, function(f) min(f$cvm))
)
cat("SHH: CV summary across alpha values\n")
print(cv_summary_shh)

##########################################################
#10) CV comparison curves (same y-scale)
##########################################################
fits_to_plot_shh <- list(
  ridge = cv_fits_shh[["alpha_0"]],
  enet  = cv_fits_shh[["alpha_0.5"]],
  lasso = cv_fits_shh[["alpha_1"]]
)
yl_shh <- range(unlist(lapply(fits_to_plot_shh, function(f) f$cvm)), na.rm = TRUE)

par(mfrow = c(1,3), mar = c(4,4,2,1), oma = c(0,0,3,0))
plot(fits_to_plot_shh$ridge, main = "", ylim = yl_shh); mtext("Alpha = 0 (ridge)", side=3, line=-1, cex=0.9)
plot(fits_to_plot_shh$enet,  main = "", ylim = yl_shh); mtext("Alpha = 0.5 (elastic net)", side=3, line=-1, cex=0.9)
plot(fits_to_plot_shh$lasso, main = "", ylim = yl_shh); mtext("Alpha = 1 (lasso)", side=3, line=-1, cex=0.9)
mtext("SHH: 10-fold CV curves across alpha (ridge / elastic net / lasso)", outer=TRUE, cex=1.1, line=1)
par(mfrow = c(1,1))

##########################################################
#11) Elastic net (alpha=0.5): extract genes, refit Cox
##########################################################
best_fit_shh_enet <- cv_fits_shh[["alpha_0.5"]]

lambda_min_shh <- best_fit_shh_enet$lambda.min
lambda_1se_shh <- best_fit_shh_enet$lambda.1se
lambda_min_shh; lambda_1se_shh

coef_matrix_shh <- coef(best_fit_shh_enet, s = "lambda.min")
selected_coefs_shh <- coef_matrix_shh[coef_matrix_shh[,1] != 0, , drop = FALSE]
selected_genes_shh <- rownames(selected_coefs_shh)

cat("SHH: Selected", length(selected_genes_shh), "genes at lambda.min\n")
if (length(selected_genes_shh) == 0) stop("SHH: No genes selected at lambda.min for alpha=0.5")

###############################################
#12) Refit unpenalised Cox on selected genes
###############################################
cox_data_shh <- data.frame(
  time  = dat_glm_shh$OS..years.,
  event = dat_glm_shh$Dead,
  x_glm_shh[, selected_genes_shh, drop = FALSE],
  check.names = FALSE
)

cox_model_shh <- coxph(Surv(time, event) ~ ., data = cox_data_shh)
print(summary(cox_model_shh))

################################################################################
#13) SHH: Risk score (lp) from refit Cox + validation + basic checks
################################################################################

# Coefficients from refit Cox model
final_coefs_shh <- coef(cox_model_shh)

# Risk score (linear predictor) from predict()
risk_score_shh_lp <- drop(predict(cox_model_shh, type = "lp"))

# Manual risk score using the same design matrix as coxph (validation)
X_shh <- model.matrix(cox_model_shh)
risk_score_shh_manual <- drop(X_shh[, names(final_coefs_shh), drop = FALSE] %*% final_coefs_shh)

cat("SHH cor(manual, predict lp) =", cor(risk_score_shh_manual, risk_score_shh_lp), "\n")
summary(risk_score_shh_manual - risk_score_shh_lp)

########################
# Density plot + rug + N
########################
d <- density(risk_score_shh_lp)

plot(d,
     main = "SHH risk score distribution",
     xlab = "Risk score (linear predictor)",
     ylab = "Density",
     ylim = c(0, 0.15),
     yaxt = "n")

axis(2, at = seq(0, 0.15, by = 0.05), las = 1)
rug(risk_score_shh_lp)
mtext(paste0("N = ", length(risk_score_shh_lp)), side = 1, line = 4)

# Score-only Cox model (tests association of score with survival)
summary(coxph(Surv(time, event) ~ risk_score_shh_lp, data = cox_data_shh))


################################################################################
# Prognostic Gene Analysis for SHH Medulloblastoma
################################################################################
#
# Purpose: Identify top prognostic genes in SHH-subtype medulloblastoma using
#          elastic net Cox regression, then validate with univariable Cox models
#          and visualise results as a forest plot.
#
# Output:  - Forest plot (PNG)
#          - Results table (CSV)
#          - Saved data objects (RDS)
################################################################################

#########################
# Configuration
#########################

# Use pre-existing SHH objects from workspace
clinical_shh <- dat_glm_shh
rownames(clinical_shh) <- clinical_shh$Study_ID  
expr_shh     <- x_glm_shh

# Column names in clinical data
survival_col <- "OS..years."
event_col    <- "Dead"

# Model parameters
elastic_net_alpha <- 0.5
n_top_genes       <- 10
random_seed       <- 1


#########################
# Validate inputs
#########################

# Check required columns exist
required_cols <- c(survival_col, event_col)
missing_cols  <- setdiff(required_cols, colnames(clinical_shh))

if (length(missing_cols) > 0) {
  stop("Missing columns in clinical data: ", paste(missing_cols, collapse = ", "))
}

# Ensure alignment between clinical and expression data
if (!identical(rownames(clinical_shh), rownames(expr_shh))) {
  common <- intersect(rownames(clinical_shh), rownames(expr_shh))
  clinical_shh <- clinical_shh[common, , drop = FALSE]
  expr_shh     <- expr_shh[common, , drop = FALSE]
  cat("Aligned to", length(common), "common samples\n")
}

cat("Input validation passed\n")
cat("  SHH samples:", nrow(clinical_shh), "\n")
cat("  Genes:", ncol(expr_shh), "\n\n")


#################################################
# Fit elastic net Cox model
#################################################

surv_object <- Surv(
  time  = clinical_shh[[survival_col]],
  event = clinical_shh[[event_col]]
)

cat("Fitting elastic net Cox model (alpha =", elastic_net_alpha, ")...\n")

set.seed(random_seed)
cv_fit <- cv.glmnet(
  x      = as.matrix(expr_shh),
  y      = surv_object,
  family = "cox",
  alpha  = elastic_net_alpha,
  nfolds = 10
)

cat("  lambda.min:", round(cv_fit$lambda.min, 4), "\n")
cat("  lambda.1se:", round(cv_fit$lambda.1se, 4), "\n")


#################################################
# Extract top genes by coefficient magnitude
#################################################

coefficients <- as.matrix(coef(cv_fit, s = "lambda.min"))
nonzero_coef <- coefficients[coefficients[, 1] != 0, , drop = FALSE]

cat("\nNon-zero coefficients:", nrow(nonzero_coef), "\n")

if (nrow(nonzero_coef) == 0) {
  stop("No genes selected by elastic net. Try adjusting alpha or lambda.")
}

# Rank by absolute coefficient value
ranked_genes <- rownames(nonzero_coef)[order(abs(nonzero_coef[, 1]), decreasing = TRUE)]
top_genes    <- head(ranked_genes, n_top_genes)


#################################################
# Map Ensembl IDs to gene symbols
#################################################

map_to_symbols <- function(ensembl_ids) {
  # Strip version numbers (e.g., ENSG00000123456.7 -> ENSG00000123456)
  clean_ids <- gsub("\\..*", "", ensembl_ids)
  
  symbols <- mapIds(
    org.Hs.eg.db,
    keys      = clean_ids,
    column    = "SYMBOL",
    keytype   = "ENSEMBL",
    multiVals = "first"
  )
  
  # Fall back to Ensembl ID if no symbol found
  symbols <- as.character(symbols)
  symbols[is.na(symbols)] <- ensembl_ids[is.na(symbols)]
  
  return(symbols)
}

gene_symbols <- map_to_symbols(top_genes)


#################################################
# Univariable Cox regression for each top gene
#################################################

run_univariable_cox <- function(gene_id, gene_symbol) {
  
  analysis_df <- data.frame(
    time       = clinical_shh[[survival_col]],
    status     = clinical_shh[[event_col]],
    expression = as.numeric(expr_shh[, gene_id])
  )
  analysis_df <- na.omit(analysis_df)
  
  cox_model   <- coxph(Surv(time, status) ~ expression, data = analysis_df)
  cox_summary <- summary(cox_model)
  
  data.frame(
    gene_id     = gene_id,
    symbol      = gene_symbol,
    glmnet_beta = nonzero_coef[gene_id, 1],
    HR          = cox_summary$coef[1, "exp(coef)"],
    lower_95    = cox_summary$conf.int[1, "lower .95"],
    upper_95    = cox_summary$conf.int[1, "upper .95"],
    p_value     = cox_summary$coef[1, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
}

cat("Running univariable Cox models for top", length(top_genes), "genes...\n")

results <- do.call(rbind, Map(run_univariable_cox, top_genes, gene_symbols))
results$p_adjusted <- p.adjust(results$p_value, method = "BH")
results <- results[order(results$HR), ]

cat("\nResults:\n")
print(results, row.names = FALSE)


#################################################
# Forest plot SHH Cox PH top genes — gene names in italics
#################################################

cat("\nGenerating forest plot...\n")

png("SHH_forest_plot.png", width = 10, height = 8, units = "in", res = 150)

par(mar = c(5, 10, 4, 8))

y_positions <- seq_len(nrow(results))
x_range     <- range(c(results$lower_95, results$upper_95, 1), na.rm = TRUE)

plot(
  results$HR, y_positions,
  xlim = x_range,
  xlab = "Hazard Ratio (95% CI)",
  ylab = "",
  yaxt = "n",
  log  = "x",
  pch  = 16,
  cex  = 1.2,
  main = sprintf("SHH Medulloblastoma: Top %d Prognostic Genes (Elastic Net \u03b1=%.1f)",
                 nrow(results), elastic_net_alpha)
)

# Confidence interval whiskers
arrows(
  results$lower_95, y_positions,
  results$upper_95, y_positions,
  angle = 90, code = 3, length = 0.04
)

# Gene labels on y-axis — italic (font = 3)
axis(2, at = y_positions, labels = results$symbol,
     las = 1, cex.axis = 0.9, font = 3)

# Reference line at HR = 1
abline(v = 1, lty = 2, col = "grey40")

# HR and p-value annotations
text(
  results$upper_95, y_positions,
  labels = sprintf("HR = %.2f (p = %.1e)", results$HR, results$p_value),
  pos = 4, cex = 0.7, xpd = NA
)

dev.off()

##############################################################################
# Stratify SHH patients into risk tertiles (Low / Mid / High)
##############################################################################

rs_shh <- as.numeric(risk_score_shh_lp)
cuts_shh <- quantile(rs_shh, probs = c(1/3, 2/3), na.rm = TRUE, type = 7)

dat_glm_shh$risk_score   <- rs_shh
dat_glm_shh$risk_tertile <- cut(
  rs_shh,
  breaks         = c(-Inf, cuts_shh[1], cuts_shh[2], Inf),
  labels         = c("Low", "Mid", "High"),
  include.lowest = TRUE,
  right          = TRUE
)
stopifnot("Risk score length doesn't match dat_glm_shh rows" = length(rs_shh) == nrow(dat_glm_shh))

message("SHH tertile sizes:");       print(table(dat_glm_shh$risk_tertile))
message("Cutpoints (33rd / 67th):"); print(cuts_shh)
message("Score range per tertile:"); print(tapply(dat_glm_shh$risk_score, dat_glm_shh$risk_tertile, range))
message("Event counts by tertile:"); print(with(dat_glm_shh, table(risk_tertile, Dead)))

###########################################################################
# KM curves, one per risk tertile — SHH cohort + p-value
###########################################################################
km_fit_shh  <- survfit(Surv(OS..years., Dead) ~ risk_tertile, data = dat_glm_shh)
lr_test_shh <- survdiff(Surv(OS..years., Dead) ~ risk_tertile, data = dat_glm_shh)
lr_p_shh    <- 1 - pchisq(lr_test_shh$chisq, df = length(lr_test_shh$n) - 1)
cindex_shh  <- 0.932

ggsurvplot(
  km_fit_shh,
  data              = dat_glm_shh,
  palette           = c("#2166ac", "#f4a582", "#d6604d"),
  conf.int          = TRUE,
  conf.int.alpha    = 0.12,
  risk.table        = TRUE,
  risk.table.height = 0.25,
  xlab              = "Time (years)",
  ylab              = "Overall survival probability",
  title             = "SHH cohort: KM by risk tertile",
  legend.title      = "Risk group",
  legend.labs       = c("Low", "Mid", "High"),
  pval              = FALSE,
  break.time.by     = 2,
  ggtheme           = theme_classic(base_size = 13),
  xlim              = c(0, 10)
) |>
  (\(p) {
    p$plot <- p$plot +
      annotate("text",
               x     = 0.3, y = 0.20,
               label = paste0("p ", format.pval(lr_p_shh, digits = 2, eps = 0.001)),
               size  = 4, hjust = 0) +
      annotate("text",
               x     = 0.3, y = 0.12,
               label = paste0("C-index = ", cindex_shh, " (SE = 0.021)"),
               size  = 4, hjust = 0)
    p
  })()

###########################################################################
# SHH cohort — Clinical annotation figure
# Page 1: Oncoprint ordered by ascending Cox risk score
# Page 2: Ridge plots of risk score by clinical variable
###########################################################################

# ── 1. Merge MYC and MYCN amplification into dat_glm_shh
dat_glm_shh$MYC  <- ifelse(dat_glm_shh$Study_ID %in% myc_amp_ids,  1, 0)
dat_glm_shh$MYCN <- ifelse(dat_glm_shh$Study_ID %in% mycn_amp_ids, 1, 0)

cat("MYC amplification in SHH:\n");  print(table(dat_glm_shh$MYC,  useNA = "ifany"))
cat("MYCN amplification in SHH:\n"); print(table(dat_glm_shh$MYCN, useNA = "ifany"))

# ── 2. Prepare data for oncoprint
dat_plot      <- dat_glm_shh[order(dat_glm_shh$risk_score), ]
dat_plot$Dead <- as.numeric(dat_plot$Dead)

# ── 3. Clean NAs
dat_plot$Gender    <- ifelse(is.na(dat_plot$Gender),    "Unknown", as.character(dat_plot$Gender))
dat_plot$histology <- ifelse(is.na(dat_plot$histology) | dat_plot$histology == "", "Unknown", dat_plot$histology)
dat_plot$Met.status..1.Met..0.M0. <- ifelse(
  is.na(dat_plot$Met.status..1.Met..0.M0.), "Unknown",
  as.character(dat_plot$Met.status..1.Met..0.M0.)
)
dat_plot$MYC  <- as.character(dat_plot$MYC)
dat_plot$MYCN <- as.character(dat_plot$MYCN)

# ── 4. Colour palettes
tertile_col_fun <- c("Low" = "#2166ac", "Mid" = "#f4a582", "High" = "#d6604d")
gender_col_fun  <- c("M" = "white", "F" = "black", "Unknown" = "lightgrey")
hist_col_fun    <- c("Classic" = "white", "Desmoplastic" = "#fc8d59", "LCA" = "black", "MBEN" = "red4", "Unknown" = "lightgrey")
mstage_col_fun  <- c("0" = "white", "1" = "black", "Unknown" = "lightgrey")
mycn_col_fun    <- c("0" = "white", "1" = "#ff7f00")

# ── 5. Age dot annotation (capped at 20)
Age_capped <- pmin(dat_plot$Age, 20)

age_dot_anno <- HeatmapAnnotation(
  Age = anno_points(
    Age_capped,
    ylim = c(0, 20),
    pch = ifelse(is.na(dat_plot$Age), 4, 20),
    gp = gpar(col = ifelse(is.na(dat_plot$Age), "red", "grey")),
    axis_param = list(at = c(0, 5, 10, 15, 20)),
    baseline = 3,
    baseline_gp = gpar(col = "grey", lty = 2, lwd = 1)
  ),
  annotation_name_side = "left",
  annotation_height = unit(2, "cm")
)

# ── 6. OS dot annotation (capped at 10 years)
OS <- pmin(dat_plot$OS..years., 10)
pchSetup <- ifelse(dat_plot$Dead == 1 & dat_plot$OS..years. < 10, 4,
                   ifelse(dat_plot$Dead == 0 & dat_plot$OS..years. < 10, 20, 2))

OS_dot_anno <- HeatmapAnnotation(
  OS = anno_points(
    OS,
    pch = pchSetup,
    gp = gpar(col = ifelse(dat_plot$Dead == 1, "firebrick3", "lightgrey")),
    ylim = c(0, 10),
    axis_param = list(at = c(0, 5, 10)),
    baseline_gp = gpar(col = "lightgrey", lty = 2, lwd = 1)
  ),
  annotation_name_side = "left",
  annotation_height = unit(2, "cm")
)

# ── 7. Categorical bar annotations — Subtype removed
cat_anno <- HeatmapAnnotation(
  `Risk tertile` = anno_simple(as.character(dat_plot$risk_tertile), col = tertile_col_fun, border = TRUE),
  Gender         = anno_simple(dat_plot$Gender,                     col = gender_col_fun,  border = TRUE),
  Histology      = anno_simple(dat_plot$histology,                  col = hist_col_fun,    border = TRUE),
  `M-stage`      = anno_simple(dat_plot$Met.status..1.Met..0.M0.,  col = mstage_col_fun,  border = TRUE),
  MYCN           = anno_simple(dat_plot$MYCN,                       col = mycn_col_fun,    border = TRUE),
  annotation_name_side = "left"
)

# ── 8. Base heatmap
rs_matrix <- matrix(dat_plot$risk_score, nrow = 1)

ht_shh <- Heatmap(
  rs_matrix,
  name  = "Risk score",
  col   = colorRamp2(
    c(min(dat_plot$risk_score), 0, max(dat_plot$risk_score)),
    c("#2166ac", "white", "#d6604d")
  ),
  border              = TRUE,
  show_column_names   = FALSE,
  cluster_columns     = FALSE,
  cluster_rows        = FALSE,
  row_names_side      = "left",
  row_names_gp        = gpar(fontsize = 10),
  height              = unit(0.6, "cm"),
  show_heatmap_legend = TRUE,
  top_annotation      = c(cat_anno, age_dot_anno, OS_dot_anno)
)

# ── 9. Manual legends — lgd_subtype removed
lgd_tertile <- Legend(labels = c("Low", "Mid", "High"),
                      title = "Risk tertile",
                      legend_gp = gpar(fill = c("#2166ac", "#f4a582", "#d6604d")), border = TRUE)
lgd_gender  <- Legend(labels = c("Male", "Female", "Unknown"),
                      title = "Gender",
                      legend_gp = gpar(fill = c("white", "black", "lightgrey")), border = TRUE)
lgd_hist    <- Legend(labels = c("Classic", "Desmoplastic", "LCA", "MBEN", "Unknown"),
                      title = "Histology",
                      legend_gp = gpar(fill = c("white", "#fc8d59", "black", "red4", "lightgrey")), border = TRUE)
lgd_mstage  <- Legend(labels = c("M0", "M+", "Unknown"),
                      title = "M-stage",
                      legend_gp = gpar(fill = c("white", "black", "lightgrey")), border = TRUE)
lgd_mycn    <- Legend(labels = c("Not amplified", "Amplified"),
                      title = "MYCN",
                      legend_gp = gpar(fill = c("white", "#ff7f00")), border = TRUE)

# ── 10. Ridge plots
fmt_p <- function(p) {
  if (p < 2.2e-16)    "p < 2.2e-16"
  else if (p < 0.001) paste0("p = ", formatC(p, format = "e", digits = 2))
  else                paste0("p = ", round(p, 3))
}

rs_range    <- range(dat_glm_shh$risk_score, na.rm = TRUE)
ridge_theme <- theme_bw() + theme(legend.position = "none", panel.grid = element_blank())

# A - Histology
toPlot_A <- data.frame(risk_score = dat_glm_shh$risk_score, group = dat_glm_shh$histology) %>%
  filter(!is.na(group), !group %in% c("Unknown", ""))
p_A <- fmt_p(anova(aov(risk_score ~ group, data = toPlot_A))$`Pr(>F)`[1])
hist_cols <- c("Classic" = "grey80", "Desmoplastic" = "#fc8d59", "LCA" = "black", "MBEN" = "red4")

A <- ggplot(toPlot_A, aes(risk_score, group)) +
  geom_density_ridges(aes(fill = group), scale = 1.5, rel_min_height = 0.001, alpha = 0.45) +
  scale_fill_manual(values = hist_cols) +
  xlim(rs_range) + ridge_theme +
  labs(x = "Risk score", y = "Histology") +
  annotate("text", x = rs_range[2] - 0.1, y = length(unique(toPlot_A$group)) + 1.3, label = p_A, hjust = 1, size = 3.5)

# B - M-stage
toPlot_B <- data.frame(risk_score = dat_glm_shh$risk_score,
                       group = as.character(dat_glm_shh$Met.status..1.Met..0.M0.)) %>%
  filter(!is.na(group), group != "Unknown")
p_B <- fmt_p(t.test(risk_score ~ group, data = toPlot_B)$p.value)

B <- ggplot(toPlot_B, aes(risk_score, group)) +
  geom_density_ridges(aes(fill = group), scale = 1.5, rel_min_height = 0.001, alpha = 0.45) +
  scale_fill_manual(values = c("0" = "grey80", "1" = "black")) +
  scale_y_discrete(labels = c("0" = "M0", "1" = "M+")) +
  xlim(rs_range) + ridge_theme +
  labs(x = "Risk score", y = "M-stage") +
  annotate("text", x = rs_range[2] - 0.1, y = 3.3, label = p_B, hjust = 1, size = 3.5)

# C - Gender
toPlot_C <- data.frame(risk_score = dat_glm_shh$risk_score, group = dat_glm_shh$Gender) %>%
  filter(!is.na(group), group != "Unknown")
p_C <- fmt_p(t.test(risk_score ~ group, data = toPlot_C)$p.value)

C <- ggplot(toPlot_C, aes(risk_score, group)) +
  geom_density_ridges(aes(fill = group), scale = 1.5, rel_min_height = 0.001, alpha = 0.45) +
  scale_fill_manual(values = c("F" = "#d73027", "M" = "#4575b4")) +
  xlim(rs_range) + ridge_theme +
  labs(x = "Risk score", y = "Gender") +
  annotate("text", x = rs_range[2] - 0.1, y = 3.3, label = p_C, hjust = 1, size = 3.5)

# D - MYCN amplification
toPlot_D <- data.frame(risk_score = dat_glm_shh$risk_score,
                       group = as.character(dat_glm_shh$MYCN)) %>%
  filter(!is.na(group))
p_D <- fmt_p(t.test(risk_score ~ group, data = toPlot_D)$p.value)

D <- ggplot(toPlot_D, aes(risk_score, group)) +
  geom_density_ridges(aes(fill = group), scale = 1.5, rel_min_height = 0.001, alpha = 0.45) +
  scale_fill_manual(values = c("0" = "grey80", "1" = "#ff7f00")) +
  scale_y_discrete(labels = c("0" = "Not amplified", "1" = "Amplified")) +
  xlim(rs_range) + ridge_theme +
  labs(x = "Risk score", y = "MYCN") +
  annotate("text", x = rs_range[2] - 0.1, y = 3.3, label = p_D, hjust = 1, size = 3.5)

# ── 11. Combined PDF
pdf("SHH_clinical_annotation.pdf", width = 14, height = 8)

draw(
  ht_shh,
  heatmap_legend_side    = "right",
  annotation_legend_side = "right",
  annotation_legend_list = list(lgd_tertile, lgd_gender, lgd_hist, lgd_mstage, lgd_mycn),
  padding = unit(c(2, 30, 2, 2), "mm")
)

decorate_annotation("Age", {
  grid.lines(x = unit(c(0, 1), "npc"), y = unit(c(3, 3), "native"),
             gp = gpar(col = "black", lty = "dotted", lwd = 1))
})

decorate_annotation("OS", {
  grid.lines(x = unit(c(0, 1), "npc"), y = unit(c(5, 5), "native"),
             gp = gpar(col = "black", lty = "dotted", lwd = 1))
})

grid.newpage()

print(ggarrange(A, B, C, D,
                ncol = 2, nrow = 2,
                labels = c("A", "B", "C", "D")))

dev.off()

########################################################################################################################
########################################################################################################################
#
#    8. METHYLOMIC VS TRANSCRIPTOMIC RISK SCORE COMPARISON
#       Ewan and Edward: scatter plots for G3/4 and SHH
#
########################################################################################################################
########################################################################################################################

# 1) Load data
methyl_g34  <- read.csv("EdGroup3_4_alpha0.5_RiskScores.csv")
methyl_shh  <- read.csv("EdSHH_alpha1_RiskScores.csv")
metadata    <- read.csv("taylor_pheno_forDan.csv")

# 2) Tidy transcriptomic scores
transcr_g34 <- read.csv("G3G4_predictor.csv", row.names = 1) %>%
  rownames_to_column("Study_ID") %>%
  setNames(c("Study_ID", "transcr_score"))

transcr_shh <- read.csv("SHH_predictor.csv", row.names = 1) %>%
  rownames_to_column("Study_ID") %>%
  setNames(c("Study_ID", "transcr_score"))

# 3) Build metadata lookup: GSM -> Study_ID + Subgroup
meta_lookup <- metadata %>%
  dplyr::select(title, Study_ID, Subgroup) %>%
  dplyr::rename(GSM = title)

# 4) Strip array suffix from methylomic SampleIDs
methyl_g34 <- methyl_g34 %>% mutate(GSM = sub("_.*", "", SampleID))
methyl_shh <- methyl_shh %>% mutate(GSM = sub("_.*", "", SampleID))

# 5) Merge
merged_g34 <- methyl_g34 %>%
  dplyr::rename(methyl_score = RiskScore, methyl_risk_group = RiskGroup) %>%
  left_join(meta_lookup, by = "GSM") %>%
  left_join(transcr_g34, by = "Study_ID") %>%
  drop_na(methyl_score, transcr_score)

merged_shh <- methyl_shh %>%
  dplyr::rename(methyl_score = RiskScore, methyl_risk_group = RiskGroup) %>%
  left_join(meta_lookup, by = "GSM") %>%
  left_join(transcr_shh, by = "Study_ID") %>%
  drop_na(methyl_score, transcr_score)

# 6) Check matched samples
cat("G3/4 matched samples:", nrow(merged_g34), "\n")
cat("SHH  matched samples:", nrow(merged_shh), "\n")

################################################################################
# G3/4 SCATTER PLOT
################################################################################

ct_g34 <- cor.test(merged_g34$methyl_score, merged_g34$transcr_score, method = "pearson")
r_g34  <- unname(ct_g34$estimate)
p_g34  <- ct_g34$p.value

p1 <- ggplot(merged_g34, aes(x = methyl_score, y = transcr_score, colour = Subgroup)) +
  geom_point(alpha = 0.75, size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey40",
              linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = Inf, y = Inf,
           label = paste0(
             "r = ", round(r_g34, 3),
             "\np = ", format.pval(p_g34, digits = 2, eps = 1e-3),
             "\nn = ", nrow(merged_g34)
           ),
           hjust = 1.1, vjust = 1.5, size = 4) +
  scale_colour_manual(values = c("Group3" = "#FFD700", "Group4" = "#008000")) +
  labs(title  = "G3/4: Methylomic vs Transcriptomic Risk Scores",
       x      = "Methylomic Risk Score",
       y      = "Transcriptomic Risk Score",
       colour = "Subgroup") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

################################################################################
# SHH SCATTER PLOT
################################################################################

ct_shh <- cor.test(merged_shh$methyl_score, merged_shh$transcr_score, method = "pearson")
r_shh  <- unname(ct_shh$estimate)
p_shh  <- ct_shh$p.value

p2 <- ggplot(merged_shh, aes(x = methyl_score, y = transcr_score, colour = Subgroup)) +
  geom_point(alpha = 0.75, size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey40",
              linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = Inf, y = Inf,
           label = paste0(
             "r = ", round(r_shh, 3),
             "\np = ", format.pval(p_shh, digits = 2, eps = 1e-3),
             "\nn = ", nrow(merged_shh)
           ),
           hjust = 1.1, vjust = 1.5, size = 4) +
  scale_colour_manual(values = c("SHH" = "#FF6666")) +
  labs(title  = "SHH: Methylomic vs Transcriptomic Risk Scores",
       x      = "Methylomic Risk Score",
       y      = "Transcriptomic Risk Score",
       colour = "Subgroup") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

print(p1 + p2)

################################################################################
# SHH SCATTER PLOT
################################################################################

ct_shh <- cor.test(merged_shh$methyl_score, merged_shh$transcr_score, method = "pearson")
r_shh  <- unname(ct_shh$estimate)
p_shh  <- ct_shh$p.value

p2 <- ggplot(merged_shh, aes(x = methyl_score, y = transcr_score, colour = Subgroup)) +
  geom_point(alpha = 0.75, size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey40",
              linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = Inf, y = Inf,
           label = paste0(
             "r = ", round(r_shh, 3),
             "\np = ", format.pval(p_shh, digits = 2, eps = 1e-3),
             "\nn = ", nrow(merged_shh)
           ),
           hjust = 1.1, vjust = 1.5, size = 4) +
  scale_colour_manual(values = c("SHH" = "#FF6666")) +
  labs(
    title  = "SHH: Methylomic vs Transcriptomic Risk Scores",
    x      = "Methylomic Risk Score",
    y      = "Transcriptomic Risk Score",
    colour = "Subgroup"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

print(p2)


########################################################################################################################
########################################################################################################################
#
#    9. RNA-SEQ CROSS-PLATFORM VALIDATION (E-MTAB-10767)
#       Microarray-derived risk signatures projected onto RNA-seq cohort
#       Subgroups: Group 3/4 and SHH
#
########################################################################################################################
########################################################################################################################

# ── 1. Load the RNA-seq expression matrix ─────────────────────────────────────
expr <- fread("MB.vsd.txt", data.table = FALSE)

gene_ids  <- expr[[1]]
expr[[1]] <- NULL
expr[]    <- lapply(expr, function(x) as.numeric(x))
expr_mat  <- as.matrix(expr)
rownames(expr_mat) <- gene_ids

rownames(expr_mat) <- sub("_.*$",   "", rownames(expr_mat))
rownames(expr_mat) <- sub("\\..*$", "", rownames(expr_mat))

if (any(duplicated(rownames(expr_mat)))) {
  expr_dt  <- as.data.table(expr_mat, keep.rownames = "ENSG")
  expr_dt  <- expr_dt[, lapply(.SD, mean, na.rm = TRUE), by = ENSG]
  expr_mat <- as.matrix(expr_dt[, -1])
  rownames(expr_mat) <- expr_dt$ENSG
}

cat("Expression matrix:", nrow(expr_mat), "genes x", ncol(expr_mat), "samples\n")

# ── 2. Load clinical metadata ─────────────────────────────────────────────────
sdrf        <- fread("E-MTAB-10767.sdrf (1).txt", data.table = FALSE)
names(sdrf) <- make.unique(names(sdrf))

pheno <- sdrf %>%
  transmute(
    SampleID   = as.character(`Source Name`),
    Subgroup   = as.character(`Characteristics[subgroup]`),
    time       = as.numeric(`Characteristics[os time]`),
    status_raw = as.character(`Characteristics[os status]`)
  ) %>%
  mutate(event = as.numeric(status_raw)) %>%
  filter(!is.na(SampleID), !is.na(time), !is.na(event), time > 0)

# ── 3. Subset to G3/4 and build x_glm ────────────────────────────────────────
dat_g34        <- pheno %>% filter(Subgroup %in% c("Grp3", "Grp4", "Grp3/Grp4"))
common_ids_g34 <- dat_g34$SampleID[dat_g34$SampleID %in% colnames(expr_mat)]
dat_glm        <- dat_g34[match(common_ids_g34, dat_g34$SampleID), ]

x_glm <- t(expr_mat[, common_ids_g34, drop = FALSE])
rownames(x_glm) <- common_ids_g34
colnames(x_glm) <- rownames(expr_mat)

sd_vec <- apply(x_glm, 2, sd, na.rm = TRUE)
x_glm  <- x_glm[, sd_vec > 0, drop = FALSE]
cat("G3/4 samples:", nrow(x_glm), "\n")
cat("Genes after removing zero variance:", ncol(x_glm), "\n")

y_g34 <- with(dat_glm, Surv(time, event))


###########################################################################
# SECTION 1: INITIAL EXPLORATION (RETAINED FOR REFERENCE)
# ─────────────────────────────────────────────────────────────────────────
# This block was used to explore whether predict.coxph() could be applied
# directly across platforms. It was retained to document the initial
# approach and illustrate the scale mismatch between microarray and RNA-seq
# risk score distributions. The density plot shows both distributions
# overlaid — the shift confirms platforms cannot share cutpoints directly.
# The manual dot product approach (Section 2) supersedes this.
###########################################################################

# Microarray risk scores (linear predictor from training model)
risk_score_discovery.g34.microarray <- predict(cox_model, type = "lp")

# Attempted direct application to RNA-seq via predict.coxph()
all(names(final_coefs) %in% colnames(x_glm))
risk_score_discovery.g34.rnaseq <- predict(cox_model, type = "lp",
                                           newdata = data.frame(x_glm))

# Distribution comparison -- illustrates platform scale mismatch
plot(xlim = c(-12, 6), density(risk_score_discovery.g34.microarray))
lines(col = "cornflowerblue", density(risk_score_discovery.g34.rnaseq))


###########################################################################
# SECTION 2: G3/4 VALIDATION — MICROARRAY SIGNATURE PROJECTED ONTO RNA-SEQ
###########################################################################

# ── 2.1 Gene overlap check ────────────────────────────────────────────────────
genes_in_model <- names(final_coefs)

overlap <- intersect(genes_in_model, colnames(x_glm))
missing <- setdiff(genes_in_model, colnames(x_glm))

cat("Genes in model:", length(genes_in_model), "\n")
cat("Found in RNA-seq:", length(overlap), "\n")
cat("Missing:", length(missing), "\n")

# ── 2.2 Compute RNA-seq risk scores (manual dot product) ─────────────────────
coefs_matched <- final_coefs[overlap]
x_glm_matched <- x_glm[, overlap, drop = FALSE]

risk_score_discovery.g34.rnaseq <- as.numeric(x_glm_matched %*% coefs_matched)
names(risk_score_discovery.g34.rnaseq) <- rownames(x_glm)

# ── 2.3 Build RNA-seq G3/4 clinical metadata ──────────────────────────────────
sdrf        <- fread("E-MTAB-10767.sdrf (1).txt", data.table = FALSE)
names(sdrf) <- make.unique(names(sdrf))

pheno_rnaseq <- sdrf %>%
  transmute(
    SampleID   = as.character(`Source Name`),
    Subgroup   = as.character(`Characteristics[subgroup]`),
    time       = as.numeric(`Characteristics[os time]`),
    status_raw = as.character(`Characteristics[os status]`)
  ) %>%
  mutate(event = as.numeric(status_raw)) %>%
  filter(Subgroup %in% c("Grp3", "Grp4"),
         !is.na(time), !is.na(event), time > 0,
         SampleID %in% names(risk_score_discovery.g34.rnaseq))

# ── 2.4 Attach risk scores and RNA-seq-derived tertile groups ─────────────────
# NOTE: Tertile cutpoints derived from RNA-seq distribution, not microarray.
# Platforms operate on different numerical scales -- applying microarray
# cutpoints to RNA-seq scores assigns all samples to a single group.
pheno_rnaseq$risk_score <- risk_score_discovery.g34.rnaseq[pheno_rnaseq$SampleID]

tertile_cuts_rnaseq <- quantile(risk_score_discovery.g34.rnaseq, probs = c(1/3, 2/3))

pheno_rnaseq$risk_group <- cut(
  pheno_rnaseq$risk_score,
  breaks = c(-Inf, tertile_cuts_rnaseq, Inf),
  labels = c("Low", "Mid", "High")
)

table(pheno_rnaseq$risk_group)

# ── 2.5 Kaplan-Meier validation plot ─────────────────────────────────────────
km_fit <- survfit(Surv(time, event) ~ risk_group, data = pheno_rnaseq)

ggsurvplot(km_fit, data = pheno_rnaseq,
           pval             = TRUE,
           pval.coord       = c(0.5, 0.15),
           conf.int         = FALSE,
           risk.table       = TRUE,
           xlim             = c(0, 12),
           break.time.by    = 2,
           palette          = c("cornflowerblue", "orange", "red"),
           legend.labs      = c("Low", "Mid", "High"),
           legend.title     = "Risk group",
           xlab             = "Time (years)",
           ylab             = "Overall survival probability",
           title            = "RNA-seq G3/4 \u2014 Microarray-derived risk score validation",
           ggtheme          = theme_classic(),
           risk.table.col   = "strata",
           risk.table.title = "Number at risk",
           risk.table.y.text = TRUE,
           tables.theme     = theme_cleantable(),
           fontsize         = 4)

# ── 2.6 Build risk_out_rnaseq for oncoprint ───────────────────────────────────
risk_out_rnaseq <- pheno_rnaseq %>%
  mutate(
    risk_lp = risk_score,
    tertile = risk_group
  ) %>%
  dplyr::select(SampleID, risk_lp, event, time, tertile)


###########################################################################
# SECTION 3: G3/4 CLINICAL ANNOTATION FIGURE
# Oncoprint (Page 1) + Ridge plots (Page 2)
# Input: risk_out_rnaseq (microarray-derived risk scores projected onto RNA-seq)
###########################################################################

names(sdrf) <- make.unique(names(sdrf))

# ── 3.1 Standardise tertile labels ───────────────────────────────────────────
risk_out_rnaseq$tertile <- as.character(risk_out_rnaseq$tertile)
risk_out_rnaseq$tertile[risk_out_rnaseq$tertile == "Mid"] <- "Intermediate"
risk_out_rnaseq$tertile <- factor(risk_out_rnaseq$tertile,
                                  levels = c("Low", "Intermediate", "High"))

# ── 3.2 Pull clinical metadata from SDRF ─────────────────────────────────────
clinical_rnaseq <- sdrf %>%
  transmute(
    SampleID   = as.character(`Source Name`),
    Subgroup   = as.character(`Characteristics[subgroup]`),
    Sex_raw    = as.character(`Characteristics[sex]`),
    Age_raw    = as.character(`Characteristics[age]`),
    Mstage_raw = as.character(`Characteristics[m+]`),
    MYC_raw    = as.character(`Characteristics[myc amplification]`),
    MYCN_raw   = as.character(`Characteristics[mycn amplification]`)
  ) %>%
  filter(Subgroup %in% c("Grp3", "Grp4", "Grp3/Grp4")) %>%
  mutate(
    Gender = ifelse(tolower(trimws(Sex_raw)) %in% c("m","male"),   "M",
                    ifelse(tolower(trimws(Sex_raw)) %in% c("f","female"), "F", "Unknown")),
    Age = local({
      x <- tolower(trimws(as.character(Age_raw)))
      ifelse(grepl("^0|0.*3|infant|toddler", x),                   "0-3",
             ifelse(grepl("3.*16|3-16|child|pediatric|paediatric", x),    "3-16",
                    ifelse(grepl(">.*16|adult|>16", x),                          ">16",
                           ifelse(x == "" | x == "na" | is.na(x),                       "Unknown", x))))
    }),
    Mstage = ifelse(tolower(trimws(Mstage_raw)) %in% c("true","1","yes","y"), "1",
                    ifelse(tolower(trimws(Mstage_raw)) %in% c("false","0","no","n"), "0", "Unknown")),
    MYC    = ifelse(tolower(trimws(MYC_raw))  %in% c("true","1","yes","y"),   "1", "0"),
    MYCN   = ifelse(tolower(trimws(MYCN_raw)) %in% c("true","1","yes","y"),   "1", "0")
  )

# ── 3.3 Merge and order by ascending risk score ───────────────────────────────
dat_plot <- risk_out_rnaseq %>%
  left_join(clinical_rnaseq, by = "SampleID") %>%
  arrange(risk_lp)

cat("Samples:", nrow(dat_plot), "\n")
cat("Raw Age_raw examples:\n"); print(head(unique(clinical_rnaseq$Age_raw), 10))
cat("Parsed Age examples:\n");  print(head(dat_plot$Age[!is.na(dat_plot$Age)], 10))

# ── 3.4 Clean NAs ────────────────────────────────────────────────────────────
dat_plot$Gender   <- ifelse(is.na(dat_plot$Gender)   | dat_plot$Gender   == "", "Unknown", dat_plot$Gender)
dat_plot$Subgroup <- ifelse(is.na(dat_plot$Subgroup) | dat_plot$Subgroup == "", "Unknown", dat_plot$Subgroup)
dat_plot$Mstage   <- ifelse(is.na(dat_plot$Mstage)   | dat_plot$Mstage   == "", "Unknown", dat_plot$Mstage)
dat_plot$Subgroup <- dplyr::recode(dat_plot$Subgroup,
                                   "Grp3"      = "Group 3",
                                   "Grp4"      = "Group 4",
                                   "Grp3/Grp4" = "Unknown",
                                   "Unknown"   = "Unknown")
dat_plot$MYC     <- as.character(dat_plot$MYC)
dat_plot$MYCN    <- as.character(dat_plot$MYCN)
dat_plot$tertile <- ifelse(is.na(dat_plot$tertile) | dat_plot$tertile == "", "Unknown",
                           as.character(dat_plot$tertile))

# ── 3.5 Colour palettes ───────────────────────────────────────────────────────
tertile_col_fun  <- c("Low" = "#2166ac", "Intermediate" = "#f4a582", "High" = "#d6604d", "Unknown" = "lightgrey")
gender_col_fun   <- c("M" = "white",     "F" = "black",              "Unknown" = "lightgrey")
subgroup_col_fun <- c("Group 3" = "#FFD700", "Group 4" = "#008000",  "Unknown" = "lightgrey")
mstage_col_fun   <- c("0" = "white",     "1" = "black",              "Unknown" = "lightgrey")
myc_col_fun      <- c("0" = "white",     "1" = "#e31a1c")
mycn_col_fun     <- c("0" = "white",     "1" = "#ff7f00")

# ── 3.6 Age group ─────────────────────────────────────────────────────────────
dat_plot$AgeGroup <- ifelse(is.na(dat_plot$Age) | dat_plot$Age == "", "Unknown", dat_plot$Age)
dat_plot <- dat_plot %>% filter(!AgeGroup %in% c(">16", "Unknown"))
cat("Samples after excluding >16 age group:", nrow(dat_plot), "\n")
age_col_fun <- c("0-3" = "#ffffcc", "3-16" = "#41b6c4")

# ── 3.7 OS dot annotation (capped at 10 years) ───────────────────────────────
OS       <- pmin(dat_plot$time, 10)
pchSetup <- ifelse(dat_plot$event == 1 & dat_plot$time <  10,  4,
                   ifelse(dat_plot$event == 0 & dat_plot$time <  10, 20, 2))

OS_dot_anno <- HeatmapAnnotation(
  OS = anno_points(
    OS,
    pch         = pchSetup,
    gp          = gpar(col = ifelse(dat_plot$event == 1, "firebrick3", "lightgrey")),
    ylim        = c(0, 10),
    axis_param  = list(at = c(0, 5, 10)),
    baseline_gp = gpar(col = "lightgrey", lty = 2, lwd = 1)
  ),
  annotation_name_side = "left",
  annotation_height    = unit(2, "cm")
)

# ── 3.8 Categorical bar annotations ──────────────────────────────────────────
cat_anno <- HeatmapAnnotation(
  `Risk tertile` = anno_simple(dat_plot$tertile,   col = tertile_col_fun,  border = TRUE),
  Gender         = anno_simple(dat_plot$Gender,    col = gender_col_fun,   border = TRUE),
  `Age group`    = anno_simple(dat_plot$AgeGroup,  col = age_col_fun,      border = TRUE),
  Subgroup       = anno_simple(dat_plot$Subgroup,  col = subgroup_col_fun, border = TRUE),
  `M-stage`      = anno_simple(dat_plot$Mstage,    col = mstage_col_fun,   border = TRUE),
  MYC            = anno_simple(dat_plot$MYC,       col = myc_col_fun,      border = TRUE),
  MYCN           = anno_simple(dat_plot$MYCN,      col = mycn_col_fun,     border = TRUE),
  annotation_name_side = "left"
)

# ── 3.9 Risk score heatmap ────────────────────────────────────────────────────
rs_matrix <- matrix(dat_plot$risk_lp, nrow = 1)
colnames(rs_matrix) <- dat_plot$SampleID

ht_rnaseq <- Heatmap(
  rs_matrix,
  name  = "Risk score",
  col   = colorRamp2(
    c(min(dat_plot$risk_lp, na.rm = TRUE), 0, max(dat_plot$risk_lp, na.rm = TRUE)),
    c("#2166ac", "white", "#d6604d")
  ),
  border              = TRUE,
  show_column_names   = FALSE,
  cluster_columns     = FALSE,
  cluster_rows        = FALSE,
  row_names_side      = "left",
  row_names_gp        = gpar(fontsize = 10),
  height              = unit(0.6, "cm"),
  show_heatmap_legend = TRUE,
  top_annotation      = c(cat_anno, OS_dot_anno)
)

# ── 3.10 Manual legends ───────────────────────────────────────────────────────
lgd_tertile  <- Legend(labels = c("Low", "Intermediate", "High", "Unknown"),
                       title  = "Risk tertile",
                       legend_gp = gpar(fill = c("#2166ac","#f4a582","#d6604d","lightgrey")),
                       border = TRUE)
lgd_subgroup <- Legend(labels = c("Group 3", "Group 4", "Unknown"),
                       title  = "Subgroup",
                       legend_gp = gpar(fill = c("#FFD700","#008000","lightgrey")),
                       border = TRUE)
lgd_gender   <- Legend(labels = c("Male", "Female", "Unknown"),
                       title  = "Gender",
                       legend_gp = gpar(fill = c("white","black","lightgrey")),
                       border = TRUE)
lgd_age      <- Legend(labels = c("0-3", "3-16"),
                       title  = "Age group",
                       legend_gp = gpar(fill = c("#ffffcc","#41b6c4")),
                       border = TRUE)
lgd_mstage   <- Legend(labels = c("M0", "M+", "Unknown"),
                       title  = "M-stage",
                       legend_gp = gpar(fill = c("white","black","lightgrey")),
                       border = TRUE)
lgd_myc      <- Legend(labels = c("Not amplified", "Amplified"),
                       title  = "MYC",
                       legend_gp = gpar(fill = c("white","#e31a1c")),
                       border = TRUE)
lgd_mycn     <- Legend(labels = c("Not amplified", "Amplified"),
                       title  = "MYCN",
                       legend_gp = gpar(fill = c("white","#ff7f00")),
                       border = TRUE)

# ── 3.11 Ridge plots ──────────────────────────────────────────────────────────
fmt_p <- function(p) {
  if (is.na(p))    return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  if (p < 0.001)   return(paste0("p = ", formatC(p, format = "e", digits = 2)))
  paste0("p = ", round(p, 3))
}

rs_range    <- range(dat_plot$risk_lp, na.rm = TRUE)
ridge_theme <- theme_bw() + theme(legend.position = "none", panel.grid = element_blank())

# A — Subgroup
toPlot_A <- data.frame(risk_score = dat_plot$risk_lp, group = dat_plot$Subgroup) %>%
  filter(!is.na(group), group %in% c("Group 3", "Group 4"))
p_A <- fmt_p(t.test(risk_score ~ group, data = toPlot_A)$p.value)
A <- ggplot(toPlot_A, aes(risk_score, group)) +
  geom_density_ridges(aes(fill = group), scale = 1.5, rel_min_height = 0.001, alpha = 0.45) +
  scale_fill_manual(values = c("Group 3" = "#FFD700", "Group 4" = "#008000")) +
  xlim(rs_range) + ridge_theme + labs(x = "Risk score", y = "Subgroup") +
  annotate("text", x = rs_range[2] - 0.1, y = Inf, label = p_A, hjust = 1, vjust = 1.5, size = 3.5)

# B — M-stage
toPlot_B <- data.frame(risk_score = dat_plot$risk_lp, group = dat_plot$Mstage) %>%
  filter(!is.na(group), group != "Unknown")
p_B <- fmt_p(t.test(risk_score ~ group, data = toPlot_B)$p.value)
B <- ggplot(toPlot_B, aes(risk_score, group)) +
  geom_density_ridges(aes(fill = group), scale = 1.5, rel_min_height = 0.001, alpha = 0.45) +
  scale_fill_manual(values = c("0" = "grey80", "1" = "black")) +
  scale_y_discrete(labels = c("0" = "M0", "1" = "M+")) +
  xlim(rs_range) + ridge_theme + labs(x = "Risk score", y = "M-stage") +
  annotate("text", x = rs_range[2] - 0.1, y = Inf, label = p_B, hjust = 1, vjust = 1.5, size = 3.5)

# C — Gender
toPlot_C <- data.frame(risk_score = dat_plot$risk_lp, group = dat_plot$Gender) %>%
  filter(!is.na(group), group != "Unknown")
p_C <- fmt_p(t.test(risk_score ~ group, data = toPlot_C)$p.value)
C <- ggplot(toPlot_C, aes(risk_score, group)) +
  geom_density_ridges(aes(fill = group), scale = 1.5, rel_min_height = 0.001, alpha = 0.45) +
  scale_fill_manual(values = c("F" = "#d73027", "M" = "#4575b4")) +
  scale_y_discrete(labels = c("F" = "Female", "M" = "Male")) +
  xlim(rs_range) + ridge_theme + labs(x = "Risk score", y = "Gender") +
  annotate("text", x = rs_range[2] - 0.1, y = Inf, label = p_C, hjust = 1, vjust = 1.5, size = 3.5)

# D — MYC
toPlot_D <- data.frame(risk_score = dat_plot$risk_lp, group = as.character(dat_plot$MYC)) %>%
  filter(!is.na(group))
p_D <- fmt_p(t.test(risk_score ~ group, data = toPlot_D)$p.value)
D <- ggplot(toPlot_D, aes(risk_score, group)) +
  geom_density_ridges(aes(fill = group), scale = 1.5, rel_min_height = 0.001, alpha = 0.45) +
  scale_fill_manual(values = c("0" = "grey80", "1" = "#e31a1c")) +
  scale_y_discrete(labels = c("0" = "Not amplified", "1" = "Amplified")) +
  xlim(rs_range) + ridge_theme + labs(x = "Risk score", y = "MYC") +
  annotate("text", x = rs_range[2] - 0.1, y = Inf, label = p_D, hjust = 1, vjust = 1.5, size = 3.5)

# E — MYCN
toPlot_E <- data.frame(risk_score = dat_plot$risk_lp, group = as.character(dat_plot$MYCN)) %>%
  filter(!is.na(group))
p_E <- fmt_p(t.test(risk_score ~ group, data = toPlot_E)$p.value)
E_plot <- ggplot(toPlot_E, aes(risk_score, group)) +
  geom_density_ridges(aes(fill = group), scale = 1.5, rel_min_height = 0.001, alpha = 0.45) +
  scale_fill_manual(values = c("0" = "grey80", "1" = "#ff7f00")) +
  scale_y_discrete(labels = c("0" = "Not amplified", "1" = "Amplified")) +
  xlim(rs_range) + ridge_theme + labs(x = "Risk score", y = "MYCN") +
  annotate("text", x = rs_range[2] - 0.1, y = Inf, label = p_E, hjust = 1, vjust = 1.5, size = 3.5)

# F — Age group
toPlot_F <- data.frame(risk_score = dat_plot$risk_lp, group = dat_plot$AgeGroup) %>%
  filter(!is.na(group), group %in% c("0-3", "3-16")) %>%
  mutate(group = factor(group, levels = c("0-3", "3-16")))
p_F <- fmt_p(t.test(risk_score ~ group, data = toPlot_F)$p.value)
F_plot <- ggplot(toPlot_F, aes(risk_score, group)) +
  geom_density_ridges(aes(fill = group), scale = 1.5, rel_min_height = 0.001, alpha = 0.45) +
  scale_fill_manual(values = c("0-3" = "#ffffcc", "3-16" = "#41b6c4")) +
  xlim(rs_range) + ridge_theme + labs(x = "Risk score", y = "Age group") +
  annotate("text", x = rs_range[2] - 0.1, y = Inf, label = p_F, hjust = 1, vjust = 1.5, size = 3.5)

# ── 3.12 Export PDF ───────────────────────────────────────────────────────────
pdf("RNAseq_G3G4_validation_oncoprint.pdf", width = 14, height = 8)

draw(
  ht_rnaseq,
  heatmap_legend_side    = "right",
  annotation_legend_side = "right",
  annotation_legend_list = list(lgd_tertile, lgd_subgroup, lgd_gender,
                                lgd_age, lgd_mstage, lgd_myc, lgd_mycn),
  padding = unit(c(2, 30, 2, 2), "mm")
)

decorate_annotation("OS", {
  grid.lines(x = unit(c(0, 1), "npc"), y = unit(c(5, 5), "native"),
             gp = gpar(col = "black", lty = "dotted", lwd = 1))
})

grid.newpage()

print(ggarrange(A, B, C, D, E_plot, F_plot,
                ncol = 2, nrow = 3,
                labels = c("A", "B", "C", "D", "E", "F")))

dev.off()
cat("Done \u2014 RNAseq_G3G4_validation_oncoprint.pdf written.\n")


###########################################################################
# SECTION 4: SHH VALIDATION — MICROARRAY SIGNATURE PROJECTED ONTO RNA-SEQ
# NOTE: No OS/PFS survival data in E-MTAB-10767 SDRF for SHH samples.
# KM validation not possible. Clinical annotation figure only.
# NOTE: MYC excluded — no amplified cases in SHH cohort (n = 0/66).
###########################################################################

# ── 4.1 Build SHH RNA-seq expression matrix ───────────────────────────────────
shh_samples <- sdrf %>%
  filter(`Characteristics[subgroup]` == "SHH") %>%
  pull(`Source Name`) %>%
  as.character()

shh_in_mat <- intersect(shh_samples, colnames(expr_mat))
x_glm_shh  <- t(expr_mat[, shh_in_mat, drop = FALSE])

cat("x_glm_shh dimensions:", nrow(x_glm_shh), "samples x", ncol(x_glm_shh), "genes\n")

# ── 4.2 Gene overlap check ────────────────────────────────────────────────────
genes_in_model_shh <- names(final_coefs_shh)

overlap_shh <- intersect(genes_in_model_shh, colnames(x_glm_shh))
missing_shh <- setdiff(genes_in_model_shh, colnames(x_glm_shh))

cat("Genes in SHH model:", length(genes_in_model_shh), "\n")
cat("Found in RNA-seq:", length(overlap_shh), "\n")
cat("Missing:", length(missing_shh), "\n")

# ── 4.3 Compute RNA-seq SHH risk scores (manual dot product) ─────────────────
coefs_matched_shh <- final_coefs_shh[overlap_shh]
x_glm_shh_matched <- x_glm_shh[, overlap_shh, drop = FALSE]

risk_score_shh.rnaseq <- as.numeric(x_glm_shh_matched %*% coefs_matched_shh)
names(risk_score_shh.rnaseq) <- rownames(x_glm_shh)

# ── 4.4 Build risk_out_shh with tertiles ──────────────────────────────────────
tertile_cuts_shh <- quantile(risk_score_shh.rnaseq, probs = c(1/3, 2/3))

risk_out_shh <- data.frame(
  SampleID = names(risk_score_shh.rnaseq),
  risk_lp  = risk_score_shh.rnaseq
) %>%
  mutate(
    tertile = cut(risk_lp,
                  breaks = c(-Inf, tertile_cuts_shh, Inf),
                  labels = c("Low", "Intermediate", "High"))
  )

table(risk_out_shh$tertile)


###########################################################################
# SECTION 5: SHH CLINICAL ANNOTATION FIGURE
# Oncoprint (Page 1) + Ridge plots (Page 2)
###########################################################################

names(sdrf) <- make.unique(names(sdrf))

# ── 5.1 Pull clinical metadata from SDRF ─────────────────────────────────────
clinical_shh <- sdrf %>%
  transmute(
    SampleID   = as.character(`Source Name`),
    Subgroup   = as.character(`Characteristics[subgroup]`),
    Sex_raw    = as.character(`Characteristics[sex]`),
    Age_raw    = as.character(`Characteristics[age]`),
    Mstage_raw = as.character(`Characteristics[m+]`),
    MYCN_raw   = as.character(`Characteristics[mycn amplification]`)
  ) %>%
  filter(Subgroup == "SHH") %>%
  mutate(
    Gender = ifelse(tolower(trimws(Sex_raw)) %in% c("m","male"),   "M",
                    ifelse(tolower(trimws(Sex_raw)) %in% c("f","female"), "F", "Unknown")),
    Age = local({
      x <- tolower(trimws(as.character(Age_raw)))
      # NOTE: SHH SDRF uses different age string formats to G3/4 -- requires
      # explicit pattern matching for "not available", "over 16", "0 to 3" etc.
      ifelse(grepl("not available|na|unknown", x),                  "Unknown",
             ifelse(grepl("over 16|>16|> 16|adult",  x),                  ">16",
                    ifelse(grepl("0 to 3|0-3|infant|toddler", x),                "0-3",
                           ifelse(grepl("3 to 16|3-16|child|pediatric|paediatric", x),  "3-16",
                                  "Unknown"))))
    }),
    Mstage = ifelse(tolower(trimws(Mstage_raw)) %in% c("true","1","yes","y"), "1",
                    ifelse(tolower(trimws(Mstage_raw)) %in% c("false","0","no","n"), "0", "Unknown")),
    MYCN   = ifelse(tolower(trimws(MYCN_raw)) %in% c("true","1","yes","y"),   "1", "0")
  )

# ── 5.2 Merge and order by ascending risk score ───────────────────────────────
dat_plot_shh <- risk_out_shh %>%
  left_join(clinical_shh, by = "SampleID") %>%
  arrange(risk_lp)

cat("Samples:", nrow(dat_plot_shh), "\n")

# ── 5.3 Clean NAs ─────────────────────────────────────────────────────────────
dat_plot_shh$Gender   <- ifelse(is.na(dat_plot_shh$Gender)   | dat_plot_shh$Gender   == "", "Unknown", dat_plot_shh$Gender)
dat_plot_shh$Subgroup <- ifelse(is.na(dat_plot_shh$Subgroup) | dat_plot_shh$Subgroup == "", "Unknown", dat_plot_shh$Subgroup)
dat_plot_shh$Mstage   <- ifelse(is.na(dat_plot_shh$Mstage)   | dat_plot_shh$Mstage   == "", "Unknown", dat_plot_shh$Mstage)
dat_plot_shh$MYCN     <- as.character(dat_plot_shh$MYCN)
dat_plot_shh$tertile  <- ifelse(is.na(dat_plot_shh$tertile)  | dat_plot_shh$tertile  == "", "Unknown",
                                as.character(dat_plot_shh$tertile))

# ── 5.4 Colour palettes ───────────────────────────────────────────────────────
tertile_col_fun <- c("Low" = "#2166ac", "Intermediate" = "#f4a582", "High" = "#d6604d", "Unknown" = "lightgrey")
gender_col_fun  <- c("M" = "white",     "F" = "black",              "Unknown" = "lightgrey")
mstage_col_fun  <- c("0" = "white",     "1" = "black",              "Unknown" = "lightgrey")
mycn_col_fun    <- c("0" = "white",     "1" = "#ff7f00")
age_col_fun     <- c("0-3" = "#ffffcc", "3-16" = "#41b6c4")

# ── 5.5 Age group — explicit parsing for SHH SDRF format ─────────────────────
dat_plot_shh$AgeGroup <- local({
  x <- tolower(trimws(as.character(dat_plot_shh$Age)))
  ifelse(grepl("not available|na|unknown", x),                  "Unknown",
         ifelse(grepl("over 16|>16|> 16|adult",  x),                  ">16",
                ifelse(grepl("0 to 3|0-3|infant|toddler", x),                "0-3",
                       ifelse(grepl("3 to 16|3-16|child|pediatric|paediatric", x),  "3-16",
                              "Unknown"))))
})

dat_plot_shh <- dat_plot_shh %>% filter(!AgeGroup %in% c(">16", "Unknown"))
cat("Samples after excluding >16 and unknown age:", nrow(dat_plot_shh), "\n")
table(dat_plot_shh$AgeGroup)

# ── 5.6 Categorical bar annotations ──────────────────────────────────────────
cat_anno_shh <- HeatmapAnnotation(
  `Risk tertile` = anno_simple(dat_plot_shh$tertile,   col = tertile_col_fun,  border = TRUE),
  Gender         = anno_simple(dat_plot_shh$Gender,    col = gender_col_fun,   border = TRUE),
  `Age group`    = anno_simple(dat_plot_shh$AgeGroup,  col = age_col_fun,      border = TRUE),
  `M-stage`      = anno_simple(dat_plot_shh$Mstage,    col = mstage_col_fun,   border = TRUE),
  MYCN           = anno_simple(dat_plot_shh$MYCN,      col = mycn_col_fun,     border = TRUE),
  annotation_name_side = "left"
)

# ── 5.7 Risk score heatmap ────────────────────────────────────────────────────
rs_matrix_shh <- matrix(dat_plot_shh$risk_lp, nrow = 1)
colnames(rs_matrix_shh) <- dat_plot_shh$SampleID

ht_shh <- Heatmap(
  rs_matrix_shh,
  name  = "Risk score",
  col   = colorRamp2(
    c(min(dat_plot_shh$risk_lp, na.rm = TRUE),
      median(dat_plot_shh$risk_lp, na.rm = TRUE),
      max(dat_plot_shh$risk_lp, na.rm = TRUE)),
    c("#2166ac", "white", "#d6604d")
  ),
  border              = TRUE,
  show_column_names   = FALSE,
  cluster_columns     = FALSE,
  cluster_rows        = FALSE,
  row_names_side      = "left",
  row_names_gp        = gpar(fontsize = 10),
  height              = unit(0.6, "cm"),
  show_heatmap_legend = TRUE,
  top_annotation      = cat_anno_shh
)

# ── 5.8 Manual legends ────────────────────────────────────────────────────────
lgd_tertile <- Legend(labels = c("Low", "Intermediate", "High", "Unknown"),
                      title  = "Risk tertile",
                      legend_gp = gpar(fill = c("#2166ac","#f4a582","#d6604d","lightgrey")),
                      border = TRUE)
lgd_gender  <- Legend(labels = c("Male", "Female", "Unknown"),
                      title  = "Gender",
                      legend_gp = gpar(fill = c("white","black","lightgrey")),
                      border = TRUE)
lgd_age     <- Legend(labels = c("0-3", "3-16"),
                      title  = "Age group",
                      legend_gp = gpar(fill = c("#ffffcc","#41b6c4")),
                      border = TRUE)
lgd_mstage  <- Legend(labels = c("M0", "M+", "Unknown"),
                      title  = "M-stage",
                      legend_gp = gpar(fill = c("white","black","lightgrey")),
                      border = TRUE)
lgd_mycn    <- Legend(labels = c("Not amplified", "Amplified"),
                      title  = "MYCN",
                      legend_gp = gpar(fill = c("white","#ff7f00")),
                      border = TRUE)

# ── 5.9 Ridge plots ───────────────────────────────────────────────────────────
rs_range    <- range(dat_plot_shh$risk_lp, na.rm = TRUE)
ridge_theme <- theme_bw() + theme(legend.position = "none", panel.grid = element_blank())

# A — M-stage
toPlot_A <- data.frame(risk_score = dat_plot_shh$risk_lp, group = dat_plot_shh$Mstage) %>%
  filter(!is.na(group), group != "Unknown")
p_A <- fmt_p(t.test(risk_score ~ group, data = toPlot_A)$p.value)
A <- ggplot(toPlot_A, aes(risk_score, group)) +
  geom_density_ridges(aes(fill = group), scale = 1.5, rel_min_height = 0.001, alpha = 0.45) +
  scale_fill_manual(values = c("0" = "grey80", "1" = "black")) +
  scale_y_discrete(labels = c("0" = "M0", "1" = "M+")) +
  xlim(rs_range) + ridge_theme + labs(x = "Risk score", y = "M-stage") +
  annotate("text", x = rs_range[2] - 0.1, y = Inf, label = p_A, hjust = 1, vjust = 1.5, size = 3.5)

# B — Gender
toPlot_B <- data.frame(risk_score = dat_plot_shh$risk_lp, group = dat_plot_shh$Gender) %>%
  filter(!is.na(group), group != "Unknown")
p_B <- fmt_p(t.test(risk_score ~ group, data = toPlot_B)$p.value)
B <- ggplot(toPlot_B, aes(risk_score, group)) +
  geom_density_ridges(aes(fill = group), scale = 1.5, rel_min_height = 0.001, alpha = 0.45) +
  scale_fill_manual(values = c("F" = "#d73027", "M" = "#4575b4")) +
  scale_y_discrete(labels = c("F" = "Female", "M" = "Male")) +
  xlim(rs_range) + ridge_theme + labs(x = "Risk score", y = "Gender") +
  annotate("text", x = rs_range[2] - 0.1, y = Inf, label = p_B, hjust = 1, vjust = 1.5, size = 3.5)

# C — MYCN
toPlot_C <- data.frame(risk_score = dat_plot_shh$risk_lp,
                       group = as.character(dat_plot_shh$MYCN)) %>%
  filter(!is.na(group))
p_C <- fmt_p(t.test(risk_score ~ group, data = toPlot_C)$p.value)
C <- ggplot(toPlot_C, aes(risk_score, group)) +
  geom_density_ridges(aes(fill = group), scale = 1.5, rel_min_height = 0.001, alpha = 0.45) +
  scale_fill_manual(values = c("0" = "grey80", "1" = "#ff7f00")) +
  scale_y_discrete(labels = c("0" = "Not amplified", "1" = "Amplified")) +
  xlim(rs_range) + ridge_theme + labs(x = "Risk score", y = "MYCN") +
  annotate("text", x = rs_range[2] - 0.1, y = Inf, label = p_C, hjust = 1, vjust = 1.5, size = 3.5)

# D — Age group
toPlot_D <- data.frame(risk_score = dat_plot_shh$risk_lp, group = dat_plot_shh$AgeGroup) %>%
  filter(!is.na(group), group %in% c("0-3", "3-16")) %>%
  mutate(group = factor(group, levels = c("0-3", "3-16")))
p_D <- fmt_p(t.test(risk_score ~ group, data = toPlot_D)$p.value)
D <- ggplot(toPlot_D, aes(risk_score, group)) +
  geom_density_ridges(aes(fill = group), scale = 1.5, rel_min_height = 0.001, alpha = 0.45) +
  scale_fill_manual(values = c("0-3" = "#ffffcc", "3-16" = "#41b6c4")) +
  xlim(rs_range) + ridge_theme + labs(x = "Risk score", y = "Age group") +
  annotate("text", x = rs_range[2] - 0.1, y = Inf, label = p_D, hjust = 1, vjust = 1.5, size = 3.5)

# ── 5.10 Export PDF ───────────────────────────────────────────────────────────
pdf("RNAseq_SHH_validation_oncoprint.pdf", width = 14, height = 8)

draw(
  ht_shh,
  heatmap_legend_side    = "right",
  annotation_legend_side = "right",
  annotation_legend_list = list(lgd_tertile, lgd_gender,
                                lgd_age, lgd_mstage, lgd_mycn),
  padding = unit(c(2, 30, 2, 2), "mm")
)

grid.newpage()

print(ggarrange(A, B, C, D,
                ncol = 2, nrow = 2,
                labels = c("A", "B", "C", "D")))

dev.off()
cat("Done \u2014 RNAseq_SHH_validation_oncoprint.pdf written.\n")




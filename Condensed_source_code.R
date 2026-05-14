################################################################################
#
# MEDULLOBLASTOMA SURVIVAL ANALYSIS
# Dataset: GSE85217 (Taylor Lab, 763 pediatric MB samples)
#
# Student: Ewan McGibbon
# Date:   7/4/2026
#
################################################################################

# =============================================================================
# 0. SETUP
# =============================================================================

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
library(ggrepel)
library(gridExtra)

########################################################################################################################
#
#   SECTION 1: BUILD THE ANALYSIS MATRIX
#   Input:  GSE85217 expression file
#   Output: expr_top (2000 most variable autosomal genes x 763 samples)
#
########################################################################################################################

# -----------------------------------------------------------------------------
# 1.1 Load microarray expression data
# -----------------------------------------------------------------------------

expr_micro <- fread("GSE85217_M_exp_763_MB_SubtypeStudy_TaylorLab.txt")
expr_df_micro <- as.data.frame(expr_micro)

expr_mat_micro <- as.matrix(expr_df_micro[, 6:ncol(expr_df_micro)])
rownames(expr_mat_micro) <- expr_df_micro$EnsemblGeneID_from_ensemblv77

dim(expr_mat_micro)

# -----------------------------------------------------------------------------
# 1.2 Remove sex-chromosome genes (X/Y)
# -----------------------------------------------------------------------------

mirrors <- c("useast", "uswest", "asia")
ensembl <- NULL
for (mirror in mirrors) {
  ensembl <- tryCatch(
    useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = mirror),
    error = function(e) { message("Mirror '", mirror, "' failed"); NULL }
  )
  if (!is.null(ensembl)) { message("Connected via: ", mirror); break }
}
if (is.null(ensembl)) stop("All Ensembl mirrors failed.")

annot_chr <- getBM(
  attributes = c("ensembl_gene_id", "chromosome_name"),
  filters    = "ensembl_gene_id",
  values     = rownames(expr_mat_micro),
  mart       = ensembl
)

ens_keep      <- unique(annot_chr$ensembl_gene_id[annot_chr$chromosome_name %in% as.character(1:22)])
expr_noXY     <- expr_mat_micro[rownames(expr_mat_micro) %in% ens_keep, ]
message(nrow(expr_noXY), " autosomal genes retained")

# -----------------------------------------------------------------------------
# 1.3 Select top 2000 most variable genes
# -----------------------------------------------------------------------------

gene_sd     <- apply(expr_noXY, 1, sd, na.rm = TRUE)
expr_top    <- expr_noXY[order(gene_sd, decreasing = TRUE)[1:2000], ]
dim(expr_top)

########################################################################################################################
#
#   SECTION 2: PCA + UMAP
#   Input:  expr_top, meta_sub
#   Output: Figure 1 — PCA + UMAP coloured by subgroup
#
########################################################################################################################

set.seed(123)

# -----------------------------------------------------------------------------
# 2.1 Load metadata
# -----------------------------------------------------------------------------

meta <- read.csv("Taylor_Pheno_forDan.csv", stringsAsFactors = FALSE)

overlaps        <- sapply(meta, function(col) length(intersect(col, colnames(expr_top))))
sample_col_name <- names(overlaps)[which.max(overlaps)]

common_samples <- intersect(colnames(expr_top), meta[[sample_col_name]])
expr_top_sub   <- expr_top[, common_samples]
meta_sub       <- meta[match(common_samples, meta[[sample_col_name]]), ]
stopifnot(all(colnames(expr_top_sub) == meta_sub[[sample_col_name]]))

# Consensus subgroup colours
subgroup_cols <- c("WNT" = "#0000FF", "SHH" = "#FF0000", "Group3" = "#FFFF00", "Group4" = "#008000")

# -----------------------------------------------------------------------------
# 2.2 PCA
# -----------------------------------------------------------------------------

pca_res2 <- prcomp(t(expr_top_sub), scale. = TRUE)
pca_df2  <- data.frame(
  Sample   = rownames(pca_res2$x),
  PC1      = pca_res2$x[, 1],
  PC2      = pca_res2$x[, 2],
  Subgroup = meta_sub$Subgroup
)

p_pca <- ggplot(pca_df2, aes(PC1, PC2, colour = Subgroup)) +
  geom_point(size = 1.8, alpha = 0.85) +
  scale_colour_manual(values = subgroup_cols) +
  labs(title = "PCA\n(top 2000 autosomal genes)", x = "PC1", y = "PC2", tag = "A") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

# -----------------------------------------------------------------------------
# 2.3 UMAP
# -----------------------------------------------------------------------------

umap_res2 <- umap(t(expr_top_sub), n_neighbors = 15, min_dist = 0.3, metric = "euclidean")
umap_df2  <- data.frame(
  Sample   = colnames(expr_top_sub),
  UMAP1    = umap_res2[, 1],
  UMAP2    = umap_res2[, 2],
  Subgroup = meta_sub$Subgroup
)

p_umap <- ggplot(umap_df2, aes(UMAP1, UMAP2, colour = Subgroup)) +
  geom_point(size = 1.8, alpha = 0.85) +
  scale_colour_manual(values = subgroup_cols) +
  labs(title = "UMAP\n(top 2000 autosomal genes)", x = "UMAP1", y = "UMAP2", tag = "B") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

# Figure 1
fig1 <- (p_pca | p_umap) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")
print(fig1)

########################################################################################################################
#
#   SECTION 3: KAPLAN-MEIER SURVIVAL ANALYSIS
#   Input:  meta_sub
#   Output: Figure 2 — four-panel KM (subgroup, metastasis, histology, age)
#
########################################################################################################################

library(gridExtra)

# -----------------------------------------------------------------------------
# 3.1 Prepare clinical variables
# -----------------------------------------------------------------------------

meta_sub$Metastasis_status <- factor(
  meta_sub$Met.status..1.Met..0.M0.,
  levels = c("0", "1"),
  labels = c("M0 / non-metastatic", "Metastatic")
)

raw_met <- as.character(meta_sub$Met.status..1.Met..0.M0.)
raw_met[raw_met %in% c("0","0.0","0 ","0.00")] <- "M0"
raw_met[raw_met %in% c("1","1.0","1 ","1.00")] <- "Met"
meta_sub$Met_status <- relevel(factor(raw_met), ref = "M0")

meta_sub$histology <- relevel(factor(meta_sub$histology), ref = "MBEN")
meta_sub$Subgroup  <- factor(meta_sub$Subgroup, levels = c("WNT","SHH","Group3","Group4"))

meta_sub$AgeGroup <- case_when(
  as.numeric(meta_sub$Age) < 3                                   ~ "0-3",
  as.numeric(meta_sub$Age) >= 3 & as.numeric(meta_sub$Age) < 10 ~ "4-10",
  as.numeric(meta_sub$Age) >= 10                                 ~ "10+",
  TRUE                                                            ~ NA_character_
)
meta_sub$AgeGroup <- relevel(factor(meta_sub$AgeGroup, levels = c("0-3","4-10","10+")), ref = "4-10")

cat("Age group distribution:\n")
print(table(meta_sub$AgeGroup, useNA = "always"))

# -----------------------------------------------------------------------------
# 3.2 Fit KM models
# -----------------------------------------------------------------------------

fit_subgroup <- survfit(Surv(OS..years., Dead) ~ Subgroup,          data = meta_sub)
fit_met      <- survfit(Surv(OS..years., Dead) ~ Metastasis_status, data = meta_sub)
fit_hist     <- survfit(Surv(OS..years., Dead) ~ histology,         data = meta_sub)
fit_age      <- survfit(Surv(OS..years., Dead) ~ AgeGroup,          data = meta_sub)

# -----------------------------------------------------------------------------
# 3.3 Generate KM plots
# -----------------------------------------------------------------------------

km_subgroup <- ggsurvplot(fit_subgroup, data = meta_sub,
                          risk.table = FALSE, pval = TRUE, conf.int = FALSE,
                          legend.title = "", legend.labs = c("WNT","SHH","Group3","Group4"),
                          legend = "top", palette = subgroup_cols,
                          xlab = "Time (years)", ylab = "Survival probability",
                          title = "Subgroup", xlim = c(0,10), break.time.by = 2,
                          ggtheme = theme_classic(base_size = 11))

km_met <- ggsurvplot(fit_met, data = meta_sub,
                     risk.table = FALSE, pval = TRUE, conf.int = FALSE,
                     legend.title = "", legend.labs = c("Non-metastatic","Metastatic"),
                     legend = "top",
                     xlab = "Time (years)", ylab = "Survival probability",
                     title = "Metastasis", xlim = c(0,10), break.time.by = 2,
                     ggtheme = theme_classic(base_size = 11))

# legend.labs match factor levels: MBEN, Classic, Desmoplastic, LCA
km_hist <- ggsurvplot(fit_hist, data = meta_sub,
                      risk.table = FALSE, pval = TRUE, conf.int = FALSE,
                      legend.title = "", legend.labs = c("MBEN","Classic","Desmoplastic","LCA"),
                      legend = "top", palette = c("#c51b7d","#878787","#762a83","#1b7837"),
                      xlab = "Time (years)", ylab = "Survival probability",
                      title = "Histology", xlim = c(0,10), break.time.by = 2,
                      ggtheme = theme_classic(base_size = 11))

km_age <- ggsurvplot(fit_age, data = meta_sub,
                     risk.table = FALSE, pval = TRUE, conf.int = FALSE,
                     legend.title = "", legend.labs = c("0-3 years","4-10 years","10+ years"),
                     legend = "top",
                     xlab = "Time (years)", ylab = "Survival probability",
                     title = "Age group", xlim = c(0,10), break.time.by = 2,
                     ggtheme = theme_classic(base_size = 11))

# -----------------------------------------------------------------------------
# 3.4 Figure 2 — four-panel layout
# -----------------------------------------------------------------------------

fig2 <- arrangeGrob(
  km_subgroup$plot, km_met$plot,
  km_hist$plot,     km_age$plot,
  ncol = 2, nrow = 2
)
grid::grid.draw(fig2)

########################################################################################################################
#
#   SECTION 4: COX REGRESSION — CLINICAL PREDICTORS
#   Input:  meta_sub, meta_g34, meta_shh
#   Output: Univariable Cox models
#
########################################################################################################################

# -----------------------------------------------------------------------------
# 4.1 Full cohort Cox models
# -----------------------------------------------------------------------------

cox_hist     <- coxph(Surv(OS..years., Dead) ~ histology,  data = meta_sub)
cox_subgroup <- coxph(Surv(OS..years., Dead) ~ Subgroup,   data = meta_sub)
cox_met      <- coxph(Surv(OS..years., Dead) ~ Met_status, data = meta_sub)
cox_age      <- coxph(Surv(OS..years., Dead) ~ AgeGroup,   data = meta_sub)

summary(cox_hist); summary(cox_subgroup); summary(cox_met); summary(cox_age)

# -----------------------------------------------------------------------------
# 4.2 Split into subgroups
# -----------------------------------------------------------------------------

meta_g34 <- meta_sub[meta_sub$Subgroup %in% c("Group3","Group4"), ]
meta_shh <- meta_sub[meta_sub$Subgroup == "SHH", ]

meta_g34$age_group  <- factor(meta_g34$AgeGroup, levels = c("4-10","0-3"))
meta_shh$age_group  <- factor(meta_shh$AgeGroup, levels = c("4-10","0-3"))

meta_g34$AgeGroup   <- relevel(factor(as.character(meta_g34$AgeGroup)),   ref = "4-10")
meta_g34$Met_status <- relevel(factor(as.character(meta_g34$Met_status)), ref = "M0")
meta_g34$histology  <- relevel(factor(as.character(meta_g34$histology)),  ref = "MBEN")

meta_shh$AgeGroup   <- relevel(factor(as.character(meta_shh$AgeGroup)),   ref = "4-10")
meta_shh$Met_status <- relevel(factor(as.character(meta_shh$Met_status)), ref = "M0")
meta_shh$histology  <- relevel(factor(as.character(meta_shh$histology)),  ref = "MBEN")

cox_age_g34  <- coxph(Surv(OS..years., Dead) ~ AgeGroup,   data = meta_g34)
cox_hist_g34 <- coxph(Surv(OS..years., Dead) ~ histology,  data = meta_g34)
cox_met_g34  <- coxph(Surv(OS..years., Dead) ~ Met_status, data = meta_g34)

cox_age_shh  <- coxph(Surv(OS..years., Dead) ~ AgeGroup,   data = meta_shh)
cox_hist_shh <- coxph(Surv(OS..years., Dead) ~ histology,  data = meta_shh)
cox_met_shh  <- coxph(Surv(OS..years., Dead) ~ Met_status, data = meta_shh)

########################################################################################################################
#
#   SECTION 5: SUBGROUP-SPECIFIC SURVIVAL ANALYSIS
#   Input:  meta_g34, meta_shh, myc/mycn CSVs
#   Output: KM panels for metastasis, age, MYC/MYCN amplification
#
########################################################################################################################

# -----------------------------------------------------------------------------
# 5.1 MYC/MYCN amplification status
# -----------------------------------------------------------------------------

myc_ids  <- read.csv("mycTayDems.csv",  stringsAsFactors = FALSE)
mycn_ids <- read.csv("mycnTayDems.csv", stringsAsFactors = FALSE)

vec_myc_ids  <- myc_ids$Study_ID
vec_mycn_ids <- mycn_ids$Study_ID
myc_amp_ids  <- vec_myc_ids
mycn_amp_ids <- vec_mycn_ids

meta_sub$MYC_group <- "None"
meta_sub$MYC_group[meta_sub$Study_ID %in% vec_myc_ids]  <- "MYC-amplified"
meta_sub$MYC_group[meta_sub$Study_ID %in% vec_mycn_ids] <- "MYCN-amplified"
meta_sub$MYC_group <- factor(meta_sub$MYC_group, levels = c("None","MYC-amplified","MYCN-amplified"))

meta_g34$MYC  <- factor(ifelse(meta_g34$Study_ID %in% vec_myc_ids,  "Amplified","Not amplified"), levels = c("Not amplified","Amplified"))
meta_g34$MYCN <- factor(ifelse(meta_g34$Study_ID %in% vec_mycn_ids, "Amplified","Not amplified"), levels = c("Not amplified","Amplified"))
meta_shh$MYCN <- factor(ifelse(meta_shh$Study_ID %in% vec_mycn_ids, "Amplified","Not amplified"), levels = c("Not amplified","Amplified"))

# -----------------------------------------------------------------------------
# 5.2 Full cohort: MYC/MYCN KM
# -----------------------------------------------------------------------------

fit_myc_all <- survfit(Surv(OS..years., Dead) ~ MYC_group, data = meta_sub)
km_myc_all  <- ggsurvplot(fit_myc_all, data = meta_sub,
                          title = "Overall survival by MYC/MYCN amplification status",
                          legend.title = "Amplification", legend.labs = c("None","MYC-amplified","MYCN-amplified"),
                          risk.table = TRUE, risk.table.height = 0.28, xlim = c(0,10),
                          xlab = "Overall survival (years)", ylab = "Survival probability",
                          pval = TRUE, pval.coord = c(0.3, 0.15),
                          palette = c("grey60","#d73027","#fd8d3c"),
                          ggtheme = theme_classic(base_size = 11), tables.theme = theme_cleantable())
print(km_myc_all)

# -----------------------------------------------------------------------------
# 5.3 G3/4: MYC and MYCN side-by-side
# -----------------------------------------------------------------------------

fit_myc_g34  <- survfit(Surv(OS..years., Dead) ~ MYC,  data = meta_g34)
fit_mycn_g34 <- survfit(Surv(OS..years., Dead) ~ MYCN, data = meta_g34)

km_myc_g34 <- ggsurvplot(fit_myc_g34, data = meta_g34,
                         title = "Group 3/4: MYC amplification", legend.title = "MYC",
                         legend.labs = c("Not amplified","Amplified"),
                         risk.table = TRUE, risk.table.height = 0.25, xlim = c(0,10),
                         xlab = "Overall survival (years)", ylab = "Survival probability",
                         pval = TRUE, pval.coord = c(0.3, 0.15), palette = c("grey60","#d73027"),
                         ggtheme = theme_classic(base_size = 11), tables.theme = theme_cleantable())

km_mycn_g34 <- ggsurvplot(fit_mycn_g34, data = meta_g34,
                          title = "Group 3/4: MYCN amplification", legend.title = "MYCN",
                          legend.labs = c("Not amplified","Amplified"),
                          risk.table = TRUE, risk.table.height = 0.25, xlim = c(0,10),
                          xlab = "Overall survival (years)", ylab = "Survival probability",
                          pval = TRUE, pval.coord = c(0.3, 0.15), palette = c("grey60","#fd8d3c"),
                          ggtheme = theme_classic(base_size = 11), tables.theme = theme_cleantable())

arrange_ggsurvplots(list(km_myc_g34, km_mycn_g34), ncol = 2, nrow = 1)

# -----------------------------------------------------------------------------
# 5.4 SHH: MYCN only
# -----------------------------------------------------------------------------

fit_mycn_shh <- survfit(Surv(OS..years., Dead) ~ MYCN, data = meta_shh)
km_mycn_shh  <- ggsurvplot(fit_mycn_shh, data = meta_shh,
                           title = "SHH: MYCN amplification", legend.title = "MYCN",
                           legend.labs = c("Not amplified","Amplified"),
                           risk.table = TRUE, risk.table.height = 0.25, xlim = c(0,10),
                           xlab = "Overall survival (years)", ylab = "Survival probability",
                           pval = TRUE, pval.coord = c(0.3, 0.15), palette = c("grey60","#fd8d3c"),
                           ggtheme = theme_classic(base_size = 11), tables.theme = theme_cleantable())
print(km_mycn_shh)

# -----------------------------------------------------------------------------
# 5.5 KM panels: metastasis and age group (appendix)
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

fit_age_g34 <- survfit(Surv(OS..years., Dead) ~ age_group, data = meta_g34)
km_age_g34  <- ggsurvplot(fit_age_g34, data = meta_g34,
                          title = "Group 3/4: age group", legend.title = "Age group (years)",
                          legend.labs = c("4-10","0-3"), risk.table = TRUE, risk.table.height = 0.25,
                          xlim = c(0,10), xlab = "Overall survival (years)", ylab = "Survival probability",
                          pval = TRUE, pval.coord = c(0.3,0.15), palette = c("#F08080","#20B2AA"),
                          ggtheme = theme_classic(base_size = 11), tables.theme = theme_cleantable())

fit_age_shh <- survfit(Surv(OS..years., Dead) ~ age_group, data = meta_shh)
km_age_shh  <- ggsurvplot(fit_age_shh, data = meta_shh,
                          title = "SHH: age group", legend.title = "Age group (years)",
                          legend.labs = c("4-10","0-3"), risk.table = TRUE, risk.table.height = 0.25,
                          xlim = c(0,10), xlab = "Overall survival (years)", ylab = "Survival probability",
                          pval = TRUE, pval.coord = c(0.3,0.15), palette = c("#F08080","#20B2AA"),
                          ggtheme = theme_classic(base_size = 11), tables.theme = theme_cleantable())

# Appendix figure
arrange_ggsurvplots(list(km_meta_g34, km_meta_shh, km_age_g34, km_age_shh), ncol = 2, nrow = 2)

# -----------------------------------------------------------------------------
# 5.6 KM panels: metastasis 
# -----------------------------------------------------------------------------

library(cowplot)

# Extract just the plot component and wipe the title cleanly
get_plot <- function(km_obj, label) {
  km_obj$plot$labels$title <- NULL
  km_obj$plot + 
    theme(plot.title = element_text(size = 11, face = "bold", hjust = 0)) +
    ggtitle(label)
}

plot_A <- get_plot(km_myc_all,  "A  Overall cohort: MYC/MYCN amplification")
plot_B <- get_plot(km_myc_g34,  "B  Group 3/4: MYC amplification")
plot_C <- get_plot(km_mycn_g34, "C  Group 3/4: MYCN amplification")
plot_D <- get_plot(km_mycn_shh, "D  SHH: MYCN amplification")

fig_myc_mycn <- plot_grid(
  plot_A, plot_B,
  plot_C, plot_D,
  ncol = 2, nrow = 2,
  align = "hv"
)

title <- ggdraw() + 
  draw_label(
    "MYC and MYCN amplification and overall survival in medulloblastoma",
    fontface = "bold", size = 13, x = 0.5, hjust = 0.5
  )

plot_grid(title, fig_myc_mycn, ncol = 1, rel_heights = c(0.05, 1))

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
#
#   SECTION 6: GENOME-WIDE COX SCAN + VOLCANO/FOREST PLOTS
#   Input:  expr_micro (microarray), meta_g34, meta_shh
#   Output: results_g34, results_shh, volcano plots, forest plots
#
########################################################################################################################

# -----------------------------------------------------------------------------
# 6.1 Build expression matrix aligned to subgroups
# -----------------------------------------------------------------------------

gene_ids_micro   <- expr_df_micro$EnsemblGeneID_from_ensemblv77
hgnc_syms_micro  <- expr_df_micro$HGNC_symbol_from_ensemblv77
sample_cols_micro <- grep("^MB_SubtypeStudy", names(expr_df_micro), value = TRUE)

expr_mat_scan    <- as.matrix(expr_df_micro[, sample_cols_micro, drop = FALSE])
rownames(expr_mat_scan) <- gene_ids_micro

id_col <- "Study_ID"

# G3/4
ids_g34        <- meta_g34[[id_col]]
common_ids_g34 <- intersect(ids_g34, colnames(expr_mat_scan))
expr_g34_aln   <- expr_mat_scan[, common_ids_g34]
meta_g34_aln   <- meta_g34[match(common_ids_g34, ids_g34), ]
keep_g34       <- !is.na(meta_g34_aln$OS..years.) & !is.na(meta_g34_aln$Dead)
meta_g34_cc    <- meta_g34_aln[keep_g34, ]
expr_g34_cc    <- expr_g34_aln[, keep_g34]
time_g34       <- meta_g34_cc$OS..years.
status_g34     <- meta_g34_cc$Dead

# SHH
ids_shh        <- meta_shh[[id_col]]
common_ids_shh <- intersect(ids_shh, colnames(expr_mat_scan))
expr_shh_aln   <- expr_mat_scan[, common_ids_shh]
meta_shh_aln   <- meta_shh[match(common_ids_shh, ids_shh), ]
keep_shh       <- !is.na(meta_shh_aln$OS..years.) & !is.na(meta_shh_aln$Dead)
meta_shh_cc    <- meta_shh_aln[keep_shh, ]
expr_shh_cc    <- expr_shh_aln[, keep_shh]
time_shh       <- meta_shh_cc$OS..years.
status_shh     <- meta_shh_cc$Dead

# -----------------------------------------------------------------------------
# 6.2 Genome-wide Cox scan: G3/4
# -----------------------------------------------------------------------------

cox_g34_list <- apply(expr_g34_cc, 1, function(row_expr) {
  df  <- data.frame(time = time_g34, status = status_g34, expr = as.numeric(row_expr))
  fit <- coxph(Surv(time, status) ~ expr, data = df)
  s   <- summary(fit)
  c(coef = s$coef[1,"coef"], se = s$coef[1,"se(coef)"],
    HR = s$coef[1,"exp(coef)"], p = s$coef[1,"Pr(>|z|)"])
})

results_g34 <- as.data.frame(t(cox_g34_list))
rownames(results_g34) <- NULL
results_g34$ensembl_gene_id <- gene_ids_micro
results_g34$hgnc_symbol     <- hgnc_syms_micro
results_g34[, c("coef","se","HR","p")] <- lapply(results_g34[, c("coef","se","HR","p")], as.numeric)

# -----------------------------------------------------------------------------
# 6.3 Genome-wide Cox scan: SHH
# -----------------------------------------------------------------------------

cox_shh_list <- apply(expr_shh_cc, 1, function(row_expr) {
  df  <- data.frame(time = time_shh, status = status_shh, expr = as.numeric(row_expr))
  fit <- coxph(Surv(time, status) ~ expr, data = df)
  s   <- summary(fit)
  c(coef = s$coef[1,"coef"], se = s$coef[1,"se(coef)"],
    HR = s$coef[1,"exp(coef)"], p = s$coef[1,"Pr(>|z|)"])
})

results_shh <- as.data.frame(t(cox_shh_list))
rownames(results_shh) <- NULL
results_shh$ensembl_gene_id <- gene_ids_micro
results_shh$hgnc_symbol     <- hgnc_syms_micro
results_shh[, c("coef","se","HR","p")] <- lapply(results_shh[, c("coef","se","HR","p")], as.numeric)

# -----------------------------------------------------------------------------
# 6.4 Volcano plots
# -----------------------------------------------------------------------------

make_volcano <- function(results, title_str, top_n = 15) {
  dat <- results
  dat$log2HR  <- log2(dat$HR)
  dat$neglogp <- -log10(dat$p)
  dat$label   <- ifelse(is.na(dat$hgnc_symbol) | dat$hgnc_symbol == "",
                        dat$ensembl_gene_id, dat$hgnc_symbol)
  top_genes      <- head(dat[order(dat$p), ], top_n)
  dat$show_label <- dat$label %in% top_genes$label
  dat$category   <- "Not significant"
  dat$category[dat$p < 0.05 & dat$log2HR > 0] <- "Risk"
  dat$category[dat$p < 0.05 & dat$log2HR < 0] <- "Protective"
  dat$category <- factor(dat$category, levels = c("Risk","Protective","Not significant"))
  
  ggplot(dat, aes(x = log2HR, y = neglogp, colour = category)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_text_repel(data = dat[dat$show_label, ], aes(label = label),
                    size = 3, fontface = "italic", max.overlaps = 20, box.padding = 0.4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = c("Risk" = "#d73027", "Protective" = "#4575b4", "Not significant" = "grey70")) +
    labs(x = "log2 Hazard Ratio", y = "-log10(p-value)", title = title_str, colour = NULL) +
    theme_classic(base_size = 12) +
    theme(legend.position = "top")
}

volcano_g34 <- make_volcano(results_g34, "Group 3/4: genome-wide Cox results")
volcano_shh <- make_volcano(results_shh, "SHH: genome-wide Cox results")

# -----------------------------------------------------------------------------
# 6.5 Genome-wide forest plots (top 10 by p-value)
# -----------------------------------------------------------------------------

make_forest_plot <- function(results, title_str, topN = 10) {
  dat <- head(results[order(results$p), ], topN)
  dat$HR_low  <- exp(dat$coef - 1.96 * dat$se)
  dat$HR_high <- exp(dat$coef + 1.96 * dat$se)
  dat$label   <- ifelse(is.na(dat$hgnc_symbol) | dat$hgnc_symbol == "",
                        dat$ensembl_gene_id, dat$hgnc_symbol)
  dat$direction <- ifelse(dat$HR > 1, "Risk", "Protective")
  dat$direction <- factor(dat$direction, levels = c("Risk","Protective"))
  dat$label     <- factor(dat$label, levels = rev(dat$label[order(dat$HR)]))
  
  ggplot(dat, aes(x = label, y = HR, colour = direction)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = HR_low, ymax = HR_high), width = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    coord_flip() +
    scale_y_log10() +
    scale_colour_manual(values = c("Risk" = "#d73027", "Protective" = "#4575b4")) +
    labs(x = "Gene", y = "Hazard ratio (log scale)", title = title_str, colour = NULL) +
    theme_classic(base_size = 12) +
    theme(axis.text.y = element_text(face = "italic"), legend.position = "top")
}

# Volcano + forest side by side for each subgroup
volcano_g34 + make_forest_plot(results_g34, "Group 3/4: Top 10 prognostic genes")
volcano_shh + make_forest_plot(results_shh, "SHH: Top 10 prognostic genes")

########################################################################################################################
#
#   SECTION 7: PENALISED COX REGRESSION — GROUP 3/4
#   Input:  expr_micro, meta (microarray)
#   Output: cv_fits_g34, results_enet_g34, risk scores, KM, oncoprint
#
########################################################################################################################

# -----------------------------------------------------------------------------
# 7.1 Prepare data
# -----------------------------------------------------------------------------

dat <- meta[!is.na(meta$OS..years.) & !is.na(meta$Dead), ]
dat$Dead <- as.integer(dat$Dead)
stopifnot(all(dat$Dead %in% c(0,1)))

sample_cols_glm <- grep("^MB_SubtypeStudy_", names(expr_df_micro), value = TRUE)
X_gxs <- as.matrix(sapply(expr_df_micro[, sample_cols_glm, drop = FALSE], as.numeric))
rownames(X_gxs) <- gene_ids_micro
x_all_micro     <- t(X_gxs)
rownames(x_all_micro) <- sample_cols_glm

dat_g34_micro        <- subset(dat, Subgroup %in% c("Group3","Group4"))
common_ids_g34_micro <- intersect(dat_g34_micro[[id_col]], rownames(x_all_micro))
stopifnot(length(common_ids_g34_micro) > 10)

dat_glm_g34 <- dat_g34_micro[match(common_ids_g34_micro, dat_g34_micro[[id_col]]), ]
x_glm_g34   <- x_all_micro[common_ids_g34_micro, , drop = FALSE]
x_glm_g34   <- x_glm_g34[, colSums(is.na(x_glm_g34)) == 0, drop = FALSE]

cat("G3/4 samples:", nrow(dat_glm_g34), "| Genes:", ncol(x_glm_g34), "\n")

y_g34  <- with(dat_glm_g34, Surv(OS..years., Dead))
set.seed(1)
K      <- 10
foldid <- sample(rep(1:K, length.out = nrow(x_glm_g34)))
alphas <- c(0, 0.5, 1)

# -----------------------------------------------------------------------------
# 7.2 Cross-validation across alpha values
# -----------------------------------------------------------------------------

cat("G3/4: Starting cv.glmnet |", format(Sys.time()), "\n")
cv_fits_g34 <- pblapply(alphas, function(a) {
  cv.glmnet(x_glm_g34, y_g34, family = "cox", alpha = a,
            nfolds = K, foldid = foldid, type.measure = "deviance", trace.it = 1)
})
names(cv_fits_g34) <- paste0("alpha_", alphas)
cat("G3/4: Finished |", format(Sys.time()), "\n")

cv_summary_g34 <- data.frame(
  alpha      = alphas,
  lambda_min = sapply(cv_fits_g34, function(f) f$lambda.min),
  lambda_1se = sapply(cv_fits_g34, function(f) f$lambda.1se),
  cvm_min    = sapply(cv_fits_g34, function(f) min(f$cvm))
)
print(cv_summary_g34)

# CV curves
names(cv_fits_g34) <- c("alpha_0", "alpha_0.5", "alpha_1")

# Then run the plot
yl <- range(unlist(lapply(cv_fits_g34, function(f) f$cvm)), na.rm = TRUE)
par(mfrow = c(1,3), mar = c(4,4,2,1), oma = c(0,0,3,0))
plot(cv_fits_g34[["alpha_0"]],   main = "", ylim = yl); mtext("Alpha = 0 (ridge)",         side=3, line=-1, cex=0.9)
plot(cv_fits_g34[["alpha_0.5"]], main = "", ylim = yl); mtext("Alpha = 0.5 (elastic net)", side=3, line=-1, cex=0.9)
plot(cv_fits_g34[["alpha_1"]],   main = "", ylim = yl); mtext("Alpha = 1 (lasso)",         side=3, line=-1, cex=0.9)
mtext("G3/4: 10-fold CV curves", outer = TRUE, cex = 1.1, line = 1)
par(mfrow = c(1,1))


# -----------------------------------------------------------------------------
# 7.3 Extract genes and refit Cox (alpha = 0.5, lambda.min)
# -----------------------------------------------------------------------------

best_fit_g34   <- cv_fits_g34[["alpha_0.5"]]
coef_mat_g34   <- coef(best_fit_g34, s = "lambda.min")
selected_g34   <- rownames(coef_mat_g34)[coef_mat_g34[,1] != 0]
cat("Selected", length(selected_g34), "genes\n")

y2        <- as.character(y_g34)
OS_Status <- rep(1, length(y2)); OS_Status[grep("\\+", y2)] <- 0
OS_Time   <- as.numeric(gsub("\\+","", y2))

cox_data_g34    <- data.frame(time = OS_Time, event = OS_Status, x_glm_g34[, selected_g34, drop = FALSE])
cox_formula_g34 <- as.formula(paste("Surv(time, event) ~", paste(selected_g34, collapse = " + ")))
cox_model_g34   <- coxph(cox_formula_g34, data = cox_data_g34)
summary(cox_model_g34)

final_coefs_g34   <- coef(cox_model_g34)
risk_score_g34_lp <- predict(cox_model_g34, type = "lp")
write.csv(as.matrix(x_glm_g34[, names(final_coefs_g34)]) %*% final_coefs_g34, file = "G3G4_predictor.csv")

# -----------------------------------------------------------------------------
# 7.4 Elastic net forest plot (top 10 by coefficient magnitude)
# -----------------------------------------------------------------------------

get_gene_symbols <- function(ids) {
  clean <- gsub("\\..*", "", ids)
  syms  <- mapIds(org.Hs.eg.db, keys = clean, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  syms  <- as.character(syms)
  syms[is.na(syms)] <- ids[is.na(syms)]
  syms
}

lambda_idx  <- which.min(abs(best_fit_g34$lambda - best_fit_g34$lambda.min))
betas_g34   <- as.matrix(best_fit_g34$glmnet.fit$beta[, lambda_idx, drop = FALSE])
betas_g34   <- betas_g34[betas_g34[,1] != 0, , drop = FALSE]
top_g34     <- head(rownames(betas_g34)[order(abs(betas_g34[,1]), decreasing = TRUE)], 10)
names_g34   <- get_gene_symbols(top_g34)

keep_cc <- !is.na(dat_glm_g34$OS..years.) & !is.na(dat_glm_g34$Dead)

fit_enet_g34 <- function(gene_id, gene_name) {
  df <- data.frame(
    time       = dat_glm_g34$OS..years.[keep_cc],
    status     = as.integer(dat_glm_g34$Dead[keep_cc]),
    expression = x_glm_g34[keep_cc, gene_id]
  )
  fit <- coxph(Surv(time, status) ~ expression, data = df)
  s   <- summary(fit)
  data.frame(symbol = gene_name, HR = s$coef[1,"exp(coef)"],
             lower_95 = s$conf.int[1,"lower .95"], upper_95 = s$conf.int[1,"upper .95"],
             pval = s$coef[1,"Pr(>|z|)"], stringsAsFactors = FALSE)
}

results_enet_g34 <- do.call(rbind, Map(fit_enet_g34, top_g34, names_g34))
results_enet_g34$padj <- p.adjust(results_enet_g34$pval, method = "BH")
results_enet_g34 <- results_enet_g34[order(results_enet_g34$HR), ]
write.csv(results_enet_g34, "G34_top10_forest_df.csv", row.names = FALSE)

make_enet_forest <- function(results, title_str) {
  dat <- results
  dat$direction <- factor(ifelse(dat$HR > 1, "Risk", "Protective"), levels = c("Risk","Protective"))
  dat$symbol    <- factor(dat$symbol, levels = dat$symbol[order(dat$HR)])
  ggplot(dat, aes(x = HR, y = symbol, colour = direction)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = lower_95, xmax = upper_95), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
    scale_x_log10() +
    scale_colour_manual(values = c("Risk" = "#d73027", "Protective" = "#4575b4")) +
    labs(x = "Hazard ratio (log scale)", y = "Gene", title = title_str, colour = NULL) +
    theme_classic(base_size = 12) +
    theme(axis.text.y = element_text(face = "italic"), legend.position = "top")
}

forest_enet_g34 <- make_enet_forest(results_enet_g34, "G3/4: Top 10 prognostic genes (elastic net)")
print(forest_enet_g34)

# -----------------------------------------------------------------------------
# 7.5 Risk tertiles and KM
# -----------------------------------------------------------------------------

rs_g34   <- as.numeric(risk_score_g34_lp)
cuts_g34 <- quantile(rs_g34, probs = c(1/3, 2/3), na.rm = TRUE, type = 7)

dat_glm_g34$risk_score   <- rs_g34
dat_glm_g34$risk_tertile <- cut(rs_g34, breaks = c(-Inf, cuts_g34, Inf),
                                labels = c("Low","Mid","High"), include.lowest = TRUE)

message("G3/4 tertile sizes:"); print(table(dat_glm_g34$risk_tertile))

km_fit_g34  <- survfit(Surv(OS..years., Dead) ~ risk_tertile, data = dat_glm_g34)
lr_test_g34 <- survdiff(Surv(OS..years., Dead) ~ risk_tertile, data = dat_glm_g34)
lr_p_g34    <- 1 - pchisq(lr_test_g34$chisq, df = length(lr_test_g34$n) - 1)

ggsurvplot(km_fit_g34, data = dat_glm_g34,
           palette = c("#2166ac","#f4a582","#d6604d"), conf.int = TRUE, conf.int.alpha = 0.12,
           risk.table = TRUE, risk.table.height = 0.25,
           xlab = "Time (years)", ylab = "Overall survival probability",
           title = "Kaplan-Meier survival by risk tertile \u2014 G3/4 cohort",
           legend.title = "Risk group", legend.labs = c("Low","Mid","High"),
           pval = FALSE, ggtheme = theme_classic(base_size = 13), xlim = c(0,10)) |>
  (\(p) { p$plot <- p$plot + annotate("text", x = max(dat_glm_g34$OS..years., na.rm=TRUE)*0.65,
                                      y = 0.92, label = paste0("Log-rank p = ", format.pval(lr_p_g34, digits=2, eps=0.001)),
                                      size = 4, hjust = 0); p })()

# -----------------------------------------------------------------------------
# 7.6 Clinical annotation figure (oncoprint + ridge plots)
# -----------------------------------------------------------------------------

dat_glm_g34$MYC  <- ifelse(dat_glm_g34$Study_ID %in% myc_amp_ids,  1, 0)
dat_glm_g34$MYCN <- ifelse(dat_glm_g34$Study_ID %in% mycn_amp_ids, 1, 0)

dat_plot      <- dat_glm_g34[order(dat_glm_g34$risk_score), ]
dat_plot$Dead <- as.numeric(dat_plot$Dead)
dat_plot$Gender    <- ifelse(is.na(dat_plot$Gender),    "Unknown", as.character(dat_plot$Gender))
dat_plot$histology <- ifelse(is.na(dat_plot$histology) | dat_plot$histology == "", "Unknown", dat_plot$histology)
dat_plot$Subgroup  <- ifelse(is.na(dat_plot$Subgroup),  "Unknown", as.character(dat_plot$Subgroup))
dat_plot$Met.status..1.Met..0.M0. <- ifelse(is.na(dat_plot$Met.status..1.Met..0.M0.), "Unknown", as.character(dat_plot$Met.status..1.Met..0.M0.))
dat_plot$MYC  <- as.character(dat_plot$MYC)
dat_plot$MYCN <- as.character(dat_plot$MYCN)

tertile_col_fun  <- c("Low"="#2166ac","Mid"="#f4a582","High"="#d6604d")
gender_col_fun   <- c("M"="white","F"="black","Unknown"="lightgrey")
subgroup_col_fun <- c("Group3"="#FFFF00","Group4"="#008000","Unknown"="lightgrey")
hist_col_fun     <- c("Classic"="white","Desmoplastic"="#fc8d59","LCA"="black","MBEN"="red4","Unknown"="lightgrey")
mstage_col_fun   <- c("0"="white","1"="black","Unknown"="lightgrey")
myc_col_fun      <- c("0"="white","1"="#e31a1c")
mycn_col_fun     <- c("0"="white","1"="#ff7f00")

age_dot_anno <- HeatmapAnnotation(
  Age = anno_points(pmin(dat_plot$Age, 20), ylim = c(0,20),
                    pch = ifelse(is.na(dat_plot$Age), 4, 20),
                    gp  = gpar(col = ifelse(is.na(dat_plot$Age), "red", "grey")),
                    axis_param = list(at = c(0,5,10,15,20)), baseline = 3,
                    baseline_gp = gpar(col = "grey", lty = 2, lwd = 1)),
  annotation_name_side = "left", annotation_height = unit(2, "cm"))

OS_g34 <- pmin(dat_plot$OS..years., 10)
pch_g34 <- ifelse(dat_plot$Dead==1 & dat_plot$OS..years.<10, 4,
                  ifelse(dat_plot$Dead==0 & dat_plot$OS..years.<10, 20, 2))
OS_dot_anno <- HeatmapAnnotation(
  OS = anno_points(OS_g34, pch = pch_g34,
                   gp = gpar(col = ifelse(dat_plot$Dead==1, "firebrick3", "lightgrey")),
                   ylim = c(0,10), axis_param = list(at = c(0,5,10)),
                   baseline_gp = gpar(col = "lightgrey", lty = 2, lwd = 1)),
  annotation_name_side = "left", annotation_height = unit(2, "cm"))

cat_anno <- HeatmapAnnotation(
  `Risk tertile` = anno_simple(as.character(dat_plot$risk_tertile), col = tertile_col_fun,  border = TRUE),
  Gender         = anno_simple(dat_plot$Gender,                     col = gender_col_fun,   border = TRUE),
  Subgroup       = anno_simple(dat_plot$Subgroup,                   col = subgroup_col_fun, border = TRUE),
  Histology      = anno_simple(dat_plot$histology,                  col = hist_col_fun,     border = TRUE),
  `M-stage`      = anno_simple(dat_plot$Met.status..1.Met..0.M0.,  col = mstage_col_fun,   border = TRUE),
  MYC            = anno_simple(dat_plot$MYC,                        col = myc_col_fun,      border = TRUE),
  MYCN           = anno_simple(dat_plot$MYCN,                       col = mycn_col_fun,     border = TRUE),
  annotation_name_side = "left")

ht_g34 <- Heatmap(matrix(dat_plot$risk_score, nrow=1), name = "Risk score",
                  col = colorRamp2(c(min(dat_plot$risk_score), 0, max(dat_plot$risk_score)), c("#2166ac","white","#d6604d")),
                  border = TRUE, show_column_names = FALSE, cluster_columns = FALSE, cluster_rows = FALSE,
                  row_names_side = "left", height = unit(0.6, "cm"), show_heatmap_legend = TRUE,
                  top_annotation = c(cat_anno, age_dot_anno, OS_dot_anno))

lgd_tertile  <- Legend(labels=c("Low","Mid","High"), title="Risk tertile", legend_gp=gpar(fill=c("#2166ac","#f4a582","#d6604d")), border=TRUE)
lgd_subgroup <- Legend(labels=c("Group3","Group4","Unknown"), title="Subgroup", legend_gp=gpar(fill=c("#FFFF00","#008000","lightgrey")), border=TRUE)
lgd_gender   <- Legend(labels=c("Male","Female","Unknown"), title="Gender", legend_gp=gpar(fill=c("white","black","lightgrey")), border=TRUE)
lgd_hist     <- Legend(labels=c("Classic","Desmoplastic","LCA","MBEN","Unknown"), title="Histology", legend_gp=gpar(fill=c("white","#fc8d59","black","red4","lightgrey")), border=TRUE)
lgd_mstage   <- Legend(labels=c("M0","M+","Unknown"), title="M-stage", legend_gp=gpar(fill=c("white","black","lightgrey")), border=TRUE)
lgd_myc      <- Legend(labels=c("Not amplified","Amplified"), title="MYC",  legend_gp=gpar(fill=c("white","#e31a1c")), border=TRUE)
lgd_mycn     <- Legend(labels=c("Not amplified","Amplified"), title="MYCN", legend_gp=gpar(fill=c("white","#ff7f00")), border=TRUE)

fmt_p <- function(p) {
  if (p < 2.2e-16) "p < 2.2e-16"
  else if (p < 0.001) paste0("p = ", formatC(p, format="e", digits=2))
  else paste0("p = ", round(p, 3))
}

rs_range_g34 <- range(dat_glm_g34$risk_score, na.rm=TRUE)
ridge_theme  <- theme_bw() + theme(legend.position="none", panel.grid=element_blank())

make_ridge <- function(df, x_col, group_col, fill_vals, x_lab, group_label, rs_range, test_fn = "t") {
  toPlot <- data.frame(risk_score = df[[x_col]], group = as.character(df[[group_col]])) %>%
    filter(!is.na(group), group != "Unknown", group != "")
  toPlot$group <- factor(toPlot$group)
  
  p_val <- if (test_fn == "t") {
    fmt_p(t.test(risk_score ~ group, data = toPlot)$p.value)
  } else {
    fmt_p(anova(aov(risk_score ~ group, data = toPlot))$`Pr(>F)`[1])
  }
  
  ggplot(toPlot, aes(x = risk_score, y = group)) +
    geom_density_ridges(aes(fill = group), scale = 1.5, rel_min_height = 0.001, alpha = 0.45) +
    scale_fill_manual(values = fill_vals) +
    xlim(rs_range) + ridge_theme +
    labs(x = x_lab, y = group_label) +
    annotate("text", x = rs_range[2] - 0.1,
             y = length(levels(toPlot$group)) + 1.3,
             label = p_val, hjust = 1, size = 3.5)
}
A <- make_ridge(dat_glm_g34, "risk_score", "Subgroup",                c("Group3"="#FFFF00","Group4"="#008000"),        "Risk score", "Subgroup",   rs_range_g34, "t")
B <- make_ridge(dat_glm_g34, "risk_score", "histology",               c("Classic"="grey80","Desmoplastic"="#fc8d59","LCA"="black","MBEN"="red4"), "Risk score", "Histology", rs_range_g34, "anova")
C <- make_ridge(dat_glm_g34, "risk_score", "Met.status..1.Met..0.M0.", c("0"="grey80","1"="black"),                    "Risk score", "M-stage",   rs_range_g34, "t")
D <- make_ridge(dat_glm_g34, "risk_score", "Gender",                  c("F"="#d73027","M"="#4575b4"),                  "Risk score", "Gender",    rs_range_g34, "t")
E <- make_ridge(dat_glm_g34, "risk_score", "MYC",                     c("0"="grey80","1"="#e31a1c"),                   "Risk score", "MYC",       rs_range_g34, "t")

toPlot_F <- data.frame(risk_score = dat_glm_g34$risk_score, group = as.character(dat_glm_g34$MYCN)) %>% filter(!is.na(group))
p_F      <- fmt_p(t.test(risk_score ~ group, data = toPlot_F)$p.value)
F_plot   <- ggplot(toPlot_F, aes(risk_score, group)) +
  geom_density_ridges(aes(fill=group), scale=1.5, rel_min_height=0.001, alpha=0.45) +
  scale_fill_manual(values=c("0"="grey80","1"="#ff7f00")) +
  scale_y_discrete(labels=c("0"="Not amplified","1"="Amplified")) +
  xlim(rs_range_g34) + ridge_theme + labs(x="Risk score", y="MYCN") +
  annotate("text", x=rs_range_g34[2]-0.1, y=3.3, label=p_F, hjust=1, size=3.5)

pdf("G34_clinical_annotation.pdf", width=14, height=8)
draw(ht_g34, heatmap_legend_side="right", annotation_legend_side="right",
     annotation_legend_list=list(lgd_tertile, lgd_subgroup, lgd_gender, lgd_hist, lgd_mstage, lgd_myc, lgd_mycn),
     padding=unit(c(2,30,2,2),"mm"),
     column_title="Group 3/4 medulloblastoma: transcriptomic risk score and clinical annotation",
     column_title_gp=gpar(fontsize=13, fontface="bold"))
decorate_annotation("Age", { grid.lines(x=unit(c(0,1),"npc"), y=unit(c(3,3),"native"), gp=gpar(col="black",lty="dotted",lwd=1)) })
decorate_annotation("OS",  { grid.lines(x=unit(c(0,1),"npc"), y=unit(c(5,5),"native"), gp=gpar(col="black",lty="dotted",lwd=1)) })
grid.newpage()
print(ggarrange(A, B, C, D, E, F_plot, ncol=2, nrow=3, labels=c("A","B","C","D","E","F")))
dev.off()

########################################################################################################################
#
#   SECTION 8: PENALISED COX REGRESSION — SHH
#   Input:  expr_micro, meta (microarray)
#   Output: cv_fits_shh, results_enet_shh, risk scores, KM, oncoprint
#
########################################################################################################################

# -----------------------------------------------------------------------------
# 8.1 Prepare data
# -----------------------------------------------------------------------------

dat_shh_micro        <- subset(dat, Subgroup == "SHH")
common_ids_shh_micro <- intersect(dat_shh_micro[[id_col]], rownames(x_all_micro))
stopifnot(length(common_ids_shh_micro) > 10)

dat_glm_shh <- dat_shh_micro[match(common_ids_shh_micro, dat_shh_micro[[id_col]]), ]
x_glm_shh   <- x_all_micro[common_ids_shh_micro, , drop = FALSE]
x_glm_shh   <- x_glm_shh[, colSums(is.na(x_glm_shh)) == 0, drop = FALSE]
stopifnot(identical(dat_glm_shh[[id_col]], rownames(x_glm_shh)))

cat("SHH samples:", nrow(dat_glm_shh), "| Genes:", ncol(x_glm_shh), "\n")

y_shh      <- with(dat_glm_shh, Surv(OS..years., Dead))
set.seed(1)
foldid_shh <- sample(rep(1:K, length.out = nrow(x_glm_shh)))

# -----------------------------------------------------------------------------
# 8.2 Cross-validation across alpha values
# -----------------------------------------------------------------------------

cat("SHH: Starting cv.glmnet |", format(Sys.time()), "\n")
cv_fits_shh <- pblapply(alphas, function(a) {
  cv.glmnet(x_glm_shh, y_shh, family="cox", alpha=a,
            nfolds=K, foldid=foldid_shh, type.measure="deviance", trace.it=1)
})
names(cv_fits_shh) <- paste0("alpha_", alphas)
cat("SHH: Finished |", format(Sys.time()), "\n")

cv_summary_shh <- data.frame(
  alpha      = alphas,
  lambda_min = sapply(cv_fits_shh, function(f) f$lambda.min),
  lambda_1se = sapply(cv_fits_shh, function(f) f$lambda.1se),
  cvm_min    = sapply(cv_fits_shh, function(f) min(f$cvm))
)
print(cv_summary_shh)

yl_shh <- range(unlist(lapply(cv_fits_shh, function(f) f$cvm)), na.rm=TRUE)
par(mfrow=c(1,3), mar=c(4,4,2,1), oma=c(0,0,3,0))
plot(cv_fits_shh[["alpha_0"]],   main="", ylim=yl_shh); mtext("Alpha = 0 (ridge)",         side=3, line=-1, cex=0.9)
plot(cv_fits_shh[["alpha_0.5"]], main="", ylim=yl_shh); mtext("Alpha = 0.5 (elastic net)", side=3, line=-1, cex=0.9)
plot(cv_fits_shh[["alpha_1"]],   main="", ylim=yl_shh); mtext("Alpha = 1 (lasso)",         side=3, line=-1, cex=0.9)
mtext("SHH: 10-fold CV curves", outer=TRUE, cex=1.1, line=1)
par(mfrow=c(1,1))

# -----------------------------------------------------------------------------
# 8.3 Extract genes and refit Cox (alpha = 0.5, lambda.min)
# -----------------------------------------------------------------------------

best_fit_shh  <- cv_fits_shh[["alpha_0.5"]]
coef_mat_shh  <- coef(best_fit_shh, s = "lambda.min")
selected_shh  <- rownames(coef_mat_shh)[coef_mat_shh[,1] != 0]
cat("SHH: Selected", length(selected_shh), "genes\n")
if (length(selected_shh) == 0) stop("SHH: No genes selected at lambda.min")

cox_data_shh    <- data.frame(time=dat_glm_shh$OS..years., event=dat_glm_shh$Dead, x_glm_shh[, selected_shh, drop=FALSE], check.names=FALSE)
cox_model_shh   <- coxph(Surv(time, event) ~ ., data = cox_data_shh)
summary(cox_model_shh)

final_coefs_shh   <- coef(cox_model_shh)
risk_score_shh_lp <- drop(predict(cox_model_shh, type = "lp"))
write.csv(as.matrix(x_glm_shh[, names(final_coefs_shh)]) %*% final_coefs_shh, file = "SHH_predictor.csv")

# -----------------------------------------------------------------------------
# 8.4 Elastic net forest plot (top 10 by coefficient magnitude)
# -----------------------------------------------------------------------------

lambda_idx_shh <- which.min(abs(best_fit_shh$lambda - best_fit_shh$lambda.min))
betas_shh      <- as.matrix(best_fit_shh$glmnet.fit$beta[, lambda_idx_shh, drop=FALSE])
betas_shh      <- betas_shh[betas_shh[,1] != 0, , drop=FALSE]
top_shh        <- head(rownames(betas_shh)[order(abs(betas_shh[,1]), decreasing=TRUE)], 10)
names_shh      <- get_gene_symbols(top_shh)

fit_enet_shh <- function(gene_id, gene_name) {
  df <- data.frame(time=dat_glm_shh$OS..years., status=dat_glm_shh$Dead, expression=as.numeric(x_glm_shh[, gene_id]))
  df <- na.omit(df)
  fit <- coxph(Surv(time, status) ~ expression, data = df)
  s   <- summary(fit)
  data.frame(symbol=gene_name, HR=s$coef[1,"exp(coef)"],
             lower_95=s$conf.int[1,"lower .95"], upper_95=s$conf.int[1,"upper .95"],
             pval=s$coef[1,"Pr(>|z|)"], stringsAsFactors=FALSE)
}

results_enet_shh <- do.call(rbind, Map(fit_enet_shh, top_shh, names_shh))
results_enet_shh$padj <- p.adjust(results_enet_shh$pval, method="BH")
results_enet_shh <- results_enet_shh[order(results_enet_shh$HR), ]
write.csv(results_enet_shh, "SHH_top10_forest_df.csv", row.names=FALSE)

forest_enet_shh <- make_enet_forest(results_enet_shh, "SHH: Top 10 prognostic genes (elastic net)")
print(forest_enet_shh)

# Side-by-side elastic net forest plots
forest_enet_g34 + forest_enet_shh

# -----------------------------------------------------------------------------
# 8.5 Risk tertiles and KM
# -----------------------------------------------------------------------------

rs_shh   <- as.numeric(risk_score_shh_lp)
cuts_shh <- quantile(rs_shh, probs=c(1/3,2/3), na.rm=TRUE, type=7)

dat_glm_shh$risk_score   <- rs_shh
dat_glm_shh$risk_tertile <- cut(rs_shh, breaks=c(-Inf,cuts_shh,Inf),
                                labels=c("Low","Mid","High"), include.lowest=TRUE)

message("SHH tertile sizes:"); print(table(dat_glm_shh$risk_tertile))

km_fit_shh  <- survfit(Surv(OS..years., Dead) ~ risk_tertile, data=dat_glm_shh)
lr_test_shh <- survdiff(Surv(OS..years., Dead) ~ risk_tertile, data=dat_glm_shh)
lr_p_shh    <- 1 - pchisq(lr_test_shh$chisq, df=length(lr_test_shh$n)-1)

ggsurvplot(km_fit_shh, data=dat_glm_shh,
           palette=c("#2166ac","#f4a582","#d6604d"), conf.int=TRUE, conf.int.alpha=0.12,
           risk.table=TRUE, risk.table.height=0.25,
           xlab="Time (years)", ylab="Overall survival probability",
           title="Kaplan-Meier survival by risk tertile \u2014 SHH cohort",
           legend.title="Risk group", legend.labs=c("Low","Mid","High"),
           pval=FALSE, ggtheme=theme_classic(base_size=13), xlim=c(0,10)) |>
  (\(p) { p$plot <- p$plot + annotate("text", x=0.5, y=0.15,
                                      label=paste0("Log-rank p = ", format.pval(lr_p_shh, digits=2, eps=0.001)),
                                      size=4, hjust=0); p })()

# -----------------------------------------------------------------------------
# 8.6 Clinical annotation figure (oncoprint + ridge plots)
# -----------------------------------------------------------------------------

dat_glm_shh$MYC  <- ifelse(dat_glm_shh$Study_ID %in% myc_amp_ids,  1, 0)
dat_glm_shh$MYCN <- ifelse(dat_glm_shh$Study_ID %in% mycn_amp_ids, 1, 0)

dat_plot_shh      <- dat_glm_shh[order(dat_glm_shh$risk_score), ]
dat_plot_shh$Dead <- as.numeric(dat_plot_shh$Dead)
dat_plot_shh$Gender    <- ifelse(is.na(dat_plot_shh$Gender)    | dat_plot_shh$Gender    == "", "Unknown", as.character(dat_plot_shh$Gender))
dat_plot_shh$histology <- ifelse(is.na(dat_plot_shh$histology) | dat_plot_shh$histology == "", "Unknown", dat_plot_shh$histology)
dat_plot_shh$Subtype   <- ifelse(is.na(dat_plot_shh$Subtype)   | dat_plot_shh$Subtype   == "", "Unknown", dat_plot_shh$Subtype)
dat_plot_shh$Met.status..1.Met..0.M0. <- ifelse(is.na(dat_plot_shh$Met.status..1.Met..0.M0.), "Unknown", as.character(dat_plot_shh$Met.status..1.Met..0.M0.))
dat_plot_shh$MYC  <- as.character(dat_plot_shh$MYC)
dat_plot_shh$MYCN <- as.character(dat_plot_shh$MYCN)

tertile_col_shh <- c("Low"="#2166ac","Mid"="#f4a582","High"="#d6604d")
subtype_col_fun <- c("SHH_alpha"="#7b2d8b","SHH_beta"="#e7298a","SHH_gamma"="#66a61e","SHH_delta"="#e6ab02","Unknown"="lightgrey")
mycn_col_fun2   <- c("0"="white","1"="#ff7f00")

age_dot_shh <- HeatmapAnnotation(
  Age = anno_points(pmin(dat_plot_shh$Age,20), ylim=c(0,20),
                    pch=ifelse(is.na(dat_plot_shh$Age),4,20),
                    gp=gpar(col=ifelse(is.na(dat_plot_shh$Age),"red","grey")),
                    axis_param=list(at=c(0,5,10,15,20)), baseline=3, baseline_gp=gpar(col="grey",lty=2,lwd=1)),
  annotation_name_side="left", annotation_height=unit(2,"cm"))

OS_shh  <- pmin(dat_plot_shh$OS..years., 10)
pch_shh <- ifelse(dat_plot_shh$Dead==1 & dat_plot_shh$OS..years.<10, 4,
                  ifelse(dat_plot_shh$Dead==0 & dat_plot_shh$OS..years.<10, 20, 2))
OS_dot_shh <- HeatmapAnnotation(
  OS = anno_points(OS_shh, pch=pch_shh,
                   gp=gpar(col=ifelse(dat_plot_shh$Dead==1,"firebrick3","lightgrey")),
                   ylim=c(0,10), axis_param=list(at=c(0,5,10)), baseline_gp=gpar(col="lightgrey",lty=2,lwd=1)),
  annotation_name_side="left", annotation_height=unit(2,"cm"))

cat_anno_shh <- HeatmapAnnotation(
  `Risk tertile` = anno_simple(as.character(dat_plot_shh$risk_tertile), col=tertile_col_shh,  border=TRUE),
  Gender         = anno_simple(dat_plot_shh$Gender,                     col=gender_col_fun,   border=TRUE),
  Subtype        = anno_simple(dat_plot_shh$Subtype,                    col=subtype_col_fun,  border=TRUE),
  Histology      = anno_simple(dat_plot_shh$histology,                  col=hist_col_fun,     border=TRUE),
  `M-stage`      = anno_simple(dat_plot_shh$Met.status..1.Met..0.M0.,  col=mstage_col_fun,   border=TRUE),
  MYCN           = anno_simple(dat_plot_shh$MYCN,                       col=mycn_col_fun2,    border=TRUE),
  annotation_name_side="left")

ht_shh <- Heatmap(matrix(dat_plot_shh$risk_score, nrow=1), name="Risk score",
                  col=colorRamp2(c(min(dat_plot_shh$risk_score),0,max(dat_plot_shh$risk_score)), c("#2166ac","white","#d6604d")),
                  border=TRUE, show_column_names=FALSE, cluster_columns=FALSE, cluster_rows=FALSE,
                  row_names_side="left", height=unit(0.6,"cm"), show_heatmap_legend=TRUE,
                  top_annotation=c(cat_anno_shh, age_dot_shh, OS_dot_shh))

lgd_tertile_shh <- Legend(labels=c("Low","Mid","High"), title="Risk tertile", legend_gp=gpar(fill=c("#2166ac","#f4a582","#d6604d")), border=TRUE)
lgd_subtype     <- Legend(labels=c("SHH_alpha","SHH_beta","SHH_gamma","SHH_delta","Unknown"), title="Subtype", legend_gp=gpar(fill=c("#7b2d8b","#e7298a","#66a61e","#e6ab02","lightgrey")), border=TRUE)
lgd_mycn_shh    <- Legend(labels=c("Not amplified","Amplified"), title="MYCN", legend_gp=gpar(fill=c("white","#ff7f00")), border=TRUE)

rs_range_shh <- range(dat_glm_shh$risk_score, na.rm=TRUE)

A_shh <- make_ridge(dat_glm_shh, "risk_score", "histology",               c("Classic"="grey80","Desmoplastic"="#fc8d59","LCA"="black","MBEN"="red4"), "Risk score", "Histology", rs_range_shh, "anova")
B_shh <- make_ridge(dat_glm_shh, "risk_score", "Met.status..1.Met..0.M0.", c("0"="grey80","1"="black"),                    "Risk score", "M-stage",   rs_range_shh, "t")
C_shh <- make_ridge(dat_glm_shh, "risk_score", "Gender",                  c("F"="#d73027","M"="#4575b4"),                  "Risk score", "Gender",    rs_range_shh, "t")

toPlot_D_shh <- data.frame(risk_score=dat_glm_shh$risk_score, group=as.character(dat_glm_shh$MYCN)) %>% filter(!is.na(group))
p_D_shh      <- fmt_p(t.test(risk_score ~ group, data=toPlot_D_shh)$p.value)
D_shh <- ggplot(toPlot_D_shh, aes(risk_score, group)) +
  geom_density_ridges(aes(fill=group), scale=1.5, rel_min_height=0.001, alpha=0.45) +
  scale_fill_manual(values=c("0"="grey80","1"="#ff7f00")) +
  scale_y_discrete(labels=c("0"="Not amplified","1"="Amplified")) +
  xlim(rs_range_shh) + ridge_theme + labs(x="Risk score", y="MYCN") +
  annotate("text", x=rs_range_shh[2]-0.1, y=3.3, label=p_D_shh, hjust=1, size=3.5)

toPlot_E_shh <- data.frame(risk_score=dat_glm_shh$risk_score, group=dat_glm_shh$Subtype) %>%
  filter(!is.na(group), group!="Unknown") %>%
  mutate(group=factor(group, levels=c("SHH_alpha","SHH_beta","SHH_gamma","SHH_delta")))
p_E_shh <- fmt_p(anova(aov(risk_score ~ group, data=toPlot_E_shh))$`Pr(>F)`[1])
E_shh <- ggplot(toPlot_E_shh, aes(risk_score, group)) +
  geom_density_ridges(aes(fill=group), scale=1.5, rel_min_height=0.001, alpha=0.45) +
  scale_fill_manual(values=c("SHH_alpha"="#7b2d8b","SHH_beta"="#e7298a","SHH_gamma"="#66a61e","SHH_delta"="#e6ab02")) +
  xlim(rs_range_shh) + ridge_theme + labs(x="Risk score", y="Subtype") +
  annotate("text", x=rs_range_shh[2]-0.1, y=5.3, label=p_E_shh, hjust=1, size=3.5)

pdf("SHH_clinical_annotation.pdf", width=14, height=8)
draw(ht_shh, heatmap_legend_side="right", annotation_legend_side="right",
     annotation_legend_list=list(lgd_tertile_shh, lgd_subtype, lgd_gender, lgd_hist, lgd_mstage, lgd_mycn_shh),
     padding=unit(c(2,30,2,2),"mm"),
     column_title="SHH medulloblastoma: transcriptomic risk score and clinical annotation",
     column_title_gp=gpar(fontsize=13, fontface="bold"))
decorate_annotation("Age", { grid.lines(x=unit(c(0,1),"npc"), y=unit(c(3,3),"native"), gp=gpar(col="black",lty="dotted",lwd=1)) })
decorate_annotation("OS",  { grid.lines(x=unit(c(0,1),"npc"), y=unit(c(5,5),"native"), gp=gpar(col="black",lty="dotted",lwd=1)) })
grid.newpage()
print(ggarrange(A_shh, B_shh, C_shh, D_shh, E_shh, ncol=2, nrow=3, labels=c("A","B","C","D","E")))
dev.off()

########################################################################################################################
#
#   SECTION 9: METHYLOMIC VS TRANSCRIPTOMIC RISK SCORE COMPARISON
#   Input:  Ed's methylomic scores, G3G4_predictor.csv, SHH_predictor.csv
#   Output: Scatter plots with Pearson correlation
#
########################################################################################################################

methyl_g34  <- read.csv("EdGroup3_4_alpha0.5_RiskScores.csv")
methyl_shh  <- read.csv("EdSHH_alpha1_RiskScores.csv")

transcr_g34 <- read.csv("G3G4_predictor.csv", row.names=1) %>%
  rownames_to_column("Study_ID") %>%
  setNames(c("Study_ID", "transcr_score"))

transcr_shh <- read.csv("SHH_predictor.csv", row.names=1) %>%
  rownames_to_column("Study_ID") %>%
  setNames(c("Study_ID", "transcr_score"))

metadata    <- read.csv("taylor_pheno_forDan.csv")

meta_lookup <- metadata %>%
  dplyr::select(title, Study_ID, Subgroup) %>%
  dplyr::rename(GSM = title)

methyl_g34  <- methyl_g34 %>% mutate(GSM = sub("_.*", "", SampleID))
methyl_shh  <- methyl_shh  %>% mutate(GSM = sub("_.*", "", SampleID))

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

cat("G3/4 matched:", nrow(merged_g34), "| SHH matched:", nrow(merged_shh), "\n")

ct_g34 <- cor.test(merged_g34$methyl_score, merged_g34$transcr_score, method = "pearson")
p1 <- ggplot(merged_g34, aes(x = methyl_score, y = transcr_score, colour = Subgroup)) +
  geom_point(alpha = 0.75, size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey40", linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 4,
           label = paste0("r = ", round(ct_g34$estimate, 3),
                          "\np = ", format.pval(ct_g34$p.value, digits = 2, eps = 1e-3),
                          "\nn = ", nrow(merged_g34))) +
  scale_colour_manual(values = c("Group3" = "#FFD700", "Group4" = "#008000")) +
  labs(title = "G3/4: Methylomic vs Transcriptomic Risk Scores",
       x = "Methylomic Risk Score", y = "Transcriptomic Risk Score") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

ct_shh <- cor.test(merged_shh$methyl_score, merged_shh$transcr_score, method = "pearson")
p2 <- ggplot(merged_shh, aes(x = methyl_score, y = transcr_score, colour = Subgroup)) +
  geom_point(alpha = 0.75, size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey40", linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 4,
           label = paste0("r = ", round(ct_shh$estimate, 3),
                          "\np = ", format.pval(ct_shh$p.value, digits = 2, eps = 1e-3),
                          "\nn = ", nrow(merged_shh))) +
  scale_colour_manual(values = c("SHH" = "#FF6666")) +
  labs(title = "SHH: Methylomic vs Transcriptomic Risk Scores",
       x = "Methylomic Risk Score", y = "Transcriptomic Risk Score") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

print(p1 + p2)

########################################################################################################################
#
#   SECTION 10: RNA-SEQ CROSS-PLATFORM VALIDATION (E-MTAB-10767)
#   Input:  MB.vsd.txt, E-MTAB-10767.sdrf, microarray Cox models
#   Output: KM validation plots, oncoprints for G3/4 and SHH
#
########################################################################################################################

# -----------------------------------------------------------------------------
# 10.1 Load RNA-seq expression matrix
# -----------------------------------------------------------------------------

expr_rna      <- fread("MB.vsd.txt", data.table = FALSE)
gene_ids_rna  <- expr_rna[[1]]
expr_rna[[1]] <- NULL
expr_rna[]    <- lapply(expr_rna, as.numeric)
expr_mat_rna  <- as.matrix(expr_rna)
rownames(expr_mat_rna) <- gene_ids_rna
rownames(expr_mat_rna) <- sub("_.*$",   "", rownames(expr_mat_rna))
rownames(expr_mat_rna) <- sub("\\..*$", "", rownames(expr_mat_rna))

if (any(duplicated(rownames(expr_mat_rna)))) {
  expr_dt      <- as.data.table(expr_mat_rna, keep.rownames = "ENSG")
  expr_dt      <- expr_dt[, lapply(.SD, mean, na.rm = TRUE), by = ENSG]
  expr_mat_rna <- as.matrix(expr_dt[, -1])
  rownames(expr_mat_rna) <- expr_dt$ENSG
}
cat("RNA-seq:", nrow(expr_mat_rna), "genes x", ncol(expr_mat_rna), "samples\n")

# -----------------------------------------------------------------------------
# 10.2 Load SDRF metadata
# -----------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------
# 10.3 G3/4 validation — project microarray signature onto RNA-seq
# -----------------------------------------------------------------------------

dat_g34_rna    <- pheno %>% filter(Subgroup %in% c("Grp3", "Grp4", "Grp3/Grp4"))
common_rna_g34 <- dat_g34_rna$SampleID[dat_g34_rna$SampleID %in% colnames(expr_mat_rna)]

x_rna_g34 <- t(expr_mat_rna[, common_rna_g34, drop = FALSE])
rownames(x_rna_g34) <- common_rna_g34
colnames(x_rna_g34) <- rownames(expr_mat_rna)
x_rna_g34 <- x_rna_g34[, apply(x_rna_g34, 2, sd, na.rm = TRUE) > 0, drop = FALSE]

overlap_g34       <- intersect(names(final_coefs_g34), colnames(x_rna_g34))
rs_rna_g34        <- as.numeric(x_rna_g34[, overlap_g34] %*% final_coefs_g34[overlap_g34])
names(rs_rna_g34) <- rownames(x_rna_g34)

cat("G3/4 RNA-seq — genes in model:", length(names(final_coefs_g34)),
    "| found in RNA-seq:", length(overlap_g34), "\n")

pheno_rna_g34 <- sdrf %>%
  transmute(
    SampleID   = as.character(`Source Name`),
    Subgroup   = as.character(`Characteristics[subgroup]`),
    time       = as.numeric(`Characteristics[os time]`),
    status_raw = as.character(`Characteristics[os status]`)
  ) %>%
  mutate(event = as.numeric(status_raw)) %>%
  filter(Subgroup %in% c("Grp3", "Grp4"),
         !is.na(time), !is.na(event), time > 0,
         SampleID %in% names(rs_rna_g34))

pheno_rna_g34$risk_score <- rs_rna_g34[pheno_rna_g34$SampleID]
cuts_rna_g34             <- quantile(rs_rna_g34, probs = c(1/3, 2/3))
pheno_rna_g34$risk_group <- cut(pheno_rna_g34$risk_score,
                                breaks = c(-Inf, cuts_rna_g34, Inf),
                                labels = c("Low", "Mid", "High"))

# -----------------------------------------------------------------------------
# 10.4 G3/4 KM validation plot
# -----------------------------------------------------------------------------

km_rna_g34 <- survfit(Surv(time, event) ~ risk_group, data = pheno_rna_g34)
ggsurvplot(km_rna_g34, data = pheno_rna_g34,
           pval = TRUE, pval.coord = c(0.5, 0.15), conf.int = FALSE, risk.table = TRUE,
           xlim = c(0, 12), break.time.by = 2,
           palette      = c("#2166ac", "#f4a582", "#d6604d"),
           legend.labs  = c("Low", "Mid", "High"), legend.title = "Risk group",
           xlab         = "Time (years)", ylab = "Overall survival probability",
           title        = "RNA-seq G3/4 \u2014 Microarray-derived risk score validation",
           ggtheme      = theme_classic(), tables.theme = theme_cleantable(), fontsize = 4)

# -----------------------------------------------------------------------------
# 10.5 Build risk_out_rnaseq for G3/4 oncoprint
# -----------------------------------------------------------------------------

risk_out_rnaseq <- pheno_rna_g34 %>%
  mutate(risk_lp = risk_score, tertile = risk_group) %>%
  dplyr::select(SampleID, risk_lp, event, time, tertile)

risk_out_rnaseq$tertile <- as.character(risk_out_rnaseq$tertile)
risk_out_rnaseq$tertile[risk_out_rnaseq$tertile == "Mid"] <- "Intermediate"
risk_out_rnaseq$tertile <- factor(risk_out_rnaseq$tertile,
                                  levels = c("Low", "Intermediate", "High"))

# -----------------------------------------------------------------------------
# 10.6 G3/4 clinical metadata from SDRF
# -----------------------------------------------------------------------------

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
      ifelse(grepl("^0|0.*3|infant|toddler", x),                        "0-3",
             ifelse(grepl("3.*16|3-16|child|pediatric|paediatric", x),         "3-16",
                    ifelse(grepl(">.*16|adult|>16", x),                               ">16",
                           ifelse(x == "" | x == "na" | is.na(x),                            "Unknown", x))))
    }),
    Mstage = ifelse(tolower(trimws(Mstage_raw)) %in% c("true","1","yes","y"), "1",
                    ifelse(tolower(trimws(Mstage_raw)) %in% c("false","0","no","n"), "0", "Unknown")),
    MYC  = ifelse(tolower(trimws(MYC_raw))  %in% c("true","1","yes","y"), "1", "0"),
    MYCN = ifelse(tolower(trimws(MYCN_raw)) %in% c("true","1","yes","y"), "1", "0")
  )

# -----------------------------------------------------------------------------
# 10.7 G3/4 oncoprint — merge, clean, build heatmap
# -----------------------------------------------------------------------------

dat_plot <- risk_out_rnaseq %>%
  left_join(clinical_rnaseq, by = "SampleID") %>%
  arrange(risk_lp)

dat_plot$Gender   <- ifelse(is.na(dat_plot$Gender)   | dat_plot$Gender   == "", "Unknown", dat_plot$Gender)
dat_plot$Subgroup <- ifelse(is.na(dat_plot$Subgroup) | dat_plot$Subgroup == "", "Unknown", dat_plot$Subgroup)
dat_plot$Mstage   <- ifelse(is.na(dat_plot$Mstage)   | dat_plot$Mstage   == "", "Unknown", dat_plot$Mstage)
dat_plot$Subgroup <- dplyr::recode(dat_plot$Subgroup,
                                   "Grp3" = "Group 3", "Grp4" = "Group 4",
                                   "Grp3/Grp4" = "Unknown", "Unknown" = "Unknown")
dat_plot$MYC     <- as.character(dat_plot$MYC)
dat_plot$MYCN    <- as.character(dat_plot$MYCN)
dat_plot$tertile <- ifelse(is.na(dat_plot$tertile) | dat_plot$tertile == "", "Unknown",
                           as.character(dat_plot$tertile))
dat_plot$AgeGroup <- ifelse(is.na(dat_plot$Age) | dat_plot$Age == "", "Unknown", dat_plot$Age)
dat_plot <- dat_plot %>% filter(!AgeGroup %in% c(">16", "Unknown"))
cat("G3/4 RNA-seq oncoprint samples:", nrow(dat_plot), "\n")

# Colour palettes
tertile_col_rna  <- c("Low"="#2166ac","Intermediate"="#f4a582","High"="#d6604d","Unknown"="lightgrey")
gender_col_rna   <- c("M"="white","F"="black","Unknown"="lightgrey")
subgroup_col_rna <- c("Group 3"="#FFD700","Group 4"="#008000","Unknown"="lightgrey")
mstage_col_rna   <- c("0"="white","1"="black","Unknown"="lightgrey")
myc_col_rna      <- c("0"="white","1"="#e31a1c")
mycn_col_rna     <- c("0"="white","1"="#ff7f00")
age_col_rna      <- c("0-3"="#ffffcc","3-16"="#41b6c4")

OS_rna    <- pmin(dat_plot$time, 10)
pch_rna   <- ifelse(dat_plot$event==1 & dat_plot$time<10, 4,
                    ifelse(dat_plot$event==0 & dat_plot$time<10, 20, 2))

OS_dot_rna <- HeatmapAnnotation(
  OS = anno_points(OS_rna, pch = pch_rna,
                   gp          = gpar(col = ifelse(dat_plot$event==1, "firebrick3", "lightgrey")),
                   ylim        = c(0,10), axis_param = list(at = c(0,5,10)),
                   baseline_gp = gpar(col = "lightgrey", lty = 2, lwd = 1)),
  annotation_name_side = "left", annotation_height = unit(2, "cm"))

cat_anno_rna <- HeatmapAnnotation(
  `Risk tertile` = anno_simple(dat_plot$tertile,   col = tertile_col_rna,  border = TRUE),
  Gender         = anno_simple(dat_plot$Gender,    col = gender_col_rna,   border = TRUE),
  `Age group`    = anno_simple(dat_plot$AgeGroup,  col = age_col_rna,      border = TRUE),
  Subgroup       = anno_simple(dat_plot$Subgroup,  col = subgroup_col_rna, border = TRUE),
  `M-stage`      = anno_simple(dat_plot$Mstage,    col = mstage_col_rna,   border = TRUE),
  MYC            = anno_simple(dat_plot$MYC,       col = myc_col_rna,      border = TRUE),
  MYCN           = anno_simple(dat_plot$MYCN,      col = mycn_col_rna,     border = TRUE),
  annotation_name_side = "left")

rs_mat_rna <- matrix(dat_plot$risk_lp, nrow = 1)
colnames(rs_mat_rna) <- dat_plot$SampleID

ht_rnaseq <- Heatmap(rs_mat_rna, name = "Risk score",
                     col = colorRamp2(c(min(dat_plot$risk_lp, na.rm=TRUE), 0, max(dat_plot$risk_lp, na.rm=TRUE)),
                                      c("#2166ac","white","#d6604d")),
                     border = TRUE, show_column_names = FALSE, cluster_columns = FALSE, cluster_rows = FALSE,
                     row_names_side = "left", height = unit(0.6,"cm"), show_heatmap_legend = TRUE,
                     top_annotation = c(cat_anno_rna, OS_dot_rna))

lgd_tertile_rna  <- Legend(labels=c("Low","Intermediate","High","Unknown"), title="Risk tertile",
                           legend_gp=gpar(fill=c("#2166ac","#f4a582","#d6604d","lightgrey")), border=TRUE)
lgd_subgroup_rna <- Legend(labels=c("Group 3","Group 4","Unknown"), title="Subgroup",
                           legend_gp=gpar(fill=c("#FFD700","#008000","lightgrey")), border=TRUE)
lgd_gender_rna   <- Legend(labels=c("Male","Female","Unknown"), title="Gender",
                           legend_gp=gpar(fill=c("white","black","lightgrey")), border=TRUE)
lgd_age_rna      <- Legend(labels=c("0-3","3-16"), title="Age group",
                           legend_gp=gpar(fill=c("#ffffcc","#41b6c4")), border=TRUE)
lgd_mstage_rna   <- Legend(labels=c("M0","M+","Unknown"), title="M-stage",
                           legend_gp=gpar(fill=c("white","black","lightgrey")), border=TRUE)
lgd_myc_rna      <- Legend(labels=c("Not amplified","Amplified"), title="MYC",
                           legend_gp=gpar(fill=c("white","#e31a1c")), border=TRUE)
lgd_mycn_rna     <- Legend(labels=c("Not amplified","Amplified"), title="MYCN",
                           legend_gp=gpar(fill=c("white","#ff7f00")), border=TRUE)

# Ridge plots
fmt_p <- function(p) {
  if (is.na(p))     return("p = NA")
  if (p < 2.2e-16)  return("p < 2.2e-16")
  if (p < 0.001)    return(paste0("p = ", formatC(p, format="e", digits=2)))
  paste0("p = ", round(p, 3))
}

rs_range_rna <- range(dat_plot$risk_lp, na.rm = TRUE)
ridge_theme  <- theme_bw() + theme(legend.position = "none", panel.grid = element_blank())

toPlot_A <- data.frame(risk_score=dat_plot$risk_lp, group=dat_plot$Subgroup) %>% filter(group %in% c("Group 3","Group 4"))
p_A <- fmt_p(t.test(risk_score ~ group, data=toPlot_A)$p.value)
A <- ggplot(toPlot_A, aes(risk_score, group)) +
  geom_density_ridges(aes(fill=group), scale=1.5, rel_min_height=0.001, alpha=0.45) +
  scale_fill_manual(values=c("Group 3"="#FFD700","Group 4"="#008000")) +
  xlim(rs_range_rna) + ridge_theme + labs(x="Risk score", y="Subgroup") +
  annotate("text", x=rs_range_rna[2]-0.1, y=Inf, label=p_A, hjust=1, vjust=1.5, size=3.5)

toPlot_B <- data.frame(risk_score=dat_plot$risk_lp, group=dat_plot$Mstage) %>% filter(group != "Unknown")
p_B <- fmt_p(t.test(risk_score ~ group, data=toPlot_B)$p.value)
B <- ggplot(toPlot_B, aes(risk_score, group)) +
  geom_density_ridges(aes(fill=group), scale=1.5, rel_min_height=0.001, alpha=0.45) +
  scale_fill_manual(values=c("0"="grey80","1"="black")) +
  scale_y_discrete(labels=c("0"="M0","1"="M+")) +
  xlim(rs_range_rna) + ridge_theme + labs(x="Risk score", y="M-stage") +
  annotate("text", x=rs_range_rna[2]-0.1, y=Inf, label=p_B, hjust=1, vjust=1.5, size=3.5)

toPlot_C <- data.frame(risk_score=dat_plot$risk_lp, group=dat_plot$Gender) %>% filter(group != "Unknown")
p_C <- fmt_p(t.test(risk_score ~ group, data=toPlot_C)$p.value)
C <- ggplot(toPlot_C, aes(risk_score, group)) +
  geom_density_ridges(aes(fill=group), scale=1.5, rel_min_height=0.001, alpha=0.45) +
  scale_fill_manual(values=c("F"="#d73027","M"="#4575b4")) +
  scale_y_discrete(labels=c("F"="Female","M"="Male")) +
  xlim(rs_range_rna) + ridge_theme + labs(x="Risk score", y="Gender") +
  annotate("text", x=rs_range_rna[2]-0.1, y=Inf, label=p_C, hjust=1, vjust=1.5, size=3.5)

toPlot_D <- data.frame(risk_score=dat_plot$risk_lp, group=as.character(dat_plot$MYC)) %>% filter(!is.na(group))
p_D <- fmt_p(t.test(risk_score ~ group, data=toPlot_D)$p.value)
D <- ggplot(toPlot_D, aes(risk_score, group)) +
  geom_density_ridges(aes(fill=group), scale=1.5, rel_min_height=0.001, alpha=0.45) +
  scale_fill_manual(values=c("0"="grey80","1"="#e31a1c")) +
  scale_y_discrete(labels=c("0"="Not amplified","1"="Amplified")) +
  xlim(rs_range_rna) + ridge_theme + labs(x="Risk score", y="MYC") +
  annotate("text", x=rs_range_rna[2]-0.1, y=Inf, label=p_D, hjust=1, vjust=1.5, size=3.5)

toPlot_E <- data.frame(risk_score=dat_plot$risk_lp, group=as.character(dat_plot$MYCN)) %>% filter(!is.na(group))
p_E <- fmt_p(t.test(risk_score ~ group, data=toPlot_E)$p.value)
E_plot <- ggplot(toPlot_E, aes(risk_score, group)) +
  geom_density_ridges(aes(fill=group), scale=1.5, rel_min_height=0.001, alpha=0.45) +
  scale_fill_manual(values=c("0"="grey80","1"="#ff7f00")) +
  scale_y_discrete(labels=c("0"="Not amplified","1"="Amplified")) +
  xlim(rs_range_rna) + ridge_theme + labs(x="Risk score", y="MYCN") +
  annotate("text", x=rs_range_rna[2]-0.1, y=Inf, label=p_E, hjust=1, vjust=1.5, size=3.5)

toPlot_F <- data.frame(risk_score=dat_plot$risk_lp, group=dat_plot$AgeGroup) %>%
  filter(group %in% c("0-3","3-16")) %>% mutate(group=factor(group, levels=c("0-3","3-16")))
p_F <- fmt_p(t.test(risk_score ~ group, data=toPlot_F)$p.value)
F_plot <- ggplot(toPlot_F, aes(risk_score, group)) +
  geom_density_ridges(aes(fill=group), scale=1.5, rel_min_height=0.001, alpha=0.45) +
  scale_fill_manual(values=c("0-3"="#ffffcc","3-16"="#41b6c4")) +
  xlim(rs_range_rna) + ridge_theme + labs(x="Risk score", y="Age group") +
  annotate("text", x=rs_range_rna[2]-0.1, y=Inf, label=p_F, hjust=1, vjust=1.5, size=3.5)

# Export PDF
pdf("RNAseq_G3G4_validation_oncoprint.pdf", width=14, height=8)
draw(ht_rnaseq,
     heatmap_legend_side    = "right",
     annotation_legend_side = "right",
     annotation_legend_list = list(lgd_tertile_rna, lgd_subgroup_rna, lgd_gender_rna,
                                   lgd_age_rna, lgd_mstage_rna, lgd_myc_rna, lgd_mycn_rna),
     padding         = unit(c(2,30,2,2),"mm"),
     column_title    = "Group 3/4 medulloblastoma: RNA-seq validation of transcriptomic risk score",
     column_title_gp = gpar(fontsize=13, fontface="bold"))
decorate_annotation("OS", {
  grid.lines(x=unit(c(0,1),"npc"), y=unit(c(5,5),"native"),
             gp=gpar(col="black", lty="dotted", lwd=1))
})
grid.newpage()
print(ggarrange(A, B, C, D, E_plot, F_plot, ncol=2, nrow=3, labels=c("A","B","C","D","E","F")))
dev.off()
cat("Done \u2014 RNAseq_G3G4_validation_oncoprint.pdf written.\n")

########################################################################################################################
#
#   SECTION 10.8: SHH RNA-SEQ VALIDATION
#   NOTE: No OS data in E-MTAB-10767 for SHH — KM validation not possible.
#   NOTE: MYC excluded — no amplified cases in SHH cohort (n = 0/66).
#   Output: Clinical annotation oncoprint only
#
########################################################################################################################

# -----------------------------------------------------------------------------
# 10.8 Project SHH microarray signature onto RNA-seq
# -----------------------------------------------------------------------------

shh_samples_rna   <- sdrf %>% filter(`Characteristics[subgroup]` == "SHH") %>% pull(`Source Name`) %>% as.character()
shh_in_rna        <- intersect(shh_samples_rna, colnames(expr_mat_rna))
x_rna_shh         <- t(expr_mat_rna[, shh_in_rna, drop = FALSE])

overlap_shh       <- intersect(names(final_coefs_shh), colnames(x_rna_shh))
rs_rna_shh        <- as.numeric(x_rna_shh[, overlap_shh] %*% final_coefs_shh[overlap_shh])
names(rs_rna_shh) <- rownames(x_rna_shh)

cat("SHH RNA-seq — genes in model:", length(names(final_coefs_shh)),
    "| found in RNA-seq:", length(overlap_shh), "\n")

cuts_rna_shh     <- quantile(rs_rna_shh, probs = c(1/3, 2/3))
risk_out_shh_rna <- data.frame(SampleID = names(rs_rna_shh), risk_lp = rs_rna_shh) %>%
  mutate(tertile = cut(risk_lp, breaks = c(-Inf, cuts_rna_shh, Inf),
                       labels = c("Low", "Intermediate", "High")))

cat("SHH RNA-seq tertile sizes:\n"); print(table(risk_out_shh_rna$tertile))

# -----------------------------------------------------------------------------
# 10.9 SHH clinical metadata from SDRF
# -----------------------------------------------------------------------------

clinical_shh_rna <- sdrf %>%
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
      ifelse(grepl("not available|na|unknown", x),                         "Unknown",
             ifelse(grepl("over 16|>16|> 16|adult", x),                          ">16",
                    ifelse(grepl("0 to 3|0-3|infant|toddler", x),                       "0-3",
                           ifelse(grepl("3 to 16|3-16|child|pediatric|paediatric", x),         "3-16", "Unknown"))))
    }),
    Mstage = ifelse(tolower(trimws(Mstage_raw)) %in% c("true","1","yes","y"), "1",
                    ifelse(tolower(trimws(Mstage_raw)) %in% c("false","0","no","n"), "0", "Unknown")),
    MYCN = ifelse(tolower(trimws(MYCN_raw)) %in% c("true","1","yes","y"), "1", "0")
  )

# -----------------------------------------------------------------------------
# 10.10 SHH oncoprint — merge, clean, build heatmap
# -----------------------------------------------------------------------------

dat_plot_shh <- risk_out_shh_rna %>%
  left_join(clinical_shh_rna, by = "SampleID") %>%
  arrange(risk_lp)

dat_plot_shh$Gender  <- ifelse(is.na(dat_plot_shh$Gender)  | dat_plot_shh$Gender  == "", "Unknown", dat_plot_shh$Gender)
dat_plot_shh$Mstage  <- ifelse(is.na(dat_plot_shh$Mstage)  | dat_plot_shh$Mstage  == "", "Unknown", dat_plot_shh$Mstage)
dat_plot_shh$MYCN    <- as.character(dat_plot_shh$MYCN)
dat_plot_shh$tertile <- ifelse(is.na(dat_plot_shh$tertile) | dat_plot_shh$tertile == "", "Unknown",
                               as.character(dat_plot_shh$tertile))

dat_plot_shh$AgeGroup <- local({
  x <- tolower(trimws(as.character(dat_plot_shh$Age)))
  ifelse(grepl("not available|na|unknown", x),                         "Unknown",
         ifelse(grepl("over 16|>16|> 16|adult", x),                          ">16",
                ifelse(grepl("0 to 3|0-3|infant|toddler", x),                       "0-3",
                       ifelse(grepl("3 to 16|3-16|child|pediatric|paediatric", x),         "3-16", "Unknown"))))
})

dat_plot_shh <- dat_plot_shh %>% filter(!AgeGroup %in% c(">16", "Unknown"))
cat("SHH RNA-seq oncoprint samples:", nrow(dat_plot_shh), "\n")

# Colour palettes
tertile_col_shh_rna <- c("Low"="#2166ac","Intermediate"="#f4a582","High"="#d6604d","Unknown"="lightgrey")
gender_col_shh_rna  <- c("M"="white","F"="black","Unknown"="lightgrey")
mstage_col_shh_rna  <- c("0"="white","1"="black","Unknown"="lightgrey")
mycn_col_shh_rna    <- c("0"="white","1"="#ff7f00")
age_col_shh_rna     <- c("0-3"="#ffffcc","3-16"="#41b6c4")

cat_anno_shh_rna <- HeatmapAnnotation(
  `Risk tertile` = anno_simple(dat_plot_shh$tertile,   col = tertile_col_shh_rna, border = TRUE),
  Gender         = anno_simple(dat_plot_shh$Gender,    col = gender_col_shh_rna,  border = TRUE),
  `Age group`    = anno_simple(dat_plot_shh$AgeGroup,  col = age_col_shh_rna,     border = TRUE),
  `M-stage`      = anno_simple(dat_plot_shh$Mstage,    col = mstage_col_shh_rna,  border = TRUE),
  MYCN           = anno_simple(dat_plot_shh$MYCN,      col = mycn_col_shh_rna,    border = TRUE),
  annotation_name_side = "left")

rs_mat_shh <- matrix(dat_plot_shh$risk_lp, nrow = 1)
colnames(rs_mat_shh) <- dat_plot_shh$SampleID

ht_shh_rna <- Heatmap(rs_mat_shh, name = "Risk score",
                      col = colorRamp2(
                        c(min(dat_plot_shh$risk_lp, na.rm=TRUE), median(dat_plot_shh$risk_lp, na.rm=TRUE), max(dat_plot_shh$risk_lp, na.rm=TRUE)),
                        c("#2166ac","white","#d6604d")),
                      border = TRUE, show_column_names = FALSE, cluster_columns = FALSE, cluster_rows = FALSE,
                      row_names_side = "left", height = unit(0.6,"cm"), show_heatmap_legend = TRUE,
                      top_annotation = cat_anno_shh_rna)

lgd_tertile_shh_rna <- Legend(labels=c("Low","Intermediate","High","Unknown"), title="Risk tertile",
                              legend_gp=gpar(fill=c("#2166ac","#f4a582","#d6604d","lightgrey")), border=TRUE)
lgd_gender_shh_rna  <- Legend(labels=c("Male","Female","Unknown"), title="Gender",
                              legend_gp=gpar(fill=c("white","black","lightgrey")), border=TRUE)
lgd_age_shh_rna     <- Legend(labels=c("0-3","3-16"), title="Age group",
                              legend_gp=gpar(fill=c("#ffffcc","#41b6c4")), border=TRUE)
lgd_mstage_shh_rna  <- Legend(labels=c("M0","M+","Unknown"), title="M-stage",
                              legend_gp=gpar(fill=c("white","black","lightgrey")), border=TRUE)
lgd_mycn_shh_rna    <- Legend(labels=c("Not amplified","Amplified"), title="MYCN",
                              legend_gp=gpar(fill=c("white","#ff7f00")), border=TRUE)

# Ridge plots
rs_range_shh_rna <- range(dat_plot_shh$risk_lp, na.rm = TRUE)

toPlot_A_shh <- data.frame(risk_score=dat_plot_shh$risk_lp, group=dat_plot_shh$Mstage) %>% filter(group != "Unknown")
p_A_shh <- fmt_p(t.test(risk_score ~ group, data=toPlot_A_shh)$p.value)
A_shh <- ggplot(toPlot_A_shh, aes(risk_score, group)) +
  geom_density_ridges(aes(fill=group), scale=1.5, rel_min_height=0.001, alpha=0.45) +
  scale_fill_manual(values=c("0"="grey80","1"="black")) +
  scale_y_discrete(labels=c("0"="M0","1"="M+")) +
  xlim(rs_range_shh_rna) + ridge_theme + labs(x="Risk score", y="M-stage") +
  annotate("text", x=rs_range_shh_rna[2]-0.1, y=Inf, label=p_A_shh, hjust=1, vjust=1.5, size=3.5)

toPlot_B_shh <- data.frame(risk_score=dat_plot_shh$risk_lp, group=dat_plot_shh$Gender) %>% filter(group != "Unknown")
p_B_shh <- fmt_p(t.test(risk_score ~ group, data=toPlot_B_shh)$p.value)
B_shh <- ggplot(toPlot_B_shh, aes(risk_score, group)) +
  geom_density_ridges(aes(fill=group), scale=1.5, rel_min_height=0.001, alpha=0.45) +
  scale_fill_manual(values=c("F"="#d73027","M"="#4575b4")) +
  scale_y_discrete(labels=c("F"="Female","M"="Male")) +
  xlim(rs_range_shh_rna) + ridge_theme + labs(x="Risk score", y="Gender") +
  annotate("text", x=rs_range_shh_rna[2]-0.1, y=Inf, label=p_B_shh, hjust=1, vjust=1.5, size=3.5)

toPlot_C_shh <- data.frame(risk_score=dat_plot_shh$risk_lp, group=as.character(dat_plot_shh$MYCN)) %>% filter(!is.na(group))
p_C_shh <- fmt_p(t.test(risk_score ~ group, data=toPlot_C_shh)$p.value)
C_shh <- ggplot(toPlot_C_shh, aes(risk_score, group)) +
  geom_density_ridges(aes(fill=group), scale=1.5, rel_min_height=0.001, alpha=0.45) +
  scale_fill_manual(values=c("0"="grey80","1"="#ff7f00")) +
  scale_y_discrete(labels=c("0"="Not amplified","1"="Amplified")) +
  xlim(rs_range_shh_rna) + ridge_theme + labs(x="Risk score", y="MYCN") +
  annotate("text", x=rs_range_shh_rna[2]-0.1, y=Inf, label=p_C_shh, hjust=1, vjust=1.5, size=3.5)

toPlot_D_shh <- data.frame(risk_score=dat_plot_shh$risk_lp, group=dat_plot_shh$AgeGroup) %>%
  filter(group %in% c("0-3","3-16")) %>% mutate(group=factor(group, levels=c("0-3","3-16")))
p_D_shh <- fmt_p(t.test(risk_score ~ group, data=toPlot_D_shh)$p.value)
D_shh <- ggplot(toPlot_D_shh, aes(risk_score, group)) +
  geom_density_ridges(aes(fill=group), scale=1.5, rel_min_height=0.001, alpha=0.45) +
  scale_fill_manual(values=c("0-3"="#ffffcc","3-16"="#41b6c4")) +
  xlim(rs_range_shh_rna) + ridge_theme + labs(x="Risk score", y="Age group") +
  annotate("text", x=rs_range_shh_rna[2]-0.1, y=Inf, label=p_D_shh, hjust=1, vjust=1.5, size=3.5)

# Export PDF
pdf("RNAseq_SHH_validation_oncoprint.pdf", width=14, height=8)
draw(ht_shh_rna,
     heatmap_legend_side    = "right",
     annotation_legend_side = "right",
     annotation_legend_list = list(lgd_tertile_shh_rna, lgd_gender_shh_rna,
                                   lgd_age_shh_rna, lgd_mstage_shh_rna, lgd_mycn_shh_rna),
     padding         = unit(c(2,30,2,2),"mm"),
     column_title    = "SHH medulloblastoma: RNA-seq clinical annotation (no survival data available)",
     column_title_gp = gpar(fontsize=13, fontface="bold"))
grid.newpage()
print(ggarrange(A_shh, B_shh, C_shh, D_shh, ncol=2, nrow=2, labels=c("A","B","C","D")))
dev.off()
cat("Done \u2014 RNAseq_SHH_validation_oncoprint.pdf written.\n")


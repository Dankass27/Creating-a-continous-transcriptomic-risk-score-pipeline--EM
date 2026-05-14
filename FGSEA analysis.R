################################################################################
# Part of my original source code September 2025 - May 2026
# MEDULLOBLASTOMA — GENE SET ENRICHMENT ANALYSIS (fgsea)
# Dataset: GSE85217 (Taylor Lab, 763 pediatric MB samples)
#
# Student: Ewan McGibbon
# Date:    7/4/2026
#
# Description:
#   fgsea analysis for Group 3/4 and SHH medulloblastoma subgroups.
#   Genes ranked by Pearson correlation with transcriptomic risk score.
#
# REQUIRES: MB_main.R to have been run first.
# The following objects must be present in the workspace:
#   - dat_glm          : G3/4 clinical data with risk_score column
#   - dat_glm_shh      : SHH clinical data with risk_score column
#   - expr_noXY        : full expression matrix (genes x samples, sex chrs removed)
#   - expr_shh_cc      : SHH expression matrix (genes x samples)
#
################################################################################

library(fgsea)
library(biomaRt)
library(ggplot2)
library(dplyr)
library(patchwork)

# ── Load MSigDB Gene Sets ──────────────────────────────────────────────────────
# All 9 gene sets loaded
# Unhash and swap into fgsea() to test against a different collection

pathways_hallmark <- gmtPathways("h.all.v2024.1.Hs.symbols.gmt")     # Hallmark (50 pathways)
pathways_c1       <- gmtPathways("c1.all.v2024.1.Hs.symbols.gmt")    # Chromosomal position
pathways_c2       <- gmtPathways("c2.all.v2024.1.Hs.symbols.gmt")    # Curated (KEGG, Reactome)
pathways_c6       <- gmtPathways("c6.all.v2024.1.Hs.symbols.gmt")    # Oncogenic signatures
# pathways_c3     <- gmtPathways("c3.all.v2024.1.Hs.symbols.gmt")    # Regulatory targets
# pathways_c4     <- gmtPathways("c4.all.v2024.1.Hs.symbols.gmt")    # Computational cancer
# pathways_c5     <- gmtPathways("c5.all.v2024.1.Hs.symbols.gmt")    # Gene Ontology
# pathways_c7     <- gmtPathways("c7.all.v2024.1.Hs.symbols.gmt")    # Immunological
# pathways_c8     <- gmtPathways("c8.all.v2024.1.Hs.symbols.gmt")    # Cell type signatures

# ── BioMart connection ─────────────────────────────────────────────────────────
mart <- tryCatch(
  useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast"),
  error = function(e) tryCatch(
    useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "uswest"),
    error = function(e) tryCatch(
      useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "asia"),
      error = function(e) useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
    )
  )
)


########################################################################################################################
########################################################################################################################
#
#                                   1. fgsea Analysis — G3/4 Medulloblastoma
#
########################################################################################################################
########################################################################################################################
#
# Pre-processing: Build ranked gene list from expression vs risk score
# Analysis:       Run fgsea against MSigDB gene sets
#
########################################################################################################################

# ── Pre-processing ─────────────────────────────────────────────────────────────

# Use full transcriptome (sex chromosomes removed); subset to G3/4 samples
expr_mat_g34 <- expr_noXY[, colnames(expr_noXY) %in% dat_glm$Study_ID]

anno <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = rownames(expr_mat_g34),
  mart       = mart
)
anno         <- anno[!duplicated(anno$ensembl_gene_id) & anno$hgnc_symbol != "", ]
idx          <- match(rownames(expr_mat_g34), anno$ensembl_gene_id)
gene_symbols <- ifelse(is.na(idx), NA, anno$hgnc_symbol[idx])

common_pts   <- intersect(colnames(expr_mat_g34), dat_glm$Study_ID)
cat("G3/4 patients in common:", length(common_pts), "\n")

expr_matched <- apply(expr_mat_g34[, common_pts], 2, as.numeric)
risk_matched <- dat_glm$risk_score[match(common_pts, dat_glm$Study_ID)]

cors         <- apply(expr_matched, 1, function(x) {
  cor(x, risk_matched, method = "pearson", use = "complete.obs")
})
names(cors) <- gene_symbols
cors        <- cors[!is.na(cors) & !is.na(names(cors)) & names(cors) != ""]
cors        <- cors[!duplicated(names(cors))]
ranks_fgsea <- sort(cors, decreasing = TRUE)

# ── fgsea runs ─────────────────────────────────────────────────────────────────

##  ↓↓ SWAP HERE ↓↓
# active_pathways <- pathways_c2
# active_label    <- "C2"
##  ↑↑ SWAP HERE ↑↑

## Run fgsea
set.seed(42)
fgsea_res <- fgsea(
  pathways = active_pathways,
  stats    = ranks_fgsea,
  minSize  = 15,
  maxSize  = 500
)

# Order by adjusted p-value
fgsea_res <- fgsea_res[order(fgsea_res$padj), ]
head(fgsea_res, 20)

# Significant pathways (FDR < 0.05)
sig_res <- fgsea_res[fgsea_res$padj < 0.05, ]
cat("Significant", active_label, "pathways:", nrow(sig_res), "\n")

# ── Plot 1: Summary table ──────────────────────────────────────────────────────
top_up    <- sig_res[sig_res$NES > 0, ]$pathway[1:10]
top_down  <- sig_res[sig_res$NES < 0, ]$pathway[1:10]
top_paths <- c(top_up, rev(top_down))

plotGseaTable(
  active_pathways[top_paths],
  stats     = ranks_fgsea,
  fgseaRes  = fgsea_res,
  gseaParam = 0.5
)

# ── Bubble plot runs ───────────────────────────────────────────────────────────
set.seed(42); fgsea_res_g34_h  <- fgsea(pathways = pathways_hallmark, stats = ranks_fgsea, minSize = 15, maxSize = 500)
set.seed(42); fgsea_res_g34_c1 <- fgsea(pathways = pathways_c1,       stats = ranks_fgsea, minSize = 15, maxSize = 500)
set.seed(42); fgsea_res_g34_c2 <- fgsea(pathways = pathways_c2,       stats = ranks_fgsea, minSize = 15, maxSize = 500)
set.seed(42); fgsea_res_g34_c6 <- fgsea(pathways = pathways_c6,       stats = ranks_fgsea, minSize = 15, maxSize = 500)


########################################################################################################################
########################################################################################################################
#
#                                   2. fgsea Analysis — SHH Medulloblastoma
#
########################################################################################################################
########################################################################################################################
#
# Pre-processing: Build ranked gene list from expression vs risk score
# Analysis:       Run fgsea against MSigDB gene sets
#
########################################################################################################################

# ── Pre-processing ─────────────────────────────────────────────────────────────

expr_mat_shh <- expr_shh_cc

anno_shh <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = rownames(expr_mat_shh),
  mart       = mart
)
anno_shh         <- anno_shh[!duplicated(anno_shh$ensembl_gene_id) & anno_shh$hgnc_symbol != "", ]
idx_shh          <- match(rownames(expr_mat_shh), anno_shh$ensembl_gene_id)
gene_symbols_shh <- ifelse(is.na(idx_shh), NA, anno_shh$hgnc_symbol[idx_shh])

common_pts_shh   <- intersect(colnames(expr_mat_shh), dat_glm_shh$Study_ID)
cat("SHH patients in common:", length(common_pts_shh), "\n")

expr_matched_shh <- apply(expr_mat_shh[, common_pts_shh], 2, as.numeric)
risk_matched_shh <- dat_glm_shh$risk_score[match(common_pts_shh, dat_glm_shh$Study_ID)]

cors_shh         <- apply(expr_matched_shh, 1, function(x) {
  cor(x, risk_matched_shh, method = "pearson", use = "complete.obs")
})
names(cors_shh) <- gene_symbols_shh
cors_shh        <- cors_shh[!is.na(cors_shh) & !is.na(names(cors_shh)) & names(cors_shh) != ""]
cors_shh        <- cors_shh[!duplicated(names(cors_shh))]
ranks_fgsea_shh <- sort(cors_shh, decreasing = TRUE)

# ── fgsea runs ─────────────────────────────────────────────────────────────────

# ↓↓ SWAP HERE ↓↓
active_pathways_shh <- pathways_hallmark
active_label_shh    <- "hallmark"
# ↑↑ SWAP HERE ↑↑

set.seed(42)
fgsea_res_shh <- fgsea(
  pathways = active_pathways_shh,
  stats    = ranks_fgsea_shh,
  minSize  = 15,
  maxSize  = 500
)

# Order by adjusted p-value
fgsea_res_shh <- fgsea_res_shh[order(fgsea_res_shh$padj), ]
head(fgsea_res_shh, 20)

# Significant pathways (FDR < 0.05)
sig_res_shh <- fgsea_res_shh[fgsea_res_shh$padj < 0.05, ]
cat("Significant", active_label_shh, "pathways:", nrow(sig_res_shh), "\n")

# ── Plot: Summary table ────────────────────────────────────────────────────────
top_up_shh    <- sig_res_shh[sig_res_shh$NES > 0, ]$pathway[1:10]
top_down_shh  <- sig_res_shh[sig_res_shh$NES < 0, ]$pathway[1:10]
top_paths_shh <- c(top_up_shh, rev(top_down_shh))

pdf("SHH_fgsea_hallmark_table.pdf", width = 20, height = 8)
par(cex = 0.6)
plotGseaTable(
  active_pathways_shh[top_paths_shh],
  stats     = ranks_fgsea_shh,
  fgseaRes  = fgsea_res_shh,
  gseaParam = 0.5
)
par(cex = 1)
dev.off()

# ── Bubble plot runs ───────────────────────────────────────────────────────────
set.seed(42); fgsea_res_shh_h  <- fgsea(pathways = pathways_hallmark, stats = ranks_fgsea_shh, minSize = 15, maxSize = 500)
set.seed(42); fgsea_res_shh_c1 <- fgsea(pathways = pathways_c1,       stats = ranks_fgsea_shh, minSize = 15, maxSize = 500)
set.seed(42); fgsea_res_shh_c2 <- fgsea(pathways = pathways_c2,       stats = ranks_fgsea_shh, minSize = 15, maxSize = 500)
set.seed(42); fgsea_res_shh_c6 <- fgsea(pathways = pathways_c6,       stats = ranks_fgsea_shh, minSize = 15, maxSize = 500)


########################################################################################################################
# Bubble plots — G3/4 and SHH
########################################################################################################################

# ── Helper ─────────────────────────────────────────────────────────────────────
extract_top <- function(fgsea_res, collection_label, n_top = 10) {
  fgsea_res %>%
    as.data.frame() %>%
    filter(padj < 0.05) %>%
    arrange(padj) %>%
    head(n_top) %>%
    mutate(
      Collection = collection_label,
      pathway    = gsub("HALLMARK_|REACTOME_|KEGG_|GOBP_|GOCC_|GOMF_", "", pathway),
      pathway    = gsub("_", " ", pathway),
      pathway    = stringr::str_to_title(pathway)
    )
}

# ── G3/4 bubble plot ──────────────────────────────────────────────────────────
g34_bubble_df <- bind_rows(
  extract_top(fgsea_res_g34_h,  "Hallmark",       10),
  extract_top(fgsea_res_g34_c1, "C1 Chromosomal",  5),
  extract_top(fgsea_res_g34_c2, "C2 Curated",     10),
  extract_top(fgsea_res_g34_c6, "C6 Oncogenic",    5)
) %>%
  mutate(
    pathway      = factor(pathway, levels = rev(unique(pathway))),
    neglog10padj = -log10(padj),
    Direction    = ifelse(NES > 0, "Enriched", "Depleted")
  )

p_g34_bubble <- ggplot(g34_bubble_df, aes(
  x      = NES,
  y      = pathway,
  size   = size,
  colour = neglog10padj
)) +
  geom_point(alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  scale_colour_gradient(low = "#fee8c8", high = "#b30000", name = expression(-log[10](p[adj]))) +
  scale_size_continuous(range = c(3, 10), name = "Gene set size") +
  facet_wrap(~ Collection, scales = "free_y", ncol = 1) +
  labs(
    title    = "Group 3/4 - fgsea Enrichment Summary",
    subtitle = "Top significant pathways per MSigDB collection (FDR < 0.05)",
    x        = "Normalised Enrichment Score (NES)",
    y        = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "#2c3e50", colour = NA),
    strip.text       = element_text(colour = "white", face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text.y      = element_text(size = 8)
  )

# ── SHH bubble plot ───────────────────────────────────────────────────────────
shh_bubble_df <- bind_rows(
  extract_top(fgsea_res_shh_h,  "Hallmark",       10),
  extract_top(fgsea_res_shh_c1, "C1 Chromosomal",  5),
  extract_top(fgsea_res_shh_c2, "C2 Curated",     10),
  extract_top(fgsea_res_shh_c6, "C6 Oncogenic",    5)
) %>%
  mutate(
    pathway      = factor(pathway, levels = rev(unique(pathway))),
    neglog10padj = -log10(padj),
    Direction    = ifelse(NES > 0, "Enriched", "Depleted")
  )

p_shh_bubble <- ggplot(shh_bubble_df, aes(
  x      = NES,
  y      = pathway,
  size   = size,
  colour = neglog10padj
)) +
  geom_point(alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  scale_colour_gradient(low = "#e8f4f8", high = "#08519c", name = expression(-log[10](p[adj]))) +
  scale_size_continuous(range = c(3, 10), name = "Gene set size") +
  facet_wrap(~ Collection, scales = "free_y", ncol = 1) +
  labs(
    title    = "SHH - fgsea Enrichment Summary",
    subtitle = "Top significant pathways per MSigDB collection (FDR < 0.05)",
    x        = "Normalised Enrichment Score (NES)",
    y        = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "#2c3e50", colour = NA),
    strip.text       = element_text(colour = "white", face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text.y      = element_text(size = 8)
  )


########################################################################################################################
# Combined figures — bubble plot + individual enrichment plots
########################################################################################################################

# ── G3/4 enrichment plots ──────────────────────────────────────────────────────

# MYC targets V1 — top hit, canonical G3/4 high-risk biology
p_myc_g34 <- plotEnrichment(pathways_hallmark[["HALLMARK_MYC_TARGETS_V1"]], ranks_fgsea) +
  labs(title = "HALLMARK_MYC_TARGETS_V1", subtitle = "Group 3/4") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 9), plot.subtitle = element_text(size = 8))

# E2F targets — cell cycle, downstream of MYC
p_e2f_g34 <- plotEnrichment(pathways_hallmark[["HALLMARK_E2F_TARGETS"]], ranks_fgsea) +
  labs(title = "HALLMARK_E2F_TARGETS", subtitle = "Group 3/4") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 9), plot.subtitle = element_text(size = 8))

# Chr8q24 — MYC locus, chromosomal amplification
p_chr8q24 <- plotEnrichment(pathways_c1[["chr8q24"]], ranks_fgsea) +
  labs(title = "chr8q24 (MYC locus)", subtitle = "Group 3/4") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 9), plot.subtitle = element_text(size = 8))

p_g34_combined <- p_g34_bubble / (p_myc_g34 | p_e2f_g34 | p_chr8q24) +
  plot_layout(heights = c(3, 1)) +
  plot_annotation(
    title = "Group 3/4 Medulloblastoma - fgsea Analysis",
    theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
  )

ggsave("G34_fgsea_combined.pdf", p_g34_combined, width = 12, height = 20)

# ── SHH enrichment plots ───────────────────────────────────────────────────────

# mTORC1 — top hit, SHH pathway activates mTORC1 downstream of SMO/GLI
p_mtorc1_shh <- plotEnrichment(pathways_hallmark[["HALLMARK_MTORC1_SIGNALING"]], ranks_fgsea_shh) +
  labs(title = "HALLMARK_MTORC1_SIGNALING", subtitle = "SHH") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 9), plot.subtitle = element_text(size = 8))

# MYC targets V1 — proliferative high-risk signature
p_myc_shh <- plotEnrichment(pathways_hallmark[["HALLMARK_MYC_TARGETS_V1"]], ranks_fgsea_shh) +
  labs(title = "HALLMARK_MYC_TARGETS_V1", subtitle = "SHH") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 9), plot.subtitle = element_text(size = 8))

# Chr17p13 — TP53 locus, negatively enriched; consistent with TP53 mutation in oncoprint
p_chr17p13 <- plotEnrichment(pathways_c1[["chr17p13"]], ranks_fgsea_shh) +
  labs(title = "chr17p13 (TP53 locus)", subtitle = "SHH") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 9), plot.subtitle = element_text(size = 8))

p_shh_combined <- p_shh_bubble / (p_mtorc1_shh | p_myc_shh | p_chr17p13) +
  plot_layout(heights = c(3, 1)) +
  plot_annotation(
    title = "SHH Medulloblastoma - fgsea Analysis",
    theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
  )

ggsave("SHH_fgsea_combined.pdf", p_shh_combined, width = 12, height = 20)

# ── Appendix: fgsea summary tables — G3/4 and SHH ────────────────────────────
pdf("Appendix_fgsea_tables.pdf", width = 20, height = 8)

# G3/4 — Hallmark
par(cex = 0.6)
sig_h <- fgsea_res_g34_h[fgsea_res_g34_h$padj < 0.05, ]
top_paths_h <- c(sig_h[sig_h$NES > 0, ]$pathway[1:10], rev(sig_h[sig_h$NES < 0, ]$pathway[1:10]))
plotGseaTable(pathways_hallmark[top_paths_h], stats = ranks_fgsea, fgseaRes = fgsea_res_g34_h, gseaParam = 0.5)
title(main = "Group 3/4 - Hallmark", outer = TRUE, line = -1)

# G3/4 — C2 Curated
sig_c2 <- fgsea_res_g34_c2[fgsea_res_g34_c2$padj < 0.05, ]
top_paths_c2 <- c(sig_c2[sig_c2$NES > 0, ]$pathway[1:10], rev(sig_c2[sig_c2$NES < 0, ]$pathway[1:10]))
plotGseaTable(pathways_c2[top_paths_c2], stats = ranks_fgsea, fgseaRes = fgsea_res_g34_c2, gseaParam = 0.5)
title(main = "Group 3/4 - C2 Curated", outer = TRUE, line = -1)

# G3/4 — C6 Oncogenic
sig_c6 <- fgsea_res_g34_c6[fgsea_res_g34_c6$padj < 0.05, ]
top_paths_c6 <- c(sig_c6[sig_c6$NES > 0, ]$pathway[1:10], rev(sig_c6[sig_c6$NES < 0, ]$pathway[1:10]))
plotGseaTable(pathways_c6[top_paths_c6], stats = ranks_fgsea, fgseaRes = fgsea_res_g34_c6, gseaParam = 0.5)
title(main = "Group 3/4 - C6 Oncogenic", outer = TRUE, line = -1)

# SHH — Hallmark
sig_shh_h <- fgsea_res_shh_h[fgsea_res_shh_h$padj < 0.05, ]
top_paths_shh_h <- c(sig_shh_h[sig_shh_h$NES > 0, ]$pathway[1:10], rev(sig_shh_h[sig_shh_h$NES < 0, ]$pathway[1:10]))
plotGseaTable(pathways_hallmark[top_paths_shh_h], stats = ranks_fgsea_shh, fgseaRes = fgsea_res_shh_h, gseaParam = 0.5)
title(main = "SHH - Hallmark", outer = TRUE, line = -1)

# SHH — C2 Curated
sig_shh_c2 <- fgsea_res_shh_c2[fgsea_res_shh_c2$padj < 0.05, ]
top_paths_shh_c2 <- c(sig_shh_c2[sig_shh_c2$NES > 0, ]$pathway[1:10], rev(sig_shh_c2[sig_shh_c2$NES < 0, ]$pathway[1:10]))
plotGseaTable(pathways_c2[top_paths_shh_c2], stats = ranks_fgsea_shh, fgseaRes = fgsea_res_shh_c2, gseaParam = 0.5)
title(main = "SHH - C2 Curated", outer = TRUE, line = -1)

# SHH — C6 Oncogenic
sig_shh_c6 <- fgsea_res_shh_c6[fgsea_res_shh_c6$padj < 0.05, ]
top_paths_shh_c6 <- c(sig_shh_c6[sig_shh_c6$NES > 0, ]$pathway[1:10], rev(sig_shh_c6[sig_shh_c6$NES < 0, ]$pathway[1:10]))
plotGseaTable(pathways_c6[top_paths_shh_c6], stats = ranks_fgsea_shh, fgseaRes = fgsea_res_shh_c6, gseaParam = 0.5)
title(main = "SHH - C6 Oncogenic", outer = TRUE, line = -1)

par(cex = 1)
dev.off()

# ── Appendix: fgsea summary tables — G3/4 and SHH ────────────────────────────
pdf("Appendix_fgsea_tables.pdf", width = 20, height = 8)
par(cex = 0.6)

# G3/4 — Hallmark
sig_h <- fgsea_res_g34_h[fgsea_res_g34_h$padj < 0.05, ]
top_paths_h <- c(sig_h[sig_h$NES > 0, ]$pathway[1:10], rev(sig_h[sig_h$NES < 0, ]$pathway[1:10]))
plotGseaTable(pathways_hallmark[top_paths_h], stats = ranks_fgsea, fgseaRes = fgsea_res_g34_h, gseaParam = 0.5)

# G3/4 — C2 Curated
sig_c2 <- fgsea_res_g34_c2[fgsea_res_g34_c2$padj < 0.05, ]
top_paths_c2 <- c(sig_c2[sig_c2$NES > 0, ]$pathway[1:10], rev(sig_c2[sig_c2$NES < 0, ]$pathway[1:10]))
plotGseaTable(pathways_c2[top_paths_c2], stats = ranks_fgsea, fgseaRes = fgsea_res_g34_c2, gseaParam = 0.5)

# G3/4 — C6 Oncogenic
sig_c6 <- fgsea_res_g34_c6[fgsea_res_g34_c6$padj < 0.05, ]
top_paths_c6 <- c(sig_c6[sig_c6$NES > 0, ]$pathway[1:10], rev(sig_c6[sig_c6$NES < 0, ]$pathway[1:10]))
plotGseaTable(pathways_c6[top_paths_c6], stats = ranks_fgsea, fgseaRes = fgsea_res_g34_c6, gseaParam = 0.5)

# SHH — Hallmark
sig_shh_h <- fgsea_res_shh_h[fgsea_res_shh_h$padj < 0.05, ]
top_paths_shh_h <- c(sig_shh_h[sig_shh_h$NES > 0, ]$pathway[1:10], rev(sig_shh_h[sig_shh_h$NES < 0, ]$pathway[1:10]))
plotGseaTable(pathways_hallmark[top_paths_shh_h], stats = ranks_fgsea_shh, fgseaRes = fgsea_res_shh_h, gseaParam = 0.5)

# SHH — C2 Curated
sig_shh_c2 <- fgsea_res_shh_c2[fgsea_res_shh_c2$padj < 0.05, ]
top_paths_shh_c2 <- c(sig_shh_c2[sig_shh_c2$NES > 0, ]$pathway[1:10], rev(sig_shh_c2[sig_shh_c2$NES < 0, ]$pathway[1:10]))
plotGseaTable(pathways_c2[top_paths_shh_c2], stats = ranks_fgsea_shh, fgseaRes = fgsea_res_shh_c2, gseaParam = 0.5)

# SHH — C6 Oncogenic
sig_shh_c6 <- fgsea_res_shh_c6[fgsea_res_shh_c6$padj < 0.05, ]
top_paths_shh_c6 <- c(sig_shh_c6[sig_shh_c6$NES > 0, ]$pathway[1:10], rev(sig_shh_c6[sig_shh_c6$NES < 0, ]$pathway[1:10]))
plotGseaTable(pathways_c6[top_paths_shh_c6], stats = ranks_fgsea_shh, fgseaRes = fgsea_res_shh_c6, gseaParam = 0.5)

par(cex = 1)
dev.off()

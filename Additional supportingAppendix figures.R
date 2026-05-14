################################################################################
#
# MEDULLOBLASTOMA PROJECT — APPENDIX FIGURES
# Dataset: GSE85217
#
# Student: Ewan McGibbon
# Date: April 2026
#
# NOTE: Requires MB_main.R to have been run first
#
################################################################################

library(survival)
library(survminer)
library(cowplot)
library(ggplot2)
library(ggridges)
library(dplyr)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)

########################################################################################################################
#                                   SUPPLEMENTARY FIGURE 4
#                          Molecular subtype KM panels
########################################################################################################################

# -----------------------------------------------------------------------------
# Group 3 only
# -----------------------------------------------------------------------------

g3_data <- meta_sub[meta_sub$Subgroup == "Group3", ]
g3_data$Subtype <- factor(dplyr::recode(g3_data$Subtype,
                                        "Group3_alpha" = "Group 3-\u03b1",
                                        "Group3_beta"  = "Group 3-\u03b2",
                                        "Group3_gamma" = "Group 3-\u03b3"),
                          levels = c("Group 3-\u03b1","Group 3-\u03b2","Group 3-\u03b3"))x

fit_subtype_g3 <- survfit(Surv(OS..years., Dead) ~ Subtype, data = g3_data)
km_subtype_g3  <- ggsurvplot(fit_subtype_g3, data = g3_data,
                             title = "Group 3: molecular subtypes",
                             legend.title = "Subtype", legend = "top",
                             risk.table = FALSE, xlim = c(0,5), break.time.by = 1,
                             xlab = "Overall survival (years)", ylab = "Survival probability",
                             pval = TRUE, pval.coord = c(0.3,0.15),
                             palette = c("#FFD700","#FFA500","#FF6600"),
                             ggtheme = theme_classic(base_size = 11))

km_subtype_g3$plot <- km_subtype_g3$plot +
  scale_color_manual(values = c("#FFD700","#FFA500","#FF6600"),
                     labels = c("Group 3-\u03b1","Group 3-\u03b2","Group 3-\u03b3"))

# -----------------------------------------------------------------------------
# Group 4 only
# -----------------------------------------------------------------------------

g4_data <- meta_sub[meta_sub$Subgroup == "Group4", ]
g4_data$Subtype <- factor(dplyr::recode(g4_data$Subtype,
                                        "Group4_alpha" = "Group 4-\u03b1",
                                        "Group4_beta"  = "Group 4-\u03b2",
                                        "Group4_gamma" = "Group 4-\u03b3"),
                          levels = c("Group 4-\u03b1","Group 4-\u03b2","Group 4-\u03b3"))

fit_subtype_g4 <- survfit(Surv(OS..years., Dead) ~ Subtype, data = g4_data)
km_subtype_g4  <- ggsurvplot(fit_subtype_g4, data = g4_data,
                             title = "Group 4: molecular subtypes",
                             legend.title = "Subtype", legend = "top",
                             risk.table = FALSE, xlim = c(0,5), break.time.by = 1,
                             xlab = "Overall survival (years)", ylab = "Survival probability",
                             pval = TRUE, pval.coord = c(0.3,0.15),
                             palette = c("#006400","#228B22","#90EE90"),
                             ggtheme = theme_classic(base_size = 11))

km_subtype_g4$plot <- km_subtype_g4$plot +
  scale_color_manual(values = c("#006400","#228B22","#90EE90"),
                     labels = c("Group 4-\u03b1","Group 4-\u03b2","Group 4-\u03b3"))

# -----------------------------------------------------------------------------
# SHH only
# -----------------------------------------------------------------------------

shh_data <- meta_sub[meta_sub$Subgroup == "SHH", ]
shh_data$Subtype <- factor(dplyr::recode(shh_data$Subtype,
                                         "SHH_alpha" = "SHH-\u03b1",
                                         "SHH_beta"  = "SHH-\u03b2",
                                         "SHH_gamma" = "SHH-\u03b3",
                                         "SHH_delta" = "SHH-\u03b4"),
                           levels = c("SHH-\u03b1","SHH-\u03b2","SHH-\u03b3","SHH-\u03b4"))

fit_subtype_shh <- survfit(Surv(OS..years., Dead) ~ Subtype, data = shh_data)
km_subtype_shh  <- ggsurvplot(fit_subtype_shh, data = shh_data,
                              title = "SHH: molecular subtypes",
                              legend.title = "Subtype", legend = "top",
                              risk.table = FALSE, xlim = c(0,5), break.time.by = 1,
                              xlab = "Overall survival (years)", ylab = "Survival probability",
                              pval = TRUE, pval.coord = c(0.3,0.15),
                              ggtheme = theme_classic(base_size = 11))

km_subtype_shh$plot <- km_subtype_shh$plot +
  scale_color_discrete(labels = c("SHH-\u03b1","SHH-\u03b2","SHH-\u03b3","SHH-\u03b4"))

# -----------------------------------------------------------------------------
# Three panel figure
# -----------------------------------------------------------------------------

fig_supp_subtypes <- plot_grid(
  km_subtype_g3$plot, km_subtype_g4$plot, km_subtype_shh$plot,
  ncol = 3, nrow = 1,
  align = "hv",
  labels = c("A","B","C")
)

print(fig_supp_subtypes)

########################################################################################################################
#                                   SUPPLEMENTARY FIGURE 5
#                          Statistical power for subtype survival analyses
########################################################################################################################

# -----------------------------------------------------------------------------
# Log-rank power function
# d = total events, HR = hazard ratio between groups, alpha = significance level
# -----------------------------------------------------------------------------

power_logrank <- function(d, HR, alpha = 0.05) {
  theta <- log(HR)
  z     <- sqrt(d/4) * abs(theta)
  power <- pnorm(z - qnorm(1 - alpha/2)) + pnorm(-z - qnorm(1 - alpha/2))
  return(round(power, 3))
}

# -----------------------------------------------------------------------------
# Event counts from present analysis
# -----------------------------------------------------------------------------

n_g3      <- sum(!is.na(g3_data$OS..years.))
events_g3 <- sum(g3_data$Dead, na.rm = TRUE)

n_g4      <- sum(!is.na(g4_data$OS..years.))
events_g4 <- sum(g4_data$Dead, na.rm = TRUE)

n_shh      <- sum(!is.na(shh_data$OS..years.))
events_shh <- sum(shh_data$Dead, na.rm = TRUE)

cat("Group 3 — n:", n_g3, "| events:", events_g3, "\n")
cat("Group 4 — n:", n_g4, "| events:", events_g4, "\n")
cat("SHH     — n:", n_shh, "| events:", events_shh, "\n")

# -----------------------------------------------------------------------------
# Power estimates — present vs estimated complete dataset (HR = 2.0 assumed)
# -----------------------------------------------------------------------------

power_df <- data.frame(
  Subgroup = rep(c("Group 3", "Group 4", "SHH"), each = 2),
  Dataset  = rep(c("Present analysis", "Estimated complete"), 3),
  Power    = c(
    power_logrank(events_g3,  2.0), power_logrank(round(events_g3  * (nrow(g3_data)/n_g3)),  2.0),
    power_logrank(events_g4,  2.0), power_logrank(round(events_g4  * (nrow(g4_data)/n_g4)),  2.0),
    power_logrank(events_shh, 2.0), power_logrank(round(events_shh * (nrow(shh_data)/n_shh)), 2.0)
  )
)

power_df$Subgroup <- factor(power_df$Subgroup, levels = c("Group 3", "Group 4", "SHH"))
power_df$Dataset  <- factor(power_df$Dataset,  levels = c("Present analysis", "Estimated complete"))

# -----------------------------------------------------------------------------
# Power bar chart — Supplementary Figure 5
# -----------------------------------------------------------------------------

fig_supp_power <- ggplot(power_df, aes(x = Subgroup, y = Power, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "red", linewidth = 0.8) +
  annotate("text", x = 0.5, y = 0.82, label = "80% power threshold",
           hjust = 0, size = 3.5, colour = "red") +
  scale_fill_manual(values = c("Present analysis" = "#4575b4",
                               "Estimated complete" = "#d73027")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Subgroup", y = "Estimated power",
       title = "Statistical power for subtype survival analysis",
       fill = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")

print(fig_supp_power)

########################################################################################################################
#                                   SUPPLEMENTARY FIGURE 7
#                    G3/4 KM curves stratified by risk tertile within clinical subgroups
########################################################################################################################

tertile_pal <- c("Low" = "#2166ac", "Mid" = "#f4a582", "High" = "#d6604d")
km_theme    <- theme_classic(base_size = 11)

km_tertile <- function(data, title) {
  data$risk_tertile <- factor(data$risk_tertile, levels = c("Low", "Mid", "High"))
  fit <- survfit(Surv(OS..years., Dead) ~ risk_tertile, data = data)
  lr  <- survdiff(Surv(OS..years., Dead) ~ risk_tertile, data = data)
  p   <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
  
  ggsurvplot(
    fit, data = data,
    palette       = tertile_pal,
    conf.int      = FALSE,
    risk.table    = FALSE,
    legend.title  = "Risk tertile",
    legend.labs   = c("Low", "Mid", "High"),
    legend        = "top",
    xlim          = c(0, 10),
    break.time.by = 2,
    xlab          = "Overall survival (years)",
    ylab          = "Survival probability",
    title         = title,
    pval          = FALSE,
    ggtheme       = km_theme
  )$plot +
    annotate("text", x = 0.3, y = 0.1,
             label = paste0("Log-rank p = ", format.pval(p, digits = 2, eps = 0.001)),
             hjust = 0, size = 3.5)
}

A <- km_tertile(dat_glm_micro[dat_glm_micro$Subgroup == "Group3", ],           "Group 3")
B <- km_tertile(dat_glm_micro[dat_glm_micro$Subgroup == "Group4", ],           "Group 4")
C <- km_tertile(dat_glm_micro[dat_glm_micro$MYC == 1, ],                       "MYC amplified")
D <- km_tertile(dat_glm_micro[dat_glm_micro$MYCN == 1, ],                      "MYCN amplified")
E <- km_tertile(dat_glm_micro[dat_glm_micro$histology == "Classic", ],         "Classic histology")
F_plot <- km_tertile(
  dat_glm_micro[as.character(dat_glm_micro$Met.status..1.Met..0.M0.) == "1", ], "Metastatic (M+)")

quartz(type = "pdf", file = "G34_tertile_KM_by_variable.pdf", width = 14, height = 9)

print(plot_grid(A, B, C, D, E, F_plot,
                ncol = 3, nrow = 2,
                labels = c("A", "B", "C", "D", "E", "F"),
                label_size = 12))

dev.off()

########################################################################################################################
#                                   SUPPLEMENTARY FIGURE 8
#                    G3/4 risk score distribution by molecular subtype
########################################################################################################################

toPlot <- data.frame(
  risk_score = dat_glm_micro$risk_score,
  Subtype    = dat_glm_micro$Subtype,
  Subgroup   = dat_glm_micro$Subgroup
) %>%
  filter(!is.na(Subtype), Subtype != "Unknown") %>%
  mutate(Subtype = factor(dplyr::recode(Subtype, !!!subtype_labels),
                          levels = rev(names(subtype_bar_cols))))

p_anova <- anova(aov(risk_score ~ Subtype, data = toPlot))$`Pr(>F)`[1]
rs_range <- range(toPlot$risk_score, na.rm = TRUE)

# A — Ridge plot
p_ridge <- ggplot(toPlot, aes(x = risk_score, y = Subtype, fill = Subtype)) +
  geom_density_ridges(scale = 1.5, rel_min_height = 0.001, alpha = 0.75) +
  scale_fill_manual(values = subtype_bar_cols) +
  xlim(rs_range) +
  annotate("text", x = rs_range[2] - 0.1, y = 0.6,
           label = fmt_p(p_anova), hjust = 1, size = 3.5) +
  labs(x = "Risk score", y = "Subtype",
       title = "Risk score distribution by molecular subtype — Group 3/4") +
  theme_bw() +
  theme(legend.position = "none", panel.grid = element_blank())

# B — Boxplot with jitter
p_box <- ggplot(toPlot, aes(x = Subtype, y = risk_score, fill = Subtype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 0.8) +
  scale_fill_manual(values = subtype_bar_cols) +
  labs(x = "Subtype", y = "Risk score",
       title = "Risk score by molecular subtype — Group 3/4") +
  theme_bw() +
  theme(legend.position = "none", panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1))

quartz(type = "pdf", file = "G34_risk_score_by_subtype.pdf", width = 12, height = 6)

print(ggarrange(p_ridge, p_box,
                ncol = 2, nrow = 1,
                labels = c("A", "B")))

dev.off()






###########################################################################
# SHH cohort — Risk score by molecular subtype (separate figure)
###########################################################################

toPlot_shh <- data.frame(
  risk_score = dat_glm_shh$risk_score,
  Subtype    = dat_glm_shh$Subtype
) %>%
  filter(!is.na(Subtype), Subtype != "Unknown") %>%
  mutate(Subtype = factor(dplyr::recode(Subtype, !!!shh_subtype_labels),
                          levels = rev(names(shh_subtype_bar_cols))))

p_anova_shh <- anova(aov(risk_score ~ Subtype, data = toPlot_shh))$`Pr(>F)`[1]
rs_range_shh <- range(toPlot_shh$risk_score, na.rm = TRUE)

# A — Ridge plot
p_ridge_shh <- ggplot(toPlot_shh, aes(x = risk_score, y = Subtype, fill = Subtype)) +
  geom_density_ridges(scale = 1.5, rel_min_height = 0.001, alpha = 0.75) +
  scale_fill_manual(values = shh_subtype_bar_cols) +
  xlim(rs_range_shh) +
  annotate("text", x = rs_range_shh[2] - 0.1, y = 0.6,
           label = fmt_p(p_anova_shh), hjust = 1, size = 3.5) +
  labs(x = "Risk score", y = "Subtype",
       title = "Risk score distribution by molecular subtype — SHH") +
  theme_bw() +
  theme(legend.position = "none", panel.grid = element_blank())

# B — Boxplot with jitter
p_box_shh <- ggplot(toPlot_shh, aes(x = Subtype, y = risk_score, fill = Subtype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 0.8) +
  scale_fill_manual(values = shh_subtype_bar_cols) +
  labs(x = "Subtype", y = "Risk score",
       title = "Risk score by molecular subtype — SHH") +
  theme_bw() +
  theme(legend.position = "none", panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1))

quartz(type = "pdf", file = "SHH_risk_score_by_subtype.pdf", width = 12, height = 6)

print(ggarrange(p_ridge_shh, p_box_shh,
                ncol = 2, nrow = 1,
                labels = c("A", "B")))

dev.off()

# ── Shared subtype labels and colours ─────────────────────────────────────────
subtype_col_fun <- c(
  "Group3_alpha" = "#FFD700", "Group3_beta"  = "#FFA500", "Group3_gamma" = "#FF6600",
  "Group4_alpha" = "#006400", "Group4_beta"  = "#228B22", "Group4_gamma" = "#90EE90",
  "Unknown"      = "lightgrey"
)

subtype_labels <- c(
  "Group3_alpha" = "Group 3-\u03b1", "Group3_beta"  = "Group 3-\u03b2", "Group3_gamma" = "Group 3-\u03b3",
  "Group4_alpha" = "Group 4-\u03b1", "Group4_beta"  = "Group 4-\u03b2", "Group4_gamma" = "Group 4-\u03b3"
)

subtype_bar_cols <- c(
  "Group 3-\u03b1" = "#FFD700", "Group 3-\u03b2" = "#FFA500", "Group 3-\u03b3" = "#FF6600",
  "Group 4-\u03b1" = "#006400", "Group 4-\u03b2" = "#228B22", "Group 4-\u03b3" = "#90EE90"
)

# ── 1. Merge MYC and MYCN amplification ───────────────────────────────────────
dat_glm_g34$MYC  <- ifelse(dat_glm_g34$Study_ID %in% myc_amp_ids,  1, 0)
dat_glm_g34$MYCN <- ifelse(dat_glm_g34$Study_ID %in% mycn_amp_ids, 1, 0)

cat("MYC amplification in G3/4:\n");  print(table(dat_glm_g34$MYC,  useNA = "ifany"))
cat("MYCN amplification in G3/4:\n"); print(table(dat_glm_g34$MYCN, useNA = "ifany"))

# ── 2. Prepare data for oncoprint ─────────────────────────────────────────────
dat_plot      <- dat_glm_g34[order(dat_glm_g34$risk_score), ]
dat_plot$Dead <- as.numeric(dat_plot$Dead)

# ── 3. Clean NAs ──────────────────────────────────────────────────────────────
dat_plot$Gender    <- ifelse(is.na(dat_plot$Gender),    "Unknown", as.character(dat_plot$Gender))
dat_plot$histology <- ifelse(is.na(dat_plot$histology) | dat_plot$histology == "", "Unknown", dat_plot$histology)
dat_plot$Subgroup  <- ifelse(is.na(dat_plot$Subgroup),  "Unknown", as.character(dat_plot$Subgroup))
dat_plot$Subtype   <- ifelse(is.na(dat_plot$Subtype)   | dat_plot$Subtype   == "", "Unknown", as.character(dat_plot$Subtype))
dat_plot$Met.status..1.Met..0.M0. <- ifelse(
  is.na(dat_plot$Met.status..1.Met..0.M0.), "Unknown",
  as.character(dat_plot$Met.status..1.Met..0.M0.)
)
dat_plot$MYC         <- as.character(dat_plot$MYC)
dat_plot$MYCN        <- as.character(dat_plot$MYCN)
dat_plot$Dead_status <- ifelse(is.na(dat_plot$Dead), "Unknown",
                               ifelse(dat_plot$Dead == 1, "Dead", "Censored"))

# ── 4. Colour palettes ────────────────────────────────────────────────────────
tertile_col_fun  <- c("Low" = "#2166ac", "Mid" = "#f4a582", "High" = "#d6604d")
gender_col_fun   <- c("M" = "white",     "F" = "black",     "Unknown" = "lightgrey")
subgroup_col_fun <- c("Group3" = "#FFFF00", "Group4" = "#008000", "Unknown" = "lightgrey")
hist_col_fun     <- c("Classic" = "white", "Desmoplastic" = "#fc8d59",
                      "LCA" = "black",     "MBEN" = "red4", "Unknown" = "lightgrey")
mstage_col_fun   <- c("0" = "white", "1" = "black", "Unknown" = "lightgrey")
myc_col_fun      <- c("0" = "white", "1" = "#e31a1c")
mycn_col_fun     <- c("0" = "white", "1" = "#ff7f00")
dead_col_fun     <- c("Dead" = "black", "Censored" = "white", "Unknown" = "lightgrey")

# ── 5. Age dot annotation ─────────────────────────────────────────────────────
Age_capped <- pmin(dat_plot$Age, 20)

age_dot_anno <- HeatmapAnnotation(
  Age = anno_points(
    Age_capped, ylim = c(0, 20),
    pch = ifelse(is.na(dat_plot$Age), 4, 20),
    gp = gpar(col = ifelse(is.na(dat_plot$Age), "red", "grey")),
    axis_param = list(at = c(0, 5, 10, 15, 20)),
    baseline = 3, baseline_gp = gpar(col = "grey", lty = 2, lwd = 1)
  ),
  annotation_name_side = "left", annotation_height = unit(3.2, "cm")
)

# ── 6. OS dot annotation ──────────────────────────────────────────────────────
OS <- pmin(dat_plot$OS..years., 10)
pchSetup <- ifelse(dat_plot$Dead == 1 & dat_plot$OS..years. <  10,  4,
                   ifelse(dat_plot$Dead == 0 & dat_plot$OS..years. <  10, 20, 2))

OS_dot_anno <- HeatmapAnnotation(
  OS = anno_points(
    OS, pch = pchSetup,
    gp = gpar(col = ifelse(dat_plot$Dead == 1, "firebrick3", "lightgrey")),
    ylim = c(0, 10), axis_param = list(at = c(0, 5, 10)),
    baseline_gp = gpar(col = "lightgrey", lty = 2, lwd = 1)
  ),
  annotation_name_side = "left", annotation_height = unit(3, "cm")
)

# ── 7. Categorical bar annotations ───────────────────────────────────────────
cat_anno <- HeatmapAnnotation(
  `Risk tertile` = anno_simple(as.character(dat_plot$risk_tertile), col = tertile_col_fun,  border = TRUE),
  Gender         = anno_simple(dat_plot$Gender,                     col = gender_col_fun,   border = TRUE),
  Subgroup       = anno_simple(dat_plot$Subgroup,                   col = subgroup_col_fun, border = TRUE),
  Subtype        = anno_simple(dat_plot$Subtype,                    col = subtype_col_fun,  border = TRUE),
  Histology      = anno_simple(dat_plot$histology,                  col = hist_col_fun,     border = TRUE),
  `M-stage`      = anno_simple(dat_plot$Met.status..1.Met..0.M0.,  col = mstage_col_fun,   border = TRUE),
  MYC            = anno_simple(dat_plot$MYC,                        col = myc_col_fun,      border = TRUE),
  MYCN           = anno_simple(dat_plot$MYCN,                       col = mycn_col_fun,     border = TRUE),
  Outcome        = anno_simple(dat_plot$Dead_status,                col = dead_col_fun,     border = TRUE),
  annotation_name_side = "left"
)

# ── 8. Base heatmap ───────────────────────────────────────────────────────────
rs_matrix <- matrix(dat_plot$risk_score, nrow = 1)

ht_g34 <- Heatmap(
  rs_matrix, name = "Risk score",
  col = colorRamp2(
    c(min(dat_plot$risk_score), 0, max(dat_plot$risk_score)),
    c("#2166ac", "white", "#d6604d")
  ),
  border = TRUE, show_column_names = FALSE,
  cluster_columns = FALSE, cluster_rows = FALSE,
  row_names_side = "left", row_names_gp = gpar(fontsize = 10),
  height = unit(0.6, "cm"), show_heatmap_legend = TRUE,
  top_annotation = c(cat_anno, age_dot_anno, OS_dot_anno)
)

# ── 9. Manual legends ────────────────────────────────────────────────────────
lgd_tertile  <- Legend(labels = c("Low", "Mid", "High"), title = "Risk tertile",
                       legend_gp = gpar(fill = c("#2166ac", "#f4a582", "#d6604d")), border = TRUE)
lgd_subgroup <- Legend(labels = c("Group3", "Group4", "Unknown"), title = "Subgroup",
                       legend_gp = gpar(fill = c("#FFFF00", "#008000", "lightgrey")), border = TRUE)
lgd_subtype  <- Legend(labels = c("Group 3-\u03b1", "Group 3-\u03b2", "Group 3-\u03b3",
                                  "Group 4-\u03b1", "Group 4-\u03b2", "Group 4-\u03b3", "Unknown"),
                       title = "Subtype",
                       legend_gp = gpar(fill = c("#FFD700", "#FFA500", "#FF6600",
                                                 "#006400", "#228B22", "#90EE90", "lightgrey")), border = TRUE)
lgd_gender   <- Legend(labels = c("Male", "Female", "Unknown"), title = "Gender",
                       legend_gp = gpar(fill = c("white", "black", "lightgrey")), border = TRUE)
lgd_hist     <- Legend(labels = c("Classic", "Desmoplastic", "LCA", "MBEN", "Unknown"), title = "Histology",
                       legend_gp = gpar(fill = c("white", "#fc8d59", "black", "red4", "lightgrey")), border = TRUE)
lgd_mstage   <- Legend(labels = c("M0", "M+", "Unknown"), title = "M-stage",
                       legend_gp = gpar(fill = c("white", "black", "lightgrey")), border = TRUE)
lgd_myc      <- Legend(labels = c("Not amplified", "Amplified"), title = "MYC",
                       legend_gp = gpar(fill = c("white", "#e31a1c")), border = TRUE)
lgd_mycn     <- Legend(labels = c("Not amplified", "Amplified"), title = "MYCN",
                       legend_gp = gpar(fill = c("white", "#ff7f00")), border = TRUE)
lgd_outcome  <- Legend(labels = c("Dead", "Censored", "Unknown"), title = "Outcome",
                       legend_gp = gpar(fill = c("black", "white", "lightgrey")), border = TRUE)

# ── 10. Ridge plots ───────────────────────────────────────────────────────────
fmt_p <- function(p) {
  if (p < 2.2e-16)    "p < 2.2e-16"
  else if (p < 0.001) paste0("p = ", formatC(p, format = "e", digits = 2))
  else                paste0("p = ", round(p, 3))
}

rs_range    <- range(dat_glm_g34$risk_score, na.rm = TRUE)
ridge_theme <- theme_bw() + theme(legend.position = "none", panel.grid = element_blank())

toPlot_A <- data.frame(risk_score = dat_glm_g34$risk_score, group = dat_glm_g34$Subgroup) %>%
  filter(!is.na(group), group != "Unknown")
p_A <- fmt_p(t.test(risk_score ~ group, data = toPlot_A)$p.value)
A <- ggplot(toPlot_A, aes(risk_score, group)) +
  geom_density_ridges(aes(fill = group), scale = 1.5, rel_min_height = 0.001, alpha = 0.45) +
  scale_fill_manual(values = c("Group3" = "#FFFF00", "Group4" = "#008000")) +
  xlim(rs_range) + ridge_theme + labs(x = "Risk score", y = "Subgroup") +
  annotate("text", x = rs_range[2] - 0.1, y = 3.3, label = p_A, hjust = 1, size = 3.5)

toPlot_B <- data.frame(risk_score = dat_glm_g34$risk_score, group = dat_glm_g34$histology) %>%
  filter(!is.na(group), !group %in% c("Unknown", ""))
p_B <- fmt_p(anova(aov(risk_score ~ group, data = toPlot_B))$`Pr(>F)`[1])
hist_cols <- c("Classic" = "grey80", "Desmoplastic" = "#fc8d59", "LCA" = "black", "MBEN" = "red4")
B <- ggplot(toPlot_B, aes(risk_score, group)) +
  geom_density_ridges(aes(fill = group), scale = 1.5, rel_min_height = 0.001, alpha = 0.45) +
  scale_fill_manual(values = hist_cols) +
  xlim(rs_range) + ridge_theme + labs(x = "Risk score", y = "Histology") +
  annotate("text", x = rs_range[2] - 0.1, y = 4.3, label = p_B, hjust = 1, size = 3.5)

toPlot_C <- data.frame(risk_score = dat_glm_g34$risk_score,
                       group = as.character(dat_glm_g34$Met.status..1.Met..0.M0.)) %>%
  filter(!is.na(group), group != "Unknown")
p_C <- fmt_p(t.test(risk_score ~ group, data = toPlot_C)$p.value)
C <- ggplot(toPlot_C, aes(risk_score, group)) +
  geom_density_ridges(aes(fill = group), scale = 1.5, rel_min_height = 0.001, alpha = 0.45) +
  scale_fill_manual(values = c("0" = "grey80", "1" = "black")) +
  scale_y_discrete(labels = c("0" = "M0", "1" = "M+")) +
  xlim(rs_range) + ridge_theme + labs(x = "Risk score", y = "M-stage") +
  annotate("text", x = rs_range[2] - 0.1, y = 3.3, label = p_C, hjust = 1, size = 3.5)

toPlot_D <- data.frame(risk_score = dat_glm_g34$risk_score, group = dat_glm_g34$Gender) %>%
  filter(!is.na(group), group != "Unknown")
p_D <- fmt_p(t.test(risk_score ~ group, data = toPlot_D)$p.value)
D <- ggplot(toPlot_D, aes(risk_score, group)) +
  geom_density_ridges(aes(fill = group), scale = 1.5, rel_min_height = 0.001, alpha = 0.45) +
  scale_fill_manual(values = c("F" = "#d73027", "M" = "#4575b4")) +
  xlim(rs_range) + ridge_theme + labs(x = "Risk score", y = "Gender") +
  annotate("text", x = rs_range[2] - 0.1, y = 3.3, label = p_D, hjust = 1, size = 3.5)

toPlot_E <- data.frame(risk_score = dat_glm_g34$risk_score,
                       group = as.character(dat_glm_g34$MYC)) %>%
  filter(!is.na(group))
p_E <- fmt_p(t.test(risk_score ~ group, data = toPlot_E)$p.value)
E <- ggplot(toPlot_E, aes(risk_score, group)) +
  geom_density_ridges(aes(fill = group), scale = 1.5, rel_min_height = 0.001, alpha = 0.45) +
  scale_fill_manual(values = c("0" = "grey80", "1" = "#e31a1c")) +
  scale_y_discrete(labels = c("0" = "Not amplified", "1" = "Amplified")) +
  xlim(rs_range) + ridge_theme + labs(x = "Risk score", y = "MYC") +
  annotate("text", x = rs_range[2] - 0.1, y = 3.3, label = p_E, hjust = 1, size = 3.5)

toPlot_F <- data.frame(risk_score = dat_glm_g34$risk_score,
                       group = as.character(dat_glm_g34$MYCN)) %>%
  filter(!is.na(group))
p_F <- fmt_p(t.test(risk_score ~ group, data = toPlot_F)$p.value)
F_plot <- ggplot(toPlot_F, aes(risk_score, group)) +
  geom_density_ridges(aes(fill = group), scale = 1.5, rel_min_height = 0.001, alpha = 0.45) +
  scale_fill_manual(values = c("0" = "grey80", "1" = "#ff7f00")) +
  scale_y_discrete(labels = c("0" = "Not amplified", "1" = "Amplified")) +
  xlim(rs_range) + ridge_theme + labs(x = "Risk score", y = "MYCN") +
  annotate("text", x = rs_range[2] - 0.1, y = 3.3, label = p_F, hjust = 1, size = 3.5)

# ── 11. Combined PDF ──────────────────────────────────────────────────────────
pdf("G34_clinical_annotation_w_subtypes.pdf", width = 14, height = 8)

draw(
  ht_g34,
  heatmap_legend_side    = "right",
  annotation_legend_side = "right",
  annotation_legend_list = list(lgd_tertile, lgd_subgroup, lgd_subtype, lgd_gender,
                                lgd_hist, lgd_mstage, lgd_myc, lgd_mycn, lgd_outcome),
  padding = unit(c(2, 30, 2, 2), "mm")
)

decorate_heatmap_body("Risk score", {
  grid.lines(x = unit(c(1/3, 1/3), "npc"), y = unit(c(0, 1), "npc"),
             gp = gpar(col = "black", lty = "dashed", lwd = 1.5))
  grid.lines(x = unit(c(2/3, 2/3), "npc"), y = unit(c(0, 1), "npc"),
             gp = gpar(col = "black", lty = "dashed", lwd = 1.5))
})

decorate_annotation("Age", {
  grid.lines(x = unit(c(0, 1), "npc"), y = unit(c(3, 3), "native"),
             gp = gpar(col = "black", lty = "dotted", lwd = 1))
})

decorate_annotation("OS", {
  grid.lines(x = unit(c(0, 1), "npc"), y = unit(c(5, 5), "native"),
             gp = gpar(col = "black", lty = "dotted", lwd = 1))
})

grid.newpage()
print(ggarrange(A, B, C, D, E, F_plot, ncol = 2, nrow = 3, labels = c("A","B","C","D","E","F")))
dev.off()

########################################################################################################################
#                                   SUPPLEMENTARY FIGURE 9
#                    SHH oncoprint — ordered by ascending Cox risk score
#                    + Ridge plots of risk score by clinical variable
########################################################################################################################

# ── Shared subtype labels and colours ─────────────────────────────────────────
shh_subtype_labels <- c(
  "SHH_alpha" = "SHH-\u03b1",
  "SHH_beta"  = "SHH-\u03b2",
  "SHH_gamma" = "SHH-\u03b3",
  "SHH_delta" = "SHH-\u03b4"
)

shh_subtype_col_fun <- c(
  "SHH_alpha" = "#7b2d8b",
  "SHH_beta"  = "#e7298a",
  "SHH_gamma" = "#66a61e",
  "SHH_delta" = "#e6ab02",
  "Unknown"   = "lightgrey"
)

shh_subtype_bar_cols <- c(
  "SHH-\u03b1" = "#7b2d8b",
  "SHH-\u03b2" = "#e7298a",
  "SHH-\u03b3" = "#66a61e",
  "SHH-\u03b4" = "#e6ab02"
)

# ── 1. Merge MYC and MYCN amplification into dat_glm_shh ─────────────────────
dat_glm_shh$MYC  <- ifelse(dat_glm_shh$Study_ID %in% myc_amp_ids,  1, 0)
dat_glm_shh$MYCN <- ifelse(dat_glm_shh$Study_ID %in% mycn_amp_ids, 1, 0)

cat("MYC amplification in SHH:\n");  print(table(dat_glm_shh$MYC,  useNA = "ifany"))
cat("MYCN amplification in SHH:\n"); print(table(dat_glm_shh$MYCN, useNA = "ifany"))

# ── 2. Prepare data for oncoprint ─────────────────────────────────────────────
dat_plot      <- dat_glm_shh[order(dat_glm_shh$risk_score), ]
dat_plot$Dead <- as.numeric(dat_plot$Dead)

# ── 3. Clean NAs ─────────────────────────────────────────────────────────────
dat_plot$Gender    <- ifelse(is.na(dat_plot$Gender),    "Unknown", as.character(dat_plot$Gender))
dat_plot$histology <- ifelse(is.na(dat_plot$histology) | dat_plot$histology == "", "Unknown", dat_plot$histology)
dat_plot$Subtype   <- ifelse(is.na(dat_plot$Subtype)   | dat_plot$Subtype   == "", "Unknown", dat_plot$Subtype)
dat_plot$Met.status..1.Met..0.M0. <- ifelse(
  is.na(dat_plot$Met.status..1.Met..0.M0.), "Unknown",
  as.character(dat_plot$Met.status..1.Met..0.M0.)
)
dat_plot$MYC         <- as.character(dat_plot$MYC)
dat_plot$MYCN        <- as.character(dat_plot$MYCN)
dat_plot$Dead_status <- ifelse(is.na(dat_plot$Dead), "Unknown",
                               ifelse(dat_plot$Dead == 1, "Dead", "Censored"))

# ── 4. Colour palettes ────────────────────────────────────────────────────────
tertile_col_fun <- c("Low" = "#2166ac", "Mid" = "#f4a582", "High" = "#d6604d")
gender_col_fun  <- c("M" = "white", "F" = "black", "Unknown" = "lightgrey")
hist_col_fun    <- c("Classic" = "white", "Desmoplastic" = "#fc8d59",
                     "LCA" = "black", "MBEN" = "red4", "Unknown" = "lightgrey")
mstage_col_fun  <- c("0" = "white", "1" = "black", "Unknown" = "lightgrey")
mycn_col_fun    <- c("0" = "white", "1" = "#ff7f00")
dead_col_fun    <- c("Dead" = "black", "Censored" = "white", "Unknown" = "lightgrey")

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
  annotation_height = unit(3.2, "cm")
)

# ── 6. OS dot annotation (capped at 10 years) ────────────────────────────────
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
  annotation_height = unit(3, "cm")
)

# ── 7. Categorical bar annotations ───────────────────────────────────────────
cat_anno <- HeatmapAnnotation(
  `Risk tertile` = anno_simple(as.character(dat_plot$risk_tertile), col = tertile_col_fun,     border = TRUE),
  Gender         = anno_simple(dat_plot$Gender,                     col = gender_col_fun,      border = TRUE),
  Subtype        = anno_simple(dat_plot$Subtype,                    col = shh_subtype_col_fun, border = TRUE),
  Histology      = anno_simple(dat_plot$histology,                  col = hist_col_fun,        border = TRUE),
  `M-stage`      = anno_simple(dat_plot$Met.status..1.Met..0.M0.,  col = mstage_col_fun,      border = TRUE),
  MYCN           = anno_simple(dat_plot$MYCN,                       col = mycn_col_fun,        border = TRUE),
  Outcome        = anno_simple(dat_plot$Dead_status,                col = dead_col_fun,        border = TRUE),
  annotation_name_side = "left"
)

# ── 8. Base heatmap ───────────────────────────────────────────────────────────
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

# ── 9. Manual legends ────────────────────────────────────────────────────────
lgd_tertile <- Legend(labels = c("Low", "Mid", "High"),
                      title = "Risk tertile",
                      legend_gp = gpar(fill = c("#2166ac", "#f4a582", "#d6604d")), border = TRUE)
lgd_subtype <- Legend(labels = c("SHH-\u03b1", "SHH-\u03b2", "SHH-\u03b3", "SHH-\u03b4", "Unknown"),
                      title = "Subtype",
                      legend_gp = gpar(fill = c("#7b2d8b", "#e7298a", "#66a61e", "#e6ab02", "lightgrey")), border = TRUE)
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
lgd_outcome <- Legend(labels = c("Dead", "Censored", "Unknown"),
                      title = "Outcome",
                      legend_gp = gpar(fill = c("black", "white", "lightgrey")), border = TRUE)

# ── 10. Ridge plots ───────────────────────────────────────────────────────────
fmt_p <- function(p) {
  if (p < 2.2e-16)    "p < 2.2e-16"
  else if (p < 0.001) paste0("p = ", formatC(p, format = "e", digits = 2))
  else                paste0("p = ", round(p, 3))
}

rs_range    <- range(dat_glm_shh$risk_score, na.rm = TRUE)
ridge_theme <- theme_bw() + theme(legend.position = "none", panel.grid = element_blank())

# A — Histology
toPlot_A <- data.frame(risk_score = dat_glm_shh$risk_score, group = dat_glm_shh$histology) %>%
  filter(!is.na(group), !group %in% c("Unknown", ""))
p_A <- fmt_p(anova(aov(risk_score ~ group, data = toPlot_A))$`Pr(>F)`[1])
hist_cols <- c("Classic" = "grey80", "Desmoplastic" = "#fc8d59", "LCA" = "black", "MBEN" = "red4")

A <- ggplot(toPlot_A, aes(risk_score, group)) +
  geom_density_ridges(aes(fill = group), scale = 1.5, rel_min_height = 0.001, alpha = 0.45) +
  scale_fill_manual(values = hist_cols) +
  xlim(rs_range) + ridge_theme +
  labs(x = "Risk score", y = "Histology") +
  annotate("text", x = rs_range[2] - 0.1, y = length(unique(toPlot_A$group)) + 1.3,
           label = p_A, hjust = 1, size = 3.5)

# B — M-stage
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

# C — Gender
toPlot_C <- data.frame(risk_score = dat_glm_shh$risk_score, group = dat_glm_shh$Gender) %>%
  filter(!is.na(group), group != "Unknown")
p_C <- fmt_p(t.test(risk_score ~ group, data = toPlot_C)$p.value)

C <- ggplot(toPlot_C, aes(risk_score, group)) +
  geom_density_ridges(aes(fill = group), scale = 1.5, rel_min_height = 0.001, alpha = 0.45) +
  scale_fill_manual(values = c("F" = "#d73027", "M" = "#4575b4")) +
  xlim(rs_range) + ridge_theme +
  labs(x = "Risk score", y = "Gender") +
  annotate("text", x = rs_range[2] - 0.1, y = 3.3, label = p_C, hjust = 1, size = 3.5)

# D — MYCN amplification
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

# ── 11. Combined PDF — page 1: oncoprint, page 2: ridge plots ─────────────────
quartz(type = "pdf", file = "SHH_clinical_annotation.pdf", width = 16, height = 10)

draw(
  ht_shh,
  heatmap_legend_side    = "right",
  annotation_legend_side = "right",
  annotation_legend_list = list(lgd_tertile, lgd_subtype, lgd_gender,
                                lgd_hist, lgd_mstage, lgd_mycn, lgd_outcome),
  padding = unit(c(2, 30, 2, 2), "mm")
)

decorate_heatmap_body("Risk score", {
  grid.lines(x = unit(c(1/3, 1/3), "npc"),
             y = unit(c(0, 1), "npc"),
             gp = gpar(col = "black", lty = "dashed", lwd = 1.5))
  grid.lines(x = unit(c(2/3, 2/3), "npc"),
             y = unit(c(0, 1), "npc"),
             gp = gpar(col = "black", lty = "dashed", lwd = 1.5))
})

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

###########################################################################
# SHH cohort — Risk score by molecular subtype (separate figure)
###########################################################################

toPlot_shh <- data.frame(
  risk_score = dat_glm_shh$risk_score,
  Subtype    = dat_glm_shh$Subtype
) %>%
  filter(!is.na(Subtype), Subtype != "Unknown") %>%
  mutate(Subtype = factor(dplyr::recode(Subtype, !!!shh_subtype_labels),
                          levels = rev(names(shh_subtype_bar_cols))))

p_anova_shh  <- anova(aov(risk_score ~ Subtype, data = toPlot_shh))$`Pr(>F)`[1]
rs_range_shh <- range(toPlot_shh$risk_score, na.rm = TRUE)

# A — Ridge plot
p_ridge_shh <- ggplot(toPlot_shh, aes(x = risk_score, y = Subtype, fill = Subtype)) +
  geom_density_ridges(scale = 1.5, rel_min_height = 0.001, alpha = 0.75) +
  scale_fill_manual(values = shh_subtype_bar_cols) +
  xlim(rs_range_shh) +
  annotate("text", x = rs_range_shh[2] - 0.1, y = 0.6,
           label = fmt_p(p_anova_shh), hjust = 1, size = 3.5) +
  labs(x = "Risk score", y = "Subtype",
       title = "Risk score distribution by molecular subtype — SHH") +
  theme_bw() +
  theme(legend.position = "none", panel.grid = element_blank())

# B — Boxplot with jitter
p_box_shh <- ggplot(toPlot_shh, aes(x = Subtype, y = risk_score, fill = Subtype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 0.8) +
  scale_fill_manual(values = shh_subtype_bar_cols) +
  labs(x = "Subtype", y = "Risk score",
       title = "Risk score by molecular subtype — SHH") +
  theme_bw() +
  theme(legend.position = "none", panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1))

quartz(type = "pdf", file = "SHH_risk_score_by_subtype.pdf", width = 12, height = 6)

print(ggarrange(p_ridge_shh, p_box_shh,
                ncol = 2, nrow = 1,
                labels = c("A", "B")))

dev.off()


###########################################################################
# Supplementary Table — Univariable Cox regression by molecular subtype
# Dataset: GSE85217 (Cavalli microarray cohort)
# Input: meta_g34 and meta_shh from MB_main.R
###########################################################################

# ── G3/4 subtypes ─────────────────────────────────────────────────────────────
meta_g34$Subtype <- relevel(factor(as.character(meta_g34$Subtype),
                                   levels = c("Group3_alpha","Group3_beta","Group3_gamma",
                                              "Group4_alpha","Group4_beta","Group4_gamma")),
                            ref = "Group3_alpha")
cox_subtype_g34 <- coxph(Surv(OS..years., Dead) ~ Subtype, data = meta_g34)

# ── SHH subtypes ──────────────────────────────────────────────────────────────
meta_shh$Subtype <- relevel(factor(as.character(meta_shh$Subtype),
                                   levels = c("SHH_alpha","SHH_beta","SHH_gamma","SHH_delta")),
                            ref = "SHH_alpha")
cox_subtype_shh <- coxph(Surv(OS..years., Dead) ~ Subtype, data = meta_shh)

# ── Format combined table ─────────────────────────────────────────────────────
fmt_hr <- function(model) {
  s    <- summary(model)
  coef <- s$coefficients
  ci   <- s$conf.int
  data.frame(
    Characteristic = rownames(coef),
    HR_CI = paste0(sprintf("%.2f", ci[, "exp(coef)"]),
                   " (",
                   sprintf("%.2f", ci[, "lower .95"]),
                   "\u2013",
                   sprintf("%.2f", ci[, "upper .95"]),
                   ")"),
    p_value = ifelse(coef[, "Pr(>|z|)"] < 0.001, "<0.001",
                     sprintf("%.3f", coef[, "Pr(>|z|)"])),
    stringsAsFactors = FALSE
  )
}

clean_labels <- function(df, prefix) {
  df$Characteristic <- gsub(prefix, "", df$Characteristic)
  df[, c("Characteristic", "HR_CI", "p_value")]
}

g34_tbl <- clean_labels(fmt_hr(cox_subtype_g34), "Subtype")
shh_tbl <- clean_labels(fmt_hr(cox_subtype_shh), "Subtype")

# Clean subtype labels to Greek
g34_tbl$Characteristic <- dplyr::recode(g34_tbl$Characteristic,
                                        "Group3_beta"  = "Group 3-\u03b2",
                                        "Group3_gamma" = "Group 3-\u03b3",
                                        "Group4_alpha" = "Group 4-\u03b1",
                                        "Group4_beta"  = "Group 4-\u03b2",
                                        "Group4_gamma" = "Group 4-\u03b3"
)
shh_tbl$Characteristic <- dplyr::recode(shh_tbl$Characteristic,
                                        "SHH_beta"  = "SHH-\u03b2",
                                        "SHH_gamma" = "SHH-\u03b3",
                                        "SHH_delta" = "SHH-\u03b4"
)

ref_g34    <- data.frame(Characteristic = "Group 3-\u03b1 (Reference)", HR_CI = "\u2014", p_value = "\u2014")
ref_shh    <- data.frame(Characteristic = "SHH-\u03b1 (Reference)",     HR_CI = "\u2014", p_value = "\u2014")
header_g34 <- data.frame(Characteristic = "Group 3/4 (n=377, events=118)", HR_CI = "", p_value = "")
header_shh <- data.frame(Characteristic = "SHH (n=172, events=38)",        HR_CI = "", p_value = "")

g34_tbl    <- rbind(ref_g34, g34_tbl)
shh_tbl    <- rbind(ref_shh, shh_tbl)
supp_table <- rbind(header_g34, g34_tbl, header_shh, shh_tbl)
colnames(supp_table) <- c("Characteristic", "HR (95% CI)", "p-value")

# ── Export as Word document ───────────────────────────────────────────────────
ft <- flextable(supp_table) %>%
  bold(i = c(1, 9), bold = TRUE) %>%
  bold(i = c(3, 7), j = 3, bold = TRUE) %>%
  set_header_labels(
    Characteristic = "Characteristic",
    `HR (95% CI)`  = "HR (95% CI)",
    `p-value`      = "p-value"
  ) %>%
  autofit() %>%
  theme_booktabs()

doc <- read_docx() %>%
  body_add_par("Supplementary Table. Univariable Cox regression by molecular subtype — GSE85217",
               style = "heading 2") %>%
  body_add_flextable(ft) %>%
  body_add_par(paste0("G3/4 overall model: likelihood ratio p = ",
                      sprintf("%.3f", summary(cox_subtype_g34)$logtest["pvalue"]),
                      ". SHH overall model: likelihood ratio p = ",
                      sprintf("%.3f", summary(cox_subtype_shh)$logtest["pvalue"]),
                      ". Bold p-values indicate statistical significance (p<0.05)."),
               style = "Normal")

print(doc, target = "Supplementary_Table_Subtype_Cox.docx")
cat("Saved: Supplementary_Table_Subtype_Cox.docx\n")


###########################################################################
# Supplementary Figure — PCA and UMAP coloured by molecular subtype
# Group 3, Group 4 and SHH separately (3x2 panel)
# Input: expr_top_sub and meta_sub from MB_main.R
###########################################################################

library(uwot)
library(ggplot2)
library(ggpubr)

# ── Shared subtype colour palettes ────────────────────────────────────────────
g3_subtype_cols <- c(
  "Group3_alpha" = "#FFD700",
  "Group3_beta"  = "#FFA500",
  "Group3_gamma" = "#FF6600"
)
g3_subtype_labels <- c(
  "Group3_alpha" = "Group 3-\u03b1",
  "Group3_beta"  = "Group 3-\u03b2",
  "Group3_gamma" = "Group 3-\u03b3"
)

g4_subtype_cols <- c(
  "Group4_alpha" = "#006400",
  "Group4_beta"  = "#228B22",
  "Group4_gamma" = "#90EE90"
)
g4_subtype_labels <- c(
  "Group4_alpha" = "Group 4-\u03b1",
  "Group4_beta"  = "Group 4-\u03b2",
  "Group4_gamma" = "Group 4-\u03b3"
)

shh_subtype_cols <- c(
  "SHH_alpha" = "#7b2d8b",
  "SHH_beta"  = "#e7298a",
  "SHH_gamma" = "#66a61e",
  "SHH_delta" = "#e6ab02"
)
shh_subtype_labels <- c(
  "SHH_alpha" = "SHH-\u03b1",
  "SHH_beta"  = "SHH-\u03b2",
  "SHH_gamma" = "SHH-\u03b3",
  "SHH_delta" = "SHH-\u03b4"
)

# ── Helper: build expression matrix and aligned metadata ──────────────────────
build_expr_meta <- function(subgroup_filter) {
  meta_filt <- meta_sub[meta_sub$Subgroup %in% subgroup_filter &
                          !is.na(meta_sub$Subtype) &
                          meta_sub$Subtype != "", ]
  common    <- intersect(meta_filt$ExpressionID, colnames(expr_top_sub))
  expr_filt <- expr_top_sub[, common, drop = FALSE]
  meta_filt <- meta_filt[match(common, meta_filt$ExpressionID), ]
  stopifnot(all(colnames(expr_filt) == meta_filt$ExpressionID))
  list(expr = expr_filt, meta = meta_filt)
}

# ── Helper: PCA plot ──────────────────────────────────────────────────────────
make_pca_plot <- function(expr_mat, meta_df, subtype_labels, subtype_cols, title) {
  set.seed(123)
  pca     <- prcomp(t(expr_mat), scale. = TRUE)
  var_exp <- round(summary(pca)$importance[2, 1:2] * 100, 1)
  
  df <- data.frame(
    PC1     = pca$x[, 1],
    PC2     = pca$x[, 2],
    Subtype = dplyr::recode(meta_df$Subtype, !!!subtype_labels)
  ) %>% filter(!is.na(Subtype))
  df$Subtype <- factor(df$Subtype, levels = subtype_labels)
  
  ggplot(df, aes(PC1, PC2, colour = Subtype)) +
    geom_point(size = 1.8, alpha = 0.8) +
    scale_colour_manual(values = setNames(subtype_cols, subtype_labels)) +
    labs(x = paste0("PC1 (", var_exp[1], "%)"),
         y = paste0("PC2 (", var_exp[2], "%)"),
         title = title, colour = "Subtype") +
    theme_classic(base_size = 11) +
    theme(legend.key.size = unit(0.4, "cm"),
          legend.text = element_text(size = 8))
}

# ── Helper: UMAP plot ─────────────────────────────────────────────────────────
make_umap_plot <- function(expr_mat, meta_df, subtype_labels, subtype_cols, title) {
  set.seed(123)
  umap_res <- umap(t(expr_mat), n_neighbors = 15, min_dist = 0.3, metric = "euclidean")
  
  df <- data.frame(
    UMAP1   = umap_res[, 1],
    UMAP2   = umap_res[, 2],
    Subtype = dplyr::recode(meta_df$Subtype, !!!subtype_labels)
  ) %>% filter(!is.na(Subtype))
  df$Subtype <- factor(df$Subtype, levels = subtype_labels)
  
  ggplot(df, aes(UMAP1, UMAP2, colour = Subtype)) +
    geom_point(size = 1.8, alpha = 0.8) +
    scale_colour_manual(values = setNames(subtype_cols, subtype_labels)) +
    labs(x = "UMAP1", y = "UMAP2",
         title = title, colour = "Subtype") +
    theme_classic(base_size = 11) +
    theme(legend.key.size = unit(0.4, "cm"),
          legend.text = element_text(size = 8))
}

# ── 1. Group 3 ────────────────────────────────────────────────────────────────
g3 <- build_expr_meta("Group3")
cat("Group 3 samples:", ncol(g3$expr), "\n")

A <- make_pca_plot(g3$expr, g3$meta, g3_subtype_labels, g3_subtype_cols,
                   "Group 3: PCA by molecular subtype")
B <- make_umap_plot(g3$expr, g3$meta, g3_subtype_labels, g3_subtype_cols,
                    "Group 3: UMAP by molecular subtype")

# ── 2. Group 4 ────────────────────────────────────────────────────────────────
g4 <- build_expr_meta("Group4")
cat("Group 4 samples:", ncol(g4$expr), "\n")

C <- make_pca_plot(g4$expr, g4$meta, g4_subtype_labels, g4_subtype_cols,
                   "Group 4: PCA by molecular subtype")
D <- make_umap_plot(g4$expr, g4$meta, g4_subtype_labels, g4_subtype_cols,
                    "Group 4: UMAP by molecular subtype")

# ── 3. SHH ───────────────────────────────────────────────────────────────────
shh <- build_expr_meta("SHH")
cat("SHH samples:", ncol(shh$expr), "\n")

E <- make_pca_plot(shh$expr, shh$meta, shh_subtype_labels, shh_subtype_cols,
                   "SHH: PCA by molecular subtype")
F_plot <- make_umap_plot(shh$expr, shh$meta, shh_subtype_labels, shh_subtype_cols,
                         "SHH: UMAP by molecular subtype")

# ── 4. Combined PDF ───────────────────────────────────────────────────────────
quartz(type = "pdf", file = "Supp_PCA_UMAP_subtypes.pdf", width = 14, height = 14)

print(ggarrange(A, B, C, D, E, F_plot,
                ncol = 2, nrow = 3,
                labels = c("A", "B", "C", "D", "E", "F")))

dev.off()
cat("Saved: Supp_PCA_UMAP_subtypes.pdf\n")






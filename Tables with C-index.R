################################################################################
#
# MEDULLOBLASTOMA PROJECT — TABLES
# Dataset: GSE85217 + E-MTAB-10767
#
# Student: Ewan McGibbon
# Date:    20/2/2026
#
################################################################################

# =============================================================================
# 0. SETUP
# =============================================================================

library(flextable)
library(officer)
library(dplyr)
library(tidyverse)

########################################################################################################################
#                                   SHARED HELPER FUNCTIONS
########################################################################################################################

fmt_p <- function(p) {
  ifelse(p < 0.001, "<0.001", formatC(p, digits = 3, format = "f"))
}

fmt_hr_ci <- function(hr, lo, hi) {
  paste0(formatC(hr, digits = 2, format = "f"),
         " (", formatC(lo, digits = 2, format = "f"),
         "\u2013", formatC(hi, digits = 2, format = "f"), ")")
}

fmt_n_pct <- function(vec, level) {
  denom <- sum(!is.na(vec))
  n     <- sum(vec == level, na.rm = TRUE)
  if (denom == 0) return("0 (0%)")
  paste0(n, " (", round(100 * n / denom), "%)")
}

fmt_median_range <- function(x, digits = 1) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("\u2014")
  paste0(round(median(x), digits), " (",
         round(min(x), digits), "\u2013",
         round(max(x), digits), ")")
}

miss <- function(x) sum(is.na(x))

simple_ft <- function(df, col_header, header_rows = NULL,
                      sig_rows = NULL, italic_rows = NULL) {
  ft <- flextable::flextable(df)
  ft <- flextable::set_header_labels(ft, values = setNames(col_header, names(df)))
  n_cols     <- ncol(df)
  col_widths <- if (n_cols == 2) c(3.8, 2.0) else c(3.2, 1.8, 1.0)
  thin_line  <- officer::fp_border(color = "#BBBBBB", width = 0.5)
  thick_line <- officer::fp_border(color = "#000000", width = 1.0)
  ft <- ft %>%
    flextable::border_remove() %>%
    flextable::hline_top(border = thick_line, part = "header") %>%
    flextable::hline_bottom(border = thick_line, part = "header") %>%
    flextable::border(border.bottom = thin_line, part = "body") %>%
    flextable::hline_bottom(border = thick_line, part = "body") %>%
    flextable::align(j = 1, align = "left",  part = "all") %>%
    flextable::align(j = seq(2, n_cols), align = "right", part = "all") %>%
    flextable::bold(part = "header") %>%
    flextable::italic(part = "header") %>%
    flextable::fontsize(size = 11, part = "all") %>%
    flextable::font(fontname = "Times New Roman", part = "all") %>%
    flextable::padding(padding.top = 4, padding.bottom = 4,
                       padding.left = 6, padding.right = 6, part = "all")
  all_rows    <- seq_len(nrow(df))
  sub_rows    <- if (!is.null(header_rows)) setdiff(all_rows, header_rows) else all_rows
  nodata_rows <- grep("No data", df[[1]])
  sub_rows    <- setdiff(sub_rows, nodata_rows)
  if (length(sub_rows) > 0) {
    ft <- flextable::padding(ft, i = sub_rows, j = 1,
                             padding.left = 20, part = "body")
  }
  if (!is.null(header_rows)) {
    ft <- flextable::padding(ft, i = header_rows, j = 1,
                             padding.top = 10, padding.left = 6, part = "body")
  }
  for (j in seq_len(n_cols)) {
    ft <- flextable::width(ft, j = j, width = col_widths[j])
  }
  if (!is.null(header_rows)) {
    ft <- flextable::bold(ft, i = header_rows, part = "body")
  }
  if (!is.null(italic_rows)) {
    ft <- flextable::italic(ft, i = italic_rows, part = "body")
  }
  if (!is.null(sig_rows)) {
    ft <- flextable::bold(ft, i = sig_rows, j = ncol(df), part = "body")
  }
  ft
}

save_table <- function(title, ft, footnote, filename) {
  officer::read_docx() %>%
    officer::body_add_par(title, style = "Normal") %>%
    flextable::body_add_flextable(ft) %>%
    officer::body_add_par(footnote, style = "Normal") %>%
    print(target = filename)
  message("Saved: ", filename)
}

########################################################################################################################
#                                   PREPARE DATA
########################################################################################################################

cavalli <- meta %>%
  transmute(
    sex      = case_when(Gender == "M" ~ "Male", Gender == "F" ~ "Female",
                         TRUE ~ NA_character_),
    age      = as.numeric(Age),
    subgroup = case_when(Subgroup == "WNT" ~ "WNT", Subgroup == "SHH" ~ "SHH",
                         Subgroup == "Group3" ~ "Group 3",
                         Subgroup == "Group4" ~ "Group 4",
                         TRUE ~ NA_character_),
    mol_subtype = case_when(
      Subtype == "WNT_alpha" ~ "WNT-alpha", Subtype == "WNT_beta" ~ "WNT-beta",
      Subtype == "SHH_alpha" ~ "SHH-alpha", Subtype == "SHH_beta" ~ "SHH-beta",
      Subtype == "SHH_gamma" ~ "SHH-gamma", Subtype == "SHH_delta" ~ "SHH-delta",
      Subtype == "Group3_alpha" ~ "Group 3-alpha",
      Subtype == "Group3_beta"  ~ "Group 3-beta",
      Subtype == "Group3_gamma" ~ "Group 3-gamma",
      Subtype == "Group4_alpha" ~ "Group 4-alpha",
      Subtype == "Group4_beta"  ~ "Group 4-beta",
      Subtype == "Group4_gamma" ~ "Group 4-gamma",
      TRUE ~ NA_character_),
    histology = case_when(
      meta$histology == "Classic"      ~ "Classic",
      meta$histology == "LCA"          ~ "LCA",
      meta$histology == "Desmoplastic" ~ "Desmoplastic",
      meta$histology == "MBEN"         ~ "MBEN",
      TRUE                             ~ NA_character_),
    met      = case_when(Met.status..1.Met..0.M0. == 1 ~ "M+",
                         Met.status..1.Met..0.M0. == 0 ~ "M0",
                         TRUE ~ NA_character_),
    survival = case_when(Dead == 1 ~ "Dead", Dead == 0 ~ "Alive",
                         TRUE ~ NA_character_),
    followup = as.numeric(OS..years.),
    MYC  = ifelse(Study_ID %in% myc_amp_ids,  "Present", "Absent"),
    MYCN = ifelse(Study_ID %in% mycn_amp_ids, "Present", "Absent")
  )

rnaseq <- sdrf %>%
  transmute(
    sex = case_when(
      `Characteristics[sex]` == "male"   ~ "Male",
      `Characteristics[sex]` == "female" ~ "Female",
      TRUE ~ NA_character_),
    age_band = case_when(
      `Characteristics[age]` == "0 to 3"              ~ "0\u20133 years",
      `Characteristics[age]` == "3 to 16"             ~ "3\u201316 years",
      `Characteristics[age]` %in% c("> 16","over 16") ~ ">16 years",
      TRUE ~ NA_character_),
    subgroup = case_when(
      `Characteristics[subgroup]` == "WNT"       ~ "WNT",
      `Characteristics[subgroup]` == "SHH"       ~ "SHH",
      `Characteristics[subgroup]` == "Grp3"      ~ "Group 3",
      `Characteristics[subgroup]` == "Grp4"      ~ "Group 4",
      `Characteristics[subgroup]` == "Grp3/Grp4" ~ "Group 3/4",
      `Characteristics[subgroup]` == "MB-NOS"    ~ "MB-NOS",
      TRUE ~ NA_character_),
    g34_subtype = case_when(
      `Characteristics[grp3/4 subtype]` %in%
        c("I","II","III","IV","V","VI","VII","VIII") ~
        `Characteristics[grp3/4 subtype]`,
      TRUE ~ NA_character_),
    g34_category = case_when(
      `Characteristics[g3g4 category]` == "High.G3" ~ "High G3",
      `Characteristics[g3g4 category]` == "Low.G3"  ~ "Low G3",
      `Characteristics[g3g4 category]` == "G3.5"    ~ "G3.5",
      `Characteristics[g3g4 category]` == "High.G4" ~ "High G4",
      `Characteristics[g3g4 category]` == "Low.G4"  ~ "Low G4",
      TRUE ~ NA_character_),
    histology = case_when(
      `Characteristics[large cell/anaplastic]` == TRUE  ~ "LCA",
      `Characteristics[desmoplastic/nodular]`  == TRUE  ~ "DN",
      `Characteristics[large cell/anaplastic]` == FALSE &
        `Characteristics[desmoplastic/nodular]` == FALSE ~ "Classic",
      TRUE ~ NA_character_),
    met      = case_when(`Characteristics[m+]` == TRUE  ~ "M+",
                         `Characteristics[m+]` == FALSE ~ "M0",
                         TRUE ~ NA_character_),
    resection = case_when(
      `Characteristics[sub-total resection]` == TRUE  ~ "STR",
      `Characteristics[sub-total resection]` == FALSE ~ "GTR",
      TRUE ~ NA_character_),
    survival = case_when(
      `Characteristics[os status]` == 1 ~ "Dead",
      `Characteristics[os status]` == 0 ~ "Alive",
      TRUE ~ NA_character_),
    followup = suppressWarnings(as.numeric(`Characteristics[os time]`)),
    MYC  = case_when(`Characteristics[myc amplification]`  == TRUE  ~ "Present",
                     `Characteristics[myc amplification]`  == FALSE ~ "Absent",
                     TRUE ~ NA_character_),
    MYCN = case_when(`Characteristics[mycn amplification]` == TRUE  ~ "Present",
                     `Characteristics[mycn amplification]` == FALSE ~ "Absent",
                     TRUE ~ NA_character_),
    TP53 = case_when(`Characteristics[tp53 mutation]`      == TRUE  ~ "Present",
                     `Characteristics[tp53 mutation]`      == FALSE ~ "Absent",
                     TRUE ~ NA_character_),
    CTNNB1 = case_when(`Characteristics[ctnnb1 mutation]`  == TRUE  ~ "Present",
                       `Characteristics[ctnnb1 mutation]`  == FALSE ~ "Absent",
                       TRUE ~ NA_character_),
    TERT = case_when(`Characteristics[tert mutation]`      == TRUE  ~ "Present",
                     `Characteristics[tert mutation]`      == FALSE ~ "Absent",
                     TRUE ~ NA_character_)
  )

rnaseq$followup[rnaseq$followup > 50] <- NA

cat("Cavalli histology:\n")
print(table(cavalli$histology, useNA = "always"))

########################################################################################################################
#                                   TABLE 1 — CAVALLI BASELINE
########################################################################################################################

age_grp <- case_when(is.na(cavalli$age) ~ NA_character_,
                     cavalli$age < 3 ~ "<3 years", TRUE ~ "\u22653 years")

t1 <- tribble(
  ~Characteristic,                              ~`Cavalli microarray (n = 763)`,
  "Sex",                                         "",
  "  Male",                                      fmt_n_pct(cavalli$sex, "Male"),
  "  Female",                                    fmt_n_pct(cavalli$sex, "Female"),
  paste0("  No data (", miss(cavalli$sex), ")"), "",
  "Age at diagnosis",                            "",
  "  Median, years (range)",                     fmt_median_range(cavalli$age),
  "  <3 years",                                  fmt_n_pct(age_grp, "<3 years"),
  "  \u22653 years",                             fmt_n_pct(age_grp, "\u22653 years"),
  paste0("  No data (", miss(cavalli$age), ")"), "",
  "Molecular subgroup",                          "",
  "  WNT",                                       fmt_n_pct(cavalli$subgroup, "WNT"),
  "  SHH",                                       fmt_n_pct(cavalli$subgroup, "SHH"),
  "  Group 3",                                   fmt_n_pct(cavalli$subgroup, "Group 3"),
  "  Group 4",                                   fmt_n_pct(cavalli$subgroup, "Group 4"),
  "Molecular subtype",                           "",
  "  WNT-alpha",                                 fmt_n_pct(cavalli$mol_subtype, "WNT-alpha"),
  "  WNT-beta",                                  fmt_n_pct(cavalli$mol_subtype, "WNT-beta"),
  "  SHH-alpha",                                 fmt_n_pct(cavalli$mol_subtype, "SHH-alpha"),
  "  SHH-beta",                                  fmt_n_pct(cavalli$mol_subtype, "SHH-beta"),
  "  SHH-gamma",                                 fmt_n_pct(cavalli$mol_subtype, "SHH-gamma"),
  "  SHH-delta",                                 fmt_n_pct(cavalli$mol_subtype, "SHH-delta"),
  "  Group 3-alpha",                             fmt_n_pct(cavalli$mol_subtype, "Group 3-alpha"),
  "  Group 3-beta",                              fmt_n_pct(cavalli$mol_subtype, "Group 3-beta"),
  "  Group 3-gamma",                             fmt_n_pct(cavalli$mol_subtype, "Group 3-gamma"),
  "  Group 4-alpha",                             fmt_n_pct(cavalli$mol_subtype, "Group 4-alpha"),
  "  Group 4-beta",                              fmt_n_pct(cavalli$mol_subtype, "Group 4-beta"),
  "  Group 4-gamma",                             fmt_n_pct(cavalli$mol_subtype, "Group 4-gamma"),
  "Histopathological variant",                   "",
  "  Classic",                                   fmt_n_pct(cavalli$histology, "Classic"),
  "  LCA",                                       fmt_n_pct(cavalli$histology, "LCA"),
  "  Desmoplastic",                              fmt_n_pct(cavalli$histology, "Desmoplastic"),
  "  MBEN",                                      fmt_n_pct(cavalli$histology, "MBEN"),
  paste0("  No data (", miss(cavalli$histology), ")"), "",
  "Metastatic stage",                            "",
  "  M0",                                        fmt_n_pct(cavalli$met, "M0"),
  "  M+",                                        fmt_n_pct(cavalli$met, "M+"),
  paste0("  No data (", miss(cavalli$met), ")"), "",
  "Survival status",                             "",
  "  Alive",                                     fmt_n_pct(cavalli$survival, "Alive"),
  "  Dead",                                      fmt_n_pct(cavalli$survival, "Dead"),
  paste0("  No data (", miss(cavalli$survival), ")"), "",
  "Follow-up",                                   "",
  "  Median, years (range)",                     fmt_median_range(cavalli$followup),
  paste0("  No data (", miss(cavalli$followup), ")"), "",
  "MYC amplification",                           "",
  "  Present",                                   fmt_n_pct(cavalli$MYC,  "Present"),
  "  Absent",                                    fmt_n_pct(cavalli$MYC,  "Absent"),
  "MYCN amplification",                          "",
  "  Present",                                   fmt_n_pct(cavalli$MYCN, "Present"),
  "  Absent",                                    fmt_n_pct(cavalli$MYCN, "Absent")
)

t1_hdrs <- which(t1$`Cavalli microarray (n = 763)` == "" &
                   !grepl("No data", t1$Characteristic))

ft1 <- simple_ft(t1,
                 col_header  = c("Characteristic", "Cavalli microarray (n = 763)"),
                 header_rows = t1_hdrs)

########################################################################################################################
#                                   TABLE 2 — RNA-SEQ BASELINE
########################################################################################################################

t2 <- tribble(
  ~Characteristic,                                   ~`E-MTAB-10767 RNA-seq (n = 331)`,
  "Sex",                                              "",
  "  Male",                                           fmt_n_pct(rnaseq$sex, "Male"),
  "  Female",                                         fmt_n_pct(rnaseq$sex, "Female"),
  paste0("  No data (", miss(rnaseq$sex), ")"),       "",
  "Age at diagnosis (banded)",                        "",
  "  0\u20133 years",                                 fmt_n_pct(rnaseq$age_band, "0\u20133 years"),
  "  3\u201316 years",                                fmt_n_pct(rnaseq$age_band, "3\u201316 years"),
  "  >16 years",                                      fmt_n_pct(rnaseq$age_band, ">16 years"),
  paste0("  No data (", miss(rnaseq$age_band), ")"),  "",
  "Molecular subgroup",                               "",
  "  WNT",                                            fmt_n_pct(rnaseq$subgroup, "WNT"),
  "  SHH",                                            fmt_n_pct(rnaseq$subgroup, "SHH"),
  "  Group 3",                                        fmt_n_pct(rnaseq$subgroup, "Group 3"),
  "  Group 4",                                        fmt_n_pct(rnaseq$subgroup, "Group 4"),
  "  Group 3/4",                                      fmt_n_pct(rnaseq$subgroup, "Group 3/4"),
  "  MB-NOS",                                         fmt_n_pct(rnaseq$subgroup, "MB-NOS"),
  "Group 3/4 subtype",                                "",
  "  I",    fmt_n_pct(rnaseq$g34_subtype, "I"),
  "  II",   fmt_n_pct(rnaseq$g34_subtype, "II"),
  "  III",  fmt_n_pct(rnaseq$g34_subtype, "III"),
  "  IV",   fmt_n_pct(rnaseq$g34_subtype, "IV"),
  "  V",    fmt_n_pct(rnaseq$g34_subtype, "V"),
  "  VI",   fmt_n_pct(rnaseq$g34_subtype, "VI"),
  "  VII",  fmt_n_pct(rnaseq$g34_subtype, "VII"),
  "  VIII", fmt_n_pct(rnaseq$g34_subtype, "VIII"),
  paste0("  No data (", miss(rnaseq$g34_subtype), ")"),   "",
  "Group 3/4 risk category",                          "",
  "  High G3",  fmt_n_pct(rnaseq$g34_category, "High G3"),
  "  Low G3",   fmt_n_pct(rnaseq$g34_category, "Low G3"),
  "  G3.5",     fmt_n_pct(rnaseq$g34_category, "G3.5"),
  "  High G4",  fmt_n_pct(rnaseq$g34_category, "High G4"),
  "  Low G4",   fmt_n_pct(rnaseq$g34_category, "Low G4"),
  paste0("  No data (", miss(rnaseq$g34_category), ")"),  "",
  "Histopathological variant",                        "",
  "  Classic",  fmt_n_pct(rnaseq$histology, "Classic"),
  "  LCA",      fmt_n_pct(rnaseq$histology, "LCA"),
  "  DN",       fmt_n_pct(rnaseq$histology, "DN"),
  paste0("  No data (", miss(rnaseq$histology), ")"), "",
  "Metastatic stage",                                 "",
  "  M0",       fmt_n_pct(rnaseq$met, "M0"),
  "  M+",       fmt_n_pct(rnaseq$met, "M+"),
  paste0("  No data (", miss(rnaseq$met), ")"),       "",
  "Resection",                                        "",
  "  GTR",      fmt_n_pct(rnaseq$resection, "GTR"),
  "  STR",      fmt_n_pct(rnaseq$resection, "STR"),
  paste0("  No data (", miss(rnaseq$resection), ")"), "",
  "Survival status",                                  "",
  "  Alive",    fmt_n_pct(rnaseq$survival, "Alive"),
  "  Dead",     fmt_n_pct(rnaseq$survival, "Dead"),
  paste0("  No data (", miss(rnaseq$survival), ")"),  "",
  "Follow-up",                                        "",
  "  Median, years (range)", fmt_median_range(rnaseq$followup),
  paste0("  No data (", miss(rnaseq$followup), ")"),  "",
  "MYC amplification",                                "",
  "  Present",  fmt_n_pct(rnaseq$MYC,  "Present"),
  "  Absent",   fmt_n_pct(rnaseq$MYC,  "Absent"),
  "MYCN amplification",                               "",
  "  Present",  fmt_n_pct(rnaseq$MYCN, "Present"),
  "  Absent",   fmt_n_pct(rnaseq$MYCN, "Absent"),
  "TP53 mutation",                                    "",
  "  Present",  fmt_n_pct(rnaseq$TP53, "Present"),
  "  Absent",   fmt_n_pct(rnaseq$TP53, "Absent"),
  "CTNNB1 mutation",                                  "",
  "  Present",  fmt_n_pct(rnaseq$CTNNB1, "Present"),
  "  Absent",   fmt_n_pct(rnaseq$CTNNB1, "Absent"),
  "TERT mutation",                                    "",
  "  Present",  fmt_n_pct(rnaseq$TERT, "Present"),
  "  Absent",   fmt_n_pct(rnaseq$TERT, "Absent")
)

t2_hdrs <- which(t2$`E-MTAB-10767 RNA-seq (n = 331)` == "" &
                   !grepl("No data", t2$Characteristic))

ft2 <- simple_ft(t2,
                 col_header  = c("Characteristic",
                                 "E-MTAB-10767 RNA-seq (n = 331)"),
                 header_rows = t2_hdrs)

########################################################################################################################
#                                   TABLE 3 — G3/4 UNIVARIABLE COX
########################################################################################################################

# Refit histology model with Classic as reference
meta_g34$histology <- relevel(factor(as.character(meta_g34$histology)), ref = "Classic")
cox_hist_g34 <- coxph(Surv(OS..years., Dead) ~ histology, data = meta_g34)

# Extract coefficients
age_g34   <- extract_cox(cox_age_g34,      "AgeGroup0-3")
age10_g34 <- extract_cox(cox_age_g34,      "AgeGroup10+")
hd_g34    <- extract_cox(cox_hist_g34,     "histologyDesmoplastic")
hl_g34    <- extract_cox(cox_hist_g34,     "histologyLCA")
hm_g34    <- extract_cox(cox_hist_g34,     "histologyMBEN")
met_g34   <- extract_cox(cox_met_g34,      "Met_statusMet")
sg_g34    <- extract_cox(cox_subgroup_g34, "SubgroupGroup4")
myc_g34   <- extract_cox(cox_myc_g34,      "MYCAmplified")

t3 <- tribble(
  ~Characteristic,                                             ~`HR (95% CI)`, ~`p-value`,
  "Age at diagnosis",                                           "",             "",
  "  4-10 years",                                               "Reference",    "\u2014",
  "  <3 years",                                                 age_g34$hr,     age_g34$p,
  "  \u226510 years",                                           age10_g34$hr,   age10_g34$p,
  paste0("  No data (", nrow(meta_g34) - cox_age_g34$n, ")"),  "",             "",
  "Histopathological variant",                                  "",             "",
  "  Classic",                                                  "Reference",    "\u2014",
  "  Desmoplastic",                                             hd_g34$hr,      hd_g34$p,
  "  LCA",                                                      hl_g34$hr,      hl_g34$p,
  "  MBEN",                                                     hm_g34$hr,      hm_g34$p,
  paste0("  No data (", nrow(meta_g34) - cox_hist_g34$n, ")"), "",             "",
  "Metastatic stage",                                           "",             "",
  "  M0",                                                       "Reference",    "\u2014",
  "  M+",                                                       met_g34$hr,     met_g34$p,
  paste0("  No data (", nrow(meta_g34) - cox_met_g34$n, ")"),  "",             "",
  "Intra-subgroup",                                             "",             "",
  "  Group 3",                                                  "Reference",    "\u2014",
  "  Group 4",                                                  sg_g34$hr,      sg_g34$p,
  "MYC amplification",                                          "",             "",
  "  Not amplified",                                            "Reference",    "\u2014",
  "  Amplified",                                                myc_g34$hr,     myc_g34$p
)

t3_hdrs <- which(t3$`HR (95% CI)` == "" & !grepl("No data", t3$Characteristic))
t3_refs  <- which(t3$`HR (95% CI)` == "Reference")
t3_sigs  <- which(suppressWarnings(as.numeric(gsub("<", "0", t3$`p-value`))) < 0.05)

ft3 <- simple_ft(t3,
                 col_header  = c(paste0("Characteristic\n(n = ",
                                        nrow(meta_g34), ", events = ",
                                        sum(meta_g34$Dead, na.rm = TRUE), ")"),
                                 "HR (95% CI)", "p-value"),
                 header_rows = t3_hdrs,
                 italic_rows = t3_refs,
                 sig_rows    = t3_sigs)

save_table(
  title    = "Table 3. Univariable Cox regression of clinical variables \u2014 Group 3/4 medulloblastoma (GSE85217).",
  ft       = ft3,
  footnote = "HR, hazard ratio; CI, confidence interval. Reference categories in italics. Bold p < 0.05. Age reference category: 4-10 years. Histology reference category: Classic.",
  filename = "Table3_G34_Cox.docx"
)
########################################################################################################################
#                                   TABLE 4 — SHH UNIVARIABLE COX
########################################################################################################################

if (!exists("cox_mycn_shh")) {
  cox_mycn_shh <- coxph(Surv(OS..years., Dead) ~ MYCN, data = meta_shh)
}

if (!exists("cox_hist_shh") ||
    cox_hist_shh$n == nrow(meta_shh[!is.na(meta_shh$OS..years.) & !is.na(meta_shh$Dead), ])) {
  cox_hist_shh <- coxph(Surv(OS..years., Dead) ~ histology,
                        data = meta_shh[meta_shh$histology != "MBEN" &
                                          !is.na(meta_shh$histology), ])
}

age_shh   <- extract_cox(cox_age_shh,  "AgeGroup0-3")
age10_shh <- extract_cox(cox_age_shh,  "AgeGroup10+")
hd_shh    <- extract_cox(cox_hist_shh, "histologyDesmoplastic")
hl_shh    <- extract_cox(cox_hist_shh, "histologyLCA")
met_shh   <- extract_cox(cox_met_shh,  "Met_statusMet")
mycn_shh  <- extract_cox(cox_mycn_shh, "MYCNAmplified")

t4 <- tribble(
  ~Characteristic,                                             ~`HR (95% CI)`,  ~`p-value`,
  "Age at diagnosis",                                           "",               "",
  "  4-10 years",                                               "Reference",      "\u2014",
  "  <3 years",                                                 age_shh$hr,       age_shh$p,
  "  \u226510 years",                                           age10_shh$hr,     age10_shh$p,
  paste0("  No data (", nrow(meta_shh) - cox_age_shh$n, ")"),  "",               "",
  "Histopathological variant",                                  "",               "",
  "  Classic",                                                  "Reference",      "\u2014",
  "  Desmoplastic",                                             hd_shh$hr,        hd_shh$p,
  "  LCA",                                                      hl_shh$hr,        hl_shh$p,
  paste0("  No data (", nrow(meta_shh) - cox_hist_shh$n, ")"), "",               "",
  "Metastatic stage",                                           "",               "",
  "  M0",                                                       "Reference",      "\u2014",
  "  M+",                                                       met_shh$hr,       met_shh$p,
  paste0("  No data (", nrow(meta_shh) - cox_met_shh$n, ")"),  "",               "",
  "MYCN amplification",                                         "",               "",
  "  Not amplified",                                            "Reference",      "\u2014",
  "  Amplified",                                                mycn_shh$hr,      mycn_shh$p
)

t4_hdrs <- which(t4$`HR (95% CI)` == "" & !grepl("No data", t4$Characteristic))
t4_refs  <- which(t4$`HR (95% CI)` == "Reference")
t4_sigs  <- which(suppressWarnings(as.numeric(gsub("<", "0", t4$`p-value`))) < 0.05)

ft4 <- simple_ft(t4,
                 col_header  = c(paste0("Characteristic\n(n = ",
                                        nrow(meta_shh), ", events = ",
                                        sum(meta_shh$Dead, na.rm = TRUE), ")"),
                                 "HR (95% CI)", "p-value"),
                 header_rows = t4_hdrs,
                 italic_rows = t4_refs,
                 sig_rows    = t4_sigs)

########################################################################################################################
#                                   TABLE 5 — MODEL PERFORMANCE
########################################################################################################################

conc_g34 <- concordance(cox_model_g34)
conc_shh <- concordance(cox_model_shh)
conc_rna <- concordance(Surv(time, event) ~ risk_lp, data = risk_out_rnaseq)
km_rna   <- survdiff(Surv(time, event) ~ tertile, data = risk_out_rnaseq)

t5 <- tribble(
  ~Cohort,                         ~Subgroup,                ~N,
  ~Events,                         ~`C-index (SE)`,          ~`Log-rank p`,
  
  "Cavalli microarray (GSE85217)", "Group 3/4",
  as.character(cox_model_g34$n),   as.character(cox_model_g34$nevent),
  paste0(formatC(conc_g34$concordance, digits = 3, format = "f"),
         " (", formatC(sqrt(conc_g34$var), digits = 3, format = "f"), ")"),
  fmt_p(pchisq(lr_test_g34$chisq, df = length(lr_test_g34$n) - 1, lower.tail = FALSE)),
  
  "E-MTAB-10767 RNA-seq",          "Group 3/4 (validation)",
  as.character(nrow(risk_out_rnaseq)), as.character(sum(risk_out_rnaseq$event)),
  paste0(formatC(1 - conc_rna$concordance, digits = 3, format = "f"),
         " (", formatC(sqrt(conc_rna$var), digits = 3, format = "f"), ")\u2020"),
  fmt_p(pchisq(km_rna$chisq, df = 2, lower.tail = FALSE)),
  
  "Cavalli microarray (GSE85217)", "SHH",
  as.character(cox_model_shh$n),   as.character(cox_model_shh$nevent),
  paste0(formatC(conc_shh$concordance, digits = 3, format = "f"),
         " (", formatC(sqrt(conc_shh$var), digits = 3, format = "f"), ")"),
  fmt_p(pchisq(lr_test_shh$chisq, df = length(lr_test_shh$n) - 1, lower.tail = FALSE)),
  
  "E-MTAB-10767 RNA-seq",          "SHH (validation)",
  as.character(nrow(risk_out_shh_rna)), "No OS data",
  "N/A",                            "N/A"
)

thin_line  <- officer::fp_border(color = "#BBBBBB", width = 0.5)
thick_line <- officer::fp_border(color = "#000000", width = 1.0)

ft5 <- flextable::flextable(t5) %>%
  flextable::bold(part = "header") %>%
  flextable::italic(part = "header") %>%
  flextable::border_remove() %>%
  flextable::hline_top(border = thick_line, part = "header") %>%
  flextable::hline_bottom(border = thick_line, part = "header") %>%
  flextable::border(border.bottom = thin_line, part = "body") %>%
  flextable::hline_bottom(border = thick_line, part = "body") %>%
  flextable::hline(i = 2, border = officer::fp_border(color = "#555555", width = 0.8), part = "body") %>%
  flextable::align(j = 1, align = "left",  part = "all") %>%
  flextable::align(j = 2, align = "left",  part = "all") %>%
  flextable::align(j = seq(3, 6), align = "right", part = "all") %>%
  flextable::color(i = 4, color = "#888888", part = "body") %>%
  flextable::fontsize(size = 11, part = "all") %>%
  flextable::font(fontname = "Times New Roman", part = "all") %>%
  flextable::padding(padding.top = 4, padding.bottom = 4,
                     padding.left = 6, padding.right = 6, part = "all") %>%
  flextable::width(j = 1, width = 2.0) %>%
  flextable::width(j = 2, width = 2.0) %>%
  flextable::width(j = 3, width = 0.4) %>%
  flextable::width(j = 4, width = 0.7) %>%
  flextable::width(j = 5, width = 1.4) %>%
  flextable::width(j = 6, width = 0.9)

########################################################################################################################
#                                   EXPORT
########################################################################################################################

save_table(
  title    = "Table 1. Clinical and molecular characteristics of the Cavalli microarray cohort (GSE85217).",
  ft       = ft1,
  footnote = "Abbreviations: LCA, large cell/anaplastic; DN, desmoplastic/nodular; MBEN, medulloblastoma with extensive nodularity; M0, non-metastatic; M+, metastatic. Categorical variables: n (% of non-missing). Continuous variables: median (range). Missing data counts shown in parentheses. MYC and MYCN amplification status derived from Cavalli et al. (Nature, 2017) supplementary annotation.",
  filename = "Table1_Cavalli_baseline.docx"
)

save_table(
  title    = "Table 2. Clinical and molecular characteristics of the E-MTAB-10767 RNA-seq cohort.",
  ft       = ft2,
  footnote = "Abbreviations: LCA, large cell/anaplastic; DN, desmoplastic/nodular; GTR, gross total resection; STR, sub-total resection; M0, non-metastatic; M+, metastatic. Categorical variables: n (% of non-missing). Age reported as banded categories. Histology inferred from binary SDRF flags. Group 3/4 subtypes and risk categories per Cavalli et al. (Nature, 2017).",
  filename = "Table2_RNAseq_baseline.docx"
)

save_table(
  title    = "Table 3. Univariable Cox regression of clinical variables \u2014 Group 3/4 medulloblastoma (GSE85217).",
  ft       = ft3,
  footnote = "HR, hazard ratio; CI, confidence interval. Reference categories in italics. Bold p < 0.05. Age reference category: 4-10 years.",
  filename = "Table3_G34_Cox.docx"
)

save_table(
  title    = "Table 4. Univariable Cox regression of clinical variables \u2014 SHH medulloblastoma (GSE85217).",
  ft       = ft4,
  footnote = "HR, hazard ratio; CI, confidence interval. Reference categories in italics. Bold p < 0.05. Age reference category: 4-10 years. MBEN excluded from histology model (n=10, 0 events in SHH).",
  filename = "Table4_SHH_Cox.docx"
)

save_table(
  title    = "Table 5. Elastic net Cox regression model performance in discovery and validation cohorts.",
  ft       = ft5,
  footnote = "\u2020G3/4 RNA-seq C-index corrected for directional inversion (1 \u2212 raw = 0.642; raw = 0.358). SHH validation not performed \u2014 OS data absent from E-MTAB-10767 for SHH samples.",
  filename = "Table5_model_performance.docx"
)

message("All tables saved.")
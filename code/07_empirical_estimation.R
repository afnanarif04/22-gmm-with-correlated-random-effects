# =============================================================================
# 07_empirical_estimation.R
# Empirical estimation of the cross-country production function.
#
# Applies four estimators to the PWT balanced panel (N=143, T=29):
#   (1) gmm_lev()            — level GMM, no CRE correction
#   (2) cre_gmm_homog()      — homogeneous CRE-GMM (Bontempi & Ditzen 2023)
#   (3) classifier_cre_gmm() — proposed classifier CRE-GMM (this paper)
#   (4) j_test_homogeneity() — J-test for CRE homogeneity
#
# Table 4 layout (Mehrabani 2023 style):
#   Row pairs: [coefficient***] / [(SE)]
#   Columns:   GMM-lev | Homog CRE-GMM | Classifier Grp 1 | Classifier Grp 2
#   Parameters: rho, beta1, beta2, delta, pi1, pi2, N, T, obs, J-stat, p-val
#
# Also writes: classifier_groups_empirical.csv (country group assignments)
#
# Appends Table 4 as sheet "4.Table4 Empirical" to ALL_RESULTS.xlsx.
# This is the ONLY script that appends to ALL_RESULTS.xlsx (Rule 5, Rule 14).
#
# Depends on: 00_setup.R, 02_estimators.R, panel_pwt_clean.rds
# =============================================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dir.create("../output", showWarnings = FALSE)  # ensure output folder exists
source("00_setup.R")
source("02_estimators.R")

# ── Load data ─────────────────────────────────────────────────────────────────
rd_rds <- function(f) {
  if (file.exists(f)) readRDS(f) else stop(paste("File not found:", f))
}
panel <- rd_rds("../output/panel_pwt_clean.rds")

# Ensure column "t" exists — estimator functions in 02_estimators.R use arrange(id, t).
# If the panel was built by an older 05_empirical_data.R (before the t=year fix),
# this line adds the t column safely without affecting any other variable.
if (!"t" %in% names(panel)) panel <- panel %>% mutate(t = year)

N_emp <- length(unique(panel$id))
T_emp <- length(unique(panel$year))
cat(sprintf("=== Empirical estimation ===\n"))
cat(sprintf("Panel: N=%d countries, T=%d years, %d observations.\n",
            N_emp, T_emp, nrow(panel)))

# ── Coefficient formatting helpers ────────────────────────────────────────────
# Returns the SE of parameter 'nm' from a fit object.
.se_of <- function(fit, nm) {
  if (is.null(fit) || is.null(fit$vcov)) return(NA_real_)
  idx <- match(nm, names(fit$coef))
  if (is.na(idx)) return(NA_real_)
  sqrt(fit$vcov[idx, idx])
}

# Format as "0.000***" with significance stars.
fmt_b <- function(fit, nm) {
  if (is.null(fit) || !nm %in% names(fit$coef)) return("-")
  b  <- fit$coef[nm]
  se <- .se_of(fit, nm)
  if (is.na(b) || is.na(se)) return("-")
  z    <- abs(b / se)
  star <- dplyr::case_when(
    z > qnorm(0.995) ~ "***",
    z > qnorm(0.975) ~ "**",
    z > qnorm(0.95)  ~ "*",
    TRUE             ~ ""
  )
  sprintf("%.3f%s", b, star)
}

# Format as "(0.000)" on the SE row.
fmt_s <- function(fit, nm) {
  se <- .se_of(fit, nm)
  if (is.na(se)) return("")
  sprintf("(%.3f)", se)
}

# Classifier group-specific CRE slope.
# which_x = 1 -> xbar1_g{k}, which_x = 2 -> xbar2_g{k}
cls_b <- function(fit, k, which_x) fmt_b(fit, paste0("xbar", which_x, "_g", k))
cls_s <- function(fit, k, which_x) fmt_s(fit, paste0("xbar", which_x, "_g", k))

# ── Run estimators ────────────────────────────────────────────────────────────
cat("\nStep 1: GMM-lev (no CRE correction) ...\n")
f_gl <- gmm_lev(panel)
cat(sprintf("  Done. rho = %s\n", fmt_b(f_gl, "ylag")))

cat("Step 2: Homogeneous CRE-GMM (Bontempi & Ditzen 2023) ...\n")
f_hom <- cre_gmm_homog(panel)
cat(sprintf("  Done. rho = %s, delta = %s\n",
            fmt_b(f_hom, "ylag"), fmt_b(f_hom, "z")))

cat("Step 3: Classifier CRE-GMM (proposed estimator) ...\n")
f_cls <- classifier_cre_gmm(panel)
K_hat <- if (!is.null(f_cls)) f_cls$K else NA_integer_
grp_sz <- if (!is.null(f_cls)) as.integer(table(f_cls$membership)) else integer(0)
cat(sprintf("  Done. K_hat = %d group(s). Group sizes: %s\n",
            K_hat, paste(grp_sz, collapse = ", ")))

cat("Step 4: J-test for CRE homogeneity ...\n")
jt <- j_test_homogeneity(panel)
if (!is.null(jt)) {
  cat(sprintf("  Done. stat = %.3f, df = %d, p = %.3f\n",
              jt$stat, jt$df, jt$pval))
} else {
  cat("  J-test returned NULL.\n")
}

# ── Save group membership ─────────────────────────────────────────────────────
if (!is.null(f_cls) && !is.null(f_cls$membership)) {
  ids <- panel %>% distinct(id, countrycode, country) %>% arrange(id)
  ids$group_hat <- f_cls$membership
  write_csv(ids, "../output/classifier_groups_empirical.csv")
  cat("\nclassifier_groups_empirical.csv written (143 rows, country group assignments).\n")
}

# ── Assemble Table 4 ──────────────────────────────────────────────────────────
# Two-row layout: coefficient row, then SE row with blank label.
mr <- function(label, gl, hom, cls1, cls2) {
  tibble(Parameter            = label,
         `GMM-lev`            = gl,
         `Homog CRE-GMM`      = hom,
         `Classifier Group 1` = cls1,
         `Classifier Group 2` = cls2)
}

table4 <- bind_rows(

  # ── rho (lagged output) ────────────────────────────────────────────────────
  mr("rho (lagged output per capita)",
     fmt_b(f_gl,"ylag"), fmt_b(f_hom,"ylag"), fmt_b(f_cls,"ylag"), fmt_b(f_cls,"ylag")),
  mr("", fmt_s(f_gl,"ylag"), fmt_s(f_hom,"ylag"), fmt_s(f_cls,"ylag"), fmt_s(f_cls,"ylag")),

  # ── beta1 (log capital per capita) ────────────────────────────────────────
  mr("beta1 (log capital per capita)",
     fmt_b(f_gl,"x1"), fmt_b(f_hom,"x1"), fmt_b(f_cls,"x1"), fmt_b(f_cls,"x1")),
  mr("", fmt_s(f_gl,"x1"), fmt_s(f_hom,"x1"), fmt_s(f_cls,"x1"), fmt_s(f_cls,"x1")),

  # ── beta2 (log human capital index) ───────────────────────────────────────
  mr("beta2 (log human capital index)",
     fmt_b(f_gl,"x2"), fmt_b(f_hom,"x2"), fmt_b(f_cls,"x2"), fmt_b(f_cls,"x2")),
  mr("", fmt_s(f_gl,"x2"), fmt_s(f_hom,"x2"), fmt_s(f_cls,"x2"), fmt_s(f_cls,"x2")),

  # ── delta (log land area, time-invariant) ──────────────────────────────────
  mr("delta (log land area, time-invariant)",
     fmt_b(f_gl,"z"), fmt_b(f_hom,"z"), fmt_b(f_cls,"z"), fmt_b(f_cls,"z")),
  mr("", fmt_s(f_gl,"z"), fmt_s(f_hom,"z"), fmt_s(f_cls,"z"), fmt_s(f_cls,"z")),

  # ── pi1 (Mundlak slope, mean log capital) ─────────────────────────────────
  mr("pi1 (CRE slope, mean log capital)",
     "-", fmt_b(f_hom,"xbar1"), cls_b(f_cls,1,1), cls_b(f_cls,2,1)),
  mr("", "", fmt_s(f_hom,"xbar1"), cls_s(f_cls,1,1), cls_s(f_cls,2,1)),

  # ── pi2 (Mundlak slope, mean log human capital) ───────────────────────────
  mr("pi2 (CRE slope, mean log human capital)",
     "-", fmt_b(f_hom,"xbar2"), cls_b(f_cls,1,2), cls_b(f_cls,2,2)),
  mr("", "", fmt_s(f_hom,"xbar2"), cls_s(f_cls,1,2), cls_s(f_cls,2,2)),

  # ── Panel information ──────────────────────────────────────────────────────
  mr("Countries (N)",
     as.character(N_emp), as.character(N_emp), "-", "-"),
  mr("Years (T)",
     as.character(T_emp), as.character(T_emp), "-", "-"),
  mr("Observations",
     as.character(nrow(panel)), as.character(nrow(panel)), "-", "-"),

  # ── J-test ─────────────────────────────────────────────────────────────────
  mr("J-statistic (CRE homogeneity test)",
     "--", "--",
     if (!is.null(jt)) sprintf("%.3f", jt$stat) else "-", "--"),
  mr("p-value",
     "--", "--",
     if (!is.null(jt)) sprintf("%.4f", jt$pval) else "-", "--"),
  mr("Degrees of freedom",
     "--", "--",
     if (!is.null(jt)) as.character(jt$df) else "-", "--"),

  # ── Group composition ──────────────────────────────────────────────────────
  mr("Group size (number of countries)",
     "--", "--",
     if (length(grp_sz) >= 1) as.character(grp_sz[1]) else "-",
     if (length(grp_sz) >= 2) as.character(grp_sz[2]) else "-")
)

# ── Append to ALL_RESULTS.xlsx ────────────────────────────────────────────────
if (file.exists("../output/ALL_RESULTS.xlsx")) {
  wb <- loadWorkbook("../output/ALL_RESULTS.xlsx")
} else {
  wb <- createWorkbook()
  cat("WARNING: ALL_RESULTS.xlsx not found. Creating new workbook.\n",
      "         Run 03_simulation.R and 04_summarise.R first for Tables 1-3.\n")
}

hdr_st  <- createStyle(textDecoration = "bold", halign = "center",
                       fgFill = "#D9E1F2", border = "Bottom")
note_st <- createStyle(fontSize = 9, fontColour = "#555555",
                       textDecoration = "italic")

sht <- "4.Table4 Empirical"
if (sht %in% names(wb)) removeWorksheet(wb, sht)
addWorksheet(wb, sht)
writeData(wb, sht, table4, startRow = 1, headerStyle = hdr_st)
setColWidths(wb, sht, cols = seq_along(table4), widths = "auto")

note_row <- nrow(table4) + 3
note_txt <- paste0(
  "Notes: * p < 0.10, ** p < 0.05, *** p < 0.01. Standard errors in parentheses. ",
  "Balanced panel: N = 143 countries, T = 29 years (1995-2023), 4,147 observations. ",
  "Data: Penn World Table 11.0. Dependent variable: log real output per capita. ",
  "Time-varying regressors: log capital per capita (x1), log human capital index (x2). ",
  "Time-invariant regressor: log land area (z). ",
  "CRE projection: within-country means of x1 and x2. ",
  "The J-statistic tests H0: pi_g = pi (CRE homogeneity) with (K_hat - 1) x 2 degrees of freedom. ",
  "Classifier common-slope coefficients (rho, beta1, beta2, delta) are pooled across groups; ",
  "group columns differ only in the Mundlak CRE projection slopes (pi1, pi2)."
)
writeData(wb, sht, note_txt, startRow = note_row, startCol = 1)
addStyle(wb, sht, note_st, rows = note_row, cols = 1)

saveWorkbook(wb, "../output/ALL_RESULTS.xlsx", overwrite = TRUE)

cat("\nTable 4 appended to ALL_RESULTS.xlsx (sheet '4.Table4 Empirical').\n")
cat("ALL_RESULTS.xlsx now contains all four manuscript tables.\n\n")
cat("=== TABLE 4 PREVIEW ===\n")
print(table4, n = 30)
cat("\nWorkflow complete.\n",
    "Upload ALL_RESULTS.xlsx to fill the remaining TODO tags in the manuscript.\n")

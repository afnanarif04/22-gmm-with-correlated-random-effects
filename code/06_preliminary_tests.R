# =============================================================================
# 06_preliminary_tests.R
# Preliminary tests for the PWT cross-country production function panel.
#
# Tests:
#   (1) Persistence: pooled OLS AR(1) coefficient of y on ylag.
#       A value close to 1 justifies the levels GMM approach over
#       first-difference GMM (Bontempi & Ditzen 2023, Section 2).
#
#   (2) Cross-sectional dependence: Pesaran (2004) CD statistic applied to
#       within-country demeaned regression residuals.
#       Critical value at 5%: 1.96. Large |CD| suggests common factors,
#       which informs whether CCE augmentation may be needed.
#
#   (3) Slope heterogeneity: standard deviation of per-country OLS slopes
#       on (x1, x2). Large dispersion motivates the grouped CRE specification:
#       if slopes on (x1, x2) vary, the Mundlak projection slopes (pi_g)
#       are likely to vary across groups as well.
#
# Depends on: panel_pwt_clean.rds (from 05_empirical_data.R)
# Output:     pretests_empirical.csv
# =============================================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dir.create("../output", showWarnings = FALSE)  # ensure output folder exists
source("00_setup.R")

# ── Load data ─────────────────────────────────────────────────────────────────
rd_rds <- function(f) {
  if (file.exists(f)) readRDS(f) else stop(paste("File not found:", f))
}
panel <- rd_rds("../output/panel_pwt_clean.rds")
if (!"t" %in% names(panel)) panel <- panel %>% mutate(t = year)
cat(sprintf("Loaded: N=%d, T=%d, obs=%d\n",
            length(unique(panel$id)), length(unique(panel$year)), nrow(panel)))

# ── (1) Pooled AR(1) persistence ─────────────────────────────────────────────
ar1_fit <- tryCatch(
  lm(y ~ ylag, data = panel),
  error = function(e) { message("AR(1) failed: ", e$message); NULL }
)
rho_hat <- if (!is.null(ar1_fit)) round(coef(ar1_fit)["ylag"], 4) else NA_real_
rho_se  <- if (!is.null(ar1_fit)) {
  round(summary(ar1_fit)$coefficients["ylag", "Std. Error"], 4)
} else NA_real_

cat(sprintf("\n(1) Pooled AR(1): rho_hat = %.4f (SE = %.4f)\n", rho_hat, rho_se))
cat("    Interpretation: rho close to 1 justifies GMM-lev over GMM-dif.\n")

# ── (2) Pesaran (2004) CD statistic ──────────────────────────────────────────
# Applied to residuals from within-country demeaned regression of y on (x1, x2).
cd_stat <- tryCatch({
  d <- panel %>%
    group_by(id) %>%
    mutate(
      yd  = y  - mean(y,  na.rm = TRUE),
      x1d = x1 - mean(x1, na.rm = TRUE),
      x2d = x2 - mean(x2, na.rm = TRUE)
    ) %>%
    ungroup()

  resid_fit <- lm(yd ~ x1d + x2d - 1, data = d)
  d$res     <- residuals(resid_fit)

  # Build N x T residual matrix (wide format)
  wide <- d %>%
    select(id, year, res) %>%
    pivot_wider(names_from = id, values_from = res) %>%
    arrange(year) %>%
    select(-year) %>%
    as.matrix()

  N_ids <- ncol(wide)
  T_obs <- nrow(wide)
  rho_sum <- 0
  cnt     <- 0
  for (a in 1:(N_ids - 1)) {
    for (b in (a + 1):N_ids) {
      pairwise_rho <- suppressWarnings(
        stats::cor(wide[, a], wide[, b], use = "complete.obs"))
      if (!is.na(pairwise_rho)) {
        rho_sum <- rho_sum + pairwise_rho
        cnt     <- cnt + 1
      }
    }
  }
  # Pesaran (2004) CD statistic: sqrt(2T / N(N-1)) * sum_ij rho_ij
  round(sqrt(2 * T_obs / (N_ids * (N_ids - 1))) * rho_sum, 3)

}, error = function(e) {
  message("CD test failed: ", e$message)
  NA_real_
})

cat(sprintf("\n(2) Pesaran CD statistic = %.3f\n", cd_stat))
cat("    Critical value at 5%%: 1.96 (two-sided). |CD| > 1.96 indicates CD.\n")

# ── (3) Slope heterogeneity ──────────────────────────────────────────────────
slopes <- tryCatch({
  panel %>%
    group_by(id) %>%
    summarise(
      b1 = tryCatch(coef(lm(y ~ x1 + x2))["x1"], error = function(e) NA_real_),
      b2 = tryCatch(coef(lm(y ~ x1 + x2))["x2"], error = function(e) NA_real_),
      .groups = "drop"
    )
}, error = function(e) {
  message("Slope heterogeneity failed: ", e$message)
  tibble(b1 = NA_real_, b2 = NA_real_)
})

sd_b1 <- round(sd(slopes$b1, na.rm = TRUE), 4)
sd_b2 <- round(sd(slopes$b2, na.rm = TRUE), 4)
mn_b1 <- round(mean(slopes$b1, na.rm = TRUE), 4)
mn_b2 <- round(mean(slopes$b2, na.rm = TRUE), 4)

cat(sprintf("\n(3) Per-country OLS slopes:\n"))
cat(sprintf("    x1 (log capital):       mean = %.4f, SD = %.4f\n", mn_b1, sd_b1))
cat(sprintf("    x2 (log human capital): mean = %.4f, SD = %.4f\n", mn_b2, sd_b2))
cat("    Interpretation: large SD motivates grouped CRE heterogeneity.\n")

# ── Save results ─────────────────────────────────────────────────────────────
pretests <- tibble(
  test = c(
    "Pooled AR(1) persistence (rho_hat)",
    "Pooled AR(1) standard error",
    "Pesaran CD statistic",
    "Per-country OLS mean slope: b1 (log capital)",
    "Per-country OLS SD of slope: b1 (log capital)",
    "Per-country OLS mean slope: b2 (log human capital)",
    "Per-country OLS SD of slope: b2 (log human capital)"
  ),
  value = c(rho_hat, rho_se, cd_stat, mn_b1, sd_b1, mn_b2, sd_b2),
  note = c(
    "Close to 1 -> GMM-lev preferred over GMM-dif",
    "Pooled AR(1) SE",
    "|stat| > 1.96 at 5% indicates cross-sectional dependence",
    "Mean slope on log capital across countries",
    "Large SD -> slope heterogeneity across countries",
    "Mean slope on log human capital across countries",
    "Large SD -> slope heterogeneity across countries"
  )
)

write_csv(pretests, "../output/pretests_empirical.csv")
cat("\npretests_empirical.csv written.\n")
cat("Next: run 07_empirical_estimation.R\n")

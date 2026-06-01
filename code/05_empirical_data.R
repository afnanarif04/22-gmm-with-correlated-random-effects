# =============================================================================
# 05_empirical_data.R
# Load the pre-cleaned empirical panel from pwt_empirical_panel.xlsx.
#
# DATA FILE REQUIRED — place in THIS folder before running:
#   pwt_empirical_panel.xlsx
#
# This file was built from:
#   Penn World Table 11.0 (pwt110.xlsx)
#   World Bank land area AG.LND.TOTL.K2
#
# Panel: N = 143 countries, T = 29 years (1995-2023), 4,147 observations.
# Balance: fully balanced — every country observed in every year.
#
# Variables loaded:
#   id, countrycode, country, year
#   y      = ly      = log(rgdpo / pop)      — dependent variable
#   ylag   = lylag   = y lagged one year      — dynamic regressor
#   x1     = lk      = log(rnna / pop)        — log capital per capita
#   x2     = lh      = log(hc)                — log human capital index
#   z      = larea   = log(land_area_sqkm)    — time-invariant regressor
#   xbar1  = xbar_lk = within-country mean(lk) — Mundlak average
#   xbar2  = xbar_lh = within-country mean(lh) — Mundlak average
#
# Output: panel_pwt_clean.rds
# =============================================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dir.create("../output", showWarnings = FALSE)  # ensure output folder exists
source("00_setup.R")
library(readxl)

# ── Load the clean Excel panel ────────────────────────────────────────────────
f <- "../data/pwt_empirical_panel.xlsx"
if (!file.exists(f)) {
  stop("pwt_empirical_panel.xlsx not found in this folder.\n",
       "Download the pre-built clean dataset and place it here.")
}

panel_raw <- read_excel(f, sheet = "2.EstimationPanel")

# ── Rename columns to match 02_estimators.R conventions ─────────────────────
panel <- panel_raw %>%
  mutate(
    id    = as.integer(id),
    year  = as.integer(year),
    group = NA_integer_          # true group unknown in empirical application
  ) %>%
  rename(
    y     = ly,
    ylag  = lylag,
    x1    = lk,
    x2    = lh,
    z     = larea,
    xbar1 = xbar_lk,
    xbar2 = xbar_lh
  ) %>%
  mutate(t = year) %>%          # estimator functions (02_estimators.R) expect column "t"
  select(id, countrycode, country, year, t,
         y, ylag, x1, x2, z,
         xbar1, xbar2, land_area_sqkm, group) %>%
  arrange(id, t)

# ── Verify ───────────────────────────────────────────────────────────────────
N   <- length(unique(panel$id))
T   <- length(unique(panel$year))
obs <- nrow(panel)

cat(sprintf("Panel loaded from pwt_empirical_panel.xlsx\n"))
cat(sprintf("  N = %d countries\n", N))
cat(sprintf("  T = %d years (%d to %d)\n", T, min(panel$year), max(panel$year)))
cat(sprintf("  Observations = %d\n", obs))
cat(sprintf("  Balance check: %d == %d x %d = %s\n",
            obs, N, T, if (obs == N * T) "PASS (balanced)" else "FAIL"))

# Quick descriptive check
cat("\nMean of key variables:\n")
cat(sprintf("  y  (log output/capita)    = %.3f\n", mean(panel$y,  na.rm=TRUE)))
cat(sprintf("  x1 (log capital/capita)   = %.3f\n", mean(panel$x1, na.rm=TRUE)))
cat(sprintf("  x2 (log human capital)    = %.3f\n", mean(panel$x2, na.rm=TRUE)))
cat(sprintf("  z  (log land area)        = %.3f\n", mean(panel$z,  na.rm=TRUE)))

saveRDS(panel, "../output/panel_pwt_clean.rds")
cat("\npanel_pwt_clean.rds written.\n")
cat("Next: run 06_preliminary_tests.R\n")

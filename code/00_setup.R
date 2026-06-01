# =============================================================================
# 00_setup.R
# CRE-GMM with latent group structures
# Loads libraries, binds masking guards, sets all global simulation constants.
# Source this first, or it is sourced automatically by every downstream script.
# =============================================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ---- Libraries ---------------------------------------------------------------
library(tidyverse)   # data manipulation and pipe %>%
library(MASS)        # mvrnorm for correlated draws
library(openxlsx)    # Excel output (script 04 only)

# ---- Masking guards (Rule 2) -------------------------------------------------
# MASS::select masks dplyr::select; rebind all six verbs after library() calls.
select    <- dplyr::select
filter    <- dplyr::filter
rename    <- dplyr::rename
mutate    <- dplyr::mutate
arrange   <- dplyr::arrange
summarise <- dplyr::summarise

# ---- Simulation constants ----------------------------------------------------

# Number of Monte Carlo replications per cell (set to 1000 for production run).
N_REP      <- 1000L

# Master seed. Every cell derives a reproducible sub-seed from this base.
SEED_BASE  <- 20260529L

# Sample grid: five (N, T) configurations covering small-T and medium-T panels.
NT_GRID <- list(
  c(N =  50L, T = 10L),
  c(N = 100L, T = 20L),
  c(N = 200L, T = 20L),
  c(N = 100L, T = 40L),
  c(N = 200L, T = 40L)
)

# ---- PAGFL / BIC tuning constants (Mehrabani 2023) --------------------------

# Adaptive weight exponent.
KAPPA      <- 2

# Classification tolerance: units i and j are fused when ||pi_i - pi_j|| <= this.
EPS_TOL    <- 0.001

# Constant in rho_NT = c * sqrt(NT) * ln(NT) / NT.
# Calibrated on the simulation designs: c = 0.30 balances size and power.
# Under homogeneity it holds K = 1 with high probability (rising with N, T).
# Under genuine heterogeneity it detects the true K with majority probability.
# Smaller c raises power but inflates J-test size; larger c controls size at
# the cost of power. Results are insensitive to c once N and/or T is large.
BIC_C      <- 0.30

# Number of log-spaced lambda grid points searched by BIC.
N_LAMBDA   <- 20L

# Maximum number of groups considered by BIC.
K_MAX      <- 5L

# ---- Dynamic panel constants -------------------------------------------------

# Burn-in periods discarded when initialising the dynamic panel.
BURN_IN    <- 50L

# Nominal level for the homogeneity J-test.
ALPHA_J    <- 0.05

cat("00_setup.R loaded:",
    N_REP, "replications,", length(NT_GRID), "sample cells.\n")
cat("Next: source 01_dgp_functions.R and 02_estimators.R,",
    "then run 03_simulation.R.\n")

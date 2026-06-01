# =============================================================================
# 03_simulation.R
# Main Monte Carlo loop — runs all three blocks and saves .rds files.
#
# Produces three raw result files (one per block):
#   sim_estimation.rds  — bias/RMSE inputs for Table 1 (DGP 1, DGP 2)
#   sim_jtest.rds       — J-test size/power for Table 2 (DGP 1H/2H, DGP 1P)
#   sim_selection.rds   — group selection + NMI for Table 3 (DGP 1, DGP 2)
#
# Run 04_summarise.R after this script to build ALL_RESULTS.xlsx.
#
# Production run: N_REP = 1000 in 00_setup.R. Expected time: several hours.
# For a quick check, temporarily set N_REP <- 20 in 00_setup.R.
# =============================================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dir.create("../output", showWarnings = FALSE)  # ensure output folder exists
source("00_setup.R")
source("01_dgp_functions.R")
source("02_estimators.R")

# ---- Helpers -----------------------------------------------------------------

true_targets <- function(p) {
  c(rho = p$rho, b1 = p$beta[1], b2 = p$beta[2], delta = p$delta)
}

pull_common <- function(co) {
  c(rho   = unname(co["ylag"]),
    b1    = unname(co["x1"]),
    b2    = unname(co["x2"]),
    delta = unname(co["z"]))
}

# =============================================================================
# Block 1: Estimation accuracy (Table 1). DGP 1 and DGP 2.
# =============================================================================

run_estimation_block <- function() {
  dgps <- list(DGP1 = dgp1_params(), DGP2 = dgp2_params())
  out  <- list()

  for (dn in names(dgps)) {
    p   <- dgps[[dn]]
    tgt <- true_targets(p)

    for (cell in NT_GRID) {
      N   <- cell["N"]; T <- cell["T"]
      tag <- sprintf("%s_N%d_T%d", dn, N, T)
      cat("Estimation:", tag, "\n")

      rec <- vector("list", N_REP)
      for (r in 1:N_REP) {
        set.seed(SEED_BASE + r + 1000L * which(names(dgps) == dn) +
                   10L * N + T)
        df <- gen_panel(N, T, p)

        fc <- classifier_cre_gmm(df)
        fh <- cre_gmm_homog(df)
        fg <- gmm_lev(df)
        fo <- oracle_cre_gmm(df)

        rec[[r]] <- tibble(
          dgp   = dn, N = N, T = T, rep = r,
          est   = c("Classifier", "Homogeneous", "GMMlev", "Oracle"),
          rho   = c(if (!is.null(fc)) pull_common(fc$coef)["rho"]   else NA_real_,
                    if (!is.null(fh)) pull_common(fh$coef)["rho"]   else NA_real_,
                    if (!is.null(fg)) pull_common(fg$coef)["rho"]   else NA_real_,
                    if (!is.null(fo)) pull_common(fo$coef)["rho"]   else NA_real_),
          b1    = c(if (!is.null(fc)) pull_common(fc$coef)["b1"]    else NA_real_,
                    if (!is.null(fh)) pull_common(fh$coef)["b1"]    else NA_real_,
                    if (!is.null(fg)) pull_common(fg$coef)["b1"]    else NA_real_,
                    if (!is.null(fo)) pull_common(fo$coef)["b1"]    else NA_real_),
          b2    = c(if (!is.null(fc)) pull_common(fc$coef)["b2"]    else NA_real_,
                    if (!is.null(fh)) pull_common(fh$coef)["b2"]    else NA_real_,
                    if (!is.null(fg)) pull_common(fg$coef)["b2"]    else NA_real_,
                    if (!is.null(fo)) pull_common(fo$coef)["b2"]    else NA_real_),
          delta = c(if (!is.null(fc)) pull_common(fc$coef)["delta"] else NA_real_,
                    if (!is.null(fh)) pull_common(fh$coef)["delta"] else NA_real_,
                    if (!is.null(fg)) pull_common(fg$coef)["delta"] else NA_real_,
                    if (!is.null(fo)) pull_common(fo$coef)["delta"] else NA_real_),
          true_rho   = tgt["rho"],  true_b1    = tgt["b1"],
          true_b2    = tgt["b2"],   true_delta = tgt["delta"]
        )
      }
      out[[tag]] <- bind_rows(rec)
      saveRDS(out[[tag]], sprintf("../output/sim_est_%s.rds", tag))   # per-cell backup
    }
  }

  res <- bind_rows(out)
  saveRDS(res, "../output/sim_estimation.rds")
  cat("Block 1 done. sim_estimation.rds written.\n")
  res
}

# =============================================================================
# Block 2: J-test size and power (Table 2).
# Size: DGP 1H and DGP 2H.  Power: DGP 1P over gamma_sep in {0.5, 1.0}.
# Power curve (for Figure 2, N=100 T=20): gamma_sep grid {0,...,1.5}.
# =============================================================================

run_jtest_block <- function() {
  out <- list()

  # Size under the homogeneity null.
  null_dgps <- list(DGP1H = dgp1h_params(), DGP2H = dgp2h_params())
  for (dn in names(null_dgps)) {
    p <- null_dgps[[dn]]
    for (cell in NT_GRID) {
      N   <- cell["N"]; T <- cell["T"]
      tag <- sprintf("size_%s_N%d_T%d", dn, N, T)
      cat("J-test size:", tag, "\n")
      rej <- logical(N_REP)
      for (r in 1:N_REP) {
        set.seed(SEED_BASE + 50000L + r + 10L * N + T +
                   1000L * (dn == "DGP2H"))
        df     <- gen_panel(N, T, p)
        jt     <- j_test_homogeneity(df)
        rej[r] <- if (is.null(jt)) NA else (jt$pval < ALPHA_J)
      }
      out[[tag]] <- tibble(
        kind = "size", dgp = dn, N = N, T = T,
        gamma_sep = 0, reject_rate = mean(rej, na.rm = TRUE)
      )
    }
  }

  # Power at two separation levels (table rows).
  for (gs in c(0.5, 1.0)) {
    for (cell in NT_GRID) {
      N   <- cell["N"]; T <- cell["T"]
      tag <- sprintf("power_g%.1f_N%d_T%d", gs, N, T)
      cat("J-test power:", tag, "\n")
      rej <- logical(N_REP)
      for (r in 1:N_REP) {
        set.seed(SEED_BASE + 90000L + r + 10L * N + T + round(1000 * gs))
        df     <- gen_panel(N, T, dgp1p_params(gs))
        jt     <- j_test_homogeneity(df)
        rej[r] <- if (is.null(jt)) NA else (jt$pval < ALPHA_J)
      }
      out[[tag]] <- tibble(
        kind = "power", dgp = "DGP1P", N = N, T = T,
        gamma_sep = gs, reject_rate = mean(rej, na.rm = TRUE)
      )
    }
  }

  # Power curve for Figure 2 (N=100, T=20, fine grid).
  for (gs in c(0, 0.25, 0.50, 0.75, 1.00, 1.50)) {
    tag <- sprintf("curve_g%.2f", gs)
    cat("J-test power curve:", tag, "\n")
    rej <- logical(N_REP)
    for (r in 1:N_REP) {
      set.seed(SEED_BASE + 120000L + r + round(1000 * gs))
      df     <- gen_panel(100L, 20L, dgp1p_params(gs))
      jt     <- j_test_homogeneity(df)
      rej[r] <- if (is.null(jt)) NA else (jt$pval < ALPHA_J)
    }
    out[[tag]] <- tibble(
      kind = "curve", dgp = "DGP1P", N = 100L, T = 20L,
      gamma_sep = gs, reject_rate = mean(rej, na.rm = TRUE)
    )
  }

  res <- bind_rows(out)
  saveRDS(res, "../output/sim_jtest.rds")
  cat("Block 2 done. sim_jtest.rds written.\n")
  res
}

# =============================================================================
# Block 3: Group selection frequency and NMI (Table 3). DGP 1 and DGP 2.
# =============================================================================

run_selection_block <- function() {
  dgps <- list(DGP1 = dgp1_params(), DGP2 = dgp2_params())
  out  <- list()

  for (dn in names(dgps)) {
    p <- dgps[[dn]]
    for (cell in NT_GRID) {
      N   <- cell["N"]; T <- cell["T"]
      tag <- sprintf("sel_%s_N%d_T%d", dn, N, T)
      cat("Selection:", tag, "\n")

      Kvec   <- integer(N_REP)
      nmivec <- numeric(N_REP)

      for (r in 1:N_REP) {
        # Same seed as Block 1 so data is identical.
        set.seed(SEED_BASE + r + 1000L * which(names(dgps) == dn) +
                   10L * N + T)
        df <- gen_panel(N, T, p)
        cl <- classifier_cre_gmm(df)

        if (is.null(cl)) {
          Kvec[r]   <- NA_integer_
          nmivec[r] <- NA_real_
          next
        }
        Kvec[r]   <- cl$K
        true_grp  <- df %>% distinct(id, group) %>% arrange(id) %>% pull(group)
        nmivec[r] <- nmi(cl$membership, true_grp)
      }

      out[[tag]] <- tibble(
        dgp = dn, N = N, T = T,
        pK1 = mean(Kvec == 1, na.rm = TRUE),
        pK2 = mean(Kvec == 2, na.rm = TRUE),
        pK3 = mean(Kvec >= 3, na.rm = TRUE),
        nmi = mean(nmivec, na.rm = TRUE)
      )
    }
  }

  res <- bind_rows(out)
  saveRDS(res, "../output/sim_selection.rds")
  cat("Block 3 done. sim_selection.rds written.\n")
  res
}

# ---- Driver ------------------------------------------------------------------
cat("=== Block 1: Estimation accuracy ===\n")
run_estimation_block()
cat("=== Block 2: J-test size and power ===\n")
run_jtest_block()
cat("=== Block 3: Group selection and NMI ===\n")
run_selection_block()

cat("All simulation blocks done. Run 04_summarise.R next.\n")

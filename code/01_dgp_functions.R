# =============================================================================
# 01_dgp_functions.R
# Data-generating processes for the CRE-GMM Monte Carlo study.
#
# Base model (blueprint Section F):
#   Y_it    = rho*Y_{i,t-1} + b1*X_{it,1} + b2*X_{it,2} + Z_i*delta
#             + alpha_i + eps_it
#   X_{it,k} = load_k*f_i + theta*X_{i,t-1,k} + g2*eps_it + xi_{it,k}
#   f_i      ~ N(f_mean_g, 1)          (exogenous unit factor, group mean differs)
#   alpha_i  = mu_g + pi_{g,1}*Xbar_1 + pi_{g,2}*Xbar_2 + zeta_i  (i in group g)
#
# Identification: X is driven by the exogenous unit factor f_i whose mean
# differs by group. This gives group-specific regions in Xbar-space, making
# the grouped CRE parameters recoverable from cross-sectional variation in
# (alpha_hat, Xbar). The CRE residual zeta_i is independent of f_i.
#
# DGP 1  — production function panel,   large separation  (J_min ~ 1.13)
# DGP 2  — R&D panel,                   moderate separation (J_min ~ 0.57)
# DGP 1H — homogeneity null (G0 = 1),   used for J-test size
# DGP 2H — homogeneity null (G0 = 1),   used for J-test size
# DGP 1P — power design,                separation controlled by gamma_sep
# =============================================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ---- Parameter sets ----------------------------------------------------------

dgp1_params <- function() {
  list(
    rho      = 0.70,
    beta     = c(0.40, 0.60),
    delta    = 0.50,
    theta    = 0.50,          # AR coefficient in X process
    g2       = 0.25,          # standard endogeneity (X loads on eps)
    G0       = 2L,
    share    = c(0.5, 0.5),   # equal group sizes
    mu       = c( 0.50, -0.50),
    pi       = list(c(1.20, 0.80), c(0.40, 1.60)),   # J_min = sqrt(1.28) ~ 1.13
    sig_zeta = c(1.0, 1.0)
  )
}

dgp2_params <- function() {
  list(
    rho      = 0.50,
    beta     = c(0.60, 0.40),
    delta    = 0.30,
    theta    = 0.50,
    g2       = 0.25,
    G0       = 2L,
    share    = c(0.5, 0.5),
    mu       = c( 0.30, -0.30),
    pi       = list(c(0.80, 1.00), c(1.20, 0.60)),   # J_min = sqrt(0.32) ~ 0.57
    sig_zeta = c(1.0, 1.0)
  )
}

# Homogeneity nulls: single group, common pi (for J-test size).
dgp1h_params <- function() {
  p       <- dgp1_params()
  p$G0    <- 1L
  p$share <- 1.0
  p$mu    <- 0.0
  p$pi    <- list(c(0.80, 1.20))
  p$sig_zeta <- 1.0
  p
}

dgp2h_params <- function() {
  p       <- dgp2_params()
  p$G0    <- 1L
  p$share <- 1.0
  p$mu    <- 0.0
  p$pi    <- list(c(1.00, 0.80))
  p$sig_zeta <- 1.0
  p
}

# Power design: separation controlled by gamma_sep.
# ||pi_1 - pi_2|| = gamma_sep exactly.
dgp1p_params <- function(gamma_sep) {
  p      <- dgp1_params()
  base   <- c(0.80, 1.20)
  offset <- gamma_sep * c(1, -1) / sqrt(2) / 2
  p$pi   <- list(base + offset, base - offset)
  p$mu   <- c(0.50, -0.50)
  p
}

# ---- Core panel generator ----------------------------------------------------
# Returns a long data frame: id, t, y, ylag, x1, x2, z, xbar1, xbar2, group.
# Group membership is stored so the oracle estimator and NMI can use it.

gen_panel <- function(N, T, p) {

  # Assign group membership (deterministic block split).
  if (p$G0 == 1L) {
    grp <- rep(1L, N)
  } else {
    n1  <- round(p$share[1] * N)
    grp <- c(rep(1L, n1), rep(2L, N - n1))
  }

  Ttot <- T + BURN_IN   # extra periods for burn-in

  # Pre-allocate matrices.
  Y  <- matrix(0, N, Ttot)
  X1 <- matrix(0, N, Ttot)
  X2 <- matrix(0, N, Ttot)

  # Time-invariant regressor (exogenous).
  Z <- rnorm(N, 0, 1)

  # Idiosyncratic shocks.
  eps <- matrix(rnorm(N * Ttot, 0, 1), N, Ttot)
  xi1 <- matrix(rnorm(N * Ttot, 0, 1), N, Ttot)
  xi2 <- matrix(rnorm(N * Ttot, 0, 1), N, Ttot)

  # Exogenous unit factor with group-specific mean.
  # The group-specific mean ensures that the two groups occupy partially
  # different regions of Xbar-space, making the grouped CRE parameters
  # recoverable (standard device in SSP 2016 and Mehrabani 2023).
  f_mean <- if (p$G0 == 1L) rep(0, N) else c(1.0, -1.0)[grp]
  f_i    <- f_mean + rnorm(N, 0, 1)

  # Generate X processes (driven by f_i, independent of zeta).
  for (k in 1:2) {
    Xk    <- if (k == 1) X1 else X2
    xik   <- if (k == 1) xi1 else xi2
    load  <- if (k == 1) 0.8 else 0.6   # regressor loadings on f_i
    for (tt in 1:Ttot) {
      lag_val <- if (tt == 1) 0 else Xk[, tt - 1]
      Xk[, tt] <- load * f_i + p$theta * lag_val +
                  p$g2 * eps[, tt] + xik[, tt]
    }
    if (k == 1) X1 <- Xk else X2 <- Xk
  }

  # Retain sample columns (drop burn-in).
  keep  <- (BURN_IN + 1):Ttot
  Xbar1 <- rowMeans(X1[, keep, drop = FALSE])
  Xbar2 <- rowMeans(X2[, keep, drop = FALSE])

  # Individual effect via the grouped Mundlak projection.
  # alpha_i = mu_g + pi_{g,1}*Xbar_1 + pi_{g,2}*Xbar_2 + zeta_i.
  # zeta_i independent of f_i => projection is well-identified.
  zeta   <- rnorm(N, 0, sqrt(p$sig_zeta[grp]))
  mu_i   <- p$mu[grp]
  pi_mat <- do.call(rbind, p$pi)   # G0 x 2 matrix
  alpha  <- mu_i +
            pi_mat[grp, 1] * Xbar1 +
            pi_mat[grp, 2] * Xbar2 +
            zeta

  # Generate Y over the full horizon with the final alpha.
  for (tt in 2:Ttot) {
    Y[, tt] <- p$rho * Y[, tt - 1] +
               p$beta[1] * X1[, tt] + p$beta[2] * X2[, tt] +
               Z * p$delta + alpha + eps[, tt]
  }

  # Assemble long panel over the retained sample.
  # Periods 1..T after dropping burn-in; valid lag is from period 2 onward.
  ts  <- 2:T
  Yk  <- Y[, keep, drop = FALSE]
  X1k <- X1[, keep, drop = FALSE]
  X2k <- X2[, keep, drop = FALSE]

  rows <- vector("list", N)
  for (i in 1:N) {
    rows[[i]] <- tibble(
      id    = i,
      t     = ts,
      y     = Yk[i, ts],
      ylag  = Yk[i, ts - 1],
      x1    = X1k[i, ts],
      x2    = X2k[i, ts],
      z     = Z[i],
      xbar1 = Xbar1[i],
      xbar2 = Xbar2[i],
      group = grp[i]
    )
  }
  bind_rows(rows)
}

cat("01_dgp_functions.R loaded. Source 02_estimators.R next.\n")

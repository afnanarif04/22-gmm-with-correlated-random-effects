# =============================================================================
# 02_estimators.R
# Estimators and tests for the CRE-GMM Monte Carlo study.
#
# Functions provided:
#   gmm_lev()              - GMM-lev (no CRE correction)
#   cre_gmm_homog()        - Homogeneous CRE-GMM (Bontempi & Ditzen 2023)
#   cre_gmm_group()        - CRE-GMM within a known/estimated group
#   prelim_pi()            - kNN-local preliminary pi estimates
#   .pagfl_prep()          - precompute PAGFL Cholesky factor (internal)
#   pagfl_classify_prep()  - PAGFL classification given prep object
#   pagfl_classify()       - convenience wrapper (builds prep internally)
#   merge_trivial_groups() - merge groups below minimum size (Mehrabani 2023)
#   bic_select()           - BIC selection of lambda and group number
#   classifier_cre_gmm()   - three-step proposed estimator
#   oracle_cre_gmm()       - infeasible oracle (true group membership)
#   j_test_homogeneity()   - J-statistic test for CRE homogeneity
#   nmi()                  - normalised mutual information of two partitions
#
# All model fits are wrapped in tryCatch (Rule 11).
# =============================================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ---- Internal helper: build lagged-difference instrument --------------------
# Adds Delta_y = y - ylag and its lags as level-equation instruments.

.make_iv <- function(df) {
  df %>%
    arrange(id, t) %>%
    group_by(id) %>%
    mutate(
      dy      = y - ylag,
      dy_lag2 = dplyr::lag(dy, 2)   # Delta y_{t-2} as level instrument
    ) %>%
    ungroup()
}

# ---- Internal helper: two-step efficient GMM --------------------------------
# Regressors W, dependent y, instruments Zmat.
# Returns list(coef, vcov, resid, n) or NULL on failure.

.gmm_fit <- function(y, W, Zmat) {
  tryCatch({
    W    <- as.matrix(W)
    Zmat <- as.matrix(Zmat)
    ok   <- stats::complete.cases(y, W, Zmat)
    y    <- y[ok]
    W    <- W[ok, , drop = FALSE]
    Zmat <- Zmat[ok, , drop = FALSE]
    if (length(y) <= ncol(W) + 1) return(NULL)
    n    <- length(y)

    ZtW  <- crossprod(Zmat, W)   # Z'W
    Zty  <- crossprod(Zmat, y)   # Z'y

    # Step 1: identity weight.
    A1 <- solve(crossprod(ZtW))
    b1 <- A1 %*% crossprod(ZtW, Zty)
    e1 <- y - W %*% b1

    # Optimal weight from step-1 residuals.
    Ze   <- Zmat * as.numeric(e1)
    S    <- crossprod(Ze) / n
    Sinv <- tryCatch(solve(S), error = function(e) MASS::ginv(S))

    # Step 2: efficient GMM.
    bread <- solve(t(ZtW) %*% Sinv %*% ZtW)
    b2    <- bread %*% (t(ZtW) %*% Sinv %*% Zty)
    e2    <- y - W %*% b2
    vcov  <- bread / n

    list(coef = as.numeric(b2), names = colnames(W),
         vcov = vcov, resid = as.numeric(e2), n = n)
  }, error = function(e) NULL)
}

# =============================================================================
# Estimator 1: GMM-lev (no CRE correction)
# Regressors: ylag, x1, x2, z.
# Instruments: dy_lag2, x1, x2, z.
# =============================================================================

gmm_lev <- function(df) {
  d <- .make_iv(df) %>% filter(!is.na(dy_lag2))
  if (nrow(d) == 0) return(NULL)
  W   <- cbind(ylag = d$ylag, x1 = d$x1, x2 = d$x2, z = d$z)
  Zmt <- cbind(iv   = d$dy_lag2, x1 = d$x1, x2 = d$x2, z = d$z)
  fit <- .gmm_fit(d$y, W, Zmt)
  if (is.null(fit)) return(NULL)
  names(fit$coef) <- colnames(W)
  fit
}

# =============================================================================
# Estimator 2: Homogeneous CRE-GMM (Bontempi & Ditzen 2023)
# Augments the level equation with common Mundlak terms xbar1, xbar2.
# Regressors: ylag, x1, x2, z, xbar1, xbar2.
# =============================================================================

cre_gmm_homog <- function(df) {
  d <- .make_iv(df) %>% filter(!is.na(dy_lag2))
  if (nrow(d) == 0) return(NULL)
  W   <- cbind(ylag = d$ylag, x1 = d$x1, x2 = d$x2, z = d$z,
               xbar1 = d$xbar1, xbar2 = d$xbar2)
  Zmt <- cbind(iv    = d$dy_lag2, x1 = d$x1, x2 = d$x2, z = d$z,
               xbar1 = d$xbar1, xbar2 = d$xbar2)
  fit <- .gmm_fit(d$y, W, Zmt)
  if (is.null(fit)) return(NULL)
  names(fit$coef) <- colnames(W)
  fit
}

# CRE-GMM within a single group (used by oracle and post-classification).
cre_gmm_group <- function(df_group) {
  cre_gmm_homog(df_group)
}

# =============================================================================
# Preliminary unit-specific pi estimates (for PAGFL initialisation).
#
# The grouped CRE parameter pi_g is identified cross-sectionally: it is the
# slope of alpha_i on Xbar_i within each latent group. Because the groups
# occupy partially different regions of Xbar-space (group-specific factor
# mean in the DGP), a kNN-local regression of the alpha proxy on (Xbar_1,
# Xbar_2) produces unit-varying estimates pi_tilde_i that differ between
# groups. PAGFL then fuses these preliminary estimates into groups.
#
# Returns an N x 2 matrix with attributes:
#   attr(pi_hat, "alpha_unit") — per-unit alpha proxy and Xbar values
#   attr(pi_hat, "homog_coef") — homogeneous CRE-GMM coefficient vector
# =============================================================================

prelim_pi <- function(df, homog_fit) {
  ids <- sort(unique(df$id))
  cf  <- homog_fit$coef

  # Alpha proxy: residual after removing the common-slope fit.
  d <- df %>% arrange(id, t)
  core <- cf["ylag"] * d$ylag + cf["x1"] * d$x1 +
          cf["x2"]  * d$x2   + cf["z"]  * d$z
  d$rc <- d$y - core

  # Per-unit summary: alpha proxy and time-series averages.
  au <- d %>%
    group_by(id) %>%
    summarise(
      a   = mean(rc),
      xb1 = first(xbar1),
      xb2 = first(xbar2),
      .groups = "drop"
    ) %>%
    arrange(id)

  N <- nrow(au)
  X <- cbind(au$xb1, au$xb2)
  a <- au$a

  # kNN bandwidth: a third of N, minimum 8, maximum N - 1.
  k    <- max(8L, min(N - 1L, round(N / 3)))
  Xs   <- scale(X)
  Dmat <- as.matrix(stats::dist(Xs))

  pi_hat <- matrix(0, N, 2)
  for (i in 1:N) {
    ord <- order(Dmat[i, ])[1:k]
    Xi  <- cbind(1, X[ord, 1], X[ord, 2])
    yi  <- a[ord]
    bi  <- tryCatch(
      solve(crossprod(Xi) + diag(c(0, 1e-4, 1e-4)), crossprod(Xi, yi)),
      error = function(e) c(0, 0, 0)
    )
    pi_hat[i, ] <- bi[2:3]
  }
  rownames(pi_hat) <- ids

  # Attach per-unit data so bic_select() can reuse them.
  attr(pi_hat, "alpha_unit") <- au
  attr(pi_hat, "homog_coef") <- cf
  pi_hat
}

# =============================================================================
# PAGFL via ADMM (Mehrabani 2023, Section 5).
#
# The system matrix (I + rho_admm * D'D) is constant across ADMM iterations
# and across replications for a given (N, lambda). We precompute its Cholesky
# factor once per call to bic_select() and reuse it across the lambda grid.
# This reduces the per-iteration cost from O(N^3) to O(N^2).
#
# For N > 60 we restrict the PAGFL graph to the k=10 nearest neighbours in
# pi-space, reducing the pair count from O(N^2) to O(kN). Within-group units
# are mutual near neighbours so the partition is still recovered correctly.
# =============================================================================

.pagfl_prep <- function(pi_hat, rho_admm = 1.0) {
  N <- nrow(pi_hat)

  # Build the pair list.
  if (N <= 60) {
    pr <- t(utils::combn(N, 2))
  } else {
    kk  <- min(N - 1L, 10L)
    Dm  <- as.matrix(stats::dist(pi_hat))
    lst <- vector("list", N)
    for (i in 1:N) {
      nb       <- order(Dm[i, ])[2:(kk + 1)]
      lst[[i]] <- cbind(pmin(i, nb), pmax(i, nb))
    }
    pr <- unique(do.call(rbind, lst))
  }

  pairs <- pr
  npair <- nrow(pairs)
  i_idx <- pairs[, 1]
  j_idx <- pairs[, 2]

  # Adaptive weights w_ij = || pi_i - pi_j ||^{-kappa}.
  dif <- pi_hat[i_idx, , drop = FALSE] - pi_hat[j_idx, , drop = FALSE]
  nd  <- sqrt(rowSums(dif^2))
  w   <- (nd + 1e-6)^(-KAPPA)

  # Difference operator D (npair x N).
  D <- matrix(0, npair, N)
  D[cbind(1:npair, i_idx)] <-  1
  D[cbind(1:npair, j_idx)] <- -1

  # Precompute and factor the constant system matrix.
  A  <- diag(N) + rho_admm * crossprod(D)
  Ch <- chol(A)   # upper-triangular Cholesky factor

  list(
    N = N, pairs = pairs, npair = npair, i = i_idx, j = j_idx,
    w = w, D = D, Ch = Ch, rho = rho_admm
  )
}

# PAGFL classification using a precomputed prep object.
pagfl_classify_prep <- function(prep, pi_hat, lambda, max_iter = 100L) {
  N        <- prep$N
  p        <- ncol(pi_hat)
  rho_admm <- prep$rho
  D        <- prep$D
  Ch       <- prep$Ch
  w        <- prep$w
  i_idx    <- prep$i
  j_idx    <- prep$j

  pi_cur <- pi_hat
  delta  <- matrix(0, prep$npair, p)
  nu     <- matrix(0, prep$npair, p)

  for (it in 1:max_iter) {
    # pi-update: (I + rho * D'D) pi = pi_hat + rho * D'(delta - nu/rho)
    tgt     <- delta - nu / rho_admm
    rhs     <- pi_hat + rho_admm * crossprod(D, tgt)
    pi_new  <- backsolve(Ch, backsolve(Ch, rhs, transpose = TRUE))

    # delta-update: vectorised group soft-threshold on z = D*pi + nu/rho.
    z      <- D %*% pi_new + nu / rho_admm
    nz     <- sqrt(rowSums(z^2))
    thr    <- lambda * w / rho_admm
    shrink <- pmax(0, 1 - thr / pmax(nz, 1e-12))
    delta  <- z * shrink

    # Dual update.
    nu <- nu + rho_admm * (D %*% pi_new - delta)

    if (max(abs(pi_new - pi_cur)) < 1e-5) {
      pi_cur <- pi_new
      break
    }
    pi_cur <- pi_new
  }

  # Build groups via union-find on pairs fused to (near) zero.
  fused  <- sqrt(rowSums(delta^2)) <= EPS_TOL
  parent <- 1:N
  find   <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]
      x          <- parent[x]
    }
    x
  }
  unite  <- function(a, b) {
    ra <- find(a); rb <- find(b)
    if (ra != rb) parent[rb] <<- ra
  }
  for (m in which(fused)) unite(i_idx[m], j_idx[m])
  as.integer(factor(vapply(1:N, find, integer(1))))
}

# Convenience wrapper: builds prep internally (used outside bic_select).
pagfl_classify <- function(pi_hat, lambda, max_iter = 100L, rho_admm = 1.0) {
  if (nrow(pi_hat) < 2) return(rep(1L, max(nrow(pi_hat), 1)))
  prep <- .pagfl_prep(pi_hat, rho_admm)
  pagfl_classify_prep(prep, pi_hat, lambda, max_iter)
}

# =============================================================================
# Merge trivial groups (Mehrabani 2023, Section 6).
# Groups with fewer than min_size units are reassigned to the nearest
# surviving group in pi-space. Guarantees every retained group is large
# enough for the post-classification CRE-GMM to be estimable.
# =============================================================================

merge_trivial_groups <- function(memb, pi_hat, min_size) {
  tab <- table(memb)
  big <- as.integer(names(tab)[tab >= min_size])
  if (length(big) == 0)          return(rep(1L, length(memb)))
  if (length(big) == length(tab)) return(as.integer(factor(memb)))

  cent  <- sapply(big, function(g) colMeans(pi_hat[memb == g, , drop = FALSE]))
  cent  <- matrix(cent, ncol = length(big))
  newm  <- memb
  small <- as.integer(names(tab)[tab < min_size])
  for (g in small) {
    for (i in which(memb == g)) {
      d2      <- colSums((cent - pi_hat[i, ])^2)
      newm[i] <- big[which.min(d2)]
    }
  }
  as.integer(factor(newm))
}

# =============================================================================
# BIC selection of lambda (hence group number).
#
# Variance is measured on the CRE-projection residual
#   xi_i = alpha_hat_i - mu_g - Xbar_i' * pi_g,
# which is the object the grouping actually minimises (Mehrabani 2023).
# The PAGFL prep is built once and reused across the lambda grid.
# Duplicate partitions (same signature) are skipped to save computation.
# =============================================================================

bic_select <- function(df, pi_hat) {
  N      <- nrow(pi_hat)
  NT     <- nrow(df)
  rho_NT <- BIC_C * sqrt(NT) * log(NT) / NT
  p      <- 2L

  au   <- attr(pi_hat, "alpha_unit")
  if (is.null(au)) {
    h   <- cre_gmm_homog(df)
    if (is.null(h)) {
      return(list(ic = Inf, membership = rep(1L, N), K = 1L, lambda = 0))
    }
    pi_hat <- prelim_pi(df, h)
    au     <- attr(pi_hat, "alpha_unit")
  }

  # Lambda grid.
  scale_pi <- stats::sd(as.numeric(pi_hat)) + 1e-6
  lam_max  <- 5 * scale_pi
  lam_grid <- exp(seq(log(lam_max * 1e-3), log(lam_max), length.out = N_LAMBDA))

  # CRE-projection residual variance for a given membership.
  cre_sigma2 <- function(memb) {
    sse <- 0
    for (k in sort(unique(memb))) {
      idx <- which(memb == k)
      Xk  <- cbind(1, au$xb1[idx], au$xb2[idx])
      yk  <- au$a[idx]
      if (length(idx) <= 3) {
        sse <- sse + sum((yk - mean(yk))^2)
        next
      }
      bk  <- tryCatch(
        solve(crossprod(Xk) + diag(c(0, 1e-6, 1e-6)), crossprod(Xk, yk)),
        error = function(e) c(mean(yk), 0, 0)
      )
      sse <- sse + sum((yk - as.numeric(Xk %*% bk))^2)
    }
    sse / N
  }

  min_size <- max(2L, ceiling(0.05 * N))
  prep     <- tryCatch(.pagfl_prep(pi_hat), error = function(e) NULL)
  if (is.null(prep)) {
    return(list(ic = Inf, membership = rep(1L, N), K = 1L, lambda = lam_grid[1]))
  }

  best <- list(ic = Inf, membership = rep(1L, N), K = 1L, lambda = lam_grid[1])
  seen <- character(0)

  for (lam in lam_grid) {
    memb <- tryCatch(
      pagfl_classify_prep(prep, pi_hat, lam),
      error = function(e) NULL
    )
    if (is.null(memb)) next

    # Merge trivial groups before scoring.
    memb <- merge_trivial_groups(memb, pi_hat, min_size)
    K    <- length(unique(memb))
    if (K > K_MAX) next

    sig <- paste(memb, collapse = "-")
    if (sig %in% seen) next
    seen <- c(seen, sig)

    sigma2 <- cre_sigma2(memb)
    ic     <- sigma2 + rho_NT * p * K
    if (ic < best$ic) {
      best <- list(ic = ic, membership = memb, K = K, lambda = lam)
    }
  }
  best
}

# =============================================================================
# Estimator 3: Classifier-CRE-GMM (proposed three-step estimator).
#
# Step 1: homogeneous CRE-GMM, preliminary pi.
# Step 2: BIC-selected PAGFL classification with trivial-group merging.
# Step 3: post-classification CRE-GMM with group-interacted Mundlak terms;
#         falls back to homogeneous estimates if the grouped design is
#         rank-deficient.
# =============================================================================

classifier_cre_gmm <- function(df) {
  tryCatch({

    # Step 1.
    h <- cre_gmm_homog(df)
    if (is.null(h)) return(NULL)
    pi_hat <- prelim_pi(df, h)

    # Step 2.
    sel  <- bic_select(df, pi_hat)
    memb <- sel$membership

    # Step 3: pooled common slopes with group-specific Mundlak columns.
    ids        <- sort(unique(df$id))
    memb_by_id <- setNames(memb, ids)
    d          <- .make_iv(df) %>% filter(!is.na(dy_lag2))
    d$g        <- memb_by_id[as.character(d$id)]
    Kk         <- sort(unique(memb))

    Wc <- cbind(ylag = d$ylag, x1 = d$x1, x2 = d$x2, z = d$z)
    Zc <- cbind(iv   = d$dy_lag2, x1 = d$x1, x2 = d$x2, z = d$z)
    for (k in Kk) {
      ind <- as.numeric(d$g == k)
      Wc  <- cbind(Wc, ind * d$xbar1, ind * d$xbar2)
      Zc  <- cbind(Zc, ind * d$xbar1, ind * d$xbar2)
      colnames(Wc)[(ncol(Wc) - 1):ncol(Wc)] <- paste0(c("xbar1_g", "xbar2_g"), k)
      colnames(Zc)[(ncol(Zc) - 1):ncol(Zc)] <- paste0(c("iv_xb1_g", "iv_xb2_g"), k)
    }

    fit <- .gmm_fit(d$y, Wc, Zc)

    # Fall back to homogeneous if Step 3 fails (rank-deficient design).
    if (is.null(fit)) {
      return(list(
        coef = h$coef, vcov = h$vcov,
        membership = rep(1L, length(ids)), K = 1L, lambda = sel$lambda
      ))
    }

    names(fit$coef) <- colnames(Wc)
    list(
      coef = fit$coef, vcov = fit$vcov,
      membership = memb, K = sel$K, lambda = sel$lambda
    )

  }, error = function(e) {
    message("  classifier_cre_gmm failed: ", e$message)
    NULL
  })
}

# =============================================================================
# Estimator 4: Oracle CRE-GMM (infeasible — uses TRUE group membership).
# Identical pooled common-slope specification to the classifier, but the
# group partition is the true one stored in the df$group column.
# =============================================================================

oracle_cre_gmm <- function(df) {
  tryCatch({
    d  <- .make_iv(df) %>% filter(!is.na(dy_lag2))
    Kk <- sort(unique(d$group))
    Wc <- cbind(ylag = d$ylag, x1 = d$x1, x2 = d$x2, z = d$z)
    Zc <- cbind(iv   = d$dy_lag2, x1 = d$x1, x2 = d$x2, z = d$z)
    for (k in Kk) {
      ind <- as.numeric(d$group == k)
      Wc  <- cbind(Wc, ind * d$xbar1, ind * d$xbar2)
      Zc  <- cbind(Zc, ind * d$xbar1, ind * d$xbar2)
      colnames(Wc)[(ncol(Wc) - 1):ncol(Wc)] <- paste0(c("xbar1_g", "xbar2_g"), k)
    }
    fit <- .gmm_fit(d$y, Wc, Zc)
    if (is.null(fit)) return(NULL)
    names(fit$coef) <- colnames(Wc)
    fit
  }, error = function(e) NULL)
}

# =============================================================================
# J-test for CRE homogeneity.
# H0: pi_g = pi for all g.  H1: pi_l != pi_k for some l != k.
# Compares restricted (homogeneous) and unrestricted (classifier) GMM criteria.
# Returns list(stat, df, pval) or NULL.
# =============================================================================

j_test_homogeneity <- function(df) {
  tryCatch({
    h  <- cre_gmm_homog(df)
    if (is.null(h)) return(NULL)
    cl <- classifier_cre_gmm(df)
    if (is.null(cl)) return(NULL)

    K <- cl$K
    if (K < 2L) {
      return(list(stat = 0, df = 0L, pval = 1))
    }

    # Restricted residuals (homogeneous).
    sse_r  <- sum(h$resid^2)
    sig2   <- mean(h$resid^2)

    # Unrestricted residuals (classifier).
    d          <- .make_iv(df) %>% filter(!is.na(dy_lag2))
    ids        <- sort(unique(df$id))
    memb_by_id <- setNames(cl$membership, ids)
    d$g        <- memb_by_id[as.character(d$id)]
    cf         <- cl$coef

    pred <- cf["ylag"] * d$ylag + cf["x1"] * d$x1 +
            cf["x2"]  * d$x2   + cf["z"]  * d$z
    for (k in sort(unique(d$g))) {
      ind <- as.numeric(d$g == k)
      c1  <- cf[paste0("xbar1_g", k)]
      c2  <- cf[paste0("xbar2_g", k)]
      pred <- pred + ind * (c1 * d$xbar1 + c2 * d$xbar2)
    }
    sse_u <- sum((d$y - pred)^2)

    stat  <- max((sse_r - sse_u) / sig2, 0)
    dfree <- (K - 1L) * 2L   # (G0 - 1) * p
    pval  <- stats::pchisq(stat, df = dfree, lower.tail = FALSE)

    list(stat = stat, df = dfree, pval = pval)

  }, error = function(e) NULL)
}

# =============================================================================
# Normalised mutual information between two integer partitions.
# =============================================================================

nmi <- function(a, b) {
  tryCatch({
    n   <- length(a)
    ta  <- table(a); tb <- table(b); tab <- table(a, b)
    Hab <- 0
    for (i in seq_along(ta)) {
      for (j in seq_along(tb)) {
        nij <- tab[i, j]
        if (nij == 0) next
        Hab <- Hab + (nij / n) * log((nij / n) / ((ta[i] / n) * (tb[j] / n)))
      }
    }
    Ha <- -sum((ta / n) * log(ta / n))
    Hb <- -sum((tb / n) * log(tb / n))
    if (Ha == 0 || Hb == 0) return(1)
    as.numeric(Hab / sqrt(Ha * Hb))
  }, error = function(e) NA_real_)
}

cat("02_estimators.R loaded. Run 03_simulation.R next.\n")

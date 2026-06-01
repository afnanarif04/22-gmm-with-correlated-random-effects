# Replication Package — Classifier CRE-GMM for Dynamic Panels

## About This Research

This repository contains the replication materials for a study on estimation
of dynamic panel data models that contain both time-varying regressors and
time-invariant covariates, when the individual effects are correlated with
the regressors through a latent group structure.

### The Problem

Level generalised method of moments estimation in dynamic panels requires
an instrument for the individual-specific effect. The correlated random
effects approach of Bontempi and Ditzen (2023) handles this by projecting
the individual effect onto the time-series means of the regressors (the
Mundlak 1978 correction). This works well when all units share the same
projection slope. When the projection slopes differ across latent groups
of units, imposing a common slope produces inconsistent estimates of the
coefficients on time-invariant covariates — the variables that level
estimation is specifically designed to identify.

### What This Paper Does

The paper proposes the **Classifier CRE-GMM estimator**, a three-step procedure:

1. **Step 1 (Initialisation):** Estimate the model under the homogeneous
   correlated random effects restriction (Bontempi & Ditzen 2023). Use the
   residuals to construct a preliminary unit-specific projection coefficient
   for each unit, via a local regression of the individual effect proxy on
   the within-unit means of the time-varying regressors.

2. **Step 2 (Classification):** Apply the pairwise adaptive group fused-LASSO
   of Mehrabani (2023) to the preliminary projection coefficients. This penalises
   differences between unit-pairs and shrinks them together, recovering the
   latent group partition. A Bayesian information criterion selects the number
   of groups automatically.

3. **Step 3 (Post-classification estimation):** Apply correlated random effects
   GMM-in-levels within each estimated group, using group-specific Mundlak
   projection slopes. This estimator is shown to achieve oracle efficiency:
   it has the same asymptotic distribution as the infeasible estimator that
   observes the true group partition, and strictly lower asymptotic variance
   for the time-invariant covariate coefficient than the homogeneous estimator.

A J-test for correlated random effects homogeneity with a chi-squared limiting
distribution and a consistent Bayesian information criterion complete the
inferential toolkit.

---

## Simulation Study Summary

The Monte Carlo study evaluates finite-sample performance across two designs
differing in the degree of group separation:

- **Design 1 (production function panel):** Two groups with projection vectors
  `pi_1 = (1.20, 0.80)` and `pi_2 = (0.40, 1.60)`, giving minimum separation
  `J_min ≈ 1.13` (large separation).

- **Design 2 (R&D panel):** Two groups with `pi_1 = (0.80, 1.00)` and
  `pi_2 = (1.20, 0.60)`, giving `J_min ≈ 0.57` (moderate separation).

The sample grid is `(N, T) ∈ {(50, 10), (100, 20), (200, 20), (100, 40), (200, 40)}`.
One thousand Monte Carlo replications are used for each cell.

**Key findings from the simulation:**

| Result | Value |
|--------|-------|
| Classifier RMSE for δ at N=200, T=40 (Design 1) | 0.145 |
| Homogeneous CRE-GMM RMSE for δ at same cell | 0.146 |
| Level GMM (no CRE correction) RMSE for δ at same cell | 1.779 |
| Reduction in RMSE: Classifier vs Level GMM | > 90 per cent |
| J-test empirical size at N=200, T=20 (nominal 5%) | 0.05 (Design 1H), 0.03 (Design 2H) |
| J-test power at γ_sep = 1.0, N=100, T=40 | 0.70 |
| BIC selects K=2 (correct) at N=200, T=20 — Design 1 | 43% of replications |
| NMI of recovered partition at same cell | 0.60 |

---

## Empirical Application Summary

The empirical application estimates a log-linearised Cobb-Douglas production
function using Penn World Table 11.0 (Feenstra, Inklaar & Timmer 2015):

| Dimension | Value |
|-----------|-------|
| Countries (N) | 143 |
| Years (T) | 29 (1995–2023) |
| Observations | 4,147 (fully balanced) |

**Dependent variable:** log real output per capita (`ly = log(rgdpo / pop)`)

**Time-varying regressors:** log capital stock per capita (`lk = log(rnna / pop)`)
and log human capital index (`lh = log(hc)`)

**Time-invariant regressor:** log land area (`larea = log(land_area_sqkm)`)

**Preliminary tests:**

| Test | Result |
|------|--------|
| Pooled AR(1) persistence | ρ̂ = 0.9955 (near unit root — justifies GMM-lev) |
| Pesaran (2021) CD statistic | 59.09 (strong cross-sectional dependence) |
| Cross-country slope SD for lk | 2.32 (large heterogeneity) |
| Cross-country slope SD for lh | 4.56 (large heterogeneity) |

**Estimation results:**

The information criterion selects a single group (K̂ = 1) for this panel,
and the J-test for correlated random effects homogeneity does not reject
at conventional significance levels. The homogeneous CRE-GMM estimator
therefore characterises the production function for this sample. The
near-unit-root persistence motivates the levels-based approach throughout.

---

## Repository Structure

```
README.md                  — this file
LICENSE                    — MIT licence
data/
  pwt_empirical_panel.xlsx — pre-processed estimation panel (N=143, T=29)
  README_data.md           — variable definitions and download instructions
code/
  00_setup.R               — libraries, constants, simulation grid
  01_dgp_functions.R       — data-generating processes (DGPs 1, 2, 1H, 2H, 1P)
  02_estimators.R          — all estimators: GMM-lev, CRE-GMM, Classifier, Oracle,
                             PAGFL, BIC, J-test, NMI
  03_simulation.R          — main Monte Carlo loop (Tables 1–3)
  04_summarise.R           — aggregates simulation output to ALL_RESULTS.xlsx
  05_empirical_data.R      — loads pwt_empirical_panel.xlsx → panel_pwt_clean.rds
  06_preliminary_tests.R   — AR(1) persistence, Pesaran CD, slope heterogeneity
  07_empirical_estimation.R — four estimators on PWT panel, Table 4
output/
  .gitkeep                 — placeholder; populated when code runs
```

---

## Software Requirements

- R version 4.3.0 or later (tested on R 4.5.2, Windows 11)
- RStudio (recommended — scripts use `rstudioapi` for portable `setwd`)

### Required R Packages

```r
install.packages(c(
  "tidyverse",   # dplyr, tidyr, tibble, purrr, readr, ggplot2
  "MASS",        # mvrnorm, ginv
  "openxlsx",    # createWorkbook, writeData, saveWorkbook, loadWorkbook
  "readxl"       # read_excel (empirical scripts only)
))
```

| Package   | Version tested |
|-----------|----------------|
| tidyverse | 2.0.0          |
| MASS      | 7.3-60         |
| openxlsx  | 4.2.5          |
| readxl    | 1.4.3          |

---

## Data

The file `data/pwt_empirical_panel.xlsx` contains the pre-processed estimation
panel. It was constructed from publicly available sources. See
`data/README_data.md` for the full variable dictionary, construction steps,
and exact download instructions.

| Data | Source | Direct URL |
|------|--------|-----------|
| Real GDP, capital, human capital | Penn World Table 11.0 | https://www.rug.nl/ggdc/productivity/pwt/ |
| Land area (sq. km) | World Bank AG.LND.TOTL.K2 | https://data.worldbank.org/indicator/AG.LND.TOTL.K2 |

---

## Replication Instructions

### Simulation (Tables 1–3)

Run the following scripts **in order** from RStudio (open each file and
press **Source**, or run `source("script.R")` from the console):

```
code/00_setup.R              → loads libraries and sets constants
code/01_dgp_functions.R      → defines all DGP functions
code/02_estimators.R         → defines all estimator functions
code/03_simulation.R         → runs 1,000 replications × 5 cells × 2 DGPs
                               (all 3 blocks: estimation, J-test, selection)
code/04_summarise.R          → reads the .rds output and writes
                               output/ALL_RESULTS.xlsx (Sheets 1–3)
```

**Estimated runtime:** 8–10 hours for the full 1,000-rep production run
on a single core. To run a quick check, open `00_setup.R` and temporarily
set `N_REP <- 20L` before sourcing.

The random seed is fixed via `SEED_BASE <- 20260529L`. Results are
reproducible across runs on the same platform.

### Empirical Estimation (Table 4)

After the simulation is complete (or independently):

```
code/05_empirical_data.R     → loads data/pwt_empirical_panel.xlsx
                               → writes output/panel_pwt_clean.rds
code/06_preliminary_tests.R  → AR(1) persistence, Pesaran CD, slope SD
                               → writes output/pretests_empirical.csv
code/07_empirical_estimation.R → runs GMM-lev, CRE-GMM, Classifier, J-test
                               → appends Sheet 4 to output/ALL_RESULTS.xlsx
                               → writes output/classifier_groups_empirical.csv
```

**Estimated runtime:** approximately 10–20 minutes.

### All Output Files

After running all scripts, the `output/` folder will contain:

| File | Contents |
|------|----------|
| `sim_estimation.rds` | Raw estimation results (Block 1) |
| `sim_jtest.rds` | Raw J-test results (Block 2) |
| `sim_selection.rds` | Raw selection/NMI results (Block 3) |
| `ALL_RESULTS.xlsx` | Tables 1–4 in the paper (all four sheets) |
| `panel_pwt_clean.rds` | Cleaned empirical panel |
| `pretests_empirical.csv` | Preliminary test statistics |
| `classifier_groups_empirical.csv` | Country-level group assignments |

---

## Notes

- Results may vary slightly across platforms and R versions due to
  floating-point arithmetic and random number generator implementation.
- The BIC penalty constant `BIC_C = 0.30` is calibrated on the simulation
  designs. Underfitting probability falls as N and T grow. Overfitting is
  penalised by the `rho_NT` rate.
- For the N=200 cells, the pairwise PAGFL graph is restricted to the 10
  nearest neighbours in projection-coefficient space, reducing the pair
  count from O(N²) to O(kN) while preserving classification accuracy.
- The empirical panel is fully balanced: every country appears in every year.

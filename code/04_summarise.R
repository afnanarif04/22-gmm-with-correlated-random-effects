# =============================================================================
# 04_summarise.R
# Aggregate the three raw simulation .rds files into the manuscript tables
# and write the SINGLE Excel workbook ALL_RESULTS.xlsx (Rule 5).
#
# Reads:   sim_estimation.rds, sim_jtest.rds, sim_selection.rds
# Writes:  ALL_RESULTS.xlsx (sheets 1 to 3)
#
# This is the ONLY script that writes Excel (Rule 5).
# After running 07_empirical_estimation.R the fourth sheet is appended.
# =============================================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dir.create("../output", showWarnings = FALSE)  # ensure output folder exists
source("00_setup.R")

# Safe RDS reader — returns empty tibble if file not yet produced.
rd_rds <- function(f) if (file.exists(f)) readRDS(f) else tibble()

est <- rd_rds("../output/sim_estimation.rds")
jt  <- rd_rds("../output/sim_jtest.rds")
sel <- rd_rds("../output/sim_selection.rds")

# ---- Table 1: Bias and RMSE (3 decimal places, Mehrabani 2023 convention) ---

build_table1 <- function(est) {
  if (nrow(est) == 0) return(tibble())
  est %>%
    pivot_longer(c(rho, b1, b2, delta), names_to = "param", values_to = "estv") %>%
    mutate(truev = dplyr::case_when(
      param == "rho"   ~ true_rho,
      param == "b1"    ~ true_b1,
      param == "b2"    ~ true_b2,
      param == "delta" ~ true_delta
    )) %>%
    group_by(dgp, N, T, est, param) %>%
    summarise(
      bias  = round(mean(estv - truev, na.rm = TRUE), 3),
      rmse  = round(sqrt(mean((estv - truev)^2, na.rm = TRUE)), 3),
      n_ok  = sum(!is.na(estv)),
      .groups = "drop"
    ) %>%
    arrange(
      dgp, N, T,
      factor(est,   levels = c("Classifier", "Homogeneous", "GMMlev", "Oracle")),
      factor(param, levels = c("rho", "b1", "b2", "delta"))
    )
}

# ---- Table 2: J-test size and power (2 decimal places) ---------------------

build_table2 <- function(jt) {
  if (nrow(jt) == 0) return(tibble())
  jt %>%
    mutate(reject_rate = round(reject_rate, 2)) %>%
    arrange(kind, dgp, N, T, gamma_sep)
}

# ---- Table 3: Group selection and NMI (2 decimal places) -------------------

build_table3 <- function(sel) {
  if (nrow(sel) == 0) return(tibble())
  sel %>%
    mutate(across(c(pK1, pK2, pK3, nmi), ~round(.x, 2))) %>%
    arrange(dgp, N, T)
}

t1 <- build_table1(est)
t2 <- build_table2(jt)
t3 <- build_table3(sel)

# ---- Write ALL_RESULTS.xlsx --------------------------------------------------

wb     <- createWorkbook()
hdr_st <- createStyle(textDecoration = "bold", halign = "center",
                      fgFill = "#D9E1F2", border = "Bottom")
note_st <- createStyle(fontSize = 9, fontColour = "#555555",
                       textDecoration = "italic")

add_sheet <- function(wb, nm, data, note = NULL) {
  addWorksheet(wb, nm)
  if (nrow(data) == 0) {
    writeData(wb, nm, "No data. Run 03_simulation.R first.")
    return(invisible())
  }
  writeData(wb, nm, data, startRow = 1, headerStyle = hdr_st)
  setColWidths(wb, nm, cols = seq_along(data), widths = "auto")
  if (!is.null(note)) {
    nr <- nrow(data) + 3
    writeData(wb, nm, note, startRow = nr, startCol = 1)
    addStyle(wb, nm, note_st, rows = nr, cols = 1)
  }
}

add_sheet(wb, "1.Table1 Bias-RMSE", t1,
  "Table 1. Bias and RMSE of the four estimators. 1000 replications. Bias = mean(estimate - true); RMSE = square root of mean squared error.")
add_sheet(wb, "2.Table2 J-test", t2,
  "Table 2. J-test for CRE homogeneity at 5 per cent nominal level. kind=size under H0; kind=power under H1; kind=curve is the power curve at N=100, T=20.")
add_sheet(wb, "3.Table3 Selection", t3,
  "Table 3. Group selection frequency and normalised mutual information. NMI = 1.00 indicates perfect classification.")

saveWorkbook(wb, "../output/ALL_RESULTS.xlsx", overwrite = TRUE)

cat("ALL_RESULTS.xlsx written with Tables 1-3.\n")
cat("Next: run 05_empirical_data.R (after placing pwt1001.xlsx in this folder).\n")

#!/usr/bin/env Rscript

# =====================================================================
# run_all.R  --  master reproduction driver
# =====================================================================
# Generating-Element Maximum Entropy for Non-Gaussian Uncertainty
# Evaluation (Zabolotnii, 2026).
#
# Runs every verification experiment in order and writes all figures and
# per-script console logs into outputs/.  Each script is sourced in its own
# fresh environment so global state never leaks between experiments.  Every
# script sets its own seed (20260612 for the headline single-realization
# runs, 1..20 for the seed-replication studies), so results are fully
# deterministic.
#
# Usage:
#   Rscript run_all.R
#
# Output:
#   outputs/                          figures (PDF)
#   outputs/logs/<script>.log         full console transcript per script
#   outputs/run_all_summary.txt       timing + status summary
# =====================================================================

# ---- locate the repository root (dir containing this script) ----
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(file_arg) == 1 && nzchar(file_arg)) {
  ROOT <- normalizePath(dirname(file_arg))
} else {
  ROOT <- normalizePath(getwd())   # fallback: assume cwd is repo root
}

R_DIR   <- file.path(ROOT, "R")
OUT_DIR <- file.path(ROOT, "outputs")
LOG_DIR <- file.path(OUT_DIR, "logs")
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

# Figure-writing scripts read this; force all figures into outputs/.
Sys.setenv(GENELEMENT_OUT_DIR = OUT_DIR)

# ---- experiment order ----
# parity_matched_patp.R is run LAST on purpose: it (re)writes the final
# paper variants of fig_cauchy_logpdf.pdf and fig_mixture_fits.pdf (PM-PATP
# + T-MaxEnt), superseding the intermediate versions written by the two
# ablation scripts.
SCRIPTS <- c(
  "patp_maxent_simulation.R",   # CORE solver + Exp 1 (Cauchy) & Exp 2 (mixture), Form-B
  "ablation_cauchy_metrics.R",  # Exp 1 extended per-alpha metrics (Tab 1 discussion)
  "ablation_mixture.R",         # Exp 2 baselines + T-MaxEnt freq ablation (Tab 2, Tab 3)
  "ablation_ga_noise.R",        # Exp 4 GA MaxEnt-vs-MC fitness (Tab mv)
  "extra_ablations.R",          # linear-map control + 20-seed replication
  "opmm_alpha_cauchy.R",        # oPMMalpha variance-optimal selection, Cauchy (revision: reviewer Q1/Q6)
  "opmm_alpha_mixture.R",       # oPMMalpha + CV log-score, mixture control (revision: reviewer Q1/Q6)
  "free_exponent_baseline.R",   # FM-MEM free-exponent baseline (design justification)
  "log_rational_element.R",     # Exp 3 generating-element comparison (Tab genelement, fig)
  "log_rational_unbounded.R",   # unbounded closed-form matched fit (revision: reviewer Q3/W3)
  "sensitivity_L_N.R",          # truncation L x sample N sensitivity grid (revision: reviewer Q2/W6)
  "fm_mem_ga.R",                # FM-MEM GA exponent search vs alpha-scan (revision: reviewer Q4/W5)
  "mv_qmc_baseline.R",          # M&V randomized-QMC baseline + per-seed designs (revision: reviewer Q5/W7)
  "tmaxent_auto_freq.R",        # automated ECF-noise-floor T-MaxEnt frequency rule (revision: reviewer Q6)
  "parity_matched_patp.R"       # PM-PATP final figures (fig:cauchy, fig:mixture) -- LAST
)

cat("=====================================================================\n")
cat(" GENERATING-ELEMENT MaxEnt -- full reproduction run\n")
cat(sprintf(" R: %s\n", R.version.string))
cat(sprintf(" repo root : %s\n", ROOT))
cat(sprintf(" outputs   : %s\n", OUT_DIR))
cat("=====================================================================\n\n")

run_one <- function(script) {
  path <- file.path(R_DIR, script)
  log  <- file.path(LOG_DIR, sub("\\.R$", ".log", script))
  cat(sprintf("[RUN ] %-28s -> logs/%s\n", script, basename(log)))
  t0 <- Sys.time()
  con <- file(log, open = "wt")
  sink(con); sink(con, type = "message")
  status <- tryCatch({
    sys.source(path, envir = new.env(parent = globalenv()))
    "OK"
  }, error = function(e) paste0("ERROR: ", conditionMessage(e)))
  sink(type = "message"); sink(); close(con)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("[%s] %-28s  %.1f s\n",
              if (status == "OK") "DONE" else "FAIL", script, secs))
  data.frame(script = script, status = status, seconds = secs,
             stringsAsFactors = FALSE)
}

t_all <- Sys.time()
summary_df <- do.call(rbind, lapply(SCRIPTS, run_one))
total_secs <- as.numeric(difftime(Sys.time(), t_all, units = "secs"))

cat("\n---------------------------------------------------------------------\n")
cat(" SUMMARY\n")
cat("---------------------------------------------------------------------\n")
for (i in seq_len(nrow(summary_df))) {
  cat(sprintf("  %-28s %-6s %6.1f s\n",
              summary_df$script[i], summary_df$status[i], summary_df$seconds[i]))
}
cat(sprintf("  %-28s %-6s %6.1f s\n", "TOTAL", "", total_secs))

figs <- list.files(OUT_DIR, pattern = "\\.pdf$", full.names = FALSE)
cat(sprintf("\n Figures written to outputs/: %s\n",
            if (length(figs)) paste(figs, collapse = ", ") else "(none)"))

summary_path <- file.path(OUT_DIR, "run_all_summary.txt")
writeLines(c(
  sprintf("Reproduction run: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("R: %s", R.version.string),
  "",
  sprintf("%-28s %-8s %s", "script", "status", "seconds"),
  apply(summary_df, 1, function(r)
    sprintf("%-28s %-8s %.1f", r[["script"]], r[["status"]], as.numeric(r[["seconds"]]))),
  sprintf("%-28s %-8s %.1f", "TOTAL", "", total_secs),
  "",
  sprintf("figures: %s", paste(figs, collapse = ", "))
), summary_path)

cat(sprintf(" Summary written to %s\n", summary_path))

if (any(summary_df$status != "OK")) {
  cat("\n*** One or more scripts FAILED -- see logs above. ***\n")
  quit(status = 1)
}
cat("\nAll experiments completed successfully.\n")

#!/usr/bin/env Rscript

# ablation_cauchy_metrics.R
# Extended ablation of Experiment 1 (standard Cauchy) from patp_maxent_simulation.R
# for the Q1 revision: per-alpha convergence status, Hessian condition number,
# dual potential, PATP exponents, and accuracy metrics of the fitted MaxEnt
# density vs the TRUE standard Cauchy law:
#   (a) Q75/Q90/Q95/Q99 quantiles from the fitted CDF on the grid + % errors,
#   (b) central-body density MSE on [-10, 10],
#   (c) Kolmogorov-Smirnov distance sup|F_fit - F_true| on grid in [-45, 45],
#   (d) true Cauchy mass outside [-50, 50] (constant reference).
# Figure: log10 density vs x on [-20, 20] for true Cauchy and fitted
# alpha = 0, 0.2, 0.5.
#
# Base R only. Solver functions copied VERBATIM from patp_maxent_simulation.R.

# -------------------------------------------------------------------------
# 1. CORE PATP-MAXENT IMPLEMENTATION (verbatim from patp_maxent_simulation.R)
# -------------------------------------------------------------------------

patp_power <- function(i, alpha) {
  # Exponent map p_i(alpha)
  A <- 1.0 / i
  B <- 4.0 - i - 3.0 / i
  C <- 2.0 * i - 4.0 + 2.0 / i
  return(A + B * alpha + C * alpha^2)
}

patp_basis <- function(x, n, alpha) {
  # Computes the PATP basis functions for a vector x
  # Returns a matrix of size length(x) x n
  # Column 1 is x (linear constraint)
  # Columns 2 to n are sign(x) * |x|^p_i(alpha) for i = 2..n
  M <- matrix(0, nrow = length(x), ncol = n)
  M[, 1] <- x
  if (n > 1) {
    for (i in 2:n) {
      p <- patp_power(i, alpha)
      # Avoid complex numbers for negative x by using знакозберігаючий степінь
      M[, i] <- sign(x) * (abs(x)^p)
    }
  }
  return(M)
}

solve_maxent <- function(Phi, target_moments, dx, max_iter = 100, tol = 1e-6) {
  # Newton-Raphson solver with backtracking line search for dual MaxEnt problem
  n <- ncol(Phi)
  lambdas <- rep(0, n)

  for (iter in 1:max_iter) {
    unnorm_pdf <- exp(as.numeric(Phi %*% lambdas))
    if (any(is.infinite(unnorm_pdf)) || any(is.na(unnorm_pdf))) {
      return(NULL)
    }

    Z <- sum(unnorm_pdf) * dx
    if (Z == 0 || is.na(Z)) return(NULL)

    pdf <- unnorm_pdf / Z
    fitted_moments <- as.numeric(t(Phi) %*% pdf * dx)
    grad <- fitted_moments - target_moments

    # Check convergence
    if (sqrt(sum(grad^2)) < tol) {
      H <- t(Phi) %*% (Phi * pdf) * dx - fitted_moments %*% t(fitted_moments)
      diag(H) <- diag(H) + 1e-8 # Small regularization for stability
      cond_num <- kappa(H, exact = TRUE)
      return(list(lambdas = lambdas, Z = Z, grad = grad, H = H, cond_num = cond_num, converged = TRUE))
    }

    H <- t(Phi) %*% (Phi * pdf) * dx - fitted_moments %*% t(fitted_moments)
    diag(H) <- diag(H) + 1e-8

    step <- tryCatch(solve(H, grad), error = function(e) {
      tryCatch(qr.solve(H, grad), error = function(e2) NULL)
    })
    if (is.null(step)) return(NULL)

    # Backtracking Line Search
    alpha_step <- 1.0
    converged_step <- FALSE
    for (ls in 1:10) {
      new_lambdas <- lambdas - alpha_step * step
      new_unnorm <- exp(as.numeric(Phi %*% new_lambdas))
      if (!any(is.infinite(new_unnorm)) && !any(is.na(new_unnorm))) {
        new_Z <- sum(new_unnorm) * dx
        if (new_Z > 0 && !is.na(new_Z)) {
          # Check potential: log(Z) - sum(lambdas * target)
          old_pot <- log(Z) - sum(lambdas * target_moments)
          new_pot <- log(new_Z) - sum(new_lambdas * target_moments)
          if (new_pot <= old_pot + 1e-4) {
            lambdas <- new_lambdas
            converged_step <- TRUE
            break
          }
        }
      }
      alpha_step <- alpha_step * 0.5
    }

    if (!converged_step) {
      lambdas <- lambdas - 0.1 * step # Small step fallback
    }
  }

  # Return state even if not fully converged
  unnorm_pdf <- exp(as.numeric(Phi %*% lambdas))
  Z <- sum(unnorm_pdf) * dx
  pdf <- unnorm_pdf / Z
  fitted_moments <- as.numeric(t(Phi) %*% pdf * dx)
  grad <- fitted_moments - target_moments
  H <- t(Phi) %*% (Phi * pdf) * dx - fitted_moments %*% t(fitted_moments)
  diag(H) <- diag(H) + 1e-8
  cond_num <- kappa(H, exact = TRUE)
  return(list(lambdas = lambdas, Z = Z, grad = grad, H = H, cond_num = cond_num, converged = FALSE))
}

# -------------------------------------------------------------------------
# 2. EXPERIMENT SETUP (identical to Experiment 1 of patp_maxent_simulation.R)
# -------------------------------------------------------------------------

set.seed(20260612)
N <- 1000
data_cauchy <- rcauchy(N, location = 0, scale = 1)

L <- 50
x_grid <- seq(-L, L, length.out = 2000)
dx <- x_grid[2] - x_grid[1]
n_moments <- 4
alpha_grid <- seq(0.0, 1.0, by = 0.1)

# True Cauchy references
q_levels <- c(0.75, 0.90, 0.95, 0.99)
q_true <- qcauchy(q_levels, 0, 1)
mass_outside <- 2 * pcauchy(-50, 0, 1)

cat("=========================================================\n")
cat("ABLATION: STANDARD CAUCHY, PATP-MAXENT, EXTENDED METRICS\n")
cat("=========================================================\n")
cat(sprintf("Seed = 20260612, N = %d, grid = [-%g, %g] with %d points, dx = %.6f\n",
            N, L, L, length(x_grid), dx))
cat(sprintf("Constraints: %d PATP moments; alpha grid = %s\n",
            n_moments, paste(sprintf("%.1f", alpha_grid), collapse = ", ")))
cat(sprintf("\nTrue quantiles (qcauchy): Q75 = %.6f, Q90 = %.6f, Q95 = %.6f, Q99 = %.6f\n",
            q_true[1], q_true[2], q_true[3], q_true[4]))
cat(sprintf("(d) True Cauchy mass OUTSIDE [-50, 50] = 2*pcauchy(-50) = %.8f (%.5f%%)\n",
            mass_outside, 100 * mass_outside))
cat("    NOTE: the fitted MaxEnt density is normalized to total mass 1 on the\n")
cat("    truncated grid, so this is the irreducible mass the truncation discards.\n")

# Subset indices for metrics
idx_body <- which(x_grid >= -10 & x_grid <= 10)   # central body [-10, 10]
idx_ks   <- which(x_grid >= -45 & x_grid <= 45)   # KS restricted to [-45, 45]
pdf_true_body <- dcauchy(x_grid[idx_body], 0, 1)
cdf_true_ks   <- pcauchy(x_grid[idx_ks], 0, 1)

# -------------------------------------------------------------------------
# 3. ALPHA SWEEP WITH FULL METRICS
# -------------------------------------------------------------------------

results <- data.frame()
fitted_pdfs <- list()   # keep fitted pdf per alpha for the figure

for (alpha in alpha_grid) {
  a_key <- sprintf("%.1f", alpha)
  p2 <- patp_power(2, alpha)
  p3 <- patp_power(3, alpha)
  p4 <- patp_power(4, alpha)

  cat(sprintf("\n--- alpha = %.1f ---\n", alpha))
  cat(sprintf("PATP exponents: p2 = %.6f, p3 = %.6f, p4 = %.6f\n", p2, p3, p4))

  # Empirical PATP moments from the SAME sample
  Phi_data <- patp_basis(data_cauchy, n_moments, alpha)
  target_moments <- colMeans(Phi_data)
  cat(sprintf("Empirical target moments: %s\n",
              paste(sprintf("%.6g", target_moments), collapse = ", ")))

  Phi_grid <- patp_basis(x_grid, n_moments, alpha)
  fit <- solve_maxent(Phi_grid, target_moments, dx)

  if (is.null(fit)) {
    cat("Solver returned NULL (overflow/singular Hessian) -- no fit available.\n")
    results <- rbind(results, data.frame(
      alpha = alpha, status = "solver-NULL",
      p2 = p2, p3 = p3, p4 = p4,
      cond_num = NA, potential = NA,
      q75 = NA, q90 = NA, q95 = NA, q99 = NA,
      err_q75_pct = NA, err_q90_pct = NA, err_q95_pct = NA, err_q99_pct = NA,
      body_mse = NA, ks_dist = NA
    ))
    next
  }

  pdf_fitted <- exp(as.numeric(Phi_grid %*% fit$lambdas)) / fit$Z
  cdf_fitted <- cumsum(pdf_fitted) * dx
  potential <- log(fit$Z) - sum(fit$lambdas * target_moments)

  # (a) Quantiles from the fitted CDF on the grid (nearest grid point, as in
  #     the original Experiment 1) + signed relative errors vs qcauchy
  q_est <- sapply(q_levels, function(p) x_grid[which.min(abs(cdf_fitted - p))])
  err_pct <- (q_est - q_true) / q_true * 100

  # (b) Central-body density MSE on [-10, 10]
  body_mse <- mean((pdf_fitted[idx_body] - pdf_true_body)^2)

  # (c) Kolmogorov-Smirnov distance on the grid restricted to [-45, 45]
  ks_dist <- max(abs(cdf_fitted[idx_ks] - cdf_true_ks))

  status <- if (fit$converged) "TRUE" else "FALSE"
  cat(sprintf("Converged = %s, kappa(H) exact = %.6e, dual potential = %.6f\n",
              status, fit$cond_num, potential))
  cat(sprintf("Final gradient norm = %.3e\n", sqrt(sum(fit$grad^2))))
  cat(sprintf("Lambdas = %s\n", paste(sprintf("%.6g", fit$lambdas), collapse = ", ")))
  cat(sprintf("(a) Q75 = %.4f (true %.4f, err %+.2f%%) | Q90 = %.4f (true %.4f, err %+.2f%%)\n",
              q_est[1], q_true[1], err_pct[1], q_est[2], q_true[2], err_pct[2]))
  cat(sprintf("    Q95 = %.4f (true %.4f, err %+.2f%%) | Q99 = %.4f (true %.4f, err %+.2f%%)\n",
              q_est[3], q_true[3], err_pct[3], q_est[4], q_true[4], err_pct[4]))
  cat(sprintf("(b) Central-body density MSE on [-10,10] = %.6e\n", body_mse))
  cat(sprintf("(c) KS distance sup|F_fit - F_true| on [-45,45] = %.6f\n", ks_dist))

  results <- rbind(results, data.frame(
    alpha = alpha, status = status,
    p2 = p2, p3 = p3, p4 = p4,
    cond_num = fit$cond_num, potential = potential,
    q75 = q_est[1], q90 = q_est[2], q95 = q_est[3], q99 = q_est[4],
    err_q75_pct = err_pct[1], err_q90_pct = err_pct[2],
    err_q95_pct = err_pct[3], err_q99_pct = err_pct[4],
    body_mse = body_mse, ks_dist = ks_dist
  ))
  fitted_pdfs[[a_key]] <- pdf_fitted
}

# -------------------------------------------------------------------------
# 4. SUMMARY TABLES
# -------------------------------------------------------------------------

cat("\n=========================================================\n")
cat("SUMMARY TABLE (full precision)\n")
cat("=========================================================\n")
print(format(results, digits = 8), row.names = FALSE)

cat("\nCONVERGENCE REGION:\n")
conv_alphas <- results$alpha[results$status == "TRUE"]
nonconv_alphas <- results$alpha[results$status == "FALSE"]
null_alphas <- results$alpha[results$status == "solver-NULL"]
cat(sprintf("  Converged (TRUE):        %s\n",
            if (length(conv_alphas)) paste(sprintf("%.1f", conv_alphas), collapse = ", ") else "none"))
cat(sprintf("  Not converged (FALSE):   %s\n",
            if (length(nonconv_alphas)) paste(sprintf("%.1f", nonconv_alphas), collapse = ", ") else "none"))
cat(sprintf("  Solver NULL (overflow):  %s\n",
            if (length(null_alphas)) paste(sprintf("%.1f", null_alphas), collapse = ", ") else "none"))

ok <- results[results$status == "TRUE", ]
if (nrow(ok) > 0) {
  cat("\nSELECTION CHECK (among converged alphas only):\n")
  cat(sprintf("  argmin dual potential : alpha = %.1f (potential = %.6f)\n",
              ok$alpha[which.min(ok$potential)], min(ok$potential)))
  cat(sprintf("  argmin body MSE       : alpha = %.1f (MSE = %.6e)\n",
              ok$alpha[which.min(ok$body_mse)], min(ok$body_mse)))
  cat(sprintf("  argmin KS distance    : alpha = %.1f (KS = %.6f)\n",
              ok$alpha[which.min(ok$ks_dist)], min(ok$ks_dist)))
  cat(sprintf("  argmin |Q95 error|    : alpha = %.1f (err = %+.2f%%)\n",
              ok$alpha[which.min(abs(ok$err_q95_pct))], ok$err_q95_pct[which.min(abs(ok$err_q95_pct))]))
  cat(sprintf("  argmin |Q99 error|    : alpha = %.1f (err = %+.2f%%)\n",
              ok$alpha[which.min(abs(ok$err_q99_pct))], ok$err_q99_pct[which.min(abs(ok$err_q99_pct))]))
}

# -------------------------------------------------------------------------
# 5. FIGURE 1: log10 density on [-20, 20]
# -------------------------------------------------------------------------

fig_dir <- Sys.getenv("GENELEMENT_OUT_DIR", unset = "outputs")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
fig_path <- file.path(fig_dir, "fig_cauchy_logpdf.pdf")

idx_plot <- which(x_grid >= -20 & x_grid <= 20)
xp <- x_grid[idx_plot]

plot_specs <- list(
  list(key = "0.0", label = expression(PATP ~ alpha == 0),   col = "#D55E00", lty = "dashed"),
  list(key = "0.2", label = expression(PATP ~ alpha == 0.2), col = "#0072B2", lty = "dotdash"),
  list(key = "0.5", label = expression(PATP ~ alpha == 0.5), col = "#009E73", lty = "dotted")
)

pdf(fig_path, width = 7, height = 4.5)
par(mar = c(4.5, 4.5, 1, 1))

ylist <- list(log10(dcauchy(xp, 0, 1)))
for (s in plot_specs) {
  if (!is.null(fitted_pdfs[[s$key]])) ylist[[length(ylist) + 1]] <- log10(fitted_pdfs[[s$key]][idx_plot])
}
yrange <- range(unlist(lapply(ylist, function(y) y[is.finite(y)])))

plot(xp, log10(dcauchy(xp, 0, 1)), type = "l", lwd = 2.5, col = "black", lty = "solid",
     xlab = "x (dimensionless)", ylab = expression(log[10] ~ "density" ~ f(x)),
     ylim = yrange, xlim = c(-20, 20))

leg_labels <- c(expression("True standard Cauchy"))
leg_cols <- c("black"); leg_ltys <- c("solid")
for (s in plot_specs) {
  if (!is.null(fitted_pdfs[[s$key]])) {
    lines(xp, log10(fitted_pdfs[[s$key]][idx_plot]), lwd = 2, col = s$col, lty = s$lty)
    leg_labels <- c(leg_labels, s$label)
    leg_cols <- c(leg_cols, s$col)
    leg_ltys <- c(leg_ltys, s$lty)
  } else {
    cat(sprintf("Figure note: alpha = %s has no fit (solver NULL), omitted from plot.\n", s$key))
  }
}
legend("topright", legend = leg_labels, col = leg_cols, lty = leg_ltys,
       lwd = 2, bty = "n", cex = 0.9)
dev.off()

cat(sprintf("\nFigure written: %s\n", fig_path))
cat("\nAblation complete.\n")

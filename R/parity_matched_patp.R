#!/usr/bin/env Rscript

# Parity-Matched PATP-MaxEnt (PM-PATP) verification.
# Motivation: the original Form-B basis {x, sign(x)|x|^{p_i}} is all-odd, so the
# fitted log-density is odd and f(x)f(-x) = 1/Z^2 identically — it cannot
# represent symmetric densities. The parity-matched variant assigns
#   phi_i(x; alpha) = |x|^{p_i(alpha)}        for even i,
#   phi_i(x; alpha) = sign(x)|x|^{p_i(alpha)} for odd i,
# which recovers the TRUE monomials x^i exactly at alpha = 1.
# Also tests T-MaxEnt on the Cauchy (CF constraints exist even when moments don't).
# Solver copied VERBATIM from patp_maxent_simulation.R.

patp_power <- function(i, alpha) {
  A <- 1.0 / i
  B <- 4.0 - i - 3.0 / i
  C <- 2.0 * i - 4.0 + 2.0 / i
  return(A + B * alpha + C * alpha^2)
}

patp_basis_pm <- function(x, n, alpha) {
  # Parity-matched PATP basis
  M <- matrix(0, nrow = length(x), ncol = n)
  M[, 1] <- x
  if (n > 1) {
    for (i in 2:n) {
      p <- patp_power(i, alpha)
      if (i %% 2 == 0) {
        M[, i] <- abs(x)^p
      } else {
        M[, i] <- sign(x) * (abs(x)^p)
      }
    }
  }
  return(M)
}

trig_basis <- function(x, S, p) {
  M <- matrix(0, nrow = length(x), ncol = 2 * S)
  for (r in 1:S) {
    M[, 2 * r - 1] <- cos(r * p * x)
    M[, 2 * r]     <- sin(r * p * x)
  }
  return(M)
}

solve_maxent <- function(Phi, target_moments, dx, max_iter = 100, tol = 1e-6) {
  n <- ncol(Phi)
  lambdas <- rep(0, n)
  for (iter in 1:max_iter) {
    unnorm_pdf <- exp(as.numeric(Phi %*% lambdas))
    if (any(is.infinite(unnorm_pdf)) || any(is.na(unnorm_pdf))) return(NULL)
    Z <- sum(unnorm_pdf) * dx
    if (Z == 0 || is.na(Z)) return(NULL)
    pdf <- unnorm_pdf / Z
    fitted_moments <- as.numeric(t(Phi) %*% pdf * dx)
    grad <- fitted_moments - target_moments
    if (sqrt(sum(grad^2)) < tol) {
      H <- t(Phi) %*% (Phi * pdf) * dx - fitted_moments %*% t(fitted_moments)
      diag(H) <- diag(H) + 1e-8
      cond_num <- kappa(H, exact = TRUE)
      return(list(lambdas = lambdas, Z = Z, grad = grad, H = H, cond_num = cond_num, converged = TRUE))
    }
    H <- t(Phi) %*% (Phi * pdf) * dx - fitted_moments %*% t(fitted_moments)
    diag(H) <- diag(H) + 1e-8
    step <- tryCatch(solve(H, grad), error = function(e) {
      tryCatch(qr.solve(H, grad), error = function(e2) NULL)
    })
    if (is.null(step)) return(NULL)
    alpha_step <- 1.0
    converged_step <- FALSE
    for (ls in 1:10) {
      new_lambdas <- lambdas - alpha_step * step
      new_unnorm <- exp(as.numeric(Phi %*% new_lambdas))
      if (!any(is.infinite(new_unnorm)) && !any(is.na(new_unnorm))) {
        new_Z <- sum(new_unnorm) * dx
        if (new_Z > 0 && !is.na(new_Z)) {
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
    if (!converged_step) lambdas <- lambdas - 0.1 * step
  }
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

fit_metrics_cauchy <- function(Phi_grid, fit, x_grid, dx) {
  pdf_fitted <- exp(as.numeric(Phi_grid %*% fit$lambdas)) / fit$Z
  cdf_fitted <- cumsum(pdf_fitted) * dx
  qs <- c(0.75, 0.90, 0.95, 0.99)
  q_est <- sapply(qs, function(q) x_grid[which.min(abs(cdf_fitted - q))])
  q_true <- qcauchy(qs)
  err_pct <- (q_est - q_true) / q_true * 100
  body_idx <- which(x_grid >= -10 & x_grid <= 10)
  body_mse <- mean((pdf_fitted[body_idx] - dcauchy(x_grid[body_idx]))^2)
  ks_idx <- which(x_grid >= -45 & x_grid <= 45)
  ks <- max(abs(cdf_fitted[ks_idx] - pcauchy(x_grid[ks_idx])))
  list(q_est = q_est, err_pct = err_pct, body_mse = body_mse, ks = ks,
       pdf = pdf_fitted)
}

# =========================================================================
# EXPERIMENT 1-PM: STANDARD CAUCHY, parity-matched PATP + T-MaxEnt
# =========================================================================
cat("=========================================================\n")
cat("EXP 1-PM: CAUCHY — PARITY-MATCHED PATP + T-MAXENT\n")
cat("=========================================================\n")
set.seed(20260612)
N <- 1000
data_cauchy <- rcauchy(N)
L <- 50
x_grid <- seq(-L, L, length.out = 2000)
dx <- x_grid[2] - x_grid[1]
n_moments <- 4

cat("\n-- PM-PATP sweep --\n")
cat("alpha | conv | kappa_H | potential | bodyMSE[-10,10] | KS[-45,45] | Q95est (err%) | Q99est (err%)\n")
res_c <- data.frame()
pdf_store <- list()
for (alpha in seq(0, 1, by = 0.1)) {
  Phi_data <- patp_basis_pm(data_cauchy, n_moments, alpha)
  tm <- colMeans(Phi_data)
  Phi_grid <- patp_basis_pm(x_grid, n_moments, alpha)
  fit <- solve_maxent(Phi_grid, tm, dx)
  if (is.null(fit)) {
    cat(sprintf("%.1f | solver-NULL (targets: %s)\n", alpha, paste(sprintf("%.4g", tm), collapse = ", ")))
    next
  }
  m <- fit_metrics_cauchy(Phi_grid, fit, x_grid, dx)
  pot <- log(fit$Z) - sum(fit$lambdas * tm)
  res_c <- rbind(res_c, data.frame(alpha = alpha, conv = fit$converged, kH = fit$cond_num,
                                   pot = pot, bmse = m$body_mse, ks = m$ks,
                                   q95 = m$q_est[3], e95 = m$err_pct[3],
                                   q99 = m$q_est[4], e99 = m$err_pct[4]))
  pdf_store[[sprintf("a%.1f", alpha)]] <- m$pdf
  cat(sprintf("%.1f | %s | %.3e | %.4f | %.4e | %.4f | %.2f (%+.1f%%) | %.2f (%+.1f%%)\n",
              alpha, fit$converged, fit$cond_num, pot, m$body_mse, m$ks,
              m$q_est[3], m$err_pct[3], m$q_est[4], m$err_pct[4]))
}
if (nrow(res_c) > 0) {
  opt <- res_c[which.min(res_c$pot), ]
  cat(sprintf("\nalpha* (min potential among returned fits) = %.1f | bodyMSE=%.4e | KS=%.4f\n",
              opt$alpha, opt$bmse, opt$ks))
  best_acc <- res_c[which.min(res_c$bmse), ]
  cat(sprintf("best body-MSE alpha = %.1f (%.4e)\n", best_acc$alpha, best_acc$bmse))
}

cat("\n-- T-MaxEnt on Cauchy (ECF constraints; CF of Cauchy = exp(-|u|)) --\n")
cat("p | S | conv | kappa_H | bodyMSE | KS | |ECF(S*p)| | Q95est(err%) | Q99est(err%)\n")
res_tc <- data.frame()
for (p_val in c(0.2, 0.3, 0.5, 1.0)) {
  for (S in c(2, 3, 4)) {
    Phi_data <- trig_basis(data_cauchy, S, p_val)
    tm <- colMeans(Phi_data)
    ecf_top <- sqrt(mean(cos(S * p_val * data_cauchy))^2 + mean(sin(S * p_val * data_cauchy))^2)
    Phi_grid <- trig_basis(x_grid, S, p_val)
    fit <- solve_maxent(Phi_grid, tm, dx)
    if (is.null(fit)) { cat(sprintf("%.1f | %d | solver-NULL\n", p_val, S)); next }
    m <- fit_metrics_cauchy(Phi_grid, fit, x_grid, dx)
    res_tc <- rbind(res_tc, data.frame(p = p_val, S = S, conv = fit$converged, kH = fit$cond_num,
                                       bmse = m$body_mse, ks = m$ks, ecf = ecf_top,
                                       e95 = m$err_pct[3], e99 = m$err_pct[4]))
    pdf_store[[sprintf("trig_p%.1f_S%d", p_val, S)]] <- m$pdf
    cat(sprintf("%.1f | %d | %s | %.3e | %.4e | %.4f | %.4f | %.2f (%+.1f%%) | %.2f (%+.1f%%)\n",
                p_val, S, fit$converged, fit$cond_num, m$body_mse, m$ks, ecf_top,
                m$q_est[3], m$err_pct[3], m$q_est[4], m$err_pct[4]))
  }
}

# =========================================================================
# EXPERIMENT 2-PM: GAUSSIAN MIXTURE, parity-matched PATP sweep
# =========================================================================
cat("\n=========================================================\n")
cat("EXP 2-PM: MIXTURE — PARITY-MATCHED PATP SWEEP\n")
cat("=========================================================\n")
set.seed(20260612)
N <- 1000
sd_true <- sqrt(0.32)
u <- runif(N)
data_mix <- ifelse(u < 0.4, rnorm(N, -1, sd_true), rnorm(N, 1, sd_true))
Lm <- 5
xg <- seq(-Lm, Lm, length.out = 1000)
dxm <- xg[2] - xg[1]
pdf_true_mix <- 0.4 * dnorm(xg, -1, sd_true) + 0.6 * dnorm(xg, 1, sd_true)
nm <- 6

cat("alpha | conv | kappa_H | potential | PDF MSE\n")
res_m <- data.frame()
mix_pdf_store <- list()
for (alpha in seq(0, 1, by = 0.1)) {
  Phi_data <- patp_basis_pm(data_mix, nm, alpha)
  tm <- colMeans(Phi_data)
  Phi_grid <- patp_basis_pm(xg, nm, alpha)
  fit <- solve_maxent(Phi_grid, tm, dxm)
  if (is.null(fit)) { cat(sprintf("%.1f | solver-NULL\n", alpha)); next }
  pdf_f <- exp(as.numeric(Phi_grid %*% fit$lambdas)) / fit$Z
  mse <- mean((pdf_f - pdf_true_mix)^2)
  pot <- log(fit$Z) - sum(fit$lambdas * tm)
  res_m <- rbind(res_m, data.frame(alpha = alpha, conv = fit$converged, kH = fit$cond_num,
                                   pot = pot, mse = mse))
  mix_pdf_store[[sprintf("a%.1f", alpha)]] <- pdf_f
  cat(sprintf("%.1f | %s | %.3e | %.4f | %.6e\n", alpha, fit$converged, fit$cond_num, pot, mse))
}
if (nrow(res_m) > 0) {
  opt_m <- res_m[which.min(res_m$pot), ]
  cat(sprintf("\nalpha* (min potential) = %.1f | MSE=%.6e | kappa_H=%.3e\n", opt_m$alpha, opt_m$mse, opt_m$kH))
  best_m <- res_m[which.min(res_m$mse), ]
  cat(sprintf("best-MSE alpha = %.1f (MSE=%.6e, kappa_H=%.3e)\n", best_m$alpha, best_m$mse, best_m$kH))
}

# T-MaxEnt reference on the mixture (p=0.5, S=3) — for figure regeneration
Phi_data_t <- trig_basis(data_mix, 3, 0.5)
fit_t <- solve_maxent(trig_basis(xg, 3, 0.5), colMeans(Phi_data_t), dxm)
pdf_t_mix <- exp(as.numeric(trig_basis(xg, 3, 0.5) %*% fit_t$lambdas)) / fit_t$Z
mse_t <- mean((pdf_t_mix - pdf_true_mix)^2)
cat(sprintf("\nT-MaxEnt (p=0.5,S=3) reference: kappa_H=%.3e, MSE=%.6e\n", fit_t$cond_num, mse_t))

# Monomial baseline (= PM-PATP alpha=1; cross-check vs ablation_mixture.R)
Phi_mono_d <- sapply(1:6, function(i) data_mix^i)
fit_mono <- solve_maxent(sapply(1:6, function(i) xg^i), colMeans(Phi_mono_d), dxm)
if (!is.null(fit_mono)) {
  pdf_mono <- exp(as.numeric(sapply(1:6, function(i) xg^i) %*% fit_mono$lambdas)) / fit_mono$Z
  cat(sprintf("Monomial x^i cross-check: kappa_H=%.3e, MSE=%.6e (expect 1.21e5 / 5.73e-4)\n",
              fit_mono$cond_num, mean((pdf_mono - pdf_true_mix)^2)))
}

# =========================================================================
# FIGURES (regenerated for the final paper variants)
# =========================================================================
OUT_DIR <- Sys.getenv("GENELEMENT_OUT_DIR", unset = "outputs")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Fig 1: Cauchy log-density, true vs PM-PATP (alpha*, plus alpha=0) vs T-MaxEnt best
if (nrow(res_c) > 0 && nrow(res_tc) > 0) {
  a_star <- res_c[which.min(res_c$pot), "alpha"]
  tc_best <- res_tc[res_tc$p == 0.2 & res_tc$S == 2, ]   # equal-budget (4 constraints)
  key_a <- sprintf("a%.1f", a_star)
  key_t <- sprintf("trig_p%.1f_S%d", tc_best$p, tc_best$S)
  pdf(file.path(OUT_DIR, "fig_cauchy_logpdf.pdf"), width = 7, height = 4.5)
  par(mar = c(4.5, 4.5, 1, 1))
  sel <- which(x_grid >= -20 & x_grid <= 20)
  plot(x_grid[sel], log10(dcauchy(x_grid[sel])), type = "l", lwd = 2.5, col = "black",
       xlab = "x", ylab = expression(log[10] ~ f(x)), ylim = c(-6, 0.2))
  lines(x_grid[sel], log10(pmax(pdf_store[["a0.0"]][sel], 1e-12)), lwd = 2, lty = 2, col = "#D55E00")
  lines(x_grid[sel], log10(pmax(pdf_store[[key_a]][sel], 1e-12)), lwd = 2, lty = 4, col = "#0072B2")
  lines(x_grid[sel], log10(pmax(pdf_store[[key_t]][sel], 1e-12)), lwd = 2, lty = 3, col = "#009E73")
  legend("topright", bty = "n", lwd = c(2.5, 2, 2, 2), lty = c(1, 2, 4, 3),
         col = c("black", "#D55E00", "#0072B2", "#009E73"),
         legend = c("True Cauchy",
                    expression(paste("PM-PATP, ", alpha, " = 0")),
                    bquote(paste("PM-PATP, ", alpha, "* = ", .(a_star))),
                    sprintf("T-MaxEnt (p=%.1f, S=%d)", tc_best$p, tc_best$S)))
  dev.off()
  cat(sprintf("\nFigure 1 written (PM-PATP alpha*=%.1f, T-MaxEnt p=%.1f S=%d)\n", a_star, tc_best$p, tc_best$S))
}

# Fig 2: Mixture density, true vs T-MaxEnt vs PM-PATP best vs monomial
if (nrow(res_m) > 0 && !is.null(fit_mono)) {
  b_m <- res_m[which.min(res_m$pot), "alpha"]  # show the potential-SELECTED member
  key_b <- sprintf("a%.1f", b_m)
  pdf(file.path(OUT_DIR, "fig_mixture_fits.pdf"), width = 7, height = 4.5)
  par(mar = c(4.5, 4.5, 1, 1))
  sel <- which(xg >= -3.5 & xg <= 3.5)
  plot(xg[sel], pdf_true_mix[sel], type = "l", lwd = 2.5, col = "black",
       xlab = "x", ylab = "f(x)", ylim = c(0, max(pdf_true_mix) * 1.15))
  lines(xg[sel], pdf_t_mix[sel], lwd = 2, lty = 2, col = "#009E73")
  lines(xg[sel], mix_pdf_store[[key_b]][sel], lwd = 2, lty = 4, col = "#0072B2")
  lines(xg[sel], pdf_mono[sel], lwd = 2, lty = 3, col = "#D55E00")
  legend("topleft", bty = "n", lwd = c(2.5, 2, 2, 2), lty = c(1, 2, 4, 3),
         col = c("black", "#009E73", "#0072B2", "#D55E00"),
         legend = c("True mixture", "T-MaxEnt (p=0.5, S=3)",
                    bquote(paste("PM-PATP, ", alpha, " = ", .(b_m))),
                    "Monomial (6 moments)"))
  dev.off()
  cat(sprintf("Figure 2 written (PM-PATP best alpha=%.1f)\n", b_m))
}

cat("\nDone.\n")

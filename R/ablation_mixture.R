#!/usr/bin/env Rscript

# Ablation & Baseline Study for Experiment 2 (Bimodal Gaussian Mixture)
# ---------------------------------------------------------------------
# Reproduces and extends Experiment 2 from patp_maxent_simulation.R /
# trig_maxent_simulation.R:
#   Mixture: 0.4*N(-1, sd=sqrt(0.32)) + 0.6*N(1, sd=sqrt(0.32))
#   seed 20260612, N = 1000 samples, grid [-5,5] with 1000 points,
#   6 constraints for all baselines.
#
# Sub-experiments:
#   (A) Baselines at 6 constraints: T-MaxEnt(p=0.5,S=3), PATP(0.9),
#       PATP(1.0), classical monomial x^i (i=1..6), Legendre P_1..P_6(x/5).
#   (B) T-MaxEnt (p,S) ablation: p in {0.3,0.4,0.5,0.7,1.0} x S in {2,3,4,5},
#       with |ECF(S*p)| to validate the design rule |f(Sp)| >= 0.05.
#   (C) Ridge sensitivity {1e-6, 1e-8, 1e-10} for T-MaxEnt(0.5,3),
#       PATP(0.9), monomial.
#
# Convention: kappa(H) is computed on the Hessian WITH the ridge added to
# its diagonal (exactly as in the existing scripts).
#
# Base R only. Fully reproducible (single seed, deterministic solver).

# -------------------------------------------------------------------------
# 1. CORE MAXENT SOLVER (copied verbatim from trig_maxent_simulation.R;
#    sole change: hardcoded 1e-8 ridge replaced by the `ridge` argument)
# -------------------------------------------------------------------------

solve_maxent <- function(Phi, target_moments, dx, max_iter = 100, tol = 1e-6,
                         ridge = 1e-8) {
  # Newton-Raphson solver with backtracking line search for dual MaxEnt
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

    # Check convergence
    if (sqrt(sum(grad^2)) < tol) {
      H <- t(Phi) %*% (Phi * pdf) * dx - fitted_moments %*% t(fitted_moments)
      diag(H) <- diag(H) + ridge # Regularization for stability
      cond_num <- kappa(H, exact = TRUE)
      return(list(lambdas = lambdas, Z = Z, grad = grad, H = H, cond_num = cond_num, converged = TRUE))
    }

    H <- t(Phi) %*% (Phi * pdf) * dx - fitted_moments %*% t(fitted_moments)
    diag(H) <- diag(H) + ridge

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
      lambdas <- lambdas - 0.1 * step
    }
  }

  unnorm_pdf <- exp(as.numeric(Phi %*% lambdas))
  Z <- sum(unnorm_pdf) * dx
  pdf <- unnorm_pdf / Z
  fitted_moments <- as.numeric(t(Phi) %*% pdf * dx)
  grad <- fitted_moments - target_moments
  H <- t(Phi) %*% (Phi * pdf) * dx - fitted_moments %*% t(fitted_moments)
  diag(H) <- diag(H) + ridge
  cond_num <- kappa(H, exact = TRUE)
  return(list(lambdas = lambdas, Z = Z, grad = grad, H = H, cond_num = cond_num, converged = FALSE))
}

# -------------------------------------------------------------------------
# 2. BASES
# -------------------------------------------------------------------------

# Trigonometric basis (verbatim from trig_maxent_simulation.R)
trig_basis <- function(x, S, p) {
  # Returns a matrix of size length(x) x 2S
  # Odd columns: cos(r * p * x)
  # Even columns: sin(r * p * x)
  M <- matrix(0, nrow = length(x), ncol = 2 * S)
  for (r in 1:S) {
    M[, 2 * r - 1] <- cos(r * p * x)
    M[, 2 * r]     <- sin(r * p * x)
  }
  return(M)
}

# PATP basis (verbatim from patp_maxent_simulation.R)
patp_power <- function(i, alpha) {
  A <- 1.0 / i
  B <- 4.0 - i - 3.0 / i
  C <- 2.0 * i - 4.0 + 2.0 / i
  return(A + B * alpha + C * alpha^2)
}

patp_basis <- function(x, n, alpha) {
  M <- matrix(0, nrow = length(x), ncol = n)
  M[, 1] <- x
  if (n > 1) {
    for (i in 2:n) {
      p <- patp_power(i, alpha)
      M[, i] <- sign(x) * (abs(x)^p)
    }
  }
  return(M)
}

# Classical monomial basis: phi_i(x) = x^i, i = 1..n (plain powers)
monomial_basis <- function(x, n) {
  M <- matrix(0, nrow = length(x), ncol = n)
  for (i in 1:n) {
    M[, i] <- x^i
  }
  return(M)
}

# Legendre orthogonal polynomial basis P_1..P_6 evaluated at t = x/L
# (rescale to [-1,1]); explicit Legendre polynomial formulas.
legendre_basis <- function(x, L) {
  t <- x / L
  M <- matrix(0, nrow = length(x), ncol = 6)
  M[, 1] <- t
  M[, 2] <- (3 * t^2 - 1) / 2
  M[, 3] <- (5 * t^3 - 3 * t) / 2
  M[, 4] <- (35 * t^4 - 30 * t^2 + 3) / 8
  M[, 5] <- (63 * t^5 - 70 * t^3 + 15 * t) / 8
  M[, 6] <- (231 * t^6 - 315 * t^4 + 105 * t^2 - 5) / 16
  return(M)
}

# Empirical characteristic function modulus |f_hat(u)| from a sample
ecf_modulus <- function(u, x) {
  sqrt(mean(cos(u * x))^2 + mean(sin(u * x))^2)
}

# -------------------------------------------------------------------------
# 3. DATA AND GRID (identical to Experiment 2)
# -------------------------------------------------------------------------

set.seed(20260612)
N <- 1000
sd_true <- sqrt(0.32)
u <- runif(N)
data_mix <- ifelse(u < 0.4, rnorm(N, -1, sd_true), rnorm(N, 1, sd_true))

L <- 5
x_grid <- seq(-L, L, length.out = 1000)
dx <- x_grid[2] - x_grid[1]
pdf_true <- 0.4 * dnorm(x_grid, -1, sd_true) + 0.6 * dnorm(x_grid, 1, sd_true)

# Generic fit-and-evaluate helper
fit_basis <- function(Phi_data, Phi_grid, ridge = 1e-8) {
  target_moments <- colMeans(Phi_data)
  fit <- solve_maxent(Phi_grid, target_moments, dx, ridge = ridge)
  if (is.null(fit)) {
    return(list(ok = FALSE))
  }
  pdf_fit <- exp(as.numeric(Phi_grid %*% fit$lambdas)) / fit$Z
  mse <- mean((pdf_fit - pdf_true)^2)
  potential <- log(fit$Z) - sum(fit$lambdas * target_moments)
  return(list(ok = TRUE, converged = fit$converged, cond_num = fit$cond_num,
              mse = mse, potential = potential, pdf = pdf_fit,
              grad_norm = sqrt(sum(fit$grad^2))))
}

cat("=========================================================\n")
cat("ABLATION & BASELINE STUDY: BIMODAL GAUSSIAN MIXTURE\n")
cat("0.4*N(-1, sd=sqrt(0.32)) + 0.6*N(1, sd=sqrt(0.32))\n")
cat(sprintf("seed = 20260612, N = %d, grid [-%d, %d] with %d points\n",
            N, L, L, length(x_grid)))
cat("kappa(H) convention: computed on Hessian WITH ridge added to diagonal\n")
cat("=========================================================\n")

# -------------------------------------------------------------------------
# 4. SUB-EXPERIMENT (A): BASELINES AT 6 CONSTRAINTS
# -------------------------------------------------------------------------

cat("\n---------------------------------------------------------\n")
cat("(A) BASELINES AT 6 CONSTRAINTS (ridge = 1e-8)\n")
cat("---------------------------------------------------------\n")

baselines <- list(
  list(name = "T-MaxEnt (p=0.5, S=3)",
       data = trig_basis(data_mix, 3, 0.5), grid = trig_basis(x_grid, 3, 0.5)),
  list(name = "PATP (alpha=0.9)",
       data = patp_basis(data_mix, 6, 0.9), grid = patp_basis(x_grid, 6, 0.9)),
  list(name = "PATP (alpha=1.0)",
       data = patp_basis(data_mix, 6, 1.0), grid = patp_basis(x_grid, 6, 1.0)),
  list(name = "Monomial x^i, i=1..6",
       data = monomial_basis(data_mix, 6), grid = monomial_basis(x_grid, 6)),
  list(name = "Legendre P_1..P_6(x/5)",
       data = legendre_basis(data_mix, L), grid = legendre_basis(x_grid, L))
)

res_A <- list()
cat(sprintf("%-24s | %-9s | %-12s | %-12s | %-12s\n",
            "Basis", "Converged", "kappa(H)", "PDF MSE", "Dual potential"))
cat(strrep("-", 84), "\n", sep = "")
for (b in baselines) {
  r <- fit_basis(b$data, b$grid)
  res_A[[b$name]] <- r
  if (!r$ok) {
    cat(sprintf("%-24s | %-9s | %-12s | %-12s | %-12s\n",
                b$name, "FAILED", "-", "-", "-"))
  } else {
    cat(sprintf("%-24s | %-9s | %.6e | %.6e | %.6f\n",
                b$name, r$converged, r$cond_num, r$mse, r$potential))
  }
}

# -------------------------------------------------------------------------
# 5. SUB-EXPERIMENT (B): T-MAXENT (p, S) ABLATION
# -------------------------------------------------------------------------

cat("\n---------------------------------------------------------\n")
cat("(B) T-MAXENT (p, S) ABLATION (ridge = 1e-8)\n")
cat("    |ECF(S*p)| = empirical CF modulus at highest harmonic frequency\n")
cat("    Design rule: |f(S*p)| >= 0.05\n")
cat("---------------------------------------------------------\n")

p_vals <- c(0.3, 0.4, 0.5, 0.7, 1.0)
S_vals <- c(2, 3, 4, 5)

cat(sprintf("%-5s | %-3s | %-6s | %-9s | %-12s | %-12s | %-10s\n",
            "p", "S", "u=S*p", "Converged", "kappa(H)", "PDF MSE", "|ECF(S*p)|"))
cat(strrep("-", 76), "\n", sep = "")
res_B <- data.frame()
for (p in p_vals) {
  for (S in S_vals) {
    u_max <- S * p
    cf_mod <- ecf_modulus(u_max, data_mix)
    r <- fit_basis(trig_basis(data_mix, S, p), trig_basis(x_grid, S, p))
    if (!r$ok) {
      cat(sprintf("%-5.1f | %-3d | %-6.2f | %-9s | %-12s | %-12s | %.6f\n",
                  p, S, u_max, "FAILED", "-", "-", cf_mod))
      res_B <- rbind(res_B, data.frame(p = p, S = S, u_max = u_max,
                                       converged = NA, cond_num = NA,
                                       mse = NA, cf_mod = cf_mod))
    } else {
      cat(sprintf("%-5.1f | %-3d | %-6.2f | %-9s | %.6e | %.6e | %.6f\n",
                  p, S, u_max, r$converged, r$cond_num, r$mse, cf_mod))
      res_B <- rbind(res_B, data.frame(p = p, S = S, u_max = u_max,
                                       converged = r$converged,
                                       cond_num = r$cond_num,
                                       mse = r$mse, cf_mod = cf_mod))
    }
  }
}

# -------------------------------------------------------------------------
# 6. SUB-EXPERIMENT (C): RIDGE SENSITIVITY
# -------------------------------------------------------------------------

cat("\n---------------------------------------------------------\n")
cat("(C) RIDGE SENSITIVITY: ridge in {1e-6, 1e-8, 1e-10}\n")
cat("    kappa(H) computed on Hessian WITH the ridge added\n")
cat("---------------------------------------------------------\n")

ridge_vals <- c(1e-6, 1e-8, 1e-10)
ridge_cases <- list(
  list(name = "T-MaxEnt (p=0.5, S=3)",
       data = trig_basis(data_mix, 3, 0.5), grid = trig_basis(x_grid, 3, 0.5)),
  list(name = "PATP (alpha=0.9)",
       data = patp_basis(data_mix, 6, 0.9), grid = patp_basis(x_grid, 6, 0.9)),
  list(name = "Monomial x^i, i=1..6",
       data = monomial_basis(data_mix, 6), grid = monomial_basis(x_grid, 6))
)

cat(sprintf("%-24s | %-7s | %-9s | %-12s | %-12s\n",
            "Basis", "Ridge", "Converged", "kappa(H)", "PDF MSE"))
cat(strrep("-", 74), "\n", sep = "")
for (cse in ridge_cases) {
  for (rv in ridge_vals) {
    r <- fit_basis(cse$data, cse$grid, ridge = rv)
    if (!r$ok) {
      cat(sprintf("%-24s | %-7.0e | %-9s | %-12s | %-12s\n",
                  cse$name, rv, "FAILED", "-", "-"))
    } else {
      cat(sprintf("%-24s | %-7.0e | %-9s | %.6e | %.6e\n",
                  cse$name, rv, r$converged, r$cond_num, r$mse))
    }
  }
}

# -------------------------------------------------------------------------
# 7. FIGURE 2: DENSITY FITS ON [-3.5, 3.5]
# -------------------------------------------------------------------------

# Best PATP = lower PDF MSE among alpha = 0.9 and alpha = 1.0
r_trig <- res_A[["T-MaxEnt (p=0.5, S=3)"]]
r_p09  <- res_A[["PATP (alpha=0.9)"]]
r_p10  <- res_A[["PATP (alpha=1.0)"]]
r_mono <- res_A[["Monomial x^i, i=1..6"]]

if (r_p09$ok && r_p10$ok) {
  if (r_p09$mse <= r_p10$mse) {
    r_patp <- r_p09; patp_label <- "PATP (alpha=0.9)"
  } else {
    r_patp <- r_p10; patp_label <- "PATP (alpha=1.0)"
  }
} else if (r_p09$ok) {
  r_patp <- r_p09; patp_label <- "PATP (alpha=0.9)"
} else {
  r_patp <- r_p10; patp_label <- "PATP (alpha=1.0)"
}
cat(sprintf("\nFigure: best PATP selected = %s (PDF MSE = %.6e)\n",
            patp_label, r_patp$mse))

fig_dir <- Sys.getenv("GENELEMENT_OUT_DIR", unset = "outputs")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
fig_path <- file.path(fig_dir, "fig_mixture_fits.pdf")

idx <- x_grid >= -3.5 & x_grid <= 3.5
fmt_e <- function(v) sprintf("%.2e", v)

pdf(fig_path, width = 7, height = 4.5)
par(mar = c(4.5, 4.5, 1, 1))
ymax <- max(pdf_true[idx], r_trig$pdf[idx], r_patp$pdf[idx], r_mono$pdf[idx])
plot(x_grid[idx], pdf_true[idx], type = "l", lwd = 2.5, col = "black", lty = 1,
     xlab = "x (dimensionless)", ylab = "Probability density f(x)",
     ylim = c(0, 1.05 * ymax), xaxs = "i")
lines(x_grid[idx], r_trig$pdf[idx], lwd = 2, lty = 2, col = "#0072B2")
lines(x_grid[idx], r_patp$pdf[idx], lwd = 2, lty = 4, col = "#D55E00")
lines(x_grid[idx], r_mono$pdf[idx], lwd = 2, lty = 3, col = "#009E73")

leg <- vector("expression", 4)
leg[[1]] <- bquote("True mixture")
leg[[2]] <- bquote("T-MaxEnt (p=0.5, S=3):" ~ kappa[H] == .(fmt_e(r_trig$cond_num)) * "," ~ "MSE" == .(fmt_e(r_trig$mse)))
leg[[3]] <- bquote(.(patp_label) * ":" ~ kappa[H] == .(fmt_e(r_patp$cond_num)) * "," ~ "MSE" == .(fmt_e(r_patp$mse)))
leg[[4]] <- bquote("Monomial" ~ x^i * ", i=1..6:" ~ kappa[H] == .(fmt_e(r_mono$cond_num)) * "," ~ "MSE" == .(fmt_e(r_mono$mse)))
legend("topleft", legend = leg,
       col = c("black", "#0072B2", "#D55E00", "#009E73"),
       lty = c(1, 2, 4, 3), lwd = c(2.5, 2, 2, 2),
       bty = "n", cex = 0.72, seg.len = 3)
dev.off()

cat(sprintf("Figure written: %s\n", fig_path))
cat("\nAblation study complete.\n")

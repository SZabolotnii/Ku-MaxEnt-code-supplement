#!/usr/bin/env Rscript

# PATP-MaxEnt Simulation and Verification Script
# Implements the Parametrically-Adaptive Transition Polynomial (PATP) as MaxEnt constraints.
# Verifies performance on: (1) Standard Cauchy (heavy tails), (2) Gaussian Mixture (multimodality).

# -------------------------------------------------------------------------
# 1. CORE PATP-MAXENT IMPLEMENTATION
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
# 2. EXPERIMENT 1: HEAVY TAILS (STANDARD CAUCHY)
# -------------------------------------------------------------------------
run_cauchy_experiment <- function() {
  cat("\n=========================================================\n")
  cat("EXPERIMENT 1: HEAVY-TAILED LAW (STANDARD CAUCHY)\n")
  cat("=========================================================\n")
  
  set.seed(20260612)
  N <- 1000
  data_cauchy <- rcauchy(N, location = 0, scale = 1)
  
  # Standard Cauchy quantiles for comparison
  q_true_95 <- qcauchy(0.95, 0, 1) # ~6.314
  q_true_99 <- qcauchy(0.99, 0, 1) # ~31.82
  cat(sprintf("True 95th Percentile: %.4f\n", q_true_95))
  cat(sprintf("True 99th Percentile: %.4f\n", q_true_99))
  
  # We construct a grid over [-50, 50] with 2000 points to capture the Cauchy range
  L <- 50
  x_grid <- seq(-L, L, length.out = 2000)
  dx <- x_grid[2] - x_grid[1]
  n_moments <- 4 # We use 4 moments (S=4)
  
  # Sweep alpha ∈ [0, 1]
  alpha_grid <- seq(0.0, 1.0, by = 0.1)
  results <- data.frame()
  
  for (alpha in alpha_grid) {
    # Compute empirical PATP moments
    Phi_data <- patp_basis(data_cauchy, n_moments, alpha)
    target_moments <- colMeans(Phi_data)
    
    # Fit MaxEnt on grid
    Phi_grid <- patp_basis(x_grid, n_moments, alpha)
    fit <- solve_maxent(Phi_grid, target_moments, dx)
    
    if (is.null(fit)) {
      cat(sprintf("Alpha = %.1f: Solver failed to converge (divergent moments)\n", alpha))
      next
    }
    
    # Calculate fitted PDF
    pdf_fitted <- exp(as.numeric(Phi_grid %*% fit$lambdas)) / fit$Z
    cdf_fitted <- cumsum(pdf_fitted) * dx
    
    # Estimate quantiles
    q_est_95 <- x_grid[which.min(abs(cdf_fitted - 0.95))]
    q_est_99 <- x_grid[which.min(abs(cdf_fitted - 0.99))]
    
    # Compute potential: log(Z) - sum(lambdas * target)
    potential <- log(fit$Z) - sum(fit$lambdas * target_moments)
    
    results <- rbind(results, data.frame(
      alpha = alpha,
      converged = fit$converged,
      cond_num = fit$cond_num,
      potential = potential,
      q95 = q_est_95,
      q99 = q_est_99,
      err_q95_pct = abs(q_est_95 - q_true_95) / q_true_95 * 100,
      err_q99_pct = abs(q_est_99 - q_true_99) / q_true_99 * 100
    ))
    
    cat(sprintf("Alpha = %.1f: Conv=%s, Cond(H)=%.2e, Potential=%.4f, Q95=%.3f (Err=%.1f%%), Q99=%.3f (Err=%.1f%%)\n",
                alpha, fit$converged, fit$cond_num, potential, 
                q_est_95, abs(q_est_95 - q_true_95) / q_true_95 * 100,
                q_est_99, abs(q_est_99 - q_true_99) / q_true_99 * 100))
  }
  
  # Find optimal alpha* (minimizes dual potential)
  opt_idx <- which.min(results$potential)
  opt_res <- results[opt_idx, ]
  cat(sprintf("\n--> Optimal Alpha* (Min Dual Potential): %.1f\n", opt_res$alpha))
  cat(sprintf("--> Optimal Q95 Estimate: %.3f (Error = %.1f%%)\n", opt_res$q95, opt_res$err_q95_pct))
  cat(sprintf("--> Optimal Q99 Estimate: %.3f (Error = %.1f%%)\n", opt_res$q99, opt_res$err_q99_pct))
  
  # Compare alpha = 0 (pure fractional) vs alpha = 1 (pure integer)
  cat("\nComparison of Regimes:\n")
  cat("Regime       | Cond(H)   | Potential | Q95 Error | Q99 Error\n")
  cat("------------------------------------------------------------\n")
  for (r_name in c("Fractal (a=0)", "Linear (a=0.5)", "Integer (a=1)")) {
    a_val <- if (r_name == "Fractal (a=0)") 0.0 else if (r_name == "Linear (a=0.5)") 0.5 else 1.0
    row <- results[abs(results$alpha - a_val) < 1e-5, ]
    if (nrow(row) > 0) {
      cat(sprintf("%-12s | %.2e | %9.4f | %9.1f%% | %9.1f%%\n",
                  r_name, row$cond_num, row$potential, row$err_q95_pct, row$err_q99_pct))
    } else {
      cat(sprintf("%-12s | FAILED    |      -    |      -    |      -\n", r_name))
    }
  }
}

# -------------------------------------------------------------------------
# 3. EXPERIMENT 2: MULTIMODAL DISTRIBUTION (GAUSSIAN MIXTURE)
# -------------------------------------------------------------------------
run_mixture_experiment <- function() {
  cat("\n=========================================================\n")
  cat("EXPERIMENT 2: MULTIMODAL DISTRIBUTION (GAUSSIAN MIXTURE)\n")
  cat("=========================================================\n")
  
  set.seed(20260612)
  N <- 1000
  # True PDF: 0.4 * N(-1, sd=sqrt(0.32)) + 0.6 * N(1, sd=sqrt(0.32))
  sd_true <- sqrt(0.32)
  u <- runif(N)
  data_mix <- ifelse(u < 0.4, rnorm(N, -1, sd_true), rnorm(N, 1, sd_true))
  
  # Set grid over [-5, 5]
  L <- 5
  x_grid <- seq(-L, L, length.out = 1000)
  dx <- x_grid[2] - x_grid[1]
  n_moments <- 6 # We use 6 moments (S=6) to capture the multimodality
  
  alpha_grid <- seq(0.0, 1.0, by = 0.1)
  results <- data.frame()
  
  # True PDF values for MSE calculation
  pdf_true <- 0.4 * dnorm(x_grid, -1, sd_true) + 0.6 * dnorm(x_grid, 1, sd_true)
  
  for (alpha in alpha_grid) {
    Phi_data <- patp_basis(data_mix, n_moments, alpha)
    target_moments <- colMeans(Phi_data)
    
    Phi_grid <- patp_basis(x_grid, n_moments, alpha)
    fit <- solve_maxent(Phi_grid, target_moments, dx)
    
    if (is.null(fit)) {
      cat(sprintf("Alpha = %.1f: Solver failed to converge\n", alpha))
      next
    }
    
    # Calculate fitted PDF
    pdf_fitted <- exp(as.numeric(Phi_grid %*% fit$lambdas)) / fit$Z
    
    # Compute PDF Mean Squared Error (MSE)
    pdf_mse <- mean((pdf_fitted - pdf_true)^2)
    potential <- log(fit$Z) - sum(fit$lambdas * target_moments)
    
    results <- rbind(results, data.frame(
      alpha = alpha,
      converged = fit$converged,
      cond_num = fit$cond_num,
      potential = potential,
      pdf_mse = pdf_mse
    ))
    
    cat(sprintf("Alpha = %.1f: Conv=%s, Cond(H)=%.2e, Potential=%.4f, PDF MSE=%.6f\n",
                alpha, fit$converged, fit$cond_num, potential, pdf_mse))
  }
  
  opt_idx <- which.min(results$potential)
  opt_res <- results[opt_idx, ]
  cat(sprintf("\n--> Optimal Alpha* (Min Dual Potential): %.1f\n", opt_res$alpha))
  cat(sprintf("--> Optimal PDF MSE: %.6f\n", opt_res$pdf_mse))
  
  # Compare alpha = 0 (pure fractional) vs alpha = 1 (pure integer)
  cat("\nComparison of Regimes:\n")
  cat("Regime       | Cond(H)   | Potential | PDF MSE\n")
  cat("-------------------------------------------------\n")
  for (r_name in c("Fractal (a=0)", "Linear (a=0.5)", "Integer (a=1)")) {
    a_val <- if (r_name == "Fractal (a=0)") 0.0 else if (r_name == "Linear (a=0.5)") 0.5 else 1.0
    row <- results[abs(results$alpha - a_val) < 1e-5, ]
    if (nrow(row) > 0) {
      cat(sprintf("%-12s | %.2e | %9.4f | %.6f\n",
                  r_name, row$cond_num, row$potential, row$pdf_mse))
    } else {
      cat(sprintf("%-12s | FAILED    |      -    |      -\n", r_name))
    }
  }
}

# -------------------------------------------------------------------------
# MAIN EXECUTION
# -------------------------------------------------------------------------
main <- function() {
  cat("=========================================================\n")
  cat("RUNNING PATP-MAXENT SIMULATION AND VERIFICATION IN R\n")
  cat("=========================================================\n")
  
  run_cauchy_experiment()
  run_mixture_experiment()
  
  cat("\nVerification complete.\n")
}

main()

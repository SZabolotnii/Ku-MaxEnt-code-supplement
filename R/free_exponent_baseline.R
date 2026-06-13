#!/usr/bin/env Rscript

# FM-MEM-style free-exponent baseline on the mixture benchmark:
# optimize the three exponents (q2,q3,q4,...,q6 for 6 constraints -> 5 free
# exponents) of a parity-matched basis WITHOUT the one-parameter map,
# by Nelder-Mead over the converged dual potential (same selection criterion
# as the alpha-scan). Multi-start. Compares cost (solver calls) and outcome
# (MSE, kappa_H) against the 11-point one-dimensional alpha scan.
# Solver copied verbatim from patp_maxent_simulation.R.

solve_maxent <- function(Phi, target_moments, dx, max_iter = 100, tol = 1e-6) {
  n <- ncol(Phi); lambdas <- rep(0, n)
  for (iter in 1:max_iter) {
    unnorm_pdf <- exp(as.numeric(Phi %*% lambdas))
    if (any(is.infinite(unnorm_pdf)) || any(is.na(unnorm_pdf))) return(NULL)
    Z <- sum(unnorm_pdf) * dx
    if (Z == 0 || is.na(Z)) return(NULL)
    pdf <- unnorm_pdf / Z
    fitted <- as.numeric(t(Phi) %*% pdf * dx)
    grad <- fitted - target_moments
    if (sqrt(sum(grad^2)) < tol) {
      H <- t(Phi) %*% (Phi * pdf) * dx - fitted %*% t(fitted)
      diag(H) <- diag(H) + 1e-8
      return(list(lambdas = lambdas, Z = Z, cond_num = kappa(H, exact = TRUE), converged = TRUE))
    }
    H <- t(Phi) %*% (Phi * pdf) * dx - fitted %*% t(fitted)
    diag(H) <- diag(H) + 1e-8
    step <- tryCatch(solve(H, grad), error = function(e) tryCatch(qr.solve(H, grad), error = function(e2) NULL))
    if (is.null(step)) return(NULL)
    a_s <- 1.0; ok <- FALSE
    for (ls in 1:10) {
      nl <- lambdas - a_s * step
      nu <- exp(as.numeric(Phi %*% nl))
      if (!any(is.infinite(nu)) && !any(is.na(nu))) {
        nZ <- sum(nu) * dx
        if (nZ > 0 && !is.na(nZ)) {
          if (log(nZ) - sum(nl * target_moments) <= log(Z) - sum(lambdas * target_moments) + 1e-4) {
            lambdas <- nl; ok <- TRUE; break
          }
        }
      }
      a_s <- a_s * 0.5
    }
    if (!ok) lambdas <- lambdas - 0.1 * step
  }
  unnorm_pdf <- exp(as.numeric(Phi %*% lambdas))
  Z <- sum(unnorm_pdf) * dx
  pdf <- unnorm_pdf / Z
  fitted <- as.numeric(t(Phi) %*% pdf * dx)
  H <- t(Phi) %*% (Phi * pdf) * dx - fitted %*% t(fitted)
  diag(H) <- diag(H) + 1e-8
  list(lambdas = lambdas, Z = Z, cond_num = kappa(H, exact = TRUE), converged = FALSE)
}

basis_free <- function(x, q) {
  # parity-matched basis with free exponent vector q = (q2,...,q6); phi1 = x
  n <- length(q) + 1
  M <- matrix(0, nrow = length(x), ncol = n)
  M[, 1] <- x
  for (i in 2:n) {
    p <- q[i - 1]
    M[, i] <- if (i %% 2 == 0) abs(x)^p else sign(x) * (abs(x)^p)
  }
  M
}

set.seed(20260612)
N <- 1000; sd_true <- sqrt(0.32)
u <- runif(N)
data_mix <- ifelse(u < 0.4, rnorm(N, -1, sd_true), rnorm(N, 1, sd_true))
xg <- seq(-5, 5, length.out = 1000); dxm <- xg[2] - xg[1]
pdf_true <- 0.4 * dnorm(xg, -1, sd_true) + 0.6 * dnorm(xg, 1, sd_true)

n_calls <- 0
objective <- function(q) {
  n_calls <<- n_calls + 1
  if (any(q <= 0.05) || any(q > 6)) return(1e6)        # box constraint
  tm <- colMeans(basis_free(data_mix, q))
  fit <- solve_maxent(basis_free(xg, q), tm, dxm)
  if (is.null(fit) || !fit$converged) return(1e6)
  log(fit$Z) - sum(fit$lambdas * tm)                    # dual potential
}

evaluate <- function(q) {
  tm <- colMeans(basis_free(data_mix, q))
  Phi_g <- basis_free(xg, q)
  fit <- solve_maxent(Phi_g, tm, dxm)
  if (is.null(fit)) return(NULL)
  pdf_f <- exp(as.numeric(Phi_g %*% fit$lambdas)) / fit$Z
  list(mse = mean((pdf_f - pdf_true)^2), kH = fit$cond_num,
       pot = log(fit$Z) - sum(fit$lambdas * tm), conv = fit$converged)
}

patp_power <- function(i, alpha) 1/i + (4 - i - 3/i)*alpha + (2*i - 4 + 2/i)*alpha^2

starts <- list(
  fractional = sapply(2:6, function(i) 1/i),                 # PATP alpha=0 exponents
  integerish = c(2, 3, 4, 5, 6) * 0.5,                       # mid-scale
  patp_star  = sapply(2:6, function(i) patp_power(i, 0.7))   # warm start at the scan optimum
)

cat("==========================================================\n")
cat("FREE-EXPONENT (FM-MEM-STYLE) BASELINE, mixture, 6 constraints\n")
cat("Nelder-Mead over 5 free exponents; objective = converged dual potential\n")
cat("==========================================================\n\n")
results <- list()
for (nm in names(starts)) {
  n_calls <<- 0
  t0 <- Sys.time()
  opt <- optim(starts[[nm]], objective, method = "Nelder-Mead",
               control = list(maxit = 300, reltol = 1e-6))
  t1 <- Sys.time()
  ev <- evaluate(opt$par)
  results[[nm]] <- list(opt = opt, ev = ev, calls = n_calls,
                        secs = as.numeric(difftime(t1, t0, units = "secs")))
  cat(sprintf("start=%-10s: exponents (%s) -> potential %.4f, PDF MSE %.4e, kappa_H %.3e, conv=%s\n",
              nm, paste(sprintf("%.3f", opt$par), collapse = ", "),
              ev$pot, ev$mse, ev$kH, ev$conv))
  cat(sprintf("               solver calls: %d, wall time: %.1f s, NM convergence code: %d\n\n",
              n_calls, results[[nm]]$secs, opt$convergence))
}

best <- names(results)[which.min(sapply(results, function(r) r$ev$pot))]
cat(sprintf("Best start by potential: %s\n", best))
cat("\nReference (1-D PATP scan, same data): 11 solver calls total;\n")
cat("selected alpha*=0.7 -> potential 1.4198, MSE 8.028e-05, kappa_H 2.431e+04\n")
cat("\nDone.\n")

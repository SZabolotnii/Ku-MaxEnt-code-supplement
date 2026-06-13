#!/usr/bin/env Rscript

# PATP-MaxEnt Ablation: GA Fitness Noise — Monte Carlo vs Analytical MaxEnt
# Extends Experiment 3 (patp_maxent_mv_opt.R). All core functions below are
# copied VERBATIM from patp_maxent_mv_opt.R (do not edit them here).
#
# Sub-experiments:
#   (A) Fitness-noise probe at fixed candidate n = (100, 50, 30)
#   (B) Seeded GA comparison, 10 runs per method, matched seeds
#   (C) Decomposition: per-seed cost pairs; MaxEnt-arm sd = GA stochasticity alone
#
# Base R only. Fully seeded for reproducibility.

# -------------------------------------------------------------------------
# 1. CORE MOMENT CONVERSION AND PROPAGATION  (verbatim from patp_maxent_mv_opt.R)
# -------------------------------------------------------------------------

raw_moments_to_cumulants <- function(mu_prime) {
  r_max <- length(mu_prime)
  kappa <- rep(0, r_max)
  for (r in 1:r_max) {
    val <- mu_prime[r]
    if (r > 1) {
      for (j in 1:(r-1)) {
        val <- val - choose(r-1, j-1) * kappa[j] * mu_prime[r-j]
      }
    }
    kappa[r] <- val
  }
  return(kappa)
}

cumulants_to_raw_moments <- function(kappa) {
  r_max <- length(kappa)
  mu_prime <- rep(0, r_max)
  for (r in 1:r_max) {
    val <- 0
    for (j in 1:r) {
      coef_term <- choose(r-1, j-1)
      prev_mu <- if (r - j == 0) 1 else mu_prime[r - j]
      val <- val + coef_term * kappa[j] * prev_mu
    }
    mu_prime[r] <- val
  }
  return(mu_prime)
}

sample_mean_raw_moments <- function(pop_raw_moments, n) {
  pop_cumulants <- raw_moments_to_cumulants(pop_raw_moments)
  r_max <- length(pop_cumulants)
  sample_cumulants <- pop_cumulants
  for (r in 1:r_max) {
    sample_cumulants[r] <- pop_cumulants[r] / (n^(r-1))
  }
  sample_raw_moments <- cumulants_to_raw_moments(sample_cumulants)
  return(sample_raw_moments)
}

product_raw_moments <- function(n1, n2, n3, a = 1000, r_max = 4) {
  # 1. W1 = mean(X1) where X1 ~ Beta(2, 5)
  alpha_b <- 2
  beta_b <- 5
  pop_raw_W1 <- rep(0, r_max)
  val <- 1
  for (r in 1:r_max) {
    val <- val * (alpha_b + r - 1) / (alpha_b + beta_b + r - 1)
    pop_raw_W1[r] <- val
  }
  sample_raw_W1 <- sample_mean_raw_moments(pop_raw_W1, n1)

  # 2. W2 = mean(X2) where X2 ~ N(100, 10)
  # kappa_1 = 100, kappa_2 = 100/n2, higher cumulants = 0
  cumulants_W2 <- rep(0, r_max)
  cumulants_W2[1] <- 100
  cumulants_W2[2] <- 100 / n2
  sample_raw_W2 <- cumulants_to_raw_moments(cumulants_W2)

  # 3. W3 = mean(X3) - 1 where X3 ~ N(1.2, 0.05)
  # kappa_1 = 0.2, kappa_2 = 0.0025/n3, higher cumulants = 0
  cumulants_W3 <- rep(0, r_max)
  cumulants_W3[1] <- 0.2
  cumulants_W3[2] <- 0.0025 / n3
  sample_raw_W3 <- cumulants_to_raw_moments(cumulants_W3)

  # 4. Raw moments of V = W1 * W2 * W3
  raw_V <- sample_raw_W1 * sample_raw_W2 * sample_raw_W3

  # 5. Raw moments of Y = a * V
  raw_Y <- raw_V
  for (r in 1:r_max) {
    raw_Y[r] <- (a^r) * raw_V[r]
  }
  return(raw_Y)
}

# -------------------------------------------------------------------------
# 2. MAXENT QUANTILE SOLVER  (verbatim from patp_maxent_mv_opt.R)
# -------------------------------------------------------------------------

solve_maxent_raw <- function(raw_moments) {
  # Fits a 4-moment classical MaxEnt distribution using the exact raw moments.
  # We first centralize and normalize the moments to mean=0, variance=1
  m1 <- raw_moments[1]
  m2 <- raw_moments[2]
  m3 <- raw_moments[3]
  m4 <- raw_moments[4]

  variance <- m2 - m1^2
  sd_val <- sqrt(variance)

  # Normalized central moments
  m_tilde_1 <- 0
  m_tilde_2 <- 1
  m_tilde_3 <- (m3 - 3*m2*m1 + 2*m1^3) / (sd_val^3)
  m_tilde_4 <- (m4 - 4*m3*m1 + 6*m2*m1^2 - 3*m1^4) / (sd_val^4)

  target_moments <- c(m_tilde_1, m_tilde_2, m_tilde_3, m_tilde_4)

  # Grid for integration (normalized space)
  x_grid <- seq(-6, 6, length.out = 1000)
  dx <- x_grid[2] - x_grid[1]

  # Monomial basis
  Phi <- matrix(0, nrow = length(x_grid), ncol = 4)
  for (i in 1:4) {
    Phi[, i] <- x_grid^i
  }

  # Newton solver
  lambdas <- rep(0, 4)
  converged <- FALSE

  for (iter in 1:50) {
    unnorm <- exp(as.numeric(Phi %*% lambdas))
    if (any(is.infinite(unnorm)) || any(is.na(unnorm))) break
    Z <- sum(unnorm) * dx
    if (Z == 0 || is.na(Z)) break

    pdf <- unnorm / Z
    fitted_moments <- as.numeric(t(Phi) %*% pdf * dx)
    grad <- fitted_moments - target_moments

    if (sqrt(sum(grad^2)) < 1e-5) {
      converged <- TRUE
      break
    }

    H <- t(Phi) %*% (Phi * pdf) * dx - fitted_moments %*% t(fitted_moments)
    diag(H) <- diag(H) + 1e-7

    step <- tryCatch(solve(H, grad), error = function(e) {
      tryCatch(qr.solve(H, grad), error = function(e2) NULL)
    })
    if (is.null(step)) break

    # Simple line search
    alpha <- 1.0
    for (ls in 1:5) {
      new_lambdas <- lambdas - alpha * step
      new_unnorm <- exp(as.numeric(Phi %*% new_lambdas))
      if (!any(is.infinite(new_unnorm)) && !any(is.na(new_unnorm))) {
        lambdas <- new_lambdas
        break
      }
      alpha <- alpha * 0.5
    }
  }

  if (!converged) return(NULL)

  # Compute CDF and extract the 5% quantile in normalized space
  pdf_final <- exp(as.numeric(Phi %*% lambdas)) / Z
  cdf_final <- cumsum(pdf_final) * dx
  q05_norm <- x_grid[which.min(abs(cdf_final - 0.05))]

  # Denormalize to get the 5% quantile in original space
  q05_orig <- m1 + q05_norm * sd_val
  return(q05_orig)
}

# -------------------------------------------------------------------------
# 3. MONTE CARLO ESTIMATOR  (verbatim from patp_maxent_mv_opt.R)
# -------------------------------------------------------------------------

mc_quantile_05 <- function(n1, n2, n3, a = 1000, N_mc = 10000) {
  # Simulates the sample means and returns the 5% quantile of Y
  # 1. W1 ~ Beta(2, 5) sample mean of size n1
  W1_samples <- replicate(N_mc, mean(rbeta(n1, 2, 5)))
  # 2. W2 ~ N(100, 10) sample mean of size n2
  W2_samples <- rnorm(N_mc, 100, 10 / sqrt(n2))
  # 3. W3 ~ N(1.2, 0.05) sample mean - 1 of size n3
  W3_samples <- rnorm(N_mc, 0.2, 0.05 / sqrt(n3))

  Y_samples <- a * W1_samples * W2_samples * W3_samples
  return(quantile(Y_samples, 0.05, names = FALSE))
}

# -------------------------------------------------------------------------
# 4. OPTIMIZATION AND GENETIC ALGORITHM  (verbatim from patp_maxent_mv_opt.R)
# -------------------------------------------------------------------------

cost_function <- function(n) {
  # Cost: n1 + 2.5 * n2 + 5 * n3
  return(n[1] + 2.5 * n[2] + 5.0 * n[3])
}

evaluate_fitness <- function(n, method, N_mc = 10000) {
  n1 <- round(n[1])
  n2 <- round(n[2])
  n3 <- round(n[3])

  # Bounds check
  if (n1 < 10 || n1 > 500 || n2 < 10 || n2 > 500 || n3 < 10 || n3 > 500) {
    return(-1e6)
  }

  cost <- cost_function(c(n1, n2, n3))

  # Compute the 5% quantile of Y
  q05 <- if (method == "MaxEnt") {
    raw_moments <- product_raw_moments(n1, n2, n3, a = 1000, r_max = 4)
    q <- solve_maxent_raw(raw_moments)
    if (is.null(q)) return(-1e6)
    q
  } else {
    mc_quantile_05(n1, n2, n3, a = 1000, N_mc = N_mc)
  }

  # Constraint: q05 >= 4800
  if (q05 < 4800) {
    # Penalize constraint violation
    penalty <- 1000 * (4800 - q05)
    return(-(cost + penalty))
  } else {
    return(-cost)
  }
}

run_genetic_algorithm <- function(method, N_mc = 10000, pop_size = 30, generations = 20) {
  # Simple integer-based Genetic Algorithm
  # Population is a matrix of size pop_size x 3
  pop <- matrix(round(runif(pop_size * 3, 10, 150)), nrow = pop_size, ncol = 3)

  for (gen in 1:generations) {
    # 1. Evaluate fitness
    fitness <- rep(0, pop_size)
    for (i in 1:pop_size) {
      fitness[i] <- evaluate_fitness(pop[i, ], method, N_mc)
    }

    # 2. Selection (Tournament)
    new_pop <- pop
    for (i in 1:pop_size) {
      idx1 <- sample(1:pop_size, 1)
      idx2 <- sample(1:pop_size, 1)
      winner <- if (fitness[idx1] > fitness[idx2]) idx1 else idx2
      new_pop[i, ] <- pop[winner, ]
    }
    pop <- new_pop

    # 3. Crossover (Uniform)
    for (i in seq(1, pop_size - 1, by = 2)) {
      if (runif(1) < 0.8) {
        mask <- runif(3) < 0.5
        child1 <- ifelse(mask, pop[i, ], pop[i+1, ])
        child2 <- ifelse(mask, pop[i+1, ], pop[i, ])
        pop[i, ] <- child1
        pop[i+1, ] <- child2
      }
    }

    # 4. Mutation (Random perturb)
    for (i in 1:pop_size) {
      if (runif(1) < 0.2) {
        gene <- sample(1:3, 1)
        pop[i, gene] <- pop[i, gene] + round(rnorm(1, 0, 15))
        # Keep inside bounds
        pop[i, gene] <- max(10, min(500, pop[i, gene]))
      }
    }
  }

  # Final evaluation
  fitness <- rep(0, pop_size)
  for (i in 1:pop_size) {
    fitness[i] <- evaluate_fitness(pop[i, ], method, N_mc)
  }

  best_idx <- which.max(fitness)
  best_sol <- round(pop[best_idx, ])
  best_cost <- cost_function(best_sol)

  # True quantile evaluation (using large Monte Carlo N = 2*10^5 to check constraint violation)
  true_q05 <- mc_quantile_05(best_sol[1], best_sol[2], best_sol[3], a = 1000, N_mc = 200000)
  violation <- true_q05 < 4800

  return(list(sol = best_sol, cost = best_cost, true_q05 = true_q05, violation = violation))
}

# =========================================================================
# ABLATION MAIN (new code)
# =========================================================================

fmt_num <- function(x, d = 4) formatC(x, format = "f", digits = d)

main <- function() {
  cat("=========================================================\n")
  cat("ABLATION: GA FITNESS NOISE — MC vs ANALYTICAL MAXENT\n")
  cat("=========================================================\n")
  cat(sprintf("R version: %s\n", R.version.string))
  cat(sprintf("Date: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

  # -----------------------------------------------------------------------
  # (A) FITNESS-NOISE PROBE at fixed candidate n = (100, 50, 30)
  # -----------------------------------------------------------------------
  cat("---------------------------------------------------------\n")
  cat("(A) FITNESS-NOISE PROBE: candidate n = (100, 50, 30)\n")
  cat("---------------------------------------------------------\n")
  n_cand <- c(100, 50, 30)
  K <- 50
  cat(sprintf("Deterministic cost of candidate: %.1f\n", cost_function(n_cand)))

  # A.1 Monte Carlo fitness (N_mc = 5000), 50 evaluations
  set.seed(20260612)
  mc_vals <- numeric(K)
  t0 <- Sys.time()
  for (k in 1:K) {
    mc_vals[k] <- evaluate_fitness(n_cand, "MC", N_mc = 5000)
  }
  t1 <- Sys.time()
  mc_total_time <- as.numeric(difftime(t1, t0, units = "secs"))
  mc_percall <- mc_total_time / K

  cat("\nMC fitness (N_mc = 5000), 50 evaluations:\n")
  cat("All 50 values:\n")
  print(round(mc_vals, 4))
  cat(sprintf("  mean = %s\n", fmt_num(mean(mc_vals))))
  cat(sprintf("  sd   = %s\n", fmt_num(sd(mc_vals))))
  cat(sprintf("  min  = %s\n", fmt_num(min(mc_vals))))
  cat(sprintf("  max  = %s\n", fmt_num(max(mc_vals))))
  cat(sprintf("  range (max - min) = %s\n", fmt_num(max(mc_vals) - min(mc_vals))))
  cat(sprintf("  total wall time for 50 calls = %.4f s\n", mc_total_time))
  cat(sprintf("  mean per-call wall time      = %.6f s\n", mc_percall))

  # A.2 Analytical MaxEnt fitness, 50 evaluations (deterministic; no RNG used)
  set.seed(20260612)  # set identically for symmetry; MaxEnt path consumes no RNG
  me_vals <- numeric(K)
  t0 <- Sys.time()
  for (k in 1:K) {
    me_vals[k] <- evaluate_fitness(n_cand, "MaxEnt")
  }
  t1 <- Sys.time()
  me_total_time <- as.numeric(difftime(t1, t0, units = "secs"))
  me_percall <- me_total_time / K

  cat("\nMaxEnt analytical fitness, 50 evaluations:\n")
  cat(sprintf("  mean = %s\n", fmt_num(mean(me_vals))))
  cat(sprintf("  sd   = %s\n", fmt_num(sd(me_vals))))
  cat(sprintf("  unique values among 50 calls = %d\n", length(unique(me_vals))))
  cat(sprintf("  total wall time for 50 calls = %.4f s\n", me_total_time))
  cat(sprintf("  mean per-call wall time      = %.6f s\n", me_percall))

  # Assertion: analytical fitness must be exactly deterministic
  me_sd <- sd(me_vals)
  if (me_sd == 0) {
    cat("  ASSERTION PASSED: sd(MaxEnt fitness over 50 calls) == 0 exactly.\n")
  } else {
    cat(sprintf("  ASSERTION FAILED: sd(MaxEnt fitness) = %.10e != 0\n", me_sd))
  }
  stopifnot(me_sd == 0)

  cat(sprintf("\nPer-call speedup (MC time / MaxEnt time): %.2fx\n",
              mc_percall / me_percall))
  cat(sprintf("(MC per-call = %.6f s, MaxEnt per-call = %.6f s)\n",
              mc_percall, me_percall))

  # A.3 Latent noise diagnostic: the underlying MC q05 estimates at the same
  # candidate. With the same seed, these are EXACTLY the q05 values drawn
  # inside the 50 MC fitness calls of A.1 (evaluate_fitness consumes RNG only
  # via mc_quantile_05). The penalty in evaluate_fitness is one-sided, so the
  # fitness is constant whenever none of the q05 draws crosses the 4800
  # threshold — the MC noise is latent in q05 itself.
  set.seed(20260612)
  q05_vals <- numeric(K)
  for (k in 1:K) {
    q05_vals[k] <- mc_quantile_05(n_cand[1], n_cand[2], n_cand[3], a = 1000, N_mc = 5000)
  }
  cat("\nLatent MC q05 estimates underlying the 50 MC fitness calls:\n")
  cat(sprintf("  mean = %s\n", fmt_num(mean(q05_vals))))
  cat(sprintf("  sd   = %s\n", fmt_num(sd(q05_vals))))
  cat(sprintf("  min  = %s\n", fmt_num(min(q05_vals))))
  cat(sprintf("  max  = %s\n", fmt_num(max(q05_vals))))
  cat(sprintf("  draws below the 4800 constraint threshold: %d / %d\n",
              sum(q05_vals < 4800), K))
  me_q05_fixed <- solve_maxent_raw(product_raw_moments(n_cand[1], n_cand[2], n_cand[3],
                                                       a = 1000, r_max = 4))
  cat(sprintf("  analytical MaxEnt q05 at same candidate  : %s (deterministic)\n",
              fmt_num(me_q05_fixed)))

  # -----------------------------------------------------------------------
  # (B) SEEDED GA COMPARISON: 10 runs per method, matched seeds 1001..1010
  # -----------------------------------------------------------------------
  cat("\n---------------------------------------------------------\n")
  cat("(B) SEEDED GA COMPARISON: 10 runs/method, seeds 1001-1010\n")
  cat("    pop_size = 30, generations = 20, MC N_mc = 5000\n")
  cat("---------------------------------------------------------\n")
  runs <- 10

  mc_costs <- numeric(runs); mc_times <- numeric(runs)
  mc_q05   <- numeric(runs); mc_viol  <- logical(runs)
  mc_sols  <- matrix(0L, nrow = runs, ncol = 3)

  cat("\n--- GA with Monte Carlo fitness (N_mc = 5000) ---\n")
  for (r in 1:runs) {
    set.seed(1000 + r)
    t_start <- Sys.time()
    res <- run_genetic_algorithm("MC", N_mc = 5000, pop_size = 30, generations = 20)
    t_end <- Sys.time()
    mc_times[r] <- as.numeric(difftime(t_end, t_start, units = "secs"))
    mc_costs[r] <- res$cost
    mc_q05[r]   <- res$true_q05
    mc_viol[r]  <- res$violation
    mc_sols[r, ] <- res$sol
    cat(sprintf("Seed %d: Sol=(%d, %d, %d), Cost=%.1f, True Q05=%.1f, Violated=%s, Time=%.2fs\n",
                1000 + r, res$sol[1], res$sol[2], res$sol[3], res$cost,
                res$true_q05, res$violation, mc_times[r]))
  }

  me_costs <- numeric(runs); me_times <- numeric(runs)
  me_q05   <- numeric(runs); me_viol  <- logical(runs)
  me_sols  <- matrix(0L, nrow = runs, ncol = 3)

  cat("\n--- GA with Analytical MaxEnt fitness ---\n")
  for (r in 1:runs) {
    set.seed(1000 + r)  # same seed as the MC arm -> identical initial population
    t_start <- Sys.time()
    res <- run_genetic_algorithm("MaxEnt", N_mc = 5000, pop_size = 30, generations = 20)
    t_end <- Sys.time()
    me_times[r] <- as.numeric(difftime(t_end, t_start, units = "secs"))
    me_costs[r] <- res$cost
    me_q05[r]   <- res$true_q05
    me_viol[r]  <- res$violation
    me_sols[r, ] <- res$sol
    cat(sprintf("Seed %d: Sol=(%d, %d, %d), Cost=%.1f, True Q05=%.1f, Violated=%s, Time=%.2fs\n",
                1000 + r, res$sol[1], res$sol[2], res$sol[3], res$cost,
                res$true_q05, res$violation, me_times[r]))
  }

  cat("\n--- SUMMARY (B) ---\n")
  cat(sprintf("Metric                  | Monte Carlo (N=5k) | Analytical MaxEnt\n"))
  cat("-----------------------------------------------------------------\n")
  cat(sprintf("Mean optimized cost     | %18.2f | %17.2f\n", mean(mc_costs), mean(me_costs)))
  cat(sprintf("SD optimized cost       | %18.2f | %17.2f\n", sd(mc_costs), sd(me_costs)))
  cat(sprintf("Min / Max cost          | %8.1f / %7.1f | %7.1f / %7.1f\n",
              min(mc_costs), max(mc_costs), min(me_costs), max(me_costs)))
  cat(sprintf("Mean wall time (s)      | %18.2f | %17.2f\n", mean(mc_times), mean(me_times)))
  cat(sprintf("SD wall time (s)        | %18.2f | %17.2f\n", sd(mc_times), sd(me_times)))
  cat(sprintf("Constraint violations   | %12d / %d     | %11d / %d\n",
              sum(mc_viol), runs, sum(me_viol), runs))
  cat(sprintf("Whole-GA speedup        | mean(time_MC)/mean(time_MaxEnt) = %.2fx\n",
              mean(mc_times) / mean(me_times)))

  # -----------------------------------------------------------------------
  # (C) DECOMPOSITION: per-seed cost pairs (identical initial populations)
  # -----------------------------------------------------------------------
  cat("\n---------------------------------------------------------\n")
  cat("(C) DECOMPOSITION: per-seed cost pairs (MC vs MaxEnt)\n")
  cat("---------------------------------------------------------\n")
  cat("Seed  | Cost(MC)  | Cost(MaxEnt) | Diff(MC - MaxEnt)\n")
  cat("----------------------------------------------------\n")
  for (r in 1:runs) {
    cat(sprintf("%d  | %9.1f | %12.1f | %+10.1f\n",
                1000 + r, mc_costs[r], me_costs[r], mc_costs[r] - me_costs[r]))
  }
  cat("----------------------------------------------------\n")
  cat(sprintf("Mean per-seed diff (MC - MaxEnt): %+.2f\n", mean(mc_costs - me_costs)))
  cat(sprintf("SD   per-seed diff              : %.2f\n", sd(mc_costs - me_costs)))
  cat("\nInterpretation aid:\n")
  cat(sprintf("  SD within MaxEnt arm (fitness deterministic) = %.2f  <- GA stochasticity alone\n",
              sd(me_costs)))
  cat(sprintf("  SD within MC arm (fitness noisy + GA random) = %.2f\n", sd(mc_costs)))
  cat(sprintf("  At matched seeds the initial populations are identical;\n"))
  cat(sprintf("  per-seed cost differences reflect the fitness evaluator\n"))
  cat(sprintf("  (and downstream divergence of the GA RNG stream it induces).\n"))

  cat("\nDone.\n")
}

main()

#!/usr/bin/env Rscript

# =========================================================================
# FM-MEM (GA-optimized free exponents) vs the PATP one-parameter scan
# Reviewer Q4 / Weakness W5.
# =========================================================================
# The reviewer asks for a head-to-head against a recent fractional-moment
# MaxEnt that JOINTLY optimizes exponents and multipliers (GA-based), reporting
# accuracy AND runtime/robustness. We implement exactly that and compare on the
# bimodal mixture (the case where exponent selection is the lever):
#
#   (1) Monomial MaxEnt (alpha=1)                  -- classical baseline.
#   (2) PATP alpha-scan: 11 configs tied to one    -- OURS. ~11 inner solves,
#       scalar alpha, Gamma-selected.                 deterministic.
#   (3) FM-MEM-GA (Gamma objective): outer GA over  -- SOTA-style. The realizable
#       5 free parity-matched exponents, inner         selector both methods can
#       MaxEnt for the multipliers, select by the      use (no ground truth).
#       same dual-potential criterion.
#   (4) FM-MEM-GA (oracle MSE objective): GA driven -- UPPER BOUND on what free
#       by the true reconstruction MSE.                exponents can buy at all.
#
# Honest question: does paying ~25x the solver calls for a full GA exponent
# search beat the 11-point scan in ACCURACY? (4) bounds the best achievable.
# Solver/basis copied verbatim. Base R only; self-contained GA (no deps).
# =========================================================================

patp_power <- function(i, alpha) 1.0/i + (4.0 - i - 3.0/i)*alpha + (2.0*i - 4.0 + 2.0/i)*alpha^2
# basis with arbitrary per-column exponents (col 1 = x; parity matched by column)
basis_free <- function(x, exps) {           # exps: length-(n-1) vector for cols 2..n
  n <- length(exps) + 1
  M <- matrix(0, nrow = length(x), ncol = n); M[, 1] <- x
  for (i in 2:n) { p <- exps[i-1]
    M[, i] <- if (i %% 2 == 0) abs(x)^p else sign(x) * (abs(x)^p) }
  M
}
basis_patp <- function(x, n, alpha) basis_free(x, sapply(2:n, function(i) patp_power(i, alpha)))
solve_maxent <- function(Phi, target_moments, dx, max_iter = 100, tol = 1e-6) {
  n <- ncol(Phi); lambdas <- rep(0, n)
  for (iter in 1:max_iter) {
    unnorm_pdf <- exp(as.numeric(Phi %*% lambdas))
    if (any(is.infinite(unnorm_pdf)) || any(is.na(unnorm_pdf))) return(NULL)
    Z <- sum(unnorm_pdf) * dx; if (Z == 0 || is.na(Z)) return(NULL)
    pdf <- unnorm_pdf / Z; fitted <- as.numeric(t(Phi) %*% pdf * dx)
    grad <- fitted - target_moments
    if (sqrt(sum(grad^2)) < tol)
      return(list(lambdas = lambdas, Z = Z, converged = TRUE))
    H <- t(Phi) %*% (Phi * pdf) * dx - fitted %*% t(fitted); diag(H) <- diag(H) + 1e-8
    step <- tryCatch(solve(H, grad), error = function(e) tryCatch(qr.solve(H, grad), error = function(e2) NULL))
    if (is.null(step)) return(NULL)
    a_s <- 1.0; ok <- FALSE
    for (ls in 1:10) {
      nl <- lambdas - a_s * step; nu <- exp(as.numeric(Phi %*% nl))
      if (!any(is.infinite(nu)) && !any(is.na(nu))) { nZ <- sum(nu) * dx
        if (nZ > 0 && !is.na(nZ) && log(nZ) - sum(nl*target_moments) <= log(Z) - sum(lambdas*target_moments) + 1e-4) {
          lambdas <- nl; ok <- TRUE; break } }
      a_s <- a_s * 0.5
    }
    if (!ok) lambdas <- lambdas - 0.1 * step
  }
  list(converged = FALSE)
}

# ---- problem: bimodal mixture (Table 2 setup) ----
N <- 1000; sd_true <- sqrt(0.32); n <- 6
xg <- seq(-5, 5, length.out = 1000); dx <- xg[2] - xg[1]
ptrue <- 0.4*dnorm(xg,-1,sd_true) + 0.6*dnorm(xg,1,sd_true)

# fit a given exponent vector; returns list(mse, Gamma, conv) and counts one solve
SOLVES <- new.env(); SOLVES$n <- 0
fit_exps <- function(d, exps) {
  SOLVES$n <- SOLVES$n + 1
  tg <- colMeans(basis_free(d, exps)); Pg <- basis_free(xg, exps)
  f <- solve_maxent(Pg, tg, dx)
  if (is.null(f) || !f$converged) return(list(conv = FALSE, mse = Inf, Gamma = Inf))
  pdf <- exp(as.numeric(Pg %*% f$lambdas)) / f$Z
  list(conv = TRUE, mse = mean((pdf - ptrue)^2), Gamma = log(f$Z) - sum(f$lambdas * tg))
}

# ---- self-contained real-coded GA over the 5 free exponents (cols 2..6) ----
ga_optimize <- function(d, objective, pop = 16, gen = 12, lo = 0.1, hi = 3.0, seed = 1) {
  set.seed(1000 + seed)
  D <- n - 1
  P <- matrix(runif(pop * D, lo, hi), pop, D)
  score <- function(e) { r <- fit_exps(d, e); if (!r$conv) return(Inf)
                         if (objective == "gamma") r$Gamma else r$mse }
  fit <- apply(P, 1, score)
  for (g in 1:gen) {
    newP <- P
    for (k in 1:pop) {
      a <- sample(pop, 2); b <- sample(pop, 2)               # tournament parents
      p1 <- P[a[which.min(fit[a])], ]; p2 <- P[b[which.min(fit[b])], ]
      ch <- runif(D, pmin(p1,p2) - 0.25*abs(p1-p2), pmax(p1,p2) + 0.25*abs(p1-p2))  # BLX-0.25
      m <- runif(D) < 0.2; ch[m] <- ch[m] + rnorm(sum(m), 0, 0.3)  # Gaussian mutation
      newP[k, ] <- pmin(pmax(ch, lo), hi)
    }
    newFit <- apply(newP, 1, score)
    elite <- which.min(fit)                                  # elitism
    comb <- rbind(P, newP); combFit <- c(fit, newFit)
    keep <- order(combFit)[1:pop]
    P <- comb[keep, , drop = FALSE]; fit <- combFit[keep]
  }
  best <- which.min(fit); list(exps = P[best, ], score = fit[best])
}

# ---- run all four methods over seeds ----
seeds <- 1:5
res <- data.frame()
t_scan <- t_gaG <- t_gaM <- 0
for (s in seeds) {
  set.seed(s); u <- runif(N)
  d <- ifelse(u < 0.4, rnorm(N,-1,sd_true), rnorm(N,1,sd_true))

  # (1) monomial
  rm_ <- fit_exps(d, 2:n)                                    # exponents 2,3,4,5,6
  # (2) PATP alpha-scan, Gamma-select
  SOLVES$n <- 0; t0 <- Sys.time()
  scan <- data.frame()
  for (alpha in setdiff(seq(0,1,by=0.1), 0.5)) {
    r <- fit_exps(d, sapply(2:n, function(i) patp_power(i, alpha)))
    if (r$conv) scan <- rbind(scan, data.frame(alpha = alpha, mse = r$mse, Gamma = r$Gamma))
  }
  sel <- scan[which.min(scan$Gamma), ]; scan_solves <- SOLVES$n
  t_scan <- t_scan + as.numeric(difftime(Sys.time(), t0, units = "secs"))
  # (3) FM-MEM-GA, Gamma objective
  SOLVES$n <- 0; t0 <- Sys.time()
  gG <- ga_optimize(d, "gamma", seed = s); rG <- fit_exps(d, gG$exps); gaG_solves <- SOLVES$n
  t_gaG <- t_gaG + as.numeric(difftime(Sys.time(), t0, units = "secs"))
  # (4) FM-MEM-GA, oracle MSE objective (upper bound)
  SOLVES$n <- 0; t0 <- Sys.time()
  gM <- ga_optimize(d, "mse", seed = s); rM <- fit_exps(d, gM$exps); gaM_solves <- SOLVES$n
  t_gaM <- t_gaM + as.numeric(difftime(Sys.time(), t0, units = "secs"))

  oracle <- min(scan$mse)                                   # best achievable on the scan grid
  res <- rbind(res, data.frame(seed = s,
    mono = rm_$mse, scan_sel = sel$mse, scan_oracle = oracle,
    gaG = rG$mse, gaM = rM$mse,
    scan_solves = scan_solves, gaG_solves = gaG_solves, gaM_solves = gaM_solves))
}

cat("\n=========== FM-MEM-GA vs PATP alpha-scan (mixture, 5 seeds) ===========\n")
cat("Reconstruction PDF MSE (lower better):\n\n")
cat(sprintf("%5s %10s %10s %11s %10s %10s\n", "seed","monomial","scan(Gamma)","scan-oracle","FM-GA(Gam)","FM-GA(MSE)"))
for (i in 1:nrow(res)) with(res[i,],
  cat(sprintf("%5d %10.2e %10.2e %11.2e %10.2e %10.2e\n", seed, mono, scan_sel, scan_oracle, gaG, gaM)))
cat(sprintf("\n%5s %10.2e %10.2e %11.2e %10.2e %10.2e\n", "mean",
            mean(res$mono), mean(res$scan_sel), mean(res$scan_oracle), mean(res$gaG), mean(res$gaM)))
cat(sprintf("%5s %10s %10s %11s %10s %10s  (sd)\n", "", "",
            sprintf("%.2e", sd(res$scan_sel)), sprintf("%.2e", sd(res$scan_oracle)),
            sprintf("%.2e", sd(res$gaG)), sprintf("%.2e", sd(res$gaM))))
cat(sprintf("\nMean inner MaxEnt solves per fit:  scan = %.0f,  FM-GA(Gamma) = %.0f,  FM-GA(MSE) = %.0f\n",
            mean(res$scan_solves), mean(res$gaG_solves), mean(res$gaM_solves)))
cat(sprintf("Total wall time (5 seeds):         scan = %.2fs, FM-GA(Gamma) = %.2fs, FM-GA(MSE) = %.2fs\n",
            t_scan, t_gaG, t_gaM))
cat(sprintf("Cost ratio (solves): FM-GA(Gamma)/scan = %.1fx ; FM-GA(MSE)/scan = %.1fx\n",
            mean(res$gaG_solves)/mean(res$scan_solves), mean(res$gaM_solves)/mean(res$scan_solves)))
cat("\nInterpretation: scan-oracle vs FM-GA(MSE) bounds what free exponents can buy;\n")
cat("scan(Gamma) vs FM-GA(Gamma) is the apples-to-apples realizable comparison.\n")
cat("\nDone.\n")

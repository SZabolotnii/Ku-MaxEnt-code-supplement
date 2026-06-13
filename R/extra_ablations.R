#!/usr/bin/env Rscript

# Extra ablations requested by review gates:
# (A) Linear exponent map baseline: q_i(alpha) = (1-alpha)/i + alpha*i,
#     parity-matched basis — does the PATP *quadratic* map matter,
#     or does any one-parameter fractional family do?
# (B) 20-seed replication of the mixture headline numbers (PM-PATP scan vs
#     monomial vs T-MaxEnt) — converts single-seed case study into mean +/- sd.
# (C) 20-seed replication of the Cauchy headline numbers (alpha=0 KS/bodyMSE,
#     convergence region, monomial infeasibility frequency).
# Solver copied verbatim from patp_maxent_simulation.R.

patp_power <- function(i, alpha) {
  1.0 / i + (4.0 - i - 3.0 / i) * alpha + (2.0 * i - 4.0 + 2.0 / i) * alpha^2
}
lin_power <- function(i, alpha) (1 - alpha) / i + alpha * i

basis_pm <- function(x, n, alpha, pow_fun) {
  M <- matrix(0, nrow = length(x), ncol = n)
  M[, 1] <- x
  if (n > 1) for (i in 2:n) {
    p <- pow_fun(i, alpha)
    M[, i] <- if (i %% 2 == 0) abs(x)^p else sign(x) * (abs(x)^p)
  }
  M
}

trig_basis <- function(x, S, p) {
  M <- matrix(0, nrow = length(x), ncol = 2 * S)
  for (r in 1:S) { M[, 2*r-1] <- cos(r*p*x); M[, 2*r] <- sin(r*p*x) }
  M
}

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

fit_mix <- function(data_mix, xg, dxm, pdf_true, n, alpha, pow_fun) {
  tm <- colMeans(basis_pm(data_mix, n, alpha, pow_fun))
  Phi_g <- basis_pm(xg, n, alpha, pow_fun)
  fit <- solve_maxent(Phi_g, tm, dxm)
  if (is.null(fit)) return(NULL)
  pdf_f <- exp(as.numeric(Phi_g %*% fit$lambdas)) / fit$Z
  list(mse = mean((pdf_f - pdf_true)^2),
       pot = log(fit$Z) - sum(fit$lambdas * tm),
       kH = fit$cond_num, conv = fit$converged)
}

# =========================================================================
# (A) LINEAR EXPONENT MAP vs PATP QUADRATIC MAP — mixture, seed 20260612
# =========================================================================
cat("==========================================================\n")
cat("(A) LINEAR MAP q_i(a)=(1-a)/i + a*i  vs  PATP QUADRATIC MAP\n")
cat("==========================================================\n")
set.seed(20260612)
N <- 1000; sd_true <- sqrt(0.32)
u <- runif(N)
data_mix <- ifelse(u < 0.4, rnorm(N, -1, sd_true), rnorm(N, 1, sd_true))
xg <- seq(-5, 5, length.out = 1000); dxm <- xg[2] - xg[1]
pdf_true <- 0.4 * dnorm(xg, -1, sd_true) + 0.6 * dnorm(xg, 1, sd_true)

cat("\nmap | alpha | conv | kappa_H | potential | PDF MSE\n")
for (mp in c("PATP", "LIN")) {
  pf <- if (mp == "PATP") patp_power else lin_power
  res <- data.frame()
  for (alpha in seq(0, 1, by = 0.1)) {
    f <- fit_mix(data_mix, xg, dxm, pdf_true, 6, alpha, pf)
    if (is.null(f)) { cat(sprintf("%s | %.1f | solver-NULL\n", mp, alpha)); next }
    res <- rbind(res, data.frame(alpha = alpha, conv = f$conv, kH = f$kH, pot = f$pot, mse = f$mse))
    cat(sprintf("%s | %.1f | %s | %.3e | %.4f | %.6e\n", mp, alpha, f$conv, f$kH, f$pot, f$mse))
  }
  sel <- res[which.min(res$pot), ]
  best <- res[which.min(res$mse), ]
  cat(sprintf(">> %s map: selected alpha*=%.1f (MSE=%.4e, kH=%.3e); best alpha=%.1f (MSE=%.4e)\n\n",
              mp, sel$alpha, sel$mse, sel$kH, best$alpha, best$mse))
}

# =========================================================================
# (B) 20-SEED REPLICATION — MIXTURE
# =========================================================================
cat("==========================================================\n")
cat("(B) 20-SEED REPLICATION, MIXTURE (6 constraints, [-5,5])\n")
cat("==========================================================\n")
alpha_grid <- setdiff(seq(0, 1, by = 0.1), 0.5)  # exclude degenerate point
rep_m <- data.frame()
for (s in 1:20) {
  set.seed(s)
  u <- runif(N)
  dm <- ifelse(u < 0.4, rnorm(N, -1, sd_true), rnorm(N, 1, sd_true))
  # PM-PATP scan with potential selection
  res <- data.frame()
  for (alpha in alpha_grid) {
    f <- fit_mix(dm, xg, dxm, pdf_true, 6, alpha, patp_power)
    if (!is.null(f) && f$conv) res <- rbind(res, data.frame(alpha = alpha, kH = f$kH, pot = f$pot, mse = f$mse))
  }
  if (nrow(res) == 0) next
  sel <- res[which.min(res$pot), ]
  # monomial
  fm <- fit_mix(dm, xg, dxm, pdf_true, 6, 1.0, patp_power)
  # T-MaxEnt (0.5, 3)
  tmv <- colMeans(trig_basis(dm, 3, 0.5))
  Phi_t <- trig_basis(xg, 3, 0.5)
  ft <- solve_maxent(Phi_t, tmv, dxm)
  mse_t <- if (!is.null(ft)) mean((exp(as.numeric(Phi_t %*% ft$lambdas)) / ft$Z - pdf_true)^2) else NA
  rep_m <- rbind(rep_m, data.frame(seed = s, a_sel = sel$alpha, mse_sel = sel$mse, kH_sel = sel$kH,
                                   mse_mono = fm$mse, kH_mono = fm$kH, mse_trig = mse_t,
                                   ratio = fm$mse / sel$mse))
  cat(sprintf("seed %2d: a*=%.1f, MSE_sel=%.3e, MSE_mono=%.3e, MSE_trig=%.3e, ratio mono/sel=%.2f\n",
              s, sel$alpha, sel$mse, fm$mse, mse_t, fm$mse / sel$mse))
}
cat(sprintf("\nSummary over %d seeds:\n", nrow(rep_m)))
cat(sprintf("selected alpha*: median %.2f, range [%.1f, %.1f]; mode tally:\n",
            median(rep_m$a_sel), min(rep_m$a_sel), max(rep_m$a_sel)))
print(table(rep_m$a_sel))
cat(sprintf("MSE selected : mean %.3e, sd %.3e\n", mean(rep_m$mse_sel), sd(rep_m$mse_sel)))
cat(sprintf("MSE monomial : mean %.3e, sd %.3e\n", mean(rep_m$mse_mono), sd(rep_m$mse_mono)))
cat(sprintf("MSE T-MaxEnt : mean %.3e, sd %.3e\n", mean(rep_m$mse_trig), sd(rep_m$mse_trig)))
cat(sprintf("ratio mono/selected: mean %.2f, sd %.2f, min %.2f, max %.2f\n",
            mean(rep_m$ratio), sd(rep_m$ratio), min(rep_m$ratio), max(rep_m$ratio)))
cat(sprintf("kappa_H selected: geometric mean %.3e; monomial: %.3e\n",
            exp(mean(log(rep_m$kH_sel))), exp(mean(log(rep_m$kH_mono)))))

# =========================================================================
# (C) 20-SEED REPLICATION — CAUCHY
# =========================================================================
cat("\n==========================================================\n")
cat("(C) 20-SEED REPLICATION, CAUCHY (4 constraints, [-50,50])\n")
cat("==========================================================\n")
xgc <- seq(-50, 50, length.out = 2000); dxc <- xgc[2] - xgc[1]
rep_c <- data.frame()
for (s in 1:20) {
  set.seed(s)
  dc <- rcauchy(N)
  # alpha = 0 fit
  tm0 <- colMeans(basis_pm(dc, 4, 0, patp_power))
  Phi0 <- basis_pm(xgc, 4, 0, patp_power)
  f0 <- solve_maxent(Phi0, tm0, dxc)
  ks0 <- NA; bm0 <- NA
  if (!is.null(f0) && f0$converged) {
    pdf0 <- exp(as.numeric(Phi0 %*% f0$lambdas)) / f0$Z
    cdf0 <- cumsum(pdf0) * dxc
    ki <- which(xgc >= -45 & xgc <= 45)
    ks0 <- max(abs(cdf0[ki] - pcauchy(xgc[ki])))
    bi <- which(xgc >= -10 & xgc <= 10)
    bm0 <- mean((pdf0[bi] - dcauchy(xgc[bi]))^2)
  }
  # convergence region size and monomial infeasibility
  nconv <- 0
  for (alpha in c(0, 0.1, 0.2, 0.3)) {
    tma <- colMeans(basis_pm(dc, 4, alpha, patp_power))
    fa <- solve_maxent(basis_pm(xgc, 4, alpha, patp_power), tma, dxc)
    if (!is.null(fa) && fa$converged) nconv <- nconv + 1
  }
  # Monomial 4-moment MaxEnt feasibility certificate (matches the paper, sec:exp1):
  # the empirical target of the binding constraint must be attainable on the
  # truncated support D=[-L,L]. The binding monomial is x^4, so the problem is
  # INFEASIBLE iff the empirical 4th-moment target exceeds max_{x in D} x^4 = L^4.
  # (The 2nd-moment target m2_hat > L^2 = 2500 is a weaker, secondary witness and
  #  is reported alongside for transparency; the paper's "19/20" uses the x^4 cert.)
  L <- 50; max_x4 <- L^4            # 6.25e6
  m2 <- mean(dc^2)
  m4 <- mean(dc^4)
  mono_infeasible <- (m4 > max_x4)  # paper's 4-moment attainability certificate
  rep_c <- rbind(rep_c, data.frame(seed = s, conv0 = !is.na(ks0), ks0 = ks0, bm0 = bm0,
                                   nconv03 = nconv, m2_hat = m2, m4_hat = m4,
                                   m2_witness = (m2 > L^2), mono_infeasible = mono_infeasible))
  cat(sprintf("seed %2d: a=0 conv=%s, KS=%.3f, bodyMSE=%.3e, conv in {0..0.3}: %d/4, m2_hat=%.0f, m4_hat=%.3e (mono infeasible m4>%.2e: %s)\n",
              s, !is.na(ks0), ks0, bm0, nconv, m2, m4, max_x4, mono_infeasible))
}
cat(sprintf("\nSummary over 20 seeds:\n"))
cat(sprintf("alpha=0 converged: %d/20; KS: mean %.3f, sd %.3f; bodyMSE: mean %.3e, sd %.3e\n",
            sum(rep_c$conv0), mean(rep_c$ks0, na.rm = TRUE), sd(rep_c$ks0, na.rm = TRUE),
            mean(rep_c$bm0, na.rm = TRUE), sd(rep_c$bm0, na.rm = TRUE)))
cat(sprintf("mean # converged alphas in {0,0.1,0.2,0.3}: %.2f / 4\n", mean(rep_c$nconv03)))
cat(sprintf("4-moment monomial infeasible (m4_hat > L^4 = %.2e): %d / 20 seeds\n", max_x4, sum(rep_c$mono_infeasible)))
cat(sprintf("  [secondary 2nd-moment witness m2_hat > L^2 = 2500: %d / 20 seeds]\n", sum(rep_c$m2_witness)))
cat("\nDone.\n")

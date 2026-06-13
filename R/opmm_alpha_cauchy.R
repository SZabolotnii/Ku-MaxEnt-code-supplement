#!/usr/bin/env Rscript

# =========================================================================
# oPMMalpha — variance-optimal generating-element selection, Cauchy case
# =========================================================================
# Prototype answering the reviewer's central critique (Q1/Q6, Weakness W4):
# the dual-potential (Gamma) selection rule is heuristic and mis-picks on
# heavy tails (it selects alpha*=0.3, the LEAST accurate converged member).
#
# Proposed replacement (oPMMalpha): pick the generating element that minimizes
# the asymptotic variance of the *reported* functional T, via the PMM/delta
# method, using quantities the solver already forms:
#
#   sqrt(N)(That - T) -> N(0, V(alpha)),  V(alpha) = gradT' H^{-1} Sigma H^{-1} gradT
#
#   H(alpha)     = Cov_f[phi]   (the dual Hessian)
#   Sigma(alpha) = Cov_data[phi] (sampling cov of the empirical constraints)
#   gradT        = d T / d lambda  (delta method)
#
# Mechanism: for monomial-like high-exponent members on Cauchy, Var[phi]
# explodes (Sigma large) -> V large -> the criterion rejects them; it targets
# estimator precision of the reported quantity, not fitted entropy.
#
# Solver / basis copied verbatim from extra_ablations.R (-> patp_maxent_simulation.R).
# Seed and setup match Table 1 (sec:exp1): N=1000 Cauchy, D=[-50,50], 4 constraints.
# =========================================================================

patp_power <- function(i, alpha) {
  1.0 / i + (4.0 - i - 3.0 / i) * alpha + (2.0 * i - 4.0 + 2.0 / i) * alpha^2
}
basis_pm <- function(x, n, alpha, pow_fun) {
  M <- matrix(0, nrow = length(x), ncol = n)
  M[, 1] <- x
  if (n > 1) for (i in 2:n) {
    p <- pow_fun(i, alpha)
    M[, i] <- if (i %% 2 == 0) abs(x)^p else sign(x) * (abs(x)^p)
  }
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
  list(lambdas = lambdas, Z = Z, converged = FALSE)
}

# ---- setup (matches Table 1) --------------------------------------------
set.seed(20260612)
N  <- 1000
L  <- 50
n  <- 4
xg <- seq(-L, L, length.out = 2000); dx <- xg[2] - xg[1]
d  <- rcauchy(N)
f_true <- dcauchy(xg); F_true <- pcauchy(xg)
grid_alpha <- c(0, 0.1, 0.2, 0.3)   # converged members for this seed (Table 1)

# delta-method gradients of reported functionals wrt lambda --------------
# coverage probability  T = P(|X| <= c) = F(c)-F(-c):  dT/dl_k = int_{-c}^{c}(phi_k-mu)f dx
grad_coverage <- function(Phi_g, pdf, mu_f, c) {
  w <- (xg >= -c & xg <= c)
  centered <- sweep(Phi_g, 2, mu_f)            # phi - mu_f
  as.numeric(t(centered[w, , drop = FALSE]) %*% pdf[w] * dx)
}
# upper quantile q_beta:  dq/dl_k = -(1/f(q)) int_{-L}^{q}(phi_k-mu)f dx
grad_quantile <- function(Phi_g, pdf, cdf, mu_f, beta) {
  qi <- which(cdf >= beta)[1]
  if (is.na(qi) || pdf[qi] <= 0) return(NULL)
  centered <- sweep(Phi_g, 2, mu_f)
  cumint <- apply(centered * pdf, 2, function(col) sum(col[1:qi]) * dx)
  list(g = -cumint / pdf[qi], q = xg[qi])
}

rows <- list()
for (alpha in grid_alpha) {
  Phi_s  <- basis_pm(d,  n, alpha, patp_power)   # constraints on the SAMPLE (N x n)
  target <- colMeans(Phi_s)                      # full-sample empirical targets
  Phi_g  <- basis_pm(xg, n, alpha, patp_power)   # constraints on the grid
  fit <- solve_maxent(Phi_g, target, dx)
  if (is.null(fit) || !fit$converged) next

  lam <- fit$lambdas
  pdf <- exp(as.numeric(Phi_g %*% lam)) / fit$Z
  cdf <- cumsum(pdf) * dx
  mu_f <- as.numeric(t(Phi_g) %*% pdf * dx)                       # E_f[phi]
  H    <- t(Phi_g) %*% (Phi_g * pdf) * dx - mu_f %*% t(mu_f)      # Cov_f[phi] = dual Hessian
  diag(H) <- diag(H) + 1e-8
  Hinv <- solve(H)
  Sigma <- cov(Phi_s)                                            # Cov_data[phi] (sampling cov)
  M <- Hinv %*% Sigma %*% Hinv                                   # sandwich middle

  # accuracy truth + old criterion
  bi <- which(xg >= -10 & xg <= 10); bodyMSE <- mean((pdf[bi] - f_true[bi])^2)
  ki <- which(xg >= -45 & xg <= 45); KS <- max(abs(cdf[ki] - F_true[ki]))
  Gamma <- log(fit$Z) - sum(lam * target)

  # variance of reported functionals (per the N samples)
  g_cov10 <- grad_coverage(Phi_g, pdf, mu_f, 10)
  V_cov10 <- as.numeric(t(g_cov10) %*% M %*% g_cov10) / N
  q975 <- grad_quantile(Phi_g, pdf, cdf, mu_f, 0.975)
  V_q975 <- if (is.null(q975)) NA else as.numeric(t(q975$g) %*% M %*% q975$g) / N
  V_A <- sum(diag(M)) / N                                        # A-optimality (target-free)

  rows[[as.character(alpha)]] <- data.frame(
    alpha = alpha, Gamma = Gamma, bodyMSE = bodyMSE, KS = KS,
    trSigma = sum(diag(Sigma)), V_A = V_A, V_cov10 = V_cov10, V_q975 = V_q975)
}
res <- do.call(rbind, rows)

cat("\n================ oPMMalpha on Cauchy (seed 20260612) ================\n")
print(format(res, digits = 4, scientific = TRUE), row.names = FALSE)

pick <- function(v) res$alpha[which.min(v)]
cat(sprintf("\nSelection picks:\n"))
cat(sprintf("  Gamma-min (OLD heuristic)      : alpha* = %.1f   [Gamma=%.3f]\n",
            pick(res$Gamma), min(res$Gamma)))
cat(sprintf("  body-MSE truth (most accurate) : alpha  = %.1f\n", pick(res$bodyMSE)))
cat(sprintf("  oPMMalpha V_A  (A-optimal)     : alpha* = %.1f\n", pick(res$V_A)))
cat(sprintf("  oPMMalpha V_cov10 (P|X|<10)    : alpha* = %.1f\n", pick(res$V_cov10)))
cat(sprintf("  oPMMalpha V_q975 (q_0.975)     : alpha* = %.1f\n", pick(res$V_q975)))

# normalized (best member = 1.0) for readability
norm_tbl <- data.frame(alpha = res$alpha,
                       Gamma_rank = res$Gamma - min(res$Gamma),
                       bodyMSE_x = res$bodyMSE / min(res$bodyMSE),
                       V_A_x = res$V_A / min(res$V_A),
                       V_cov10_x = res$V_cov10 / min(res$V_cov10),
                       V_q975_x = res$V_q975 / min(res$V_q975, na.rm = TRUE))
cat("\nRelative (best converged member = 1.0; Gamma shown as excess over min):\n")
print(format(norm_tbl, digits = 3), row.names = FALSE)

# =========================================================================
# 20-SEED REPLICATION: how often does each criterion select the member that
# is actually most accurate (lowest body MSE), and what accuracy penalty does
# its pick incur? Converts the single-seed result into a replication claim.
# =========================================================================
cat("\n================ 20-seed replication of the selection rule ================\n")
agree_G <- agree_VA <- agree_Vc <- 0; nseed <- 0
pen_G <- pen_VA <- pen_Vc <- c()
for (s in 1:20) {
  set.seed(s); ds <- rcauchy(N)
  tab <- list()
  for (alpha in grid_alpha) {
    Ps <- basis_pm(ds, n, alpha, patp_power); tg <- colMeans(Ps)
    Pg <- basis_pm(xg, n, alpha, patp_power)
    ft <- solve_maxent(Pg, tg, dx); if (is.null(ft) || !ft$converged) next
    lm <- ft$lambdas; pf <- exp(as.numeric(Pg %*% lm)) / ft$Z
    mf <- as.numeric(t(Pg) %*% pf * dx)
    Hh <- t(Pg) %*% (Pg * pf) * dx - mf %*% t(mf); diag(Hh) <- diag(Hh) + 1e-8
    Hi <- tryCatch(solve(Hh), error = function(e) NULL); if (is.null(Hi)) next
    Mm <- Hi %*% cov(Ps) %*% Hi
    bi <- which(xg >= -10 & xg <= 10); bmse <- mean((pf[bi] - f_true[bi])^2)
    gc <- grad_coverage(Pg, pf, mf, 10)
    tab[[as.character(alpha)]] <- data.frame(alpha = alpha,
      Gamma = log(ft$Z) - sum(lm * tg), bodyMSE = bmse,
      V_A = sum(diag(Mm)) / N, V_cov10 = as.numeric(t(gc) %*% Mm %*% gc) / N)
  }
  T <- do.call(rbind, tab); if (is.null(T) || nrow(T) < 2) next
  nseed <- nseed + 1
  truth <- T$alpha[which.min(T$bodyMSE)]; best_mse <- min(T$bodyMSE)
  aG  <- T$alpha[which.min(T$Gamma)];   agree_G  <- agree_G  + (aG  == truth)
  aVA <- T$alpha[which.min(T$V_A)];     agree_VA <- agree_VA + (aVA == truth)
  aVc <- T$alpha[which.min(T$V_cov10)]; agree_Vc <- agree_Vc + (aVc == truth)
  pen_G  <- c(pen_G,  T$bodyMSE[T$alpha == aG]  / best_mse)
  pen_VA <- c(pen_VA, T$bodyMSE[T$alpha == aVA] / best_mse)
  pen_Vc <- c(pen_Vc, T$bodyMSE[T$alpha == aVc] / best_mse)
}
cat(sprintf("Over %d seeds with >=2 converged members:\n", nseed))
cat(sprintf("  Gamma-min            picks the most-accurate member: %2d/%d  | mean body-MSE penalty %.2fx\n",
            agree_G,  nseed, mean(pen_G)))
cat(sprintf("  oPMMa V_A (A-opt)    picks the most-accurate member: %2d/%d  | mean body-MSE penalty %.2fx\n",
            agree_VA, nseed, mean(pen_VA)))
cat(sprintf("  oPMMa V_cov10 (body) picks the most-accurate member: %2d/%d  | mean body-MSE penalty %.2fx\n",
            agree_Vc, nseed, mean(pen_Vc)))
cat("\nDone.\n")

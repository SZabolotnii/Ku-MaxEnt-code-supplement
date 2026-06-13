#!/usr/bin/env Rscript

# =========================================================================
# oPMMalpha — mixture (light-tailed) control: does the variance criterion
# AGREE with the dual-potential heuristic where the heuristic already works?
# =========================================================================
# The companion Cauchy prototype (opmm_alpha_cauchy.R) shows oPMMalpha FIXES
# the Gamma-minimization mis-pick on heavy tails. This script checks the other
# side of the claim: on the bimodal Gaussian mixture (Sigma finite, light tails)
# the paper reports Gamma-selection is reliable. A principled replacement must
# NOT break that. We confirm oPMMalpha (V targeting the reported coverage
# functional) tracks the most-accurate member and agrees with Gamma here.
#
# Setup matches Table 2 / extra_ablations.R section (B):
#   N=1000, mixture 0.4*N(-1,s) + 0.6*N(+1,s), s=sqrt(0.32), D=[-5,5], 6 constraints.
# Solver/basis copied verbatim from extra_ablations.R.
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
  unnorm_pdf <- exp(as.numeric(Phi %*% lambdas)); Z <- sum(unnorm_pdf) * dx
  list(lambdas = lambdas, Z = Z, converged = FALSE)
}

# coverage-probability functional T = P(a<=X<=b); dT/dl_k = int_a^b (phi_k - mu) f dx
grad_coverage <- function(Phi_g, xg, pdf, mu_f, dx, a, b) {
  w <- (xg >= a & xg <= b)
  centered <- sweep(Phi_g, 2, mu_f)
  as.numeric(t(centered[w, , drop = FALSE]) %*% pdf[w] * dx)
}

# ---- setup (matches Table 2) --------------------------------------------
N <- 1000; sd_true <- sqrt(0.32); n <- 6
xg <- seq(-5, 5, length.out = 1000); dx <- xg[2] - xg[1]
pdf_true <- 0.4 * dnorm(xg, -1, sd_true) + 0.6 * dnorm(xg, 1, sd_true)
alpha_grid <- setdiff(seq(0, 1, by = 0.1), 0.5)   # exclude degenerate point

eval_seed <- function(seed) {
  set.seed(seed); u <- runif(N)
  dm <- ifelse(u < 0.4, rnorm(N, -1, sd_true), rnorm(N, 1, sd_true))
  half <- seq_len(N) %% 2 == 0           # deterministic 50/50 train/holdout split
  rows <- list()
  for (alpha in alpha_grid) {
    Ps <- basis_pm(dm, n, alpha, patp_power); tg <- colMeans(Ps)
    Pg <- basis_pm(xg, n, alpha, patp_power)
    ft <- solve_maxent(Pg, tg, dx); if (is.null(ft) || !ft$converged) next
    lm <- ft$lambdas; pf <- exp(as.numeric(Pg %*% lm)) / ft$Z
    mf <- as.numeric(t(Pg) %*% pf * dx)
    Hh <- t(Pg) %*% (Pg * pf) * dx - mf %*% t(mf); diag(Hh) <- diag(Hh) + 1e-8
    Hi <- tryCatch(solve(Hh), error = function(e) NULL); if (is.null(Hi)) next
    Mm <- Hi %*% cov(Ps) %*% Hi
    mse <- mean((pf - pdf_true)^2)
    gc <- grad_coverage(Pg, xg, pf, mf, dx, -3, 3)      # reported bulk-coverage functional
    # held-out predictive log-score (CV): fit on train half, score on holdout half.
    Pst <- basis_pm(dm[!half], n, alpha, patp_power); tgt <- colMeans(Pst)
    fcv <- solve_maxent(Pg, tgt, dx)
    nll <- NA
    if (!is.null(fcv) && fcv$converged) {
      Ph <- basis_pm(dm[half], n, alpha, patp_power)
      logf <- as.numeric(Ph %*% fcv$lambdas) - log(fcv$Z)   # exact log-density of fit
      nll <- -mean(logf[is.finite(logf)])                   # mean negative log-likelihood (lower=better)
    }
    rows[[as.character(alpha)]] <- data.frame(alpha = alpha,
      Gamma = log(ft$Z) - sum(lm * tg), pdfMSE = mse,
      V_cov = as.numeric(t(gc) %*% Mm %*% gc) / N, NLL = nll)
  }
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

# ---- single seed (Table 2 seed) -----------------------------------------
T0 <- eval_seed(20260612)
cat("\n================ oPMMalpha on the mixture (seed 20260612) ================\n")
print(format(T0, digits = 4, scientific = TRUE), row.names = FALSE)
cat(sprintf("\n  Gamma-min picks alpha* = %.1f ; pdf-MSE truth alpha = %.1f ; V_cov picks alpha* = %.1f ; CV-logscore picks alpha* = %.1f\n",
            T0$alpha[which.min(T0$Gamma)], T0$alpha[which.min(T0$pdfMSE)],
            T0$alpha[which.min(T0$V_cov)], T0$alpha[which.min(T0$NLL)]))

# ---- 20-seed: which criterion tracks shape accuracy on the mixture? -----
cat("\n================ 20-seed selection on the mixture (shape reconstruction) ================\n")
aG <- aVc <- aCV <- 0; nseed <- 0; pen_G <- pen_Vc <- pen_CV <- c()
for (s in 1:20) {
  T <- eval_seed(s); if (is.null(T) || nrow(T) < 2) next
  nseed <- nseed + 1
  truth <- T$alpha[which.min(T$pdfMSE)]; best <- min(T$pdfMSE)
  gG  <- T$alpha[which.min(T$Gamma)]; gVc <- T$alpha[which.min(T$V_cov)]
  Tcv <- T[is.finite(T$NLL), ]; gCV <- if (nrow(Tcv)) Tcv$alpha[which.min(Tcv$NLL)] else NA
  aG  <- aG  + (gG  == truth); aVc <- aVc + (gVc == truth); aCV <- aCV + (!is.na(gCV) && gCV == truth)
  pen_G  <- c(pen_G,  T$pdfMSE[T$alpha == gG]  / best)
  pen_Vc <- c(pen_Vc, T$pdfMSE[T$alpha == gVc] / best)
  if (!is.na(gCV)) pen_CV <- c(pen_CV, T$pdfMSE[T$alpha == gCV] / best)
}
cat(sprintf("Over %d seeds (most-accurate by pdf MSE):\n", nseed))
cat(sprintf("  Gamma-min               picks most-accurate: %2d/%d | mean penalty %.2fx\n", aG,  nseed, mean(pen_G)))
cat(sprintf("  oPMMa V_cov (saturated) picks most-accurate: %2d/%d | mean penalty %.2fx   <- WRONG target (P|X|<3 ~ 1)\n", aVc, nseed, mean(pen_Vc)))
cat(sprintf("  CV held-out log-score   picks most-accurate: %2d/%d | mean penalty %.2fx\n", aCV, nseed, mean(pen_CV)))
cat("\nDone.\n")

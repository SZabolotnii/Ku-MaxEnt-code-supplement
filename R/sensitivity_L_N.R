#!/usr/bin/env Rscript

# =========================================================================
# Sensitivity to truncation length L and sample size N (reviewer Q2 / W6)
# =========================================================================
# "How sensitive are the conclusions to L? Show feasibility, conditioning, and
#  tail/quantile errors as L varies, for Cauchy and mixture."
#
# Cauchy : fit the body member PM-PATP alpha=0 (4 constraints) on D=[-L,L];
#          report convergence rate, conditioning, body MSE, KS, q95/q99 error,
#          and the four-moment monomial infeasibility rate (m4_hat > L^4).
# Mixture: scan PM-PATP (6 constraints), Gamma-select; report conditioning,
#          PDF MSE, selected alpha*. (Light tails: expect near-insensitivity.)
#
# Fixed grid spacing dx (not fixed point count) so quantile resolution is
# comparable across L. 10 seeds per (L,N) cell. Solver/basis copied verbatim.
# =========================================================================

patp_power <- function(i, alpha) 1.0/i + (4.0 - i - 3.0/i)*alpha + (2.0*i - 4.0 + 2.0/i)*alpha^2
basis_pm <- function(x, n, alpha, pow_fun) {
  M <- matrix(0, nrow = length(x), ncol = n); M[, 1] <- x
  if (n > 1) for (i in 2:n) { p <- pow_fun(i, alpha)
    M[, i] <- if (i %% 2 == 0) abs(x)^p else sign(x) * (abs(x)^p) }
  M
}
solve_maxent <- function(Phi, target_moments, dx, max_iter = 100, tol = 1e-6) {
  n <- ncol(Phi); lambdas <- rep(0, n)
  for (iter in 1:max_iter) {
    unnorm_pdf <- exp(as.numeric(Phi %*% lambdas))
    if (any(is.infinite(unnorm_pdf)) || any(is.na(unnorm_pdf))) return(NULL)
    Z <- sum(unnorm_pdf) * dx; if (Z == 0 || is.na(Z)) return(NULL)
    pdf <- unnorm_pdf / Z; fitted <- as.numeric(t(Phi) %*% pdf * dx)
    grad <- fitted - target_moments
    if (sqrt(sum(grad^2)) < tol) {
      H <- t(Phi) %*% (Phi * pdf) * dx - fitted %*% t(fitted); diag(H) <- diag(H) + 1e-8
      return(list(lambdas = lambdas, Z = Z, cond_num = kappa(H, exact = TRUE), converged = TRUE))
    }
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
gm <- function(v) if (length(v)) exp(mean(log(v))) else NA   # geometric mean

NSEED <- 10

# ---------------- CAUCHY: vary L and N ----------------
cat("\n=================== CAUCHY sensitivity to (L, N) ===================\n")
cat("PM-PATP alpha=0, 4 constraints; dx fixed ~0.05; convergence/conditioning/body/KS/quantiles + monomial infeasibility\n\n")
cat(sprintf("%4s %5s | %5s %9s %10s %7s %8s %8s %8s\n",
            "L","N","conv","kH(gm)","bodyMSE","KS","q95err%","q99err%","mono.inf"))
for (L in c(20, 50, 100, 200)) {
  xg <- seq(-L, L, by = 0.05); dx <- xg[2]-xg[1]
  Phi0 <- basis_pm(xg, 4, 0, patp_power)
  ftrue <- dcauchy(xg); Ftrue <- pcauchy(xg)
  for (N in c(200, 1000, 5000)) {
    conv <- 0; kHs <- bms <- kss <- q95e <- q99e <- c(); monoinf <- 0
    for (s in 1:NSEED) {
      set.seed(s); d <- rcauchy(N)
      if (mean(d^4) > L^4) monoinf <- monoinf + 1
      tm <- colMeans(basis_pm(d, 4, 0, patp_power))
      f <- solve_maxent(Phi0, tm, dx)
      if (is.null(f) || !f$converged) next
      conv <- conv + 1; kHs <- c(kHs, f$cond_num)
      pdf <- exp(as.numeric(Phi0 %*% f$lambdas))/f$Z; cdf <- cumsum(pdf)*dx
      bi <- which(xg >= -10 & xg <= 10); bms <- c(bms, mean((pdf[bi]-ftrue[bi])^2))
      ki <- which(xg >= -0.9*L & xg <= 0.9*L); kss <- c(kss, max(abs(cdf[ki]-Ftrue[ki])))
      q95 <- xg[which(cdf >= 0.95)[1]]; q99 <- xg[which(cdf >= 0.99)[1]]
      q95e <- c(q95e, 100*(q95-6.314)/6.314); q99e <- c(q99e, 100*(q99-31.821)/31.821)
    }
    cat(sprintf("%4d %5d | %2d/%2d %9.1e %10.2e %7.3f %+8.0f %+8.0f %6d/%d\n",
                L, N, conv, NSEED, gm(kHs), mean(bms), mean(kss),
                mean(q95e), mean(q99e), monoinf, NSEED))
  }
  cat("\n")
}

# ---------------- MIXTURE: vary L and N ----------------
cat("=================== MIXTURE sensitivity to (L, N) ===================\n")
cat("PM-PATP scan (6 constraints), Gamma-selected; light tails -> expect near-insensitivity to L\n\n")
sd_true <- sqrt(0.32); alpha_grid <- setdiff(seq(0,1,by=0.1), 0.5)
cat(sprintf("%4s %5s | %5s %9s %10s %8s\n", "L","N","conv","kH(gm)","pdfMSE(sel)","a*(med)"))
for (L in c(3, 5, 8, 12)) {
  xg <- seq(-L, L, by = 0.01); dx <- xg[2]-xg[1]
  ptrue <- 0.4*dnorm(xg,-1,sd_true) + 0.6*dnorm(xg,1,sd_true)
  for (N in c(200, 1000, 5000)) {
    conv <- 0; kHs <- mss <- asel <- c()
    for (s in 1:NSEED) {
      set.seed(s); u <- runif(N)
      dm <- ifelse(u < 0.4, rnorm(N,-1,sd_true), rnorm(N,1,sd_true))
      best_pot <- Inf; best <- NULL; ba <- NA
      for (alpha in alpha_grid) {
        tm <- colMeans(basis_pm(dm, 6, alpha, patp_power))
        Pg <- basis_pm(xg, 6, alpha, patp_power)
        f <- solve_maxent(Pg, tm, dx)
        if (is.null(f) || !f$converged) next
        pot <- log(f$Z) - sum(f$lambdas*tm)
        if (pot < best_pot) { best_pot <- pot; pdf <- exp(as.numeric(Pg %*% f$lambdas))/f$Z
          best <- list(kH=f$cond_num, mse=mean((pdf-ptrue)^2)); ba <- alpha }
      }
      if (is.null(best)) next
      conv <- conv + 1; kHs <- c(kHs, best$kH); mss <- c(mss, best$mse); asel <- c(asel, ba)
    }
    cat(sprintf("%4d %5d | %2d/%2d %9.1e %10.2e %8.2f\n",
                L, N, conv, NSEED, gm(kHs), mean(mss), median(asel)))
  }
  cat("\n")
}
cat("Done.\n")

#!/usr/bin/env Rscript

# =========================================================================
# Unbounded log-rational element fit (reviewer Q3 / Weakness W3)
# =========================================================================
# The truncated-grid log-rational fit (log_rational_element.R) suffers from
# quantile saturation near the support edge +-L. Here we fit the SAME matched
# element on the WHOLE real line in closed form, with lambda constrained to the
# normalizable set lambda < -1/2 -- removing the truncation artifact entirely.
#
# Model:  f(x) = (1/Z) (1 + (x/s)^2)^lambda  on R,  Z = s*B(1/2, -lambda-1/2).
# This is the Student/Cauchy family: nu = -2*lambda - 1 degrees of freedom,
# exact tail slope 2*lambda, and X = (s/sqrt(nu)) * T_nu (standard Student-t).
#
# Single-constraint fit: match the one generalized moment
#   mu = E[log(1+(X/s)^2)].
# For the model, this has a CLOSED FORM (no quadrature, no grid):
#   E_f[log(1+(X/s)^2)] = psi(-lambda) - psi(-lambda - 1/2)      (digamma psi)
# Sanity: lambda=-1 (Cauchy) -> psi(1)-psi(1/2) = -gamma - (-gamma-2ln2) = 2 ln 2.
# So the fit is a 1-D root-find: solve  psi(-lambda)-psi(-lambda-1/2) = mu_hat.
#
# We then read EXACT analytic quantiles/CDF (qt/pt) -- no saturation possible.
# Seeds match the rest of the suite. Base R only.
# =========================================================================

s <- 1                                   # scale fixed a priori (Cauchy scale)
moment_of_lambda <- function(lam) digamma(-lam) - digamma(-lam - 0.5)

fit_unbounded <- function(x) {
  mu <- mean(log(1 + (x / s)^2))         # empirical generalized moment (finite for Cauchy)
  # psi(-lam)-psi(-lam-1/2) is positive and decreasing in lam on (-Inf,-1/2);
  # ranges from +Inf (lam->-1/2^-) down to 0 (lam->-Inf). Bracket and solve.
  f <- function(lam) moment_of_lambda(lam) - mu
  lam <- tryCatch(uniroot(f, lower = -200, upper = -0.5 - 1e-6, tol = 1e-10)$root,
                  error = function(e) NA)
  list(lambda = lam, mu = mu, nu = -2 * lam - 1)
}

# analytic CDF / quantile of the fitted law  X = (s/sqrt(nu)) T_nu
F_fit <- function(x, lam) { nu <- -2*lam - 1; pt(x * sqrt(nu) / s, df = nu) }
Q_fit <- function(b, lam) { nu <- -2*lam - 1; (s / sqrt(nu)) * qt(b, df = nu) }

# true-law references
targets <- list(
  Cauchy   = list(r = function(n) rcauchy(n),        F = function(x) pcauchy(x),
                  q = function(b) qcauchy(b),          slope = -2, lam_true = -1),
  Student2 = list(r = function(n) rt(n, df = 2),      F = function(x) pt(x, df = 2),
                  q = function(b) qt(b, df = 2),       slope = -3, lam_true = -1.5),
  Student3 = list(r = function(n) rt(n, df = 3),      F = function(x) pt(x, df = 3),
                  q = function(b) qt(b, df = 3),       slope = -4, lam_true = -2)
)

N <- 1000
cat("\n=========== Unbounded (closed-form) log-rational fit ===========\n")
cat("model  f(x) ~ (1+(x/s)^2)^lambda on R ;  fit:  psi(-l)-psi(-l-1/2) = mu_hat\n\n")
xref <- seq(-2000, 2000, length.out = 400001)   # fine grid for KS sup over R
for (nm in names(targets)) {
  tg <- targets[[nm]]
  set.seed(20260612); x <- tg$r(N)
  ft <- fit_unbounded(x)
  lam <- ft$lambda
  KS <- max(abs(F_fit(xref, lam) - tg$F(xref)))
  q95 <- Q_fit(0.95, lam); q99 <- Q_fit(0.99, lam)
  tq95 <- tg$q(0.95); tq99 <- tg$q(0.99)
  cat(sprintf("%-9s mu_hat=%.4f -> lambda_hat=%.4f (true %.3f), nu_hat=%.3f | tail slope %.3f (true %d)\n",
              nm, ft$mu, lam, tg$lam_true, ft$nu, 2*lam, tg$slope))
  cat(sprintf("          KS(fit,true) over R = %.4f | q95 %.3f (true %.3f, err %+.1f%%) | q99 %.3f (true %.3f, err %+.1f%%)\n\n",
              KS, q95, tq95, 100*(q95-tq95)/tq95, q99, tq99, 100*(q99-tq99)/tq99))
}

# ---- contrast with the truncated-grid fit on Cauchy (the saturation artifact) ----
cat("=========== Contrast: truncated-grid vs unbounded, Cauchy ===========\n")
# truncated grid fit (single log-rational constraint), as in log_rational_element.R
set.seed(20260612); xc <- rcauchy(N)
L <- 50; xg <- seq(-L, L, length.out = 4000); dx <- xg[2] - xg[1]
mu_hat <- mean(log(1 + xc^2))
phi <- log(1 + xg^2)                      # single constraint on the grid
# 1-D Newton on lambda for the truncated (grid-normalized) model
lam <- -1
for (it in 1:200) {
  w  <- exp(lam * phi); Z <- sum(w) * dx
  Ef <- sum(phi * w) * dx / Z
  Vf <- sum(phi^2 * w) * dx / Z - Ef^2
  step <- (Ef - mu_hat) / Vf; lam <- lam - step
  if (abs(step) < 1e-10) break
}
w <- exp(lam * phi); Z <- sum(w) * dx; pdf <- w / Z; cdf <- cumsum(pdf) * dx
ks_tr <- max(abs(cdf - pcauchy(xg)))
q95_tr <- xg[which(cdf >= 0.95)[1]]; q99_tr <- xg[which(cdf >= 0.99)[1]]
fu <- fit_unbounded(xc)
cat(sprintf("truncated [-50,50]: lambda=%.4f, slope=%.3f, KS=%.4f, q95=%.2f (true 6.31), q99=%.2f (true 31.82)\n",
            lam, 2*lam, ks_tr, q95_tr, q99_tr))
cat(sprintf("unbounded R       : lambda=%.4f, slope=%.3f, KS=%.4f, q95=%.2f, q99=%.2f  (exact analytic, no saturation)\n",
            fu$lambda, 2*fu$lambda, max(abs(F_fit(xref, fu$lambda) - pcauchy(xref))),
            Q_fit(0.95, fu$lambda), Q_fit(0.99, fu$lambda)))
cat("\nDone.\n")

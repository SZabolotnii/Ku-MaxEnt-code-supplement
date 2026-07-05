#!/usr/bin/env Rscript
# ============================================================================
# Head-to-head benchmark: Ku-MaxEnt (PATP generating-element) vs monomial-MaxEnt
# (GOPoly-accuracy proxy on known support) vs Pearson system, on the published
# multimodal test distributions and expanded-uncertainty error metric of
# Rajan et al. 2018 (IEEE Access), "Moment-Constrained Maximum Entropy Method
# for Expanded Uncertainty Evaluation".
#
# Protocol (matches Rajan 2018 Sec. III):
#   * EXACT high-order moments of each analytical test distribution (no sampling).
#   * Reconstruct PDF from n moments; estimate quantiles at tail percentile levels.
#   * Metric: expanded-uncertainty error eps = |x - x*| / |x - m1|, averaged over
#     the lower- and upper-tail percentile levels used in Rajan 2018.
#   * Competitors: Pearson (Rajan's stated best 4-moment method) + monomial MaxEnt.
#
# NOTE on the GOPoly proxy: on a fixed, known integration grid the moment-
# constrained MaxEnt density is the unique solution of a convex problem and is
# independent of the basis used to represent the same monomial-moment span.
# GOPoly (Rajan) orthogonalises the monomial basis for numerical conditioning and
# auto-estimates the support; on our KNOWN grid its converged density -- hence its
# eps accuracy -- equals that of the raw monomial-MaxEnt arm. We therefore report
# monomial-MaxEnt as the GOPoly *accuracy* proxy and flag any monomial non-
# convergence separately (that is the conditioning gap GOPoly is designed to close).
# ============================================================================

suppressWarnings(suppressMessages(library(PearsonDS)))

# ---- Ku-MaxEnt core (reused from scripts/patp_maxent_simulation.R) ----------
patp_power <- function(i, alpha) {
  A <- 1.0 / i; B <- 4.0 - i - 3.0 / i; C <- 2.0 * i - 4.0 + 2.0 / i
  A + B * alpha + C * alpha^2
}
patp_basis <- function(x, n, alpha) {
  M <- matrix(0, nrow = length(x), ncol = n)
  M[, 1] <- x
  if (n > 1) for (i in 2:n) M[, i] <- sign(x) * (abs(x)^patp_power(i, alpha))
  M
}
monomial_basis <- function(x, n) {
  M <- matrix(0, nrow = length(x), ncol = n)
  for (i in 1:n) M[, i] <- x^i
  M
}
trig_basis <- function(x, S, p) {
  # 2S trigonometric constraints: cos(r p x), sin(r p x), r=1..S
  # Correct Ku-MaxEnt element for bounded / multimodal / light-tailed targets.
  M <- matrix(0, nrow = length(x), ncol = 2 * S)
  for (r in 1:S) { M[, 2*r-1] <- cos(r*p*x); M[, 2*r] <- sin(r*p*x) }
  M
}
solve_maxent <- function(Phi, target_moments, dx, max_iter = 200, tol = 1e-7) {
  n <- ncol(Phi); lambdas <- rep(0, n)
  for (iter in 1:max_iter) {
    unnorm <- exp(as.numeric(Phi %*% lambdas))
    if (any(!is.finite(unnorm))) return(NULL)
    Z <- sum(unnorm) * dx
    if (!is.finite(Z) || Z <= 0) return(NULL)
    pdf <- unnorm / Z
    fm <- as.numeric(t(Phi) %*% pdf * dx)
    grad <- fm - target_moments
    H <- t(Phi) %*% (Phi * pdf) * dx - fm %*% t(fm); diag(H) <- diag(H) + 1e-9
    if (sqrt(sum(grad^2)) < tol) {
      return(list(lambdas = lambdas, Z = Z, converged = TRUE,
                  cond = kappa(H, exact = TRUE)))
    }
    step <- tryCatch(solve(H, grad), error = function(e)
      tryCatch(qr.solve(H, grad), error = function(e2) NULL))
    if (is.null(step)) return(NULL)
    a <- 1.0; ok <- FALSE
    for (ls in 1:20) {
      nl <- lambdas - a * step
      nu <- exp(as.numeric(Phi %*% nl))
      if (all(is.finite(nu))) {
        nZ <- sum(nu) * dx
        if (is.finite(nZ) && nZ > 0) {
          op <- log(Z) - sum(lambdas * target_moments)
          np <- log(nZ) - sum(nl * target_moments)
          if (np <= op + 1e-5) { lambdas <- nl; ok <- TRUE; break }
        }
      }
      a <- a * 0.5
    }
    if (!ok) lambdas <- lambdas - 0.05 * step
  }
  unnorm <- exp(as.numeric(Phi %*% lambdas)); Z <- sum(unnorm) * dx
  if (!is.finite(Z) || Z <= 0) return(NULL)
  pdf <- unnorm / Z; fm <- as.numeric(t(Phi) %*% pdf * dx)
  H <- t(Phi) %*% (Phi * pdf) * dx - fm %*% t(fm); diag(H) <- diag(H) + 1e-9
  list(lambdas = lambdas, Z = Z, converged = sqrt(sum((fm-target_moments)^2)) < 1e-4,
       cond = kappa(H, exact = TRUE))
}

# ---- Rajan 2018 multimodal test distributions (Appendix B, eqs 15-20) -------
# Each is given as an UNNORMALISED shape on a stated support; we normalise
# numerically, so exact leading constants are irrelevant.
norm_mix <- function(w, mu, s2) function(x) {
  out <- numeric(length(x))
  for (k in seq_along(w)) out <- out + w[k] * dnorm(x, mu[k], sqrt(s2[k]))
  out
}
dists <- list(
  D2_bimodal_wide  = list(f = norm_mix(c(.5,.5), c(-4,4),  c(16,9)),   lo=-20, hi=20),
  D3_bimodal_close = list(f = norm_mix(c(.4,.6), c(-1,1),  c(.16,.16)),lo=-4,  hi=4),
  D4_bimodal_skew  = list(f = norm_mix(c(.4,.6), c(0,3),   c(1,1)),    lo=-6,  hi=10),
  D5_bimodal_far   = list(f = norm_mix(c(.4,.6), c(0,25),  c(1,1)),    lo=-6,  hi=32),
  D1_beta_mix      = list(f = function(x){ b<-function(a,bb) dbeta(x,a,bb)
                            (b(8,1)+b(16,1)+b(64,0.5))/3 }, lo=0, hi=1),
  D6_beta_shift    = list(f = function(x){
                            ((x+1)^23*(1-x)^11 + (x+1)^11*(1-x)^23) }, lo=-1, hi=1)
)

pct_lower <- c(0.001,0.006,0.01,0.1,0.27,1,2,3,4,5,6,7,8,9,10)/100
pct_upper <- c(90,91,92,93,94,95,96,97,98,99,99.73,99.9,99.99,99.994)/100
levels_p  <- c(pct_lower, pct_upper)

# ---- helpers ---------------------------------------------------------------
standardize <- function(d, ngrid = 40000) {
  xf <- seq(d$lo, d$hi, length.out = ngrid); dxf <- xf[2]-xf[1]
  sh <- d$f(xf); sh[!is.finite(sh) | sh < 0] <- 0
  Zc <- sum(sh)*dxf; f <- sh/Zc
  mu <- sum(xf*f)*dxf; v <- sum((xf-mu)^2*f)*dxf; s <- sqrt(v)
  list(mu=mu, s=s, xf=xf, f=f, dxf=dxf)
}
true_quantiles <- function(st, p) {
  cdf <- cumsum(st$f)*st$dxf; cdf <- cdf/cdf[length(cdf)]
  zt <- approx(cdf, (st$xf-st$mu)/st$s, xout=p, ties="ordered", rule=2)$y
  zt
}
zgrid <- function(st, m=4000) {
  zlo <- (min(st$xf)-st$mu)/st$s; zhi <- (max(st$xf)-st$mu)/st$s
  seq(zlo, zhi, length.out=m)
}
z_moments_true <- function(st, basisfun) {
  # integrate basis against true density in z-space
  zf <- (st$xf-st$mu)/st$s
  # true density in z: f_z(z) = s * f_x(x); on xf grid weight = f*dxf already integrates in x,
  # and dz = dx/s, so integral g(z) f_z dz = sum(g(zf) * f * dxf)
  B <- basisfun(zf)
  as.numeric(t(B) %*% st$f * st$dxf)
}
q_from_fit <- function(zg, fit, basisfun, p) {
  Phi <- basisfun(zg); dz <- zg[2]-zg[1]
  pdf <- exp(as.numeric(Phi %*% fit$lambdas))/fit$Z
  cdf <- cumsum(pdf)*dz; cdf <- cdf/cdf[length(cdf)]
  approx(cdf, zg, xout=p, ties="ordered", rule=2)$y
}
eps_metric <- function(zt, zest) {
  # eps = |z - z*| / |z - m1|, m1(z)=0 ; average over levels (drop |zt|<0.05 to avoid /0)
  keep <- abs(zt) > 0.05 & is.finite(zest)
  mean(abs(zt[keep]-zest[keep])/abs(zt[keep])) * 100
}

# ---- main sweep ------------------------------------------------------------
alpha_grid <- seq(0,1,by=0.1)
n_moments_set <- c(4,6,8)
rows <- list()

for (dn in names(dists)) {
  st <- standardize(dists[[dn]])
  zt <- true_quantiles(st, levels_p)
  zg <- zgrid(st)

  # ----- Pearson (4 moments only) -----
  m3 <- z_moments_true(st, function(z) matrix(z^3, ncol=1))[1]
  m4 <- z_moments_true(st, function(z) matrix(z^4, ncol=1))[1]
  pear_eps <- NA
  pear <- tryCatch({
    pp <- pearsonFitM(mean=0, variance=1, skewness=m3, kurtosis=m4)
    zest <- qpearson(levels_p, params=pp)
    eps_metric(zt, zest)
  }, error=function(e) NA)
  rows[[length(rows)+1]] <- data.frame(dist=dn, method="Pearson", n=4,
                                       eps=pear, cond=NA, note="")

  for (nm in n_moments_set) {
    # ----- monomial-MaxEnt (GOPoly-accuracy proxy) -----
    bf_mono <- function(z) monomial_basis(z, nm)
    tm <- z_moments_true(st, bf_mono)
    fit <- solve_maxent(bf_mono(zg), tm, zg[2]-zg[1])
    if (!is.null(fit)) {
      zest <- q_from_fit(zg, fit, bf_mono, levels_p)
      rows[[length(rows)+1]] <- data.frame(dist=dn, method="Monomial-MaxEnt(GOPoly)",
        n=nm, eps=eps_metric(zt, zest), cond=fit$cond,
        note=ifelse(fit$converged,"","NOT-CONVERGED"))
    } else {
      rows[[length(rows)+1]] <- data.frame(dist=dn, method="Monomial-MaxEnt(GOPoly)",
        n=nm, eps=NA, cond=NA, note="FAILED")
    }

    # ----- Ku-MaxEnt (PATP, alpha-scan + dual-potential selection) -----
    best <- NULL; best_pot <- Inf; oracle_eps <- Inf; sel_eps <- NA; sel_cond <- NA
    for (al in alpha_grid) {
      bf <- (function(a){ function(z) patp_basis(z, nm, a) })(al)
      tmp <- z_moments_true(st, bf)
      f2 <- solve_maxent(bf(zg), tmp, zg[2]-zg[1])
      if (is.null(f2) || !f2$converged) next
      pot <- log(f2$Z) - sum(f2$lambdas*tmp)
      zest <- q_from_fit(zg, f2, bf, levels_p)
      e <- eps_metric(zt, zest)
      if (is.finite(e)) oracle_eps <- min(oracle_eps, e)
      if (pot < best_pot) { best_pot <- pot; sel_eps <- e; sel_cond <- f2$cond }
    }
    rows[[length(rows)+1]] <- data.frame(dist=dn, method="Ku-MaxEnt(PATP,sel)",
      n=nm, eps=sel_eps, cond=sel_cond, note="")
    rows[[length(rows)+1]] <- data.frame(dist=dn, method="Ku-MaxEnt(PATP,oracle)",
      n=nm, eps=ifelse(is.finite(oracle_eps),oracle_eps,NA), cond=NA, note="")

    # ----- Ku-MaxEnt (T-MaxEnt, correct element for multimodal); 2S=nm constraints
    S <- nm %/% 2
    p_grid <- c(0.2,0.3,0.4,0.5,0.6,0.8,1.0)
    t_best_pot <- Inf; t_sel_eps <- NA; t_sel_cond <- NA; t_oracle <- Inf
    for (pv in p_grid) {
      bf <- (function(pp){ function(z) trig_basis(z, S, pp) })(pv)
      tmp <- z_moments_true(st, bf)
      f3 <- solve_maxent(bf(zg), tmp, zg[2]-zg[1])
      if (is.null(f3) || !f3$converged) next
      pot <- log(f3$Z) - sum(f3$lambdas*tmp)
      zest <- q_from_fit(zg, f3, bf, levels_p)
      e <- eps_metric(zt, zest)
      if (is.finite(e)) t_oracle <- min(t_oracle, e)
      if (pot < t_best_pot) { t_best_pot <- pot; t_sel_eps <- e; t_sel_cond <- f3$cond }
    }
    rows[[length(rows)+1]] <- data.frame(dist=dn, method="Ku-MaxEnt(T-MaxEnt,sel)",
      n=nm, eps=t_sel_eps, cond=t_sel_cond, note="")
    rows[[length(rows)+1]] <- data.frame(dist=dn, method="Ku-MaxEnt(T-MaxEnt,oracle)",
      n=nm, eps=ifelse(is.finite(t_oracle),t_oracle,NA), cond=NA, note="")
  }
  cat(sprintf("done: %s (mu=%.3f, s=%.3f)\n", dn, st$mu, st$s))
}

res <- do.call(rbind, rows)
write.csv(res, "outputs/head_to_head/benchmark_results.csv", row.names=FALSE)

cat("\n================ RESULTS: expanded-uncertainty error eps (%) ================\n")
res$eps_r <- ifelse(is.na(res$eps), NA, round(res$eps,2))
print(res[,c("dist","method","n","eps_r","note")], row.names=FALSE)

cat("\n---- Mean eps across the 6 distributions, by method x n ----\n")
agg <- aggregate(eps ~ method + n, data=res[is.finite(res$eps),], FUN=mean)
agg$eps <- round(agg$eps,2)
print(agg[order(agg$n, agg$eps),], row.names=FALSE)

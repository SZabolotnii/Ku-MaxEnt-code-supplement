#!/usr/bin/env Rscript

# =========================================================================
# Stronger M&V baseline: randomized quasi-Monte-Carlo (RQMC) fitness
# Reviewer Q5 / Weakness W7.
# =========================================================================
# The paper compares the analytical MaxEnt M&V fitness against plain Monte
# Carlo (N=5000). The reviewer asks whether a stronger, variance-reduced /
# quasi-MC baseline closes the gap, plus per-seed final-design distributions.
#
# We add a randomized QMC estimator of q_0.05(Y): a 3-D Halton sequence
# (bases 2,3,5) with a Cranley-Patterson random shift, inverse-CDF mapped to
# the three input sample-means (W2,W3 exact normal; W1 the CLT-normal of the
# Beta(2,5) sample mean). RQMC is lower-variance than MC at equal N but still
# stochastic (the shift) and still costs N evaluations -- unlike the exactly
# deterministic, faster analytical evaluator.
#
# Outputs: (A) per-call noise + speed probe (MC vs RQMC vs MaxEnt);
#          (B) 3-arm GA with PER-SEED final designs (n1,n2,n3), cost, violation.
# Core functions copied verbatim from ablation_ga_noise.R. Base R only.
# Env: MV_GA_SEEDS (default 10) controls the GA arm; set 0 to run the probe only.
# =========================================================================

raw_moments_to_cumulants <- function(mu_prime) {
  r_max <- length(mu_prime); kappa <- rep(0, r_max)
  for (r in 1:r_max) { val <- mu_prime[r]
    if (r > 1) for (j in 1:(r-1)) val <- val - choose(r-1, j-1) * kappa[j] * mu_prime[r-j]
    kappa[r] <- val }
  kappa
}
cumulants_to_raw_moments <- function(kappa) {
  r_max <- length(kappa); mu_prime <- rep(0, r_max)
  for (r in 1:r_max) { val <- 0
    for (j in 1:r) { prev_mu <- if (r - j == 0) 1 else mu_prime[r - j]
      val <- val + choose(r-1, j-1) * kappa[j] * prev_mu }
    mu_prime[r] <- val }
  mu_prime
}
sample_mean_raw_moments <- function(pop_raw_moments, n) {
  pc <- raw_moments_to_cumulants(pop_raw_moments); sc <- pc
  for (r in seq_along(pc)) sc[r] <- pc[r] / (n^(r-1))
  cumulants_to_raw_moments(sc)
}
product_raw_moments <- function(n1, n2, n3, a = 1000, r_max = 4) {
  pop_raw_W1 <- rep(0, r_max); val <- 1
  for (r in 1:r_max) { val <- val * (2 + r - 1) / (7 + r - 1); pop_raw_W1[r] <- val }
  sW1 <- sample_mean_raw_moments(pop_raw_W1, n1)
  cW2 <- rep(0, r_max); cW2[1] <- 100; cW2[2] <- 100 / n2; sW2 <- cumulants_to_raw_moments(cW2)
  cW3 <- rep(0, r_max); cW3[1] <- 0.2; cW3[2] <- 0.0025 / n3; sW3 <- cumulants_to_raw_moments(cW3)
  raw_V <- sW1 * sW2 * sW3; raw_Y <- raw_V
  for (r in 1:r_max) raw_Y[r] <- (a^r) * raw_V[r]
  raw_Y
}
solve_maxent_raw <- function(raw_moments) {
  m1<-raw_moments[1]; m2<-raw_moments[2]; m3<-raw_moments[3]; m4<-raw_moments[4]
  sd_val <- sqrt(m2 - m1^2)
  mt3 <- (m3 - 3*m2*m1 + 2*m1^3)/(sd_val^3); mt4 <- (m4 - 4*m3*m1 + 6*m2*m1^2 - 3*m1^4)/(sd_val^4)
  target <- c(0, 1, mt3, mt4)
  xg <- seq(-6, 6, length.out = 1000); dx <- xg[2]-xg[1]
  Phi <- sapply(1:4, function(i) xg^i)
  lam <- rep(0, 4); conv <- FALSE
  for (it in 1:50) {
    un <- exp(as.numeric(Phi %*% lam)); if (any(!is.finite(un))) break
    Z <- sum(un)*dx; if (Z == 0 || is.na(Z)) break
    pdf <- un/Z; fit <- as.numeric(t(Phi) %*% pdf * dx); grad <- fit - target
    if (sqrt(sum(grad^2)) < 1e-5) { conv <- TRUE; break }
    H <- t(Phi) %*% (Phi*pdf)*dx - fit %*% t(fit); diag(H) <- diag(H) + 1e-7
    step <- tryCatch(solve(H, grad), error=function(e) tryCatch(qr.solve(H,grad), error=function(e2) NULL))
    if (is.null(step)) break
    a_s <- 1
    for (ls in 1:5) { nl <- lam - a_s*step; nu <- exp(as.numeric(Phi %*% nl))
      if (all(is.finite(nu))) { lam <- nl; break }; a_s <- a_s*0.5 }
  }
  if (!conv) return(NULL)
  pdf <- exp(as.numeric(Phi %*% lam))/Z; cdf <- cumsum(pdf)*dx
  m1 + xg[which.min(abs(cdf - 0.05))] * sd_val
}
mc_quantile_05 <- function(n1, n2, n3, a = 1000, N_mc = 10000) {
  W1 <- replicate(N_mc, mean(rbeta(n1, 2, 5)))
  W2 <- rnorm(N_mc, 100, 10/sqrt(n2)); W3 <- rnorm(N_mc, 0.2, 0.05/sqrt(n3))
  quantile(a * W1 * W2 * W3, 0.05, names = FALSE)
}

# ---- RQMC estimator: precomputed Halton + Cranley-Patterson shift ----
halton1 <- function(N, base) {
  r <- numeric(N)
  for (i in 1:N) { f <- 1; x <- 0; k <- i
    while (k > 0) { f <- f/base; x <- x + f*(k %% base); k <- k %/% base }
    r[i] <- x }
  r
}
N_QMC <- 5000
HALTON <- cbind(halton1(N_QMC, 2), halton1(N_QMC, 3), halton1(N_QMC, 5))  # computed once
varBeta <- (2*5)/((7^2)*(7+1))                                            # Var Beta(2,5) = 10/392
qmc_quantile_05 <- function(n1, n2, n3, a = 1000) {
  U <- sweep(HALTON, 2, runif(3), `+`) %% 1            # randomized (Cranley-Patterson) shift
  W1 <- qnorm(U[,1], 2/7, sqrt(varBeta/n1))            # CLT-normal of Beta sample mean
  W2 <- qnorm(U[,2], 100, 10/sqrt(n2)); W3 <- qnorm(U[,3], 0.2, 0.05/sqrt(n3))
  quantile(a * W1 * W2 * W3, 0.05, names = FALSE)
}

cost_function <- function(n) n[1] + 2.5*n[2] + 5*n[3]
evaluate_fitness <- function(n, method, N_mc = 5000) {
  n1<-round(n[1]); n2<-round(n[2]); n3<-round(n[3])
  if (n1<10||n1>500||n2<10||n2>500||n3<10||n3>500) return(-1e6)
  cost <- cost_function(c(n1,n2,n3))
  q05 <- if (method == "MaxEnt") { q <- solve_maxent_raw(product_raw_moments(n1,n2,n3)); if (is.null(q)) return(-1e6); q
         } else if (method == "QMC") qmc_quantile_05(n1,n2,n3)
         else mc_quantile_05(n1,n2,n3, N_mc = N_mc)
  if (q05 < 4800) -(cost + 1000*(4800 - q05)) else -cost
}
run_ga <- function(method, N_mc = 5000, pop = 30, gens = 20) {
  popm <- matrix(round(runif(pop*3, 10, 150)), pop, 3)
  for (g in 1:gens) {
    fit <- apply(popm, 1, evaluate_fitness, method = method, N_mc = N_mc)
    np <- popm
    for (i in 1:pop) { a<-sample(pop,1); b<-sample(pop,1); np[i,] <- popm[if (fit[a]>fit[b]) a else b, ] }
    popm <- np
    for (i in seq(1, pop-1, by=2)) if (runif(1)<0.8) { m<-runif(3)<0.5
      c1<-ifelse(m,popm[i,],popm[i+1,]); c2<-ifelse(m,popm[i+1,],popm[i,]); popm[i,]<-c1; popm[i+1,]<-c2 }
    for (i in 1:pop) if (runif(1)<0.2) { gn<-sample(1:3,1)
      popm[i,gn] <- max(10, min(500, popm[i,gn] + round(rnorm(1,0,15)))) }
  }
  fit <- apply(popm, 1, evaluate_fitness, method = method, N_mc = N_mc)
  sol <- round(popm[which.max(fit), ])
  tq <- mc_quantile_05(sol[1], sol[2], sol[3], N_mc = 200000)
  list(sol = sol, cost = cost_function(sol), true_q05 = tq, violation = tq < 4800)
}

# ---------------- (A) noise + speed probe ----------------
cat("\n=============== (A) M&V fitness probe at n=(100,50,30) ===============\n")
nc <- c(100, 50, 30); K <- 50
set.seed(20260612); t0<-Sys.time(); mc <- replicate(K, mc_quantile_05(nc[1],nc[2],nc[3],N_mc=5000)); tmc<-as.numeric(difftime(Sys.time(),t0,units="secs"))/K
set.seed(20260612); t0<-Sys.time(); qm <- replicate(K, qmc_quantile_05(nc[1],nc[2],nc[3]));            tqm<-as.numeric(difftime(Sys.time(),t0,units="secs"))/K
t0<-Sys.time(); me <- solve_maxent_raw(product_raw_moments(nc[1],nc[2],nc[3])); tme<-as.numeric(difftime(Sys.time(),t0,units="secs"))
cat(sprintf("                 mean q05    sd q05    <4800   per-call(s)   speedup vs MaxEnt\n"))
cat(sprintf("MC   (N=5000)  %9.2f %9.3f   %3d/%d   %.6f      %5.1fx\n", mean(mc), sd(mc), sum(mc<4800), K, tmc, tmc/tme))
cat(sprintf("RQMC (N=5000)  %9.2f %9.3f   %3d/%d   %.6f      %5.1fx\n", mean(qm), sd(qm), sum(qm<4800), K, tqm, tqm/tme))
cat(sprintf("MaxEnt (exact) %9.2f %9.3f   %3s     %.6f      %5s\n", me, 0, "-", tme, "1x"))
cat(sprintf("\nRQMC variance reduction vs MC: sd %.3f -> %.3f (%.1fx smaller)\n", sd(mc), sd(qm), sd(mc)/sd(qm)))

# ---------------- (B) 3-arm GA with per-seed designs ----------------
NS <- as.integer(Sys.getenv("MV_GA_SEEDS", "10"))
if (NS > 0) {
  cat(sprintf("\n=============== (B) 3-arm GA, %d matched seeds (per-seed designs) ===============\n", NS))
  arms <- c("MC", "QMC", "MaxEnt")
  store <- list()
  for (m in arms) {
    rows <- data.frame()
    for (r in 1:NS) {
      set.seed(1000 + r); t0 <- Sys.time(); res <- run_ga(m)
      rows <- rbind(rows, data.frame(seed=1000+r, n1=res$sol[1], n2=res$sol[2], n3=res$sol[3],
                    cost=res$cost, true_q05=res$true_q05, viol=res$violation,
                    time=as.numeric(difftime(Sys.time(),t0,units="secs"))))
    }
    store[[m]] <- rows
  }
  for (m in arms) {
    cat(sprintf("\n--- %s arm: per-seed final designs ---\n", m)); r <- store[[m]]
    for (i in 1:nrow(r)) cat(sprintf("  seed %d: (n1,n2,n3)=(%3d,%3d,%3d) cost=%.1f trueQ05=%.1f viol=%s\n",
        r$seed[i], r$n1[i], r$n2[i], r$n3[i], r$cost[i], r$true_q05[i], r$viol[i]))
    cat(sprintf("  cost: mean %.1f sd %.1f | violations %d/%d | mean time %.2fs\n",
        mean(r$cost), sd(r$cost), sum(r$viol), nrow(r), mean(r$time)))
  }
  cat("\n--- summary: violations and cost by arm ---\n")
  for (m in arms) { r <- store[[m]]
    cat(sprintf("  %-7s violations %d/%d | cost %.1f +/- %.1f | design sd (n1,n2,n3)=(%.0f,%.0f,%.0f)\n",
        m, sum(r$viol), nrow(r), mean(r$cost), sd(r$cost), sd(r$n1), sd(r$n2), sd(r$n3))) }
}
cat("\nDone.\n")

#!/usr/bin/env Rscript

# =========================================================================
# Automated T-MaxEnt frequency selection (reviewer Q6)
# =========================================================================
# The paper picks the trigonometric fundamental p by a hand-set rule
# |psi_hat_N(Sp)| >= 0.05, noting it is necessary but not sufficient. The
# reviewer asks for an automated rule that accounts for the O(1/sqrt(N)) ECF
# noise and reports its performance. We replace 0.05 with:
#
#   (i)  ADMISSIBILITY (N-aware noise floor): a config (p,S) is admissible iff
#        every harmonic clears the ECF noise floor by a 3-sigma margin,
#        |psi_hat_N(j p)| >= 3 / sqrt(N) for j = 1..S.  (At N=1000 this is
#        ~0.095, generalizing the fixed 0.05 to be sample-size aware.)
#   (ii) SELECTION among admissible configs by a held-out predictive log-score
#        (CV) -- the same variance/prediction logic as oPMMalpha. Amplitude
#        admissibility alone is necessary-not-sufficient (a coarse p=1.0 can
#        clear the floor yet fit poorly); the CV score rejects those.
#
# Reported: per-config diagnostics, and how often the automated rule selects a
# near-oracle config vs the fixed-0.05 default, over 10 seeds. Mixture bench.
# Solver/basis copied verbatim. Base R only.
# =========================================================================

trig_basis <- function(x, S, p) {
  M <- matrix(0, nrow = length(x), ncol = 2*S)
  for (r in 1:S) { M[, 2*r-1] <- cos(r*p*x); M[, 2*r] <- sin(r*p*x) }
  M
}
solve_maxent <- function(Phi, target_moments, dx, max_iter = 100, tol = 1e-6) {
  n <- ncol(Phi); lambdas <- rep(0, n)
  for (iter in 1:max_iter) {
    un <- exp(as.numeric(Phi %*% lambdas)); if (any(!is.finite(un))) return(NULL)
    Z <- sum(un)*dx; if (Z == 0 || is.na(Z)) return(NULL)
    pdf <- un/Z; fit <- as.numeric(t(Phi) %*% pdf * dx); grad <- fit - target_moments
    if (sqrt(sum(grad^2)) < tol) return(list(lambdas = lambdas, Z = Z, converged = TRUE))
    H <- t(Phi) %*% (Phi*pdf)*dx - fit %*% t(fit); diag(H) <- diag(H) + 1e-8
    step <- tryCatch(solve(H, grad), error=function(e) tryCatch(qr.solve(H,grad), error=function(e2) NULL))
    if (is.null(step)) return(NULL)
    a_s<-1; ok<-FALSE
    for (ls in 1:10) { nl<-lambdas-a_s*step; nu<-exp(as.numeric(Phi %*% nl))
      if (all(is.finite(nu))) { nZ<-sum(nu)*dx
        if (nZ>0 && !is.na(nZ) && log(nZ)-sum(nl*target_moments) <= log(Z)-sum(lambdas*target_moments)+1e-4) {
          lambdas<-nl; ok<-TRUE; break } }
      a_s<-a_s*0.5 }
    if (!ok) lambdas <- lambdas - 0.1*step
  }
  list(converged = FALSE)
}

N <- 1000; sd_true <- sqrt(0.32)
xg <- seq(-5, 5, length.out = 1000); dx <- xg[2]-xg[1]
ptrue <- 0.4*dnorm(xg,-1,sd_true) + 0.6*dnorm(xg,1,sd_true)
floor_3sig <- 3 / sqrt(N)                       # 3-sigma ECF noise floor
configs <- expand.grid(p = c(0.3,0.5,0.7,1.0), S = c(2,3,4,5))

eval_config <- function(d, p, S) {
  half <- seq_len(length(d)) %% 2 == 0
  Pg <- trig_basis(xg, S, p); tg <- colMeans(trig_basis(d, S, p))
  f <- solve_maxent(Pg, tg, dx)
  if (is.null(f) || !f$converged) return(NULL)
  pdf <- exp(as.numeric(Pg %*% f$lambdas))/f$Z
  mse <- mean((pdf - ptrue)^2)
  # ECF modulus at each harmonic; binding = top harmonic Sp
  ecf_mod <- sapply(1:S, function(j) Mod(mean(exp(1i * j*p * d))))
  admissible <- all(ecf_mod >= floor_3sig)
  # held-out predictive log-score (fit on train half, NLL on holdout half)
  ft <- solve_maxent(Pg, colMeans(trig_basis(d[!half], S, p)), dx)
  nll <- NA
  if (!is.null(ft) && ft$converged) {
    lf <- as.numeric(trig_basis(d[half], S, p) %*% ft$lambdas) - log(ft$Z)
    nll <- -mean(lf[is.finite(lf)])
  }
  data.frame(p=p, S=S, top_ecf=ecf_mod[S], admissible=admissible, NLL=nll, mse=mse)
}

cat(sprintf("\n=========== Automated T-MaxEnt frequency rule (mixture) ===========\n"))
cat(sprintf("ECF noise floor 1/sqrt(N) = %.4f ; 3-sigma admissibility threshold = %.4f (replaces fixed 0.05)\n\n", 1/sqrt(N), floor_3sig))

# single-seed config table
set.seed(20260612); u <- runif(N); d0 <- ifelse(u<0.4, rnorm(N,-1,sd_true), rnorm(N,1,sd_true))
tab <- do.call(rbind, lapply(1:nrow(configs), function(i) eval_config(d0, configs$p[i], configs$S[i])))
tab <- tab[order(tab$p, tab$S), ]
cat(sprintf("%4s %3s %9s %6s %9s %11s\n", "p","S","top|ECF|","adm.","CV-NLL","PDF MSE"))
for (i in 1:nrow(tab)) with(tab[i,],
  cat(sprintf("%4.1f %3d %9.3f %6s %9.4f %11.2e\n", p, S, top_ecf, ifelse(admissible,"yes","NO"), NLL, mse)))
adm <- tab[tab$admissible & is.finite(tab$NLL), ]
auto <- adm[which.min(adm$NLL), ]
oracle <- tab[which.min(tab$mse), ]
cat(sprintf("\nAutomated rule (3-sigma admissible + best CV log-score) picks (p=%.1f, S=%d): PDF MSE %.2e\n",
            auto$p, auto$S, auto$mse))
cat(sprintf("Oracle-best config                         (p=%.1f, S=%d): PDF MSE %.2e\n", oracle$p, oracle$S, oracle$mse))
cat(sprintf("Fixed-0.05 paper default                   (p=0.5, S=3): PDF MSE %.2e\n",
            tab$mse[tab$p==0.5 & tab$S==3]))

# 10-seed: how often does the automated rule land on a near-oracle config?
cat(sprintf("\n--- 10-seed replication ---\n"))
pick_auto <- pick_def <- c(); pen_auto <- pen_def <- c()
for (s in 1:10) {
  set.seed(s); u <- runif(N); d <- ifelse(u<0.4, rnorm(N,-1,sd_true), rnorm(N,1,sd_true))
  T <- do.call(rbind, lapply(1:nrow(configs), function(i) eval_config(d, configs$p[i], configs$S[i])))
  best <- min(T$mse)
  adm <- T[T$admissible & is.finite(T$NLL), ]
  a <- adm[which.min(adm$NLL), ]
  pick_auto <- c(pick_auto, sprintf("(%.1f,%d)", a$p, a$S)); pen_auto <- c(pen_auto, a$mse/best)
  defmse <- T$mse[T$p==0.5 & T$S==3]; pen_def <- c(pen_def, defmse/best)
}
cat(sprintf("Automated rule: mean PDF-MSE penalty vs oracle = %.2fx ; picks = %s\n",
            mean(pen_auto), paste(table(pick_auto), names(table(pick_auto)), collapse=", ")))
cat(sprintf("Fixed 0.05 default (p=0.5,S=3): mean penalty vs oracle = %.2fx\n", mean(pen_def)))
cat("\nDone.\n")

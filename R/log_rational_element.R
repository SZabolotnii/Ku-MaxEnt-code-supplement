#!/usr/bin/env Rscript

# Testing logarithmic / rational generating elements for MaxEnt reconstruction
# of heavy-tailed laws, against the paper's power-PATP and T-MaxEnt bases.
#
# Hypothesis: a LOG element log(1+(x/s)^2) gives f ∝ (1+(x/s)^2)^λ — the
# Student/Cauchy family — so it CAN represent algebraic (power-law) tails;
# a BOUNDED rational element 1/(1+(x/s)^2) cannot (its exponent → 0 on the
# tail, flattening it).  Solver copied VERBATIM from patp_maxent_simulation.R.

solve_maxent <- function(Phi, target_moments, dx, max_iter = 200, tol = 1e-6) {
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
  pdf <- unnorm_pdf / Z; fitted <- as.numeric(t(Phi) %*% pdf * dx)
  H <- t(Phi) %*% (Phi * pdf) * dx - fitted %*% t(fitted); diag(H) <- diag(H) + 1e-8
  list(lambdas = lambdas, Z = Z, cond_num = kappa(H, exact = TRUE), converged = FALSE)
}

# ----- basis constructors -----
patp_power <- function(i, alpha) 1/i + (4 - i - 3/i)*alpha + (2*i - 4 + 2/i)*alpha^2
basis_patp0 <- function(x) { # power-PATP alpha=0, parity-matched, 4 constraints
  M <- cbind(x, abs(x)^patp_power(2,0), sign(x)*abs(x)^patp_power(3,0), abs(x)^patp_power(4,0)); M }
basis_trig <- function(x, S=2, p=0.2) { M<-matrix(0,length(x),2*S); for(r in 1:S){M[,2*r-1]<-cos(r*p*x);M[,2*r]<-sin(r*p*x)}; M }
basis_log1 <- function(x) cbind(log(1+x^2))                                   # 1 constraint
basis_logmulti <- function(x, s=c(0.5,1,2,4)) sapply(s, function(ss) log(1+(x/ss)^2))   # 4
basis_hybrid <- function(x) cbind(x^2, log(1+x^2))                            # Gaussian core + algebraic tail, 2
basis_brat <- function(x, s=c(0.5,1,2,4)) sapply(s, function(ss) 1/(1+(x/ss)^2))        # bounded rational, 4

BASES <- list(
  "power-PATP a=0" = basis_patp0,
  "T-MaxEnt(0.2,2)" = function(x) basis_trig(x),
  "LogRat-1"        = basis_log1,
  "LogRat-multi"    = basis_logmulti,
  "Hybrid x^2+log"  = basis_hybrid,
  "BoundedRat"      = basis_brat
)

# ----- true densities -----
stable_pdf <- function(x, alpha=1.5) { # symmetric alpha-stable via CF inversion: (1/pi)∫_0^inf exp(-u^a)cos(ux)du
  du <- 0.005; u <- seq(0, 40, by=du); w <- exp(-u^alpha)
  w[1] <- w[1]/2; w[length(w)] <- w[length(w)]/2   # TRAPEZOID endpoints (rectangle rule injects du/(2pi) DC offset)
  sapply(x, function(xx) sum(w * cos(u*xx)) * du / pi)
}
stable_cdf_grid <- function(xg, pdf) cumsum(pdf)*(xg[2]-xg[1])

TARGETS <- list(
  Cauchy   = list(gen=function(N) rcauchy(N),       pdf=function(x) dcauchy(x),       cdf=function(x) pcauchy(x),       L=50, true_slope=-2),
  "t(df=2)"= list(gen=function(N) rt(N,2),          pdf=function(x) dt(x,2),          cdf=function(x) pt(x,2),          L=50, true_slope=-3),
  "SaS a=1.5"=list(gen=NULL, pdf=function(x) stable_pdf(x,1.5), cdf=NULL,             L=50, true_slope=-2.5, stable=TRUE),
  Gaussian = list(gen=function(N) rnorm(N),         pdf=function(x) dnorm(x),         cdf=function(x) pnorm(x),         L=8,  true_slope=NA)
)

tail_slope <- function(xg, pdf_fit, lo, hi) { # slope of log10 f vs log10|x| on [lo,hi] (x>0)
  idx <- which(xg >= lo & xg <= hi & pdf_fit > 1e-12)
  if (length(idx) < 5) return(NA)
  coef(lm(log10(pdf_fit[idx]) ~ log10(xg[idx])))[2]
}

set.seed(20260612); N <- 1000
# stable sample via Chambers-Mallows-Stuck (symmetric)
rstable_sym <- function(N, alpha) {
  U <- runif(N, -pi/2, pi/2); W <- rexp(N)
  sin(alpha*U)/(cos(U))^(1/alpha) * (cos(U - alpha*U)/W)^((1-alpha)/alpha)
}

cat("================================================================\n")
cat(" GENERATING-ELEMENT COMPARISON (N=1000, seed 20260612)\n")
cat(" metrics: conv | kappa_H | KS | bodyMSE | Q95 err% | Q99 err% | tail-slope (true)\n")
cat("================================================================\n")

results <- data.frame()
store <- list()
for (tname in names(TARGETS)) {
  tg <- TARGETS[[tname]]
  if (!is.null(tg$stable) && tg$stable) { set.seed(20260612); xs <- rstable_sym(N, 1.5) }
  else { set.seed(20260612); xs <- tg$gen(N) }
  L <- tg$L; xg <- seq(-L, L, length.out = 4000); dx <- xg[2]-xg[1]
  pdf_true <- tg$pdf(xg)
  if (!is.null(tg$cdf)) cdf_true <- tg$cdf(xg) else cdf_true <- stable_cdf_grid(xg, pdf_true)
  body_lo <- if (tname=="Gaussian") -4 else -5; body_hi <- -body_lo
  bidx <- which(xg>=body_lo & xg<=body_hi)
  kidx <- which(xg>=-0.9*L & xg<=0.9*L)
  q95t <- xg[which.min(abs(cdf_true-0.95))]; q99t <- xg[which.min(abs(cdf_true-0.99))]
  ts_lo <- if (tname=="Gaussian") 2 else 8; ts_hi <- 0.8*L
  cat(sprintf("\n--- %s  (true tail slope %s; true Q95=%.2f Q99=%.2f) ---\n",
              tname, ifelse(is.na(tg$true_slope),"n/a",sprintf("%.1f",tg$true_slope)), q95t, q99t))
  for (bname in names(BASES)) {
    Phi_d <- BASES[[bname]](xs); tgt <- colMeans(as.matrix(Phi_d))
    Phi_g <- as.matrix(BASES[[bname]](xg))
    fit <- solve_maxent(Phi_g, tgt, dx)
    if (is.null(fit)) { cat(sprintf("  %-16s | solver-NULL\n", bname)); next }
    pf <- exp(as.numeric(Phi_g %*% fit$lambdas))/fit$Z
    ks <- max(abs(cumsum(pf)*dx - cdf_true)[kidx])
    bmse <- mean((pf[bidx]-pdf_true[bidx])^2)
    cdf_f <- cumsum(pf)*dx
    q95 <- xg[which.min(abs(cdf_f-0.95))]; q99 <- xg[which.min(abs(cdf_f-0.99))]
    e95 <- (q95-q95t)/q95t*100; e99 <- (q99-q99t)/q99t*100
    sl <- tail_slope(xg, pf, ts_lo, ts_hi)
    results <- rbind(results, data.frame(target=tname, basis=bname, conv=fit$converged,
                     kH=fit$cond_num, ks=ks, bmse=bmse, e95=e95, e99=e99, slope=sl,
                     true_slope=tg$true_slope))
    store[[paste(tname,bname)]] <- pf
    cat(sprintf("  %-16s | %-5s | %.2e | %.4f | %.3e | %+7.1f | %+7.1f | %s\n",
                bname, fit$converged, fit$cond_num, ks, bmse, e95, e99,
                ifelse(is.na(sl),"n/a",sprintf("%.2f",sl))))
  }
}

# ----- summary: KS and tail-slope recovery for heavy targets -----
cat("\n================================================================\n")
cat(" SUMMARY: KS (lower=better) by basis x target\n")
cat("================================================================\n")
hv <- results[results$target!="Gaussian",]
tab <- tapply(results$ks, list(results$basis, results$target), function(z) z[1])
print(round(tab,4))
cat("\n LogRat-1 fitted lambda (Cauchy should be ~ -1):\n")
# refit lograt-1 on Cauchy to print lambda
set.seed(20260612); xc<-rcauchy(N); xg<-seq(-50,50,length.out=4000); dx<-xg[2]-xg[1]
fl<-solve_maxent(as.matrix(basis_log1(xg)), colMeans(as.matrix(basis_log1(xc))), dx)
cat(sprintf("   lambda_hat = %.4f  (=> tail slope 2*lambda = %.2f, true Cauchy -2)\n", fl$lambdas[1], 2*fl$lambdas[1]))

# ----- figure: Cauchy log-density, generating elements compared -----
OUT_DIR <- Sys.getenv("GENELEMENT_OUT_DIR", unset = "outputs")
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)
xg<-seq(-50,50,length.out=4000)
pdf(file.path(OUT_DIR, "fig_genelement_cauchy.pdf"), width=7, height=4.5)
par(mar=c(4.5,4.5,1,1))
sel<-which(xg>=-30 & xg<=30)
plot(xg[sel], log10(dcauchy(xg[sel])), type="l", lwd=2.5, col="black",
     xlab="x", ylab=expression(log[10]~f(x)), ylim=c(-5,0.2))
cols<-c("LogRat-1"="#0072B2","power-PATP a=0"="#D55E00","T-MaxEnt(0.2,2)"="#009E73","BoundedRat"="#CC79A7")
ltys<-c("LogRat-1"=1,"power-PATP a=0"=2,"T-MaxEnt(0.2,2)"=4,"BoundedRat"=3)
for (b in names(cols)) { k<-paste("Cauchy",b); if(!is.null(store[[k]])) lines(xg[sel], log10(pmax(store[[k]][sel],1e-12)), lwd=2, lty=ltys[b], col=cols[b]) }
legend("topright", bty="n", lwd=2, lty=c(1,ltys), col=c("black",cols),
       legend=c("True Cauchy", names(cols)))
dev.off()
cat("\nFigure written: figures/fig_genelement_cauchy.pdf\n\nDone.\n")

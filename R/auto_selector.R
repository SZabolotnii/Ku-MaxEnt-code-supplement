#!/usr/bin/env Rscript
# ============================================================================
# Automatic element selector (closes the open problem flagged in Exp 5).
# Two-stage: (1) route by Hartigan's dip test of multimodality --- significant
# dip -> trigonometric element (bounded/multimodal); else -> logarithmic element
# (a strong general unimodal element, Gaussian through Cauchy). (2) within-family
# hyperparameters by held-out log-score (T-MaxEnt frequency) or fixed scales
# (LogRat). Each element is paired with its natural support treatment: T-MaxEnt
# on a tight standardized grid, LogRat on a wide grid (heavy tails).
#
# Claim to test: the AUTO-routed fit matches the per-case ORACLE (best of the two
# elements) and never collapses to the wrong element the way a single held-out
# log-score selector does. Metric: tail-quantile error at 90/95/99% vs truth.
# ============================================================================
suppressWarnings(suppressMessages({library(diptest); library(moments)}))

solve_maxent <- function(Phi, tm, dx, max_iter=250, tol=1e-7){ Phi<-as.matrix(Phi); lam<-rep(0,ncol(Phi))
  for(it in 1:max_iter){ u<-exp(as.numeric(Phi%*%lam)); if(any(!is.finite(u))) return(NULL)
    Z<-sum(u)*dx; if(!is.finite(Z)||Z<=0) return(NULL); pdf<-u/Z; fm<-as.numeric(t(Phi)%*%pdf*dx); g<-fm-tm
    H<-t(Phi)%*%(Phi*pdf)*dx-fm%*%t(fm); diag(H)<-diag(H)+1e-9
    if(sqrt(sum(g^2))<tol) return(list(lam=lam,Z=Z,conv=TRUE))
    st<-tryCatch(solve(H,g),error=function(e) tryCatch(qr.solve(H,g),error=function(e2) NULL)); if(is.null(st)) return(NULL)
    a<-1; ok<-FALSE; for(ls in 1:20){ nl<-lam-a*st; nu<-exp(as.numeric(Phi%*%nl))
      if(all(is.finite(nu))){ nZ<-sum(nu)*dx; if(is.finite(nZ)&&nZ>0){ if(log(nZ)-sum(nl*tm)<=log(Z)-sum(lam*tm)+1e-5){ lam<-nl; ok<-TRUE; break } } }; a<-a*0.5 }
    if(!ok) lam<-lam-0.05*st }
  u<-exp(as.numeric(Phi%*%lam)); Z<-sum(u)*dx; if(!is.finite(Z)||Z<=0) return(NULL)
  pdf<-u/Z; fm<-as.numeric(t(Phi)%*%pdf*dx); list(lam=lam,Z=Z,conv=sqrt(sum((fm-tm)^2))<1e-4) }
trig_bf <- function(p,S=3) function(x){ M<-matrix(0,length(x),2*S); for(r in 1:S){M[,2*r-1]<-cos(r*p*x);M[,2*r]<-sin(r*p*x)}; M }
logm_bf <- function(x,s=c(0.5,1,2,4)) sapply(s,function(ss) log(1+(x/ss)^2))

qlev<-c(.90,.95,.99)
# robust standardization (median / MAD): finite and stable for both infinite-variance
# heavy tails (Cauchy MAD=1) and small-scale real data (returns MAD~0.01).
rstd <- function(x){ m<-median(x); sc<-mad(x); if(!is.finite(sc)||sc==0) sc<-sd(x); list(m=m,sc=sc,z=(x-m)/sc) }
# --- fit one element on data x, return quantile estimates at qlev (original units) ---
fit_lograt <- function(x){ r<-rstd(x); z<-r$z; L<-50; zg<-seq(-L,L,length.out=4000); dz<-zg[2]-zg[1]
  zw<-z[abs(z)<=L]; f<-solve_maxent(logm_bf(zg),colMeans(logm_bf(zw)),dz); if(is.null(f)||!f$conv) return(NULL)
  pf<-exp(as.numeric(logm_bf(zg)%*%f$lam))/f$Z; cdf<-cumsum(pf)*dz; cdf<-cdf/cdf[length(cdf)]
  r$m + r$sc*approx(cdf,zg,xout=qlev,ties="ordered",rule=2)$y }
fit_trig <- function(x){ r<-rstd(x); mu<-r$m; s<-r$sc; z<-r$z; L<-max(abs(z))+1; zg<-seq(-L,L,length.out=3000); dz<-zg[2]-zg[1]
  ntr<-floor(.7*length(z)); tr<-z[1:ntr]; te<-z[(ntr+1):length(z)]
  best<-NULL; bsc<--Inf
  for(p in c(0.2,0.3,0.4,0.5,0.6,0.8)){ bf<-trig_bf(p); f<-solve_maxent(bf(zg),colMeans(bf(tr)),dz)
    if(is.null(f)||!f$conv) next; sc<-mean((as.numeric(bf(te)%*%f$lam)-log(f$Z))); if(is.finite(sc)&&sc>bsc){bsc<-sc;best<-list(f=f,p=p)} }
  if(is.null(best)) return(NULL); bf<-trig_bf(best$p); f<-solve_maxent(bf(zg),colMeans(bf(z)),dz); if(is.null(f)||!f$conv) return(NULL)
  pf<-exp(as.numeric(bf(zg)%*%f$lam))/f$Z; cdf<-cumsum(pf)*dz; cdf<-cdf/cdf[length(cdf)]
  qz<-approx(cdf,zg,xout=qlev,ties="ordered",rule=2)$y; mu+s*qz }
qerr <- function(qest,qtrue) if(is.null(qest)) NA else mean(abs(qest-qtrue)/abs(qtrue))*100

# targets with true quantiles
sas_q<-function(p){du<-.002;u<-seq(0,60,du);w<-exp(-u^1.5);w[1]<-w[1]/2;w[length(w)]<-w[length(w)]/2
  xr<-seq(-200,200,length.out=8000);pd<-sapply(xr,function(xx) sum(w*cos(u*xx))*du/pi);pd[pd<0]<-0
  cd<-cumsum(pd)*(xr[2]-xr[1]);cd<-cd/cd[length(cd)];approx(cd,xr,xout=p,ties="ordered",rule=2)$y}
targets<-list(
  Cauchy=list(r=function(N)rcauchy(N),q=qcauchy(qlev)),
  t2=list(r=function(N)rt(N,2),q=qt(qlev,2)),
  stable1.5=list(r=function(N){U<-runif(N,-pi/2,pi/2);W<-rexp(N);a<-1.5;sin(a*U)/cos(U)^(1/a)*(cos(U-a*U)/W)^((1-a)/a)},q=sas_q(qlev)),
  bimodal_close=list(r=function(N){u<-runif(N);ifelse(u<.4,rnorm(N,-1,.4),rnorm(N,1,.4))},
    q={xr<-seq(-6,6,length.out=20000);pd<-.4*dnorm(xr,-1,.4)+.6*dnorm(xr,1,.4);cd<-cumsum(pd)*(xr[2]-xr[1]);cd<-cd/cd[length(cd)];approx(cd,xr,xout=qlev,ties="ordered",rule=2)$y}),
  Gaussian=list(r=function(N)rnorm(N),q=qnorm(qlev)))

cat("Auto-selector vs fixed elements vs oracle: tail-quantile error (%) at 90/95/99, mean/route over 8 seeds\n\n")
cat(sprintf("%-14s %8s %8s %8s %8s   %s\n","target","AUTO","fix-Log","fix-Trig","ORACLE","route(of 8)"))
for(tn in names(targets)){
  eA<-eL<-eT<-eO<-c(); rt<-c()
  for(s in 1:8){ set.seed(20260800+s); x<-targets[[tn]]$r(1500); qt<-targets[[tn]]$q
    route<- if(dip.test(x)$p.value<0.05) "T" else "Log"; rt<-c(rt,route)
    qL<-fit_lograt(x); qT<-fit_trig(x)
    qA<-if(route=="T") qT else qL
    eA<-c(eA,qerr(qA,qt)); eL<-c(eL,qerr(qL,qt)); eT<-c(eT,qerr(qT,qt))
    both<-c(qerr(qL,qt),qerr(qT,qt)); both<-both[is.finite(both)]; eO<-c(eO, if(length(both)) min(both) else NA) }
  tab<-table(rt)
  cat(sprintf("%-14s %8.1f %8.1f %8.1f %8.1f   %s\n",tn,mean(eA,na.rm=T),mean(eL,na.rm=T),mean(eT,na.rm=T),mean(eO,na.rm=T),
    paste(sprintf("%s:%d",names(tab),as.integer(tab)),collapse=" ")))
}
# real data
R<-diff(log(EuStockMarkets))
cat("\n-- real EuStockMarkets (validate vs empirical quantiles) --\n")
for(j in 1:4){ x<-as.numeric(R[,j]); qe<-as.numeric(quantile(x,qlev,type=7))
  route<- if(dip.test(x)$p.value<0.05) "T" else "Log"; qL<-fit_lograt(x); qT<-fit_trig(x)
  qA<-if(route=="T") qT else qL
  cat(sprintf("%-14s AUTO=%.1f (route %s)  fix-Log=%.1f  fix-Trig=%.1f\n",
    colnames(R)[j],qerr(qA,qe),route,qerr(qL,qe),qerr(qT,qe))) }

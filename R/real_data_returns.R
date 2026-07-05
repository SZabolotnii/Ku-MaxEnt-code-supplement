#!/usr/bin/env Rscript
# ============================================================================
# Real-data expanded-uncertainty demonstration.
# Data: EuStockMarkets (base R) daily log-returns of 4 European indices
#       (DAX, SMI, CAC, FTSE) -- real measured, moderately heavy-tailed series
#       (excess kurtosis 2.4-6.3), finite variance, near-symmetric.
# Task: reconstruct the return distribution from a few moments and estimate the
#       tail quantiles (the coverage bounds of the expanded uncertainty), then
#       validate against the EMPIRICAL quantiles of the full sample (n=1859).
# Methods: Ku-MaxEnt matched elements (LogRat-multi, PM-PATP) vs Pearson vs
#          monomial-MaxEnt (GOPoly-accuracy proxy). Standardised (mean/sd).
# Metric: mean |q_est - q_emp| / |q_emp| over the 1/5/10/90/95/99% levels.
# ============================================================================
suppressWarnings(suppressMessages({library(PearsonDS); library(moments)}))

solve_maxent <- function(Phi, tm, dx, max_iter=250, tol=1e-7){
  Phi<-as.matrix(Phi); lambdas<-rep(0,ncol(Phi))
  for(iter in 1:max_iter){ u<-exp(as.numeric(Phi%*%lambdas)); if(any(!is.finite(u))) return(NULL)
    Z<-sum(u)*dx; if(!is.finite(Z)||Z<=0) return(NULL)
    pdf<-u/Z; fm<-as.numeric(t(Phi)%*%pdf*dx); grad<-fm-tm
    H<-t(Phi)%*%(Phi*pdf)*dx-fm%*%t(fm); diag(H)<-diag(H)+1e-9
    if(sqrt(sum(grad^2))<tol) return(list(lambdas=lambdas,Z=Z,converged=TRUE,cond=kappa(H,exact=TRUE)))
    step<-tryCatch(solve(H,grad),error=function(e) tryCatch(qr.solve(H,grad),error=function(e2) NULL)); if(is.null(step)) return(NULL)
    a<-1; ok<-FALSE
    for(ls in 1:20){ nl<-lambdas-a*step; nu<-exp(as.numeric(Phi%*%nl))
      if(all(is.finite(nu))){ nZ<-sum(nu)*dx; if(is.finite(nZ)&&nZ>0){
        if(log(nZ)-sum(nl*tm)<=log(Z)-sum(lambdas*tm)+1e-5){ lambdas<-nl; ok<-TRUE; break } } }; a<-a*0.5 }
    if(!ok) lambdas<-lambdas-0.05*step }
  u<-exp(as.numeric(Phi%*%lambdas)); Z<-sum(u)*dx; if(!is.finite(Z)||Z<=0) return(NULL)
  pdf<-u/Z; fm<-as.numeric(t(Phi)%*%pdf*dx); H<-t(Phi)%*%(Phi*pdf)*dx-fm%*%t(fm); diag(H)<-diag(H)+1e-9
  list(lambdas=lambdas,Z=Z,converged=sqrt(sum((fm-tm)^2))<1e-4,cond=kappa(H,exact=TRUE)) }
patp_power <- function(i,alpha) 1/i + (4-i-3/i)*alpha + (2*i-4+2/i)*alpha^2
basis_logmulti <- function(x,s=c(0.5,1,2,4)) sapply(s,function(ss) log(1+(x/ss)^2))
basis_pmpatp   <- function(x) cbind(x, abs(x)^patp_power(2,0), sign(x)*abs(x)^patp_power(3,0), abs(x)^patp_power(4,0))
monomial_basis <- function(x,n){ M<-matrix(0,length(x),n); for(i in 1:n) M[,i]<-x^i; M }

lev <- c(0.01,0.05,0.10,0.90,0.95,0.99)
R <- diff(log(EuStockMarkets)); series <- colnames(R)
L<-15; zg<-seq(-L,L,length.out=4000); dx<-zg[2]-zg[1]
qfit_err <- function(fit,bf,z,qemp){ pf<-exp(as.numeric(as.matrix(bf(zg))%*%fit$lambdas))/fit$Z
  cdf<-cumsum(pf)*dx; cdf<-cdf/cdf[length(cdf)]; qe<-approx(cdf,zg,xout=lev,ties="ordered",rule=2)$y
  mean(abs(qe-qemp)/abs(qemp))*100 }

rows<-list()
for(sname in series){
  x0<-as.numeric(R[,sname]); mu<-mean(x0); s<-sd(x0); z<-(x0-mu)/s
  qemp<-as.numeric(quantile(z, lev, type=7)); exk<-round(kurtosis(z)-3,1)
  # Pearson
  pe<-tryCatch({ pp<-pearsonFitM(moments=c(mean=0,variance=1,skewness=skewness(z),kurtosis=kurtosis(z)))
    qe<-qpearson(lev,params=pp); mean(abs(qe-qemp)/abs(qemp))*100 }, error=function(e) NA)
  rows[[length(rows)+1]]<-data.frame(series=sname,exkurt=exk,method="Pearson",q_err=round(pe,1))
  # monomial n=4
  fm<-solve_maxent(monomial_basis(zg,4), colMeans(monomial_basis(z,4)), dx)
  rows[[length(rows)+1]]<-data.frame(series=sname,exkurt=exk,method="Monomial-MaxEnt(GOPoly)",
    q_err=if(!is.null(fm)&&fm$converged) round(qfit_err(fm,function(z)monomial_basis(z,4),z,qemp),1) else NA)
  # Ku PM-PATP
  fp<-solve_maxent(basis_pmpatp(zg), colMeans(basis_pmpatp(z)), dx)
  rows[[length(rows)+1]]<-data.frame(series=sname,exkurt=exk,method="Ku-MaxEnt(PM-PATP)",
    q_err=if(!is.null(fp)&&fp$converged) round(qfit_err(fp,basis_pmpatp,z,qemp),1) else NA)
  # Ku LogRat
  fl<-solve_maxent(basis_logmulti(zg), colMeans(basis_logmulti(z)), dx)
  rows[[length(rows)+1]]<-data.frame(series=sname,exkurt=exk,method="Ku-MaxEnt(LogRat)",
    q_err=if(!is.null(fl)&&fl$converged) round(qfit_err(fl,basis_logmulti,z,qemp),1) else NA)
}
res<-do.call(rbind,rows)
write.csv(res,"outputs/head_to_head/real_data_returns_results.csv",row.names=FALSE)
cat("\n===== REAL DATA: EuStockMarkets log-returns =====\n")
cat("tail-quantile error vs empirical, mean over 1/5/10/90/95/99% (%), lower=better\n\n")
ord<-c("Ku-MaxEnt(LogRat)","Ku-MaxEnt(PM-PATP)","Pearson","Monomial-MaxEnt(GOPoly)")
w<-reshape(res, idvar=c("series","exkurt"), timevar="method", direction="wide")
colnames(w)<-sub("q_err.","",colnames(w))
print(w[,c("series","exkurt",ord)], row.names=FALSE)
cat("\nMean across the 4 series:\n")
ag<-aggregate(q_err~method, data=res, FUN=function(v) round(mean(v,na.rm=TRUE),1))
print(ag[match(ord,ag$method),], row.names=FALSE)

#!/usr/bin/env Rscript
# ============================================================================
# Heavy-tailed head-to-head (CORRECTED): Ku-MaxEnt matched elements
# (logarithmic LogRat-multi + parity-matched PATP) vs Pearson system vs
# monomial-MaxEnt (GOPoly-accuracy proxy). Multi-seed, empirical moments.
#
# Correction over v1: v1 used the raw all-ODD sign-preserving PATP basis, which
# by parity cannot represent a symmetric decaying density (it piles mass at the
# truncation edge -> garbage quantiles). The matched Ku element for symmetric
# algebraic tails is the LOGARITHMIC one, log(1+(x/s)^2), giving f ~ (1+(x/s)^2)^L
# (the Student/Cauchy family); the parity-matched PATP mixes even & odd terms.
#
# Metric: KS on the modelling window + tail-quantile errors at 90/95/99% (within
# window). Reported as mean over seeds with convergence/feasibility rate.
# ============================================================================
suppressWarnings(suppressMessages({library(PearsonDS); library(moments)}))

solve_maxent <- function(Phi, tm, dx, max_iter=250, tol=1e-7){
  Phi<-as.matrix(Phi); lambdas<-rep(0,ncol(Phi))
  for(iter in 1:max_iter){
    u<-exp(as.numeric(Phi%*%lambdas)); if(any(!is.finite(u))) return(NULL)
    Z<-sum(u)*dx; if(!is.finite(Z)||Z<=0) return(NULL)
    pdf<-u/Z; fm<-as.numeric(t(Phi)%*%pdf*dx); grad<-fm-tm
    H<-t(Phi)%*%(Phi*pdf)*dx-fm%*%t(fm); diag(H)<-diag(H)+1e-9
    if(sqrt(sum(grad^2))<tol) return(list(lambdas=lambdas,Z=Z,converged=TRUE,cond=kappa(H,exact=TRUE)))
    step<-tryCatch(solve(H,grad),error=function(e) tryCatch(qr.solve(H,grad),error=function(e2) NULL))
    if(is.null(step)) return(NULL)
    a<-1; ok<-FALSE
    for(ls in 1:20){ nl<-lambdas-a*step; nu<-exp(as.numeric(Phi%*%nl))
      if(all(is.finite(nu))){ nZ<-sum(nu)*dx; if(is.finite(nZ)&&nZ>0){
        if(log(nZ)-sum(nl*tm) <= log(Z)-sum(lambdas*tm)+1e-5){ lambdas<-nl; ok<-TRUE; break } } }
      a<-a*0.5 }
    if(!ok) lambdas<-lambdas-0.05*step
  }
  u<-exp(as.numeric(Phi%*%lambdas)); Z<-sum(u)*dx; if(!is.finite(Z)||Z<=0) return(NULL)
  pdf<-u/Z; fm<-as.numeric(t(Phi)%*%pdf*dx); H<-t(Phi)%*%(Phi*pdf)*dx-fm%*%t(fm); diag(H)<-diag(H)+1e-9
  list(lambdas=lambdas,Z=Z,converged=sqrt(sum((fm-tm)^2))<1e-4,cond=kappa(H,exact=TRUE))
}
patp_power <- function(i,alpha) 1/i + (4-i-3/i)*alpha + (2*i-4+2/i)*alpha^2
# Ku elements
basis_logmulti <- function(x, s=c(0.5,1,2,4)) sapply(s, function(ss) log(1+(x/ss)^2))     # matched heavy-tail
basis_pmpatp   <- function(x) cbind(x, abs(x)^patp_power(2,0), sign(x)*abs(x)^patp_power(3,0), abs(x)^patp_power(4,0)) # parity-matched
monomial_basis <- function(x,n){ M<-matrix(0,length(x),n); for(i in 1:n) M[,i]<-x^i; M }   # GOPoly proxy

targets <- list(
  Cauchy     = list(r=function(N) rcauchy(N),        q=function(p) qcauchy(p)),
  Student_t3 = list(r=function(N) rt(N,3),           q=function(p) qt(p,3)),
  Student_t2 = list(r=function(N) rt(N,2),           q=function(p) qt(p,2)),
  Stable1p5  = list(r=function(N){ U<-runif(N,-pi/2,pi/2); W<-rexp(N); a<-1.5
                     sin(a*U)/(cos(U))^(1/a)*(cos(U-a*U)/W)^((1-a)/a) },
                    q=NULL)  # symmetric stable; true quantiles from fine reference below
)
# reference quantiles for SaS(1.5) via CF inversion cdf on a fine grid
sas_q <- function(p){ du<-0.002; u<-seq(0,60,du); w<-exp(-u^1.5); w[1]<-w[1]/2; w[length(w)]<-w[length(w)]/2
  xr<-seq(-200,200,length.out=8000); pdf<-sapply(xr,function(xx) sum(w*cos(u*xx))*du/pi)
  pdf[pdf<0]<-0; cdf<-cumsum(pdf)*(xr[2]-xr[1]); cdf<-cdf/cdf[length(cdf)]
  approx(cdf,xr,xout=p,ties="ordered",rule=2)$y }

qlev <- c(0.90,0.95,0.99); plev<-c(qlev, 1-qlev)
L<-50; xg<-seq(-L,L,length.out=4000); dx<-xg[2]-xg[1]
qfit <- function(fit,bf){ pf<-exp(as.numeric(as.matrix(bf(xg))%*%fit$lambdas))/fit$Z
  cdf<-cumsum(pf)*dx; cdf<-cdf/cdf[length(cdf)]; list(cdf=cdf,pf=pf) }
qerr <- function(cdf, qtrue){ qe<-approx(cdf,xg,xout=qlev,ties="ordered",rule=2)$y
  mean(abs(qe-qtrue)/abs(qtrue))*100 }
ks_body <- function(cdf, cdftrue){ idx<-which(xg>=-0.9*L&xg<=0.9*L); max(abs(cdf-cdftrue)[idx]) }

N<-1000; seeds<-1:10; rows<-list()
for(tn in names(targets)){
  qtrue <- if(is.null(targets[[tn]]$q)) sas_q(qlev) else targets[[tn]]$q(qlev)
  # true CDF on grid
  if(tn=="Stable1p5"){ du<-0.002; u<-seq(0,60,du); w<-exp(-u^1.5); w[1]<-w[1]/2; w[length(w)]<-w[length(w)]/2
    pdft<-sapply(xg,function(xx) sum(w*cos(u*xx))*du/pi); pdft[pdft<0]<-0; cdft<-cumsum(pdft)*dx; cdft<-cdft/cdft[length(cdft)]
  } else if(tn=="Cauchy"){ cdft<-pcauchy(xg) } else { df<-as.integer(sub("Student_t","",tn)); cdft<-pt(xg,df) }
  for(s in seeds){
    set.seed(20260700+s); x<-targets[[tn]]$r(N); xw<-x[abs(x)<=L]
    # Pearson
    pe<-tryCatch({ pp<-pearsonFitM(moments=c(mean=mean(x),variance=var(x),skewness=skewness(x),kurtosis=kurtosis(x)))
      qe<-qpearson(qlev,params=pp); mean(abs(qe-qtrue)/abs(qtrue))*100 }, error=function(e) NA)
    rows[[length(rows)+1]]<-data.frame(target=tn,seed=s,method="Pearson",q_err=pe,ks=NA,conv=!is.na(pe))
    # monomial-MaxEnt n=4 (GOPoly proxy)
    fm<-solve_maxent(monomial_basis(xg,4), colMeans(monomial_basis(xw,4)), dx)
    if(!is.null(fm)&&fm$converged){ o<-qfit(fm,function(z)monomial_basis(z,4))
      rows[[length(rows)+1]]<-data.frame(target=tn,seed=s,method="Monomial-MaxEnt(GOPoly)",q_err=qerr(o$cdf,qtrue),ks=ks_body(o$cdf,cdft),conv=TRUE)
    } else rows[[length(rows)+1]]<-data.frame(target=tn,seed=s,method="Monomial-MaxEnt(GOPoly)",q_err=NA,ks=NA,conv=FALSE)
    # Ku parity-matched PATP
    fp<-solve_maxent(basis_pmpatp(xg), colMeans(basis_pmpatp(xw)), dx)
    if(!is.null(fp)&&fp$converged){ o<-qfit(fp,basis_pmpatp)
      rows[[length(rows)+1]]<-data.frame(target=tn,seed=s,method="Ku-MaxEnt(PM-PATP)",q_err=qerr(o$cdf,qtrue),ks=ks_body(o$cdf,cdft),conv=TRUE)
    } else rows[[length(rows)+1]]<-data.frame(target=tn,seed=s,method="Ku-MaxEnt(PM-PATP)",q_err=NA,ks=NA,conv=FALSE)
    # Ku LogRat-multi (matched heavy-tail element)
    fl<-solve_maxent(basis_logmulti(xg), colMeans(basis_logmulti(xw)), dx)
    if(!is.null(fl)&&fl$converged){ o<-qfit(fl,basis_logmulti)
      rows[[length(rows)+1]]<-data.frame(target=tn,seed=s,method="Ku-MaxEnt(LogRat)",q_err=qerr(o$cdf,qtrue),ks=ks_body(o$cdf,cdft),conv=TRUE)
    } else rows[[length(rows)+1]]<-data.frame(target=tn,seed=s,method="Ku-MaxEnt(LogRat)",q_err=NA,ks=NA,conv=FALSE)
  }
  cat("done:",tn,"\n")
}
res<-do.call(rbind,rows)
write.csv(res,"outputs/head_to_head/heavy_tailed_v2_results.csv",row.names=FALSE)
cat("\n===== HEAVY-TAILED (corrected): mean over 10 seeds =====\n")
cat("q_err = mean |q_est-q_true|/|q_true| at 90/95/99% (%); ks = body KS; conv = feasible-rate\n\n")
for(tn in names(targets)){ cat("--",tn,"--\n"); sub<-res[res$target==tn,]
  a<-aggregate(cbind(q_err,ks,conv)~method, data=sub, FUN=function(v) mean(v,na.rm=TRUE), na.action=na.pass)
  a$q_err<-round(a$q_err,1); a$ks<-round(a$ks,4); a$conv<-round(a$conv,2)
  ord<-c("Ku-MaxEnt(LogRat)","Ku-MaxEnt(PM-PATP)","Pearson","Monomial-MaxEnt(GOPoly)")
  print(a[match(ord,a$method),],row.names=FALSE); cat("\n") }

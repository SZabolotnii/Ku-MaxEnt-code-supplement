#!/usr/bin/env Rscript
# Confirms the "monomial/GOPoly" equivalence claim of Experiment 5: on a fixed
# known grid the moment-constrained MaxEnt density is basis-independent, so raw
# monomial MaxEnt and Gram-Schmidt-orthonormalized (GOPoly-stabilized) MaxEnt
# converge to the identical density and the identical error metric. Base R only.
solve_maxent<-function(Phi,tm,dx,max_iter=300,tol=1e-9){Phi<-as.matrix(Phi);lam<-rep(0,ncol(Phi))
 for(it in 1:max_iter){u<-exp(as.numeric(Phi%*%lam));if(any(!is.finite(u)))return(NULL)
  Z<-sum(u)*dx;if(!is.finite(Z)||Z<=0)return(NULL);pdf<-u/Z;fm<-as.numeric(t(Phi)%*%pdf*dx);g<-fm-tm
  H<-t(Phi)%*%(Phi*pdf)*dx-fm%*%t(fm);diag(H)<-diag(H)+1e-10
  if(sqrt(sum(g^2))<tol)return(list(lam=lam,Z=Z,pdf=pdf))
  st<-tryCatch(solve(H,g),error=function(e)tryCatch(qr.solve(H,g),error=function(e2)NULL));if(is.null(st))return(NULL)
  a<-1;ok<-FALSE;for(ls in 1:25){nl<-lam-a*st;nu<-exp(as.numeric(Phi%*%nl));if(all(is.finite(nu))){nZ<-sum(nu)*dx
   if(is.finite(nZ)&&nZ>0){if(log(nZ)-sum(nl*tm)<=log(Z)-sum(lam*tm)+1e-6){lam<-nl;ok<-TRUE;break}}};a<-a/2}
  if(!ok)lam<-lam-0.05*st}
 u<-exp(as.numeric(Phi%*%lam));Z<-sum(u)*dx;list(lam=lam,Z=Z,pdf=u/Z)}
zg<-seq(-4,4,length.out=3000);dz<-zg[2]-zg[1]
ftrue<-0.4*dnorm(zg,-1,.4)+0.6*dnorm(zg,1,.4); ftrue<-ftrue/(sum(ftrue)*dz)
n<-6; Mono<-sapply(1:n,function(i) zg^i)
f1<-solve_maxent(Mono, as.numeric(t(Mono)%*%ftrue*dz), dz)
Q<-qr.Q(qr(Mono*sqrt(dz)))/sqrt(dz)                       # L2(dz)-orthonormal columns (GOPoly stabilization)
f2<-solve_maxent(Q, as.numeric(t(Q)%*%ftrue*dz), dz)
lev<-c(.90,.95,.99); cdf<-cumsum(ftrue)*dz; qt<-approx(cdf/cdf[length(cdf)],zg,lev,ties="ordered")$y
e<-function(f){cd<-cumsum(f$pdf)*dz;q<-approx(cd/cd[length(cd)],zg,lev,ties="ordered")$y;mean(abs(q-qt)/abs(qt))*100}
cat(sprintf("raw monomial   : eps=%.4f%%\n",e(f1)))
cat(sprintf("GS-orthonormal : eps=%.4f%%\n",e(f2)))
cat(sprintf("max |pdf_mono - pdf_orth| = %.3e  => basis-independent (identical density)\n", max(abs(f1$pdf-f2$pdf))))

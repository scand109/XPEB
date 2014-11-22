logistic=function(x) 1/(1+exp(-x))
ilogistic=function(y) -log(1/y-1)




makeCluster=function(dummy, marker.info, ld.trim.bp=2e5, flank=1e4) {
  ### dummy=1 if the marker exceeded the threshold
  ### marker.info has format (snpid, chr, bp);
  bump.indx=which(dummy==1);
  
  sp=diff(marker.info[bump.indx,3]);
  sp.chr=diff(marker.info[bump.indx,2]);
  loci.cut=c(0, which(sp>=ld.trim.bp| sp.chr!=0), length(bump.indx));

  window.mat=matrix(ncol=6,nrow=length(loci.cut)-1);
  rm.indx=list();
  for (j in 1:(length(loci.cut)-1)) {
    a=bump.indx[loci.cut[j]+1]; b=bump.indx[loci.cut[j+1]];
    left=min(which(marker.info[,2]==marker.info[a,2] & abs(marker.info[a,3]-marker.info[,3])<flank))
    right=max(which(marker.info[,2]==marker.info[b,2] & abs(marker.info[b,3]-marker.info[,3])<flank))
    ##window.mat[j,]=c(bump.indx[loci.cut[j]+1], bump.indx[loci.cut[j+1]], marker.info[bump.indx[loci.cut[j]+1],2],marker.info[bump.indx[loci.cut[j+1]],2], marker.info[bump.indx[loci.cut[j]+1],3], marker.info[bump.indx[loci.cut[j+1]],3])
    window.mat[j,]=c(left, right, marker.info[left,2], marker.info[right,2], marker.info[left,3], marker.info[right,3])
    
    rm.indx[[j]]=window.mat[j,1]: window.mat[j,2];  #bump.indx[loci.cut[j]+1]:bump.indx[loci.cut[j+1]]
    }
  return(list(window.mat=window.mat, rm.indx=unlist(rm.indx)))
}



grid.count=function(chi1.base, grid.endpoints) {
  grid.length=length(grid.endpoints)-1
  tabulate(as.integer(cut(chi1.base,grid.endpoints)),nbins=grid.length)
}




runif.int=function(n,max.int,min.int=1,replace=TRUE) { ## uniform on min.int..max.int
  x=sample.int(max.int-min.int+1,n,replace=replace)
  return(x+min.int-1)
}



ddirichlet=function(x,alpha,log.p=FALSE) {
  if (log.p == FALSE) return(exp(ddirichlet(x,alpha,log.p=TRUE)))
  else {
    lgamma(sum(alpha))-sum(lgamma(alpha))+sum((alpha-1)*log(x))
  }
}


## density for x in (0,1] for mean parameter 0<tau<=1/2 (or equiv. 0<=alpha<1); expected value tau=(1-alpha)/(2-alpha)
d.zero.one.pow=function(x,tau,alpha=tau.to.alpha(tau),log.p=FALSE) {
  if (log.p==FALSE) (1-alpha)*x^(-alpha)
  else log(1-alpha)-alpha*log(x)
}
p.zero.one.pow=function(x,tau,alpha=tau.to.alpha(tau)) {
  x^(1-alpha)
}
q.zero.one.pow=function(p,tau,alpha=tau.to.alpha(tau)) {
  p^(1/(1-alpha))
}
tau.to.alpha=function(tau) (1-2*tau)/(1-tau) ## for tau in (0,1/2] the mean parameter, convert to alpha scale



## These two functions are for reallocating mass p1 and p2; p1,p2>=0, p1+p2<=1
## the relative mass is converted to a Cauchy inverse-cdf scale where it can be adjusted neatly
## then mapped back

qD=qnorm
pD=pnorm
dD=dnorm

#pre-vectorization
p1p2.to.q=function(p1,p2) { ## use the more numerically stable way to compute
  if (p1<p2) qD(p1/(p1+p2))
  else -qD(p2/(p1+p2))
}
vectorized.p1p2.to.q=function(p1,p2) { ## use the more numerically stable way to compute
  v=qD(p1/(p1+p2))
  v2=-qD(p2/(p1+p2))
  v[p2>p1]=v2[p2>p1]
  return(v)
}
get.p1p2=function(q,psum) {
  c(pD(q)*psum, pD(-q)*psum)
}




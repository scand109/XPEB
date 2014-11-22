p11.barrier.scale=.01
#p11.constraint=.9  #now set in main script



logit2=function(x) {
	(1/(2+exp(x)+exp(-x)))	 ### deriviative
}


negllik=function(beta, y, density.delta.target, density.null.target) {
	p11=logistic(beta[1]+beta[2]); ## P(T=1 | B=1)
        ## include a barrier to enforce p11<p11.constraint; barrier will have derivative 1 at p11.constraint-p11.barrier.scale, rapidly increasing to the right and decreasing to the left
        if (p11>=p11.constraint) return(Inf)
        barrier=p11.barrier.scale^2/(p11.constraint-p11)
	p10=logistic(beta[1]); ## P(T=1| B=0)
	v=p11*density.delta.target*y + p10*(1-y)*density.delta.target+density.null.target	
	return( - sum(log(v)) + barrier)
}

dnegllik=function(beta, y, density.delta.target, density.null.target) {
	p11=logistic(beta[1]+beta[2]); ## P(T=1 | B=1)
	p10=logistic(beta[1]); ## P(T=1| B=0)
        dbarrierdp11=((p11.constraint-p11)/p11.barrier.scale)^(-2)
        if (p11>=p11.constraint) dbarrierdp11=0 ## this condition should not happen with the bfgs_optim method, but could with the optim one, I suppose.
        ddty=density.delta.target*y
        v=p11*ddty + p10*(1-y)*density.delta.target+density.null.target
        ##r=-sum(log(v))
        dp11db1=logit2(beta[1]+beta[2])
        dp11db2=dp11db1
        dp10db1=logit2(beta[1]);
        dvdb1=dp11db1*ddty + dp10db1*(1-y)*density.delta.target
        dvdb2=dp11db2*ddty
        drdb1=-sum(dvdb1/v)+dbarrierdp11*dp11db1
        drdb2=-sum(dvdb2/v)+dbarrierdp11*dp11db2
        
	return( c(drdb1,drdb2) )
}

eb.modelFitting=function( density.null.target, density.alt.target, post.alt.base) {
  ##beta0=c(-8, 6)
  beta0=c(-7, 7)
  #beta0=c(-8.3,9.3)
	
  density.delta.target=density.alt.target-density.null.target

  gamma.to.beta=function(gamma){c(gamma[1],exp(gamma[2]))}
  beta.to.gamma=function(beta){c(beta[1],log(beta[2]))}
  	
  opt=bfgs_options()
  opt$initialstep=.01
  opt$growstep=1.2
  opt$trace=100
  opt$f=function(gamma) negllik(gamma.to.beta(gamma),y=post.alt.base, density.delta.target=density.delta.target, density.null.target=density.null.target)
  opt$df=function(gamma) c(1,exp(gamma[2]))*dnegllik(gamma.to.beta(gamma),y=post.alt.base, density.delta.target=density.delta.target, density.null.target=density.null.target)
  state=bfgs_init(beta.to.gamma(beta0),opt)
  res=bfgs_run(state,opt)
  state=res$state
  beta1.1=gamma.to.beta(state$x)
  
  o2=optim(beta.to.gamma(beta0),opt$f,gr=opt$df,method='BFGS')
  beta1.2=gamma.to.beta(o2$par)

  o3=optim(beta.to.gamma(beta0),opt$f)
  beta1.3=gamma.to.beta(o3$par)

  beta1=beta1.1;method=1;value=state$f1
  if (o2$value < value) {beta1=beta1.2;method=2;value=o2$value}
  if (o3$value < value) {beta1=beta1.3;method=3;value=o3$value}
  
  
  
  p11=logistic(beta1[1]+beta1[2]); ## P(T=1 | B=1)
  p10=logistic(beta1[1]); ## P(T=1| B=0)

  tmp1=density.alt.target*p11*post.alt.base+density.alt.target*p10*(1-post.alt.base)
  tmp0=density.null.target*(1-p11)*post.alt.base+density.null.target*(1-p10)*(1-post.alt.base)
  post.prob.alt=tmp1/(tmp1+tmp0) # conditioned on chi1.base and chi1.target
  #print(sprintf("Max likelihood a: %g", -negllik(beta=beta1, y=post.alt.base, density.delta.target=density.delta.target, density.null.target=density.null.target)))
  #print(sprintf("max log likelihood =%g", sum(log(tmp1+tmp0))))
##  return(list(betas=beta1, post.prob.alt.target=post.prob.alt, method=method));
  return(list(betas=beta1, post.prob.alt.target=post.prob.alt, method=method, nl=function(beta) negllik(beta,y=post.alt.base, density.delta.target=density.delta.target, density.null.target=density.null.target)));
}

analyzeArun5=function(dat,n.base,n.target,tau=1e-6,grid.length=500,data.frac=1,alt.basis.size=25,plot.style='hua',diagnostics=FALSE,n.iter=1e6,block.scale=1,noise.scale=.2,block.scale.shape=c(.02,.1,.2),burn.in=5000,skip.initial.bases=0) {

  noise.scale.block=block.scale*block.scale.shape
  epsilon=1E-20

  if (min(n.base,n.target)<=100) stop('sample size too small')
  trans.base=make.transforms(transform.b=n.base,J=alt.basis.size,grid.length=grid.length+1)
  trans.target=make.transforms(transform.b=n.target,J=alt.basis.size,grid.length=grid.length+1)

  chi1.target=dat[,1]
  chi1.base=dat[,2]
  x.target=trans.target$psi(chi1.target)
  x.base=trans.base$psi(chi1.base)
  i.target=pmin(pmax(ceiling(x.target*grid.length),1),grid.length)
  i.base=pmin(pmax(ceiling(x.base*grid.length),1),grid.length)
  ct.target=tabulate(i.target,nbins=grid.length)
  ct.base=tabulate(i.base,nbins=grid.length)

  make.grid.basis=function(alt.basis) {
    min.val=0
    grid.endpoints=seq(0,1,length=grid.length+1)
    grid.basis=matrix(NA,grid.length,alt.basis$J+1)
    grid.basis[,1]=pmax(diff(alt.basis$G0(grid.endpoints)),min.val)
    for (j in 1:alt.basis$J) {
      grid.basis[,j+1]=pmax(diff(alt.basis$Gl[[j]](grid.endpoints)),min.val)
    }
    grid.basis=t(t(grid.basis)/colSums(grid.basis))
  }
  grid.basis.base=make.grid.basis(trans.base)
  grid.basis.target=make.grid.basis(trans.target)
  B.base=grid.basis.base[,-(1:(1+skip.initial.bases))]
  B.target=grid.basis.target[,-(1:(1+skip.initial.bases))]
  f0.base=grid.basis.base[,1]
  f0.target=grid.basis.target[,1]

  J=alt.basis.size-skip.initial.bases
  ## q.base, q.target, q.alpha1, alt.mix (of length J)
  theta0=c(qD(1E-5),qD(1E-5),0, rep(1,J)/J)

  lpost=function(theta) {
    qblock=theta[1:3]
    pblock=pD(qblock)
    alt.mix=theta[-(1:3)]
    l0=sum(dD(qblock,log=TRUE))
    alpha.base=pblock[1]
    alpha.target=pblock[2]
    alpha1=pblock[3]
    if (alpha1<=0 || alpha.base>.01 ||alpha.target>.01) return(-Inf)
    #if (alpha1<=0) return(-Inf)
    f1.base=B.base %*% alt.mix
    f1.target=B.target %*% alt.mix
    m.base=alpha.base*f1.base+(1-alpha.base)*f0.base
    m.target=alpha.target*f1.target+(1-alpha.target)*f0.target
    l1=sum(ct.base*log(m.base))+sum(ct.target*log(m.target))
    l2=d.zero.one.pow(alpha.base,tau,log.p=TRUE)+d.zero.one.pow(alpha.target,tau,log.p=TRUE)
    ##l2=dexp(alpha.base,rate=1/tau,log=TRUE)+dexp(alpha.target,rate=1/tau,log=TRUE)
    l3=ddirichlet(alt.mix,rep(alpha1,J),log.p=TRUE)
    l=l0+data.frac*l1+l2+l3
    if (is.na(l) || l==Inf) l=-Inf
    return(l)
  }

  mix.proposal=function(theta) {
    qblock=theta[1:3]
    pi=theta[-(1:3)]
    ww=runif.int(1,2)
    if (ww==1) {
      iijj=runif.int(2,min.int=1,max.int=length(pi),replace=FALSE)
      ii=iijj[1];jj=iijj[2]
    } else {
      ii=runif.int(1,min.int=1,max.int=length(pi)-1)
      jj=ii+1
    }
    
    L=epsilon;R=min(1,pi[ii]+pi[jj])
    
    psum=pi[ii]+pi[jj]
    q=p1p2.to.q(pi[ii],pi[jj])
    qprime=rcauchy(1,q,noise.scale)
    p1p2prime=get.p1p2(qprime,psum)
    lpropratio=dD(qprime,log=TRUE)-dD(q,log=TRUE)
    pi[ii]=p1p2prime[1];pi[jj]=p1p2prime[2]
    pi=pmax(pi,epsilon);pi=pi/sum(pi);
    
    return(list(theta=c(qblock,pi),pi=pi,lpropratio=lpropratio,ii=ii,jj=jj))
  }

  qblock.proposal=function(theta) {
    qblock=theta[1:3]
    pi=theta[-(1:3)]

    ii=runif.int(1,length(qblock))
    q=qblock[ii]
    qp=rcauchy(1,qblock,noise.scale.block[ii])
    qblock[ii]=qp
    lpropratio=0
    return(list(theta=c(qblock,pi),lpropratio=lpropratio,ii=ii))
  }

  if (diagnostics) theta.hist=matrix(NA,n.iter,length(theta0)) else theta.hist=NULL

  theta=theta0
  lprob1=lpost(theta)
  pblock.sum=rep(0,3)
  alt.mix.sum=rep(0,J)
  
  for (kk in 1:n.iter) {
    pick=runif.int(1,5)
    if (pick==1) {
      prop=qblock.proposal(theta)
    } else {
      prop=mix.proposal(theta)
    }
    lprobp=lpost(prop$theta)
    
    delt=lprobp-lprob1+prop$lpropratio
    step=FALSE
    if (!is.na(delt) && (delt>0 || runif(1)<exp(delt))) {
      theta=prop$theta
      lprob1=lprobp
      step=TRUE
    }
    if (diagnostics) {
      theta.hist[kk,]=theta
      if (FALSE && kk%%1000==1) {
        print(theta)
      }
    }
    if (kk>burn.in) {
      alt.mix.sum=alt.mix.sum+theta[-(1:3)]
      qblock=theta[1:3]
      pblock=pD(qblock)
      pblock.sum=pblock.sum+pblock
    }
  }

  n.eff=n.iter-burn.in
  pblock=pblock.sum/n.eff
  alt.mix=alt.mix.sum/n.eff
  alpha.base=pblock[1] ## pi_0
  alpha.target=pblock[2] ## pi^'_0
  alpha1=pblock[3] ## alpha
  
  f1.base=B.base %*% alt.mix
  f1.target=B.target %*% alt.mix
  m.base=alpha.base*f1.base+(1-alpha.base)*f0.base
  m.target=alpha.target*f1.target+(1-alpha.target)*f0.target

  PR.base=alpha.base*f1.base/m.base
  PR.target=alpha.target*f1.target/m.target
  post.alt.base=PR.base[i.base]
  post.alt.target=PR.target[i.target]

  density.null.target=f0.target[i.target]
  density.alt.target=f1.target[i.target]
  
  eb=eb.modelFitting( density.null.target, density.alt.target, post.alt.base)

  beta=eb$beta
  p1=logistic(beta[1]+beta[2]); ## P(T=1 | B=1)
  p0=logistic(beta[1]); ## P(T=1| B=0)

  #direct.posterior=function(f1.base,f0.base,i.base,alpha.base) {
  #  m.base=alpha.base*f1.base+(1-alpha.base)*f0.base
  #  PR.base=alpha.base*f1.base/m.base
  #  return(PR.base[i.base])
  #}
  
  
  ##print(sum(eb$post.prob.alt.target>.95))
  ##browser()
  
  return(list(theta.hist=theta.hist,pblock=pblock,alt.mix=alt.mix,alpha1=alpha1,alpha.base=alpha.base,alpha.target=alpha.target,
              B.base=B.base, B.target=B.target, f0.base=f0.base, f0.target=f0.target, f1.base=f1.base, f1.target=f1.target, m.base=m.base, m.target=m.target, i.base=i.base, i.target=i.target, ct.base=ct.base, ct.target=ct.target,
              eb=eb, post.alt.base=post.alt.base, post.alt.target=post.alt.target, trans.base=trans.base, trans.target=trans.target, n.base=n.base, n.target=n.target,p1=p1,p0=p0, PR.base=PR.base, PR.target=PR.target, grid.length=grid.length
              ))
}

eb.predict.postprob=function(dat, res) {
  chi1.target=dat[,1];
  chi1.base=dat[,2]
  
  betas=res$eb$betas
  trans.base=res$trans.base
  trans.target=res$trans.target
  grid.length=res$grid.length
  x.target=trans.target$psi(chi1.target)
  x.base=trans.base$psi(chi1.base)
  i.target=pmin(pmax(ceiling(x.target*grid.length),1),grid.length)
  i.base=pmin(pmax(ceiling(x.base*grid.length),1),grid.length)
  
  density.null.target=res$f0.target[i.target]
  density.alt.target=res$f1.target[i.target]

  ######
  alpha.base=res$alpha.base ## pi_0
  alpha.target=res$alpha.target ## pi^'_0

  f0.base=res$f0.base
  f1.base=res$f1.base
  m.base=alpha.base*f1.base+(1-alpha.base)*f0.base
  PR.base=alpha.base*f1.base/m.base
  
  post.alt.base=PR.base[i.base] ## recalculated PR.base incorporates changes to alpha.base

  p11=logistic(betas[1]+betas[2]); ## P(T=1 | B=1)
  p10=logistic(betas[1]); ## P(T=1| B=0)

  ## handle the complete data
  tmp1=density.alt.target*p11*post.alt.base+density.alt.target*p10*(1-post.alt.base)
  tmp0=density.null.target*(1-p11)*post.alt.base+density.null.target*(1-p10)*(1-post.alt.base)
  post.prob.alt=tmp1/(tmp1+tmp0) # conditioned on chi1.base and chi1.target

  ## target only
  ix=is.na(i.base) & !is.na(i.target)
  post.prob.alt[ix]=res$PR.target[i.target[ix]]

  ## base only
  ix=!is.na(i.base) & is.na(i.target)
  bb=res$PR.base[i.base[ix]] ## post.prob.alt.base[ix]
  post.prob.alt[ix]=bb*res$p1+(1-bb)*res$p0
  
  return(post.prob.alt)
}

plot.shared.deconvolve=function(dat,res,xseq=seq(0,40,.25),yseq=seq(0,40,.25)) {
  X=outer(xseq,yseq,function(x,y) x)
  Y=outer(xseq,yseq,function(x,y) y)
  post.prob.alt.grid=eb.predict.postprob(cbind(as.vector(X),as.vector(Y)), res)
  Ar=matrix(post.prob.alt.grid,nrow(X),ncol(X))
  plot(c(0,40),c(0,40),type='n',xlab='target chi-square stat', ylab='base chi-square stat')
  contour(xseq,yseq,Ar,levels=c(.5,.75,.9,.95,.99,.999),add=TRUE)
  chi1.target=dat[,1];
  chi1.base=dat[,2]
  xx1=pmin(chi1.target,40)
  xx1[xx1==40]=jitter(xx1[xx1==40])
  yy1=pmin(chi1.base,40)
  yy1[yy1==40]=jitter(yy1[yy1==40])
  points(xx1,yy1,pch='.',col='blue')
  post.prob.alt.target=eb.predict.postprob(dat, res)
  ix=post.prob.alt.target>.95
  points(xx1[ix],yy1[ix],pch=8,col='blue')
}




eb.predict.postprob.target.only=function(dat, res) {
  chi1.target=dat[,1];
  trans.target=res$trans.target
  grid.length=res$grid.length
  x.target=trans.target$psi(chi1.target)
  i.target=pmin(pmax(ceiling(x.target*grid.length),1),grid.length)
  post.prob.alt=res$PR.target[i.target]
  return(post.prob.alt)
}


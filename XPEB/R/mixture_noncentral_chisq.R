make.transforms=function(transform.b=1E5,transform.c=9,J=10,grid.length=500,subdivisions=1E2) {

ipsi=function(x) transform.b*(exp(transform.c*x^2)-1)/(exp(transform.c)-1)
dipsi=function(x) 2*transform.b*transform.c*x/(exp(transform.c)-1)*exp(transform.c*x^2)
psi=function(y) sqrt(log(y/transform.b*(exp(transform.c)-1)+1)/transform.c)
dpsi=function(y) {
  a=y/transform.b*(exp(transform.c)-1)+1
  b=log(a)/transform.c
  1/2/sqrt(b)*(exp(transform.c)-1)/transform.b/transform.c/a
}

#ipsi=function(x) transform.b*(exp(transform.c*x)-1)/(exp(transform.c)-1)
#dipsi=function(x) transform.b*transform.c/(exp(transform.c)-1)*exp(transform.c*x)
#psi=function(y) log(y/transform.b*(exp(transform.c)-1)+1)/transform.c
#dpsi=function(y) (exp(transform.c)-1)/transform.b/transform.c/(y/transform.b*(exp(transform.c)-1)+1)

g0=function(y) {
	dipsi(y)*dchisq(ipsi(y),1)
}

G0=function(y) pchisq(ipsi(y),1)
G1=function(y,theta) pchisq(ipsi(y),1,ncp=ipsi(theta))
G=Vectorize(function(y,alpha,beta){
	f=function(theta) G1(y,theta)*dbeta(theta,alpha,beta)
	integrate(f,0,1,subdivisions=subdivisions)$value
})

make.interpG=function(alpha,beta,G,val0=0,val1=1) {
	ygrid=seq(0,1,length=grid.length)[-c(1,grid.length)]
	val=G(ygrid,alpha,beta)
	yy=c(0,ygrid,1)
	vv=c(val0,val,val1)
	GG=approxfun(yy,vv)
}

make.interpg=function(GG,delta=1E-5) {
	ygrid=seq(0,1,length=grid.length)
	ygrid[grid.length]=1-delta
	ygrid[1]=delta
	valL=GG(ygrid-delta)
	valR=GG(ygrid+delta)
	vv=(valR-valL)/(2*delta)
	yy=seq(0,1,length=grid.length)
	gg=approxfun(yy,vv)
}



make.basis.interps=function(J=10) {
	Gl=list()
	gl=list()
	for (j in 1:J) {
		GG=make.interpG(j,J+1-j,G)
		gg=make.interpg(GG)
		Gl[[j]]=GG
		gl[[j]]=gg
	}
	return(list(Gl=Gl,gl=gl,G0=G0,g0=g0,J=J,psi=psi,ipsi=ipsi,dipsi=dipsi,dpsi=dpsi))
}




show.noncentrality.distribs=function(J=10,cdf=FALSE, ...) {
	v=exp(seq(log(1/10000),log(1E6),length=5001))
	##B=make.basis.interps(J)
	for (j in 1:J) {
		if (cdf) vv=pbeta(psi(v),j,J+1-j) else vv=dbeta(psi(v),j,J+1-j)*dpsi(v)
		if (j==1) plot(v,vv,type='l',col=j,lty=j,xlab='x',ylab=if (cdf) 'cdf of noncentrality distributions' else 'density of noncentrality distribution',...) else {
			lines(v,vv,col=j,lty=j)
			}
		}
	##abline(v=0,col='red')
}

#test.make.basis.interps()
#test.make.basis.interps(J=10,cdf=FALSE)
#show.mixed.noncentral.chisq.basis()
#show.mixed.noncentral.chisq.basis(cdf=TRUE)
#show.mixed.noncentral.chisq.basis(v=exp(seq(log(1/1000),log(1E6),length=5001)),log='xy',ylim=c(1E-6,10))
#show.noncentrality.distribs(log='x')
#show.noncentrality.distribs(log='xy')
#show.noncentrality.distribs(cdf=TRUE,log='x',ylim=c(0,1))
#show.noncentrality.distribs(cdf=FALSE,ylim=c(0,1),xlim=c(0,80))
#show.noncentrality.distribs(cdf=TRUE,ylim=c(0,1),xlim=c(0,80))
#v=seq(0,1,.1);rbind(v,ipsi(v))

return( make.basis.interps(J=J) )

}

show.mixed.noncentral.chisq.basis=function(B,v=seq(0,20,length=1001),cdf=TRUE,rescale=FALSE,J=B$J,...) {
  ##B=make.basis.interps(J)
  
  for (j in 1:J) {
    if (!cdf) {
      vv=B$gl[[j]](B$psi(v))*B$dpsi(v)
      if (rescale) vv=vv/max(vv[-1])
    } else {
      vv=B$Gl[[j]](B$psi(v))
    }
    if (j==1) plot(v,vv,type='l',col=j,lty=j,...) else {
      lines(v,vv,col=j,lty=j)
    }
  }
  if (!cdf) {
    vv=B$g0(B$psi(v))*B$dpsi(v)
    if (rescale) vv=vv/max(vv[-1])
    lines(v,vv,col='red')
  } else {
    vv=B$G0(B$psi(v))
    lines(v,vv,col='red')
  }
}

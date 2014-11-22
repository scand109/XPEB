#opt=bfgs_options();
#  Obtain the default options structure for bfgs operations.
#   see edit('bfgs_options') for details on the settings.
#  Mandatory: set opt$f and opt$df to the objective and gradient functions

bfgs_options=function() {

opt=list()
opt$trace=0; # 0: quiet
             # >=1: trace bfgs iterations
             # >=2: trace bfgs changes
             # >=3: trace bfgs internal steps
opt$diagnose=1; # 0: quiet
                # >=1: print convergence diagnostics at termination
opt$onfigure=0; # plot on figure #, or 0 for no plot.

opt$stop_tolerance=1E-9; # stop if the norm of the change in x is
                         # less than this
opt$exhaustion=1E-6; # stop if the change in change in objective is
                    # smaller than this (on asinh scale)

opt$step=function(x,s) x+s; # controls how the gradient-based step s
                     # modifies x

opt$initialstep=1.0; # used for the trial step size on the first
                   # step; just to help get the inital order of
                   # magnitude right. Results should adapt promptly anyway
opt$growstep=2.0; # factor by which steps are allowed to lengthen
                  # on each iteration (>1)
opt$maxstep=Inf; # use if you see steps bigger than you desire
opt$maxiter=10000; # stop for exhaustion 
opt$H0scale=1.0; # scale of initial guess at Hessian
opt$mu_min=0.0; # fudge if BFGS fails to maintain pos.def.; no fudge
            # should be necessary now.
opt$descentconst=0.0001; # small positive. (0< .. <1)
                         #  Explanation: for small enough steps,
                         #  descent should be as fast as the
                         #  gradient predicts, but even for larger
                         #  steps we want to retain some fraction
                         #  of that efficacy or our step is too big
return(opt)
}

#work alike for norm(x,'F'), except work on vectors
fnorm=function(x,dummy) sqrt(sum(x*x))

# state=bfgs_init(x,opt,d)
# x: the initial value for the optimizaton
# opt: a bfgs_options structure; use bfgs_options to get the default
bfgs_init=function(x,opt) {
state=list()
state$x=x;
state$f1=opt$f(x); state$fcallct=1;
state$df1=opt$df(x);  state$dfcallct=1;
state$H=opt$H0scale*diag(rep(1,length(state$df1)));
state$laststep=opt$initialstep/opt$growstep; 
# ((the divide by opt$growstep is because the code will grow it back before the
# first call))
return(state)
}

# Run a sequence of bfgs steps
# maxiter is optional override
#function [state f1v laststepv fcallctv dfcallctv descent]=bfgs_run(state,opt,maxiter=opt$maxiter)
bfgs_run=function(state,opt,maxiter=opt$maxiter) {

# pull state out
x=state$x;
f1=state$f1;
df1=state$df1;
H=state$H;
laststep=state$laststep;
fcallct=state$fcallct;
dfcallct=state$dfcallct;

#Record our rate of progress for plotting/analysis
f1v=matrix(0,maxiter+1,1);
laststepv=matrix(0,maxiter+1,1);
fcallctv=matrix(0,maxiter+1,1);
dfcallctv=matrix(0,maxiter+1,1);
f1v[1]=f1;
laststepv[1]=laststep;
fcallctv[1]=fcallct;
dfcallctv[1]=dfcallct;

if (opt$trace>=1) {
    #print('iter: objective (retries: stepsize)');
    #print(sprintf('0: %g',f1));
  }
for (kk in 1:maxiter) {
   # the stepsize heuristic used here is that largest step we want to allow is at most the largest allowable step or 110% of the last step's size.
   # this means we waste less time trying to take steps that are larger than the accuracy of our current Hessian estimate will allow
   # and similarly don't waste time taking steps that are unnecessarily small.
   b=bfgs_step_obj(x,f1,df1,H, min(opt$growstep*laststep,opt$maxstep), opt$mu_min, opt$f, opt$df, opt$step,fcallct,dfcallct,opt$descentconst,opt$trace);
   x=b$xn;f1=b$fn;df1=b$dfn;H=b$Hn;descent=b$descent;k=b$k;s=b$s1;fcallct=b$fcallct;dfcallct=b$dfcallct;

   laststep=fnorm(s,'F');
   H=(H+t(H))/2; # fix numerical asymmetry. (necessary?)
   f1v[kk+1]=f1;
   laststepv[kk+1]=laststep;
   fcallctv[kk+1]=fcallct;
   dfcallctv[kk+1]=dfcallct;
   if (opt$trace>=1) {
       #print(sprintf('%d: %g (%d:%g)',kk,f1,k,laststep));
       #if (opt$trace>=100) print(x)
     }
   if (laststep<opt$stop_tolerance) {break;}
   if (kk>=21 && asinh(f1)-asinh(f1v[kk-20])>-opt$exhaustion) {
     break;
   }
 }

# trim back the records
f1v=f1v[1:kk+1];
laststepv=laststepv[1:kk+1];
fcallctv=fcallctv[1:kk+1];
dfcallctv=dfcallctv[1:kk+1];

# wrap state back in
state$x=x;
state$f1=f1;
state$df1=df1;
state$H=H;
state$laststep=laststep;
state$fcallct=fcallct;
state$dfcallct=dfcallct;

if (descent==0 && opt$diagnose>=1) {
    #print('Convergence diagnostics:')
    #print(sprintf('f calls: %d',fcallct));
    #print(sprintf('df calls: %d',dfcallct));
    # check Hessian
    e=eigen(H)
    m=min(e$values);
    if (m<=0) {
        #print(sprintf('error: Hessian is not positive definite. Minimum eigen value is %g',m));
    } else {
        #print(sprintf('Hessian okay. Minimum e.value is: m=%g',m));
      }
    # check that gradient is an ascent direction
    sc=.0001;
    fv=c(state$f1);
    scv=c(0);
    ascent=0;
    for (j in 1:5) { #8 
        s=state$df1*sc;
        #try
          f2=opt$f(opt$step(state$x,s));
        #catch
        #  f2=NaN;

        if (ascent==0 && f2>f1) {
            ascent=j; # record the index of first ascent
          }
        if (ascent==0 && f2<f1) {
            ascent=-j; # record the negative index of first descent
          }
        fv=c(fv, f2);
        scv=c(scv, sc);
        sc=sc*10;
      }
    if (ascent>0) {
        #print(sprintf('Gradient ok. Ascent direction.'));
      }
    if (ascent==0) {
        #print(sprintf('Flat at bottom.'));
      }
    if (ascent<0) {
        #print('error: Gradient step is not an ascent. Error in gradient (or energy or step) calculation?');
      }
    #print(rbind(scv,fv));
  }

#kkv=0:kk;
#if (opt$onfigure>0)
#callctv=fcallctv+dfcallctv;
#figure(opt$onfigure);clf;
##semilogy(callctv, f1v); # valid if f1v>0
#plot(kkv, asinh(f1v)/log(10));
#hold on;
##semilogy(callctv, laststepv,'g:');
#plot(kkv, asinh(laststepv)/log(10),'g:');
#hold off;

out=list(state=state, f1v=f1v, laststepv=laststepv, fcallctv=fcallctv, dfcallctv=dfcallctv, descent=descent)
return(out)
}

#######################

#function [xn,fn,dfn,Hn,descent,k,s1,fcallct,dfcallct]=bfgs_step_obj(x,f1,df1,H, maxstep,mu_min,f,df,step,fcallct,dfcallct,descentconst,trace)
bfgs_step_obj=function(x,f1,df1,H, maxstep,mu_min,f,df,step,fcallct,dfcallct,descentconst,trace) {

maxdescent=-maxstep*fnorm(df1,'F'); # this is as large as the drop could
                               # be by going with s1=-df1*(maxstep/fnorm(df1,'F'))
                               # according to the secant approx.
if (maxdescent>-1E-6) {# if the small maxstep precludes progress
    if (trace>=2) {
        #print('bfgs: enabling larger step');
      }
    maxstep=Inf;
  }

targetstep=maxstep;
ts=truststep(df1,H,maxstep,mu_min);
s1=ts$y1;mu=ts$mu
# inspiring citation: Powell, Some global convergence properties of a variable metric algorithm for minimization
# http://books.google.com/books?id=3jK9ZlXeYdAC&lpg=PA53&ots=7kIbi4eqBx&dq=Some#20global#20convergence#20properties#20of#20a#20variable#20metric#20algorithm#20for#20minimization&lr&pg=PA57#v=onepage&q&f=false

# Check that it's a descent step...
for (k in 1:20) {
  xn=step(x,s1); # Call step to update x; in general, we allow x to
                 # have different size than s, but morally this is xn=x+s1;
  fn=f(xn); fcallct=fcallct+1;
  if (fn-f1<descentconst*t(s1) %*% df1) {
      if (trace>=3) {
          #print(sprintf('bfgs: descent: stepsize %g',fnorm(s1,'F')))
        }
      break;
    }
  if (trace>=3) {
          #print(sprintf('bfgs:       : stepsize %g',fnorm(s1,'F')))
        }
  targetstep=min(targetstep,fnorm(s1,'F'))/2;
  ts=truststep(df1,H,targetstep,mu);
  s1=ts$y1;mu=ts$mu
}
    
if (fn>=f1) {
  Hn=H;
  descent=0;
  s1=0*s1;
  xn=x;
  dfn=df1;
  if (trace>=2) {
      #print('bfgs: failure to descend');
    }
  res=list(xn=xn,fn=fn,dfn=dfn,Hn=Hn,descent=descent,k=k,s1=s1,fcallct=fcallct,dfcallct=dfcallct)
  return(res)
}

descent=1;

dfn=df(xn); dfcallct=dfcallct+1;

if (any(is.nan(dfn))) {
	#print('error: bfgs got nan\'s in df computation at xn:')
	#print(xn)
	stop('error in df')
}

g1=dfn-df1;
b=as.numeric(1/(t(g1) %*% s1));

# Only apply the BFGS update if b, below is positive; otherwise it fails to maintain pos.def.
# In the traditional BFGS, one does a full line-search, so this is ?assured??

#BFGS update
if (b>0) {
  c1=as.numeric(-1/(t(s1) %*% H %*% s1));
  Hn=H + b * g1 %*% t(g1) + c1*(H %*% s1) %*% (t(s1) %*% H);
} else {
  if (trace>=2) {
          #print('bfgs: skipping BFGS update')
  }
  Hn=H;
}


# Confirm secant condition:
#   Hn * s1 == g
# and postive definicy:
#   [ev,evv]=eig(Hn)
# If we'd had Hn earlier we'd propose moving to:
#  Hn \ -df1
# instead of:
#  H \ -df1
## n: [s1, mu]=truststep(df1,Hn,maxstep,0)
## o: [s1, mu]=truststep(df1,H,maxstep,0)

res=list(xn=xn,fn=fn,dfn=dfn,Hn=Hn,descent=descent,k=k,s1=s1,fcallct=fcallct,dfcallct=dfcallct)
return(res)
}
#####################

#function [y1,mu]=truststep(df1,H0,maxstep,mu)
truststep=function(df1,H0,maxstep,mu) {

tiny=1E-9;
#mu=0;

I=diag(rep(1,(length(df1))));

y1=solve((H0+mu*I), (-df1));
l1=sqrt(sum(y1^2));
if (l1<maxstep) {
  return(list(y1=y1,mu=mu))
}

phi1=1/maxstep-1/l1;

ct=0;
while (1) {
  # Do a newton update to find a mu-regularization value to get the length under maxstep (but close)
  # conceptually, we treat this as finding a zero of a function phi and use newton steps
  # (c.f. Lange p187 Trust Regions)
  # Here, phi(mu)=1/maxstep-1/fnorm(y(mu),'F'), where y(mu)=(H0+mu*I)\(-df1);

  dphi=as.numeric(-t(y1) %*% (solve( (H0+mu*I), y1)) / l1^3);
  if (dphi>0) {
      stop('H0 fails to be positive definite')
    }
  ct=ct+1;
  
  mustep=-phi1/(dphi);

  # Ensure "descent" (in terms of abs(phi1))
  for (k in 0:10) {
  
    mun=mu+mustep/2^k;
    if (mun<0) {
      mun=mu*(1-2^(-k));
    }
  
    if (ct>10) {
      #print(sprintf('ct=%d, k=%d, mu=%d, phi1=%d, dmu=%d',ct,k,mu,phi1,dphi));
      if (ct>100) {
	stop('convergence too slow in truststep!');
      }
    }
    y1n=solve((H0+mun*I), (-df1));
    l1n=sqrt(sum(y1n^2));
    if (l1n<1.01*maxstep && l1n>.99*maxstep) {
      y1=y1n;
      mu=mun;
      return(list(y1=y1,mu=mu))
    }
    phi1n=1/maxstep-1/l1n;
    if (abs(phi1n)<abs(phi1)) {
      break;
    }
  }
  if (abs(phi1n)>abs(phi1)) {
    #print('descent failed in truststep');
  }
  
  y1=y1n; #accept the descent
  l1=l1n;
  phi1=phi1n;
  #@add a stoping condition when mu just won't move much any more...
  if (abs(mun-mu)<tiny) {
    mu=mun;
    if (l1>1.01*maxstep) {
        #print('truststep failed to find a short enough step')
        #print(sprintf('ct=%d, k=%d, mu=%d, phi1=%d, dmu=%d',ct,k,mu,phi1,dphi));
      }
    return(list(y1=y1,mu=mu))
  }
  mu=mun; 
}
return(list(y1=y1,mu=mu))
}

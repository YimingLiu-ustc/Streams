module  Gauss_Err_Func
  Implicit None
  Private gammq , gammp , gser , gcf , gammln
contains

  Real FUNCTION fc_erf(x)
    REAL x
    if(x < 0.0)then
      fc_erf=-gammp(.5,x**2)
    else
      fc_erf=gammp(.5,x**2)
    endif
  END FUNCTION fc_erf

  Real FUNCTION fc_erfc(x)
    REAL x
    if(x < 0.)then
      fc_erfc=1.+gammp(.5,x**2)
    else
      fc_erfc=gammq(.5,x**2)
    endif
  END FUNCTION fc_erfc

  Real FUNCTION gammq(a,x)
    REAL a,x
    REAL gammcf,gamser,gln
    if(x<0.0 .or. a<=0.0) return
    if(x<a+1.)then
      call gser(gamser,a,x,gln)
      gammq=1.-gamser
    else
      call gcf(gammcf,a,x,gln)
      gammq=gammcf
    endif
  End Function gammq

  Real Function gammp(a,x)
    REAL a,x
    REAL gammcf,gamser,gln
    if(x<0..or.a<=0.) return !'bad arguments in gammp'
    if(x<a+1.)then
      call gser(gamser,a,x,gln)
      gammp=gamser
    else
      call gcf(gammcf,a,x,gln)
      gammp=1.-gammcf
    endif
  End Function gammp

  Subroutine gser(gamser,a,x,gln)
    INTEGER ITMAX
    REAL a,gamser,gln,x,EPS
    PARAMETER (ITMAX=100,EPS=3.e-7)
    INTEGER n
    REAL ap,del,sum
    gln=gammln(a)
    if(x<=0.)then
      if(x<0.) return ! 'x < 0 in gser'
      gamser=0.
      return
    endif
    ap=a
    sum=1./a
    del=sum
    do 11 n=1,ITMAX
      ap=ap+1.
      del=del*x/ap
      sum=sum+del
      if(abs(del)<abs(sum)*EPS)goto 1
11  continue
    return
1   gamser=sum*exp(-x+a*log(x)-gln)
  End Subroutine gser

  SUBROUTINE gcf(gammcf,a,x,gln)
    INTEGER ITMAX
    REAL a,gammcf,gln,x,EPS,FPMIN
    PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
    INTEGER i
    REAL an,b,c,d,del,h
    gln=gammln(a)
    b=x+1.-a
    c=1./FPMIN
    d=1./b
    h=d
    do 11 i=1,ITMAX
      an=-i*(i-a)
      b=b+2.
      d=an*d+b
      if(abs(d)<FPMIN)d=FPMIN
      c=b+an/c
      if(abs(c)<FPMIN)c=FPMIN
      d=1./d
      del=d*c
      h=h*del
      if(abs(del-1.)<EPS)goto 1
11  continue
    return ! 'a too large, ITMAX too small in gcf'
1   gammcf=exp(-x+a*log(x)-gln)*h
  END Subroutine gcf
    
  Real Function gammln(xx)
    REAL xx
    INTEGER j
    DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
    SAVE cof,stp
    DATA cof,stp/76.18009172947146d0,-86.50532032941677d0, &
   24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, &
   -.5395239384953d-5,2.5066282746310005d0/
    x=xx
    y=x
    tmp=x+5.5d0
    tmp=(x+0.5d0)*log(tmp)-tmp
    ser=1.000000000190015d0
    do 11 j=1,6
      y=y+1.d0
      ser=ser+cof(j)/y
11  continue
    gammln=tmp+log(stp*ser/x)
  End Function gammln
End Module Gauss_Err_Func

subroutine init_rt
!
! initialize Rayleigh-Taylor instability case
!
  use mod_streams
  use Gauss_Err_Func
  implicit none
!
  integer :: i,j,k,l,m
  real(mykind) :: dx, perda,delx,disamplify,disomega,disreducespeed,disshape,defecteddis,Epsilon,disturbedx
  real(mykind) :: rhoplus,rhominus,tempA,tempB,TempP
  disamplify=0.1
  delx = 0.01
  disomega=1._mykind
  disreducespeed=0.95
  Epsilon=0.1**10._mykind
  perda=0.01
  do k=1-ng,nz+ng
    do j=1-ng,ny+ng
      do i=1-ng,nx+ng
        disshape=disamplify*cos(2._mykind*pi*disomega*z(k))
        defecteddis=disshape*(disreducespeed**(abs(x(i)-disshape)/(abs(disamplify)+Epsilon))**3._mykind)
        disturbedx=x(i)-defecteddis
        rhoplus=(1+At)*EXP(-(1+At)/Fr*disturbedx)*(1+fc_erf(real(disturbedx/delx)))/2._mykind
        rhominus=(1-At)*EXP(-(1-At)/Fr*disturbedx)*(1-fc_erf(real(disturbedx/delx)))/2._mykind
        w(1,i,j,k)=rhoplus+rhominus
        tempA=(1._mykind+At)/Fr
        tempB=(1._mykind-At)/Fr
        TempP=0.5*Exp(-tempA*disturbedx)+0.5*Exp(-tempB*disturbedx)+&
        0.5*Exp(-tempA*disturbedx)*fc_erf(real(disturbedx/delx))-0.5*Exp(-tempB*disturbedx)*fc_erf(real(disturbedx/delx))+&  
        0.5*Exp((tempB*delx/2._mykind)**2._mykind)*(fc_erf(real(disturbedx/delx+tempB*delx/2._mykind))-fc_erf(real(tempB*delx/2._mykind)))-&   
        0.5*Exp((tempA*delx/2._mykind)**2._mykind)*(fc_erf(real(disturbedx/delx+tempA*delx/2._mykind))-fc_erf(real(tempA*delx/2._mykind)))
        w(5,i,j,k)=TempP*gm
      enddo
    enddo
  enddo

! rhou initialization
  do k=1-ng,nz+ng
    do j=1-ng,ny+ng
      do i=1-ng,nx+ng
        do m=2,4
          w(m,i,j,k) =0._mykind
        enddo
      enddo
    enddo
  enddo
! rhoe=p*gm initialization



end subroutine init_rt

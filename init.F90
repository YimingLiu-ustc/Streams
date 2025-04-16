subroutine init
!
! Providing initial conditions
!
  use mod_streams
  use Gauss_Err_Func
  implicit none
!
  integer :: i,j,k,m
  real(mykind) :: perda,delx,disamplify,disomega,disreducespeed,disshape,defecteddis,Epsilon,disturbedx
  real(mykind) :: rhoplus,rhominus,tempA,tempB,TempP
!
!
  disamplify=0.1
  delx = 0.01
  disomega=1._mykind
  disreducespeed=0.95
  Epsilon=0.1**10
  perda=0.01
!
!
  dfupdated      = .true.    ! 数字滤波更新标志
  ncyc0          = 0         ! number of cycles is set to zero
  istore         = 1         ! index of the first solution to store
  istore_restart = 1         ! index of the first solution to store restart
  telaps0        = 0._mykind ! time counter is set to zero
  stat_io        = 'rewind'  ! I/O status
!
! Freestream initialization (just to avoid 0 in double ghost nodes)

  do k=1-ng,nz+ng
    do j=1-ng,ny+ng
      do i=1-ng,nx+ng
        do m=1,nv
          w(m,i,j,k) = winf(m)
        enddo
      enddo
    enddo
  enddo

  ducros = .false.
  if (tresduc<=0._mykind) ducros = .true.
!
  if (idiski>=1) then !续算
!
! Reading solution from file
!
    if(io_type == 1) call readrst_serial  ! serial 串行 I/O
    if(io_type == 2) call readrst         ! parallel 并行 I/O
    if (iflow==-1) then
    elseif (iflow==0) then
      if (idiski > 1) call readstat1d
      call generatewmean_channel
    else
      if(io_type == 1) then
        call readdf_serial
        if (idiski > 1) call readstat2d_serial
      elseif(io_type == 2) then
        call readdf
        if (idiski > 1) call readstat2d
      endif
      call generatewmean
      if (masterproc) call target_reystress
      call mpi_bcast(amat_df,9*ny,mpi_prec,0,mp_cart,iermpi)
    endif
!
    istore          = int(telaps0/dtsave) + 1
    istore_restart  = int(telaps0/dtsave_restart) + 1
    stat_io         = 'append'
!
    if (idiski==1) then
      w_av    = 0._mykind
      w_av_1d = 0._mykind
      itav    = 0
    endif
!
  else ! start from scratch
!
! Setting the initial conditions
!
    w_av    = 0._mykind
    w_av_1d = 0._mykind
    itav    = 0
!
! rho1bc and rho2bc boundary data for temperature stratification
  do k=1,nzmax
    do j=1,nymax
      i=1
      disshape=disamplify*cos(2._mykind*pi*disomega*zg(k))
      defecteddis=disshape*(disreducespeed**(abs(xg(i)-disshape)/(abs(disamplify)+Epsilon))**3._mykind)
      disturbedx=xg(i)-defecteddis
      rhoplus=(1+At)*EXP(-(1+At)/Fr*disturbedx)*(1+fc_erf(real(disturbedx/delx)))/2._mykind
      rhominus=(1-At)*EXP(-(1-At)/Fr*disturbedx)*(1-fc_erf(real(disturbedx/delx)))/2._mykind
      rho1bc(j,k)=rhoplus+rhominus
      tempA=(1._mykind+At)/Fr
      tempB=(1._mykind-At)/Fr
      TempP=0.5*Exp(-tempA*disturbedx)+0.5*Exp(-tempB*disturbedx)+&
      0.5*Exp(-tempA*disturbedx)*fc_erf(real(disturbedx/delx))-0.5*Exp(-tempB*disturbedx)*fc_erf(real(disturbedx/delx))+&  
      0.5*Exp((tempB*delx/2._mykind)**2._mykind)*(fc_erf(real(disturbedx/delx+tempB*delx/2._mykind))-fc_erf(real(tempB*delx/2._mykind)))-&   
      0.5*Exp((tempA*delx/2._mykind)**2._mykind)*(fc_erf(real(disturbedx/delx+tempA*delx/2._mykind))-fc_erf(real(tempA*delx/2._mykind)))
      T1bc(j,k)=TempP/rho1bc(j,k)
    enddo
  enddo

  do k=1,nzmax
    do j=1,nymax
      i=nxmax
      disshape=disamplify*cos(2._mykind*pi*disomega*zg(k))
      defecteddis=disshape*(disreducespeed**(abs(xg(i)-disshape)/(abs(disamplify)+Epsilon))**3._mykind)
      disturbedx=xg(i)-defecteddis
      rhoplus=(1+At)*EXP(-(1+At)/Fr*disturbedx)*(1+fc_erf(real(disturbedx/delx)))/2._mykind
      rhominus=(1-At)*EXP(-(1-At)/Fr*disturbedx)*(1-fc_erf(real(disturbedx/delx)))/2._mykind
      rho2bc(j,k)=rhoplus+rhominus
      tempA=(1._mykind+At)/Fr
      tempB=(1._mykind-At)/Fr
      TempP=0.5*Exp(-tempA*disturbedx)+0.5*Exp(-tempB*disturbedx)+&
      0.5*Exp(-tempA*disturbedx)*fc_erf(real(disturbedx/delx))-0.5*Exp(-tempB*disturbedx)*fc_erf(real(disturbedx/delx))+&  
      0.5*Exp((tempB*delx/2._mykind)**2._mykind)*(fc_erf(real(disturbedx/delx+tempB*delx/2._mykind))-fc_erf(real(tempB*delx/2._mykind)))-&   
      0.5*Exp((tempA*delx/2._mykind)**2._mykind)*(fc_erf(real(disturbedx/delx+tempA*delx/2._mykind))-fc_erf(real(tempA*delx/2._mykind)))
      T2bc(j,k)=TempP/rho2bc(j,k)
    enddo
  enddo
!
    if (iflow==-1) then
      call init_windtunnel
    elseif (iflow==0) then
      call generatewmean_channel
      call init_channel
    elseif (iflow==3) then
      call init_rt
    else
      call generatewmean
      if (masterproc) call target_reystress
      call mpi_bcast(amat_df,9*ny,mpi_prec,0,mp_cart,iermpi)
      call initurb
    endif
!
  endif ! End of definition of the initial flowfield
!
end subroutine init

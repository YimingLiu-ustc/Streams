subroutine bcwall_adiabatic(ilat)
!
! Apply wall boundary conditions(adiabatic)
!
  use mod_streams
  use Gauss_Err_Func
  implicit none
!
  integer :: i,j,k,l,ilat
  real(mykind) :: rho,uu,vv,ww,qq,pp,tt,rhoe
  real(mykind), dimension(ng), device :: Tj
!
  if (ilat==1) then     ! left side
    !$cuf kernel do(2) <<<*,*>>>
    do k=1,nz
      do j=1,ny
        w_gpu(1,j,k,1)= rho1bc_gpu(j,k)
        !w_gpu(1,j,k,1) = 1-At
        w_gpu(1,j,k,2) = 0._mykind
        w_gpu(1,j,k,3) = 0._mykind
        w_gpu(1,j,k,4) = 0._mykind
        do l=1,ng
          rho  = w_gpu(1+l,j,k,1)
          uu   = w_gpu(1+l,j,k,2)/w_gpu(1+l,j,k,1)
          vv   = w_gpu(1+l,j,k,3)/w_gpu(1+l,j,k,1)
          ww   = w_gpu(1+l,j,k,4)/w_gpu(1+l,j,k,1)
          rhoe = w_gpu(1+l,j,k,5)
          qq   = 0.5_mykind*(uu*uu+vv*vv+ww*ww)
          pp   = gm1*(rhoe-rho*qq)
          tt   = pp/rho
          Tj(l)=tt
          !tt   = 2._mykind*twall-tt
          !rho  = pp/tt
          w_gpu(1-l,j,k,1) =  rho
          w_gpu(1-l,j,k,2) = -rho*uu
          w_gpu(1-l,j,k,3) = -rho*vv
          w_gpu(1-l,j,k,4) = -rho*ww
          w_gpu(1-l,j,k,5) = pp*gm+qq*rho
          ! w_gpu(1-l,j,k,1) = rho
          ! w_gpu(1-l,j,k,2) = rho*uu
          ! w_gpu(1-l,j,k,3) = rho*vv
          ! w_gpu(1-l,j,k,4) = rho*ww
          ! w_gpu(1-l,j,k,5) = pp*gm+qq*rho
        enddo
        w_gpu(1,j,k,5) = w_gpu(1,j,k,1)*gm*(4._mykind*Tj(1)-Tj(2))/3._mykind
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
  elseif (ilat==2) then ! right side
    !$cuf kernel do(2) <<<*,*>>>
    do k=1,nz
      do j=1,ny
        !  w_gpu(nx,j,k,1)=rho2bc_gpu(j,k)
        w_gpu(nx,j,k,1) = 1+At
        w_gpu(nx,j,k,2) = 0._mykind
        w_gpu(nx,j,k,3) = 0._mykind
        w_gpu(nx,j,k,4) = 0._mykind
        do l=1,ng
          rho  = w_gpu(nx-l,j,k,1)
          uu   = w_gpu(nx-l,j,k,2)/w_gpu(nx-l,j,k,1)
          vv   = w_gpu(nx-l,j,k,3)/w_gpu(nx-l,j,k,1)
          ww   = w_gpu(nx-l,j,k,4)/w_gpu(nx-l,j,k,1)
          rhoe = w_gpu(nx-l,j,k,5)
          qq   = 0.5_mykind*(uu*uu+vv*vv+ww*ww)
          pp   = gm1*(rhoe-rho*qq)
          tt   = pp/rho
          !tt   = 2._mykind*twall-tt
          !rho  = pp/tt
          Tj(l)=tt
          w_gpu(nx+l,j,k,1) =  rho
          w_gpu(nx+l,j,k,2) = -rho*uu
          w_gpu(nx+l,j,k,3)  = -rho*vv
          w_gpu(nx+l,j,k,4)  = -rho*ww
          w_gpu(nx+l,j,k,5)  = pp*gm+qq*rho
        enddo
        w_gpu(nx,j,k,5) = w_gpu(nx,j,k,1)*gm*(4._mykind*Tj(1)-Tj(2))/3._mykind
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
  elseif (ilat==3) then ! lower side
  elseif (ilat==4) then  ! upper side
  elseif (ilat==5) then  ! back side
  elseif (ilat==6) then  ! fore side
  endif
!
end subroutine bcwall_adiabatic
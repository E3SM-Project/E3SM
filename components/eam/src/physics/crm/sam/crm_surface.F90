module crm_surface_mod
  implicit none

contains

  subroutine crm_surface(ncrms,bflx)
    use vars
    use params
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd), intent (in) :: bflx(ncrms)
    real(crm_rknd) u_h0, tau00, tmp
    integer i,j,icrm

    !--------------------------------------------------------
    !$acc parallel loop async(asyncid)
    do icrm = 1 , ncrms
      uhl(icrm) = uhl(icrm) + dtn*utend(icrm,1)
      vhl(icrm) = vhl(icrm) + dtn*vtend(icrm,1)
      taux0(icrm) = 0.
      tauy0(icrm) = 0.
    end do
    !$acc parallel loop collapse(3) async(asyncid)
    do j=1,ny
      do i=1,nx
        do icrm = 1 , ncrms
          u_h0 = max(real(1.,crm_rknd),sqrt((0.5D0*(u(icrm,i+1,j,1)+u(icrm,i,j,1))+ug)**2+(0.5D0*(v(icrm,i,j+YES3D,1)+v(icrm,i,j,1))+vg)**2))
          tau00 = rho(icrm,1) * diag_ustar(z(icrm,1),bflx(icrm),u_h0,z0(icrm))**2
          fluxbu(icrm,i,j) = -(0.5D0*(u(icrm,i+1,j,1)+u(icrm,i,j,1))+ug-uhl(icrm))/u_h0*tau00
          fluxbv(icrm,i,j) = -(0.5D0*(v(icrm,i,j+YES3D,1)+v(icrm,i,j,1))+vg-vhl(icrm))/u_h0*tau00
          tmp = fluxbu(icrm,i,j)/dble(nx*ny)
          !$acc atomic update
          taux0(icrm) = taux0(icrm) + tmp
          tmp = fluxbv(icrm,i,j)/dble(nx*ny)
          !$acc atomic update
          tauy0(icrm) = tauy0(icrm) + tmp
        end do
      end do
    end do
  end subroutine crm_surface

  ! ----------------------------------------------------------------------
  !
  ! DISCLAIMER : this code appears to be correct but has not been
  !              very thouroughly tested. If you do notice any
  !              anomalous behaviour then please contact Andy and/or
  !              Bjorn
  !
  ! Function diag_ustar:  returns value of ustar using the below
  ! similarity functions and a specified buoyancy flux (bflx) given in
  ! kinematic units
  !
  ! phi_m (zeta > 0) =  (1 + am * zeta)
  ! phi_m (zeta < 0) =  (1 - bm * zeta)^(-1/4)
  !
  ! where zeta = z/lmo and lmo = (theta_rev/g*vonk) * (ustar^2/tstar)
  !
  ! Ref: Businger, 1973, Turbulent Transfer in the Atmospheric Surface
  ! Layer, in Workshop on Micormeteorology, pages 67-100.
  !
  ! Code writen March, 1999 by Bjorn Stevens
  !
  ! Code corrected 8th June 1999 (obukhov length was wrong way up,
  ! so now used as reciprocal of obukhov length)
  real(crm_rknd) function diag_ustar(z,bflx,wnd,z0)
    use params, only: crm_rknd
    !$acc routine seq
    implicit none
    real(crm_rknd), parameter      :: vonk =  0.4D0   ! von Karmans constant
    real(crm_rknd), parameter      :: g    = 9.81D0   ! gravitational acceleration
    real(crm_rknd), parameter      :: am   =  4.8D0   !   "          "         "
    real(crm_rknd), parameter      :: bm   = 19.3D0   !   "          "         "
    real(crm_rknd), parameter      :: eps  = 1.D-10 ! non-zero, small number
    real(crm_rknd), intent (in)    :: z             ! height where u locates
    real(crm_rknd), intent (in)    :: bflx          ! surface buoyancy flux (m^2/s^3)
    real(crm_rknd), intent (in)    :: wnd           ! wind speed at z
    real(crm_rknd), intent (in)    :: z0            ! momentum roughness height
    integer :: iterate
    real(crm_rknd)    :: lnz, klnz, c1, x, psi1, zeta, rlmo, ustar
    lnz   = log(z/z0)
    klnz  = vonk/lnz
    c1    = 3.14159D0/2.D0 - 3.D0*log(2.D0)
    ustar =  wnd*klnz
    if (bflx /= 0.0D0) then
      do iterate=1,8
        rlmo   = -bflx * vonk/(ustar**3 + eps)   !reciprocal of
        !obukhov length
        zeta  = min(real(1.,crm_rknd),z*rlmo)
        if (zeta > 0.D0) then
          ustar =  vonk*wnd  /(lnz + am*zeta)
        else
          x     = sqrt( sqrt( 1.0D0 - bm*zeta ) )
          psi1  = 2.D0*log(1.0D0+x) + log(1.0D0+x*x) - 2.D0*atan(x) + c1
          ustar = wnd*vonk/(lnz - psi1)
        end if
      end do
    end if
    diag_ustar = ustar
  end function diag_ustar
  ! ----------------------------------------------------------------------

  real(crm_rknd) function z0_est(z,bflx,wnd,ustar)
    use params, only: crm_rknd
    ! Compute z0 from buoyancy flux, wind, and friction velocity
    ! 2004, Marat Khairoutdinov
    implicit none
    !$acc routine seq
    real(crm_rknd), parameter      :: vonk =  0.4D0   ! von Karmans constant
    real(crm_rknd), parameter      :: g    = 9.81D0   ! gravitational acceleration
    real(crm_rknd), parameter      :: am   =  4.8D0   !   "          "         "
    real(crm_rknd), parameter      :: bm   = 19.3D0   !   "          "         "
    real(crm_rknd), parameter      :: eps  = 1.D-10 ! non-zero, small number
    real(crm_rknd), intent (in)    :: z             ! height where u locates
    real(crm_rknd), intent (in)    :: bflx          ! surface buoyancy flux (m^2/s^3)
    real(crm_rknd), intent (in)    :: wnd           ! wind speed at z
    real(crm_rknd), intent (in)    :: ustar         ! friction velocity
    real(crm_rknd)    :: lnz, klnz, c1, x, psi1, zeta, rlmo
    c1    = 3.14159D0/2.D0 - 3.D0*log(2.D0)
    rlmo   = -bflx*vonk/(ustar**3+eps)   !reciprocal of
    zeta   = min(real(1.,crm_rknd),z*rlmo)
    if (zeta >= 0.D0) then
      psi1 = -am*zeta
    else
      x     = sqrt( sqrt( 1.0D0 - bm*zeta ) )
      psi1  = 2.D0*log(1.0D0+x) + log(1.0D0+x*x) - 2.D0*atan(x) + c1
    end if
    lnz = max(real(0.,crm_rknd),vonk*wnd/(ustar + eps) + psi1)
    z0_est = z*exp(-lnz)
  end function z0_est
  ! ----------------------------------------------------------------------

end module crm_surface_mod

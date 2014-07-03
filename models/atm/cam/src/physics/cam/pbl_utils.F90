module pbl_utils
!-----------------------------------------------------------------------!
! Module to hold PBL-related subprograms that may be used with multiple !
! different vertical diffusion schemes.                                 !
!                                                                       !
! Public subroutines:                                                   !
!     
!     calc_obklen                                                       !
!                                                                       !
!------------------ History --------------------------------------------!
! Created: Apr. 2012, by S. Santos                                      !
!-----------------------------------------------------------------------!

use shr_kind_mod, only: r8 => shr_kind_r8

implicit none
private

! Public Procedures
!----------------------------------------------------------------------!
! Excepting the initialization procedure, these are elemental
! procedures, so they can accept scalars or any dimension of array as
! arguments, as long as all arguments have the same number of
! elements.
public pbl_utils_init
public calc_ustar
public calc_obklen
public virtem

real(r8), parameter :: ustar_min = 0.01_r8

real(r8) :: g         ! acceleration of gravity
real(r8) :: vk        ! Von Karman's constant
real(r8) :: cpair     ! specific heat of dry air
real(r8) :: rair      ! gas constant for dry air
real(r8) :: zvir      ! rh2o/rair - 1

contains

subroutine pbl_utils_init(g_in,vk_in,cpair_in,rair_in,zvir_in)

  !-----------------------------------------------------------------------!
  ! Purpose: Set constants to be used in calls to later functions         !
  !-----------------------------------------------------------------------!

  real(r8), intent(in) :: g_in       ! acceleration of gravity
  real(r8), intent(in) :: vk_in      ! Von Karman's constant
  real(r8), intent(in) :: cpair_in   ! specific heat of dry air
  real(r8), intent(in) :: rair_in    ! gas constant for dry air
  real(r8), intent(in) :: zvir_in    ! rh2o/rair - 1

  g = g_in
  vk = vk_in
  cpair = cpair_in
  rair = rair_in
  zvir = zvir_in

end subroutine pbl_utils_init

elemental subroutine calc_ustar( t,    pmid, taux, tauy, &
                                 rrho, ustar)

  !-----------------------------------------------------------------------!
  ! Purpose: Calculate ustar and bottom level density (necessary for      !
  !  Obukhov length calculation).                                         !
  !-----------------------------------------------------------------------!

  real(r8), intent(in) :: t         ! surface temperature
  real(r8), intent(in) :: pmid      ! midpoint pressure (bottom level)
  real(r8), intent(in) :: taux      ! surface u stress [N/m2]
  real(r8), intent(in) :: tauy      ! surface v stress [N/m2]

  real(r8), intent(out) :: rrho     ! 1./bottom level density
  real(r8), intent(out) :: ustar    ! surface friction velocity [m/s]

  rrho = rair * t / pmid
  ustar = max( sqrt( sqrt(taux**2 + tauy**2)*rrho ), ustar_min )
  
end subroutine calc_ustar

elemental subroutine calc_obklen( ths,  thvs, qflx, shflx, rrho, ustar, &
                                  khfs, kqfs, kbfs, obklen)

  !-----------------------------------------------------------------------!
  ! Purpose: Calculate Obukhov length and kinematic fluxes.               !
  !-----------------------------------------------------------------------!

  real(r8), intent(in)  :: ths           ! potential temperature at surface [K]
  real(r8), intent(in)  :: thvs          ! virtual potential temperature at surface
  real(r8), intent(in)  :: qflx          ! water vapor flux (kg/m2/s)
  real(r8), intent(in)  :: shflx         ! surface heat flux (W/m2)

  real(r8), intent(in)  :: rrho          ! 1./bottom level density [ m3/kg ]
  real(r8), intent(in)  :: ustar         ! Surface friction velocity [ m/s ]
  
  real(r8), intent(out) :: khfs          ! sfc kinematic heat flux [mK/s]
  real(r8), intent(out) :: kqfs          ! sfc kinematic water vapor flux [m/s]
  real(r8), intent(out) :: kbfs          ! sfc kinematic buoyancy flux [m^2/s^3]
  real(r8), intent(out) :: obklen        ! Obukhov length
  
  ! Need kinematic fluxes for Obukhov:
  khfs = shflx*rrho/cpair
  kqfs = qflx*rrho
  kbfs = khfs + zvir*ths*kqfs
  
  ! Compute Obukhov length:
  obklen = -thvs * ustar**3 / (g*vk*(kbfs + sign(1.e-10_r8,kbfs)))

end subroutine calc_obklen

elemental real(r8) function virtem(t,q)

  !-----------------------------------------------------------------------!
  ! Purpose: Calculate virtual temperature from temperature and specific  !
  !  humidity.                                                            !
  !-----------------------------------------------------------------------!

  real(r8), intent(in) :: t, q

  virtem = t * (1.0_r8 + zvir*q)

end function virtem

end module pbl_utils

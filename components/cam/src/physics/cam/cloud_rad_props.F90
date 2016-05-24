module cloud_rad_props
 ! a shell for cloud radiative properties.  Implement on rrtmg side but not for cam radiative transfer
use ppgrid,           only: pcols, pver, pverp
use physics_types,    only: physics_state

use radconstants,     only: nswbands, nlwbands


implicit none
private
save

public :: &
   cloud_rad_props_init !,     &
!   cloud_rad_props_get_sw, & ! return SW optical props of total bulk aerosols
!   cloud_rad_props_get_lw    ! return LW optical props of total bulk aerosols

!==============================================================================
contains
!==============================================================================

subroutine cloud_rad_props_init()

    return

end subroutine cloud_rad_props_init

!==============================================================================

!subroutine cloud_rad_props_get_sw(state,  &
!!                                  tau, tau_w, tau_w_g, tau_w_f,&
!                                  diagnosticindex, oldliq, oldice)
!   ! Arguments
!   type(physics_state), intent(in) :: state
!   
!   integer, optional,   intent(in) :: diagnosticindex      ! index (if present) to radiation diagnostic information
!
!   real(r8), intent(out) :: tau    (nswbands,pcols,pver) ! aerosol extinction optical depth
!   real(r8), intent(out) :: tau_w  (nswbands,pcols,pver) ! aerosol single scattering albedo * tau
!   real(r8), intent(out) :: tau_w_g(nswbands,pcols,pver) ! aerosol assymetry parameter * tau * w
!   real(r8), intent(out) :: tau_w_f(nswbands,pcols,pver) ! aerosol forward scattered fraction * tau * w
!
!   logical, optional, intent(in) :: oldliq,oldice
!
!
!  return
!end subroutine cloud_rad_props_get_sw
!!==============================================================================
!
!subroutine cloud_rad_props_get_lw(state,  cld_abs_od, diagnosticindex, oldliq, oldice, oldcloud)
!   ! Arguments
!   type(physics_state), intent(in)  :: state
!   
!   real(r8),            intent(out) :: cld_abs_od(nlwbands,pcols,pver) ! [fraction] absorption optical depth, per layer
!   integer, optional,   intent(in)  :: diagnosticindex
!   logical, optional,   intent(in)  :: oldliq  ! use old liquid optics
!   logical, optional,   intent(in)  :: oldice  ! use old ice optics
!   logical, optional,   intent(in)  :: oldcloud  ! use old optics for both (b4b)
!
!
!   return
!end subroutine cloud_rad_props_get_lw
!
end module cloud_rad_props

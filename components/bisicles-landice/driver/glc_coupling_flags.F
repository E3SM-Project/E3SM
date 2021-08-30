!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
!
! !MODULE: glc_coupling_flags - determine coupling flags
!
module glc_coupling_flags

! !DESCRIPTION:
!
! This module determines various coupling flags
!
! !REVISION HISTORY:
!  Author: Bill Sacks

! !USES:

  use glc_kinds_mod
  use glc_constants, only: stdout, zero_gcm_fluxes
  use glc_exit_mod

  implicit none
  private
  save

! !PUBLIC MEMBER FUNCTIONS:

  public :: has_ocn_coupling
  public :: has_ice_coupling

!EOP

!***********************************************************************

contains

!***********************************************************************
!BOP
! !IROUTINE: has_ocn_coupling
! !INTERFACE:

  logical function has_ocn_coupling()

! !DESCRIPTION:
! Returns true if glc has coupling to ocn

! !USES:

    use glc_route_ice_runoff, only: ice_needs_ocean_coupling

!EOP
!-----------------------------------------------------------------------

! Local variables

    logical :: liq_to_ocean
    logical :: ice_to_ocean

!-----------------------------------------------------------------------

    if (zero_gcm_fluxes) then
       has_ocn_coupling = .false.
    else
       ! For now, liquid runoff is always sent to the ocean
       liq_to_ocean = .true.
       
       ice_to_ocean = ice_needs_ocean_coupling()
       
       has_ocn_coupling = (liq_to_ocean .or. ice_to_ocean)
    end if

  end function has_ocn_coupling

!***********************************************************************

!***********************************************************************
!BOP
! !IROUTINE: has_ice_coupling
! !INTERFACE:

  logical function has_ice_coupling()

! !DESCRIPTION:
! Returns true if glc has coupling to ice

! !USES:

    use glc_route_ice_runoff, only: ice_needs_sea_ice_coupling

!EOP
!-----------------------------------------------------------------------

    if (zero_gcm_fluxes) then
       has_ice_coupling = .false.
    else
       has_ice_coupling = ice_needs_sea_ice_coupling()
    end if

  end function has_ice_coupling

!***********************************************************************

end module glc_coupling_flags

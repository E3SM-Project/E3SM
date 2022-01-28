!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
!
! !MODULE: glc_route_ice_runoff - route ice runoff to ocean or sea ice
!
module glc_route_ice_runoff

! !DESCRIPTION:
!
! This module handles the routing of the solid ice runoff flux (i.e., calving) to either
! the ocean or sea ice, depending on a control flag.
!
! !REVISION HISTORY:
!  Author: Bill Sacks

! !USES:

  use glc_kinds_mod
  use glc_constants, only: stdout
  use glc_exit_mod

  implicit none
  private
  save

! !PUBLIC MEMBER FUNCTIONS:

  public :: set_routing
  public :: route_ice_runoff
  public :: ice_needs_ocean_coupling
  public :: ice_needs_sea_ice_coupling

!EOP

  ! Define possible routing settings
  integer(int_kind), parameter :: ROUTING_NULL       = 0
  integer(int_kind), parameter :: ROUTING_TO_OCEAN   = 1
  integer(int_kind), parameter :: ROUTING_TO_SEA_ICE = 2

  ! the routing used for this run
  integer(int_kind) :: routing = ROUTING_NULL

!***********************************************************************

contains

!***********************************************************************
!BOP
! !IROUTINE: set_routing
! !INTERFACE:

  subroutine set_routing(routing_code)

! !DESCRIPTION:
! Sets the routing type, given a routing code.
! Possible routing codes are:
!  'ocn': all ice runoff goes to ocean
!  'ice': all ice runoff goes to sea ice

! !INPUT PARAMETERS:

    character(len=*), intent(in) :: routing_code    ! name of the destination for ice runoff

!EOP
!-----------------------------------------------------------------------

! Local variables

    character(len=*), parameter :: subname = 'set_routing'

!-----------------------------------------------------------------------
    
    select case (routing_code)
    case ('ocn')
       routing = ROUTING_TO_OCEAN
    case ('ice')
       routing = ROUTING_TO_SEA_ICE
    case default
       write (stdout,*) subname, ' ERROR: Unknown routing: ', trim(routing_code)
       call exit_glc(sigAbort, ' ')
    end select

  end subroutine set_routing


!***********************************************************************
!BOP
! !IROUTINE: route_ice_runoff
! !INTERFACE:

  subroutine route_ice_runoff(rofi, rofi_to_ocn, rofi_to_ice)

! !DESCRIPTION:
! Routes solid ice runoff to the appropriate destination(s).
! Assumes that set_routing has already been called

! !INPUT PARAMETERS:

    real(r8), intent(in)  :: rofi          ! total solid ice runoff computed by glc for one grid cell

! !OUTPUT PARAMETERS:

    real(r8), intent(out) :: rofi_to_ocn   ! ice runoff to send to ocean
    real(r8), intent(out) :: rofi_to_ice   ! ice runoff to send to sea ice

!EOP
!-----------------------------------------------------------------------

! Local variables

    character(len=*), parameter :: subname = 'route_ice_runoff'

!-----------------------------------------------------------------------

    select case (routing)
    case (ROUTING_TO_OCEAN)
       rofi_to_ocn = rofi
       rofi_to_ice = 0._r8
    case (ROUTING_TO_SEA_ICE)
       rofi_to_ocn = 0._r8
       rofi_to_ice = rofi
    case default
       write (stdout,*) subname, ' ERROR: Unknown routing: ', routing
       call exit_glc(sigAbort, ' ')
    end select

  end subroutine route_ice_runoff

!***********************************************************************

!***********************************************************************
!BOP
! !IROUTINE: ice_needs_ocean_coupling
! !INTERFACE:

  logical function ice_needs_ocean_coupling()

! !DESCRIPTION:
! Returns true if the ice runoff requires ocn coupling, false otherwise

!EOP
!-----------------------------------------------------------------------

! Local variables

    character(len=*), parameter :: subname = 'ice_needs_ocean_coupling'

!-----------------------------------------------------------------------

    select case (routing)
    case (ROUTING_TO_OCEAN)
       ice_needs_ocean_coupling = .true.
    case (ROUTING_TO_SEA_ICE)
       ice_needs_ocean_coupling = .false.
    case default
       write (stdout,*) subname, ' ERROR: Unknown routing: ', routing
       call exit_glc(sigAbort, ' ')
    end select

  end function ice_needs_ocean_coupling

!***********************************************************************

!***********************************************************************
!BOP
! !IROUTINE: ice_needs_sea_ice_coupling
! !INTERFACE:

  logical function ice_needs_sea_ice_coupling()

! !DESCRIPTION:
! Returns true if the ice runoff requires ice coupling, false otherwise

!EOP
!-----------------------------------------------------------------------

! Local variables

    character(len=*), parameter :: subname = 'ice_needs_sea_ice_coupling'

!-----------------------------------------------------------------------

    select case (routing)
    case (ROUTING_TO_OCEAN)
       ice_needs_sea_ice_coupling = .false.
    case (ROUTING_TO_SEA_ICE)
       ice_needs_sea_ice_coupling = .true.
    case default
       write (stdout,*) subname, ' ERROR: Unknown routing: ', routing
       call exit_glc(sigAbort, ' ')
    end select

  end function ice_needs_sea_ice_coupling

!***********************************************************************

end module glc_route_ice_runoff

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

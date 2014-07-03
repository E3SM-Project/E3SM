module ionosphere

!---------------------------------------------------------------------------------
! Dummy interface for waccmx/ionosphere module.
!---------------------------------------------------------------------------------

use shr_kind_mod,   only : r8 => shr_kind_r8            ! Real kind to declare variables
use physics_types,  only : physics_state, physics_ptend !Structures containing physics state and tendency variables
use physics_buffer, only : physics_buffer_desc
use abortutils,     only : endrun

implicit none
save
private

public ionos_init     ! Initialization
public ionos_register ! Registration of ionosphere variables in pbuf physics buffer
public ionos_intr     ! interface to actual ionosphere simulation routines

!==============================================================================
contains
!==============================================================================

subroutine ionos_init()
  
   call endrun('ionos_init: dummy interface should not be called')

end subroutine ionos_init

!==============================================================================     

subroutine ionos_register()

   call endrun('ionos_register: dummy interface should not be called')

end subroutine ionos_register

!==============================================================================

subroutine ionos_intr(state, ptend, pbuf, ztodt)

   type(physics_state), intent(in)    :: state       ! physics state structure
   type(physics_ptend), intent(inout) :: ptend       ! parameterization tendency structure
   type(physics_buffer_desc), pointer :: pbuf(:)     ! physics buffer
   real(r8),            intent(in)    :: ztodt       ! Two times model timestep (2 delta-t)

   call endrun('ionos_intr: dummy interface should not be called')

end subroutine ionos_intr

!===============================================================================

end module ionosphere

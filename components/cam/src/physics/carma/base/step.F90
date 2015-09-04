! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  This routine performs all calculations necessary to take one timestep.
!!
!! @author McKie
!! @version Oct-1995
subroutine step(carma, cstate, rc)

	! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod

  implicit none

  type(carma_type), intent(in)         :: carma   !! the carma object
  type(carmastate_type), intent(inout) :: cstate  !! the carma state object
  integer, intent(inout)               :: rc      !! return code, negative indicates failure
  

  ! Iterate over each column. Each of these columns should be independent, so
  ! the work for each column could be done by a different thread.
  
	! Do pre-timestep processing
	if (rc >= 0) call prestep(carma, cstate, rc)

	! Update model state at new time
	if (rc >= 0) call newstate(carma, cstate, rc)

  ! Return to caller with one timestep taken
  return
end

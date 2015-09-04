! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This determines the maximum particle concentration for each group in each
!! gridbox. This can be used to make calculations more efficient by skipping
!! calculations when concentrations are low
!!
!! @author Chuck Bardeen
!! @version Nov 2009
subroutine maxconc(carma, cstate, iz, rc)

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
  integer, intent(in)                  :: iz      !! z index
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

	! Locals
	integer       :: igrp
	integer       :: iep


  ! Find maximum particle concentration for each spatial grid box
  !  (in units of cm^-3)
  do igrp = 1,NGROUP
    iep = ienconc(igrp)

    pconmax(iz,igrp) = maxval(pc(iz,:,iep))

    pconmax(iz,igrp) = pconmax(iz,igrp) &
                            / xmet(iz)        &
                            / ymet(iz)        &
                            / zmet(iz)
  enddo  ! igrp

  return
end

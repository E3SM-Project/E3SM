#include "carma_globaer.h"

!! This routine calculates new particle concentrations from coagulation
!! microphysical processes.
!!
!!  The basic form from which the solution is derived is:
!!
!!    ( new_value - old_value ) / dtime = source_term - loss_rate*new_value
!!
!! This routine derived from psolve.f code, in which particle concentrations
!! due to coagulation were formerly included, before the relatively slow
!! coagulation calcs were separated from the other microphysical processes
!! so that time splitting could be applied to these fast & slow calcs.
!!
!! @author Bill McKie
!! @version Sep-1997
subroutine csolve(carma, cstate, ibin, ielem, rc)

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
  integer, intent(in)                  :: ibin    !! bin index
  integer, intent(in)                  :: ielem   !! element index
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Local Variables
  integer                        :: igroup
  real(kind=f)                   :: xyzmet(NZ)
  real(kind=f)                   :: ppd(NZ)
  real(kind=f)                   :: pls(NZ)


  ! Define current group & particle number concentration element indices
  igroup = igelem(ielem)         ! particle group

  ! Metric scaling factor
  xyzmet = xmet(:) * ymet(:) * zmet(:)

  ! Compute total production rate due to coagulation
  ppd = coagpe(:,ibin,ielem) / xyzmet(:)

  ! Compute total loss rate due to coagulation
  pls = coaglg(:,ibin,igroup) / xyzmet(:)

  ! Update net particle number concentration during current timestep
  ! due to production and loss rates for coagulation
  pc(:,ibin,ielem) = ( pc(:,ibin,ielem) &
                       + dtime * ppd(:) ) &
                       /  ( ONE + pls(:) * dtime )

  return
end

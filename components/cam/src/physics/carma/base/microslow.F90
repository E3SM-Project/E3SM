! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  This routine drives the potentially slower microphysics calculations.
!!
!!  Originally part of microphy.  Now in this separate routine to allow
!!  time splitting of coagulation at a different timestep size from
!!  other microphysical calcs.
!!
!! @author McKie
!! @version Sep-1997
subroutine microslow(carma, cstate, rc)

  ! carma types defs
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

  ! Local Declarations
  integer   :: ibin
  integer   :: ielem


  
  !  Set production terms and loss rates due to slow microphysics
  !  processes (coagulation) to zero.
  coagpe(:,:,:) = 0._f
  coaglg(:,:,:) = 0._f

  ! Calculate (implicit) particle loss rates for coagulation.
  call coagl(carma, cstate, rc)

  ! Calculate particle production terms and solve for particle 
  ! concentrations at end of time step.
  !
  ! NOTE: The order of elements required by CARMA to work with the
  ! element loop first is: if you have a group that is both a source
  ! and product of coagulation, then it needs to come after the
  ! other group that participates in that coagulation in the element
  ! table. For example, icoag(2,1) = 1 will not work, but
  ! icoag(2,1) = 2 should work.
  do ielem = 1,NELEM
    do ibin = 1,NBIN
      call coagp(carma, cstate, ibin, ielem, rc)
      call csolve(carma, cstate, ibin, ielem, rc)
    enddo
  enddo
  
  ! Return to caller with new particle concentrations.
  return
end

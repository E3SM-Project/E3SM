! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine zeroes the fast microphysics sinks and sources, 
!! at one spatial point per call.
!!
!! @author Andy Ackerman
!! @version Oct-1997
subroutine zeromicro(carma, cstate, iz, rc)

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
  integer, intent(in)                  :: iz      !! vertical index
  integer, intent(inout)               :: rc      !! return code, negative indicates failure
  

  ! Set production terms and loss rates due to nucleation, growth,
  ! and evaporation to zero.  Also set index of smallest bin nuceleated
  ! during time step equal to <NBIN> first time through spatial loop.
  
  if (do_grow) then

    phprod         = 0._f
    rlprod         = 0._f
    dtpart(iz,:,:) = 0._f

    if (NGAS > 0) gasprod(:) = 0._f

    rhompe(:, :)  = 0._f
    rnucpe(:,:)   = 0._f
    growpe(:,:)   = 0._f
    evappe(:,:)   = 0._f
    rnuclg(:,:,:) = 0._f
    growlg(:,:)   = 0._f
    evaplg(:,:)   = 0._f

  end if

  ! Return to caller with fast microphysics sinks and sources zeroed.
  return
end

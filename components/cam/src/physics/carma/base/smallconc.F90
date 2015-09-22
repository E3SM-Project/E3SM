! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine ensures limits all particle concentrations in a grid box
!! to SMALL_PC.  In bins where this limitation results in the particle 
!! concentration changing, the core mass fraction and second moment fraction 
!! are set to <FIX_COREF>. 
!!
!! @author Andy Ackerman
!! @version Oct-1997
subroutine smallconc(carma, cstate, iz, ibin, ielem, rc)

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
  integer, intent(in)                  :: ibin    !! bin index
  integer, intent(in)                  :: ielem   !! element index
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

	! Locals
	integer       :: ig
	integer       :: ip
	real(kind=f)  :: small_val


  ig = igelem(ielem)
  ip = ienconc(ig)


  ! Element is particle concentration
  if (ielem == ip) then
     pc(iz,ibin,ielem) = max(pc(iz,ibin,ielem), SMALL_PC)
  else
  
    ! Element is core mass
    if ((itype(ielem) .eq. I_COREMASS) .or. (itype(ielem) .eq. I_VOLCORE)) then
      small_val = SMALL_PC * rmass(ibin,ig) * FIX_COREF
    
    ! Element is core second moment
    elseif (itype(ielem) .eq. I_CORE2MOM) then
      small_val = SMALL_PC * (rmass(ibin,ig) * FIX_COREF)**2
    end if

    ! Reset if either the particle concentration or the element mass are too small.
    if ((pc(iz,ibin,ip) <= SMALL_PC) .or. (pc(iz,ibin,ielem) < small_val)) then
      pc(iz,ibin,ielem) = small_val
    endif
  endif

  ! Return to caller with particle concentrations limited to SMALL_PC
  return
end

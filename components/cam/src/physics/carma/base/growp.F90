! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine calculates particle source terms due to growth <growpe>
!! for one particle size bin at one spatial grid point per call.
!!
!! @author Andy Ackerman
!! @version Dec-1995
subroutine growp(carma, cstate, iz, ibin, ielem, rc)

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

  ! Local declarations
  integer                              :: igroup  ! group index
  integer                              :: iepart
  

  ! Define group & particle # concentration indices for current element
  igroup = igelem(ielem)      ! target particle group 
  iepart = ienconc(igroup)	  ! target particle number concentration element

  ! Calculate production terms due to condensational growth <growpe>
  ! only if group to which element belongs grows.
  if( igrowgas(iepart) .ne. 0 .and. ibin .ne. 1 )then

    ! Bypass calculation if few droplets are present 
    if( pconmax(iz,igroup) .gt. FEW_PC )then
      growpe(ibin,ielem) = pc(iz,ibin-1,ielem) * growlg(ibin-1,igroup) 
    endif
  endif

  ! Return to caller with growth production terms evaluated.
  return
end

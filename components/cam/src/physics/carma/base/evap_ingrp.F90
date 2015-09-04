! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine calculates particle source terms <evappe> of droplets
!! evaporating within a particle group.
!!
!! Distinct evaporation of cores has not been treated.
!!
!! @author Andy Ackerman
!! @version Aug-2001
subroutine evap_ingrp(carma,cstate,iz,ibin,ig,ip,rc)

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
  integer, intent(in)                  :: ig      !! group index
  integer, intent(in)                  :: ip
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Local declarations
  integer                              :: ie  
  integer                              :: isub


  ! For a single group, the core mass fraction is 0.
  cmf(ibin,ig) = 0.0_f
  
  ! The smallest bin cannot be a source to smaller bins in same group
  if( ibin .eq. 1 )then
    return
  endif

  ! Evaluate evaporation source term <evappe> for all elements in group
  do isub = 1, nelemg(ig)
    ie = ip + isub - 1
    evappe(ibin-1,ie) = evappe(ibin-1,ie) + &
      pc(iz,ibin,ie)*evaplg(ibin,ig)
  enddo

  return
end

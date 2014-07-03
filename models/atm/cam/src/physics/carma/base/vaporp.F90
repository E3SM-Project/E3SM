! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  This routine calculates the vapor pressure for all gases at one altitude.
!!
!!  <pvapl> and <pvapi> are vapor pressures in units of [dyne/cm^2]
!!
!!  Uses temperature <t> as input.
!!
!! @author Andy Ackerman
!! @version Dec-1995
subroutine vaporp(carma, cstate, iz, igas, rc)

  !     types
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
  integer, intent(in)                  :: igas    !! gas index
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Each gas should have a vapor pressure routine specified for it.
  !
  ! As new gases are supported, this table should be expanded with new entries for
  ! the appropriate vapor pressure rotuines.
  select case(ivaprtn(igas))

    case (I_VAPRTN_H2O_BUCK1981)
      call vaporp_h2o_buck1981(carma, cstate, iz, rc, pvapl(iz, igas), pvapi(iz, igas))
    
    case(I_VAPRTN_H2O_MURPHY2005)
      call vaporp_h2o_murphy2005(carma, cstate, iz, rc, pvapl(iz, igas), pvapi(iz, igas))
    
    case(I_VAPRTN_H2O_GOFF1946)
      call vaporp_h2o_goff1946(carma, cstate, iz, rc, pvapl(iz, igas), pvapi(iz, igas))
    
    case(I_VAPRTN_H2SO4_AYERS1980)
      call vaporp_h2so4_ayers1980(carma, cstate, iz, rc, pvapl(iz, igas), pvapi(iz, igas))

    case default
      if (do_print) write(LUNOPRT,*) "vaporp:: ERROR - Unknown vapor pressure routine  (", ivaprtn(igas), &
        ") for gas (", igas, ")."
      rc = RC_ERROR
      return
  end select

  ! Return to caller with vapor pressures evaluated.
  return
end

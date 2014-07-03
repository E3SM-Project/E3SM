! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! Calculates the vapor pressure of water vapor over liquid water and ice according
!! to the parameterization of Murphy & Koop [2005].
!!
!! NOTE: <pvapl> and <pvapi> are vapor pressures in units of [dyne/cm^2]
!!
!! @author  Chuck Bardeen
!! @version May-2009
subroutine vaporp_h2o_murphy2005(carma, cstate, iz, rc, pvap_liq, pvap_ice)

	!     types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod

	implicit none

  type(carma_type), intent(in)         :: carma     !! the carma object
  type(carmastate_type), intent(inout) :: cstate    !! the carma state object
  integer, intent(in)                  :: iz        !! z index
  real(kind=f), intent(out)            :: pvap_liq  !! vapor pressure wrt liquid [dyne/cm2]
  real(kind=f), intent(out)            :: pvap_ice  !! vapor pressure wrt ice [dyne[cm2]
  integer, intent(inout)               :: rc        !! return code, negative indicates failure

  ! Local declarations
  real(kind=f)            :: tt
  

  ! Saturation vapor pressure over liquid water and water ice from
  !   Murphy and Koop, Quart. J. Roy. Meteo. Soc., 131, 1539-1565, [2005].
  tt = t(iz)
           
  pvap_liq = 10.0_f * exp(54.842763_f - (6763.22_f / tt) - (4.210_f * log(tt)) + (0.000367_f * tt) + &
             (tanh(0.0415_f * (tt - 218.8_f)) * &
              (53.878_f - (1331.22_f / tt) - (9.44523_f * log(tt)) + 0.014025_f * tt)))

  pvap_ice = 10.0_f * exp(9.550426_f - (5723.265_f / tt) + (3.53068_f * log(tt)) - (0.00728332_f * tt))

  ! Check to see whether temperature is ouside range of validity for the parameterization.
  !
  ! pvapl is defined for 123 < T < 332 K
  ! pvapi is defined for T > 110 K
  !
  ! NOTE: Don't stop the simulation if the limits are exceeded.
!  if ((t(iz) .le. 123.0_f) .or. (t(iz) .ge. 332.0_f)) then
!    if (do_print) write(LUNOPRT,*) 'vaporp_h2o_murphy2005::WARNING - Temperature', t(iz), &
!         ' out of range at iz = ', iz, "lat=", lat, "lon=", lon
!    rc = RC_WARNING
!  endif

  ! Return to caller with vapor pressures evaluated.
  return
end

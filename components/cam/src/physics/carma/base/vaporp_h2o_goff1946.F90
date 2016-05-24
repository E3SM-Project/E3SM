! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! Calculates the vapor pressure of water vapor over liquid water and ice according
!! to the parameterization of Goff & Gratch [1946] as used in CAM (wv_saturation.F90).
!!
!! NOTE: <pvapl> and <pvapi> are vapor pressures in units of [dyne/cm^2]
!!
!! @author  Chuck Bardeen
!! @version Dec-2010
subroutine vaporp_h2o_goff1946(carma, cstate, iz, rc, pvap_liq, pvap_ice)

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
  !   Goff and Gatch, [1946].
  tt = t(iz)
           
  pvap_liq = 10.0_f * 10._f**(-7.90298_f * (373.16_f / tt - 1._f) + &
             5.02808_f * log10(373.16_f / tt) - &
             1.3816e-7_f * (10._f**(11.344_f * (1._f - tt / 373.16_f)) - 1._f) + &
             8.1328e-3_f * (10._f**(-3.49149_f * (373.16_f / tt - 1._f)) - 1._f) + &
             log10(1013.246_f)) * 100._f                      

  pvap_ice = 10.0_f * 10._f**(-9.09718_f * (273.16_f / tt - 1._f) - 3.56654_f * &
             log10(273.16_f / tt) + 0.876793_f * (1._f - tt / 273.16_f) + &
             log10(6.1071_f)) * 100._f

  ! Check to see whether temperature is ouside range of validity for the parameterization.
  !
  ! pvapl is defined for -50 C < T < 102 C ,  Gibbons [1990]
  ! pvapi is defined for T > -100 C
  !
  ! NOTE: Don't stop the simulation if the limits are exceeded.
  if ((t(iz) .le. 173.0_f) .or. (t(iz) .ge. 375.0_f)) then
!    if (do_print) then
!       write(LUNOPRT,*) 'vaporp_h2o_goff1946::WARNING - Temperature', t(iz), &
!            ' out of range at iz = ', iz, "lat=", lat, "lon=", lon
!    end if
    rc = RC_WARNING
  endif

  ! Return to caller with vapor pressures evaluated.
  return
end

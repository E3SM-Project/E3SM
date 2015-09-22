! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! Calculates the vapor pressure of water vapor over liquid water and ice according
!! to the parameterization of Buck [1981].
!!
!! NOTE: <pvapl> and <pvapi> are vapor pressures in units of [dyne/cm^2]
!!
subroutine vaporp_h2o_buck1981(carma, cstate, iz, rc, pvap_liq, pvap_ice)

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
  real(kind=f), intent(out)            :: pvap_liq  !! vapor pressure wrt liquid
  real(kind=f), intent(out)            :: pvap_ice  !! vapor pressure wrt ice
  integer, intent(inout)               :: rc        !! return code, negative indicates failure

  ! Local declarations
  
  !  Define coefficients in Buck's formulation for saturation vapor pressures
  !  Table 2
  !
  ! Ice: valid temperature interval -80 - 0 C
  real(kind=f), parameter :: BAI = 6.1115e2_f
  real(kind=f), parameter :: BBI = 23.036_f
  real(kind=f), parameter :: BCI = 279.82_f
  real(kind=f), parameter :: BDI = 333.7_f

  ! Liquid: valid temperature interval -40 - +50 C
  real(kind=f), parameter :: BAL = 6.1121e2_f
  real(kind=f), parameter :: BBL = 18.729_f
  real(kind=f), parameter :: BCL = 257.87_f
  real(kind=f), parameter :: BDL = 227.3_f
  
  real(kind=f)            :: tt


  ! Saturation vapor pressure over liquid water and water ice from
  !   Buck [J. Atmos. Sci., 20, 1527, 1981]
  tt = t(iz) - 273.16_f
  
  pvap_liq = BAL * exp( (BBL - tt/BDL)*tt / (tt + BCL) )
  pvap_ice = BAI * exp( (BBI - tt/BDI)*tt / (tt + BCI) )

  ! Check to see whether temperature is ouside range of validity for the parameterization.
  !
  ! NOTE: Don't stop the simulation if the limits are exceeded.
  if (pvap_liq .le. 1.e-13_f) then
    if (do_print) write(LUNOPRT,*) 'vaporp_buck1981::WARNING - Temperature (', t(iz), ') too small for iz = ', iz
    rc = RC_WARNING
  endif

  ! Return to caller with vapor pressures evaluated.
  return
end

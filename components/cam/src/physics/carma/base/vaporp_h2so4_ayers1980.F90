! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  Calculates the vapor pressure for sulfuric acid.
!!
!!  <pvap_liq> and <pvap_ice> are vapor pressures in units of [dyne/cm^2]
!!
!!  Created   Dec-1995  (Ackerman) 
!!  Modified  Sep-1997  (McKie)
!!  Modified Jul-2001 (Mills)
!!
!!  NOTE: To calculate vapor pressure of H2SO4 water vapor pressure (pvapl(iz, igash2o))
!!  should be calculated before this calculation.
!!
!! @author Mike Mills, Tianyi Fan
!! @version Feb-2011
subroutine vaporp_H2SO4_Ayers1980(carma, cstate, iz, rc, pvap_liq, pvap_ice)
!     types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod
  use sulfate_utils
  
  implicit none

  type(carma_type), intent(in)         :: carma     !! the carma object
  type(carmastate_type), intent(inout) :: cstate    !! the carma state object
  integer, intent(in)                  :: iz        !! z index
  real(kind=f), intent(out)            :: pvap_liq  !! vapor pressure wrt liquid [dyne/cm2]
  real(kind=f), intent(out)            :: pvap_ice  !! vapor pressure wrt ice [dyne[cm2]
  integer, intent(inout)               :: rc        !! return code, negative indicates failure

  ! Local declarations
  real(kind=f)            :: gc_cgs                          ! water vapor mass concentration [g/cm3]
  real(kind=f)            :: fk1, fk4, fk4_1, fk4_2 
  real(kind=f)            :: factor_kulm                     ! Kulmala correction terms
  real(kind=f)            :: en, temp
  real(kind=f)            :: sulfeq
  real(kind=f), parameter :: t0_kulm     = 340._f            !  T0 set in the low end of the Ayers measurement range (338-445K)
  real(kind=f), parameter :: t_crit_kulm = 905._f            !  Critical temperature = 1.5 * Boiling point
  real(kind=f), parameter :: fk0         = -10156._f / t0_kulm + 16.259_f   !  Log(Kulmala correction factor)
  real(kind=f), parameter :: fk2         = 1._f / t0_kulm
  real(kind=f), parameter :: fk3         = 0.38_f / (t_crit_kulm - t0_kulm)    
  

  ! Saturation vapor pressure of sulfuric acid
  !  
  ! Don't allow saturation vapor pressure to underflow at very low temperatures
  temp=max(t(iz),140._f)
  
  ! Convert water vapor concentration to g/cm3:
  gc_cgs = gc(iz, igash2o) / (xmet(iz) * ymet(iz) * zmet(iz))
  
  ! Compute the sulfate composition based on Hanson parameterization
  ! to temperature and water vapor concentration.
  wtpct(iz) = wtpct_tabaz(carma, temp, gc_cgs, pvapl(iz, igash2o), rc)

  ! Parameterized fit to Giauque's (1959) enthalpies v. wt %:
  en = 4.184_f * (23624.8_f - 1.14208e8_f / ((wtpct(iz) - 105.318_f)**2 + 4798.69_f))
  en = max(en, 0.0_f)
 
  ! Ayers' (1980) fit to sulfuric acid equilibrium vapor pressure:
  ! (Remember this is the log)
  ! SULFEQ=-10156/Temp+16.259-En/(8.3143*Temp)
  !
  ! Kulmala correction (J. CHEM. PHYS. V.93, No.1, 1 July 1990)
  fk1   = -1._f / temp
  fk4_1 = log(t0_kulm / temp)
  fk4_2 = t0_kulm / temp
  fk4   = 1.0_f + fk4_1 - fk4_2
  factor_kulm = fk1 + fk2 + fk3 * fk4

  ! This is for pure H2SO4
  sulfeq = fk0 + 10156._f * factor_kulm
           
  ! Adjust for WTPCT composition:
  sulfeq = sulfeq - en / (8.3143_f * temp)
      
  ! REMEMBER TO TAKE THE EXPONENTIAL!	
  sulfeq = exp(sulfeq)
  
  ! BUT this is in Atmospheres.  Convert ==> dynes/cm2
  pvap_liq = sulfeq * 1.01325e6_f  
  pvap_ice = sulfeq * 1.01325e6_f 
  
  return
end

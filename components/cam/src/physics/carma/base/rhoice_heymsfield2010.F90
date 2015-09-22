! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine calculates the effective ice densities for each bin, based upon
!! the parameterization of Heymsfield et al. [2010].
!!
!! @author   Chuck Bardeen
!! @ version March 2010
!!
!! @see CARMAELEMENT_Create
subroutine rhoice_heymsfield2010(carma, rhoice, igroup, regime, rho, aratelem, rc)

  ! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod

  implicit none

  type(carma_type), intent(in)         :: carma       !! the carma object
  real(kind=f), intent(in)             :: rhoice      !! ice density(g/cm3)
  integer, intent(in)                  :: igroup      !! group index
  character(len=4), intent(in)         :: regime      !! crystal regime [warm | cold | conv]
  real(kind=f), intent(out)            :: rho(NBIN)   !! crystal density per bin (g/cm3)
  real(kind=f), intent(out)            :: aratelem(NBIN)  !! projected area ratio ()
  integer, intent(inout)               :: rc          !! return code, negative indicates failure

  ! Local declarations
  integer                              :: ibin        ! bin index
  real(kind=f)                         :: a           ! scalar coefficient from Heysfield and Schmitt [2010]
  real(kind=f), parameter              :: b = 2.1_f   ! exponential coefficient from Heysfield and Schmitt [2010]
  real(kind=f)                         :: rbin        ! predicated crystal radius (cm)
  real(kind=f)                         :: dmax        ! maximum diameter
  real(kind=f)                         :: totalmass   ! bin mass
  
1 format(/,'rhoice_heymsfield2010::ERROR - unknown ice regime (', a, ').')

  
  rc = RC_OK

! Figure out the 'a' coefficient.
  if (regime == "deep") then
    a = 1.10e-2_f
  else if (regime == "conv") then
    a = 6.33e-3_f
  else if (regime == "cold") then
    a = 5.74e-3_f
  else if (regime == "avg") then
    a = 5.28e-3_f
  else if (regime == "synp") then
    a = 4.22e-3_f
  else if (regime == "warm") then
    a = 3.79e-3_f
  else
    if (do_print) write(LUNOPRT,1) regime
    rc = RC_ERROR
    return
  end if

  ! Get the starting mass for the first bin and the volume ratio from the CARMA_GROUP. This
  ! call is used before initialization has happened, so the bin structure hasn't been
  ! determined yet.
  
  do ibin = 1, NBIN
  
    ! Determine the total mass of the particle.
    !
    ! NOTE: This needs to match the logic in setupbins.F90, so that the ice density
    ! and radii will be determined properly.
    totalmass = rmassmin(igroup) * (rmrat(igroup)**(ibin-1))
    
    ! Determine the radius of the particle from Heymsfield et al. [2010].
    !
    ! m(D) = a * D ^ b (all in cgs units)
    rbin = ((totalmass / a) ** (1._f/b)) / 2._f
    
    ! Determine the density of an equivalent sphere.
    rho(ibin) = totalmass / ((4._f / 3._f ) * PI * (rbin ** 3._f))
    
    ! Don't let the density be larger than the bulk density of ice. This
    ! will happen for r < ~ 50 um in the parameterization, but this is
    ! not physical.
    rho(ibin) = min(rho(ibin), rhoice)
    
    ! Determine the area ratio based on the formulation given in Schmitt and Heymsfield
    ! [2009].
    dmax = 2._f * rbin
    
    if (dmax <= 200.e-4_f) then
      aratelem(ibin) = exp(-38._f * dmax)
    else
      aratelem(ibin) = 0.16_f * (dmax ** (-0.27_f))
    end if

  end do

end subroutine rhoice_heymsfield2010

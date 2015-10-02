! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  Calculates particle production rates due to heterogeneous
!!  nucleation <rnuclg>:
!!
!!  This was moved from sulfnuc to make the code more manageable.
!!
!!  @author Mike Mills, Chuck Bardeen
!!  @version Jun-2013
subroutine sulfhetnucrate(carma, cstate, iz, igroup, nucbin, h2o, h2so4, beta1, beta2, ftry, rstar, nucrate, rc)
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod
  use sulfate_utils
  
  implicit none
  
  type(carma_type), intent(in)         :: carma       !! the carma object
  type(carmastate_type), intent(inout) :: cstate      !! the carma state object
  integer, intent(in)                  :: iz          !! level index
  integer, intent(in)                  :: igroup      !! group index
  integer, intent(in)                  :: nucbin      !! bin in which nucleation occurs
  real(kind=f), intent(in)             :: h2o         !! H2O concentrations in molec/cm3 
  real(kind=f), intent(in)             :: h2so4       !! H2SO4 concentrations in molec/cm3
  real(kind=f), intent(in)             :: beta1
  real(kind=f), intent(in)             :: beta2
  real(kind=f), intent(in)             :: ftry
  real(kind=f), intent(in)             :: rstar       !! critical radius (cm)
  real(kind=f), intent(out)            :: nucrate     !! nucleation rate #/x/y/z/s
  integer, intent(inout)               :: rc          !! return code, negative indicates failure

  !  Local declarations     
  real(kind=f)                         :: cnucl
  real(kind=f)                         :: chom
  real(kind=f)                         :: expc
  real(kind=f)                         :: chet
  real(kind=f)                         :: xm
  real(kind=f)                         :: xm1
  real(kind=f)                         :: fxm
  real(kind=f)                         :: fv2
  real(kind=f)                         :: fu2
  real(kind=f)                         :: fv3
  real(kind=f)                         :: fv4
  real(kind=f)                         :: v1
  real(kind=f)                         :: fv1
  real(kind=f)                         :: ftry1
  real(kind=f)                         :: rarea
  real(kind=f)                         :: gg
  real(kind=f)                         :: FM = cos(50._f * DEG2RAD)  ! cos(contact angle)

  ! Heterogeneous nucleation which depends on r
  cnucl = 4._f * PI * rstar**(2._f)
  chom  = h2so4 * h2o * beta1 * cnucl
  expc  = 2.4e-16_f * exp(4.51872e+11_f / RGAS / t(iz))
  chet  = chom * expc * beta2     

  xm    = r(nucbin, igroup) / rstar                                                    
                 
  if (xm .lt. 1._f) then
    fxm = sqrt(1._f - 2._f * FM * xm + xm**(2._f))                                    
    fv2 = (xm - FM) / fxm                                                 
    fu2 = (1._f - xm * FM) / fxm                                   
    fv3 = (2._f + fv2) * xm**3._f  * (fv2 - 1._f)**(2._f)                                
    fv4 = 3._f * FM    * xm**2._f  * (fv2 - 1._f)                                        
  else 
    xm1 = 1._f / xm    
    fxm = sqrt(1._f - 2._f * FM * xm1 + xm1**2._f)                                    
    fu2 = (xm1 - FM) / fxm                                                
    fv2 = (1._f - xm1 * FM) / fxm                                   
    v1  = (FM**(2._f) - 1._f) / (fv2 + 1._f) / fxm**(2._f)
    fv3 = (2._f + fv2) * xm1 * v1**2._f                              
    fv4 = 3._f * FM * v1                                        
  endif

  fv1   = 0.5_f * (1._f + fu2**3._f + fv3 + fv4)                        

  ftry1 = ftry * fv1                                              
!  ftry1 = ftry * fh                                              
  if (ftry1 .lt. -1000._f) then
    nucrate = 0._f
  else
    
    rarea = 4._f * PI * r(nucbin, igroup)**2._f ! surface area per nucleus
    gg    = exp(ftry1)

    ! Calculate heterogeneous nucleation rate [embryos/s]
    ! NOTE: for [embryos/gridpoint/s], multipy rnuclg by pc [nuclei/gridpoint]
    nucrate = chet * gg * rarea  ! embryos/s
  end if
  
  return
end subroutine sulfhetnucrate
     

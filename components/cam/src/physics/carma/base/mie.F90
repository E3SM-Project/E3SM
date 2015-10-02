! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! There are several different algorithms that can be used to solve
!! a mie calculation for the optical properties of particles. This
!! routine provides a generic front end to these different mie
!! routines.
!!
!! Current methods are:
!!  miess - Original CARMA code, from Toon and Ackerman, supports core/shell
!!  bhmie - Homogeneous sphere, from Bohren and Huffman, handles wider range of parameters
!!
!! @author Chuck Bardeen
!! @version 2011
subroutine mie(carma, miertn, radius, wavelength, nmonomer, fractaldim, rmonomer, falpha_in, m, lqext, lqsca, lasym, rc)

  ! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carma_mod
  use fractal_meanfield_mod

  implicit none

  type(carma_type), intent(in)         :: carma         !! the carma object
  integer, intent(in)                  :: miertn        !! mie routine enumeration
  real(kind=f), intent(in)             :: radius        !! radius (cm)
  real(kind=f), intent(in)             :: wavelength    !! wavelength (cm)
  real(kind=f), intent(in)             :: nmonomer      !! number of monomers per aggregate [fractal particles only]
  real(kind=f), intent(in)             :: fractaldim    !! fractal dimension [fractal particles only]
  real(kind=f), intent(in)             :: rmonomer      !! monomer size (units?) [fractal particles only]
  real(kind=f), intent(in)             :: falpha_in        !! packing coefficient [fractal particles only]
  complex(kind=f), intent(in)          :: m             !! refractive index particle
  real(kind=f), intent(out)            :: lqext         !! EFFICIENCY FACTOR FOR EXTINCTION
  real(kind=f), intent(out)            :: lqsca         !! EFFICIENCY FACTOR FOR SCATTERING
  real(kind=f), intent(out)            :: lasym         !! asymmetry factor
  integer, intent(inout)               :: rc            !! return code, negative indicates failure
  

  integer, parameter                 :: nang     = 10   ! Number of angles
    
  real(kind=f)                       :: theta(IT)
  real(kind=f)                       :: wvno 
  real(kind=f)                       :: rfr 
  real(kind=f)                       :: rfi
  real(kind=f)                       :: x 
  real(kind=f)                       :: qback 
  real(kind=f)                       :: ctbrqs 
  complex(kind=f)                    :: s1(2*nang-1)
  complex(kind=f)                    :: s2(2*nang-1)
  real(kind=f)                       :: rmonomer_out      
  real(kind=f)                       :: fractaldim_out

  ! Calculate the wave number.
  wvno = 2._f * PI / wavelength
 
  ! Select the appropriate routine.
  if (miertn == I_MIERTN_TOON1981) then

    ! We only care about the forward direction.
    theta(:) = 0.0_f
    
    rfr = real(m)
    rfi = aimag(m)
    
    call miess(carma, &
               radius, &
               rfr, &
               rfi, &
               theta, &
               1, &
               lqext, &
               lqsca, &
               qback,&
               ctbrqs, &
               0.0_f, &
               rfr, &
               rfi, &
               wvno, &
               rc)
               
    lasym = ctbrqs / lqsca

  else if (miertn == I_MIERTN_BOHREN1983) then

    x = radius * wvno

    call bhmie(carma, &
               x, &
               m, &
               nang, &
               s1, &
               s2, &
               lqext, &
               lqsca, &
               qback, &
               lasym, &
               rc)

  else if (miertn == I_MIERTN_BOTET1997) then

    rfr = real(m)
    rfi = aimag(m)

    if (radius .le. rmonomer) then
      rmonomer_out = radius
      fractaldim_out = 3.0_f 
    else
      rmonomer_out = rmonomer
      fractaldim_out = fractaldim
    end if

    call fractal_meanfield(carma, &              !! carma object
                           wavelength*1.0e4_f, &   !! lambda in microns
                           rfi, &                !! imaginary index of refraction
                           rfr, &                !! real index of refraction
                           nmonomer, &           !! number of monomers
                           falpha_in, &          !! packing coefficient
                           fractaldim_out, &     !! fractal dimension
                           rmonomer_out, &       !! monomer size
                           1.0_f, &              !! xv,"set to 1"
                           0.0_f, &              !! angle, set to 0
                           lqext, &              !! extinction efficiency
                           lqsca, &              !! scattering efficiency
                           lasym, &              !! asymmetry parameter
                           rc)
                           
  else
    if (do_print) write(LUNOPRT, *) "mie::Unknown Mie routine specified."
    rc = RC_ERROR
  end if
  
  ! The mie code isn't perfect, so don't let it return values that aren't
  ! physical.
  lqext = max(lqext, 0._f)
  lqsca = max(0._f, min(lqext, lqsca))
  lasym = max(-1.0_f, min(1.0_f, lasym))
  
  return
end subroutine mie

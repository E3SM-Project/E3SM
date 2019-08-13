! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"
  

!!-----------------------------------------------------------------------
!! 
!! Purpose: Calculating the surface resistance, using the PBL parameter.
!!
!! Method: Zhang(2001), Atmospheric Environment
!!
!!
!! @author Tianyi Fan
!! @version Nov-2010
!!

subroutine calcrs(carma, cstate, ustar, tmp, radi, cc, vfall, rs, landidx, rc)
  use carma_precision_mod
  use carmastate_mod
  use carma_enums_mod
  use carma_types_mod
  use carma_mod
  use carma_constants_mod, only: BK, PI, GRAV
!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
!
! Arguments:
!
  type(carma_type), intent(in)          :: carma    !! the carma object
  type(carmastate_type), intent(in)     :: cstate   !! the carma state object
  real(kind=f), intent(in)              :: ustar    !! friction velocity [cm/s]
  real(kind=f), intent(in)              :: tmp      !! temperature [K]
  real(kind=f), intent(in)              :: radi     !! radius of the constitutent [cm]
  real(kind=f), intent(in)              :: cc       !! slip correction factor  
  real(kind=f), intent(in)              :: vfall    !! gravitational settling velocity, [cm/s]
  real(kind=f), intent(out)             :: rs       !! surface resistance [s/cm]
  integer,      intent(in)              :: landidx  !! landscape index, 1=land, 2=ocean, 3=sea ice
  integer, intent(inout)                :: rc       !! return code, negative indicates failure
  
  

! Local variables
  real(kind=f)            :: ebrn                   ! Brownian diffusion collection efficiency
  real(kind=f)            :: eimp                   ! Impaction collection efficiency
  real(kind=f)            :: eint                   ! Interception collection efficiency                 
  real(kind=f)            :: db                     ! Brownian diffusivity
  real(kind=f)            :: sc                     ! Schmidt number
  real(kind=f)            :: st                     ! Stokes number 
  real(kind=f)            :: rhoadry                ! dry air density [g/cm3] 
  real(kind=f)            :: eta                    ! kinematic viscosity of air [cm2/s] 
  real(kind=f), parameter :: xkar = 0.4_f           ! Von Karman's constant     
  real(kind=f), parameter :: eps0 = 3._f            ! empirical constant for rs, 3.0 in [Zhang, 2001], 1.0 in [Seinfeld and Pandis]

  ! exponent in the eb dependence on sc, 2/3 in [Seinfeld and Pandis, 1998], 1/2 in [Lewis and Schwartz, 2004]
  real(kind=f)            :: lam

  integer                 :: ibot
  
  if (igridv .eq. I_CART) then
    ibot = 1
  else
    ibot = NZ
  end if   
  
  ! Unit conversion   
  rhoadry = rhoa(ibot) / zmet(ibot) / xmet(ibot) / ymet(ibot)   ! [g/cm3] 
  eta = rmu(ibot) / rhoadry                                     ! rmu, aerodynamic viscosity of air [g/cm/s]

  if (landidx .eq. 1)  then
    lam = 2._f / 3._f
  else
    lam = 1._f / 2._f
  end if   
   
  ! Surface Resistance = Brownian + Impaction + Interception
  
  ! ** Brownian diffusion
  db = BK  * tmp * cc / (6._f * PI * rmu(ibot)  * radi)   ! [cm2/s]
    
  sc = eta / db    ! [-]
  ebrn = sc**(-lam)
 
  ! ** Impaction 
  st = vfall * ustar**2 / (GRAV * eta)   ! [-]
   
  ! [Slinn, 1982]   
  ! eimp = 10. ** (-3._f/st)     

  ! [Peters and Eiden, 1992]
  eimp = (st / (0.8_f + st))**2
!  eimp = max(eimp, 1.e-10_f)
  
  ! ** Interception
  !
  ! NOTE: Interception is not currently considered for ocean and ice.
  if (landidx .eq. 1) then
!    eint = 0.3_f * (0.01_f * radi * 1.e-2_f / (radi * 1.e-2_f + 1.e-5_f) + 0.99_f *radi*1.e-2_f / (radi*1.e-2_f + 8.e-4_f))
    eint = 0.3_f * (0.01_f * radi / (radi + 1.e-3_f) + 0.99_f * radi / (radi + 8.e-2_f))
  else
    eint = 0._f
  end if
  
  if (ustar > 0._f) then       
    rs = 1._f / (eps0 * ustar * (ebrn + eimp + eint ))   ! [s/cm]
  else
    rs = 0._f
  end if
     
  return
end subroutine calcrs

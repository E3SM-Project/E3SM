
module mo_constants

  use shr_kind_mod,  only : r8 => shr_kind_r8
  use physconst,      only : gravit,&
                             rearth,&
                             avogad,&
                             boltz,&
                             pi,&
                             rgrav=>rga

  implicit none

  save

  real(r8), parameter ::  dayspy   = 365._r8             ! days per year
  real(r8), parameter ::  avogadro = avogad*1.e-3_r8  ! Avogadro numb - molecules/mole
  real(r8), parameter ::  boltzmann  = boltz*1.e7_r8 ! erg/k

  real(r8), parameter ::  d2r  = pi/180._r8               ! radians to degrees
  real(r8), parameter ::  r2d  = 180._r8/pi               ! degrees to radians 

contains

  subroutine mo_constants_inti

    implicit none

  end subroutine mo_constants_inti

end module mo_constants

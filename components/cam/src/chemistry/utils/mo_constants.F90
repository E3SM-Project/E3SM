
module mo_constants

  use shr_kind_mod,  only : r8 => shr_kind_r8
  use physconst,      only : gravity => gravit,&
                             rearth,&
                             avogadro_kmol => avogad,&
                             rgas_kmol => r_universal,&
                             boltz,&
                             pi,&
                             rgrav=>rga,&
                             rhoh2o

  implicit none

  save

  real(r8), parameter ::  dayspy   = 365._r8                ! days per year
  ! The following are converted from kmol to mol.
  real(r8), parameter ::  avogadro = avogadro_kmol*1.e-3_r8 ! Avogadro numb - molecules/mole
  real(r8), parameter ::  rgas = rgas_kmol*1.e-3_r8         ! Gas constant (J/K/mol)

  ! Call out cgs units via explicit naming.
  real(r8), parameter ::  rgas_cgs    = rgas*1.e7_r8        ! erg/K/mol
  real(r8), parameter ::  boltz_cgs   = boltz*1.e7_r8       ! erg/K
  real(r8), parameter ::  rhoh2o_cgs  = rhoh2o*1.e-3_r8     ! g/ml

  ! Not a compile-time constant.
  real(r8) :: gravity_cgs = huge(1._r8)                     ! cm/s^2

  real(r8), parameter ::  d2r  = pi/180._r8                 ! radians to degrees
  real(r8), parameter ::  r2d  = 180._r8/pi                 ! degrees to radians 

  real(r8), parameter :: seasalt_density = 2.2e+3_r8        ! [kg m-3] Aerosol density
  real(r8), parameter :: dust_density    = 2.5e+3_r8        ! [kg m-3] Aerosol density
contains

  subroutine mo_constants_inti

    gravity_cgs = gravity*1.e2_r8

  end subroutine mo_constants_inti

end module mo_constants

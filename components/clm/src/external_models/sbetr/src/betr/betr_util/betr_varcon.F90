module betr_varcon
  use bshr_kind_mod  , only : r8 => shr_kind_r8
  use bshr_const_mod , only : SHR_CONST_G,SHR_CONST_STEBOL,SHR_CONST_KARMAN,  &
                           SHR_CONST_RWV,SHR_CONST_RDAIR,SHR_CONST_CPFW,      &
                           SHR_CONST_CPICE,SHR_CONST_CPDAIR,SHR_CONST_LATVAP, &
                           SHR_CONST_LATSUB,SHR_CONST_LATICE,SHR_CONST_RHOFW, &
                           SHR_CONST_RHOICE,SHR_CONST_TKFRZ,SHR_CONST_REARTH, &
                           SHR_CONST_PDB, SHR_CONST_PI, SHR_CONST_CDAY,       &
                           SHR_CONST_RGAS
  implicit none
  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  character(len=16), parameter :: bnamec  = 'column'       ! name of columns
  character(len=16), parameter :: bnamep  = 'pft'          ! name of patches
  real(r8), public,  parameter :: bspval = 1.e36_r8  ! special value for real data
  integer , public,  parameter :: bispval = -9999     ! special value for int data (keep this negative to avoid conflicts with possible valid values)

  real(r8) :: bc14ratio = 1.e-12_r8

  real(r8) :: boneatm = 1.01325e5_r8 !one standard atmospheric pressure
  real(r8) :: brpi    = SHR_CONST_PI
  real(r8) :: bgrav   = SHR_CONST_G      !gravity constant [m/s2]
  real(r8) :: bsb     = SHR_CONST_STEBOL !stefan-boltzmann constant  [W/m2/K4]
  real(r8) :: bvkc    = SHR_CONST_KARMAN !von Karman constant [-]
  real(r8) :: brwat   = SHR_CONST_RWV    !gas constant for water vapor [J/(kg K)]
  real(r8) :: brair   = SHR_CONST_RDAIR  !gas constant for dry air [J/kg/K]
  real(r8) :: broverg = SHR_CONST_RWV/SHR_CONST_G*1000._r8 !Rw/g constant = (8.3144/0.018)/(9.80616)*1000. mm/K
  real(r8) :: bcpliq  = SHR_CONST_CPFW   !Specific heat of water [J/kg-K]
  real(r8) :: bcpice  = SHR_CONST_CPICE  !Specific heat of ice [J/kg-K]
  real(r8) :: bcpair  = SHR_CONST_CPDAIR !specific heat of dry air [J/kg/K]
  real(r8) :: bhvap   = SHR_CONST_LATVAP !Latent heat of evap for water [J/kg]
  real(r8) :: bhsub   = SHR_CONST_LATSUB !Latent heat of sublimation    [J/kg]
  real(r8) :: bhfus   = SHR_CONST_LATICE !Latent heat of fusion for ice [J/kg]
  real(r8) :: bdenh2o = SHR_CONST_RHOFW  !density of liquid water [kg/m3]
  real(r8) :: bdenice = SHR_CONST_RHOICE !density of ice [kg/m3]
  real(r8) :: brgas   = SHR_CONST_RGAS   !universal gas constant [J/K/kmole]
  real(r8) :: btkair  = 0.023_r8     !thermal conductivity of air   [W/m/K]
  real(r8) :: btkice  = 2.290_r8     !thermal conductivity of ice   [W/m/K]
  real(r8) :: btkwat  = 0.57_r8       !thermal conductivity of water [W/m/K]
  real(r8) :: btfrz   = SHR_CONST_TKFRZ  !freezing temperature [K]
  real(r8), public, parameter ::  bsecspday= SHR_CONST_CDAY  ! Seconds per day

  integer, public :: betr_maxpatch_pft = 1
  integer, public :: betr_max_soilorder = 1
  integer, public :: bspinup_state = 0
end module betr_varcon

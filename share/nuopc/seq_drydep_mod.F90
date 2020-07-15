module seq_drydep_mod

  !========================================================================
  ! Module for handling dry depostion of tracers.
  ! This module is shared by land and atmosphere models for the computations of
  ! dry deposition of tracers
  !========================================================================

  use ESMF           , only : ESMF_VMGetCurrent, ESMF_VM, ESMF_VMGet
  use ESMF           , only : ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU, ESMF_SUCCESS
  use shr_sys_mod    , only : shr_sys_abort
  use shr_kind_mod   , only : r8 => shr_kind_r8, CS => SHR_KIND_CS, CX => SHR_KIND_CX
  use shr_const_mod  , only : SHR_CONST_G, SHR_CONST_RDAIR, SHR_CONST_CPDAIR, SHR_CONST_MWWV
  use shr_mpi_mod    , only : shr_mpi_bcast
  use shr_nl_mod     , only : shr_nl_find_group_name
  use shr_log_mod    , only : s_logunit => shr_log_Unit
  use shr_infnan_mod , only : shr_infnan_posinf, assignment(=)

  implicit none
  private

  ! public member functions
  public :: seq_drydep_readnl       ! Read namelist
  public :: seq_drydep_init         ! Initialization of drydep data
  public :: seq_drydep_setHCoeff    ! Calculate Henry's law coefficients

  ! private array sizes
  integer, public,  parameter :: n_species_table = 192     ! Number of species to work with
  integer, private, parameter :: maxspc = 210              ! Maximum number of species
  integer, private, parameter :: NSeas = 5                 ! Number of seasons
  integer, private, parameter :: NLUse = 11                ! Number of land-use types

  logical, private :: drydep_initialized = .false. 

  ! public data members:
  ! method specification
  character(16),public,parameter :: DD_XATM = 'xactive_atm' ! dry-dep atmosphere
  character(16),public,parameter :: DD_XLND = 'xactive_lnd' ! dry-dep land
  character(16),public,parameter :: DD_TABL = 'table'       ! dry-dep table (atm and lnd)
  character(16),public           :: drydep_method = DD_XLND ! Which option choosen

  real(r8), public, parameter :: ph     = 1.e-5_r8         ! measure of the acidity (dimensionless)

  logical, public  :: lnd_drydep                           ! If dry-dep fields passed
  integer, public  :: n_drydep = 0                         ! Number in drypdep list
  logical          :: drydep_init = .false.                ! has seq_drydep_init been called?
  character(len=CS), public, dimension(maxspc) :: drydep_list = ''   ! List of dry-dep species

  real(r8), public, allocatable, dimension(:) :: foxd      ! reactivity factor for oxidation (dimensioness)
  real(r8), public, allocatable, dimension(:) :: drat      ! ratio of molecular diffusivity (D_H2O/D_species; dimensionless)
  integer,  public, allocatable, dimension(:) :: mapping   ! mapping to species table

  ! --- Indices for each species ---
  integer,  public :: h2_ndx, ch4_ndx, co_ndx, pan_ndx, mpan_ndx, so2_ndx, o3_ndx, o3a_ndx, xpan_ndx

  !---------------------------------------------------------------------------
  ! Table 1 from Wesely, Atmos. Environment, 1989, p1293
  ! Table 2 from Sheih, microfiche PB86-218104 and Walcek, Atmos.  Environment, 1986, p949
  ! Table 3-5 compiled by P. Hess
  !
  ! index #1 : season
  !           1 -> midsummer with lush vegetation
  !           2 -> autumn with unharvested cropland
  !           3 -> late autumn after frost, no snow
  !           4 -> winter, snow on ground, and subfreezing
  !           5 -> transitional spring with partially green short annuals
  !
  ! index #2 : landuse type
  !           1 -> urban land
  !           2 -> agricultural land
  !           3 -> range land
  !           4 -> deciduous forest
  !           5 -> coniferous forest
  !           6 -> mixed forest including wetland
  !           7 -> water, both salt and fresh
  !           8 -> barren land, mostly desert
  !           9 -> nonforested wetland
  !           10 -> mixed agricultural and range land
  !           11 -> rocky open areas with low growing shrubs
  !
  ! JFL August 2000
  !---------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  ! table to parameterize the impact of soil moisture on the deposition of H2 and
  ! CO on soils (from Sanderson et al., J. Atmos. Chem., 46, 15-28, 2003).
  !---------------------------------------------------------------------------

  !--- deposition of h2 and CO on soils ---
  real(r8), parameter, public :: h2_a(NLUse) = &
                (/  0.000_r8,  0.000_r8, 0.270_r8,  0.000_r8,  0.000_r8,  &
                    0.000_r8,  0.000_r8, 0.000_r8,  0.000_r8,  0.000_r8, 0.000_r8/)
  !--- deposition of h2 and CO on soils ---
  real(r8), parameter, public :: h2_b(NLUse) = &
                (/  0.000_r8,-41.390_r8, -0.472_r8,-41.900_r8,-41.900_r8,  &
                  -41.900_r8,  0.000_r8,  0.000_r8,  0.000_r8,-41.390_r8,  0.000_r8/)
  !--- deposition of h2 and CO on soils ---
  real(r8), parameter, public :: h2_c(NLUse) = &
                (/  0.000_r8, 16.850_r8, 1.235_r8, 19.700_r8, 19.700_r8, &
                   19.700_r8,  0.000_r8, 0.000_r8,  0.000_r8, 17.700_r8, 1.000_r8/)

  !--- deposition of h2 and CO on soils
  !
  !--- ri:   Richardson number                      (dimensionless)
  !--- rlu:  Resistance of leaves in upper canopy   (s.m-1)
  !--- rac:  Aerodynamic resistance to lower canopy (s.m-1)
  !--- rgss: Ground surface resistance for SO2      (s.m-1)
  !--- rgso: Ground surface resistance for O3       (s.m-1)
  !--- rcls: Lower canopy resistance for SO2        (s.m-1)
  !--- rclo: Lower canopy resistance for O3         (s.m-1)
  !
  real(r8), public, dimension(NSeas,NLUse) :: ri, rlu, rac, rgss, rgso, rcls, rclo

  data ri  (1,1:NLUse) &
       /1.e36_r8,  60._r8, 120._r8,  70._r8, 130._r8, 100._r8,1.e36_r8,1.e36_r8,  80._r8, 100._r8, 150._r8/
  data rlu (1,1:NLUse) &
       /1.e36_r8,2000._r8,2000._r8,2000._r8,2000._r8,2000._r8,1.e36_r8,1.e36_r8,2500._r8,2000._r8,4000._r8/
  data rac (1,1:NLUse) &
       / 100._r8, 200._r8, 100._r8,2000._r8,2000._r8,2000._r8,   0._r8,   0._r8, 300._r8, 150._r8, 200._r8/
  data rgss(1,1:NLUse) &
       / 400._r8, 150._r8, 350._r8, 500._r8, 500._r8, 100._r8,   0._r8,1000._r8,  0._r8, 220._r8, 400._r8/
  data rgso(1,1:NLUse) &
       / 300._r8, 150._r8, 200._r8, 200._r8, 200._r8, 300._r8,2000._r8, 400._r8,1000._r8, 180._r8, 200._r8/
  data rcls(1,1:NLUse) &
       /1.e36_r8,2000._r8,2000._r8,2000._r8,2000._r8,2000._r8,1.e36_r8,1.e36_r8,2500._r8,2000._r8,4000._r8/
  data rclo(1,1:NLUse) &
       /1.e36_r8,1000._r8,1000._r8,1000._r8,1000._r8,1000._r8,1.e36_r8,1.e36_r8,1000._r8,1000._r8,1000._r8/

  data ri  (2,1:NLUse) &
       /1.e36_r8,1.e36_r8,1.e36_r8,1.e36_r8, 250._r8, 500._r8,1.e36_r8,1.e36_r8,1.e36_r8,1.e36_r8,1.e36_r8/
  data rlu (2,1:NLUse) &
       /1.e36_r8,9000._r8,9000._r8,9000._r8,4000._r8,8000._r8,1.e36_r8,1.e36_r8,9000._r8,9000._r8,9000._r8/
  data rac (2,1:NLUse) &
       / 100._r8, 150._r8, 100._r8,1500._r8,2000._r8,1700._r8,   0._r8,   0._r8, 200._r8, 120._r8, 140._r8/
  data rgss(2,1:NLUse) &
       / 400._r8, 200._r8, 350._r8, 500._r8, 500._r8, 100._r8,   0._r8,1000._r8,   0._r8, 300._r8, 400._r8/
  data rgso(2,1:NLUse) &
       / 300._r8, 150._r8, 200._r8, 200._r8, 200._r8, 300._r8,2000._r8, 400._r8, 800._r8, 180._r8, 200._r8/
  data rcls(2,1:NLUse) &
       /1.e36_r8,9000._r8,9000._r8,9000._r8,2000._r8,4000._r8,1.e36_r8,1.e36_r8,9000._r8,9000._r8,9000._r8/
  data rclo(2,1:NLUse) &
       /1.e36_r8, 400._r8, 400._r8, 400._r8,1000._r8, 600._r8,1.e36_r8,1.e36_r8, 400._r8, 400._r8, 400._r8/

  data ri  (3,1:NLUse) &
       /1.e36_r8,1.e36_r8,1.e36_r8,1.e36_r8, 250._r8, 500._r8,1.e36_r8,1.e36_r8,1.e36_r8,1.e36_r8,1.e36_r8/
  data rlu (3,1:NLUse) &
       /1.e36_r8,1.e36_r8,9000._r8,9000._r8,4000._r8,8000._r8,1.e36_r8,1.e36_r8,9000._r8,9000._r8,9000._r8/
  data rac (3,1:NLUse) &
       / 100._r8,  10._r8, 100._r8,1000._r8,2000._r8,1500._r8,   0._r8,   0._r8, 100._r8, 50._r8, 120._r8/
  data rgss(3,1:NLUse) &
       / 400._r8, 150._r8, 350._r8, 500._r8, 500._r8, 200._r8,   0._r8,1000._r8,   0._r8, 200._r8, 400._r8/
  data rgso(3,1:NLUse) &
       / 300._r8, 150._r8, 200._r8, 200._r8, 200._r8, 300._r8,2000._r8, 400._r8,1000._r8, 180._r8, 200._r8/
  data rcls(3,1:NLUse) &
       /1.e36_r8,1.e36_r8,9000._r8,9000._r8,3000._r8,6000._r8,1.e36_r8,1.e36_r8,9000._r8,9000._r8,9000._r8/
  data rclo(3,1:NLUse) &
       /1.e36_r8,1000._r8, 400._r8, 400._r8,1000._r8, 600._r8,1.e36_r8,1.e36_r8, 800._r8, 600._r8, 600._r8/

  data ri  (4,1:NLUse) &
       /1.e36_r8,1.e36_r8,1.e36_r8,1.e36_r8, 400._r8, 800._r8,1.e36_r8,1.e36_r8,1.e36_r8,1.e36_r8,1.e36_r8/
  data rlu (4,1:NLUse) &
       /1.e36_r8,1.e36_r8,1.e36_r8,1.e36_r8,6000._r8,9000._r8,1.e36_r8,1.e36_r8,9000._r8,9000._r8,9000._r8/
  data rac (4,1:NLUse) &
       / 100._r8,  10._r8,  10._r8,1000._r8,2000._r8,1500._r8,   0._r8,   0._r8,  50._r8,  10._r8,  50._r8/
  data rgss(4,1:NLUse) &
       / 100._r8, 100._r8, 100._r8, 100._r8, 100._r8, 100._r8,   0._r8,1000._r8, 100._r8, 100._r8,  50._r8/
  data rgso(4,1:NLUse) &
       / 600._r8,3500._r8,3500._r8,3500._r8,3500._r8,3500._r8,2000._r8, 400._r8,3500._r8,3500._r8,3500._r8/
  data rcls(4,1:NLUse) &
       /1.e36_r8,1.e36_r8,1.e36_r8,9000._r8, 200._r8, 400._r8,1.e36_r8,1.e36_r8,9000._r8,1.e36_r8,9000._r8/
  data rclo(4,1:NLUse) &
       /1.e36_r8,1000._r8,1000._r8, 400._r8,1500._r8, 600._r8,1.e36_r8,1.e36_r8, 800._r8,1000._r8, 800._r8/

  data ri  (5,1:NLUse) &
       /1.e36_r8, 120._r8, 240._r8, 140._r8, 250._r8, 190._r8,1.e36_r8,1.e36_r8, 160._r8, 200._r8, 300._r8/
  data rlu (5,1:NLUse) &
       /1.e36_r8,4000._r8,4000._r8,4000._r8,2000._r8,3000._r8,1.e36_r8,1.e36_r8,4000._r8,4000._r8,8000._r8/
  data rac (5,1:NLUse) &
       / 100._r8,  50._r8,  80._r8,1200._r8,2000._r8,1500._r8,   0._r8,   0._r8, 200._r8, 60._r8, 120._r8/
  data rgss(5,1:NLUse) &
       / 500._r8, 150._r8, 350._r8, 500._r8, 500._r8, 200._r8,   0._r8,1000._r8,   0._r8, 250._r8, 400._r8/
  data rgso(5,1:NLUse) &
       / 300._r8, 150._r8, 200._r8, 200._r8, 200._r8, 300._r8,2000._r8, 400._r8,1000._r8, 180._r8, 200._r8/
  data rcls(5,1:NLUse) &
       /1.e36_r8,4000._r8,4000._r8,4000._r8,2000._r8,3000._r8,1.e36_r8,1.e36_r8,4000._r8,4000._r8,8000._r8/
  data rclo(5,1:NLUse) &
       /1.e36_r8,1000._r8, 500._r8, 500._r8,1500._r8, 700._r8,1.e36_r8,1.e36_r8, 600._r8, 800._r8, 800._r8/

  !---------------------------------------------------------------------------
  !         ... roughness length
  !---------------------------------------------------------------------------
  real(r8), public, dimension(NSeas,NLUse) :: z0

  data z0  (1,1:NLUse) &
       /1.000_r8,0.250_r8,0.050_r8,1.000_r8,1.000_r8,1.000_r8,0.0006_r8,0.002_r8,0.150_r8,0.100_r8,0.100_r8/
  data z0  (2,1:NLUse) &
       /1.000_r8,0.100_r8,0.050_r8,1.000_r8,1.000_r8,1.000_r8,0.0006_r8,0.002_r8,0.100_r8,0.080_r8,0.080_r8/
  data z0  (3,1:NLUse) &
       /1.000_r8,0.005_r8,0.050_r8,1.000_r8,1.000_r8,1.000_r8,0.0006_r8,0.002_r8,0.100_r8,0.020_r8,0.060_r8/
  data z0  (4,1:NLUse) &
       /1.000_r8,0.001_r8,0.001_r8,1.000_r8,1.000_r8,1.000_r8,0.0006_r8,0.002_r8,0.001_r8,0.001_r8,0.040_r8/
  data z0  (5,1:NLUse) &
       /1.000_r8,0.030_r8,0.020_r8,1.000_r8,1.000_r8,1.000_r8,0.0006_r8,0.002_r8,0.010_r8,0.030_r8,0.060_r8/

  !real(r8), private, dimension(11,5), parameter :: z0xxx = reshape ( &
  ! (/   1.000,0.250,0.050,1.000,1.000,1.000,0.0006,0.002,0.150,0.100,0.100 ,  &
  !      1.000,0.100,0.050,1.000,1.000,1.000,0.0006,0.002,0.100,0.080,0.080 ,  &
  !      1.000,0.005,0.050,1.000,1.000,1.000,0.0006,0.002,0.100,0.020,0.060 ,  &
  !      1.000,0.001,0.001,1.000,1.000,1.000,0.0006,0.002,0.001,0.001,0.040 ,  &
  !      1.000,0.030,0.020,1.000,1.000,1.000,0.0006,0.002,0.010,0.030,0.060  /), (/11,5/) )

  !---------------------------------------------------------------------------
  ! public chemical data
  !---------------------------------------------------------------------------

  !--- data for foxd (reactivity factor for oxidation) ----
  real(r8), public, parameter :: dfoxd(n_species_table) = &
          (/  1._r8     & ! OX
             ,1._r8     & ! H2O2
             ,1._r8     & ! OH
             ,.1_r8     & ! HO2
             ,1.e-36_r8 & ! CO
             ,1.e-36_r8 & ! CH4
             ,1._r8     & ! CH3O2
             ,1._r8     & ! CH3OOH
             ,1._r8     & ! CH2O
             ,1._r8     & ! HCOOH
             ,0._r8     & ! NO
             ,.1_r8     & ! NO2
             ,1.e-36_r8 & ! HNO3
             ,1.e-36_r8 & ! CO2
             ,1.e-36_r8 & ! NH3
             ,.1_r8     & ! N2O5
             ,1._r8     & ! NO3
             ,1._r8     & ! CH3OH
             ,.1_r8     & ! HO2NO2
             ,1._r8     & ! O1D
             ,1.e-36_r8 & ! C2H6
             ,.1_r8     & ! C2H5O2
             ,.1_r8     & ! PO2
             ,.1_r8     & ! MACRO2
             ,.1_r8     & ! ISOPO2
             ,1.e-36_r8 & ! C4H10
             ,1._r8     & ! CH3CHO
             ,1._r8     & ! C2H5OOH
             ,1.e-36_r8 & ! C3H6
             ,1._r8     & ! POOH
             ,1.e-36_r8 & ! C2H4
             ,.1_r8     & ! PAN
             ,1._r8     & ! CH3COOOH
             ,1.e-36_r8 & ! MTERP
             ,1._r8     & ! GLYOXAL
             ,1._r8     & ! CH3COCHO
             ,1._r8     & ! GLYALD
             ,.1_r8     & ! CH3CO3
             ,1.e-36_r8 & ! C3H8
             ,.1_r8     & ! C3H7O2
             ,1._r8     & ! CH3COCH3
             ,1._r8     & ! C3H7OOH
             ,.1_r8     & ! RO2
             ,1._r8     & ! ROOH
             ,1.e-36_r8 & ! Rn
             ,1.e-36_r8 & ! ISOP
             ,1._r8     & ! MVK
             ,1._r8     & ! MACR
             ,1._r8     & ! C2H5OH 
             ,1._r8     & ! ONITR
             ,.1_r8     & ! ONIT
             ,.1_r8     & ! ISOPNO3
             ,1._r8     & ! HYDRALD
             ,1.e-36_r8 & ! HCN
             ,1.e-36_r8 & ! CH3CN
             ,1.e-36_r8 & ! SO2
             ,0.1_r8    & ! SOAGff0
             ,0.1_r8    & ! SOAGff1
             ,0.1_r8    & ! SOAGff2
             ,0.1_r8    & ! SOAGff3
             ,0.1_r8    & ! SOAGff4
             ,0.1_r8    & ! SOAGbg0
             ,0.1_r8    & ! SOAGbg1
             ,0.1_r8    & ! SOAGbg2
             ,0.1_r8    & ! SOAGbg3
             ,0.1_r8    & ! SOAGbg4
             ,0.1_r8    & ! SOAG0
             ,0.1_r8    & ! SOAG1
             ,0.1_r8    & ! SOAG2
             ,0.1_r8    & ! SOAG3
             ,0.1_r8    & ! SOAG4
             ,0.1_r8    & ! IVOC
             ,0.1_r8    & ! SVOC 
             ,0.1_r8    & ! IVOCbb
             ,0.1_r8    & ! IVOCff
             ,0.1_r8    & ! SVOCbb
             ,0.1_r8    & ! SVOCff
             ,1.e-36_r8 & ! N2O
             ,1.e-36_r8 & ! H2
             ,1.e-36_r8 & ! C2H2
             ,1._r8     & ! CH3COOH
             ,1._r8     & ! EOOH
             ,1._r8     & ! HYAC
             ,1.e-36_r8 & ! BIGENE
             ,1.e-36_r8 & ! BIGALK
             ,1._r8     & ! MEK
             ,1._r8     & ! MEKOOH
             ,1._r8     & ! MACROOH
             ,1._r8     & ! MPAN
             ,1._r8     & ! ALKNIT
             ,1._r8     & ! NOA
             ,1._r8     & ! ISOPNITA
             ,1._r8     & ! ISOPNITB
             ,1._r8     & ! ISOPNOOH
             ,1._r8     & ! NC4CHO
             ,1._r8     & ! NC4CH2OH
             ,1._r8     & ! TERPNIT
             ,1._r8     & ! NTERPOOH
             ,1._r8     & ! ALKOOH
             ,1._r8     & ! BIGALD
             ,1._r8     & ! HPALD
             ,1._r8     & ! IEPOX
             ,1._r8     & ! XOOH
             ,1._r8     & ! ISOPOOH
             ,1.e-36_r8 & ! TOLUENE
             ,1._r8     & ! CRESOL
             ,1._r8     & ! TOLOOH
             ,1.e-36_r8 & ! BENZENE
             ,1._r8     & ! PHENOL
             ,1._r8     & ! BEPOMUC
             ,1._r8     & ! PHENOOH
             ,1._r8     & ! C6H5OOH
             ,1._r8     & ! BENZOOH
             ,1._r8     & ! BIGALD1
             ,1._r8     & ! BIGALD2
             ,1._r8     & ! BIGALD3
             ,1._r8     & ! BIGALD4
             ,1._r8     & ! TEPOMUC
             ,1._r8     & ! BZOOH
             ,1._r8     & ! BZALD
             ,1._r8     & ! PBZNIT
             ,1.e-36_r8 & ! XYLENES
             ,1._r8     & ! XYLOL
             ,1._r8     & ! XYLOLOOH
             ,1._r8     & ! XYLENOOH
             ,1.e-36_r8 & ! BCARY
             ,1._r8     & ! TERPOOH
             ,1._r8     & ! TERPROD1
             ,1._r8     & ! TERPROD2
             ,1._r8     & ! TERP2OOH
             ,1.e-36_r8 & ! DMS
             ,1.e-36_r8 & ! H2SO4
             ,1._r8     & ! HONITR
             ,1._r8     & ! MACRN   
             ,1._r8     & ! MVKN
             ,1._r8     & ! ISOPN2B
             ,1._r8     & ! ISOPN3B
             ,1._r8     & ! ISOPN4D
             ,1._r8     & ! ISOPN1D
             ,1._r8     & ! ISOPNOOHD
             ,1._r8     & ! ISOPNOOHB
             ,1._r8     & ! ISOPNBNO3
             ,1._r8     & ! NO3CH2CHO
             ,1._r8     & ! HYPERACET
             ,1._r8     & ! HCOCH2OOH
             ,1._r8     & ! DHPMPAL
             ,1._r8     & ! MVKOOH
             ,1._r8     & ! ISOPOH
             ,1._r8     & ! ISOPFDN
             ,1._r8     & ! ISOPFNP
             ,1._r8     & ! INHEB
             ,1._r8     & ! HMHP
             ,1._r8     & ! HPALD1
             ,1._r8     & ! INHED
             ,1._r8     & ! HPALD4  
             ,1._r8     & ! ISOPHFP
             ,1._r8     & ! HPALDB1C
             ,1._r8     & ! HPALDB4C
             ,1._r8     & ! ICHE
             ,1._r8     & ! ISOPFDNC
             ,1._r8     & ! ISOPFNC
             ,1._r8     & ! TERPNT
             ,1._r8     & ! TERPNS
             ,1._r8     & ! TERPNT1
             ,1._r8     & ! TERPNS1
             ,1._r8     & ! TERPNPT
             ,1._r8     & ! TERPNPS
             ,1._r8     & ! TERPNPT1
             ,1._r8     & ! TERPNPS1
             ,1._r8     & ! TERPFDN
             ,1._r8     & ! SQTN
             ,1._r8     & ! TERPHFN 
             ,1._r8     & ! TERP1OOH
             ,1._r8     & ! TERPDHDP
             ,1._r8     & ! TERPF2
             ,1._r8     & ! TERPF1
             ,1._r8     & ! TERPA
             ,1._r8     & ! TERPA2 
             ,1._r8     & ! TERPK
             ,1._r8     & ! TERPAPAN
             ,1._r8     & ! TERPACID
             ,1._r8     & ! TERPA2PAN
             ,1.e-36_r8 & ! APIN
             ,1.e-36_r8 & ! BPIN
             ,1.e-36_r8 & ! LIMON
             ,1.e-36_r8 & ! MYRC
             ,1._r8     & ! TERPACID2
             ,1._r8     & ! TERPACID3
             ,1._r8     & ! TERPA3PAN
             ,1._r8     & ! TERPOOHL
             ,1._r8     & ! TERPA3
             ,1._r8     & ! TERP2AOOH
            /)

  ! PRIVATE DATA:

  Interface seq_drydep_setHCoeff                      ! overload subroutine
     Module Procedure set_hcoeff_scalar
     Module Procedure set_hcoeff_vector
  End Interface

  real(r8), private, parameter :: small_value = 1.e-36_r8          !--- smallest value to use ---

  !---------------------------------------------------------------------------
  ! private chemical data
  !---------------------------------------------------------------------------

  !--- Names of species that can work with ---
  character(len=20), public, parameter :: species_name_table(n_species_table) = &
       (/ 'OX       ' &
         ,'H2O2     ' &
         ,'OH       ' &
         ,'HO2      ' &
         ,'CO       ' &
         ,'CH4      ' &
         ,'CH3O2    ' &
         ,'CH3OOH   ' &
         ,'CH2O     ' &
         ,'HCOOH    ' &
         ,'NO       ' &
         ,'NO2      ' &
         ,'HNO3     ' &
         ,'CO2      ' &
         ,'NH3      ' &
         ,'N2O5     ' &
         ,'NO3      ' &
         ,'CH3OH    ' &
         ,'HO2NO2   ' &
         ,'O1D      ' &
         ,'C2H6     ' &
         ,'C2H5O2   ' &
         ,'PO2      ' &
         ,'MACRO2   ' &
         ,'ISOPO2   ' &
         ,'C4H10    ' &
         ,'CH3CHO   ' &
         ,'C2H5OOH  ' &
         ,'C3H6     ' &
         ,'POOH     ' &
         ,'C2H4     ' &
         ,'PAN      ' &
         ,'CH3COOOH ' &
         ,'MTERP    ' &
         ,'GLYOXAL  ' &
         ,'CH3COCHO ' &
         ,'GLYALD   ' &
         ,'CH3CO3   ' &
         ,'C3H8     ' &
         ,'C3H7O2   ' &
         ,'CH3COCH3 ' &
         ,'C3H7OOH  ' &
         ,'RO2      ' &
         ,'ROOH     ' &
         ,'Rn       ' &
         ,'ISOP     ' &
         ,'MVK      ' &
         ,'MACR     ' &
         ,'C2H5OH   ' &
         ,'ONITR    ' &
         ,'ONIT     ' &
         ,'ISOPNO3  ' &
         ,'HYDRALD  ' &
         ,'HCN      ' &
         ,'CH3CN    ' &
         ,'SO2      ' &
         ,'SOAGff0  ' &
         ,'SOAGff1  ' &
         ,'SOAGff2  ' &
         ,'SOAGff3  ' &
         ,'SOAGff4  ' &
         ,'SOAGbg0  ' &
         ,'SOAGbg1  ' &
         ,'SOAGbg2  ' &
         ,'SOAGbg3  ' &
         ,'SOAGbg4  ' &
         ,'SOAG0    ' &
         ,'SOAG1    ' &
         ,'SOAG2    ' &
         ,'SOAG3    ' &
         ,'SOAG4    ' &
         ,'IVOC     ' &
         ,'SVOC     ' &
         ,'IVOCbb   ' &
         ,'IVOCff   ' &
         ,'SVOCbb   ' &
         ,'SVOCff   ' &
         ,'N2O      ' &
         ,'H2       ' &
         ,'C2H2     ' &
         ,'CH3COOH  ' &
         ,'EOOH     ' &
         ,'HYAC     ' &
         ,'BIGENE   ' &
         ,'BIGALK   ' &
         ,'MEK      ' &
         ,'MEKOOH   ' &
         ,'MACROOH  ' &
         ,'MPAN     ' &
         ,'ALKNIT   ' &
         ,'NOA      ' &
         ,'ISOPNITA ' &
         ,'ISOPNITB ' &
         ,'ISOPNOOH ' &
         ,'NC4CHO   ' &
         ,'NC4CH2OH ' &
         ,'TERPNIT  ' &
         ,'NTERPOOH ' &
         ,'ALKOOH   ' &
         ,'BIGALD   ' &
         ,'HPALD    ' &
         ,'IEPOX    ' &
         ,'XOOH     ' &
         ,'ISOPOOH  ' &
         ,'TOLUENE  ' &
         ,'CRESOL   ' &
         ,'TOLOOH   ' &
         ,'BENZENE  ' &
         ,'PHENOL   ' &
         ,'BEPOMUC  ' &
         ,'PHENOOH  ' &
         ,'C6H5OOH  ' &
         ,'BENZOOH  ' &
         ,'BIGALD1  ' &
         ,'BIGALD2  ' &
         ,'BIGALD3  ' &
         ,'BIGALD4  ' &
         ,'TEPOMUC  ' &
         ,'BZOOH    ' &
         ,'BZALD    ' &
         ,'PBZNIT   ' &
         ,'XYLENES  ' &
         ,'XYLOL    ' &
         ,'XYLOLOOH ' &
         ,'XYLENOOH ' &
         ,'BCARY    ' &
         ,'TERPOOH  ' &
         ,'TERPROD1 ' &
         ,'TERPROD2 ' &
         ,'TERP2OOH ' &
         ,'DMS      ' &
         ,'H2SO4    ' &
         ,'HONITR   ' &
         ,'MACRN    ' &
         ,'MVKN     ' &
         ,'ISOPN2B  ' &
         ,'ISOPN3B  ' &
         ,'ISOPN4D  ' &
         ,'ISOPN1D  ' &
         ,'ISOPNOOHD' &
         ,'ISOPNOOHB' &
         ,'ISOPNBNO3' &
         ,'NO3CH2CHO' &
         ,'HYPERACET' &
         ,'HCOCH2OOH' &
         ,'DHPMPAL  ' &
         ,'MVKOOH   ' &
         ,'ISOPOH   ' &
         ,'ISOPFDN  ' &
         ,'ISOPFNP  ' &
         ,'INHEB    ' &
         ,'HMHP     ' &
         ,'HPALD1   ' &
         ,'INHED    ' &
         ,'HPALD4   ' &
         ,'ISOPHFP  ' &
         ,'HPALDB1C ' &
         ,'HPALDB4C ' &
         ,'ICHE     ' &
         ,'ISOPFDNC ' &
         ,'ISOPFNC  ' &
         ,'TERPNT   ' &
         ,'TERPNS   ' &
         ,'TERPNT1  ' &
         ,'TERPNS1  ' &
         ,'TERPNPT  ' &
         ,'TERPNPS  ' &
         ,'TERPNPT1 ' &
         ,'TERPNPS1 ' &
         ,'TERPFDN  ' &
         ,'SQTN     ' &
         ,'TERPHFN  ' &
         ,'TERP1OOH ' &
         ,'TERPDHDP ' &
         ,'TERPF2   ' &
         ,'TERPF1   ' &
         ,'TERPA    ' &
         ,'TERPA2   ' &
         ,'TERPK    ' &
         ,'TERPAPAN ' &
         ,'TERPACID ' &
         ,'TERPA2PAN' &
         ,'APIN     ' &
         ,'BPIN     ' &
         ,'LIMON    ' &
         ,'MYRC     ' &
         ,'TERPACID2' &
         ,'TERPACID3' &
         ,'TERPA3PAN' &
         ,'TERPOOHL ' &
         ,'TERPA3   ' &
         ,'TERP2AOOH' &
         /)

  !--- data for effective Henry's Law coefficient ---
  real(r8), public, parameter :: dheff(n_species_table*6) = &
            (/1.03e-02_r8, 2830._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! OX
             ,8.70e+04_r8, 7320._r8,2.2e-12_r8,-3730._r8,0._r8     ,    0._r8  & ! H2O2
             ,3.90e+01_r8, 4300._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! OH
             ,6.90e+02_r8, 5900._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! HO2
             ,9.81e-04_r8, 1650._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! CO
             ,1.41e-03_r8, 1820._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! CH4
             ,2.38e+00_r8, 5280._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! CH3O2
             ,3.00e+02_r8, 5280._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! CH3OOH
             ,3.23e+03_r8, 7100._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! CH2O
             ,8.90e+03_r8, 6100._r8,1.8e-04_r8,  -20._r8,0._r8     ,    0._r8  & ! HCOOH
             ,1.92e-03_r8, 1762._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! NO
             ,1.20e-02_r8, 2440._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! NO2
             ,2.10e+05_r8, 8700._r8,2.2e+01_r8,    0._r8,0._r8     ,    0._r8  & ! HNO3
             ,3.44e-02_r8, 2715._r8,4.3e-07_r8,-1000._r8,4.7e-11_r8,-1760._r8  & ! CO2
             ,6.02e+01_r8, 4160._r8,1.7e-05_r8,-4325._r8,1.0e-14_r8,-6716._r8  & ! NH3
             ,2.14e+00_r8, 3362._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! N2O5
             ,3.80e-02_r8,    0._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! NO3
             ,2.03e+02_r8, 5645._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! CH3OH
             ,4.00e+01_r8, 8400._r8,1.3e-06_r8,    0._r8,0._r8     ,    0._r8  & ! HO2NO2
             ,1.00e-16_r8,    0._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! O1D
             ,1.88e-03_r8, 2750._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! C2H6
             ,2.38e+00_r8, 5280._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! C2H5O2
             ,2.38e+00_r8, 5280._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! PO2
             ,2.38e+00_r8, 5280._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! MACRO2
             ,2.38e+00_r8, 5280._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! ISOPO2
             ,1.70e-03_r8,    0._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! C4H10
             ,1.29e+01_r8, 5890._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! CH3CHO
             ,3.36e+02_r8, 5995._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! C2H5OOH
             ,5.57e-03_r8, 2800._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! C3H6
             ,1.50e+06_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! POOH
             ,5.96e-03_r8, 2200._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! C2H4
             ,2.80e+00_r8, 5730._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! PAN
             ,8.37e+02_r8, 5310._r8,1.8e-04_r8,  -20._r8,0._r8     ,    0._r8  & ! CH3COOOH
             ,2.94e-02_r8, 1800._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! MTERP
             ,4.19e+05_r8, 7480._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! GLYOXAL
             ,3.50e+03_r8, 7545._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! CH3COCHO
             ,4.00e+04_r8, 4630._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! GLYALD
             ,1.00e-01_r8,    0._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! CH3CO3
             ,1.51e-03_r8, 3120._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! C3H8
             ,2.38e+00_r8, 5280._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! C3H7O2
             ,2.78e+01_r8, 5530._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! CH3COCH3
             ,3.36e+02_r8, 5995._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! C3H7OOH
             ,2.38e+00_r8, 5280._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! RO2
             ,3.36e+02_r8, 5995._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! ROOH
             ,0.00e+00_r8,    0._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! Rn
             ,3.45e-02_r8, 4400._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! ISOP
             ,4.10e+01_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! MVK
             ,6.50e+00_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! MACR
             ,1.90e+02_r8, 6500._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! C2H5OH
             ,1.44e+03_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! ONITR
             ,1.00e+03_r8, 6000._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! ONIT
             ,2.38e+00_r8, 5280._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! ISOPNO3
             ,1.10e+05_r8, 6000._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! HYDRALD
             ,9.02e+00_r8, 8258._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! HCN
             ,5.28e+01_r8, 3970._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! CH3CN
             ,1.36e+00_r8, 3100._r8,1.30e-02_r8,1960._r8,6.6e-08_r8, 1500._r8  & ! SO2
             ,1.3e+07_r8,     0._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! SOAGff0
             ,3.2e+05_r8,     0._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! SOAGff1
             ,4.0e+05_r8,     0._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! SOAGff2 
             ,1.3e+05_r8,     0._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! SOAGff3
             ,1.6e+05_r8,     0._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! SOAGff4
             ,7.9e+11_r8,     0._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! SOAGbg0
             ,6.3e+10_r8,     0._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! SOAGbg1
             ,3.2e+09_r8,     0._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! SOAGbg2
             ,6.3e+08_r8,     0._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! SOAGbg3
             ,3.2e+07_r8,     0._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! SOAGbg4
             ,4.0e+11_r8,     0._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! SOAG0
             ,3.2e+10_r8,     0._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! SOAG1
             ,1.6e+09_r8,     0._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! SOAG2
             ,3.2e+08_r8,     0._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! SOAG3
             ,1.6e+07_r8,     0._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! SOAG4
             ,1.e+03_r8,      0._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! IVOC
             ,1.e+03_r8,      0._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! SVOC
             ,1.e+03_r8,      0._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! IVOCbb
             ,1.e+03_r8,      0._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! IVOCff
             ,1.e+03_r8,      0._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! SVOCbb
             ,1.e+03_r8,      0._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! SVOCff
             ,2.42e-02_r8, 2710._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! N2O
             ,7.9e-04_r8,   530._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! H2
             ,4.14e-02_r8, 1890._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! C2H2
             ,4.1e+03_r8,  6200._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! CH3COOH
             ,1.9e+06_r8,  6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! EOOH
             ,1.46e+03_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! HYAC
             ,5.96e-03_r8, 2365._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! BIGENE
             ,1.24e-03_r8, 3010._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! BIGALK
             ,1.80e+01_r8, 5740._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! MEK
             ,6.4e+04_r8,  6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! MEKOOH
             ,4.4e+06_r8,  6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! MACROOH
             ,1.72e+00_r8, 5700._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! MPAN
             ,1.01e+00_r8, 5790._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! ALKNIT
             ,1.e+03_r8,   6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! NOA
             ,8.34e+03_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! ISOPNITA
             ,4.82e+04_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! ISOPNITB
             ,8.75e+04_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! ISOPNOOH
             ,1.46e+03_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! NC4CHO
             ,4.02e+04_r8, 9500._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! NC4CH2OH
             ,8.41e+03_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TERPNIT
             ,6.67e+04_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! NTERPOOH
             ,3.36e+02_r8, 5995._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! ALKOOH
             ,9.6e+00_r8,  6220._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! BIGALD
             ,2.30e+05_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! HPALD
             ,3.e+07_r8,   6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! IEPOX
             ,1.e+11_r8,   5995._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! XOOH
             ,3.5e+06_r8,  5995._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! ISOPOOH
             ,1.5e-01_r8,  4300._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TOLUENE
             ,5.67e+02_r8, 5800._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! CRESOL
             ,2.30e+04_r8, 5995._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TOLOOH
             ,1.8e-01_r8,  3800._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! BENZENE
             ,2.84e+03_r8, 2700._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! PHENOL
             ,3.e+07_r8,   6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! BEPOMUC
             ,1.5e+06_r8,  5995._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! PHENOOH
             ,3.36e+02_r8, 5995._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! C6H5OOH
             ,2.3e+03_r8,  5995._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! BENZOOH 
             ,1.e+05_r8,   5890._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! BIGALD1 
             ,2.9e+04_r8,  5890._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! BIGALD2
             ,2.2e+04_r8,  5890._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! BIGALD3
             ,2.2e+04_r8,  5890._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! BIGALD4
             ,2.5e+05_r8,  6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TEPOMUC
             ,3.36e+02_r8, 5995._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! BZOOH
             ,3.24e+01_r8, 6300._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! BZALD
             ,2.8e+00_r8,  5730._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! PBZNIT
             ,2.e-01_r8,   4300._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! XYLENES
             ,1.01e+03_r8, 6800._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! XYLOL
             ,1.9e+06_r8,  5995._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! XYLOLOOH
             ,3.36e+02_r8, 5995._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! XYLENOOH
             ,5.57e-03_r8, 2800._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! BCARY
             ,3.6e+06_r8,  6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TERPOOH
             ,3.92e+04_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TERPROD1
             ,7.20e+04_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TERPROD2
             ,3.36e+02_r8, 5995._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TERP2OOH
             ,5.4e-01_r8,  3460._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! DMS
             ,1.e+11_r8,   6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! H2SO4
             ,2.64e+03_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! HONITR
             ,4.14e+06_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! MACRN
             ,1.84e+05_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! MVKN
             ,8.34e+03_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! ISOPN2B
             ,8.34e+03_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! ISOPN3B
             ,4.82e+04_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! ISOPN4D
             ,4.82e+04_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! ISOPN1D
             ,9.67e+04_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! ISOPNOOHD
             ,6.61e+04_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! ISOPNOOHB
             ,8.34e+03_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! ISOPNBNO3
             ,3.39e+04_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! NO3CH2CHO
             ,1.16e+04_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! HYPERACET
             ,2.99e+04_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! HCOCH2OOH
             ,9.37e+07_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! DHPMPAL
             ,1.24e+06_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! MVKOOH
             ,8.77e+06_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! ISOPOH
             ,5.02e+08_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! ISOPFDN
             ,2.97e+11_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! ISOPFNP
             ,1.05e+05_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! INHEB
             ,1.70e+06_r8, 9870._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! HMHP
             ,2.30e+05_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! HPALD1  
             ,1.51e+05_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! INHED
             ,2.30e+05_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! HPALD4
             ,7.60e+09_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! ISOPHFP 
             ,5.43e+04_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! HPALDB1C
             ,5.43e+04_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! HPALDB4C
             ,2.09e+06_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! ICHE
             ,7.16e+09_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! ISOPFDNC
             ,1.41e+11_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! ISOPFNC
             ,8.41e+03_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TERPNT
             ,8.41e+03_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TERPNS
             ,8.55e+03_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TERPNT1
             ,8.55e+03_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TERPNS1
             ,6.67e+04_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TERPNPT
             ,6.67e+04_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TERPNPS
             ,6.78e+04_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TERPNPT1
             ,6.78e+04_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TERPNPS1
             ,1.65e+09_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TERPFDN
             ,9.04e+03_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! SQTN
             ,7.53e+11_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TERPHFN
             ,3.64e+06_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TERP1OOH
             ,3.41e+14_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TERPDHDP
             ,6.54e+01_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TERPF2
             ,4.05e+04_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TERPF1
             ,3.92e+04_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TERPA
             ,7.20e+04_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TERPA2
             ,6.39e+01_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TERPK
             ,7.94e+03_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TERPAPAN
             ,5.63e+06_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TERPACID
             ,9.59e+03_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TERPA2PAN
             ,2.94e-02_r8, 1800._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! APIN
             ,1.52e-02_r8, 4500._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! BPIN
             ,4.86e-02_r8, 4600._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! LIMON
             ,7.30e-02_r8, 2800._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! MYRC
             ,2.64e+06_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TERPACID2
             ,3.38e+09_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TERPACID3
             ,1.23e+07_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TERPA3PAN
             ,4.41e+12_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TERPOOHL
             ,1.04e+08_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TERPA3
             ,3.67e+06_r8, 6014._r8,0._r8     ,    0._r8,0._r8     ,    0._r8  & ! TERP2AOOH
            /)

  real(r8), private, parameter :: wh2o = SHR_CONST_MWWV
  real(r8), private, parameter :: mol_wgts(n_species_table) = &
       (/ 47.9981995_r8, 34.0135994_r8, 17.0067997_r8, 33.0061989_r8, 28.0104008_r8, &
          16.0405998_r8, 47.0320015_r8, 48.0393982_r8, 30.0251999_r8, 46.0246010_r8, &
          30.0061398_r8, 46.0055389_r8, 63.0123405_r8, 44.0098000_r8, 17.0289402_r8, &
          108.010483_r8, 62.0049400_r8, 32.0400009_r8, 79.0117416_r8, 15.9994001_r8, &
          30.0664005_r8, 61.0578003_r8, 91.0830002_r8, 119.093399_r8, 117.119797_r8, &
          58.1180000_r8, 44.0509987_r8, 62.0652008_r8, 42.0774002_r8, 92.0904007_r8, &
          28.0515995_r8, 121.047943_r8, 76.0497971_r8, 136.228394_r8, 58.0355988_r8, &
          72.0614014_r8, 60.0503998_r8, 75.0423965_r8, 44.0922012_r8, 75.0836029_r8, &
          58.0768013_r8, 76.0910034_r8, 89.070126_r8,  90.078067_r8,  222.000000_r8, &
          68.1141968_r8, 70.0877991_r8, 70.0877991_r8, 46.0657997_r8, 147.125946_r8, &
          119.074341_r8, 162.117935_r8, 100.112999_r8, 27.0256_r8   , 41.0524_r8   , &
          64.064800_r8,  250._r8,       250._r8,       250._r8,       250._r8,       &
          250._r8,       250._r8,       250._r8,       250._r8,       250._r8,       &
          250._r8,       250._r8,       250._r8,       250._r8,       250._r8,       &
          250._r8,       170.3_r8,      170.3_r8,      170.3_r8,       170.3_r8,     &
          170.3_r8,      170.3_r8,      44.0129_r8,    2.0148_r8,     26.0368_r8,    &
          60.0504_r8,    78.0646_r8,    74.0762_r8,    56.1032_r8,    72.1438_r8,    &
          72.1026_r8,    104.101_r8,    120.101_r8,    147.085_r8,    133.141_r8,    &
          119.074_r8,    147.126_r8,    147.126_r8,    163.125_r8,    145.111_r8,    &
          147.126_r8,    215.24_r8,     231.24_r8,     104.143_r8,    98.0982_r8,    &
          116.112_r8,    118.127_r8,    150.126_r8,    118.127_r8,    92.1362_r8,    &
          108.136_r8,    174.148_r8,    78.1104_r8,    94.1098_r8,    126.109_r8,    &
          176.122_r8,    110.109_r8,    160.122_r8,    84.0724_r8,    98.0982_r8,    &
          98.0982_r8,    112.124_r8,    140.134_r8,    124.135_r8,    106.121_r8,    &
          183.118_r8,    106.162_r8,    122.161_r8,    204.173_r8,    188.174_r8,    &
          204.343_r8,    186.241_r8,    168.227_r8,    154.201_r8,    200.226_r8,    &
          62.1324_r8,    98.0784_r8,    135.118733_r8, 149.102257_r8, 149.102257_r8, &
          147.129469_r8, 147.129469_r8, 147.129469_r8, 147.129469_r8, 163.128874_r8, &
          163.128874_r8, 147.129469_r8, 105.049617_r8, 90.078067_r8,  76.05145_r8,   &
          136.103494_r8, 120.104089_r8, 102.131897_r8, 226.141733_r8, 197.143565_r8, &
          163.128874_r8, 64.040714_r8,  116.11542_r8,  163.128874_r8, 116.11542_r8,  &
          150.130112_r8, 116.11542_r8,  116.11542_r8,  116.11542_r8,  224.125851_r8, &
          195.127684_r8, 215.246675_r8, 215.246675_r8, 215.246675_r8, 215.246675_r8, &
          231.24608_r8,  231.24608_r8,  231.24608_r8,  231.24608_r8,  294.258938_r8, &
          283.36388_r8,  265.260771_r8, 186.248507_r8, 236.262604_r8, 110.153964_r8, &
          168.233221_r8, 168.233221_r8, 154.206603_r8, 138.207199_r8, 245.229603_r8, &
          200.232031_r8, 231.202986_r8, 136.228394_r8, 136.228394_r8, 136.228394_r8, &
          136.228394_r8, 186.205413_r8, 202.204818_r8, 247.202391_r8, 218.247317_r8, &
          170.206008_r8, 186.248507_r8 /)


!===============================================================================
CONTAINS
!===============================================================================

  subroutine seq_drydep_readnl(NLFilename, drydep_nflds)

    !========================================================================
    ! reads drydep_inparm namelist and determines the number of drydep velocity
    ! fields that are sent from the land component
    !========================================================================

    character(len=*), intent(in)  :: NLFilename ! Namelist filename
    integer, intent(out)          :: drydep_nflds

    !----- local -----
    integer       :: i                ! Indices
    integer       :: unitn            ! namelist unit number
    integer       :: ierr             ! error code
    logical       :: exists           ! if file exists or not
    type(ESMF_VM) :: vm
    integer       :: localPet
    integer       :: mpicom
    integer       :: rc
    character(*),parameter :: F00   = "('(seq_drydep_read) ',8a)"
    character(*),parameter :: FI1   = "('(seq_drydep_init) ',a,I2)"
    character(*),parameter :: subName = '(seq_drydep_read) '
    !-----------------------------------------------------------------------------

    namelist /drydep_inparm/ drydep_list, drydep_method

    !-----------------------------------------------------------------------------
    ! Read namelist and figure out the drydep field list to pass
    ! First check if file exists and if not, n_drydep will be zero
    !-----------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !--- Open and read namelist ---
    if ( len_trim(NLFilename) == 0  )then
       call shr_sys_abort( subName//'ERROR: nlfilename not set' )
    end if

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return 

    call ESMF_VMGet(vm, localPet=localPet, mpiCommunicator=mpicom, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return 

    if (localPet==0) then
       inquire( file=trim(NLFileName), exist=exists)
       if ( exists ) then
          open(newunit=unitn, file=trim(NLFilename), status='old' )
          write(s_logunit,F00) 'Read in drydep_inparm namelist from: ', trim(NLFilename)
          call shr_nl_find_group_name(unitn, 'drydep_inparm', ierr)
          if (ierr == 0) then
             ! Note that ierr /= 0, no namelist is present.
             read(unitn, drydep_inparm, iostat=ierr)
             if (ierr > 0) then
                call shr_sys_abort( 'problem on read of drydep_inparm namelist in seq_drydep_readnl')
             end if
          endif
          close( unitn )
       end if
    end if
    call shr_mpi_bcast( drydep_list, mpicom )
    call shr_mpi_bcast( drydep_method, mpicom )

    do i=1,maxspc
       if(len_trim(drydep_list(i)) > 0) then
          drydep_nflds=drydep_nflds+1
       endif
    enddo

    ! set module variable
    n_drydep = drydep_nflds

    ! Make sure method is valid and determine if land is passing drydep fields
    lnd_drydep = (drydep_nflds>0 .and. drydep_method == DD_XLND)
    if (localpet==0) then
       write(s_logunit,*) 'seq_drydep_read: drydep_method: ', trim(drydep_method)
       if ( drydep_nflds == 0 )then
          write(s_logunit,F00) 'No dry deposition fields will be transfered'
       else
          write(s_logunit,FI1) 'Number of dry deposition fields transfered is ', drydep_nflds
       end if
    end if

    if ( trim(drydep_method)/=trim(DD_XATM) .and. &
         trim(drydep_method)/=trim(DD_XLND) .and. &
         trim(drydep_method)/=trim(DD_TABL) ) then
       write(s_logunit,*) 'seq_drydep_read: drydep_method : ', trim(drydep_method)
       write(s_logunit,*) 'seq_drydep_read: drydep_method must be set to : ', &
            DD_XATM,', ', DD_XLND,', or ', DD_TABL
       call shr_sys_abort('seq_drydep_read: incorrect dry deposition method specification')
    endif

    if (.not. drydep_initialized) then
       call seq_drydep_init()
    end if

  end subroutine seq_drydep_readnl

!====================================================================================

  subroutine seq_drydep_init( )

    !========================================================================
    ! Initialization of dry deposition fields
    ! reads drydep_inparm namelist and sets up CCSM driver list of fields for
    ! land-atmosphere communications.
    !========================================================================

    !----- local -----
    integer :: i, l                      ! Indices
    character(len=32) :: test_name       ! field test name

    !----- formats -----
    character(*),parameter :: subName = '(seq_drydep_init) '
    character(*),parameter :: F00   = "('(seq_drydep_init) ',8a)"

    !-----------------------------------------------------------------------------
    ! Return if this routine has already been called (e.g. cam and clm both call this)
    !-----------------------------------------------------------------------------

    !-----------------------------------------------------------------------------
    ! Allocate and fill foxd, drat and mapping as well as species indices
    !-----------------------------------------------------------------------------

    if ( n_drydep > 0 ) then

       allocate( foxd(n_drydep) )
       allocate( drat(n_drydep) )
       allocate( mapping(n_drydep) )

       ! This initializes these variables to infinity.
       foxd = shr_infnan_posinf
       drat = shr_infnan_posinf

       mapping(:) = 0

    end if

    h2_ndx=-1; ch4_ndx=-1; co_ndx=-1; mpan_ndx = -1; pan_ndx = -1; so2_ndx=-1; o3_ndx=-1; xpan_ndx=-1

    !--- Loop over drydep species that need to be worked with ---
    do i=1,n_drydep
       if ( len_trim(drydep_list(i))==0 ) exit

       test_name = drydep_list(i)

       if( trim(test_name) == 'O3' ) then
          test_name = 'OX'
       end if

       !--- Figure out if species maps to a species in the species table ---
       do l = 1,n_species_table
          if(  trim( test_name ) == trim( species_name_table(l) ) ) then
             mapping(i)  = l
             exit
          end if
       end do

       !--- If it doesn't map to a species in the species table find species close enough ---
       if( mapping(i) < 1 ) then
          select case( trim(test_name) )
          case( 'O3S', 'O3INERT' )
             test_name = 'OX'
          case( 'Pb' )
             test_name = 'HNO3'
          case( 'SOGM','SOGI','SOGT','SOGB','SOGX' )
             test_name = 'CH3OOH'
          case( 'SOA', 'SO4', 'CB1', 'CB2', 'OC1', 'OC2', 'NH4', 'SA1', 'SA2', 'SA3', 'SA4' )
             test_name = 'OX'  ! this is just a place holder. values are explicitly set below
          case( 'SOAM', 'SOAI', 'SOAT', 'SOAB', 'SOAX' )
             test_name = 'OX'  ! this is just a place holder. values are explicitly set below
          case( 'SOAGbb0' )
             test_name = 'SOAGff0'
          case( 'SOAGbb1' )
             test_name = 'SOAGff1'
          case( 'SOAGbb2' )
             test_name = 'SOAGff2'
          case( 'SOAGbb3' )
             test_name = 'SOAGff3'
          case( 'SOAGbb4' )
             test_name = 'SOAGff4'
          case( 'O3A' )
             test_name = 'OX'
          case( 'XMPAN' )
             test_name = 'MPAN'
          case( 'XPAN' )
             test_name = 'PAN'
          case( 'XNO' )
             test_name = 'NO'
          case( 'XNO2' )
             test_name = 'NO2'
          case( 'XHNO3' )
             test_name = 'HNO3'
          case( 'XONIT' )
             test_name = 'ONIT'
          case( 'XONITR' )
             test_name = 'ONITR'
          case( 'XHO2NO2')
             test_name = 'HO2NO2'
          case( 'XNH4NO3' )
             test_name = 'HNO3'
          case( 'NH4NO3' )
             test_name = 'HNO3'
          case default
             test_name = 'blank'
          end select

          !--- If found a match check the species table again ---
          if( trim(test_name) /= 'blank' ) then
             do l = 1,n_species_table
                if( trim( test_name ) == trim( species_name_table(l) ) ) then
                   mapping(i)  = l
                   exit
                end if
             end do
          else
             write(s_logunit,F00) trim(drydep_list(i)),' not in tables; will have dep vel = 0'
             call shr_sys_abort( subName//': '//trim(drydep_list(i))//' is not in tables' )
          end if
       end if

       !--- Figure out the specific species indices ---
       if ( trim(drydep_list(i)) == 'H2' )   h2_ndx   = i
       if ( trim(drydep_list(i)) == 'CO' )   co_ndx   = i
       if ( trim(drydep_list(i)) == 'CH4' )  ch4_ndx  = i
       if ( trim(drydep_list(i)) == 'MPAN' ) mpan_ndx = i
       if ( trim(drydep_list(i)) == 'PAN' )  pan_ndx  = i
       if ( trim(drydep_list(i)) == 'SO2' )  so2_ndx  = i
       if ( trim(drydep_list(i)) == 'OX' .or. trim(drydep_list(i)) == 'O3' ) o3_ndx  = i
       if ( trim(drydep_list(i)) == 'O3A' ) o3a_ndx  = i
       if ( trim(drydep_list(i)) == 'XPAN' ) xpan_ndx = i

       if( mapping(i) > 0) then
         l = mapping(i)
         foxd(i) = dfoxd(l)
         drat(i) = sqrt(mol_wgts(l)/wh2o)
       endif

    enddo

    where( rgss < 1._r8 )
       rgss = 1._r8
    endwhere

    where( rac < small_value)
       rac = small_value
    endwhere

    drydep_initialized = .true.

  end subroutine seq_drydep_init

!====================================================================================

  subroutine set_hcoeff_scalar( sfc_temp, heff )

    !========================================================================
    ! Interface to seq_drydep_setHCoeff when input is scalar
    ! wrapper routine used when surface temperature is a scalar (single column) rather
    ! than an array (multiple columns).
    !
    ! !REVISION HISTORY:
    !  2008-Nov-12 - F. Vitt - first version
    !========================================================================

    implicit none

    real(r8), intent(in)     :: sfc_temp         ! Input surface temperature
    real(r8), intent(out)    :: heff(n_drydep)   ! Output Henry's law coefficients

    !----- local -----
    real(r8) :: sfc_temp_tmp(1)    ! surface temp

    sfc_temp_tmp(:) = sfc_temp
    call set_hcoeff_vector( 1, sfc_temp_tmp, heff(:n_drydep) )

  end subroutine set_hcoeff_scalar

!====================================================================================

  subroutine set_hcoeff_vector( ncol, sfc_temp, heff )

    !========================================================================
    ! Interface to seq_drydep_setHCoeff when input is vector
    ! sets dry depositions coefficients -- used by both land and atmosphere models
    !========================================================================

    integer, intent(in)      :: ncol                  ! Input size of surface-temp vector
    real(r8), intent(in)     :: sfc_temp(ncol)        ! Surface temperature
    real(r8), intent(out)    :: heff(ncol,n_drydep)   ! Henry's law coefficients

    !----- local -----
    real(r8), parameter :: t0     = 298._r8    ! Standard Temperature
    real(r8), parameter :: ph_inv = 1._r8/ph   ! Inverse of PH
    integer  :: m, l, id       ! indices
    real(r8) :: e298           ! Henry's law coefficient @ standard temperature (298K)
    real(r8) :: dhr            ! temperature dependence of Henry's law coefficient
    real(r8) :: dk1s(ncol)     ! DK Work array 1
    real(r8) :: dk2s(ncol)     ! DK Work array 2
    real(r8) :: wrk(ncol)      ! Work array

    !----- formats -----
    character(*),parameter :: subName = '(seq_drydep_set_hcoeff) '
    character(*),parameter :: F00   = "('(seq_drydep_set_hcoeff) ',8a)"

    !-------------------------------------------------------------------------------
    ! notes:
    !-------------------------------------------------------------------------------

    wrk(:) = (t0 - sfc_temp(:))/(t0*sfc_temp(:))
    do m = 1,n_drydep
       l    = mapping(m)
       id   = 6*(l - 1)
       e298 = dheff(id+1)
       dhr  = dheff(id+2)
       heff(:,m) = e298*exp( dhr*wrk(:) )
       !--- Calculate coefficients based on the drydep tables ---
       if( dheff(id+3) /= 0._r8 .and. dheff(id+5) == 0._r8 ) then
          e298 = dheff(id+3)
          dhr  = dheff(id+4)
          dk1s(:) = e298*exp( dhr*wrk(:) )
          where( heff(:,m) /= 0._r8 )
             heff(:,m) = heff(:,m)*(1._r8 + dk1s(:)*ph_inv)
          elsewhere
             heff(:,m) = dk1s(:)*ph_inv
          endwhere
       end if
       !--- For coefficients that are non-zero AND CO2 or NH3 handle things this way ---
       if( dheff(id+5) /= 0._r8 ) then
          if( trim( drydep_list(m) ) == 'CO2' .or. trim( drydep_list(m) ) == 'NH3' &
              .or. trim( drydep_list(m) ) == 'SO2' ) then
             e298 = dheff(id+3)
             dhr  = dheff(id+4)
             dk1s(:) = e298*exp( dhr*wrk(:) )
             e298 = dheff(id+5)
             dhr  = dheff(id+6)
             dk2s(:) = e298*exp( dhr*wrk(:) )
             !--- For Carbon dioxide ---
             if( trim(drydep_list(m)) == 'CO2'.or. trim( drydep_list(m) ) == 'SO2' ) then
                heff(:,m) = heff(:,m)*(1._r8 + dk1s(:)*ph_inv*(1._r8 + dk2s(:)*ph_inv))
             !--- For NH3 ---
             else if( trim( drydep_list(m) ) == 'NH3' ) then
                heff(:,m) = heff(:,m)*(1._r8 + dk1s(:)*ph/dk2s(:))
             !--- This can't happen ---
             else
                write(s_logunit,F00) 'Bad species ',drydep_list(m)
                call shr_sys_abort( subName//'ERROR: in assigning coefficients' )
             end if
          end if
       end if
    end do

  end subroutine set_hcoeff_vector

!===============================================================================

end module seq_drydep_mod

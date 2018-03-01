module Tracer_varcon
  !-----------------------------------------------------------------------
  !BOP
  !
  ! !MODULE: Tracer_varcon
  !
  ! !DESCRIPTION: Module containing parameters and logical switches
  ! and routine to read constants from namelist for tracer
  ! transport set up.
  !
  ! NOTE(bja, 201604) Do NOT add a save statement to this module or
  ! assign values to non-parameter variables. If any value can be
  ! changed, it should be read into a class and an instance of the
  ! class passed to routines that need th evalue.
  !
  ! !USES:
  use bshr_kind_mod, only : r8 => shr_kind_r8
  use betr_constants, only : betr_string_length
  implicit none
  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  !
  ! !PUBLIC TYPES:

  public

  real(r8), parameter :: SHR_CONST_VSMOW_O18 = 2005.20e-6_R8 ! ratio of 18O/16O in Vienna Standard Mean Ocean Water (VSMOW)
  real(r8), parameter :: SHR_CONST_VSMOW_O17 = 379.9e-6_R8   ! ratio of 17O/16O in Vienna Standard Mean Ocean Water (VSMOW)
  real(r8), parameter :: SHR_CONST_VSMOW_D = 155.76e-6_R8    ! ratio of D/H in Vienna Standard Mean Ocean Water (VSMOW)
  real(r8), parameter :: SHR_CONST_VSMOW_T = 1.85e-6_R8      ! ratio of T/H in Vienna Standard Mean Ocean Water (VSMOW)

  real(r8), parameter :: catomw = 12.011_r8     ! molar mass of C atoms (g/mol)
  real(r8), parameter :: natomw = 14.007_r8     ! molar mass of N atoms (g/mol)
  real(r8), parameter :: patomw = 30.97_r8      ! molar mass of P atmos (g/mol)
  real(r8), parameter :: c13atomw = 13._r8      ! molar mass of C13 atmos (g/mol)
  real(r8), parameter :: c14atomw = 14._r8      ! molar mass of C14 atmos (g/mol)
  integer, parameter  :: bndcond_as_conc = 1 ! top boundary conditions as tracer concentration
  integer, parameter  :: bndcond_as_flux = 2 ! top boundary condition as tracer flux

  logical, parameter  :: l2ndadvsolver = .false. ! by default use 1st order solver for advection

  logical :: is_active_betr_bgc=.false.
  logical :: use_c13_betr=.false.
  logical :: use_c14_betr=.false.
  logical :: is_nitrogen_active = .true.
  logical :: is_phosphorus_active=.true.
  integer, parameter :: sorp_isotherm_linear=1
  integer, parameter :: sorp_isotherm_langmuir=2

  logical  :: advection_on, diffusion_on, reaction_on, ebullition_on
  character(len=betr_string_length), public :: reaction_method
  save

  integer, public :: betr_nlevsoi
  integer, public :: betr_nlevsno
  integer, public :: betr_nlevtrc_soil

  !
  ! NOTE(bja, 201604) Do NOT add a save statement to this module and
  ! create instances of these types. Instances should be created in a
  ! relevant class.
  !

  ! the following will likely become case dependent
!X!  type, public :: betr_respiration_parameters_type
!X!     real(r8), public :: rr_dif_scal = 1._r8     ! scaling factor for how much root respiration is diffused out into soil
!X!     real(r8), public :: mr_dif_scal = 0._r8     ! how much fraction of stem respiration is back into xylem
!X!     real(r8), public :: co2_refix_scal = 0.0_r8 ! how much fraction of co2 in the xylem is refixed in leaf
!X!  end type betr_respiration_parameters_type


!X!  type, public :: betr_soil_chem_constants_type
!X!     real(r8), public :: site_pH = 7._r8 ! pH value of the site
!X!  end type betr_soil_chem_constants_type
!X!
!X!  type, public :: betr_atm_composition_type
!X!     ! atmospheric compositions, (v/v)%
!X!     real(r8), public :: atm_n2 = 0.78084_r8
!X!     real(r8), public :: atm_o2 = 0.20946_r8
!X!     real(r8), public :: atm_ar = 0.009340_r8
!X!     real(r8), public :: atm_co2 = 379e-6_r8   ! this will be set to the value provided from co2_ppmv
!X!     real(r8), public :: atm_ch4 = 1.7e-6_r8   ! this will be set to the value provided from atmch4 if clm4me is on
!X!     real(r8), public :: atm_n2o = 3.1e-7_r8
!X!     real(r8), public :: atm_no = 4.56e-9_r8
!X!     real(r8), public :: atm_nh3 = 300.e-12_r8
!X!     real(r8), public :: atm_h2 = 0.55e-6_r8
!X!  end type betr_atm_composition_type
!X!
!X!  type, public :: betr_atm_isotope_composition_type
!X!     ! atmospheric isotopic signatures
!X!     ! the zeros will be replaced with updated value from literature searching.
!X!     real(r8), public :: atm_deld_h2 = 0._r8    ! relative to VSMOW
!X!     real(r8), public :: atm_delt_h2 = 0._r8    ! relative to VSMOW
!X!     real(r8), public :: atm_del13c_co2 = -6._r8 ! set to pre-industrial value by default, it will be used to set the value of c13ratio, PDB
!X!     real(r8), public :: atm_del13c_ch4 = 0._r8 ! relative to PDB
!X!     real(r8), public :: atm_del14c_co2 = 0._r8 ! relative to what?
!X!     real(r8), public :: atm_del14c_ch4 = 0._r8 ! relative to what?
!X!     real(r8), public :: atm_del18o_co2 = 0._r8 ! relative to VSMOW
!X!     real(r8), public :: atm_del18o_h2o = 0._r8 ! relative to VSMOW
!X!     real(r8), public :: atm_del18o_o2 = 0._r8 ! relative to VSMOW
!X!     real(r8), public :: atm_del17o_co2 = 0._r8 ! relative to VSMOW
!X!     real(r8), public :: atm_del17o_h2o = 0._r8 ! relative to VSMOW
!X!     real(r8), public :: atm_del17o_o2 = 0._r8 ! relative to VSMOW
!X!     real(r8), public :: atm_deld_ch4 = 0._r8 ! realtive to VSMOW
!X!     real(r8), public :: atm_deld_h2o = 0._r8 ! relative to VSMOW
!X!
!X!     ! true fractions of the isotopologues in the atmosphere
!X!     real(r8), public :: atm_dratio_h2, atm_tratio_h2
!X!     real(r8), public :: atm_c13rc12_co2, atm_c14rc12_co2, atm_o18ro16_co2, atm_o17ro16_co2
!X!     real(r8), public :: atm_drh_h2o, atm_tratio_h2o, atm_o18ro16_h2o, atm_o17ro16_h2o
!X!     real(r8), public :: atm_c13rc12_ch4, atm_c14rc12_ch4, atm_drh_ch4
!X!  end type betr_atm_isotope_composition_type
   public :: set_cnpbgc
   contains

!-------------------------------------------------------------------------------
    subroutine set_cnpbgc(cnpset)
    !
    !DESCRIPTION
    !set n and p switches of the betr bgc model
    implicit none
    character(len=*), intent(in) :: cnpset
    !set
    select case (trim(cnpset))
    case ('C')
      is_nitrogen_active = .false.
      is_phosphorus_active=.false.
    case ('CN')
      is_nitrogen_active = .true.
      is_phosphorus_active=.false.
    case ('CP')
      is_nitrogen_active = .false.
      is_phosphorus_active=.true.
    case ('CNP')
      is_nitrogen_active = .true.
      is_phosphorus_active=.true.
    case default
      is_nitrogen_active = .true.
      is_phosphorus_active=.true.
    end select
    end subroutine set_cnpbgc
end module Tracer_varcon

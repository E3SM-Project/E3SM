module Tracer_varcon
  !-----------------------------------------------------------------------
  !BOP
  !
  ! !MODULE: Tracer_varcon
  !
  ! !DESCRIPTION:
  ! Module containing parameters and logical switches and routine to read constants from CLM namelist for tracer transport set up.
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils  , only : endrun
  use clm_varctl  , only : iulog
  !
  ! !PUBLIC TYPES:
  implicit none
  save

  logical,  public           :: l2ndadvsolver = .false.             ! by default use 1st order solver for advection

  real(r8),public, parameter :: SHR_CONST_VSMOW_O18 = 2005.20e-6_R8 ! ratio of 18O/16O in Vienna Standard Mean Ocean Water (VSMOW)
  real(r8),public, parameter :: SHR_CONST_VSMOW_O17 = 379.9e-6_R8   ! ratio of 17O/16O in Vienna Standard Mean Ocean Water (VSMOW)
  real(r8),public, parameter :: SHR_CONST_VSMOW_D = 155.76e-6_R8    ! ratio of D/H in Vienna Standard Mean Ocean Water (VSMOW)
  real(r8),public, parameter :: SHR_CONST_VSMOW_T = 1.85e-6_R8      ! ratio of T/H in Vienna Standard Mean Ocean Water (VSMOW)

  ! underground tracer transport logical switches
  logical, public            :: ltracer_offline=.true.       ! true=> do not pass volatile tracers from/to atmosphere
  logical, public            :: ltrcunsat=.false.            ! ture=> swith on tracer transport for specified underground processes in unsaturated upland soil
  logical, public            :: ltrcsat  =.false.            ! ture=> swith on tracer transport for specified underground processes in unsaturated wetland soil
  logical, public            :: ltrclake =.false.            ! ture=> swith on tracer transport for specified underground processes in lake water and lake soil
  logical, public            :: laquadv_off =.false.         ! true=> turn off aqueous advection
  logical, public            :: lgasadv_off = .false.        ! true=> turn off gas advection
  logical, public            :: lzero_restart = .false.      ! true => start with nil tracer concentration, by default
  logical, public            :: is_online_soilchem = .false. ! true=> chemistry is done outside TracerUpdate, added for plug&play capability, say microbial model
  logical, public            :: ldsolvn_vtransport = .false. ! this is not in the namelist, and its value will be determined in SoilTracersMod, don transport?
  logical, public            :: ldsolvc_vtransport = .false. ! this is not in the namelist, and its value will be determined in SoilTracersMod, doc transport?
  logical, public            :: lco2_refix = .false.         ! true => refix co2 transported to leaf
  logical, public            :: lneut = .false.              ! true => only allow neutral molecules to go through xylem
  logical, public            :: ltracer_stem = .false.       ! true => model valatile tracer in stem
  logical, public            :: use_pH_data = .false.
  logical, public            :: licecoat = .false.           ! true => switch on ice coating for dissolved tracers, the coating is defined as the dice/h2oliq,
  ! where, dice is the change of ice content during to free-thaw cyles
  logical, public            :: is_active_betr_bgc = .false.
  logical, public            :: do_betr_leaching = .false.
  logical, public            :: liceseal = .true.       ! true => allow ice to seal the surface soil and keep the gas tracer
  real(r8),public            :: rr_dif_scal = 1._r8     ! scaling factor for how much root respiration is diffused out into soil
  real(r8),public            :: mr_dif_scal = 0._r8     ! how much fraction of stem respiration is back into xylem
  real(r8),public            :: co2_refix_scal = 0.0_r8 ! how much fraction of co2 in the xylem is refixed in leaf
  real(r8),public            :: site_pH = 7._r8         ! pH value of the site

  !  atmospheric compositions, (v/v)
  real(r8),public            :: atm_n2  = 0.78084_r8
  real(r8),public            :: atm_o2  = 0.20946_r8
  real(r8),public            :: atm_ar  = 0.009340_r8
  real(r8),public            :: atm_co2 = 379e-6_r8   !this will be set to the value provided from co2_ppmv
  real(r8),public            :: atm_ch4 = 1.7e-6_r8   !this will be set to the value provided from atmch4 if clm4me is on
  real(r8),public            :: atm_n2o = 3.1e-7_r8
  real(r8),public            :: atm_no  = 4.56e-9_r8
  real(r8),public            :: atm_nh3 = 300.e-12_r8 !
  real(r8),public            :: atm_h2 =  0.55e-6_r8

  !  atmospheric isotopic signatures
  !  the zeros will be replaced with updated value from literature searching.
  real(r8),public            :: atm_deld_h2 = 0._r8    !relative to VSMOW
  real(r8),public            :: atm_delt_h2 = 0._r8    !relative to VSMOW
  real(r8),public            :: atm_del13c_co2 =-6._r8 !set to pre-industrial value by default, it will be used to set the value of c13ratio, PDB
  real(r8),public            :: atm_del13c_ch4 = 0._r8 !relative to PDB
  real(r8),public            :: atm_del14c_co2 = 0._r8 !relative to what?
  real(r8),public            :: atm_del14c_ch4 = 0._r8 !relative to what?
  real(r8),public            :: atm_del18o_co2 = 0._r8 !relative to VSMOW
  real(r8),public            :: atm_del18o_h2o = 0._r8 !relative to VSMOW
  real(r8),public            :: atm_del18o_o2  = 0._r8 !relative to VSMOW
  real(r8),public            :: atm_del17o_co2 = 0._r8 !relative to VSMOW
  real(r8),public            :: atm_del17o_h2o = 0._r8 !relative to VSMOW
  real(r8),public            :: atm_del17o_o2  = 0._r8 !relative to VSMOW
  real(r8),public            :: atm_deld_ch4   = 0._r8 !realtive to VSMOW
  real(r8),public            :: atm_deld_h2o   = 0._r8 !relative to VSMOW

  integer, parameter, public :: bndcond_as_conc = 1    !top boundary conditions as tracer concentration
  integer, parameter, public :: bndcond_as_flux=2      !top boundary condition as tracer flux


  !true fractions of the isotopologues in the atmosphere
  real(r8),public            :: atm_dratio_h2, atm_tratio_h2
  real(r8),public            :: atm_c13rc12_co2, atm_c14rc12_co2, atm_o18ro16_co2, atm_o17ro16_co2
  real(r8),public            :: atm_drh_h2o,atm_tratio_h2o,atm_o18ro16_h2o, atm_o17ro16_h2o
  real(r8),public            :: atm_c13rc12_ch4, atm_c14rc12_ch4, atm_drh_ch4

end module Tracer_varcon

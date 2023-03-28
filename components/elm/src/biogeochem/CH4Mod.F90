module CH4Mod

   #include "shr_assert.h"
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module holding routines to calculate methane fluxes
  ! The driver averages up to gridcell, weighting by finundated, and checks for balance errors.
  ! Sources, sinks, "competition" for CH4 & O2, & transport are resolved in ch4_tran.
  !
  ! !USES:
  use shr_kind_mod       , only : r8 => shr_kind_r8
  use shr_infnan_mod     , only : nan => shr_infnan_nan
  use shr_log_mod        , only : errMsg => shr_log_errMsg
  use elm_varpar         , only : nlevsoi, ngases, nlevsno, nlevdecomp
  use elm_varcon         , only : denh2o, denice, tfrz, grav, spval, rgas, grlnd
  use elm_varcon         , only : catomw, s_con, d_con_w, d_con_g, c_h_inv, kh_theta, kh_tbase
  use landunit_varcon    , only : istdlak
  use clm_time_manager   , only : get_step_size, get_nstep
  use elm_varctl         , only : iulog, use_cn, use_lch4
  use abortutils         , only : endrun
  use decompMod          , only : bounds_type
  use SharedParamsMod    , only : ParamsShareInst
  use atm2lndType        , only : atm2lnd_type
  use CanopyStateType    , only : canopystate_type
  use EnergyFluxType     , only : energyflux_type
  use LakeStateType      , only : lakestate_type
  use lnd2atmType        , only : lnd2atm_type
  use SoilHydrologyType  , only : soilhydrology_type
  use SoilStateType      , only : soilstate_type
  use GridcellType       , only : grc_pp
  use TopounitDataType   , only : top_as  ! for topounit-level atmospheric state forcing
  use LandunitType       , only : lun_pp
  use ColumnType         , only : col_pp
  use ColumnDataType     , only : col_es, col_ws, col_wf, col_cf, col_nf
  use VegetationType     , only : veg_pp
  use VegetationDataType , only : veg_wf, veg_cs, veg_cf
  use timeinfoMod

  !
  implicit none
  save
  private

  ! Non-tunable constants
  real(r8) :: rgasm  ! J/mol.K; rgas / 1000; will be set below
  real(r8), parameter :: rgasLatm = 0.0821_r8 ! L.atm/mol.K
  !$acc declare copyin(rgasLatm)
  !$acc declare copyin(rgasm)

  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: readCH4Params
  public  :: CH4

  ! !PRIVATE MEMBER FUNCTIONS:
  private :: ch4_prod
  private :: ch4_oxid
  private :: ch4_aere
  private :: ch4_ebul
  private :: ch4_tran
  private :: ch4_annualupdate
  private :: get_jwt

  type, private :: CH4ParamsType
     ! ch4 production constants
     real(r8), pointer :: q10ch4            => null()    ! additional Q10 for methane production ABOVE the soil decomposition temperature relationship
     real(r8), pointer :: q10ch4base        => null()    ! temperature at which the effective f_ch4 actually equals the constant f_ch4
     real(r8), pointer :: f_ch4             => null()    ! ratio of CH4 production to total C mineralization
     real(r8), pointer :: rootlitfrac       => null()    ! Fraction of soil organic matter associated with roots
     real(r8), pointer :: cnscalefactor     => null()    ! scale factor on CN decomposition for assigning methane flux
     real(r8), pointer :: redoxlag          => null()    ! Number of days to lag in the calculation of finundated_lag
     real(r8), pointer :: lake_decomp_fact  => null()    ! Base decomposition rate (1/s) at 25C
     real(r8), pointer :: redoxlag_vertical => null()    ! time lag (days) to inhibit production for newly unsaturated layers
     real(r8), pointer :: pHmax             => null()    ! maximum pH for methane production(= 9._r8)
     real(r8), pointer :: pHmin             => null()    ! minimum pH for methane production(= 2.2_r8)
     real(r8), pointer :: oxinhib           => null()    ! inhibition of methane production by oxygen (m^3/mol)

     ! ch4 oxidation constants
     real(r8), pointer :: vmax_ch4_oxid     => null()     ! oxidation rate constant (= 45.e-6_r8 * 1000._r8 / 3600._r8) [mol/m3-w/s];
     real(r8), pointer :: k_m               => null()     ! Michaelis-Menten oxidation rate constant for CH4 concentration
     real(r8), pointer :: q10_ch4oxid       => null()     ! Q10 oxidation constant
     real(r8), pointer :: smp_crit          => null()     ! Critical soil moisture potential
     real(r8), pointer :: k_m_o2            => null()     ! Michaelis-Menten oxidation rate constant for O2 concentration
     real(r8), pointer :: k_m_unsat         => null()     ! Michaelis-Menten oxidation rate constant for CH4 concentration
     real(r8), pointer :: vmax_oxid_unsat   => null()     ! (= 45.e-6_r8 * 1000._r8 / 3600._r8 / 10._r8) [mol/m3-w/s]

     ! ch4 aerenchyma constants
     real(r8), pointer :: aereoxid => null()            ! fraction of methane flux entering aerenchyma rhizosphere that will be

     ! oxidized rather than emitted
     real(r8), pointer :: scale_factor_aere  => null()   ! scale factor on the aerenchyma area for sensitivity tests
     real(r8), pointer :: nongrassporosratio => null()   ! Ratio of root porosity in non-grass to grass, used for aerenchyma transport
     real(r8), pointer :: unsat_aere_ratio   => null()   ! Ratio to multiply upland vegetation aerenchyma porosity by compared to inundated systems (= 0.05_r8 / 0.3_r8)
     real(r8), pointer :: porosmin           => null()   ! minimum aerenchyma porosity (unitless)(= 0.05_r8)

     ! ch4 ebbulition constants
     real(r8), pointer :: vgc_max => null()              ! ratio of saturation pressure triggering ebullition

     ! ch4 transport constants
     real(r8), pointer :: satpow               => null() ! exponent on watsat for saturated soil solute diffusion
     real(r8), pointer :: scale_factor_gasdiff => null() ! For sensitivity tests; convection would allow this to be > 1
     real(r8), pointer :: scale_factor_liqdiff => null() ! For sensitivity tests; convection would allow this to be > 1
     real(r8), pointer :: capthick             => null() ! min thickness before assuming h2osfc is impermeable (mm) (= 100._r8)

     ! additional constants
     real(r8), pointer :: f_sat        => null()         ! volumetric soil water defining top of water table or where production is allowed (=0.95)
     real(r8), pointer :: qflxlagd     => null()         ! days to lag qflx_surf_lag in the tropics (days) ( = 30._r8)
     real(r8), pointer :: highlatfact  => null()         ! multiple of qflxlagd for high latitudes	(= 2._r8)
     real(r8), pointer :: q10lakebase  => null()         ! (K) base temperature for lake CH4 production (= 298._r8)
     real(r8), pointer :: atmch4       => null()         ! Atmospheric CH4 mixing ratio to prescribe if not provided by the atmospheric model (= 1.7e-6_r8) (mol/mol)
     real(r8), pointer :: rob          => null()         ! ratio of root length to vertical depth ("root obliquity") (= 3._r8)
  end type CH4ParamsType
  type(CH4ParamsType), public ::  CH4ParamsInst
  !$acc declare create(CH4ParamsInst)
  type, public :: ch4_type
     real(r8), pointer :: ch4_prod_depth_sat_col     (:,:) =>null()! col CH4 production rate from methanotrophs (mol/m3/s) (nlevsoi)
     real(r8), pointer :: ch4_prod_depth_unsat_col   (:,:) =>null()! col CH4 production rate from methanotrophs (mol/m3/s) (nlevsoi)
     real(r8), pointer :: ch4_prod_depth_lake_col    (:,:) =>null()! col CH4 production rate from methanotrophs (mol/m3/s) (nlevsoi)
     real(r8), pointer :: ch4_oxid_depth_sat_col     (:,:) =>null()! col CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer :: ch4_oxid_depth_unsat_col   (:,:) =>null()! col CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer :: ch4_oxid_depth_lake_col    (:,:) =>null()! col CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer :: ch4_aere_depth_sat_col     (:,:) =>null()! col CH4 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer :: ch4_aere_depth_unsat_col   (:,:) =>null()! col CH4 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer :: ch4_tran_depth_sat_col     (:,:) =>null()! col CH4 loss rate via transpiration in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer :: ch4_tran_depth_unsat_col   (:,:) =>null()! col CH4 loss rate via transpiration in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer :: ch4_ebul_depth_sat_col     (:,:) =>null()! col CH4 loss rate via ebullition in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer :: ch4_ebul_depth_unsat_col   (:,:) =>null()! col CH4 loss rate via ebullition in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer :: ch4_ebul_total_sat_col     (:)   =>null()! col Total col CH4 ebullition (mol/m2/s)
     real(r8), pointer :: ch4_ebul_total_unsat_col   (:)   =>null()! col Total col CH4 ebullition (mol/m2/s)
     real(r8), pointer :: ch4_surf_aere_sat_col      (:)   =>null()! col CH4 aerenchyma flux to atmosphere (after oxidation) (mol/m2/s)
     real(r8), pointer :: ch4_surf_aere_unsat_col    (:)   =>null()! col CH4 aerenchyma flux to atmosphere (after oxidation) (mol/m2/s)
     real(r8), pointer :: ch4_surf_ebul_sat_col      (:)   =>null()! col CH4 ebullition flux to atmosphere (after oxidation) (mol/m2/s)
     real(r8), pointer :: ch4_surf_ebul_unsat_col    (:)   =>null()! col CH4 ebullition flux to atmosphere (after oxidation) (mol/m2/s)
     real(r8), pointer :: ch4_surf_ebul_lake_col     (:)   =>null()! col CH4 ebullition flux to atmosphere (after oxidation) (mol/m2/s)
     real(r8), pointer :: co2_aere_depth_sat_col     (:,:) =>null()! col CO2 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer :: co2_aere_depth_unsat_col   (:,:) =>null()! col CO2 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer :: o2_oxid_depth_sat_col      (:,:) =>null()! col O2 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer :: o2_oxid_depth_unsat_col    (:,:) =>null()! col O2 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer :: o2_aere_depth_sat_col      (:,:) =>null()! col O2 gain rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer :: o2_aere_depth_unsat_col    (:,:) =>null()! col O2 gain rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer :: co2_decomp_depth_sat_col   (:,:) =>null()! col CO2 production during decomposition in each soil layer (nlevsoi) (mol/m3/s)
     real(r8), pointer :: co2_decomp_depth_unsat_col (:,:) =>null()! col CO2 production during decomposition in each soil layer (nlevsoi) (mol/m3/s)
     real(r8), pointer :: co2_oxid_depth_sat_col     (:,:) =>null()! col CO2 production rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer :: co2_oxid_depth_unsat_col   (:,:) =>null()! col CO2 production rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer :: conc_o2_lake_col           (:,:) =>null()! col O2 conc in each soil layer (mol/m3) (nlevsoi)
     real(r8), pointer :: conc_ch4_sat_col           (:,:) =>null()! col CH4 conc in each soil layer (mol/m3) (nlevsoi)
     real(r8), pointer :: conc_ch4_unsat_col         (:,:) =>null()! col CH4 conc in each soil layer (mol/m3) (nlevsoi)
     real(r8), pointer :: conc_ch4_lake_col          (:,:) =>null()! col CH4 conc in each soil layer (mol/m3) (nlevsoi)
     real(r8), pointer :: ch4_surf_diff_sat_col      (:)   =>null()! col CH4 surface flux (mol/m2/s)
     real(r8), pointer :: ch4_surf_diff_unsat_col    (:)   =>null()! col CH4 surface flux (mol/m2/s)
     real(r8), pointer :: ch4_surf_diff_lake_col     (:)   =>null()! col CH4 surface flux (mol/m2/s)
     real(r8), pointer :: ch4_dfsat_flux_col         (:)   =>null()! col CH4 flux to atm due to decreasing fsat (kg C/m^2/s) [+]

     real(r8), pointer :: zwt_ch4_unsat_col          (:)   =>null()! col depth of water table for unsaturated fraction (m)
     real(r8), pointer :: fsat_bef_col               (:)   =>null()! col fsat from previous timestep
     real(r8), pointer :: lake_soilc_col             (:,:) =>null()! col total soil organic matter found in level (g C / m^3) (nlevsoi)
     real(r8), pointer :: totcolch4_col              (:)   =>null()! col total methane found in soil col (g C / m^2)
     real(r8), pointer :: annsum_counter_col         (:)   =>null()! col seconds since last annual accumulator turnover
     real(r8), pointer :: tempavg_somhr_col          (:)   =>null()! col temporary average SOM heterotrophic resp. (gC/m2/s)
     real(r8), pointer :: annavg_somhr_col           (:)   =>null()! col annual average SOM heterotrophic resp. (gC/m2/s)
     real(r8), pointer :: tempavg_finrw_col          (:)   =>null()! col respiration-weighted annual average of finundated
     real(r8), pointer :: annavg_finrw_col           (:)   =>null()! col respiration-weighted annual average of finundated
     real(r8), pointer :: sif_col                    (:)   =>null()! col (unitless) ratio applied to sat. prod. to account for seasonal inundation
     real(r8), pointer :: ch4stress_unsat_col        (:,:) =>null()! col Ratio of methane available to the total per-timestep methane sinks (nlevsoi)
     real(r8), pointer :: ch4stress_sat_col          (:,:) =>null()! col Ratio of methane available to the total per-timestep methane sinks (nlevsoi)
     real(r8), pointer :: qflx_surf_lag_col          (:)   =>null()! col time-lagged surface runoff (mm H2O /s)
     real(r8), pointer :: finundated_lag_col         (:)   =>null()! col time-lagged fractional inundated area
     real(r8), pointer :: layer_sat_lag_col          (:,:) =>null()! col Lagged saturation status of soil layer in the unsaturated zone (1 = sat)
     real(r8), pointer :: zwt0_col                   (:)   =>null()! col coefficient for determining finundated (m)
     real(r8), pointer :: f0_col                     (:)   =>null()! col maximum inundated fraction for a gridcell (for methane code)
     real(r8), pointer :: p3_col                     (:)   =>null()! col coefficient for determining finundated (m)
     real(r8), pointer :: pH_col                     (:)   =>null()! col pH values for methane production
    !
     real(r8), pointer :: c_atm_grc                  (:,:) =>null()! grc atmospheric conc of CH4, O2, CO2 (mol/m3)
     real(r8), pointer :: ch4co2f_grc                (:)   =>null()! grc CO2 production from CH4 oxidation (g C/m**2/s)
     real(r8), pointer :: ch4prodg_grc               (:)   =>null()! grc average CH4 production (g C/m^2/s)
     !
     real(r8), pointer :: finundated_col             (:)  =>null() ! col fractional inundated area (excluding dedicated wetland cols)
     real(r8), pointer :: o2stress_unsat_col         (:,:)=>null() ! col Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
     real(r8), pointer :: o2stress_sat_col           (:,:)=>null() ! col Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
     real(r8), pointer :: conc_o2_sat_col            (:,:)=>null() ! col O2 conc in each soil layer (mol/m3) (nlevsoi)
     real(r8), pointer :: conc_o2_unsat_col          (:,:)=>null() ! col O2 conc in each soil layer (mol/m3) (nlevsoi)
     real(r8), pointer :: o2_decomp_depth_sat_col    (:,:)=>null() ! col O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
     real(r8), pointer :: o2_decomp_depth_unsat_col  (:,:)=>null() ! col O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
     real(r8), pointer :: ch4_surf_flux_tot_col      (:)  =>null() ! col CH4 surface flux (to atm) (kg C/m**2/s)

     real(r8), pointer :: grnd_ch4_cond_patch        (:) => null()  ! patch tracer conductance for boundary layer [m/s]
     real(r8), pointer :: grnd_ch4_cond_col          (:) => null()  ! col tracer conductance for boundary layer [m/s]

   contains

     procedure, public  :: Init
     procedure, public :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold
     procedure, public  :: Restart

  end type ch4_type
contains 

  !------------------------------------------------------------------------
  subroutine Init( this, bounds,cellorg_col )

    class(ch4_type)               :: this
    type(bounds_type), intent(in) :: bounds
    real(r8)         , intent(in) :: cellorg_col (bounds%begc:, 1:)

    call this%InitAllocate (bounds)
    if (use_lch4) then
        call this%InitHistory (bounds)
        call this%InitCold (bounds, cellorg_col)
    end if

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Allocate module variables and data structures
    !
    ! !USES:
    use shr_infnan_mod, only: nan => shr_infnan_nan, assignment(=)
    use elm_varpar    , only: nlevgrnd
    !
    ! !ARGUMENTS:
    class(ch4_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: begp, endp
    integer  :: begc, endc
    integer  :: begg, endg
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc
    begg = bounds%begg; endg = bounds%endg

    allocate(this%ch4_prod_depth_sat_col     (begc:endc,1:nlevgrnd)) ;  this%ch4_prod_depth_sat_col     (:,:) = nan
    allocate(this%ch4_prod_depth_unsat_col   (begc:endc,1:nlevgrnd)) ;  this%ch4_prod_depth_unsat_col   (:,:) = nan
    allocate(this%ch4_prod_depth_lake_col    (begc:endc,1:nlevgrnd)) ;  this%ch4_prod_depth_lake_col    (:,:) = nan
    allocate(this%ch4_oxid_depth_sat_col     (begc:endc,1:nlevgrnd)) ;  this%ch4_oxid_depth_sat_col     (:,:) = nan
    allocate(this%ch4_oxid_depth_unsat_col   (begc:endc,1:nlevgrnd)) ;  this%ch4_oxid_depth_unsat_col   (:,:) = nan
    allocate(this%ch4_oxid_depth_lake_col    (begc:endc,1:nlevgrnd)) ;  this%ch4_oxid_depth_lake_col    (:,:) = nan
    allocate(this%o2_oxid_depth_sat_col      (begc:endc,1:nlevgrnd)) ;  this%o2_oxid_depth_sat_col      (:,:) = nan
    allocate(this%o2_oxid_depth_unsat_col    (begc:endc,1:nlevgrnd)) ;  this%o2_oxid_depth_unsat_col    (:,:) = nan
    allocate(this%o2_aere_depth_sat_col      (begc:endc,1:nlevgrnd)) ;  this%o2_aere_depth_sat_col      (:,:) = nan
    allocate(this%o2_aere_depth_unsat_col    (begc:endc,1:nlevgrnd)) ;  this%o2_aere_depth_unsat_col    (:,:) = nan
    allocate(this%co2_decomp_depth_sat_col   (begc:endc,1:nlevgrnd)) ;  this%co2_decomp_depth_sat_col   (:,:) = nan
    allocate(this%co2_decomp_depth_unsat_col (begc:endc,1:nlevgrnd)) ;  this%co2_decomp_depth_unsat_col (:,:) = nan
    allocate(this%co2_oxid_depth_sat_col     (begc:endc,1:nlevgrnd)) ;  this%co2_oxid_depth_sat_col     (:,:) = nan
    allocate(this%co2_oxid_depth_unsat_col   (begc:endc,1:nlevgrnd)) ;  this%co2_oxid_depth_unsat_col   (:,:) = nan
    allocate(this%ch4_aere_depth_sat_col     (begc:endc,1:nlevgrnd)) ;  this%ch4_aere_depth_sat_col     (:,:) = nan
    allocate(this%ch4_aere_depth_unsat_col   (begc:endc,1:nlevgrnd)) ;  this%ch4_aere_depth_unsat_col   (:,:) = nan
    allocate(this%ch4_tran_depth_sat_col     (begc:endc,1:nlevgrnd)) ;  this%ch4_tran_depth_sat_col     (:,:) = nan
    allocate(this%ch4_tran_depth_unsat_col   (begc:endc,1:nlevgrnd)) ;  this%ch4_tran_depth_unsat_col   (:,:) = nan
    allocate(this%co2_aere_depth_sat_col     (begc:endc,1:nlevgrnd)) ;  this%co2_aere_depth_sat_col     (:,:) = nan
    allocate(this%co2_aere_depth_unsat_col   (begc:endc,1:nlevgrnd)) ;  this%co2_aere_depth_unsat_col   (:,:) = nan
    allocate(this%ch4_surf_aere_sat_col      (begc:endc))            ;  this%ch4_surf_aere_sat_col      (:)   = nan
    allocate(this%ch4_surf_aere_unsat_col    (begc:endc))            ;  this%ch4_surf_aere_unsat_col    (:)   = nan
    allocate(this%ch4_ebul_depth_sat_col     (begc:endc,1:nlevgrnd)) ;  this%ch4_ebul_depth_sat_col     (:,:) = nan
    allocate(this%ch4_ebul_depth_unsat_col   (begc:endc,1:nlevgrnd)) ;  this%ch4_ebul_depth_unsat_col   (:,:) = nan
    allocate(this%ch4_ebul_total_sat_col     (begc:endc))            ;  this%ch4_ebul_total_sat_col     (:)   = nan
    allocate(this%ch4_ebul_total_unsat_col   (begc:endc))            ;  this%ch4_ebul_total_unsat_col   (:)   = nan
    allocate(this%ch4_surf_ebul_sat_col      (begc:endc))            ;  this%ch4_surf_ebul_sat_col      (:)   = nan
    allocate(this%ch4_surf_ebul_unsat_col    (begc:endc))            ;  this%ch4_surf_ebul_unsat_col    (:)   = nan
    allocate(this%ch4_surf_ebul_lake_col     (begc:endc))            ;  this%ch4_surf_ebul_lake_col     (:)   = nan
    allocate(this%conc_ch4_sat_col           (begc:endc,1:nlevgrnd)) ;  this%conc_ch4_sat_col           (:,:) = spval ! detect file input
    allocate(this%conc_ch4_unsat_col         (begc:endc,1:nlevgrnd)) ;  this%conc_ch4_unsat_col         (:,:) = spval ! detect file input
    allocate(this%conc_ch4_lake_col          (begc:endc,1:nlevgrnd)) ;  this%conc_ch4_lake_col          (:,:) = nan
    allocate(this%ch4_surf_diff_sat_col      (begc:endc))            ;  this%ch4_surf_diff_sat_col      (:)   = nan
    allocate(this%ch4_surf_diff_unsat_col    (begc:endc))            ;  this%ch4_surf_diff_unsat_col    (:)   = nan
    allocate(this%ch4_surf_diff_lake_col     (begc:endc))            ;  this%ch4_surf_diff_lake_col     (:)   = nan
    allocate(this%conc_o2_lake_col           (begc:endc,1:nlevgrnd)) ;  this%conc_o2_lake_col           (:,:) = nan
    allocate(this%ch4_dfsat_flux_col         (begc:endc))            ;  this%ch4_dfsat_flux_col         (:)   = nan
    allocate(this%zwt_ch4_unsat_col          (begc:endc))            ;  this%zwt_ch4_unsat_col          (:)   = nan
    allocate(this%fsat_bef_col               (begc:endc))            ;  this%fsat_bef_col               (:)   = nan
    allocate(this%lake_soilc_col             (begc:endc,1:nlevgrnd)) ;  this%lake_soilc_col             (:,:) = spval !first time-step
    allocate(this%totcolch4_col              (begc:endc))            ;  this%totcolch4_col              (:)   = nan
    allocate(this%annsum_counter_col         (begc:endc))            ;  this%annsum_counter_col         (:)   = nan
    allocate(this%tempavg_somhr_col          (begc:endc))            ;  this%tempavg_somhr_col          (:)   = nan
    allocate(this%annavg_somhr_col           (begc:endc))            ;  this%annavg_somhr_col           (:)   = nan
    allocate(this%tempavg_finrw_col          (begc:endc))            ;  this%tempavg_finrw_col          (:)   = nan
    allocate(this%annavg_finrw_col           (begc:endc))            ;  this%annavg_finrw_col           (:)   = nan
    allocate(this%sif_col                    (begc:endc))            ;  this%sif_col                    (:)   = nan
    allocate(this%ch4stress_unsat_col        (begc:endc,1:nlevgrnd)) ;  this%ch4stress_unsat_col        (:,:) = nan
    allocate(this%ch4stress_sat_col          (begc:endc,1:nlevgrnd)) ;  this%ch4stress_sat_col          (:,:) = nan
    allocate(this%qflx_surf_lag_col          (begc:endc))            ;  this%qflx_surf_lag_col          (:)   = nan
    allocate(this%finundated_lag_col         (begc:endc))            ;  this%finundated_lag_col         (:)   = nan
    allocate(this%layer_sat_lag_col          (begc:endc,1:nlevgrnd)) ;  this%layer_sat_lag_col          (:,:) = nan
    allocate(this%zwt0_col                   (begc:endc))            ;  this%zwt0_col                   (:)   = nan
    allocate(this%f0_col                     (begc:endc))            ;  this%f0_col                     (:)   = nan
    allocate(this%p3_col                     (begc:endc))            ;  this%p3_col                     (:)   = nan
    allocate(this%pH_col                     (begc:endc))            ;  this%pH_col                     (:)   = nan
    allocate(this%ch4_surf_flux_tot_col      (begc:endc))            ;  this%ch4_surf_flux_tot_col      (:)   = nan

    allocate(this%c_atm_grc                  (begg:endg,1:ngases))   ;  this%c_atm_grc                  (:,:) = nan
    allocate(this%ch4co2f_grc                (begg:endg))            ;  this%ch4co2f_grc                (:)   = nan
    allocate(this%ch4prodg_grc               (begg:endg))            ;  this%ch4prodg_grc               (:)   = nan

    allocate(this%finundated_col             (begc:endc))            ;  this%finundated_col             (:)   = nan
    allocate(this%o2stress_unsat_col         (begc:endc,1:nlevgrnd)) ;  this%o2stress_unsat_col         (:,:) = nan
    allocate(this%o2stress_sat_col           (begc:endc,1:nlevgrnd)) ;  this%o2stress_sat_col           (:,:) = nan
    allocate(this%conc_o2_sat_col            (begc:endc,1:nlevgrnd)) ;  this%conc_o2_sat_col            (:,:) = nan
    allocate(this%conc_o2_unsat_col          (begc:endc,1:nlevgrnd)) ;  this%conc_o2_unsat_col          (:,:) = nan
    allocate(this%o2_decomp_depth_sat_col    (begc:endc,1:nlevgrnd)) ;  this%o2_decomp_depth_sat_col    (:,:) = nan
    allocate(this%o2_decomp_depth_unsat_col  (begc:endc,1:nlevgrnd)) ;  this%o2_decomp_depth_unsat_col  (:,:) = nan

    allocate(this%grnd_ch4_cond_patch        (begp:endp)) ;  this%grnd_ch4_cond_patch(:) = nan
    allocate(this%grnd_ch4_cond_col          (begc:endc)) ;  this%grnd_ch4_cond_col  (:) = nan

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !USES:
    use elm_varpar , only : nlevgrnd, nlevdecomp
    use elm_varctl , only : hist_wrtch4diag
    use histFileMod, only : hist_addfld1d, hist_addfld2d, hist_addfld_decomp
    use CH4varcon  , only : allowlakeprod
    !
    ! !ARGUMENTS:
    class(ch4_type) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    character(8)  :: vr_suffix
    character(10) :: active
    integer       :: begc,endc
    integer       :: begg,endg
    !---------------------------------------------------------------------

    begc = bounds%begc; endc = bounds%endc
    begg = bounds%begg; endg = bounds%endg

    if (nlevdecomp > 1) then
       vr_suffix = "_vr"
    else
       vr_suffix = ""
    endif

    if (hist_wrtch4diag) then
       active = "active"
    else
       active = "inactive"
    end if

    this%finundated_col(begc:endc) = spval
    call hist_addfld1d (fname='FINUNDATED', units='1', &
         avgflag='A', long_name='fractional inundated area of vegetated columns', &
         ptr_col=this%finundated_col)

    this%finundated_lag_col(begc:endc) = spval
    call hist_addfld1d (fname='FINUNDATED_LAG', units='1',  &
         avgflag='A', long_name='time-lagged inundated fraction of vegetated columns', &
         ptr_col=this%finundated_lag_col)

    this%ch4_surf_diff_sat_col(begc:endc) = spval
    call hist_addfld1d (fname='CH4_SURF_DIFF_SAT', units='mol/m2/s',  &
         avgflag='A', long_name='diffusive surface CH4 flux for inundated / lake area; (+ to atm)', &
         ptr_col=this%ch4_surf_diff_sat_col)

    this%ch4_surf_diff_unsat_col(begc:endc) = spval
    call hist_addfld1d (fname='CH4_SURF_DIFF_UNSAT', units='mol/m2/s',  &
         avgflag='A', long_name='diffusive surface CH4 flux for non-inundated area; (+ to atm)', &
         ptr_col=this%ch4_surf_diff_unsat_col)

    this%ch4_ebul_total_sat_col(begc:endc) = spval
    call hist_addfld1d (fname='CH4_EBUL_TOTAL_SAT', units='mol/m2/s',  &
         avgflag='A', long_name='ebullition surface CH4 flux; (+ to atm)', &
         ptr_col=this%ch4_ebul_total_sat_col, default='inactive')

    this%ch4_ebul_total_unsat_col(begc:endc) = spval
    call hist_addfld1d (fname='CH4_EBUL_TOTAL_UNSAT', units='mol/m2/s',  &
         avgflag='A', long_name='ebullition surface CH4 flux; (+ to atm)', &
         ptr_col=this%ch4_ebul_total_unsat_col, default='inactive')

    this%ch4_surf_ebul_sat_col(begc:endc) = spval
    call hist_addfld1d (fname='CH4_SURF_EBUL_SAT', units='mol/m2/s',  &
         avgflag='A', long_name='ebullition surface CH4 flux for inundated / lake area; (+ to atm)', &
         ptr_col=this%ch4_surf_ebul_sat_col)

    this%ch4_surf_ebul_unsat_col(begc:endc) = spval
    call hist_addfld1d (fname='CH4_SURF_EBUL_UNSAT', units='mol/m2/s',  &
         avgflag='A', long_name='ebullition surface CH4 flux for non-inundated area; (+ to atm)', &
         ptr_col=this%ch4_surf_ebul_unsat_col)

    this%ch4_surf_aere_sat_col(begc:endc) = spval
    call hist_addfld1d (fname='CH4_SURF_AERE_SAT', units='mol/m2/s',  &
         avgflag='A', long_name='aerenchyma surface CH4 flux for inundated area; (+ to atm)', &
         ptr_col=this%ch4_surf_aere_sat_col)

    this%ch4_surf_aere_unsat_col(begc:endc) = spval
    call hist_addfld1d (fname='CH4_SURF_AERE_UNSAT', units='mol/m2/s',  &
         avgflag='A', long_name='aerenchyma surface CH4 flux for non-inundated area; (+ to atm)', &
         ptr_col=this%ch4_surf_aere_unsat_col)

    this%totcolch4_col(begc:endc) = spval
    call hist_addfld1d (fname='TOTCOLCH4', units='gC/m2',  &
         avgflag='A', long_name='total belowground CH4, (0 for non-lake special landunits)', &
         ptr_col=this%totcolch4_col)

    this%conc_ch4_sat_col(begc:endc,1:nlevgrnd) = spval
    call hist_addfld2d (fname='CONC_CH4_SAT', units='mol/m3', type2d='levgrnd', &
         avgflag='A', long_name='CH4 soil Concentration for inundated / lake area', &
         ptr_col=this%conc_ch4_sat_col)

    this%conc_ch4_unsat_col(begc:endc,1:nlevgrnd) = spval
    call hist_addfld2d (fname='CONC_CH4_UNSAT', units='mol/m3', type2d='levgrnd', &
         avgflag='A', long_name='CH4 soil Concentration for non-inundated area', &
         ptr_col=this%conc_ch4_unsat_col)

    if (hist_wrtch4diag) then
       this%ch4_prod_depth_sat_col(begc:endc,1:nlevgrnd) = spval
       call hist_addfld2d (fname='CH4_PROD_DEPTH_SAT', units='mol/m3/s', type2d='levgrnd', &
            avgflag='A', long_name='CH4 soil production for inundated / lake area', &
            ptr_col=this%ch4_prod_depth_sat_col)
    end if

    if (hist_wrtch4diag) then
       this%ch4_prod_depth_unsat_col(begc:endc,1:nlevgrnd) = spval
       call hist_addfld2d (fname='CH4_PROD_DEPTH_UNSAT', units='mol/m3/s', type2d='levgrnd', &
            avgflag='A', long_name='CH4 soil production for non-inundated area', &
            ptr_col=this%ch4_prod_depth_unsat_col)
    end if

    if (hist_wrtch4diag) then
       this%ch4_oxid_depth_sat_col(begc:endc,1:nlevgrnd) = spval
       call hist_addfld2d (fname='CH4_OXID_DEPTH_SAT', units='mol/m3/s', type2d='levgrnd', &
            avgflag='A', long_name='CH4 soil oxidation for inundated / lake area', &
            ptr_col=this%ch4_oxid_depth_sat_col)
    end if

    if (hist_wrtch4diag) then
       this%ch4_oxid_depth_unsat_col(begc:endc,1:nlevgrnd) = spval
       call hist_addfld2d (fname='CH4_OXID_DEPTH_UNSAT', units='mol/m3/s', type2d='levgrnd', &
            avgflag='A', long_name='CH4 soil oxidation for non-inundated area', &
            ptr_col=this%ch4_oxid_depth_unsat_col)
    end if

    if (hist_wrtch4diag) then
       this%ch4_aere_depth_sat_col(begc:endc,1:nlevgrnd) = spval
       call hist_addfld2d (fname='CH4_AERE_DEPTH_SAT', units='mol/m3/s', type2d='levgrnd', &
            avgflag='A', long_name='CH4 soil aerenchyma loss for inundated / lake area '// &
            ' (including transpiration flux if activated)', &
            ptr_col=this%ch4_aere_depth_sat_col)
    end if

    if (hist_wrtch4diag) then
       this%ch4_aere_depth_unsat_col(begc:endc,1:nlevgrnd) = spval
       call hist_addfld2d (fname='CH4_AERE_DEPTH_UNSAT', units='mol/m3/s', type2d='levgrnd', &
            avgflag='A', long_name='CH4 soil aerenchyma loss for non-inundated area '// &
            ' (including transpiration flux if activated)', &
            ptr_col=this%ch4_aere_depth_unsat_col)
    end if

    if (hist_wrtch4diag) then
       this%o2_aere_depth_sat_col(begc:endc,1:nlevgrnd) = spval
       call hist_addfld2d (fname='O2_AERE_DEPTH_SAT', units='mol/m3/s', type2d='levgrnd', &
            avgflag='A', long_name='O2 aerenchyma diffusion into soil for inundated / lake area', &
            ptr_col=this%o2_aere_depth_sat_col)
    end if

    if (hist_wrtch4diag) then
       this%o2_aere_depth_unsat_col(begc:endc,1:nlevgrnd) = spval
       call hist_addfld2d (fname='O2_AERE_DEPTH_UNSAT', units='mol/m3/s', type2d='levgrnd', &
            avgflag='A', long_name='O2 aerenchyma diffusion into soil for non-inundated area', &
            ptr_col=this%o2_aere_depth_unsat_col)
    end if

    if (hist_wrtch4diag) then
       call hist_addfld2d (fname='O2_DECOMP_DEPTH_SAT', units='mol/m3/s', type2d='levgrnd', &
            avgflag='A', long_name='O2 consumption from HR and AR for inundated / lake area', &
            ptr_col=this%o2_decomp_depth_sat_col)
    end if

    this%o2_decomp_depth_unsat_col(begc:endc,1:nlevgrnd) = spval
    if (hist_wrtch4diag) then
       this%o2_decomp_depth_unsat_col(begc:endc,1:nlevgrnd) = spval
       call hist_addfld2d (fname='O2_DECOMP_DEPTH_UNSAT', units='mol/m3/s', type2d='levgrnd', &
            avgflag='A', long_name='O2 consumption from HR and AR for non-inundated area', &
            ptr_col=this%o2_decomp_depth_unsat_col)
    else
       call hist_addfld2d (fname='o2_decomp_depth_unsat', units='mol/m3/2', type2d='levgrnd', &
            avgflag='A', long_name='o2_decomp_depth_unsat', &
            ptr_col=this%o2_decomp_depth_unsat_col)
    end if

    if (hist_wrtch4diag) then
       this%ch4_tran_depth_sat_col(begc:endc,1:nlevgrnd) = spval
       call hist_addfld2d (fname='CH4_TRAN_DEPTH_SAT', units='mol/m3/s', type2d='levgrnd', &
            avgflag='A', long_name='CH4 soil loss from transpiration for inundated / lake area', &
            ptr_col=this%ch4_tran_depth_sat_col)
    end if

    if (hist_wrtch4diag) then
       this%ch4_tran_depth_unsat_col(begc:endc,1:nlevgrnd) = spval
       call hist_addfld2d (fname='CH4_TRAN_DEPTH_UNSAT', units='mol/m3/s', type2d='levgrnd', &
            avgflag='A', long_name='CH4 soil loss from transpiration for non-inundated area', &
            ptr_col=this%ch4_tran_depth_unsat_col)
    end if

    if (hist_wrtch4diag) then
       this%ch4_ebul_depth_sat_col(begc:endc,1:nlevgrnd) = spval
       call hist_addfld2d (fname='CH4_EBUL_DEPTH_SAT', units='mol/m3/s', type2d='levgrnd', &
            avgflag='A', long_name='CH4 soil ebullition for inundated / lake area', &
            ptr_col=this%ch4_ebul_depth_sat_col)
    end if

    if (hist_wrtch4diag) then
       this%ch4_ebul_depth_unsat_col(begc:endc,1:nlevgrnd) = spval
       call hist_addfld2d (fname='CH4_EBUL_DEPTH_UNSAT', units='mol/m3/s', type2d='levgrnd', &
            avgflag='A', long_name='CH4 soil ebullition for non-inundated area', &
            ptr_col=this%ch4_ebul_depth_unsat_col)
    end if

    if (hist_wrtch4diag) then
       this%o2stress_sat_col(begc:endc,1:nlevgrnd) = spval
       call hist_addfld2d (fname='O2STRESS_SAT', units='1', type2d='levgrnd',  &
            avgflag='A', long_name='Ratio of oxygen available to demanded for non-inundated area', &
            ptr_col=this%o2stress_sat_col)
    end if

    if (hist_wrtch4diag) then
       this%o2stress_unsat_col(begc:endc,1:nlevgrnd) = spval
       call hist_addfld2d (fname='O2STRESS_UNSAT', units='1', type2d='levgrnd',  &
            avgflag='A', long_name='Ratio of oxygen available to demanded for inundated / lake area', &
            ptr_col=this%o2stress_unsat_col)
    end if

    if (hist_wrtch4diag) then
       this%ch4stress_unsat_col(begc:endc,1:nlevgrnd) = spval
       call hist_addfld2d (fname='CH4STRESS_UNSAT', units='1', type2d='levgrnd',  &
            avgflag='A', long_name='Ratio of methane available to total potential sink for inundated / lake area', &
            ptr_col=this%ch4stress_unsat_col)
    end if

    if (hist_wrtch4diag) then
       this%ch4stress_sat_col(begc:endc,1:nlevgrnd) = spval
       call hist_addfld2d (fname='CH4STRESS_SAT', units='1', type2d='levgrnd',  &
            avgflag='A', long_name='Ratio of methane available to total potential sink for non-inundated area', &
            ptr_col=this%ch4stress_sat_col)
    end if

    if (hist_wrtch4diag .and. allowlakeprod) then
       this%ch4_prod_depth_sat_col(begc:endc,1:nlevgrnd) = spval
       call hist_addfld2d (fname='CH4_PROD_DEPTH_LAKE', units='mol/m3/s', type2d='levgrnd', &
            avgflag='A', long_name='CH4 production in each soil layer, lake col. only', &
            ptr_col=this%ch4_prod_depth_sat_col)
    end if

    if (hist_wrtch4diag .and. allowlakeprod) then
       this%conc_ch4_sat_col(begc:endc,1:nlevgrnd) = spval
       call hist_addfld2d (fname='CONC_CH4_LAKE', units='mol/m3', type2d='levgrnd', &
            avgflag='A', long_name='CH4 Concentration each soil layer, lake col. only', &
            ptr_col=this%conc_ch4_sat_col)
    end if

    if (hist_wrtch4diag .and. allowlakeprod) then
       this%conc_o2_sat_col(begc:endc,1:nlevgrnd) = spval
       call hist_addfld2d (fname='CONC_O2_LAKE', units='mol/m3', type2d='levgrnd', &
            avgflag='A', long_name='O2 Concentration each soil layer, lake col. only', &
            ptr_col=this%conc_o2_sat_col)
    end if

    if (hist_wrtch4diag .and. allowlakeprod) then
       this%ch4_surf_diff_sat_col(begc:endc) = spval
       call hist_addfld1d (fname='CH4_SURF_DIFF_LAKE', units='mol/m2/s',  &
            avgflag='A', long_name='diffusive surface CH4 flux, lake col. only (+ to atm)', &
            ptr_col=this%ch4_surf_diff_sat_col)
    end if

    if (hist_wrtch4diag .and. allowlakeprod) then
       this%ch4_surf_ebul_sat_col(begc:endc) = spval
       call hist_addfld1d (fname='CH4_SURF_EBUL_LAKE', units='mol/m2/s',  &
            avgflag='A', long_name='ebullition surface CH4 flux, lake col. only (+ to atm)', &
            ptr_col=this%ch4_surf_ebul_sat_col)
    end if

    if (hist_wrtch4diag .and. allowlakeprod) then
       this%ch4_oxid_depth_sat_col(begc:endc,1:nlevgrnd) = spval
       call hist_addfld2d (fname='CH4_OXID_DEPTH_LAKE', units='mol/m2/s', type2d='levgrnd',  &
            avgflag='A', long_name='CH4 oxidation in each soil layer, lake col. only', &
            ptr_col=this%ch4_oxid_depth_sat_col)
    end if

    if (hist_wrtch4diag) then
       this%layer_sat_lag_col(begc:endc,1:nlevgrnd) = spval
       call hist_addfld2d (fname='LAYER_SAT_LAG', units='1', type2d='levgrnd',  &
            avgflag='A', long_name='lagged saturation status of layer in unsat. zone', &
            ptr_col=this%layer_sat_lag_col)
    end if

    if (hist_wrtch4diag) then
       this%annavg_finrw_col(begc:endc) = spval
       call hist_addfld1d (fname='ANNAVG_FINRW', units='1',  &
            avgflag='A', long_name='annual average respiration-weighted FINUNDATED', &
            ptr_col=this%annavg_finrw_col)
    end if

    if (hist_wrtch4diag) then
       this%sif_col(begc:endc) = spval
       call hist_addfld1d (fname='SIF', units='1',  &
            avgflag='A', long_name='seasonal inundation factor calculated for sat. CH4 prod. (non-lake)', &
            ptr_col=this%sif_col)
    end if

    this%conc_o2_sat_col(begc:endc,1:nlevgrnd) = spval
    call hist_addfld2d (fname='CONC_O2_SAT', units='mol/m3', type2d='levgrnd', &
         avgflag='A', long_name='O2 soil Concentration for inundated / lake area', &
         ptr_col=this%conc_o2_sat_col)

    this%conc_o2_unsat_col(begc:endc,1:nlevgrnd) = spval
    call hist_addfld2d (fname='CONC_O2_UNSAT', units='mol/m3', type2d='levgrnd', &
         avgflag='A', long_name='O2 soil Concentration for non-inundated area', &
         ptr_col=this%conc_o2_unsat_col)

    this%ch4co2f_grc(begg:endg) = spval
    call hist_addfld1d (fname='FCH4TOCO2', units='gC/m2/s', &
         avgflag='A', long_name='Gridcell oxidation of CH4 to CO2', &
         ptr_lnd=this%ch4co2f_grc)

    this%ch4prodg_grc(begg:endg) = spval
    call hist_addfld1d (fname='CH4PROD', units='gC/m2/s', &
         avgflag='A', long_name='Gridcell total production of CH4', &
         ptr_lnd=this%ch4prodg_grc)

    this%ch4_dfsat_flux_col(begc:endc) = spval
    call hist_addfld1d (fname='FCH4_DFSAT', units='kgC/m2/s',  &
         avgflag='A', long_name='CH4 additional flux due to changing fsat, vegetated landunits only', &
         ptr_col=this%ch4_dfsat_flux_col)

    this%zwt_ch4_unsat_col(begc:endc) = spval
    call hist_addfld1d (fname='ZWT_CH4_UNSAT', units='m',  &
         avgflag='A', long_name='depth of water table for methane production used in non-inundated area', &
         ptr_col=this%zwt_ch4_unsat_col)

    this%qflx_surf_lag_col(begc:endc) = spval
    call hist_addfld1d (fname='QOVER_LAG', units='mm/s',  &
         avgflag='A', long_name='time-lagged surface runoff for soil columns', &
         ptr_col=this%qflx_surf_lag_col)

    if (allowlakeprod) then
       this%lake_soilc_col(begc:endc,1:nlevgrnd) = spval
       call hist_addfld2d (fname='LAKE_SOILC', units='gC/m3', type2d='levgrnd', &
            avgflag='A', long_name='Soil carbon under lakes', &
            ptr_col=this%lake_soilc_col)
    end if

    this%grnd_ch4_cond_col(begc:endc) = spval
    call hist_addfld1d (fname='WTGQ', units='m/s',  &
         avgflag='A', long_name='surface tracer conductance', &
         ptr_col=this%grnd_ch4_cond_col)

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, cellorg_col)
    !
    ! !DESCRIPTION:
    ! - Sets cold start values for time varying values.
    ! Initializes the following time varying variables:
    ! conc_ch4_sat, conc_ch4_unsat, conc_o2_sat, conc_o2_unsat,
    ! lake_soilc, o2stress, finunduated
    ! - Sets variables for ch4 code that will not be input
    ! from restart/inic file.
    ! - Sets values for inactive CH4 columns to spval so that they will
    ! not be averaged in history file.
    !
    ! !USES:
    use shr_kind_mod    , only : r8 => shr_kind_r8
    use elm_varpar      , only : nlevsoi, nlevgrnd, nlevdecomp
    use landunit_varcon , only : istsoil, istdlak, istcrop
    use elm_varctl      , only : iulog, fsurdat
    use CH4varcon       , only : allowlakeprod, usephfact, fin_use_fsat
    use spmdMod         , only : masterproc
    use fileutils       , only : getfil
    use ncdio_pio
    !
    ! !ARGUMENTS:
    class(ch4_type) :: this
    type(bounds_type) , intent(in) :: bounds
    real(r8)          , intent(in) :: cellorg_col (bounds%begc:, 1:)
    !
    ! !LOCAL VARIABLES:
    integer               :: j ,g, l,c,p ! indices
    type(file_desc_t)     :: ncid        ! netcdf id
    real(r8)     ,pointer :: zwt0_in (:) ! read in - zwt0
    real(r8)     ,pointer :: f0_in (:)   ! read in - f0
    real(r8)     ,pointer :: p3_in (:)   ! read in - p3
    real(r8)     ,pointer :: pH_in (:)   ! read in - pH
    logical               :: readvar
    character(len=256)    :: locfn
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(cellorg_col) == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))

    !----------------------------------------
    ! Initialize time constant variables
    !----------------------------------------

    allocate(zwt0_in(bounds%begg:bounds%endg))
    allocate(f0_in(bounds%begg:bounds%endg))
    allocate(p3_in(bounds%begg:bounds%endg))
    if (usephfact) allocate(ph_in(bounds%begg:bounds%endg))

    ! Methane code parameters for finundated

    if (.not. fin_use_fsat) then
       call getfil (fsurdat, locfn, 0)
       call ncd_pio_openfile (ncid, locfn, 0)
       call ncd_io(ncid=ncid, varname='ZWT0', flag='read', data=zwt0_in, dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun(msg=' ERROR: Running with CH4 Model but ZWT0 not on surfdata file'//&
               errMsg(__FILE__, __LINE__))
       end if
       call ncd_io(ncid=ncid, varname='F0', flag='read', data=f0_in, dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun(msg=' ERROR: Running with CH4 Model but F0 not on surfdata file'//&
               errMsg(__FILE__, __LINE__))
       end if
       call ncd_io(ncid=ncid, varname='P3', flag='read', data=p3_in, dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun(msg=' ERROR: Running with CH4 Model but P3 not on surfdata file'//&
               errMsg(__FILE__, __LINE__))
       end if
    end if

    ! pH factor for methane model
    if (usephfact) then
       call ncd_io(ncid=ncid, varname='PH', flag='read', data=ph_in, dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun(msg=' ERROR: CH4 pH production factor activated in ch4par_in'//&
               'but pH is not on surfdata file'//errMsg(__FILE__, __LINE__))
       end if
    end if

    do c = bounds%begc, bounds%endc
       g = col_pp%gridcell(c)

       if (.not. fin_use_fsat) then
          this%zwt0_col(c) = zwt0_in(g)
          this%f0_col(c)   = f0_in(g)
          this%p3_col(c)   = p3_in(g)
       end if
       if (usephfact) this%pH_col(c) = pH_in(g)
    end do

    deallocate(zwt0_in, f0_in, p3_in)
    if (usephfact) deallocate(pH_in)

    !----------------------------------------
    ! Initialize time varying variables
    !----------------------------------------

    if ( masterproc ) write (iulog,*) 'Setting initial data to non-spun up values for CH4 Mod'

    do c = bounds%begc,bounds%endc

       ! To detect first time-step
       this%fsat_bef_col       (c) = spval
       this%annsum_counter_col (c) = spval
       this%totcolch4_col      (c) = spval

       ! To detect first year
       this%annavg_somhr_col(c)   = spval
       this%annavg_finrw_col(c)   = spval

       ! To detect file input
       this%qflx_surf_lag_col  (c)   = spval
       this%finundated_lag_col (c)   = spval
       this%layer_sat_lag_col  (c,:) = spval
       this%conc_ch4_sat_col   (c,:) = spval
       this%conc_ch4_unsat_col (c,:) = spval
       this%conc_o2_sat_col    (c,:) = spval
       this%conc_o2_unsat_col  (c,:) = spval
       this%o2stress_sat_col   (c,:) = spval
       this%o2stress_unsat_col (c,:) = spval
       this%ch4stress_sat_col  (c,:) = spval
       this%ch4stress_unsat_col(c,:) = spval
       this%lake_soilc_col     (c,:) = spval

       ! To detect first time-step for denitrification code
       this%o2_decomp_depth_unsat_col(c,:)= spval

       l = col_pp%landunit(c)
       if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then

          this%conc_ch4_sat_col   (c,1:nlevsoi) = 0._r8
          this%conc_ch4_unsat_col (c,1:nlevsoi) = 0._r8
          this%conc_o2_sat_col    (c,1:nlevsoi) = 0._r8
          this%conc_o2_unsat_col  (c,1:nlevsoi) = 0._r8
          this%o2stress_sat_col   (c,1:nlevsoi) = 1._r8
          this%o2stress_unsat_col (c,1:nlevsoi) = 1._r8
          this%layer_sat_lag_col  (c,1:nlevsoi) = 1._r8
          this%qflx_surf_lag_col  (c)           = 0._r8
          this%finundated_lag_col (c)           = 0._r8
          this%finundated_col(c)                = 0._r8
          ! finundated will be used to calculate soil decomposition if anoxia is used
          ! Note that finundated will be overwritten with this%fsat_bef_col upon reading
          ! a restart file - either in a continuation, branch or startup spun-up case

       else if (lun_pp%itype(l) == istdlak) then

          this%conc_ch4_sat_col(c,1:nlevsoi) = 0._r8
          this%conc_o2_sat_col (c,1:nlevsoi) = 0._r8
          this%lake_soilc_col  (c,1:nlevsoi) = 580._r8 * cellorg_col(c,1:nlevsoi)

       end if

       ! Set values for all columns equal  below nlevsoi

       this%conc_ch4_sat_col           (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%conc_ch4_unsat_col         (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%conc_o2_sat_col            (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%conc_o2_unsat_col          (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%lake_soilc_col             (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%o2stress_sat_col           (c,nlevsoi+1:nlevgrnd) = 1._r8
       this%o2stress_unsat_col         (c,nlevsoi+1:nlevgrnd) = 1._r8
       this%layer_sat_lag_col          (c,nlevsoi+1:nlevgrnd) = 1._r8
       this%ch4_prod_depth_sat_col     (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%ch4_prod_depth_unsat_col   (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%ch4_prod_depth_lake_col    (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%ch4_oxid_depth_sat_col     (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%ch4_oxid_depth_unsat_col   (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%ch4_oxid_depth_lake_col    (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%o2_oxid_depth_sat_col      (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%o2_oxid_depth_unsat_col    (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%o2_decomp_depth_sat_col    (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%o2_decomp_depth_unsat_col  (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%o2_aere_depth_sat_col      (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%o2_aere_depth_unsat_col    (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%co2_decomp_depth_sat_col   (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%co2_decomp_depth_unsat_col (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%co2_oxid_depth_sat_col     (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%co2_oxid_depth_unsat_col   (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%ch4_aere_depth_sat_col     (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%ch4_aere_depth_unsat_col   (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%ch4_tran_depth_sat_col     (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%ch4_tran_depth_unsat_col   (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%co2_aere_depth_sat_col     (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%co2_aere_depth_unsat_col   (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%ch4_ebul_depth_sat_col     (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%ch4_ebul_depth_unsat_col   (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%conc_ch4_lake_col          (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%conc_o2_lake_col           (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%ch4stress_unsat_col        (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%ch4stress_sat_col          (c,nlevsoi+1:nlevgrnd) = 0._r8

       if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then

          this%conc_ch4_lake_col       (c,:) = spval
          this%conc_o2_lake_col        (c,:) = spval
          this%ch4_surf_diff_lake_col  (c)   = spval
          this%ch4_surf_ebul_lake_col  (c)   = spval
          this%ch4_prod_depth_lake_col (c,:) = spval
          this%ch4_oxid_depth_lake_col (c,:) = spval

       else if (lun_pp%itype(l) == istdlak .and. allowlakeprod) then

          this%ch4_prod_depth_unsat_col   (c,:) = spval
          this%ch4_oxid_depth_unsat_col   (c,:) = spval
          this%o2_oxid_depth_unsat_col    (c,:) = spval
          this%o2_decomp_depth_unsat_col  (c,:) = spval
          this%o2_aere_depth_unsat_col    (c,:) = spval
          this%co2_decomp_depth_unsat_col (c,:) = spval
          this%co2_oxid_depth_unsat_col   (c,:) = spval
          this%ch4_aere_depth_unsat_col   (c,:) = spval
          this%ch4_tran_depth_unsat_col   (c,:) = spval
          this%co2_aere_depth_unsat_col   (c,:) = spval
          this%ch4_surf_aere_unsat_col    (c)   = spval
          this%ch4_ebul_depth_unsat_col   (c,:) = spval
          this%ch4_ebul_total_unsat_col   (c)   = spval
          this%ch4_surf_ebul_unsat_col    (c)   = spval
          this%ch4_surf_diff_unsat_col    (c)   = spval
          this%ch4_dfsat_flux_col         (c)   = spval
          this%zwt_ch4_unsat_col          (c)   = spval
          this%sif_col                    (c)   = spval
          this%o2stress_unsat_col         (c,:) = spval
          this%ch4stress_unsat_col        (c,:) = spval
          this%finundated_col             (c)   = spval

       else  ! Inactive CH4 columns

          this%ch4_prod_depth_sat_col     (c,:) = spval
          this%ch4_prod_depth_unsat_col   (c,:) = spval
          this%ch4_prod_depth_lake_col    (c,:) = spval
          this%ch4_oxid_depth_sat_col     (c,:) = spval
          this%ch4_oxid_depth_unsat_col   (c,:) = spval
          this%ch4_oxid_depth_lake_col    (c,:) = spval
          this%o2_oxid_depth_sat_col      (c,:) = spval
          this%o2_oxid_depth_unsat_col    (c,:) = spval
          this%o2_decomp_depth_sat_col    (c,:) = spval
          this%o2_decomp_depth_unsat_col  (c,:) = spval
          this%o2_aere_depth_sat_col      (c,:) = spval
          this%o2_aere_depth_unsat_col    (c,:) = spval
          this%co2_decomp_depth_sat_col   (c,:) = spval
          this%co2_decomp_depth_unsat_col (c,:) = spval
          this%co2_oxid_depth_sat_col     (c,:) = spval
          this%co2_oxid_depth_unsat_col   (c,:) = spval
          this%ch4_aere_depth_sat_col     (c,:) = spval
          this%ch4_aere_depth_unsat_col   (c,:) = spval
          this%ch4_tran_depth_sat_col     (c,:) = spval
          this%ch4_tran_depth_unsat_col   (c,:) = spval
          this%co2_aere_depth_sat_col     (c,:) = spval
          this%co2_aere_depth_unsat_col   (c,:) = spval
          this%ch4_surf_aere_sat_col      (c)   = spval
          this%ch4_surf_aere_unsat_col    (c)   = spval
          this%ch4_ebul_depth_sat_col     (c,:) = spval
          this%ch4_ebul_depth_unsat_col   (c,:) = spval
          this%ch4_ebul_total_sat_col     (c)   = spval
          this%ch4_ebul_total_unsat_col   (c)   = spval
          this%ch4_surf_ebul_sat_col      (c)   = spval
          this%ch4_surf_ebul_unsat_col    (c)   = spval
          this%ch4_surf_ebul_lake_col     (c)   = spval
          this%ch4_surf_diff_sat_col      (c)   = spval
          this%ch4_surf_diff_unsat_col    (c)   = spval
          this%ch4_surf_diff_lake_col     (c)   = spval
          this%ch4_dfsat_flux_col         (c)   = spval
          this%zwt_ch4_unsat_col          (c)   = spval
          this%conc_ch4_lake_col          (c,:) = spval
          this%conc_o2_lake_col           (c,:) = spval
          this%sif_col                    (c)   = spval
          this%o2stress_unsat_col         (c,:) = spval
          this%o2stress_sat_col           (c,:) = spval
          this%ch4stress_unsat_col        (c,:) = spval
          this%ch4stress_sat_col          (c,:) = spval
          this%finundated_col             (c)   = spval
          this%grnd_ch4_cond_col          (c)   = spval

          ! totcolch4 Set to zero for inactive columns so that this can be used
          ! as an appropriate area-weighted gridcell average soil methane content.
          this%totcolch4_col              (c)   = 0._r8

       end if
    end do

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine Restart( this, bounds, ncid, flag )
    !
    ! !DESCRIPTION:
    ! Read/Write biogeophysics information to/from restart file.
    !
    ! !USES:
    use ncdio_pio , only : ncd_double
    use pio       , only : file_desc_t
    use decompMod , only : bounds_type
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(ch4_type) :: this
    type(bounds_type), intent(in)    :: bounds
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*),  intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    logical :: readvar      ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='CONC_O2_SAT', xtype=ncd_double, &
         dim1name='column', dim2name='levgrnd', switchdim=.true., &
         long_name='oxygen soil concentration', units='mol/m^3', &
         readvar=readvar, interpinic_flag='interp', data=this%conc_o2_sat_col)

    call restartvar(ncid=ncid, flag=flag, varname='CONC_O2_UNSAT', xtype=ncd_double, &
         dim1name='column', dim2name='levgrnd', switchdim=.true., &
         long_name='oxygen soil concentration', units='mol/m^3', &
         readvar=readvar, interpinic_flag='interp', data=this%conc_o2_unsat_col)

    call restartvar(ncid=ncid, flag=flag, varname='O2STRESS_SAT', xtype=ncd_double, &
         dim1name='column', dim2name='levgrnd', switchdim=.true., &
         long_name='oxygen stress fraction', units='', &
         readvar=readvar, interpinic_flag='interp', data=this%o2stress_sat_col)

    call restartvar(ncid=ncid, flag=flag, varname='O2STRESS_UNSAT', xtype=ncd_double, &
         dim1name='column', dim2name='levgrnd', switchdim=.true., &
         long_name='oxygen stress fraction', units='', &
         readvar=readvar, interpinic_flag='interp', data=this%o2stress_unsat_col)

    call restartvar(ncid=ncid, flag=flag, varname='O2_DECOMP_DEPTH_SAT', xtype=ncd_double, &
         dim1name='column', dim2name='levgrnd', switchdim=.true., &
         long_name='O2 consumption during decomposition', units='mol/m3/s', &
         readvar=readvar, interpinic_flag='interp', data=this%o2_decomp_depth_sat_col)

    call restartvar(ncid=ncid, flag=flag, varname='O2_DECOMP_DEPTH_UNSAT', xtype=ncd_double, &
         dim1name='column', dim2name='levgrnd', switchdim=.true., &
         long_name='O2 consumption during decomposition', units='mol/m3/s', &
         readvar=readvar, interpinic_flag='interp', data=this%o2_decomp_depth_unsat_col)

    call restartvar(ncid=ncid, flag=flag, varname='CONC_CH4_SAT', xtype=ncd_double, &
         dim1name='column', dim2name='levgrnd', switchdim=.true., &
         long_name='methane soil concentration', units='mol/m^3', &
         readvar=readvar, interpinic_flag='interp', data=this%conc_ch4_sat_col)

    call restartvar(ncid=ncid, flag=flag, varname='CONC_CH4_UNSAT', xtype=ncd_double, &
         dim1name='column', dim2name='levgrnd', switchdim=.true., &
         long_name='methane soil concentration', units='mol/m^3', &
         readvar=readvar, interpinic_flag='interp', data=this%conc_ch4_unsat_col)

    call restartvar(ncid=ncid, flag=flag, varname='LAYER_SAT_LAG', xtype=ncd_double, &
         dim1name='column', dim2name='levgrnd', switchdim=.true., &
         long_name='lagged saturation status of layer in unsat. zone', units='', &
         readvar=readvar, interpinic_flag='interp', data=this%layer_sat_lag_col)

    call restartvar(ncid=ncid, flag=flag, varname='QFLX_SURF_LAG', xtype=ncd_double, &
         dim1name='column', &
         long_name='time-lagged surface runoff', units='mm/s', &
         readvar=readvar, interpinic_flag='interp', data=this%qflx_surf_lag_col)

    call restartvar(ncid=ncid, flag=flag, varname='FINUNDATED_LAG', xtype=ncd_double, &
         dim1name='column', &
         long_name='time-lagged inundated fraction', units='', &
         readvar=readvar, interpinic_flag='interp', data=this%finundated_lag_col)

    call restartvar(ncid=ncid, flag=flag, varname='FINUNDATED', xtype=ncd_double, &
            dim1name='column', &
            long_name='inundated fraction', units='', &
            readvar=readvar, interpinic_flag='interp', data=this%fsat_bef_col)

    call restartvar(ncid=ncid, flag=flag, varname='annavg_somhr', xtype=ncd_double,  &
         dim1name='column',&
         long_name='Annual Average SOMHR',units='gC/m^2/s', &
         readvar=readvar, interpinic_flag='interp', data=this%annavg_somhr_col)

    call restartvar(ncid=ncid, flag=flag, varname='annavg_finrw', xtype=ncd_double,  &
         dim1name='column',&
         long_name='Annual Average Respiration-Weighted FINUNDATED',units='', &
         readvar=readvar, interpinic_flag='interp', data=this%annavg_finrw_col)

    call restartvar(ncid=ncid, flag=flag, varname='annsum_counter_ch4', xtype=ncd_double,  &
         dim1name='column',&
         long_name='CH4 Ann. Sum Time Counter',units='s', &
         readvar=readvar, interpinic_flag='interp', data=this%annsum_counter_col)

    call restartvar(ncid=ncid, flag=flag, varname='tempavg_somhr', xtype=ncd_double,  &
         dim1name='column',&
         long_name='Temp. Average SOMHR',units='gC/m^2/s', &
         readvar=readvar, interpinic_flag='interp', data=this%tempavg_somhr_col)

    call restartvar(ncid=ncid, flag=flag, varname='tempavg_finrw', xtype=ncd_double,  &
         dim1name='column',&
         long_name='Temp. Average Respiration-Weighted FINUNDATED',units='', &
         readvar=readvar, interpinic_flag='interp', data=this%tempavg_finrw_col)

    call restartvar(ncid=ncid, flag=flag, varname='LAKE_SOILC', xtype=ncd_double, &
         dim1name='column', dim2name='levgrnd', switchdim=.true.,&
         long_name='lake soil carbon concentration', units='g/m^3', &
         readvar=readvar, interpinic_flag='interp', data=this%lake_soilc_col)

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine readCH4Params ( ncid )
    !
    ! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use shr_infnan_mod, only: spval => shr_infnan_nan, assignment(=)
    use ncdio_pio    , only : file_desc_t,ncd_io
    use CH4varcon    , only : use_aereoxid_prog
    !
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'CH4ParamsType'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: tempr ! temporary to read in constant
    character(len=100) :: tString ! temp. var for reading
    !--------------------------------------------------------------------

    allocate(CH4ParamsInst%q10ch4           )
    allocate(CH4ParamsInst%q10ch4base       )
    allocate(CH4ParamsInst%f_ch4            )
    allocate(CH4ParamsInst%rootlitfrac      )
    allocate(CH4ParamsInst%cnscalefactor    )
    allocate(CH4ParamsInst%redoxlag         )
    allocate(CH4ParamsInst%lake_decomp_fact )
    allocate(CH4ParamsInst%redoxlag_vertical)
    allocate(CH4ParamsInst%pHmax            )
    allocate(CH4ParamsInst%pHmin            )
    allocate(CH4ParamsInst%oxinhib          )
    allocate(CH4ParamsInst%vmax_ch4_oxid    )
    allocate(CH4ParamsInst%k_m              )
    allocate(CH4ParamsInst%q10_ch4oxid      )
    allocate(CH4ParamsInst%smp_crit         )
    allocate(CH4ParamsInst%k_m_o2           )
    allocate(CH4ParamsInst%k_m_unsat        )
    allocate(CH4ParamsInst%vmax_oxid_unsat  )
    allocate(CH4ParamsInst%aereoxid         )
    allocate(CH4ParamsInst%scale_factor_aere )
    allocate(CH4ParamsInst%nongrassporosratio)
    allocate(CH4ParamsInst%unsat_aere_ratio  )
    allocate(CH4ParamsInst%porosmin          )
    allocate(CH4ParamsInst%vgc_max           )
    allocate(CH4ParamsInst%satpow              )
    allocate(CH4ParamsInst%scale_factor_gasdiff)
    allocate(CH4ParamsInst%scale_factor_liqdiff)
    allocate(CH4ParamsInst%capthick            )
    allocate(CH4ParamsInst%f_sat               )
    allocate(CH4ParamsInst%qflxlagd            )
    allocate(CH4ParamsInst%highlatfact         )
    allocate(CH4ParamsInst%q10lakebase         )
    allocate(CH4ParamsInst%atmch4              )
    allocate(CH4ParamsInst%rob                 )


    if ( .not. use_aereoxid_prog ) then
        tString='aereoxid'
        call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
        CH4ParamsInst%aereoxid=tempr
     else
        ! value should never be used.
        CH4ParamsInst%aereoxid=nan
     endif

     tString='q10ch4'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%q10ch4=tempr

     tString='q10ch4base'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%q10ch4base=tempr

     tString='f_ch4'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%f_ch4=tempr

     tString='rootlitfrac'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%rootlitfrac=tempr

     tString='cnscalefactor'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%cnscalefactor=tempr

     tString='redoxlag'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%redoxlag=tempr

     tString='lake_decomp_fact'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%lake_decomp_fact=tempr

     tString='redoxlag_vertical'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%redoxlag_vertical=tempr

     tString='pHmax'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%pHmax=tempr

     tString='pHmin'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%pHmin=tempr

     tString='vmax_ch4_oxid'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%vmax_ch4_oxid=45.e-6_r8 * 1000._r8 / 3600._r8
     ! FIX(FIX(SPM,032414),032414) can't be read off of param file.  not bfb since it is a divide
     !CH4ParamsInst%vmax_ch4_oxid=tempr

     tString='oxinhib'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%oxinhib=tempr

     tString='k_m'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%k_m= 5.e-6_r8 * 1000._r8
     ! FIX(FIX(SPM,032414),032414) can't be read off of param file.  not bfb since it is a divide
     !CH4ParamsInst%k_m=tempr

     tString='q10_ch4oxid'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%q10_ch4oxid=tempr

     tString='smp_crit'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%smp_crit=tempr

     tString='k_m_o2'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%k_m_o2  = 20.e-6_r8 * 1000._r8
     ! FIX(FIX(SPM,032414),032414) can't be read off of param file.  not bfb since it is a divide
     !CH4ParamsInst%k_m_o2=tempr

     tString='k_m_unsat'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%k_m_unsat= 5.e-6_r8 * 1000._r8 / 10._r8
     ! FIX(FIX(SPM,032414),032414) can't be read off of param file.  not bfb since it is a divide
     !CH4ParamsInst%k_m_unsat=tempr

     tString='vmax_oxid_unsat'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%vmax_oxid_unsat = 45.e-6_r8 * 1000._r8 / 3600._r8 / 10._r8
     ! FIX(FIX(SPM,032414),032414) can't be read off of param file.  not bfb since it is a divide
     !CH4ParamsInst%vmax_oxid_unsat=tempr

     tString='scale_factor_aere'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%scale_factor_aere=tempr

     tString='nongrassporosratio'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%nongrassporosratio=tempr

     tString='unsat_aere_ratio'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%unsat_aere_ratio= 0.05_r8 / 0.3_r8
     ! FIX(FIX(SPM,032414),032414) can't be read off of param file.  not bfb since it is a divide
     !CH4ParamsInst%unsat_aere_ratio=tempr

     tString='porosmin'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%porosmin=tempr

     tString='vgc_max'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%vgc_max=tempr

     tString='satpow'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%satpow=tempr

     tString='scale_factor_gasdiff'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%scale_factor_gasdiff=tempr

     tString='scale_factor_liqdiff'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%scale_factor_liqdiff=tempr

     tString='f_sat'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%f_sat=tempr

     tString='qflxlagd'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%qflxlagd=tempr

     tString='highlatfact'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%highlatfact=tempr

     tString='q10lakebase'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%q10lakebase=tempr

     tString='atmch4'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%atmch4=tempr

     tString='rob'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%rob=tempr

     tString='capthick'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CH4ParamsInst%capthick=tempr

  end subroutine readCH4Params

  !-----------------------------------------------------------------------
  subroutine CH4 (bounds, &
       num_soilc, filter_soilc, &
       num_lakec, filter_lakec, &
       num_soilp, filter_soilp, &
       lakestate_vars, canopystate_vars, soilstate_vars, soilhydrology_vars, &
       energyflux_vars, ch4_vars, lnd2atm_vars )
    !
    ! !DESCRIPTION:
    ! Driver for the methane emissions model
    !
    ! !USES:
    use subgridAveMod , only : p2c_1d_filter_parallel, p2c_2d_parallel, c2g_1d_parallel, unity 
    use elm_varpar    , only : nlevgrnd, nlevdecomp
    use pftvarcon     , only : noveg
    use CH4varcon     , only : replenishlakec, allowlakeprod, ch4offline, fin_use_fsat
    use elm_varcon    , only : secspday
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc          ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_soilc(:)    ! column filter for soil points
    integer                  , intent(in)    :: num_lakec          ! number of column lake points in column filter
    integer                  , intent(in)    :: filter_lakec(:)    ! column filter for lake points
    integer                  , intent(in)    :: num_soilp          ! number of soil points in pft filter
    integer                  , intent(in)    :: filter_soilp(:)    ! patch filter for soil points
    type(lakestate_type)     , intent(in)    :: lakestate_vars
    type(canopystate_type)   , intent(in)    :: canopystate_vars
    type(soilstate_type)     , intent(inout) :: soilstate_vars
    type(soilhydrology_type) , intent(in)    :: soilhydrology_vars
    type(energyflux_type)    , intent(inout) :: energyflux_vars
    type(ch4_type)           , intent(inout) :: ch4_vars
    type(lnd2atm_type)       , intent(inout) :: lnd2atm_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: sat                                    ! 0 = unsatured, 1 = saturated
    logical  :: lake                                   ! lake or not lake
    integer  :: j,fc,c,g,fp,p,t                        ! indices
    real(r8) :: dtime_ch4                              ! ch4 model time step (sec)
    integer  :: nstep
    integer  :: jwt(bounds%begc:bounds%endc)           ! index of the soil layer right above the water table (-)
    real(r8) :: ch4_prod_tot(bounds%begc:bounds%endc)  ! CH4 production for column (g C/m**2/s)
    real(r8) :: ch4_oxid_tot(bounds%begc:bounds%endc)  ! CH4 oxidation for column (g C/m**2/s)
    real(r8) :: nem_col(bounds%begc:bounds%endc)       ! net adjustment to atm. C flux from methane production (g C/m**2/s)
    real(r8) :: totalsat
    real(r8) :: totalunsat
    real(r8) :: dfsat
    real(r8) :: rootfraction(1:num_soilp, 1:nlevgrnd)
    real(r8) :: totsoilcolch4_bef(1:num_soilc) ! g C / m^2
    real(r8) :: totlakecolch4_bef(1:num_lakec) ! g C / m^2
    real(r8) :: errch4                                 ! g C / m^2
    real(r8) :: zwt_actual
    real(r8) :: qflxlags                               ! Time to lag qflx_surf_lag (s)
    real(r8) :: redoxlag                               ! Redox time lag
    real(r8) :: redoxlag_vertical                      ! Vertical redox lag time
    real(r8) :: atmch4                                 ! Atmospheric CH4 mixing ratio to
                                                       ! prescribe if not provided by the atmospheric model (= 1.7e-6_r8) (mol/mol)
    real(r8) :: redoxlags                              ! Redox time lag in s
    real(r8) :: redoxlags_vertical                     ! Vertical redox lag time in s
    real(r8) :: qflxlagd                               ! days to lag qflx_surf_lag in the tropics (days)
    real(r8) :: highlatfact                            ! multiple of qflxlagd for high latitudes
    integer  :: dummyfilter(1)                         ! empty filter
    integer  :: erridx
    real(r8) :: sum1,sum2,sum3
    integer  :: begg,endg,begc,endc,begp,endp
    real :: startt, stopt
    !-----------------------------------------------------------------------

    associate(                                                                 &
         dz                   =>   col_pp%dz                                 , & ! Input:  [real(r8) (:,:) ]  layer thickness (m)  (-nlevsno+1:nlevsoi)
         zi                   =>   col_pp%zi                                 , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m)
         z                    =>   col_pp%z                                  , & ! Input:  [real(r8) (:,:) ]  layer depth (m) (-nlevsno+1:nlevsoi)

         forc_t               =>   top_as%tbot                               , & ! Input:  [real(r8) (:)   ]  atmospheric temperature (Kelvin)
         forc_pbot            =>   top_as%pbot                               , & ! Input:  [real(r8) (:)   ]  atmospheric pressure (Pa)
         forc_po2             =>   top_as%po2bot                             , & ! Input:  [real(r8) (:)   ]  O2 partial pressure (Pa)
         forc_pco2            =>   top_as%pco2bot                            , & ! Input:  [real(r8) (:)   ]  CO2 partial pressure (Pa)
         forc_pch4            =>   top_as%pch4bot                            , & ! Input:  [real(r8) (:)   ]  CH4 partial pressure (Pa)

         zwt                  =>   soilhydrology_vars%zwt_col                , & ! Input:  [real(r8) (:)   ]  water table depth (m)
         zwt_perched          =>   soilhydrology_vars%zwt_perched_col        , & ! Input:  [real(r8) (:)   ]  perched water table depth (m)

         rootfr               =>   soilstate_vars%rootfr_patch               , & ! Input:  [real(r8) (:,:) ]  fraction of roots in each soil layer  (nlevgrnd)
         rootfr_col           =>   soilstate_vars%rootfr_col                 , & ! Output: [real(r8) (:,:) ]  fraction of roots in each soil layer  (nlevgrnd) (p2c)

         frac_h2osfc          =>   col_ws%frac_h2osfc            , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by surface water (0 to 1)
         qflx_surf            =>   col_wf%qflx_surf              , & ! Input:  [real(r8) (:)   ]  surface runoff (mm H2O /s)

         conc_o2_sat          =>   ch4_vars%conc_o2_sat_col                  , & ! Input:  [real(r8) (:,:) ]  O2 conc  in each soil layer (mol/m3) (nlevsoi)
         zwt0                 =>   ch4_vars%zwt0_col                         , & ! Input:  [real(r8) (:)   ]  decay factor for finundated (m)
         f0                   =>   ch4_vars%f0_col                           , & ! Input:  [real(r8) (:)   ]  maximum gridcell fractional inundated area
         p3                   =>   ch4_vars%p3_col                           , & ! Input:  [real(r8) (:)   ]  coefficient for qflx_surf_lag for finunated (s/mm)

         grnd_ch4_cond_patch  =>   ch4_vars%grnd_ch4_cond_patch              , & ! Input:  [real(r8) (:)   ]  tracer conductance for boundary layer [m/s]
         grnd_ch4_cond_col    =>   ch4_vars%grnd_ch4_cond_col                , & ! Output: [real(r8) (:)   ]  tracer conductance for boundary layer [m/s] (p2c)

         ch4_surf_diff_sat    =>   ch4_vars%ch4_surf_diff_sat_col            , & ! Output: [real(r8) (:)   ]  CH4 surface flux (mol/m2/s)
         ch4_surf_diff_unsat  =>   ch4_vars%ch4_surf_diff_unsat_col          , & ! Output: [real(r8) (:)   ]  CH4 surface flux (mol/m2/s)
         ch4_surf_diff_lake   =>   ch4_vars%ch4_surf_diff_lake_col           , & ! Output: [real(r8) (:)   ]  CH4 surface flux (mol/m2/s)
         ch4_surf_ebul_sat    =>   ch4_vars%ch4_surf_ebul_sat_col            , & ! Output: [real(r8) (:)   ]  CH4 ebullition to atmosphere (mol/m2/s)
         ch4_surf_ebul_unsat  =>   ch4_vars%ch4_surf_ebul_unsat_col          , & ! Output: [real(r8) (:)   ]  CH4 ebullition to atmosphere (mol/m2/s)
         ch4_surf_ebul_lake   =>   ch4_vars%ch4_surf_ebul_lake_col           , & ! Output: [real(r8) (:)   ]  CH4 ebullition to atmosphere (mol/m2/s)
         ch4_surf_aere_sat    =>   ch4_vars%ch4_surf_aere_sat_col            , & ! Output: [real(r8) (:)   ]  Total column CH4 aerenchyma (mol/m2/s)
         ch4_surf_aere_unsat  =>   ch4_vars%ch4_surf_aere_unsat_col          , & ! Output: [real(r8) (:)   ]  Total column CH4 aerenchyma (mol/m2/s)
         fsat_bef             =>   ch4_vars%fsat_bef_col                     , & ! Output: [real(r8) (:)   ]  finundated from previous timestep
         ch4_oxid_depth_sat   =>   ch4_vars%ch4_oxid_depth_sat_col           , & ! Output: [real(r8) (:,:) ]  CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
         ch4_oxid_depth_unsat =>   ch4_vars%ch4_oxid_depth_unsat_col         , & ! Output: [real(r8) (:,:) ]  CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
         ch4_oxid_depth_lake  =>   ch4_vars%ch4_oxid_depth_lake_col          , & ! Output: [real(r8) (:,:) ]  CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
         ch4_prod_depth_sat   =>   ch4_vars%ch4_prod_depth_sat_col           , & ! Output: [real(r8) (:,:) ]  production of CH4 in each soil layer (nlevsoi) (mol/m3/s)
         ch4_prod_depth_unsat =>   ch4_vars%ch4_prod_depth_unsat_col         , & ! Output: [real(r8) (:,:) ]  production of CH4 in each soil layer (nlevsoi) (mol/m3/s)
         ch4_prod_depth_lake  =>   ch4_vars%ch4_prod_depth_lake_col          , & ! Output: [real(r8) (:,:) ]  production of CH4 in each soil layer (nlevsoi) (mol/m3/s)
         lake_soilc           =>   ch4_vars%lake_soilc_col                   , & ! Output: [real(r8) (:,:) ]  total soil organic matter found in level (g C / m^3) (nlevsoi)
         conc_ch4_sat         =>   ch4_vars%conc_ch4_sat_col                 , & ! Output: [real(r8) (:,:) ]  CH4 conc in each soil layer (mol/m3) (nlevsoi)
         conc_ch4_unsat       =>   ch4_vars%conc_ch4_unsat_col               , & ! Output: [real(r8) (:,:) ]  CH4 conc in each soil layer (mol/m3) (nlevsoi)
         conc_ch4_lake        =>   ch4_vars%conc_ch4_lake_col                , & ! Output: [real(r8) (:,:) ]  CH4 conc in each soil layer (mol/m3) (nlevsoi)
         conc_o2_lake         =>   ch4_vars%conc_o2_lake_col                 , & ! Output: [real(r8) (:,:) ]  O2 conc  in each soil layer (mol/m3) (nlevsoi)
         ch4_dfsat_flux       =>   ch4_vars%ch4_dfsat_flux_col               , & ! Output: [real(r8) (:)   ]  CH4 flux to atm due to decreasing finundated (kg C/m^2/s) [+]
         zwt_ch4_unsat        =>   ch4_vars%zwt_ch4_unsat_col                , & ! Output: [real(r8) (:)   ]  depth of water table for unsaturated fraction (m)
         totcolch4            =>   ch4_vars%totcolch4_col                    , & ! Output: [real(r8) (:)   ]  total methane in soil column (g C / m^2)
         finundated           =>   ch4_vars%finundated_col                   , & ! Output: [real(r8) (:)   ]  fractional inundated area in soil column (excluding dedicated wetland columns)
         qflx_surf_lag        =>   ch4_vars%qflx_surf_lag_col                , & ! Output: [real(r8) (:)   ]  time-lagged surface runoff (mm H2O /s)
         finundated_lag       =>   ch4_vars%finundated_lag_col               , & ! Output: [real(r8) (:)   ]  time-lagged fractional inundated area
         layer_sat_lag        =>   ch4_vars%layer_sat_lag_col                , & ! Output: [real(r8) (:,:) ]  Lagged saturation status of soil layer in the unsaturated zone (1 = sat)
         c_atm                =>   ch4_vars%c_atm_grc                        , & ! Output: [real(r8) (:,:) ]  CH4, O2, CO2 atmospheric conc  (mol/m3)
         ch4co2f              =>   ch4_vars%ch4co2f_grc                      , & ! Output: [real(r8) (:)   ]  gridcell CO2 production from CH4 oxidation (g C/m**2/s)
         ch4prodg             =>   ch4_vars%ch4prodg_grc                     , & ! Output: [real(r8) (:)   ]  gridcell average CH4 production (g C/m^2/s)
         ch4_surf_flux_tot    =>   ch4_vars%ch4_surf_flux_tot_col            , & ! Output: [real(r8) (:)   ]  col CH4 flux to atm. (kg C/m**2/s)
         !Added to pass directly to sat, unsat routines 
         o2_decomp_depth_unsat  => ch4_vars%o2_decomp_depth_unsat_col , &
         o2stress_unsat         => ch4_vars%o2stress_unsat_col        , &
         ch4_aere_depth_unsat   => ch4_vars%ch4_aere_depth_unsat_col  , &
         ch4_tran_depth_unsat   => ch4_vars%ch4_tran_depth_unsat_col  , &
         ch4_ebul_depth_unsat   => ch4_vars%ch4_ebul_depth_unsat_col  , &
         ch4_ebul_total_unsat   => ch4_vars%ch4_ebul_total_unsat_col  , &
         co2_oxid_depth_unsat   => ch4_vars%co2_oxid_depth_unsat_col  , &
         o2_oxid_depth_unsat    => ch4_vars%o2_oxid_depth_unsat_col   , &
         o2_aere_depth_unsat    => ch4_vars%o2_aere_depth_unsat_col   , &
         co2_aere_depth_unsat    => ch4_vars%co2_aere_depth_unsat_col , &
         ch4stress_unsat        => ch4_vars%ch4stress_unsat_col       , &
         co2_decomp_depth_unsat => ch4_vars%co2_decomp_depth_unsat_col, &
         conc_o2_unsat          => ch4_vars%conc_o2_unsat_col         , &
         ! Sat variables 
         o2_decomp_depth_sat  => ch4_vars%o2_decomp_depth_sat_col  , &
         o2stress_sat         => ch4_vars%o2stress_sat_col         , &
         co2_aere_depth_sat   => ch4_vars%co2_aere_depth_sat_col   , &
         ch4_aere_depth_sat   => ch4_vars%ch4_aere_depth_sat_col   , &
         ch4_tran_depth_sat   => ch4_vars%ch4_tran_depth_sat_col   , &
         ch4_ebul_depth_sat   => ch4_vars%ch4_ebul_depth_sat_col   , &
         ch4_ebul_total_sat   => ch4_vars%ch4_ebul_total_sat_col   , &
         co2_oxid_depth_sat   => ch4_vars%co2_oxid_depth_sat_col   , &
         o2_oxid_depth_sat    => ch4_vars%o2_oxid_depth_sat_col    , &
         o2_aere_depth_sat    => ch4_vars%o2_aere_depth_sat_col    , &
         ch4stress_sat        => ch4_vars%ch4stress_sat_col        , &
         co2_decomp_depth_sat => ch4_vars%co2_decomp_depth_sat_col , &
         !
         nem_grc              =>   lnd2atm_vars%nem_grc            , & ! Output: [real(r8) (:)   ]  gridcell average net methane correction to CO2 flux (g C/m^2/s)
         redoxlag             => CH4ParamsInst%redoxlag , &
         redoxlag_vertical    => CH4ParamsInst%redoxlag_vertical , &
         atmch4               => CH4ParamsInst%atmch4      , &
         qflxlagd             => CH4ParamsInst%qflxlagd    , &
         highlatfact          => CH4ParamsInst%highlatfact   &
         )
      
      call cpu_time(startt) 
      !$acc enter data create(&
      !$acc jwt(:), &
      !$acc ch4_prod_tot(:), &
      !$acc ch4_oxid_tot(:), &
      !$acc nem_col(:), &
      !$acc rootfraction(:,:), &
      !$acc totsoilcolch4_bef(:), &
      !$acc totlakecolch4_bef(:), &
      !$acc dummyfilter(:))
      !$acc enter data create(sum1,sum2,sum3)

      begg =   bounds%begg
      endg =   bounds%endg
      begc =   bounds%begc
      endc =   bounds%endc
      begp =   bounds%begp
      endp =   bounds%endp

      dtime_ch4 = dtime_mod
      redoxlags = redoxlag*secspday ! days --> s
      redoxlags_vertical = redoxlag_vertical*secspday ! days --> s
      rgasm = rgas / 1000._r8
      lake = .false.; sat = 0 
      !$acc enter data copyin(redoxlags,redoxlags_vertical, dtime_ch4, lake, sat)
      !$acc update device(rgasm) 

      ! Initialize local fluxes to zero: necessary for columns outside the filters because averaging up to gridcell will be done
      ! ch4_surf_flux_tot(begc:endc) = 0._r8
      ! Adjustment to NEE for methane production - oxidation
      ! nem_col(begc:endc)           = 0._r8

      !$acc parallel loop independent gang vector default(present)
      do c = begc , endc 
         ch4_surf_flux_tot(c) = 0._r8 
         nem_col(c) = 0._r8
      end do 

      !$acc parallel loop independent gang vector default(present)
      do g= begg, endg
         ! using the first topounit on this gridcell for access to atmospheric forcing
         ! PET 8/10/2018 - replace this later once gas concentrations (c_atm) are also being tracked at
         ! topounit level
         t = grc_pp%topi(g)
         if (ch4offline) then
            forc_pch4(t) = atmch4*forc_pbot(t)
         else
            if (forc_pch4(t) == 0._r8) then
               !write(iulog,*) 'not using ch4offline, but methane concentration not passed from the atmosphere'//&
               !      'to land model! CLM Model is stopping.'
               !call endrun(msg=' ERROR: Methane not being passed to atmosphere'//errMsg(__FILE__, __LINE__))
               print *, ' ERROR: Methane not being passed to atmosphere'
            end if
         end if
         c_atm(g,1) =  forc_pch4(t) / rgasm / forc_t(t) ! [mol/m3 air]
         c_atm(g,2) =  forc_po2(t)  / rgasm / forc_t(t) ! [mol/m3 air]
         c_atm(g,3) =  forc_pco2(t) / rgasm / forc_t(t) ! [mol/m3 air]
      end do
      
      ! Initialize CH4 balance and calculate finundated
      !$acc parallel loop independent gang vector default(present)
      do fc = 1, num_soilc
         c = filter_soilc(fc)
         g = col_pp%gridcell(c)

         totsoilcolch4_bef(fc) = totcolch4(c)
         totcolch4(c) = 0._r8

         ! Update lagged surface runoff

         if (grc_pp%latdeg(g) < 45._r8) then
            qflxlags = qflxlagd * secspday ! 30 days
         else
            qflxlags = qflxlagd * secspday * highlatfact ! 60 days
         end if
         qflx_surf_lag(c) = qflx_surf_lag(c) * exp(-dtime_mod/qflxlags) &
              + qflx_surf(c) * (1._r8 - exp(-dtime_mod/qflxlags))

         !There may be ways to improve this for irrigated crop columns...
         if (fin_use_fsat) then
            finundated(c) = frac_h2osfc(c)
         else
            if (zwt0(c) > 0._r8) then
               if (zwt_perched(c) < z(c,nlevsoi)-1.e-5_r8 .and. zwt_perched(c) < zwt(c)) then
                  zwt_actual = zwt_perched(c)
               else
                  zwt_actual = zwt(c)
               end if
               finundated(c) = f0(c) * exp(-zwt_actual/zwt0(c)) + p3(c)*qflx_surf_lag(c)
            else
               finundated(c) = p3(c)*qflx_surf_lag(c)
            end if
         end if
         finundated(c) = max( min(finundated(c),1._r8), 0._r8)

         ! Update lagged finundated for redox calculation
         if (redoxlags > 0._r8) then
            finundated_lag(c) = finundated_lag(c) * exp(-dtime_mod/redoxlags) &
                 + finundated(c) * (1._r8 - exp(-dtime_mod/redoxlags))
         else
            finundated_lag(c) = finundated(c)
         end if

      end do

      !$acc parallel loop independent gang vector default(present)
      do fc = 1, num_lakec
         c = filter_lakec(fc)

         totlakecolch4_bef(fc) = totcolch4(c)
         totcolch4(c) = 0._r8
      end do

      ! Check to see if finundated changed since the last timestep.  If it increased, then reduce conc_ch4_sat
      ! proportionally.  If it decreased, then add flux to atm.
      !$acc parallel loop independent gang worker default(present) private(sum1)
      do fc = 1, num_soilc
         c = filter_soilc(fc)
         sum1 = 0._r8
         if (fsat_bef(c) /= spval .and. finundated(c) > fsat_bef(c)) then !Reduce conc_ch4_sat
            !$acc loop independent vector 
            do j=1,nlevsoi
               dfsat = finundated(c) - fsat_bef(c)
               conc_ch4_sat(c,j) = (fsat_bef(c)*conc_ch4_sat(c,j) + dfsat*conc_ch4_unsat(c,j)) / finundated(c)
            end do 

         else if (fsat_bef(c) /= spval .and. finundated(c) < fsat_bef(c)) then
            !$acc loop vector reduction(+:sum1)
            do j=1,nlevsoi 
               sum1 = sum1 + (fsat_bef(c) - finundated(c))*(conc_ch4_sat(c,j) - conc_ch4_unsat(c,j)) * &
                    dz(c,j) / dtime_mod * catomw / 1000._r8 ! mol --> kg
            end do
         end if
         ch4_dfsat_flux(c) = sum1 
      end do

      call cpu_time(stopt) 
      print *, "TIMING Init : ",(stopt-startt)*1.E+3,"ms"
      !!!! Begin biochemistry
      
      ! First for soil
      lake = .false.
      ! Do CH4 Annual Averages
      call ch4_annualupdate(bounds, &
           num_soilc, filter_soilc, &
           num_soilp, filter_soilp, ch4_vars, dtime_mod, 365.d0)

      ! Determine rootfr_col and also check for inactive columns

      if (nlevdecomp == 1) then

         ! Set rootfraction to spval for non-veg points, unless veg_pp%wtcol > 0.99,
         ! in which case set it equal to uniform dist.
         !$acc parallel loop independent gang vector default(present) collapse(2) 
         do j=1, nlevsoi
            do fp = 1, num_soilp
               p = filter_soilp(fp)
               c = veg_pp%column(p)

               if (veg_pp%itype(p) /= noveg) then
                  rootfraction(fp,j) = rootfr(p,j)
               else if (veg_pp%wtcol(p) < 0.99_r8) then
                  rootfraction(fp,j) = spval
               else
                  rootfraction(fp,j) = dz(c,j) / zi(c,nlevsoi)   ! Set equal to uniform distribution
               end if
            end do
         end do

         call p2c_2d_parallel(bounds, nlevgrnd, &
              rootfraction(bounds%begp:bounds%endp, :), &
              rootfr_col(bounds%begc:bounds%endc, :))

         !$acc parallel loop independent gang vector default(present) collapse(2) 
         do j=1, nlevsoi
            do fc = 1, num_soilc
               c = filter_soilc(fc)
               if (.not. col_pp%active(c)) rootfr_col(c,j) = dz(c,j) / zi(c,nlevsoi)
            end do
         end do
      end if

      ! Determine grnd_ch4_cond_col
      ! Needed to use non-filter form above so that spval would be treated properly.

      call p2c_1d_filter_parallel(num_soilc, filter_soilc, &
           grnd_ch4_cond_patch(bounds%begp:bounds%endp), &
           grnd_ch4_cond_col(bounds%begc:bounds%endc))

      ! Set the gridcell atmospheric CH4 and O2 concentrations
      !$acc parallel loop independent gang vector default(present)
      do fc = 1, num_soilc
         c = filter_soilc(fc)
         g = col_pp%gridcell(c)

         ! using the first topounit on this gridcell for access to atmospheric forcing
         ! PET 8/10/2018 - replace this later once gas concentrations (c_atm) are also being tracked at
         ! topounit level
         t = grc_pp%topi(g)
         c_atm(g,1) =  forc_pch4(t) / rgasm / forc_t(t) ! [mol/m3 air]
         c_atm(g,2) =  forc_po2(t)  / rgasm / forc_t(t) ! [mol/m3 air]
        !c_atm(g,3) =  forc_pco2(t) / rgasm / forc_t(t) ! [mol/m3 air] - Not currently used
      enddo

      !-------------------------------------------------
      ! Loop over saturated and unsaturated, non-lakes
      !------------------------------------------------

      ! do sat = 0, 1 ! 0 == unsaturated; 1 = saturated

      ! Get index of water table
      ! if (sat == 0) then ! unsaturated
      sat = 0 
      call get_jwt (bounds, num_soilc, filter_soilc, jwt(begc:endc), &
           soilstate_vars )

      !$acc parallel loop independent gang vector default(present)
      do fc = 1, num_soilc
         c = filter_soilc(fc)
         zwt_ch4_unsat(c) = zi(c,jwt(c))
      end do

      ! Update lagged saturation status of layer
      !$acc parallel loop independent gang vector default(present) collapse(2) 
      do j=1,nlevsoi
         do fc = 1, num_soilc
            c = filter_soilc(fc)

            if (j > jwt(c) .and. redoxlags_vertical > 0._r8) then ! saturated currently
               layer_sat_lag(c,j) = layer_sat_lag(c,j) * exp(-dtime_mod/redoxlags_vertical) &
                    + (1._r8 - exp(-dtime_mod/redoxlags_vertical))
            else if (redoxlags_vertical > 0._r8) then
               layer_sat_lag(c,j) = layer_sat_lag(c,j) * exp(-dtime_mod/redoxlags_vertical)
            else if (j > jwt(c)) then  ! redoxlags_vertical = 0
               layer_sat_lag(c,j) = 1._r8
            else
               layer_sat_lag(c,j) = 0._r8
            end if
         end do
      end do



      ! calculate CH4 production in each soil layer
      call ch4_prod (bounds, &
           num_soilc, filter_soilc, &
           num_soilp, filter_soilp, &
           jwt(begc:endc), sat, lake, &
           soilstate_vars, ch4_vars, dtime_mod, &
           conc_o2_unsat(begc:endc,:), ch4_prod_depth_unsat(begc:endc,:), &
           o2_decomp_depth_unsat(begc:endc,:), co2_decomp_depth_unsat(begc:endc,:) &
           )
      

      ! calculate CH4 oxidation in each soil layer
      call ch4_oxid (bounds, &
           num_soilc, filter_soilc, &
           jwt(begc:endc), sat, lake, &
           soilstate_vars, dtime_mod,&
           ch4_oxid_depth_unsat(begc:endc,:), &
           o2_oxid_depth_unsat(begc:endc,:),&
           co2_oxid_depth_unsat(begc:endc,:), &   
           conc_ch4_unsat(begc:endc,:), conc_o2_unsat(begc:endc,:), &
           o2_decomp_depth_unsat(begc:endc,:)  )

      

      ! calculate CH4 aerenchyma losses in each soil layer
      call ch4_aere (bounds, &
           num_soilc, filter_soilc, &
           num_soilp, filter_soilp, &
           jwt(begc:endc), sat, lake, &
           canopystate_vars, soilstate_vars,&
           energyflux_vars, ch4_vars,dtime_mod, &
           ch4_aere_depth_unsat(begc:endc,:), ch4_tran_depth_unsat(begc:endc,:), &
           o2_aere_depth_unsat(begc:endc,:), &
           co2_aere_depth_unsat(begc:endc,:), conc_ch4_unsat(begc:endc,:), &
           conc_o2_unsat(begc:endc,:), ch4_oxid_depth_unsat(begc:endc,:),&
           ch4_prod_depth_unsat(begc:endc,:))

      ! calculate CH4 ebullition losses in each soil layer
      call ch4_ebul (bounds, &
           num_soilc, filter_soilc, &
           jwt(begc:endc), sat, lake, &
           lakestate_vars, soilstate_vars, ch4_vars, dtime_mod, &
           ch4_ebul_depth_unsat(begc:endc,:), ch4_ebul_total_unsat(begc:endc),&
           conc_ch4_unsat(begc:endc,:), &    
           ch4_aere_depth_unsat(begc:endc,:), ch4_oxid_depth_unsat(begc:endc,:))


      ! Solve CH4 reaction/diffusion equation
      ! Competition for oxygen will occur here.
      call ch4_tran (bounds, &
           num_soilc, filter_soilc, &
           jwt(begc:endc), dtime_mod, sat, lake, &
           soilstate_vars, energyflux_vars, ch4_vars, dtime_mod,&
           o2_decomp_depth_unsat(begc:endc,:), o2stress_unsat(begc:endc,:), &
           ch4_oxid_depth_unsat(begc:endc,:) ,& 
           ch4_prod_depth_unsat(begc:endc,:), ch4_aere_depth_unsat(begc:endc,:), &
           ch4_surf_aere_unsat(begc:endc) ,&
           ch4_ebul_depth_unsat(begc:endc,:) ,&
           ch4_ebul_total_unsat(begc:endc) ,&
           ch4_surf_ebul_unsat(begc:endc) ,&   
           ch4_surf_diff_unsat(begc:endc), o2_oxid_depth_unsat(begc:endc,:), &
           o2_aere_depth_unsat(begc:endc,:), & 
           ch4stress_unsat(begc:endc,:), co2_decomp_depth_unsat(begc:endc,:), &
           conc_ch4_unsat(begc:endc,:), conc_o2_unsat(begc:endc,:) &
           )

      ! Saturated
      sat = 1 
      !$acc update device(sat)
      !$acc parallel loop independent gang vector default(present)
      do fc = 1, num_soilc
         c = filter_soilc(fc)
         jwt(c) = 0
      end do

      ! calculate CH4 production in each soil layer
      call ch4_prod (bounds, &
           num_soilc, filter_soilc, &
           num_soilp, filter_soilp, &
           jwt(begc:endc), sat, lake, &
           soilstate_vars, ch4_vars, dtime_mod, &
           conc_o2_sat(begc:endc,:), &
           ch4_prod_depth_sat(begc:endc,:),&
           o2_decomp_depth_sat(begc:endc,:), &
           co2_decomp_depth_sat(begc:endc,:) &
           )

      ! calculate CH4 oxidation in each soil layer
      call ch4_oxid (bounds, &
           num_soilc, filter_soilc, &
           jwt(begc:endc), sat, lake, &
           soilstate_vars, dtime_mod,&
           ch4_oxid_depth_sat(begc:endc,:), o2_oxid_depth_sat(begc:endc,:),&
           co2_oxid_depth_sat(begc:endc,:), &   
           conc_ch4_sat(begc:endc,:), conc_o2_sat(begc:endc,:), &
           o2_decomp_depth_sat(begc:endc,:) )

      ! calculate CH4 aerenchyma losses in each soil layer
      call ch4_aere (bounds, &
           num_soilc, filter_soilc, &
           num_soilp, filter_soilp, &
           jwt(begc:endc), sat, lake, &
           canopystate_vars, soilstate_vars, energyflux_vars, ch4_vars,dtime_mod, &
           ch4_aere_depth_sat(begc:endc,:), &
           ch4_tran_depth_sat(begc:endc,:),&
           o2_aere_depth_sat(begc:endc,:), &
           co2_aere_depth_sat(begc:endc,:), & 
           conc_ch4_sat(begc:endc,:), &
           conc_o2_sat(begc:endc,:), &
           ch4_oxid_depth_sat(begc:endc,:), &
           ch4_prod_depth_sat(begc:endc,:))
      
      ! calculate CH4 ebullition losses in each soil layer
      call ch4_ebul (bounds, &
           num_soilc, filter_soilc, &
           jwt(begc:endc), sat, lake, &
           lakestate_vars, soilstate_vars, ch4_vars, dtime_mod, &
           ch4_ebul_depth_sat(begc:endc,:), &
           ch4_ebul_total_sat(begc:endc), &
           conc_ch4_sat(begc:endc,:), &    
           ch4_aere_depth_sat(begc:endc,:), &
           ch4_oxid_depth_sat(begc:endc,:))
      
      ! Solve CH4 reaction/diffusion equation
      ! Competition for oxygen will occur here.
      call ch4_tran (bounds, &
           num_soilc, filter_soilc, &
           jwt(begc:endc), dtime_mod, sat, lake, &
           soilstate_vars, energyflux_vars, ch4_vars, dtime_mod, &
           o2_decomp_depth_sat(begc:endc,:), &
           o2stress_sat(begc:endc,:), &
           ch4_oxid_depth_sat(begc:endc,:),& 
           ch4_prod_depth_sat(begc:endc,:), &
           ch4_aere_depth_sat(begc:endc,:), &
           ch4_surf_aere_sat(begc:endc), &
           ch4_ebul_depth_sat(begc:endc,:),&
           ch4_ebul_total_sat(begc:endc),&
           ch4_surf_ebul_sat(begc:endc), &   
           ch4_surf_diff_sat(begc:endc), &
           o2_oxid_depth_sat(begc:endc,:), &
           o2_aere_depth_sat(begc:endc,:), & 
           ch4stress_sat(begc:endc,:), &
           co2_decomp_depth_sat(begc:endc,:), &
           conc_ch4_sat(begc:endc,:), &
           conc_o2_sat(begc:endc,:) &
           )

      !-------------------------------------------------
      ! Now do over lakes
      !-------------------------------------------------

      if (allowlakeprod) then
         lake = .true.
         sat = 1
         !$acc update device(lake,sat)
         !$acc parallel loop independent gang vector default(present)
         do fc = 1, num_lakec
            c = filter_lakec(fc)
            jwt(c) = 0
         end do

         ! calculate CH4 production in each lake layer
         call ch4_prod (bounds, &
              num_lakec, filter_lakec, &
              0, dummyfilter, jwt(begc:endc), sat, lake, &
              soilstate_vars, ch4_vars, dtime_mod, &
              conc_o2_sat(begc:endc,:), ch4_prod_depth_sat(begc:endc,:), &
              o2_decomp_depth_sat(begc:endc,:), co2_decomp_depth_sat(begc:endc,:) &
              )

         ! calculate CH4 oxidation in each lake layer
         call ch4_oxid (bounds, &
              num_lakec, filter_lakec, &
              jwt(begc:endc), sat, lake, soilstate_vars, dtime_mod,&
              ch4_oxid_depth_sat(begc:endc,:), o2_oxid_depth_sat(begc:endc,:), &
              co2_oxid_depth_sat(begc:endc,:), &   
              conc_ch4_sat(begc:endc,:), &
              conc_o2_sat(begc:endc,:), &
              o2_decomp_depth_sat(begc:endc,:) )

         ! calculate CH4 aerenchyma losses in each lake layer
         ! The p filter will not be used here; the relevant column vars will just be set to 0.
         call ch4_aere (bounds, &
              num_lakec, filter_lakec, &
              0, dummyfilter, jwt(begc:endc), sat, lake, &
              canopystate_vars, soilstate_vars, energyflux_vars, ch4_vars, dtime_mod, &
              ch4_aere_depth_sat(begc:endc,:), &
              ch4_tran_depth_sat(begc:endc,:), &
              o2_aere_depth_sat(begc:endc,:), &
              co2_aere_depth_sat(begc:endc,:), &
              conc_ch4_sat(begc:endc,:), &
              conc_o2_sat(begc:endc,:), &
              ch4_oxid_depth_sat(begc:endc,:), &
              ch4_prod_depth_sat(begc:endc,:))

         ! calculate CH4 ebullition losses in each lake layer
         call ch4_ebul (bounds, &
              num_lakec, filter_lakec, &
              jwt(begc:endc), sat, lake, &
              lakestate_vars, soilstate_vars,ch4_vars, dtime_mod, &
              ch4_ebul_depth_sat(begc:endc,:), &
              ch4_ebul_total_sat(begc:endc), &
              conc_ch4_sat(begc:endc,:), &    
              ch4_aere_depth_sat(begc:endc,:), &
              ch4_oxid_depth_sat(begc:endc,:))

         ! Solve CH4 reaction/diffusion equation
         ! Competition for oxygen will occur here.
         call ch4_tran (bounds, &
              num_lakec, filter_lakec, &
              jwt(begc:endc), dtime_mod, sat, lake, &
              soilstate_vars, energyflux_vars, ch4_vars, dtime_mod, &
              o2_decomp_depth_sat(begc:endc,:), &
              o2stress_sat(begc:endc,:), &
              ch4_oxid_depth_sat(begc:endc,:),& 
              ch4_prod_depth_sat(begc:endc,:), &
              ch4_aere_depth_sat(begc:endc,:), ch4_surf_aere_sat(begc:endc), &
              ch4_ebul_depth_sat(begc:endc,:),ch4_ebul_total_sat(begc:endc), &
              ch4_surf_ebul_sat(begc:endc), &   
              ch4_surf_diff_sat(begc:endc), &
              o2_oxid_depth_sat(begc:endc,:), &
              o2_aere_depth_sat(begc:endc,:), & 
              ch4stress_sat(begc:endc,:), &
              co2_decomp_depth_sat(begc:endc,:), &
              conc_ch4_sat(begc:endc,:),&
              conc_o2_sat(begc:endc,:) &
              )

      end if

      !-------------------------------------------------
      ! Average up to gridcell flux and column oxidation and production rate.
      !-------------------------------------------------

      ! First weight the soil columns by finundated.
      !$acc parallel loop independent gang worker default(present) private(totalsat, totalunsat,sum1,sum2)
      do fc = 1, num_soilc
         c = filter_soilc(fc)
         
         totalsat = ch4_surf_diff_sat(c) + ch4_surf_aere_sat(c) + ch4_surf_ebul_sat(c)
         totalunsat = ch4_surf_diff_unsat(c) + ch4_surf_aere_unsat(c) + ch4_surf_ebul_unsat(c)
         ch4_surf_flux_tot(c) = (finundated(c)*totalsat + (1._r8 - finundated(c))*totalunsat) * catomw / 1000._r8

         !Convert from mol to kg C
         ! ch4_oxid_tot and ch4_prod_tot are initialized to zero above 
         sum1 = 0._r8 
         sum2 = 0._r8
         !$acc loop vector reduction(+:sum1,sum2)
         do j=1,nlevsoi

            sum1 = sum1 + (finundated(c)*ch4_oxid_depth_sat(c,j) + &
                 (1._r8 - finundated(c))*ch4_oxid_depth_unsat(c,j))*dz(c,j) * catomw
            !Convert from mol to g C
            sum2 = sum2 + (finundated(c)*ch4_prod_depth_sat(c,j) + &
                 (1._r8 - finundated(c))*ch4_prod_depth_unsat(c,j))*dz(c,j) * catomw

         end do
         ch4_oxid_tot(c) = sum1
         ch4_prod_tot(c) = sum2  
         ! Adjustment to NEE flux to atm. for methane production
         nem_col(c) = nem_col(c) - ch4_prod_tot(c)
         ! Adjustment to NEE flux to atm. for methane oxidation
         nem_col(c) = nem_col(c) + ch4_oxid_tot(c)
      end do

      ! Correct for discrepancies in CH4 concentration from changing finundated

      !$acc parallel loop independent gang vector default(present)
      do fc = 1, num_soilc
         c = filter_soilc(fc)

         if (fsat_bef(c) /= spval) then ! not first timestep
            ch4_surf_flux_tot(c) = ch4_surf_flux_tot(c) + ch4_dfsat_flux(c)
         end if
         fsat_bef(c) = finundated(c)
      end do

      if (allowlakeprod) then
         !$acc parallel loop independent gang worker default(present) private(sum1,sum2)
         do fc = 1, num_lakec
            c = filter_lakec(fc)
            ! ch4_oxid_tot and ch4_prod_tot are initialized to zero above
            totalsat = ch4_surf_diff_sat(c) + ch4_surf_aere_sat(c) + ch4_surf_ebul_sat(c)
            ch4_surf_flux_tot(c) = totalsat*catomw / 1000._r8

            ch4_surf_diff_lake(c) = ch4_surf_diff_sat(c)
            ch4_surf_ebul_lake(c) = ch4_surf_ebul_sat(c)
            
            sum1 = 0._r8
            sum2 = 0._r8 
            !$acc loop vector reduction(+:sum1,sum2)
            do j=1,nlevsoi
               sum1 = sum1 + ch4_oxid_depth_sat(c,j)*dz(c,j)*catomw
               sum2 = sum2 + ch4_prod_depth_sat(c,j)*dz(c,j)*catomw
            end do 

            ch4_oxid_tot(c) = ch4_oxid_tot(c) + sum1
            ch4_prod_tot(c) = ch4_prod_tot(c) + sum2

            ! Adjustment to NEE flux to atm. for methane production
            if (.not. replenishlakec) then
               nem_col(c) = nem_col(c) + ch4_prod_tot(c)
               ! Here this is positive because it is actually the CO2 that comes off with the methane
               ! NOTE THIS MODE ASSUMES TRANSIENT CARBON SUPPLY FROM LAKES; COUPLED MODEL WILL NOT CONSERVE CARBON
               ! IN THIS MODE.
            else ! replenishlakec
               nem_col(c) = nem_col(c) - ch4_prod_tot(c)
               ! Keep total C constant, just shift from CO2 to methane
            end if

            ! Adjustment to NEE flux to atm. for methane oxidation
            nem_col(c) = nem_col(c) + ch4_oxid_tot(c)
         end do

         !$acc parallel loop independent gang vector collapse(2) default(present) 
         do j=1,nlevsoi 
            do fc = 1, num_lakec 
               c = filter_lakec(fc) 
               if (.not. replenishlakec) then
                  !Adjust lake_soilc for production.
                  lake_soilc(c,j) = lake_soilc(c,j) - 2._r8*ch4_prod_depth_sat(c,j)*dtime_mod*catomw
                  ! Factor of 2 is for CO2 that comes off with CH4 because of stoichiometry
               end if
               !Set lake diagnostic output variables
               ch4_prod_depth_lake(c,j) = ch4_prod_depth_sat(c,j)
               conc_ch4_lake(c,j)       = conc_ch4_sat(c,j)
               conc_o2_lake(c,j)        = conc_o2_sat(c,j)
               ch4_oxid_depth_lake(c,j) = ch4_oxid_depth_sat(c,j)
            end do 
         end do 

      end if  ! ch4_surf_flux_tot, ch4_oxid_tot, and ch4_prod_tot should be initialized to 0 above if .not. allowlakeprod

      ! Finalize CH4 balance and check for errors
      erridx = 0 
      !$acc parallel loop independent gang worker default(present) private(sum1) copy(erridx)
      do fc = 1, num_soilc
         c = filter_soilc(fc)
         sum1 = 0._r8

         !$acc loop vector reduction(+:sum1) 
         do j=1,nlevsoi
            sum1 = sum1 + (finundated(c)*conc_ch4_sat(c,j) + (1._r8-finundated(c))*conc_ch4_unsat(c,j))*dz(c,j)*catomw
            ! mol CH4 --> g C
         end do
         totcolch4(c) = totcolch4(c) + sum1 

         if (totsoilcolch4_bef(fc) /= spval) then ! not first timestep
            ! Check balance
            errch4 = totcolch4(c) - totsoilcolch4_bef(fc) - dtime_mod*(ch4_prod_tot(c) - ch4_oxid_tot(c) &
                 - ch4_surf_flux_tot(c)*1000._r8) ! kg C --> g C
            
            if (abs(errch4) > 1.e-7_r8) then ! g C / m^2 / timestep
               erridx = c
            end if
         end if
      end do 
      if(erridx > 0) then 
         c = erridx 
         g = col_pp%gridcell(erridx)
         write(iulog, *) 'CH4 Conservation Error in CH4Mod driver, nstep, c, errch4 (gC /m^2.timestep)',nstep,c
         write(iulog,*) 'Latdeg,Londeg=',grc_pp%latdeg(g),grc_pp%londeg(g)
         ! call endrun(msg=' ERROR: Methane conservation error'//errMsg(__FILE__, __LINE__))
      end if 

      if (allowlakeprod) then
         erridx = 0
         !$acc parallel loop independent gang worker default(present) private(sum1,errch4) copy(erridx)
         do fc = 1, num_lakec
            c = filter_lakec(fc)
            sum1 = 0._r8

            !$acc loop vector reduction(+:sum1) 
            do j=1,nlevsoi
               sum1 = sum1 + conc_ch4_sat(c,j)*dz(c,j)*catomw ! mol CH4 --> g C
            end do 
            totcolch4(c) = totcolch4(c) + sum1 

            if (totlakecolch4_bef(fc) /= spval) then ! not first timestep
               ! Check balance
               errch4 = totcolch4(c) - totlakecolch4_bef(fc) - dtime_mod*(ch4_prod_tot(c) - ch4_oxid_tot(c) &
                    - ch4_surf_flux_tot(c)*1000._r8) ! kg C --> g C
               if (abs(errch4) > 1.e-7_r8) then ! g C / m^2 / timestep
                  erridx = c 
               end if
            end if
         end do
         if(erridx > 0) then 
            c = erridx 
            g = col_pp%gridcell(c)
            write(iulog,*)'CH4 Conservation Error in CH4Mod driver for lake column, nstep, c, errch4 (gC/m^2.timestep)',nstep,c
            write(iulog,*) 'Latdeg,Londeg=',grc_pp%latdeg(g),grc_pp%londeg(g)
            ! call endrun(msg=' ERROR: Methane conservation error, allowlakeprod'//errMsg(__FILE__, __LINE__))
         end if 

      end if

      ! Now average up to gridcell for fluxes
      call c2g_1d_parallel( bounds, ch4_oxid_tot(begc:endc), ch4co2f(begg:endg),&
           c2l_scale_type=unity, l2g_scale_type=unity,para=.true.)

      call c2g_1d_parallel( bounds, ch4_prod_tot(begc:endc), ch4prodg(begg:endg),&
           c2l_scale_type=unity, l2g_scale_type=unity,para=.true.)

      call c2g_1d_parallel( bounds, nem_col(begc:endc), nem_grc(begg:endg), &
           c2l_scale_type=unity, l2g_scale_type=unity,para=.true.)

    !$acc exit data delete(&
    !$acc jwt(:), &
    !$acc ch4_prod_tot(:), &
    !$acc ch4_oxid_tot(:), &
    !$acc nem_col(:), &
    !$acc rootfraction(:,:), &
    !$acc totsoilcolch4_bef(:), &
    !$acc totlakecolch4_bef(:), &
    !$acc dummyfilter(:),sum1,sum2,sum3)
   !$acc exit data delete(redoxlags,redoxlags_vertical,rgasm, dtime_ch4, lake, sat)

    end associate

  end subroutine CH4

  !-----------------------------------------------------------------------
  subroutine ch4_prod (bounds, num_methc, filter_methc, num_methp, &
       filter_methp, jwt, sat, lake, soilstate_vars, ch4_vars, dtime, &
       conc_o2, ch4_prod_depth, o2_decomp_depth, co2_decomp_depth )
    !
    ! !DESCRIPTION:
    ! Production is done below the water table, based on CN heterotrophic respiration.
    ! O2 is consumed by roots & by heterotrophic aerobes.
    ! Production is done separately for sat & unsat, and is adjusted for temperature, seasonal inundation,
    ! pH (optional), & redox lag factor.
    !
    ! !USES:
    use CH4varcon          , only: usephfact, anoxicmicrosites, ch4rmcnlim
    use elm_varctl         , only: anoxia
    use elm_varpar         , only: nlevdecomp, nlevdecomp_full
    use SharedParamsMod    , only: nlev_soildecomp_standard
    use pftvarcon          , only: noveg
    !
    ! !ARGUMENTS:
    type(bounds_type)       , intent(in)    :: bounds
    integer                 , intent(in)    :: num_methc           ! number of column soil points in column filter
    integer                 , intent(in)    :: filter_methc(:)     ! column filter for soil points
    integer                 , intent(in)    :: num_methp           ! number of soil points in pft filter
    integer                 , intent(in)    :: filter_methp(:)     ! patch filter for soil points
    integer                 , intent(in)    :: jwt( bounds%begc: ) ! index of the soil layer right above the water table (-) [col]
    integer                 , intent(in)    :: sat                 ! 0 = unsaturated; 1 = saturated
    logical                 , intent(in)    :: lake                ! function called with lake filter
    type(soilstate_type)    , intent(inout) :: soilstate_vars
    type(ch4_type)          , intent(inout) :: ch4_vars
    real(r8), intent(in) :: dtime
    real(r8), intent(inout) :: conc_o2(:,:)           ! Input:  [real(r8) (:,:)]  O2 conc in each soil layer (mol/m3) (nlevsoi)
    real(r8), intent(inout) :: ch4_prod_depth(:,:)    ! Output: [real(r8) (:,:)]  production of CH4 in each soil layer (nlevsoi) (mol/m3/s)
    real(r8), intent(inout) :: o2_decomp_depth(:,:)   ! Output: [real(r8) (:,:)]  O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
    real(r8), intent(inout) :: co2_decomp_depth(:,:)  ! Output: [real(r8) (:,:)]  CO2 production during decomposition in each soil layer (nlevsoi) (mol/m3/s)
    !
    ! !LOCAL VARIABLES:
    integer  :: p,c,j,g        ! indices
    integer  :: fc             ! column index
    integer  :: fp             ! PATCH index
    real(r8) :: base_decomp    ! base rate (mol/m2/s)
    real(r8) :: q10lake          
    real(r8) :: partition_z      
    real(r8) :: pH_fact_ch4    ! pH factor in methane production
    ! Factors for methanogen temperature dependence being greater than soil aerobes
    real(r8) :: f_ch4_adj      ! Adjusted f_ch4
    real(r8) :: t_fact_ch4     ! Temperature factor calculated using additional Q10
    ! O2 limitation on decomposition and methanogenesis
    real(r8) :: seasonalfin    ! finundated in excess of respiration-weighted annual average
    ! For calculating column average (rootfrac(p,j)*rr(p,j))
    real(r8) :: rr_vr(1:num_methc, 1:nlevsoi) ! vertically resolved column-mean root respiration (g C/m^2/s)
    real(r8) :: sum1 
    !-----------------------------------------------------------------------
    associate(                                        &
         wtcol          =>    veg_pp%wtcol          , & ! Input:  [real(r8) (:)    ]  weight (relative to column)
         dz             =>    col_pp%dz             , & ! Input:  [real(r8) (:,:)  ]  layer thickness (m)  (-nlevsno+1:nlevsoi)
         z              =>    col_pp%z              , & ! Input:  [real(r8) (:,:)  ]  layer depth (m) (-nlevsno+1:nlevsoi)
         zi             =>    col_pp%zi             , & ! Input:  [real(r8) (:,:)  ]  interface level below a "z" level (m)
 
         t_soisno       =>    col_es%t_soisno       , & ! Input:  [real(r8) (:,:)  ]  soil temperature (Kelvin)  (-nlevsno+1:nlevsoi)

         h2osoi_vol     =>    col_ws%h2osoi_vol     , & ! Input:  [real(r8) (:,:)  ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]

         rr             =>    veg_cf%rr             , & ! Input:  [real(r8) (:)    ]  (gC/m2/s) root respiration (fine root MR + total root GR)
         somhr          =>    col_cf%somhr          , & ! Input:  [real(r8) (:)    ]  (gC/m2/s) soil organic matter heterotrophic respiration
         lithr          =>    col_cf%lithr          , & ! Input:  [real(r8) (:)    ]  (gC/m2/s) litter heterotrophic respiration
         hr_vr          =>    col_cf%hr_vr          , & ! Input:  [real(r8) (:,:)  ]  total vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
         o_scalar       =>    col_cf%o_scalar       , & ! Input:  [real(r8) (:,:)  ]  fraction by which decomposition is limited by anoxia
         col_rr         =>    col_cf%rr             , & ! Input:  [real(r8) (:)    ]  (gC/m2/s) root respiration (fine root MR + total root GR)
         fphr           =>    col_cf%fphr           , & ! Input:  [real(r8) (:,:)  ]  fraction of potential heterotrophic respiration

         pot_f_nit_vr   =>    col_nf%pot_f_nit_vr , & ! Input:  [real(r8) (:,:)  ]  (gN/m3/s) potential soil nitrification flux

         watsat         =>    soilstate_vars%watsat_col      , & ! Input:  [real(r8) (:,:)  ]  volumetric soil water at saturation (porosity)
         rootfr         =>    soilstate_vars%rootfr_patch    , & ! Input:  [real(r8) (:,:)  ]  fraction of roots in each soil layer  (nlevsoi)
         rootfr_col     =>    soilstate_vars%rootfr_col      , & ! Input:  [real(r8) (:,:)  ]  fraction of roots in each soil layer  (nlevsoi)

         finundated     =>    ch4_vars%finundated_col        , & ! Input:  [real(r8) (:)    ]  fractional inundated area in soil column
         pH             =>    ch4_vars%pH_col                , & ! Input:  [real(r8) (:)    ]  soil water pH
         lake_soilc     =>    ch4_vars%lake_soilc_col        , & ! Input:  [real(r8) (:,:)  ]  total soil organic matter found in level (g C / m^3) (nlevsoi)
         annavg_finrw   =>    ch4_vars%annavg_finrw_col      , & ! Input:  [real(r8) (:)    ]  respiration-weighted annual average of finundated
         finundated_lag =>    ch4_vars%finundated_lag_col    , & ! Input:  [real(r8) (:)    ]  time-lagged fractional inundated area
         layer_sat_lag  =>    ch4_vars%layer_sat_lag_col     , & ! Input:  [real(r8) (: ,:) ]  Lagged saturation status of soil layer in the unsaturated zone (1 = sat)
         sif            =>    ch4_vars%sif_col               , & ! Output: [real(r8) (:)    ]  (unitless) ratio applied to sat. prod. to account for seasonal inundation
         q10ch4           => CH4ParamsInst%q10ch4      , &  ! additional Q10 for methane production ABOVE the soil decomposition temperature relationship
         q10ch4base       => CH4ParamsInst%q10ch4base  , &  ! temperature at which the effective f_ch4 actually equals the constant f_ch4
         f_ch4            => CH4ParamsInst%f_ch4       , &  ! ratio of CH4 production to total C mineralization
         rootlitfrac      => CH4ParamsInst%rootlitfrac , &  ! Fraction of soil organic matter associated with roots
         cnscalefactor    => CH4ParamsInst%cnscalefactor, & ! scale factor on CN decomposition for assigning methane flux
         lake_decomp_fact => CH4ParamsInst%lake_decomp_fact , & ! Base decomposition rate (1/s) at 25C
         ! added by Lei Meng to account for pH influence of CH4 production
         pHmax            => CH4ParamsInst%pHmax , &
         pHmin            => CH4ParamsInst%pHmin , &
         oxinhib          => CH4ParamsInst%oxinhib , &    ! inhibition of methane production by oxygen (m^3/mol)
         q10lakebase      => CH4ParamsInst%q10lakebase ,& ! (K) base temperature for lake CH4 production
         ! Shared constant with other modules
         mino2lim => ParamsShareInst%mino2lim & ! minimum anaerobic decomposition rate as a fraction of potential aerobic rate
         )

      !$acc enter data create(&
      !$acc rr_vr(:,:), &
      !$acc base_decomp, &
      !$acc f_ch4_adj)

      q10lake = q10ch4 * 1.5_r8 ! For now, take to be the same as q10ch4 * 1.5.
      !$acc enter data copyin(q10lake)

      ! PATCH loop to calculate vertically resolved column-averaged root respiration
      if (.not. lake) then
         ! rr_vr(bounds%begc:bounds%endc,:) = spval

         !$acc parallel loop independent collapse(2) gang worker default(present) private(sum1) 
         do j=1,nlevsoi
            do fc = 1, num_methc
               c = filter_methc(fc)
               sum1 = 0._r8 
               !$acc loop vector reduction(+:sum1) 
               do p = col_pp%pfti(c), col_pp%pftf(c) 
                  if (wtcol(p) > 0._r8 .and. veg_pp%itype(p) /= noveg) then
                     sum1 = sum1 + rr(p)*rootfr(p,j)*wtcol(p)
                  end if
               end do
               rr_vr(fc,j) = sum1  
            end do
         end do
      end if



      ! column loop to partition decomposition_rate into each soil layer
      !$acc parallel loop independent gang vector default(present) collapse(2) 
      do j=1,nlevsoi
         do fc = 1, num_methc
            c = filter_methc (fc)
            g = col_pp%gridcell(c)
            partition_z = 1._r8
            base_decomp = 0.0_r8
            
            if (.not. lake) then

               if (use_cn) then
                  ! Use soil heterotrophic respiration (based on Wania)
                  base_decomp = (somhr(c)+lithr(c)) / catomw
                  ! Convert from gC to molC
                  ! Multiply base_decomp by factor accounting for lower carbon stock in seasonally inundated areas than
                  ! if it were inundated all year.
                  ! This is to reduce emissions in seasonally inundated zones, because the eq.
                  ! C-flux will be less than predicted by a non-O2-lim model
                  if (sat == 1) then
                     sif(c) = 1._r8
                     if (.not. anoxia) then
                        if (annavg_finrw(c) /= spval) then
                           seasonalfin = max(finundated(c)-annavg_finrw(c), 0._r8)
                           if (seasonalfin > 0._r8) then
                              sif(c) = (annavg_finrw(c) + mino2lim*seasonalfin) / finundated(c)
                              base_decomp = base_decomp * sif(c)
                           end if
                        end if
                     end if ! anoxia
                  end if
               else
                  #ifndef _OPENACC 
                  call endrun(msg=' ERROR: No source for decomp rate in CH4Prod.'//&
                       ' CH4 model currently requires CN.'//errMsg(__FILE__, __LINE__))
                  #endif 
               end if ! use_cn

               ! For sensitivity studies
               base_decomp = base_decomp * cnscalefactor

            else !lake

               base_decomp = lake_decomp_fact * lake_soilc(c,j) * dz(c,j) * &
                    q10lake**( (t_soisno(c,j)-q10lakebase)/10._r8) / catomw
               ! convert from g C to mol C
            end if

            ! For all landunits, prevent production or oxygen consumption when soil is at or below freezing.
            ! If using VERTSOILC, it is OK to use base_decomp as given because liquid water stress will limit decomp.
            if (t_soisno(c,j) <= tfrz .and. (nlevdecomp == 1 .or. lake)) base_decomp = 0._r8

            ! depth dependence of production either from rootfr or decomp model
            if (.not. lake) then ! use default rootfr, averaged to the column level in the ch4 driver, or vert HR
               if (nlevdecomp == 1) then ! not VERTSOILC
                  if (j <= nlev_soildecomp_standard) then  ! Top 5 levels are also used in the CLM code for establishing temperature
                     ! and moisture constraints on SOM activity
                     partition_z = rootfr_col(c,j)*rootlitfrac + (1._r8 - rootlitfrac)*dz(c,j)/zi(c,nlev_soildecomp_standard)
                  else
                     partition_z = rootfr_col(c,j)*rootlitfrac
                  end if
               else
                  if ( (somhr(c) + lithr(c)) > 0._r8) then
                     partition_z = hr_vr(c,j) * dz(c,j) / (somhr(c) + lithr(c))
                  else
                     partition_z = 1._r8
                  end if
               end if
            else ! lake
               partition_z = 1._r8
            endif

            ! Adjust f_ch4 to account for the fact that methanogens may have a higher Q10 than aerobic decomposers.
            ! Note this is crude and should ideally be applied to all anaerobic decomposition rather than just the
            ! f_ch4.
            f_ch4_adj = 1.0_r8
            if (.not. lake) then
               t_fact_ch4 = q10ch4**((t_soisno(c,j) - q10ch4base)/10._r8)
               ! Adjust f_ch4 by the ratio
               f_ch4_adj = f_ch4 * t_fact_ch4

               ! Remove CN nitrogen limitation, as methanogenesis is not N limited.
               ! Also remove (low) moisture limitation
               if (ch4rmcnlim) then
                  if (j > nlevdecomp) then
                     if (fphr(c,1) > 0._r8) then
                        f_ch4_adj = f_ch4_adj / fphr(c,1)
                     end if
                  else ! j == 1 or VERTSOILC
                     if (fphr(c,j) > 0._r8) then
                        f_ch4_adj = f_ch4_adj / fphr(c,j)
                     end if
                  end if
               end if

            else ! lake
               f_ch4_adj = 0.5_r8 ! For lakes assume no redox limitation. Production only depends on temp, soil C, and
               ! lifetime parameter.
            end if

            ! If switched on, use pH factor for production based on spatial pH data defined in surface data.
            if (.not. lake .and. usephfact .and. pH(c) >  pHmin .and.pH(c) <  pHmax) then
               pH_fact_ch4 = 10._r8**(-0.2235_r8*pH(c)*pH(c) + 2.7727_r8*pH(c) - 8.6_r8)
               ! fitted function using data from Dunfield et al. 1993
               ! Strictly less than one, with optimum at 6.5
               ! From Lei Meng
               f_ch4_adj = f_ch4_adj * pH_fact_ch4
            else
               ! if no data, then no pH effects
            end if

            ! Redox factor
            if ( (.not. lake) .and. sat == 1 .and. finundated_lag(c) < finundated(c)) then
               f_ch4_adj = f_ch4_adj * finundated_lag(c) / finundated(c)
            else if (sat == 0 .and. j > jwt(c)) then ! Assume lag in decay of alternative electron acceptors vertically
               f_ch4_adj = f_ch4_adj * layer_sat_lag(c,j)
            end if
            ! Alternative electron acceptors will be consumed first after soil is inundated.

            f_ch4_adj = min(f_ch4_adj, 0.5_r8)
            ! Must be less than 0.5 because otherwise the actual implied aerobic respiration would be negative.
            ! The total of aer. respiration + methanogenesis must remain equal to the SOMHR calculated in CN,
            ! so that the NEE is sensible. Even perfectly anaerobic conditions with no alternative
            ! electron acceptors would predict no more than 0.5 b/c some oxygen is present in organic matter.
            ! e.g. 2CH2O --> CH4 + CO2.

            ! Decomposition uses 1 mol O2 per mol CO2 produced (happens below WT also, to deplete O2 below WT)
            ! o2_decomp_depth is the demand in the absense of O2 supply limitation, in addition to autotrophic respiration.
            ! Competition will be done in ch4_oxid

            o2_decomp_depth(c,j) = base_decomp * partition_z / dz(c,j)
            if (anoxia) then
               ! Divide off o_scalar to use potential O2-unlimited HR to represent aerobe demand for oxygen competition
               if (.not. lake .and. j > nlevdecomp) then
                  if (o_scalar(c,1) > 0._r8) then
                     o2_decomp_depth(c,j) = o2_decomp_depth(c,j) / o_scalar(c,1)
                  end if
               else if (.not. lake) then ! j == 1 or VERTSOILC
                  if (o_scalar(c,j) > 0._r8) then
                     o2_decomp_depth(c,j) = o2_decomp_depth(c,j) / o_scalar(c,j)
                  end if
               end if
            end if ! anoxia

            ! Add root respiration
            if (.not. lake) then
               !o2_decomp_depth(c,j) = o2_decomp_depth(c,j) + col_rr(c)*rootfr(c,j)/catomw/dz(c,j) ! mol/m^3/s
               o2_decomp_depth(c,j) = o2_decomp_depth(c,j) + rr_vr(fc,j)/catomw/dz(c,j) ! mol/m^3/s
               ! g C/m2/s ! gC/mol O2 ! m
            end if

            ! Add oxygen demand for nitrification
            if (.not. lake .and. j<= nlevdecomp_full ) then
               o2_decomp_depth(c,j) = o2_decomp_depth(c,j) + pot_f_nit_vr(c,j) * 2.0_r8/14.0_r8
               ! g N/m^3/s           mol O2 / g N
            end if

            if (j  >  jwt(c)) then ! Below the water table so anaerobic CH4 production can occur
               ! partition decomposition to layer
               ! turn into per volume-total by dz
               ch4_prod_depth(c,j) = f_ch4_adj * base_decomp * partition_z / dz (c,j)! [mol/m3-total/s]
            else ! Above the WT
               if (anoxicmicrosites) then
                  ch4_prod_depth(c,j) = f_ch4_adj * base_decomp * partition_z / dz (c,j) &
                       / (1._r8 + oxinhib*conc_o2(c,j))
               else
                  ch4_prod_depth(c,j) = 0._r8 ! [mol/m3 total/s]
               endif ! anoxicmicrosites
            endif ! WT

         end do ! fc
      end do ! nlevsoi

      !$acc exit data delete(&
      !$acc rr_vr(:,:), &
      !$acc base_decomp, &
      !$acc f_ch4_adj,q10lake)

    end associate

  end subroutine ch4_prod

  !-----------------------------------------------------------------------
  subroutine ch4_oxid (bounds, &
       num_methc, filter_methc, &
       jwt, sat, lake, soilstate_vars, dtime,&
       ch4_oxid_depth, o2_oxid_depth, co2_oxid_depth, &   
       conc_ch4, conc_o2, o2_decomp_depth  )
    !
    ! !DESCRIPTION:
    ! Oxidation is based on double Michaelis-Mentin kinetics, and is adjusted for low soil moisture.
    ! Oxidation will be limited by available oxygen in ch4_tran.

    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in) :: bounds
    integer                , intent(in) :: num_methc           ! number of column soil points in column filter
    integer                , intent(in) :: filter_methc(:)     ! column filter for soil points
    integer                , intent(in) :: jwt( bounds%begc: ) ! index of the soil layer right above the water table (-) [col]
    integer                , intent(in) :: sat                 ! 0 = unsaturated; 1 = saturated
    logical                , intent(in) :: lake                ! function called with lake filter
    type(soilstate_type)   , intent(inout) :: soilstate_vars
    real(r8), intent(in) :: dtime                         ! land model time step (sec)
    real(r8) , intent(inout) :: ch4_oxid_depth(:,:) 
    real(r8) , intent(inout) :: o2_oxid_depth (:,:) 
    real(r8) , intent(inout) :: co2_oxid_depth(:,:) 
    real(r8) , intent(inout) :: conc_ch4(:,:) 
    real(r8) , intent(inout) :: conc_o2 (:,:) 
    real(r8) , intent(inout) :: o2_decomp_depth(:,:)
    !
    ! !LOCAL VARIABLES:
    integer :: c,j                            ! indices
    integer :: fc                             ! column index
    real(r8), parameter :: t0 = tfrz+12._r8 ! Base temperature for Q10  ! Walter, for Michigan site where the 45 M/h comes from
    real(r8):: porevol                        ! air-filled volume ratio to total soil volume
    real(r8):: h2osoi_vol_min                 ! h2osoi_vol restricted to be below watsat
    real(r8):: conc_ch4_rel                   ! concentration with respect to water volume (mol/m^3 water)
    real(r8):: conc_o2_rel                    ! concentration with respect to water volume (mol/m^3 water)
    real(r8):: oxid_a                         ! Oxidation predicted by method A (temperature & enzyme limited) (mol CH4/m3/s)
    real(r8):: smp_fact                       ! factor for reduction based on soil moisture (unitless)
    real(r8):: porewatfrac                    ! fraction of soil pore space that is filled with water
    real(r8):: k_h_cc, k_h_inv                ! see functions below for description
    real(r8):: k_m_eff                        ! effective k_m
    real(r8):: vmax_eff                       ! effective vmax
                                              ! ch4 oxidation parameters
    real(r8) :: vmax_ch4_oxid                 ! oxidation rate constant (= 45.e-6_r8 * 1000._r8 / 3600._r8) [mol/m3-w/s];
    real(r8) :: k_m 			      ! Michaelis-Menten oxidation rate constant for CH4 concentration
    real(r8) :: q10_ch4oxid                   ! Q10 oxidation constant
    real(r8) :: smp_crit                      ! Critical soil moisture potential
    real(r8) :: k_m_o2                        ! Michaelis-Menten oxidation rate constant for O2 concentration
    real(r8) :: k_m_unsat                     ! Michaelis-Menten oxidation rate constant for CH4 concentration
    real(r8) :: vmax_oxid_unsat               ! (= 45.e-6_r8 * 1000._r8 / 3600._r8 / 10._r8) [mol/m3-w/s]
    !
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes

    associate(                                          &
         h2osoi_vol => col_ws%h2osoi_vol , & ! Input:  [real(r8) (:,:)  ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]

         smp_l      => soilstate_vars%smp_l_col       , & ! Input:  [real(r8) (: ,:) ]  soil matrix potential [mm]
         watsat     => soilstate_vars%watsat_col      , & ! Input:  [real(r8) (:,:)  ]  volumetric soil water at saturation (porosity)

         t_soisno   => col_es%t_soisno   , & ! Input:  [real(r8) (:,:)  ]  soil temperature (Kelvin)  (-nlevsno+1:nlevsoi)
            ! Set oxidation parameters
         vmax_ch4_oxid   => CH4ParamsInst%vmax_ch4_oxid , &
         k_m             => CH4ParamsInst%k_m , &
         q10_ch4oxid     => CH4ParamsInst%q10_ch4oxid , &
         smp_crit        => CH4ParamsInst%smp_crit , &
         k_m_o2          => CH4ParamsInst%k_m_o2   , &
         k_m_unsat       => CH4ParamsInst%k_m_unsat, &
         vmax_oxid_unsat => CH4ParamsInst%vmax_oxid_unsat &
         )


      ! Loop to determine oxidation in each layer
      !$acc parallel loop independent gang vector default(present) collapse(2) 
      do j=1,nlevsoi
         do fc = 1, num_methc
            c = filter_methc(fc)

            if (sat == 1 .or. j > jwt(c)) then
               ! Literature (e.g. Bender & Conrad, 1992) suggests lower k_m and vmax for high-CH4-affinity methanotrophs in
               ! upland soils consuming ambient methane.
               k_m_eff = k_m
               vmax_eff = vmax_ch4_oxid
            else
               k_m_eff = k_m_unsat
               vmax_eff = vmax_oxid_unsat
            end if

            porevol = max(watsat(c,j) - h2osoi_vol(c,j), 0._r8)
            h2osoi_vol_min = min(watsat(c,j), h2osoi_vol(c,j))
            if (j <= jwt(c) .and. smp_l(c,j) < 0._r8) then
               smp_fact = exp(-smp_l(c,j)/smp_crit)
               ! Schnell & King, 1996, Figure 3
            else
               smp_fact = 1._r8
            end if

            if (j  <=  jwt(c)) then ! Above the water table
               k_h_inv = exp(-c_h_inv(1) * (1._r8 / t_soisno(c,j) - 1._r8 / kh_tbase) + log (kh_theta(1)))
               k_h_cc = t_soisno(c,j) / k_h_inv * rgasLatm ! (4.21) Wania [(mol/m3w) / (mol/m3g)]
               conc_ch4_rel = conc_ch4(c,j) / (h2osoi_vol_min + porevol/k_h_cc)

               k_h_inv = exp(-c_h_inv(2) * (1._r8 / t_soisno(c,j) - 1._r8 / kh_tbase) + log (kh_theta(2)))
               k_h_cc = t_soisno(c,j) / k_h_inv * rgasLatm ! (4.21) Wania [(mol/m3w) / (mol/m3g)]
               conc_o2_rel  = conc_o2(c,j) / (h2osoi_vol_min + porevol/k_h_cc)
            else
               conc_ch4_rel = conc_ch4(c,j) / watsat(c,j)
               conc_o2_rel  = conc_o2(c,j) / watsat(c,j)
            endif

            ![mol/m3-t/s]         [mol/m3-w/s]    [m3-w/m3-t]     [mol/m3-w]    [mol/m3-w]  [mol/m3-w]
            oxid_a              = vmax_eff     * h2osoi_vol_min* conc_ch4_rel / (k_m_eff + conc_ch4_rel) &
                 * conc_o2_rel / (k_m_o2 + conc_o2_rel) &
                 * q10_ch4oxid ** ((t_soisno(c,j) - t0) / 10._r8) * smp_fact

            ! For all landunits / levels, prevent oxidation if at or below freezing
            if (t_soisno(c,j) <= tfrz) oxid_a = 0._r8

            ch4_oxid_depth(c,j) = oxid_a
            o2_oxid_depth(c,j) = ch4_oxid_depth(c,j) * 2._r8

         end do
      end do

    end associate
  end subroutine ch4_oxid

  !-----------------------------------------------------------------------
  subroutine ch4_aere (bounds, &
       num_methc, filter_methc, &
       num_methp, filter_methp, &
       jwt, sat, lake, &
       canopystate_vars, soilstate_vars, energyflux_vars, ch4_vars, dtime, &
       ch4_aere_depth, ch4_tran_depth, o2_aere_depth , &
       co2_aere_depth, conc_ch4, conc_o2, ch4_oxid_depth, ch4_prod_depth)
    !
    ! !DESCRIPTION:
    ! Arctic c3 grass (which is often present in fens) and all vegetation in inundated areas is assumed to have
    ! some root porosity. Currently, root porosity is allowed to be different for grasses & non-grasses.
    ! CH4 diffuses out and O2 diffuses into the soil.  CH4 is also lossed via transpiration, which is both
    ! included in the "aere" variables and output separately.  In practice this value is small.
    ! By default upland veg. has small 5% porosity but this can be switched to be equal to inundated porosity.

    ! !USES:
    use elm_varcon       , only : rpi
    use pftvarcon        , only : nc3_arctic_grass, crop, nc3_nonarctic_grass, nc4_grass, noveg
    use CH4varcon        , only : transpirationloss, usefrootc, use_aereoxid_prog
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_methc           ! number of column soil points in column filter
    integer                , intent(in)    :: filter_methc(:)     ! column filter for soil points
    integer                , intent(in)    :: num_methp           ! number of soil points in pft filter
    integer                , intent(in)    :: filter_methp(:)     ! patch filter for soil points
    integer                , intent(in)    :: jwt( bounds%begc: ) ! index of the soil layer right above the water table (-) [col]
    integer                , intent(in)    :: sat                 ! 0 = unsaturated; 1 = saturated
    logical                , intent(in)    :: lake             ! function called with lake filter
    type(canopystate_type) , intent(in)    :: canopystate_vars
    type(soilstate_type)   , intent(inout) :: soilstate_vars
    type(energyflux_type)  , intent(inout) :: energyflux_vars
    type(ch4_type)         , intent(inout) :: ch4_vars
    real(r8)               , intent(in) :: dtime
    real(r8), intent(inout) :: ch4_aere_depth(:,:)
    real(r8), intent(inout) :: ch4_tran_depth(:,:)
    real(r8), intent(inout) :: o2_aere_depth(:,:)
    real(r8), intent(inout) :: co2_aere_depth(:,:)
    real(r8), intent(inout) :: ch4_oxid_depth(:,:)
    real(r8), intent(inout) :: ch4_prod_depth(:,:)
    real(r8), intent(inout) :: conc_o2(:,:)
    real(r8), intent(inout) :: conc_ch4(:,:)
    !
    ! !LOCAL VARIABLES:
    integer  :: p,c,g,j                ! indices
    integer  :: fc,fp                  ! soil filter column index
    real(r8) :: f_oxid                 ! fraction of CH4 oxidized in oxic zone around roots
    real(r8) :: diffus_aere            ! gas diffusivity through aerenchyma (m^2/s)
    real(r8) :: m_tiller
    real(r8) :: n_tiller
    real(r8) :: poros_tiller
    real(r8) :: rob                    ! root obliquity, e.g. csc of root angle relative to vertical
                                       ! (ratio of root total length to depth)
    real(r8) :: area_tiller            ! cross-sectional area of tillers (m^2/m^2)
    real(r8) :: tranloss               ! loss due to transpiration (mol / m3 /s)
    real(r8) :: aere, aeretran, oxaere ! (mol / m3 /s)
    real(r8) :: k_h_cc, k_h_inv, oxdiffus, anpp, nppratio, h2osoi_vol_min, conc_ch4_wat
    real(r8) :: aerecond               ! aerenchyma conductance (m/s)
    ! ch4 aerenchyma parameters
    real(r8) :: aereoxid               ! fraction of methane flux entering aerenchyma rhizosphere
    real(r8) :: scale_factor_aere      ! scale factor on the aerenchyma area for sensitivity tests
    real(r8) :: nongrassporosratio     ! Ratio of root porosity in non-grass to grass, used for aerenchyma transport
    real(r8) :: unsat_aere_ratio       ! Ratio to multiply upland vegetation aerenchyma porosity by compared to inundated systems (= 0.05_r8 / 0.3_r8)
    real(r8) :: porosmin               ! minimum aerenchyma porosity (unitless)(= 0.05_r8)

    real(r8), parameter :: smallnumber = 1.e-12_r8
    real(r8) :: sum1, sum2, sum3
    
    !-----------------------------------------------------------------------
    ! Enforce expected array sizes

    associate(                                                     &
         z             =>    col_pp%z                               , & ! Input:  [real(r8) (:,:)  ]  layer depth (m) (-nlevsno+1:nlevsoi)
         dz            =>    col_pp%dz                              , & ! Input:  [real(r8) (:,:)  ]  layer thickness (m)  (-nlevsno+1:nlevsoi)
         wtcol         =>    veg_pp%wtcol                           , & ! Input:  [real(r8) (:)    ]  weight (relative to column)

         elai          =>    canopystate_vars%elai_patch         , & ! Input:  [real(r8) (:)    ]  one-sided leaf area index with burying by snow

         annsum_npp    =>    veg_cf%annsum_npp    , & ! Input:  [real(r8) (:)    ]  annual sum NPP (gC/m2/yr)

         t_soisno      =>    col_es%t_soisno       , & ! Input:  [real(r8) (:,:)  ]  soil temperature (Kelvin)  (-nlevsno+1:nlevsoi)

         watsat        =>    soilstate_vars%watsat_col           , & ! Input:  [real(r8) (:,:)  ]  volumetric soil water at saturation (porosity)
         rootr         =>    soilstate_vars%rootr_patch          , & ! Input:  [real(r8) (:,:)  ]  effective fraction of roots in each soil layer  (nlevgrnd)
         rootfr        =>    soilstate_vars%rootfr_patch         , & ! Input:  [real(r8) (:,:)  ]  fraction of roots in each soil layer  (nlevsoi)

         h2osoi_vol    =>    col_ws%h2osoi_vol      , & ! Input:  [real(r8) (:,:)  ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]

         qflx_tran_veg =>    veg_wf%qflx_tran_veg  , & ! Input:  [real(r8) (:)    ]  vegetation transpiration (mm H2O/s) (+ = to atm)

         canopy_cond   =>    energyflux_vars%canopy_cond_patch   , & ! Input:  [real(r8) (:)    ]  tracer conductance for canopy [m/s]

         frootc        =>    veg_cs%frootc       , & ! Input:  [real(r8) (:)    ]  (gC/m2) fine root C

         annavg_agnpp  =>    veg_cf%annavg_agnpp  , & ! Input:  [real(r8) (:)    ]  (gC/m2/s) annual average aboveground NPP
         annavg_bgnpp  =>    veg_cf%annavg_bgnpp  , & ! Input:  [real(r8) (:)    ]  (gC/m2/s) annual average belowground NPP

         grnd_ch4_cond =>    ch4_vars%grnd_ch4_cond_patch        , & ! Input:  [real(r8) (:)    ]  tracer conductance for boundary layer [m/s]
         c_atm         =>    ch4_vars%c_atm_grc                  , & ! Input:  [real(r8) (: ,:) ]  CH4, O2, CO2 atmospheric conc  (mol/m3)

         ! Set aerenchyma parameters
         aereoxid           => CH4ParamsInst%aereoxid  ,& 
         scale_factor_aere  => CH4ParamsInst%scale_factor_aere ,&
         nongrassporosratio => CH4ParamsInst%nongrassporosratio ,&
         unsat_aere_ratio   => CH4ParamsInst%unsat_aere_ratio ,&
         porosmin           => CH4ParamsInst%porosmin ,&
         rob                => CH4ParamsInst%rob &
         )


    !$acc enter data create(&
    !$acc ch4_aere_depth(:,:), &
    !$acc ch4_tran_depth(:,:), &
    !$acc o2_aere_depth(:,:), &
    !$acc co2_aere_depth(:,:), &
    !$acc ch4_oxid_depth(:,:), &
    !$acc ch4_prod_depth(:,:), &
    !$acc conc_o2(:,:), &
    !$acc conc_ch4(:,:), &
    !$acc poros_tiller, &
    !$acc tranloss, &
    !$acc aere, &
    !$acc oxaere, &
    !$acc anpp, &
    !$acc aerecond, &
    !$acc sum1, &
    !$acc sum2, &
    !$acc sum3 )

      ! point loop to partition aerenchyma flux into each soil layer
      if (.not. lake) then
         !$acc parallel loop independent gang worker collapse(2) default(present) private(sum1,sum2,sum3,c,g)
         do j=1,nlevsoi
            do fc = 1, num_methc
               c = filter_methc(fc) 
               g = col_pp%gridcell(c)
               sum1 = 0._r8 
               sum2 = 0._r8 
               sum3 = 0._r8 
               !$acc loop vector reduction(+:sum1,sum2,sum3)
               do p = col_pp%pfti(c), col_pp%pftf(c) 
                  if(.not. veg_pp%active(p)) continue 

                  ! This parameter is poorly constrained and should be done on a PFT-specific basis...
                  diffus_aere = d_con_g(1,1)*1.e-4_r8  ! for CH4: m^2/s
                  
                  ! Calculate transpiration loss
                  if (transpirationloss .and. veg_pp%itype(p) /= noveg) then !allow tloss above WT ! .and. j > jwt(c)) then
                     ! Calculate water concentration
                     h2osoi_vol_min = min(watsat(c,j), h2osoi_vol(c,j))
                     k_h_inv = exp(-c_h_inv(1) * (1._r8 / t_soisno(c,j) - 1._r8 / kh_tbase) + log (kh_theta(1)))
                     k_h_cc = t_soisno(c,j) / k_h_inv * rgasLatm
                     conc_ch4_wat = conc_ch4(c,j) / ( (watsat(c,j)-h2osoi_vol_min)/k_h_cc + h2osoi_vol_min)
                     tranloss = conc_ch4_wat *             rootr(p,j)*qflx_tran_veg(p) / dz(c,j) / 1000._r8
                     ! mol/m3/s    mol/m3                                   mm / s         m           mm/m
                     ! Use rootr here for effective per-layer transpiration, which may not be the same as rootfr
                     tranloss = max(tranloss, 0._r8) ! in case transpiration is pathological
                  else
                     tranloss = 0._r8
                  end if

                  ! Calculate aerenchyma diffusion
                  if (j > jwt(c) .and. t_soisno(c,j) > tfrz .and. veg_pp%itype(p) /= noveg) then
                     ! Attn EK: This calculation of aerenchyma properties is very uncertain. Let's check in once all
                     ! the new components are in; if there is any tuning to be done to get a realistic global flux,
                     ! this would probably be the place.  We will have to document clearly in the Tech Note
                     ! any major changes from the Riley et al. 2011 version. (There are a few other minor ones.)

                     anpp = annsum_npp(p) ! g C / m^2/yr
                     anpp = max(anpp, 0._r8) ! NPP can be negative b/c of consumption of storage pools

                     if (annavg_agnpp(p) /= spval .and. annavg_bgnpp(p) /= spval .and. &
                           annavg_agnpp(p) > 0._r8 .and. annavg_bgnpp(p) > 0._r8) then
                        nppratio = annavg_bgnpp(p) / (annavg_agnpp(p) + annavg_bgnpp(p))
                     else
                        nppratio = 0.5_r8
                     end if

                     ! Estimate area of tillers (see Wania thesis)
                     ! m_tiller = anpp * r_leaf_root * lai ! (4.17 Wania)
                     ! m_tiller = 600._r8 * 0.5_r8 * 2._r8  ! used to be 300
                     ! Note: this calculation is based on Arctic graminoids, and should be refined for woody plants, if not
                     ! done on a PFT-specific basis.

                     if (usefrootc) then
                        m_tiller = frootc(p) ! This will yield much smaller aere area.
                     else
                        m_tiller = anpp * nppratio * elai(p)
                     end if

                     n_tiller = m_tiller / 0.22_r8

                     if (veg_pp%itype(p) == nc3_arctic_grass .or. crop(veg_pp%itype(p)) == 1 .or. &
                          veg_pp%itype(p) == nc3_nonarctic_grass .or. veg_pp%itype(p) == nc4_grass) then
                        poros_tiller = 0.3_r8  ! Colmer 2003
                     else
                        poros_tiller = 0.3_r8 * nongrassporosratio
                     end if

                     if (sat == 0) then
                        poros_tiller = poros_tiller * unsat_aere_ratio
                     end if

                     poros_tiller = max(poros_tiller, porosmin)

                     area_tiller = scale_factor_aere * n_tiller * poros_tiller * rpi * 2.9e-3_r8**2._r8 ! (m2/m2)

                     k_h_inv = exp(-c_h_inv(1) * (1._r8 / t_soisno(c,j) - 1._r8 / kh_tbase) + log (kh_theta(1))) ! (4.12) Wania (L atm/mol)
                     k_h_cc = t_soisno(c,j) / k_h_inv * rgasLatm ! (4.21) Wania [(mol/m3w) / (mol/m3g)]
                     aerecond = area_tiller * rootfr(p,j) * diffus_aere / (z(c,j)*rob)
                     ! Add in boundary layer resistance
                     aerecond = 1._r8 / (1._r8/(aerecond+smallnumber) + 1._r8/(grnd_ch4_cond(p)+smallnumber))

                     aere = aerecond * (conc_ch4(c,j)/watsat(c,j)/k_h_cc - c_atm(g,1)) / dz(c,j) ![mol/m3-total/s]
                     !ZS: Added watsat & Henry's const.
                     aere = max(aere, 0._r8) ! prevent backwards diffusion

                     ! Do oxygen diffusion into layer
                     k_h_inv = exp(-c_h_inv(2) * (1._r8 / t_soisno(c,j) - 1._r8 / kh_tbase) + log (kh_theta(2)))
                     k_h_cc = t_soisno(c,j) / k_h_inv * rgasLatm ! (4.21) Wania [(mol/m3w) / (mol/m3g)]
                     oxdiffus = diffus_aere * d_con_g(2,1) / d_con_g(1,1) ! adjust for O2:CH4 molecular diffusion
                     aerecond = area_tiller * rootfr(p,j) * oxdiffus / (z(c,j)*rob)
                     aerecond = 1._r8 / (1._r8/(aerecond+smallnumber) + 1._r8/(grnd_ch4_cond(p)+smallnumber))
                     oxaere = -aerecond *(conc_o2(c,j)/watsat(c,j)/k_h_cc - c_atm(g,2)) / dz(c,j) ![mol/m3-total/s]
                     oxaere = max(oxaere, 0._r8)
                     ! Diffusion in is positive; prevent backwards diffusion
                     if ( .not. use_aereoxid_prog ) then ! fixed aere oxid proportion; will be done in ch4_tran
                        oxaere = 0._r8
                     end if
                  else
                     aere = 0._r8
                     oxaere = 0._r8
                  end if ! veg type, below water table, & above freezing

                  ! Impose limitation based on available methane during timestep
                  ! By imposing the limitation here, don't allow aerenchyma access to methane from other Patches.
                  aeretran = min(aere+tranloss, conc_ch4(c,j)/dtime + ch4_prod_depth(c,j))
                  sum1 = sum1 + aeretran*wtcol(p) !pft weight in col.
                  sum2 = sum2 + min(tranloss, aeretran)*wtcol(p)
                  sum3 = sum3 + oxaere*wtcol(p)
               end do ! p loop
               ch4_aere_depth (c, j) = sum1 
               ch4_tran_depth (c, j) = sum2
               o2_aere_depth  (c, j) = sum3 
            end do ! c filter
         end do ! over levels
      end if ! not lake


    !$acc exit data delete(&
    !$acc ch4_aere_depth(:,:), &
    !$acc ch4_tran_depth(:,:), &
    !$acc o2_aere_depth(:,:), &
    !$acc co2_aere_depth(:,:), &
    !$acc ch4_oxid_depth(:,:), &
    !$acc ch4_prod_depth(:,:), &
    !$acc conc_o2(:,:), &
    !$acc conc_ch4(:,:), &
    !$acc poros_tiller, &
    !$acc tranloss, &
    !$acc aere, &
    !$acc oxaere, &
    !$acc anpp, &
    !$acc aerecond, &
    !$acc sum1, &
    !$acc sum2, &
    !$acc sum3 )

    end associate

  end subroutine ch4_aere

  !-----------------------------------------------------------------------
  subroutine ch4_ebul (bounds, &
       num_methc, filter_methc, &
       jwt, sat, lake, &
       lakestate_vars, soilstate_vars, ch4_vars, dtime, &
       ch4_ebul_depth, ch4_ebul_total, conc_ch4, &    
       ch4_aere_depth, ch4_oxid_depth)
    !
    ! !DESCRIPTION:
    ! Bubbling is based on temperature & pressure dependent solubility (k_h_cc),
    ! with assumed proportion of bubbles
    ! which are CH4, and assumed early nucleation at vgc_max sat (Wania).
    ! Bubbles are released to the water table surface in ch4_tran.

    ! !USES:
    use LakeCon
    !
    ! !ARGUMENTS:
    type(bounds_type)    , intent(in)    :: bounds
    integer              , intent(in)    :: num_methc           ! number of column soil points in column filter
    integer              , intent(in)    :: filter_methc(:)     ! column filter for soil points
    integer              , intent(in)    :: jwt( bounds%begc: ) ! index of the soil layer right above the water table (-) [col]
    integer              , intent(in)    :: sat                 ! 0 = unsaturated; 1 = saturated
    logical              , intent(in)    :: lake                ! function called with lake filter
    type(lakestate_type) , intent(in)    :: lakestate_vars
    type(soilstate_type) , intent(in)    :: soilstate_vars
    type(ch4_type)       , intent(inout) :: ch4_vars
    real(r8) ,intent(in)  :: dtime   ! land model time step (sec)
    real(r8), intent(inout) :: ch4_ebul_depth(:,:) ! Output: [real(r8) (:,:)]  CH4 loss rate via ebullition in each soil layer (mol/m3/s) (nlevsoi)
    real(r8), intent(inout) :: ch4_ebul_total(:)   ! Output: [real(r8) (:)]  Total column CH4 ebullition (mol/m2/s)
    real(r8), intent(inout) :: conc_ch4(:,:)       ! Output: [real(r8) (:,:)]  CH4 conc in each soil layer (mol/m3) (nlevsoi)
    real(r8), intent(in   ) :: ch4_aere_depth(:,:) ! Input:  [real(r8) (:,:)]  CH4 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
    real(r8), intent(in   ) :: ch4_oxid_depth(:,:) ! Input:  [real(r8) (:,:)]  CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
    !
    ! !LOCAL VARIABLES:
    integer :: c,j,t    ! indices
    integer :: fc       ! soil filter column index
    integer :: fp       ! soil filter pft index
    real(r8) :: vgc     ! volumetric CH4 content (m3 CH4/m3 pore air)
    real(r8) :: vgc_min ! minimum aqueous CH4 content when ebullition ceases
    real(r8) :: k_h_inv !
    real(r8) :: k_h     !
    real(r8) :: k_h_cc  !
    real(r8) :: pressure! sum atmospheric and hydrostatic pressure
    real(r8),parameter :: bubble_f = 0.57_r8 ! CH4 content in gas bubbles (Kellner et al. 2006)

    real(r8) :: ebul_timescale
    real(r8) :: vgc_max   ! ratio of saturation pressure triggering ebullition
    
    !-----------------------------------------------------------------------
    ! Enforce expected array sizes

    associate(                                                      &
         z            =>    col_pp%z                              , & ! Input:  [real(r8) (:,:) ]  soil layer depth (m)
         dz           =>    col_pp%dz                             , & ! Input:  [real(r8) (:,:) ]  layer thickness (m)  (-nlevsno+1:nlevsoi)
         zi           =>    col_pp%zi                             , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m)
         lakedepth    =>    col_pp%lakedepth                      , & ! Input:  [real(r8) (:)   ]  column lake depth (m)
         forc_pbot    =>    top_as%pbot                           , & ! Input:  [real(r8) (:)   ]  atmospheric pressure (Pa)
         t_soisno     =>    col_es%t_soisno         , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)  (-nlevsno+1:nlevsoi)
         lake_icefrac =>    lakestate_vars%lake_icefrac_col       , & ! Input:  [real(r8) (:,:) ]  mass fraction of lake layer that is frozen
         watsat       =>    soilstate_vars%watsat_col             , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)
         h2osoi_vol   =>    col_ws%h2osoi_vol        , & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
         h2osfc       =>    col_ws%h2osfc            , & ! Input:  [real(r8) (:)   ]  surface water (mm)
         frac_h2osfc  =>    col_ws%frac_h2osfc       , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by surface water (0 to 1)
         vgc_max => CH4ParamsInst%vgc_max &
         )

    !$acc enter data create(&
    !$acc ch4_ebul_depth(:,:), &
    !$acc ch4_ebul_total(:), &
    !$acc conc_ch4(:,:), &
    !$acc ch4_aere_depth(:,:), &
    !$acc ch4_oxid_depth(:,:), &
    !$acc pressure)

      vgc_min = vgc_max
      ebul_timescale = dtime ! Allow fast bubbling
      !NOTE: Seems redundant/wasteful to have these local variables 
      !$acc enter data copyin(vgc_min, ebul_timescale) 

      ! column loop to estimate ebullition CH4 flux from each soil layer
      !$acc parallel loop independent gang vector collapse(2) default(present) 
      do j=1,nlevsoi
         do fc = 1, num_methc
            c = filter_methc (fc)
            t = col_pp%topounit(c)

            if (j  >  jwt(c) .and. t_soisno(c,j) > tfrz) then ! Ebullition occurs only below the water table

               k_h_inv = exp(-c_h_inv(1) * (1._r8 / t_soisno(c,j) - 1._r8 / kh_tbase) + log (kh_theta(1))) ! (4.12 Wania) (atm.L/mol)
               k_h = 1._r8 / k_h_inv ! (mol/L.atm)
               k_h_cc = t_soisno(c,j) * k_h * rgasLatm ! (4.21) Wania [(mol/m3w) / (mol/m3g)]

               if (.not. lake) then
                  pressure = forc_pbot(t) + denh2o * grav * (z(c,j)-zi(c,jwt(c))) ! (Pa)
                  if (sat == 1 .and. frac_h2osfc(c) > 0._r8) then ! Add ponding pressure head
                     pressure = pressure + denh2o * grav * h2osfc(c)/1000._r8/frac_h2osfc(c)
                     ! mm     / mm/m
                  end if
               else
                  pressure = forc_pbot(t) + denh2o * grav * (z(c,j) + lakedepth(c))
               end if

               ! Compare partial pressure to ambient pressure.
               vgc = conc_ch4(c,j) / watsat(c,j) / k_h_cc * rgasm * t_soisno(c,j) / pressure
               ! [mol/m3t]      [m3w/m3t]   [m3g/m3w]  [Pa/(mol/m3g)]          [Pa]

               if (vgc > vgc_max * bubble_f) then ! If greater than max value, remove amount down to vgc_min
                  ch4_ebul_depth (c,j) = (vgc - vgc_min * bubble_f) * conc_ch4(c,j) / ebul_timescale
                  ! [mol/m3t/s]                                       [mol/m3t]         [s]
               else
                  ch4_ebul_depth (c,j) = 0._r8
               endif

            else ! above the water table or freezing
               ch4_ebul_depth (c,j) = 0._r8
            endif ! below the water table and not freezing

            ! Prevent ebullition from reaching the surface for frozen lakes
            if (lake .and. lake_icefrac(c,1) > 0.1_r8) ch4_ebul_depth(c,j) = 0._r8

         end do ! fc
      end do ! j


    !$acc exit data delete(&
    !$acc ch4_ebul_depth(:,:), &
    !$acc ch4_ebul_total(:), &
    !$acc conc_ch4(:,:), &
    !$acc ch4_aere_depth(:,:), &
    !$acc ch4_oxid_depth(:,:), &
    !$acc pressure)

    end associate

  end subroutine ch4_ebul

  !-----------------------------------------------------------------------
  subroutine ch4_tran (bounds, &
       num_methc, filter_methc, &
       jwt, dtime_ch4, sat, lake, &
       soilstate_vars, energyflux_vars, ch4_vars, dtime ,&
       o2_decomp_depth, o2stress, ch4_oxid_depth ,& 
       ch4_prod_depth, ch4_aere_depth, ch4_surf_aere, &
       ch4_ebul_depth,ch4_ebul_total,ch4_surf_ebul, &   
       ch4_surf_diff, o2_oxid_depth, o2_aere_depth, & 
       ch4stress, co2_decomp_depth, conc_ch4, conc_o2 &
       )
    !
    ! !DESCRIPTION:
    ! Solves the reaction & diffusion equation for the timestep.  First "competition" between processes for
    ! CH4 & O2 demand is done.  Then concentrations are apportioned into gas & liquid fractions; only the gas
    ! fraction is considered for diffusion in unsat.  Snow and lake water resistance to diffusion is added as
    ! a bulk term in the ground conductance (which is really a surface layer conductance), but concentrations
    ! are not tracked and oxidation is not allowed inside snow and lake water.
    ! Diffusivity is set based on soil texture and organic matter fraction. A Crank-Nicholson solution is used.
    ! Then CH4 diffusive flux is calculated and consistency is checked.

    ! !USES:
    use TridiagonalMod     , only : Tridiagonal_filter
    use CH4varcon          , only : ch4frzout, use_aereoxid_prog
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_methc           ! number of column soil points in column filter
    integer                , intent(in)    :: filter_methc(:)     ! column filter for soil points
    integer                , intent(in)    :: jwt( bounds%begc: ) ! index of the soil layer right above the water table (-) [col]
    integer                , intent(in)    :: sat                 ! 0 = unsaturated; 1 = saturated
    logical                , intent(in)    :: lake      ! function called with lake filter
    real(r8)               , intent(in)    :: dtime_ch4           ! time step for ch4 calculations
    type(soilstate_type)   , intent(in)    :: soilstate_vars
    type(energyflux_type)  , intent(in)    :: energyflux_vars
    type(ch4_type)         , intent(inout) :: ch4_vars
    real(r8) ,intent(in)  :: dtime                               ! land model time step (sec)
    real(r8), intent(inout) :: ch4_prod_depth   (:,:)
    real(r8), intent(inout) :: ch4_oxid_depth   (:,:)
    real(r8), intent(inout) :: ch4_aere_depth   (:,:)
    real(r8), intent(inout) :: ch4_surf_aere    (:)
    real(r8), intent(inout) :: ch4_ebul_depth   (:,:)
    real(r8), intent(inout) :: ch4_ebul_total   (:)
    real(r8), intent(inout) :: ch4_surf_ebul    (:)
    real(r8), intent(inout) :: ch4_surf_diff    (:)
    real(r8), intent(inout) :: o2_oxid_depth    (:,:)
    real(r8), intent(inout) :: o2_decomp_depth  (:,:)
    real(r8), intent(inout) :: o2_aere_depth    (:,:)
    real(r8), intent(inout) :: o2stress         (:,:)
    real(r8), intent(inout) :: ch4stress        (:,:)
    real(r8), intent(inout) :: co2_decomp_depth (:,:)
    real(r8), intent(inout) :: conc_o2          (:,:)
    real(r8), intent(inout) :: conc_ch4         (:,:)
    !
    ! !LOCAL VARIABLES:
    integer :: c,j,g,p,s,i,ll                                  ! indices
    integer :: fc                                              ! soil filter column index
    integer :: fp                                              ! soil filter pft index
    integer,parameter  :: jtop = 0                             ! top level at each column
    integer :: iter                                            ! iteration counter when dtime_ch4 < dtime
    real(r8) :: at (1:num_methc,0:nlevsoi)                     ! "a" vector for tridiagonal matrix
    real(r8) :: bt (1:num_methc,0:nlevsoi)                     ! "b" vector for tridiagonal matrix
    real(r8) :: ct (1:num_methc,0:nlevsoi)                     ! "c" vector for tridiagonal matrix
    real(r8) :: rt (1:num_methc,0:nlevsoi)                     ! "r" vector for tridiagonal solution
    real(r8) :: f_a                                            ! air-filled fraction of available pore space
    real(r8) :: diffus (1:num_methc,0:nlevsoi)                 ! diffusivity (m2/s)
    real(r8) :: k_h_inv                                        ! 1/Henry's Law Constant in Latm/mol
    real(r8) :: k_h_cc(1:num_methc,0:nlevsoi,ngases)           ! ratio of mol/m3 in liquid to mol/m3 in gas
    real(r8) :: dzj                                            !
    real(r8) :: dp1_zp1 (1:num_methc,0:nlevsoi)                ! diffusivity/delta_z for next j
    real(r8) :: dm1_zm1 (1:num_methc,0:nlevsoi)                ! diffusivity/delta_z for previous j
    real(r8) :: t_soisno_c                                     ! soil temperature (C)  (-nlevsno+1:nlevsoi)
    real(r8) :: eps                                            ! either epsilon_a or epsilon_w, depending on where in soil, wrt WT
    real(r8) :: deficit                                        ! mol CH4 /m^2 that must be subtracted from diffusive flux to atm. to make up
    ! for keeping concentrations always above zero
    real(r8) :: conc_ch4_bef(1:num_methc,1:nlevsoi)            ! concentration at the beginning of the timestep
    real(r8) :: errch4(1:num_methc)                            ! Error (Mol CH4 /m^2) [+ = too much CH4]
    real(r8) :: conc_ch4_rel(1:num_methc,0:nlevsoi)            ! Concentration per volume of air or water
    real(r8) :: conc_o2_rel(1:num_methc,0:nlevsoi)             ! Concentration per volume of air or water
    real(r8) :: conc_ch4_rel_old(1:num_methc,0:nlevsoi)        ! Concentration during last Crank-Nich. loop
    real(r8) :: h2osoi_vol_min(1:num_methc,1:nlevsoi)          ! h2osoi_vol restricted to be <= watsat
    real(r8), parameter :: smallnumber = 1.e-12_r8
    real(r8) :: snowdiff                                       ! snow diffusivity (m^2/s)
    real(r8) :: snowres                           ! Cumulative Snow resistance (s/m). Also includes
    real(r8) :: pondres                                        ! Additional resistance from ponding, up to pondmx water on top of top soil layer (s/m)
    real(r8) :: pondz                                          ! Depth of ponding (m)
    real(r8) :: ponddiff                                       ! Pondwater diffusivity (m^2/s)
    real(r8) :: spec_grnd_cond(1:num_methc,1:ngases)           ! species grnd conductance (s/m)
    real(r8) :: airfrac                                                    ! air fraction in snow
    real(r8) :: waterfrac                                                  ! water fraction in snow
    real(r8) :: icefrac                                                    ! ice fraction in snow
    real(r8) :: epsilon_t (1:num_methc,1:nlevsoi,1:ngases)     !
    real(r8) :: epsilon_t_old (1:num_methc,1:nlevsoi,1:ngases) ! epsilon_t from last time step !Currently deprecated
    real(r8) :: source (1:num_methc,1:nlevsoi,1:ngases)        ! source
    real(r8) :: source_old (1:num_methc,1:nlevsoi,1:ngases)    ! source from last time step !Currently deprecated
    real(r8) :: om_frac                                                    ! organic matter fraction
    real(r8) :: o2demand, ch4demand                                        ! mol/m^3/s
    real(r8) :: liqfrac(1:num_methc, 1:nlevsoi)
    real(r8) :: capthick                                                   ! (mm) min thickness before assuming h2osfc is impermeable
    real(r8) :: satpow                                                     ! exponent on watsat for saturated soil solute diffusion
    real(r8) :: scale_factor_gasdiff                                       ! For sensitivity tests; convection would allow this to be > 1
    real(r8) :: scale_factor_liqdiff                                       ! For sensitivity tests; convection would allow this to be > 1
    real(r8) :: organic_max                                                ! organic matter content (kg/m3) where soil is assumed to act like peat
    real(r8) :: aereoxid                                                   ! fraction of methane flux entering aerenchyma rhizosphere

    integer  :: nstep                       ! time step number
    character(len=32) :: subname='ch4_tran' ! subroutine name
    real(r8) :: sum1 
    integer  :: erridx 
    !-----------------------------------------------------------------------
    
    associate(                                                 &
         z             =>    col_pp%z                           , & ! Input:  [real(r8) (:,:) ]  soil layer depth (m)
         dz            =>    col_pp%dz                          , & ! Input:  [real(r8) (:,:) ]  layer thickness (m)  (-nlevsno+1:nlevsoi)
         zi            =>    col_pp%zi                          , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m)
         snl           =>    col_pp%snl                         , & ! Input:  [integer  (:)   ]  negative of number of snow layers

         bsw           =>    soilstate_vars%bsw_col          , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b" (nlevgrnd)
         watsat        =>    soilstate_vars%watsat_col       , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)
         cellorg       =>    soilstate_vars%cellorg_col      , & ! Input:  [real(r8) (:,:) ]  column 3D org (kg/m^3 organic matter) (nlevgrnd)

         t_soisno      =>    col_es%t_soisno   , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)  (-nlevsno+1:nlevsoi)
         t_grnd        =>    col_es%t_grnd     , & ! Input:  [real(r8) (:)   ]  ground temperature (Kelvin)
         t_h2osfc      =>    col_es%t_h2osfc   , & ! Input:  [real(r8) (:)   ]  surface water temperature

         frac_h2osfc   =>    col_ws%frac_h2osfc , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by surface water (0 to 1)
         snow_depth    =>    col_ws%snow_depth  , & ! Input:  [real(r8) (:)   ]  snow height (m)
         h2osoi_vol    =>    col_ws%h2osoi_vol  , & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
         h2osoi_liq    =>    col_ws%h2osoi_liq  , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2) [for snow & soil layers]
         h2osoi_ice    =>    col_ws%h2osoi_ice  , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2) [for snow & soil layers]
         h2osno        =>    col_ws%h2osno      , & ! Input:  [real(r8) (:)   ]  snow water (mm H2O)
         h2osfc        =>    col_ws%h2osfc      , & ! Input:  [real(r8) (:)   ]  surface water (mm)

         c_atm         =>    ch4_vars%c_atm_grc              , & ! Input:  [real(r8) (:,:) ]  CH4, O2, CO2 atmospheric conc  (mol/m3)
         grnd_ch4_cond =>    ch4_vars%grnd_ch4_cond_col      ,  & ! Output: [real(r8) (:)   ]  tracer conductance for boundary layer [m/s]

         ! Set transport parameters
         satpow               => CH4ParamsInst%satpow , &
         scale_factor_gasdiff => CH4ParamsInst%scale_factor_gasdiff,&
         scale_factor_liqdiff => CH4ParamsInst%scale_factor_liqdiff,&
         capthick             => CH4ParamsInst%capthick, &
         aereoxid             => CH4ParamsInst%aereoxid, &

         ! Set shared constant
         organic_max => ParamsShareInst%organic_max &
         )

    !$acc enter data create(&
    !$acc at(:,:), &
    !$acc bt(:,:), &
    !$acc ct(:,:), &
    !$acc rt(:,:), &
    !$acc diffus(:,:), &
    !$acc k_h_cc(:,:,:), &
    !$acc dp1_zp1(:,:), &
    !$acc dm1_zm1(:,:), &
    !$acc conc_ch4_bef(:,:), &
    !$acc errch4(:), &
    !$acc conc_ch4_rel(:,:), &
    !$acc conc_o2_rel(:,:), &
    !$acc conc_ch4_rel_old(:,:), &
    !$acc h2osoi_vol_min(:,:), &
    !$acc snowres, &
    !$acc spec_grnd_cond(:,:), &
    !$acc epsilon_t(:,:,:), &
    !$acc epsilon_t_old(:,:,:), &
    !$acc source(:,:,:), &
    !$acc source_old(:,:,:), &
    !$acc liqfrac(:,:), &
    !$acc snowdiff, &
    !$acc pondres)
   
      ! Perform competition for oxygen and methane in each soil layer if demands over the course of the timestep
      ! exceed that available. Assign to each process in proportion to the quantity demanded in the absense of
      ! the limitation.
      !$acc parallel loop independent gang vector default(present) collapse(2) 
      do j = 1,nlevsoi
         do fc = 1, num_methc
            c = filter_methc (fc)

            o2demand = o2_decomp_depth(c,j) + o2_oxid_depth(c,j) ! o2_decomp_depth includes autotrophic root respiration
            if (o2demand > 0._r8) then
               o2stress(c,j) = min((conc_o2(c,j) / dtime + o2_aere_depth(c,j)) / o2demand, 1._r8)
            else
               o2stress(c,j) = 1._r8
            end if

            ch4demand = ch4_oxid_depth(c,j) + ch4_aere_depth(c,j) + ch4_ebul_depth(c,j)
            if (ch4demand > 0._r8) then
               ch4stress(c,j) = min((conc_ch4(c,j) / dtime + ch4_prod_depth(c,j)) / ch4demand, 1._r8)
            else
               ch4stress(c,j) = 1._r8
            end if

            ! Resolve methane oxidation
            if (o2stress(c,j) < 1._r8 .or. ch4stress(c,j) < 1._r8) then
               if (ch4stress(c,j) <= o2stress(c,j)) then ! methane limited
                  if (o2stress(c,j) < 1._r8) then
                     ! Recalculate oxygen limitation
                     o2demand = o2_decomp_depth(c,j)
                     if (o2demand > 0._r8) then
                        o2stress(c,j) = min( (conc_o2(c,j) / dtime + o2_aere_depth(c,j) - ch4stress(c,j)*o2_oxid_depth(c,j) ) &
                             / o2demand, 1._r8)
                     else
                        o2stress(c,j) = 1._r8
                     end if
                  end if
                  ! Reset oxidation
                  ch4_oxid_depth(c,j) = ch4_oxid_depth(c,j) * ch4stress(c,j)
                  o2_oxid_depth(c,j) = o2_oxid_depth(c,j) * ch4stress(c,j)
               else                                      ! oxygen limited
                  if (ch4stress(c,j) < 1._r8) then
                     ! Recalculate methane limitation
                     ch4demand = ch4_aere_depth(c,j) + ch4_ebul_depth(c,j)
                     if (ch4demand > 0._r8) then
                        ch4stress(c,j) = min( (conc_ch4(c,j) / dtime + ch4_prod_depth(c,j) - &
                             o2stress(c,j)*ch4_oxid_depth(c,j)) / ch4demand, 1._r8)
                     else
                        ch4stress(c,j) = 1._r8
                     end if
                  end if
                  ! Reset oxidation
                  ch4_oxid_depth(c,j) = ch4_oxid_depth(c,j) * o2stress(c,j)
                  o2_oxid_depth(c,j) = o2_oxid_depth(c,j) * o2stress(c,j)
               end if
            end if

            ! Reset non-methanotroph demands
            ch4_aere_depth(c,j) = ch4_aere_depth(c,j) * ch4stress(c,j)
            ch4_ebul_depth(c,j) = ch4_ebul_depth(c,j) * ch4stress(c,j)
            o2_decomp_depth(c,j) = o2_decomp_depth(c,j) * o2stress(c,j)

         end do !c
      end do !j


      ! Accumulate ebullition to place in first layer above water table, or directly to atmosphere
      !$acc parallel loop independent gang worker default(present) private(sum1) 
      do fc = 1, num_methc
         c = filter_methc (fc)
         sum1 = 0._r8 
         !$acc loop vector reduction(+:sum1) 
         do j = 1,nlevsoi
            sum1 = sum1 + ch4_ebul_depth(c,j) * dz(c,j)
         enddo
         ch4_ebul_total(c) = sum1 
      enddo


      ! Set the Henry's Law coefficients
      !$acc parallel loop independent gang vector default(present) collapse(2) 
      do j = 0,nlevsoi
         do fc = 1, num_methc
            c = filter_methc (fc)

            !$acc loop seq
            do s=1,2
               if (j == 0) then
                  k_h_inv = exp(-c_h_inv(s) * (1._r8 / t_grnd(c) - 1._r8 / kh_tbase) + log (kh_theta(s)))
                  ! (4.12) Wania (L atm/mol)
                  k_h_cc(fc,j,s) = t_grnd(c) / k_h_inv * rgasLatm ! (4.21) Wania [(mol/m3w) / (mol/m3g)]
               else
                  k_h_inv = exp(-c_h_inv(s) * (1._r8 / t_soisno(c,j) - 1._r8 / kh_tbase) + log (kh_theta(s)))
                  ! (4.12) Wania (L atm/mol)
                  k_h_cc(fc,j,s) = t_soisno(c,j) / k_h_inv * rgasLatm ! (4.21) Wania [(mol/m3w) / (mol/m3g)]
               end if
            end do
         end do
      end do


      ! Set the source term for each species (no need to do j=0, since epsilon_t and source not used there)
      ! Note that because of the semi-implicit diffusion and the 30 min timestep combined with explicit
      ! sources, occasionally negative concentration will result. In this case it is brought to zero and the
      ! surface flux is adjusted to conserve. This results in some inaccuracy as compared to a shorter timestep
      ! or iterative solution.
      !$acc parallel loop independent gang vector default(present) collapse(2) 
      do j = 1,nlevsoi
         do fc = 1, num_methc
            c = filter_methc (fc)

            if ( .not. use_aereoxid_prog ) then
               ! First remove the CH4 oxidation that occurs at the base of root tissues (aere), and add to oxidation
               ch4_oxid_depth(c,j) = ch4_oxid_depth(c,j) + aereoxid * ch4_aere_depth(c,j)
               ch4_aere_depth(c,j) = ch4_aere_depth(c,j) - aereoxid * ch4_aere_depth(c,j)
            end if ! else oxygen is allowed to diffuse in via aerenchyma

            source(fc,j,1) = ch4_prod_depth(c,j) - ch4_oxid_depth(c,j) - &
                 ch4_aere_depth(c,j) - ch4_ebul_depth(c,j) ! [mol/m3-total/s]
            ! aerenchyma added to surface flux below
            ! ebul added to soil depth just above WT
#ifndef _OPENACC
            if (source(fc,j,1) + conc_ch4(c,j) / dtime < -1.e-12_r8) then
               write(iulog,*) 'Methane demands exceed methane available. Error in methane competition (mol/m^3/s), c,j:',&
                    source(c,j,1) + conc_ch4(c,j) / dtime, c, j
               g = col_pp%gridcell(c)
               write(iulog,*)'Latdeg,Londeg=',grc_pp%latdeg(g),grc_pp%londeg(g)
               call endrun(msg=' ERROR: Methane demands exceed methane available.'&
                    //errMsg(__FILE__, __LINE__))
            else if (ch4stress(c,j) < 1._r8 .and. source(fc,j,1) + conc_ch4(c,j) / dtime > 1.e-12_r8) then
               write(iulog,*) 'Methane limited, yet some left over. Error in methane competition (mol/m^3/s), c,j:',&
                    source(c,j,1) + conc_ch4(c,j) / dtime, c, j
               g = col_pp%gridcell(c)
               write(iulog,*)'Latdeg,Londeg=',grc_pp%latdeg(g),grc_pp%londeg(g)
               call endrun(msg=' ERROR: Methane limited, yet some left over.'//&
                    errMsg(__FILE__, __LINE__))
            end if
#endif
            source(fc,j,2) = -o2_oxid_depth(c,j) - o2_decomp_depth(c,j) + o2_aere_depth(c,j) ! O2 [mol/m3/s]
#ifndef _OPENACC
            if (source(fc,j,2) + conc_o2(c,j) / dtime < -1.e-12_r8) then
               write(iulog,*) 'Oxygen demands exceed oxygen available. Error in oxygen competition (mol/m^3/s), c,j:',&
                    source(c,j,2) + conc_o2(c,j) / dtime, c, j
               g = col_pp%gridcell(c)
               write(iulog,*) 'Latdeg,Londeg=',grc_pp%latdeg(g),grc_pp%londeg(g)
               call endrun(msg=' ERROR: Oxygen demands exceed oxygen available.'//&
                    errMsg(__FILE__, __LINE__) )
            else if (o2stress(c,j) < 1._r8 .and. source(fc,j,2) + conc_o2(c,j) / dtime > 1.e-12_r8) then
               write(iulog,*)  'Oxygen limited, yet some left over. Error in oxygen competition (mol/m^3/s), c,j:',&
                    source(c,j,2) + conc_o2(c,j) / dtime, c, j
               g = col_pp%gridcell(c)
               write(iulog,*)  'Latdeg,Londeg=',grc_pp%latdeg(g),grc_pp%londeg(g)
               call endrun(msg=' ERROR: Oxygen limited, yet some left over.'//errMsg(__FILE__, __LINE__))
            end if
#endif 
            conc_ch4_bef(fc,j) = conc_ch4(c,j) !For Balance Check
         enddo ! fc
      enddo ! j

      ! Accumulate aerenchyma to add directly to atmospheric flux
      !$acc parallel loop independent gang worker default(present) private(sum1) 
      do fc = 1, num_methc
         c = filter_methc (fc)
         sum1 = 0._r8 
         !$acc loop vector reduction(+:sum1)
         do j = 1,nlevsoi
            sum1 = sum1 + ch4_aere_depth(c,j) * dz(c,j)
         enddo
         ch4_surf_aere(c) = sum1 
      enddo

      ! Add in ebullition to source at depth just above WT
      !$acc parallel loop independent gang vector default(present)
      do fc = 1, num_methc
         c = filter_methc(fc)
         if (jwt(c) /= 0) then
            source(fc,jwt(c),1) = source(fc,jwt(c),1) + ch4_ebul_total(c)/dz(c,jwt(c))
         endif
      enddo ! fc

      ! Calculate concentration relative to m^3 of air or water: needed for the diffusion
      !$acc parallel loop independent gang vector default(present) collapse(2) 
      do j = 0,nlevsoi
         do fc = 1, num_methc
            c = filter_methc (fc)
            g = col_pp%gridcell(c)

            if (j == 0) then
               conc_ch4_rel(fc,j) = c_atm(g,1)
               conc_o2_rel(fc,j)  = c_atm(g,2)
            else
               h2osoi_vol_min(fc,j) = min(watsat(c,j), h2osoi_vol(c,j))
               if (ch4frzout) then
                  liqfrac(fc,j) = max(0.05_r8, (h2osoi_liq(c,j)/denh2o+smallnumber)/ &
                       (h2osoi_liq(c,j)/denh2o+h2osoi_ice(c,j)/denice+smallnumber))
               else
                  liqfrac(fc,j) = 1._r8
               end if
               if (j <= jwt(c)) then  ! Above the WT
                  !$acc loop seq
                  do s=1,2
                     epsilon_t(fc,j,s) = watsat(c,j)- (1._r8-k_h_cc(fc,j,s))*h2osoi_vol_min(fc,j)*liqfrac(fc,j)
                  end do
                  ! Partition between the liquid and gas phases. The gas phase will drive the diffusion.
               else ! Below the WT
                  !$acc loop seq
                  do s=1,2
                     epsilon_t(fc,j,s) = watsat(c,j)*liqfrac(fc,j)
                  end do
               end if
               conc_ch4_rel(fc,j) = conc_ch4(c,j)/epsilon_t(fc,j,1)
               conc_o2_rel(fc,j)  = conc_o2(c,j) /epsilon_t(fc,j,2)
            end if
         end do
      end do

      !$acc enter data create(s) 
      ! Loop over species
      do s = 1, 2 ! 1=CH4; 2=O2; 3=CO2

         !$acc update device(s) 

         ! Adjust the grnd_ch4_cond to keep it positive, and add the snow resistance & pond resistance
         !$acc parallel loop independent gang worker default(present) private(snowres)
         do fc = 1, num_methc
            c = filter_methc (fc)
            
            if (grnd_ch4_cond(c) < smallnumber .and. s==1) grnd_ch4_cond(c) = smallnumber
            ! Needed to prevent overflow when ground is frozen, e.g. for lakes
            snowres = 0._r8

            !$acc loop vector reduction(+:snowres)
            do j = -nlevsno + 1,-1

               ! Add snow resistance
               if (j >= snl(c) + 1) then
                  t_soisno_c = t_soisno(c,j) - tfrz
                  icefrac = h2osoi_ice(c,j)/denice/dz(c,j)
                  waterfrac = h2osoi_liq(c,j)/denh2o/dz(c,j)
                  airfrac = max(1._r8 - icefrac - waterfrac, 0._r8)
                  ! Calculate snow diffusivity
                  if (airfrac > 0.05_r8) then
                     f_a = airfrac / (airfrac + waterfrac)
                     eps = airfrac ! Air-filled fraction of total snow volume
                     ! Use Millington-Quirk Expression, as hydraulic properties (bsw) not available
                     snowdiff = (d_con_g(s,1) + d_con_g(s,2)*t_soisno_c) * 1.e-4_r8 * &
                          f_a**(10._r8/3._r8) / (airfrac+waterfrac)**2 * scale_factor_gasdiff
                  else !solute diffusion in water only
                     eps = waterfrac  ! Water-filled fraction of total soil volume
                     snowdiff = eps**satpow * (d_con_w(s,1) + d_con_w(s,2)*t_soisno_c + d_con_w(s,3)*t_soisno_c**2) * 1.e-9_r8 &
                          * scale_factor_liqdiff
                  end if
                  snowdiff = max(snowdiff, smallnumber)
                  snowres = snowres + dz(c,j)/snowdiff
               end if
            end do 

            ! Add pond resistance
            pondres = 0._r8

            ! First old pond formulation up to pondmx
            if (.not. lake .and. snl(c) == 0 .and. h2osoi_vol(c,1) > watsat(c,1)) then
               t_soisno_c = t_soisno(c,1) - tfrz
               if (t_soisno(c,1) <= tfrz) then
                  ponddiff = (d_con_w(s,1) + d_con_w(s,2)*t_soisno_c + d_con_w(s,3)*t_soisno_c**2) * 1.e-9_r8 &
                       * (h2osoi_liq(c,1)/denh2o+smallnumber)/ &
                       (h2osoi_liq(c,1)/denh2o+h2osoi_ice(c,1)/denice+smallnumber) &
                       * scale_factor_liqdiff
               else ! Unfrozen
                  ponddiff = (d_con_w(s,1) + d_con_w(s,2)*t_soisno_c + d_con_w(s,3)*t_soisno_c**2) * 1.e-9_r8 &
                       * scale_factor_liqdiff
               end if
               pondz = dz(c,1) * (h2osoi_vol(c,1) - watsat(c,1))
               pondres = pondz / ponddiff
            end if

            ! Now add new h2osfc form
            if (.not. lake .and. sat == 1 .and. frac_h2osfc(c) > 0._r8 .and. t_h2osfc(c) >= tfrz) then
               t_soisno_c = t_h2osfc(c) - tfrz
               ponddiff = (d_con_w(s,1) + d_con_w(s,2)*t_soisno_c + d_con_w(s,3)*t_soisno_c**2) * 1.e-9_r8 &
                    * scale_factor_liqdiff
               pondz = h2osfc(c) / 1000._r8 / frac_h2osfc(c) ! Assume all h2osfc corresponds to sat area
               ! mm      /  mm/m
               pondres = pondres + pondz / ponddiff
            else if (.not. lake .and. sat == 1 .and. frac_h2osfc(c) > 0._r8 .and. &
                 h2osfc(c)/frac_h2osfc(c) > capthick) then ! Assuming short-circuit logic will avoid FPE here.
               ! assume surface ice is impermeable
               pondres = 1/smallnumber
            end if

            spec_grnd_cond(fc,s) = 1._r8/(1._r8/grnd_ch4_cond(c) + snowres + pondres)

         end do ! fc 

         ! Determine gas diffusion and fraction of open pore (f_a)
         !$acc parallel loop independent gang vector default(present) collapse(2) 
         do j = 1,nlevsoi
            do fc = 1, num_methc
               c = filter_methc (fc)
               g = col_pp%gridcell(c)

               t_soisno_c = t_soisno(c,j) - tfrz

               if (j <= jwt(c)) then  ! Above the WT
                  f_a = 1._r8 - h2osoi_vol_min(fc,j) / watsat(c,j)
                  ! Provisionally calculate diffusivity as linear combination of the Millington-Quirk
                  ! expression in Wania (for peat) & Moldrup (for mineral soil)
                  eps =  watsat(c,j)-h2osoi_vol_min(fc,j) ! Air-filled fraction of total soil volume
                  if (organic_max > 0._r8) then
                     om_frac = min(cellorg(c,j)/organic_max, 1._r8)
                     ! Use first power, not square as in iniTimeConst
                  else
                     om_frac = 1._r8
                  end if
                  diffus(fc,j) = (d_con_g(s,1) + d_con_g(s,2)*t_soisno_c) * 1.e-4_r8 * &
                       (om_frac * f_a**(10._r8/3._r8) / watsat(c,j)**2._r8 + &
                       (1._r8-om_frac) * eps**2._r8 * f_a**(3._r8 / bsw(c,j)) ) &
                       * scale_factor_gasdiff
               else ! Below the WT use saturated diffusivity and only water in epsilon_t
                  ! Note the following is not currently corrected for the effect on diffusivity of excess ice in soil under
                  ! lakes (which is currently experimental only).
                  eps = watsat(c,j)  ! Water-filled fraction of total soil volume
                  diffus(fc,j) = eps**satpow * (d_con_w(s,1) + d_con_w(s,2)*t_soisno_c + d_con_w(s,3)*t_soisno_c**2) * 1.e-9_r8 &
                       * scale_factor_liqdiff
                  if (t_soisno(c,j)<=tfrz) then
                     diffus(fc,j) = diffus(fc,j)*(h2osoi_liq(c,j)/denh2o+smallnumber)/ &
                          (h2osoi_liq(c,j)/denh2o+h2osoi_ice(c,j)/denice+smallnumber)
                  end if
               endif ! Above/below the WT
               diffus(fc,j) = max(diffus(fc,j), smallnumber) ! Prevent overflow

            enddo ! fp
         enddo ! j

         !$acc parallel loop independent gang vector default(present) collapse(2) 
         do j = 1,nlevsoi
            do fc = 1, num_methc
               c = filter_methc (fc)

               ! Set up coefficients for tridiagonal solver.
               if (j == 1 .and. j /= jwt(c) .and. j /= jwt(c)+1) then
                  dm1_zm1(fc,j) = 1._r8/(1._r8/spec_grnd_cond(fc,s)+dz(c,j)/(diffus(fc,j)*2._r8))
                  ! replace Diffusivity / Delta_z by conductance (grnd_ch4_cond) for top layer
                  dp1_zp1(fc,j) = 2._r8/(dz(c,j)/diffus(fc,j)+dz(c,j+1)/diffus(fc,j+1))
               else if (j == 1 .and. j == jwt(c)) then
                  dm1_zm1(fc,j) = 1._r8/(1._r8/spec_grnd_cond(fc,s)+dz(c,j)/(diffus(fc,j)*2._r8))
                  ! layer resistance mult. by k_h_cc for dp1_zp1 term
                  dp1_zp1(fc,j) = 2._r8/(dz(c,j)*k_h_cc(fc,j,s)/diffus(fc,j)+dz(c,j+1)/diffus(fc,j+1))
               else if (j == 1) then ! water table at surface: multiply ground resistance by k_h_cc
                  dm1_zm1(fc,j) = 1._r8/(k_h_cc(fc,j-1,s)/spec_grnd_cond(fc,s)+dz(c,j)/(diffus(fc,j)*2._r8))
                  ! air concentration will be mult. by k_h_cc below
                  dp1_zp1(fc,j) = 2._r8/(dz(c,j)/diffus(fc,j)+dz(c,j+1)/diffus(fc,j+1))
               else if (j <= nlevsoi-1 .and. j /= jwt(c) .and. j /= jwt(c)+1) then
                  dm1_zm1(fc,j) = 2._r8/(dz(c,j)/diffus(fc,j)+dz(c,j-1)/diffus(fc,j-1))
                  dp1_zp1(fc,j) = 2._r8/(dz(c,j)/diffus(fc,j)+dz(c,j+1)/diffus(fc,j+1))
               else if (j <= nlevsoi-1 .and. j == jwt(c)) then ! layer resistance mult. by k_h_cc for dp1_zp1 term
                  dm1_zm1(fc,j) = 2._r8/(dz(c,j)/diffus(fc,j)+dz(c,j-1)/diffus(fc,j-1))
                  dp1_zp1(fc,j) = 2._r8/(dz(c,j)*k_h_cc(fc,j,s)/diffus(fc,j)+dz(c,j+1)/diffus(fc,j+1))
                  ! Concentration in layer will be mult. by k_h_cc below
               else if (j <= nlevsoi-1) then ! j==jwt+1: layer above resistance mult. by k_h_cc for dm1_zm1 term
                  dm1_zm1(fc,j) = 2._r8/(dz(c,j)/diffus(fc,j)+dz(c,j-1)*k_h_cc(fc,j-1,s)/diffus(fc,j-1))
                  ! Concentration in layer above will be mult. by k_h_cc below
                  dp1_zp1(fc,j) = 2._r8/(dz(c,j)/diffus(fc,j)+dz(c,j+1)/diffus(fc,j+1))
               else if (j /= jwt(c)+1) then ! j ==nlevsoi
                  dm1_zm1(fc,j) = 2._r8/(dz(c,j)/diffus(fc,j)+dz(c,j-1)/diffus(fc,j-1))
               else                    ! jwt == nlevsoi-1: layer above resistance mult. by k_h_cc for dm1_zm1 term
                  dm1_zm1(fc,j) = 2._r8/(dz(c,j)/diffus(fc,j)+dz(c,j-1)*k_h_cc(fc,j-1,s)/diffus(fc,j-1))
               end if
            enddo ! fp; pft
         end do ! j; nlevsoi

         ! Perform a second loop for the tridiagonal coefficients since need dp1_zp1 and dm1_z1 at each depth
         !$acc parallel loop independent gang vector default(present) collapse(2) 
         do j = 0,nlevsoi
            do fc = 1, num_methc
               c = filter_methc (fc)
               g = col_pp%gridcell(c)
               conc_ch4_rel_old(fc,j) = conc_ch4_rel(fc,j)

               if (j > 0) dzj = dz(c,j)
               if (j == 0) then ! top layer (atmosphere) doesn't change regardless of where WT is
                  at(fc,j) = 0._r8
                  bt(fc,j) = 1._r8
                  ct(fc,j) = 0._r8
                  rt(fc,j) = c_atm(g,s) ! 0th level stays at constant atmospheric conc
               elseif (j < nlevsoi .and. j == jwt(c)) then ! concentration inside needs to be mult. by k_h_cc for dp1_zp1 term
                  at(fc,j) = -0.5_r8 / dzj * dm1_zm1(fc,j)
                  bt(fc,j) = epsilon_t(fc,j,s) / dtime_ch4 + 0.5_r8 / dzj * (dp1_zp1(fc,j)*k_h_cc(fc,j,s) + dm1_zm1(fc,j))
                  ct(fc,j) = -0.5_r8 / dzj * dp1_zp1(fc,j)
               elseif (j < nlevsoi .and. j == jwt(c)+1) then
                  ! concentration above needs to be mult. by k_h_cc for dm1_zm1 term
                  at(fc,j) = -0.5_r8 / dzj * dm1_zm1(fc,j) * k_h_cc(fc,j-1,s)
                  bt(fc,j) = epsilon_t(fc,j,s) / dtime_ch4 + 0.5_r8 / dzj * (dp1_zp1(fc,j) + dm1_zm1(fc,j))
                  ct(fc,j) = -0.5_r8 / dzj * dp1_zp1(fc,j)
               elseif (j < nlevsoi) then
                  at(fc,j) = -0.5_r8 / dzj * dm1_zm1(fc,j)
                  bt(fc,j) = epsilon_t(fc,j,s) / dtime_ch4 + 0.5_r8 / dzj * (dp1_zp1(fc,j) + dm1_zm1(fc,j))
                  ct(fc,j) = -0.5_r8 / dzj * dp1_zp1(fc,j)
               else if (j == nlevsoi .and. j== jwt(c)+1) then
                  ! concentration above needs to be mult. by k_h_cc for dm1_zm1 term
                  at(fc,j) = -0.5_r8 / dzj * dm1_zm1(fc,j) * k_h_cc(fc,j-1,s)
                  bt(fc,j) = epsilon_t(fc,j,s) / dtime_ch4 + 0.5_r8 / dzj * dm1_zm1(fc,j)
                  ct(fc,j) = 0._r8
               else ! j==nlevsoi and jwt<nlevsoi-1 or jwt==nlevsoi: 0 flux at bottom
                  at(fc,j) = -0.5_r8 / dzj * dm1_zm1(fc,j)
                  bt(fc,j) = epsilon_t(fc,j,s) / dtime_ch4 + 0.5_r8 / dzj * dm1_zm1(fc,j)
                  ct(fc,j) = 0._r8
               endif
            enddo ! fp; pft
         enddo ! j; nlevsoi

         if (s == 1) then  ! CH4

            ! Set rt, since it depends on conc
            !$acc parallel loop independent gang vector default(present) collapse(2) 
            do j = 1,nlevsoi
               do fc = 1, num_methc
                  c = filter_methc (fc)

                  ! For correct balance, deprecate source_old.
                  ! The source terms are effectively constant over the timestep.
                  source_old(fc,j,s) = source(fc,j,s)
                  ! source_old could be removed later
                  epsilon_t_old(fc,j,s) = epsilon_t(fc,j,s)
                  ! epsilon_t acts like source also
                  dzj = dz(c,j)
                  if (j < nlevsoi .and. j == jwt(c)) then ! concentration inside needs to be mult. by k_h_cc for dp1_zp1 term
                     rt(fc,j) = epsilon_t_old(fc,j,s) / dtime_ch4 * conc_ch4_rel(fc,j) +           &
                          0.5_r8 / dzj * (dp1_zp1(fc,j) * (conc_ch4_rel(fc,j+1)-conc_ch4_rel(fc,j)*k_h_cc(fc,j,s)) - &
                          dm1_zm1(fc,j) * (conc_ch4_rel(fc,j)  -conc_ch4_rel(fc,j-1))) + &
                          0.5_r8 * (source(fc,j,s) + source_old(fc,j,s))
                  elseif (j < nlevsoi .and. j == jwt(c)+1) then
                     ! concentration above needs to be mult. by k_h_cc for dm1_zm1 term
                     rt(fc,j) = epsilon_t_old(fc,j,s) / dtime_ch4 * conc_ch4_rel(fc,j) +           &
                          0.5_r8 / dzj * (dp1_zp1(fc,j) * (conc_ch4_rel(fc,j+1)-conc_ch4_rel(fc,j)) - &
                          dm1_zm1(fc,j) * (conc_ch4_rel(fc,j) -conc_ch4_rel(fc,j-1)*k_h_cc(fc,j-1,s))) + &
                          0.5_r8 * (source(fc,j,s) + source_old(fc,j,s))
                  elseif (j < nlevsoi) then
                     rt(fc,j) = epsilon_t_old(fc,j,s) / dtime_ch4 * conc_ch4_rel(fc,j) +           &
                          0.5_r8 / dzj * (dp1_zp1(fc,j) * (conc_ch4_rel(fc,j+1)-conc_ch4_rel(fc,j)) - &
                          dm1_zm1(fc,j) * (conc_ch4_rel(fc,j)  -conc_ch4_rel(fc,j-1))) + &
                          0.5_r8 * (source(fc,j,s) + source_old(fc,j,s))
                  else if (j == nlevsoi .and. j== jwt(c)+1) then
                     ! concentration above needs to be mult. by k_h_cc for dm1_zm1 term
                     rt(fc,j) = epsilon_t_old(fc,j,s) / dtime_ch4 * conc_ch4_rel(fc,j) +           &
                          0.5_r8 / dzj * ( - dm1_zm1(fc,j) * (conc_ch4_rel(fc,j) -conc_ch4_rel(fc,j-1)*k_h_cc(fc,j-1,s))) + &
                          0.5_r8 * (source(fc,j,s) + source_old(fc,j,s))
                  else  !j==nlevsoi
                     rt(fc,j) = epsilon_t_old(fc,j,s) / dtime_ch4 * conc_ch4_rel(fc,j) +           &
                          0.5_r8 / dzj * ( - dm1_zm1(fc,j) * (conc_ch4_rel(fc,j)  -conc_ch4_rel(fc,j-1))) + &
                          0.5_r8 * (source(fc,j,s) + source_old(fc,j,s))
                  endif
                  epsilon_t_old(fc,j,s) = epsilon_t(fc,j,s)
                  source_old(fc,j,s) = source(fc,j,s)

               enddo ! fc; column
            enddo ! j; nlevsoi

            call Tridiagonal_filter( 0, nlevsoi, jtop, num_methc, filter_methc, &
                 at, &
                 bt, &
                 ct, &
                 rt, &
                 conc_ch4_rel)

            ! Calculate net ch4 flux to the atmosphere from the surface (+ to atm)
            !$acc parallel loop independent gang vector default(present)
            do fc = 1, num_methc
               c = filter_methc (fc)
               g = col_pp%gridcell(c)
               if (jwt(c) /= 0) then ! WT not at the surface
                  ch4_surf_diff(c) = dm1_zm1(fc,1) * ( (conc_ch4_rel(fc,1)+conc_ch4_rel_old(fc,1))/2._r8 &
                       - c_atm(g,s)) ! [mol/m2/s]
                  ch4_surf_ebul(c) = 0._r8 ! all the ebullition has already come out in the soil column (added to source)
                  ! Try adding directly to atm. to prevent destabilization of diffusion
                  !ch4_surf_ebul(c) = ch4_ebul_total(c) ! [mol/m2/s]
               else ! WT at the surface; i.e., jwt(c)==0
                  ch4_surf_diff(c) = dm1_zm1(fc,1) * ( (conc_ch4_rel(fc,1)+conc_ch4_rel_old(fc,1))/2._r8 &
                       - c_atm(g,s)*k_h_cc(fc,0,s)) ! [mol/m2/s]
                  ! atmospheric concentration gets mult. by k_h_cc as above
                  ch4_surf_ebul(c) = ch4_ebul_total(c) ! [mol/m2/s]
               endif
            enddo

            ! Ensure that concentrations stay above 0
            ! This should be done after the flux, so that the flux calculation is consistent.
            !$acc parallel loop independent gang vector collapse(2) default(present)
            do j = 1,nlevsoi
               do fc = 1, num_methc
                  c = filter_methc (fc)

                  if (conc_ch4_rel(fc,j) < 0._r8) then
                     deficit = - conc_ch4_rel(fc,j)*epsilon_t(fc,j,1)*dz(c,j)  ! Mol/m^2 added
                     #ifndef _OPENACC 
                     if (deficit > 1.e-3_r8 * scale_factor_gasdiff) then
                        if (deficit > 1.e-2_r8) then
                           write(iulog,*)  'Note: sink > source in ch4_tran, sources are changing '// &
                                ' quickly relative to diffusion timestep, and/or diffusion is rapid.'
                           g = col_pp%gridcell(c)
                           write(iulog,*) 'Latdeg,Londeg=',grc_pp%latdeg(g),grc_pp%londeg(g)
                           write(iulog,*) 'This typically occurs when there is a larger than normal diffusive flux '
                           write(iulog,*) 'If this occurs frequently, consider reducing land model (or  methane model) timestep, or reducing the max. sink per timestep in the methane model.'
                        end if
                     end if
                     #endif 
                     conc_ch4_rel(fc,j) = 0._r8
                     ! Subtract deficit
                     ch4_surf_diff(c) = ch4_surf_diff(c) - deficit/dtime_mod
                  end if
               enddo
            enddo


         elseif (s == 2) then  ! O2

            ! Set rt, since it depends on conc
            !$acc parallel loop independent collapse(2) gang vector default(present) 
            do j = 1,nlevsoi
               do fc = 1, num_methc
                  c = filter_methc (fc)

                  ! For correct balance, deprecate source_old.
                  source_old(fc,j,s) = source(fc,j,s)
                  ! source_old could be removed later
                  epsilon_t_old(fc,j,s) = epsilon_t(fc,j,s)
                  ! epsilon_t acts like source also
                  dzj     = dz(c,j)
                  if (j < nlevsoi .and. j == jwt(c)) then ! concentration inside needs to be mult. by k_h_cc for dp1_zp1 term
                     rt(fc,j) = epsilon_t_old(fc,j,s) / dtime_ch4 * conc_o2_rel(fc,j) +           &
                          0.5_r8 / dzj * (dp1_zp1(fc,j) * (conc_o2_rel(fc,j+1)-conc_o2_rel(fc,j)*k_h_cc(fc,j,s)) - &
                          dm1_zm1(fc,j) * (conc_o2_rel(fc,j)  -conc_o2_rel(fc,j-1))) + &
                          0.5_r8 * (source(fc,j,s) + source_old(fc,j,s))
                  elseif (j < nlevsoi .and. j == jwt(c)+1) then
                     ! concentration above needs to be mult. by k_h_cc for dm1_zm1 term
                     rt(fc,j) = epsilon_t_old(fc,j,s) / dtime_ch4 * conc_o2_rel(fc,j) +           &
                          0.5_r8 / dzj * (dp1_zp1(fc,j) * (conc_o2_rel(fc,j+1)-conc_o2_rel(fc,j)) - &
                          dm1_zm1(fc,j) * (conc_o2_rel(fc,j) -conc_o2_rel(fc,j-1)*k_h_cc(fc,j-1,s))) + &
                          0.5_r8 * (source(fc,j,s) + source_old(fc,j,s))
                  elseif (j < nlevsoi) then
                     rt(fc,j) = epsilon_t_old(fc,j,s) / dtime_ch4 * conc_o2_rel(fc,j) +           &
                          0.5_r8 / dzj * (dp1_zp1(fc,j) * (conc_o2_rel(fc,j+1)-conc_o2_rel(fc,j)) - &
                          dm1_zm1(fc,j) * (conc_o2_rel(fc,j)  -conc_o2_rel(fc,j-1))) + &
                          0.5_r8 * (source(fc,j,s) + source_old(fc,j,s))
                  else if (j == nlevsoi .and. j== jwt(c)+1) then
                     ! concentration above needs to be mult. by k_h_cc for dm1_zm1 term
                     rt(fc,j) = epsilon_t_old(fc,j,s) / dtime_ch4 * conc_o2_rel(fc,j) +           &
                          0.5_r8 / dzj * ( - dm1_zm1(fc,j) * (conc_o2_rel(fc,j) -conc_o2_rel(fc,j-1)*k_h_cc(fc,j-1,s))) + &
                          0.5_r8 * (source(fc,j,s) + source_old(fc,j,s))
                  else  !j==nlevsoi
                     rt(fc,j) = epsilon_t_old(fc,j,s) / dtime_ch4 * conc_o2_rel(fc,j) +           &
                          0.5_r8 / dzj * ( - dm1_zm1(fc,j) * (conc_o2_rel(fc,j)  -conc_o2_rel(fc,j-1))) + &
                          0.5_r8 * (source(fc,j,s) + source_old(fc,j,s))
                  endif
                  epsilon_t_old(fc,j,s) = epsilon_t(fc,j,s)
                  source_old(fc,j,s) = source(fc,j,s)

               enddo ! fc; column
            enddo ! j; nlevsoi

            call Tridiagonal_filter( 0, nlevsoi, jtop, &
                 num_methc, filter_methc, &
                 at, &
                 bt, &
                 ct, &
                 rt, &
                 conc_o2_rel)

            ! Ensure that concentrations stay above 0
            !$acc parallel loop independent gang vector default(present) collapse(2) 
            do j = 1,nlevsoi
               do fc = 1, num_methc
                  c = filter_methc (fc)
                  g = col_pp%gridcell(c)
                  conc_o2_rel(fc,j) = max (conc_o2_rel(fc,j), 1.e-12_r8)
                  ! In case of pathologically large aerenchyma conductance. Should be OK in general but
                  ! this will maintain stability even if a PFT with very small weight somehow has an absurd NPP or LAI.
                  ! Also, oxygen above ambient will probably bubble.
                  conc_o2_rel(fc,j) = min (conc_o2_rel(fc,j), c_atm(g,2)/epsilon_t(fc,j,2))
               enddo
            enddo

         endif  ! species

      enddo  ! species
      !$acc exit data delete(s) 

      ! Update absolute concentrations per unit volume
      !$acc parallel loop independent gang vector default(present) collapse(2) 
      do j = 1,nlevsoi ! No need to update the atm. level concentrations
         do fc = 1, num_methc
            c = filter_methc (fc)
            conc_ch4(c,j) = conc_ch4_rel(fc,j)*epsilon_t(fc,j,1)
            conc_o2(c,j)  = conc_o2_rel(fc,j) *epsilon_t(fc,j,2)
         end do
      end do

      ! Do Balance Check and absorb small
      !    discrepancy into surface flux.
      !$acc parallel loop independent gang worker default(present) private(sum1)
      do fc = 1, num_methc
         c = filter_methc (fc)
         sum1 = 0._r8
         !$acc loop vector reduction(+:sum1)
         do j = 1,nlevsoi
            sum1 = sum1 + (conc_ch4(c,j) - conc_ch4_bef(fc,j))*dz(c,j)
            sum1 = sum1 - ch4_prod_depth(c,j)*dz(c,j)*dtime_mod
            sum1 = sum1 + ch4_oxid_depth(c,j)*dz(c,j)*dtime_mod
         end do
         errch4(fc) = sum1 
      end do

      erridx = 0
      !$acc parallel loop independent gang vector copy(erridx)
      do fc = 1, num_methc
         c = filter_methc (fc)

         ! For history make sure that grnd_ch4_cond includes snow, for methane diffusivity
         grnd_ch4_cond(c) = spec_grnd_cond(fc,1)
         errch4(fc) = errch4(fc) + (ch4_surf_aere(c) + ch4_surf_ebul(c) + ch4_surf_diff(c))*dtime

         if (abs(errch4(fc)) < 1.e-8_r8) then
            ch4_surf_diff(c) = ch4_surf_diff(c) - errch4(fc)/dtime
         else ! errch4 > 1e-8 mol / m^2 / timestep
            erridx = c
         end if
      end do

      if(erridx > 0 ) then 
         g = col_pp%gridcell(erridx)
         write(iulog,*) 'CH4 Conservation Error in CH4Mod during diffusion', nstep, erridx, '(mol /m^2.timestep)'
         write(iulog,*) 'Latdeg,Londeg=',grc_pp%latdeg(g),grc_pp%londeg(g)
      end if 

    !$acc exit data delete(&
    !$acc at(:,:), &
    !$acc bt(:,:), &
    !$acc ct(:,:), &
    !$acc rt(:,:), &
    !$acc diffus(:,:), &
    !$acc k_h_cc(:,:,:), &
    !$acc dp1_zp1(:,:), &
    !$acc dm1_zm1(:,:), &
    !$acc conc_ch4_bef(:,:), &
    !$acc errch4(:), &
    !$acc conc_ch4_rel(:,:), &
    !$acc conc_o2_rel(:,:), &
    !$acc conc_ch4_rel_old(:,:), &
    !$acc h2osoi_vol_min(:,:), &
    !$acc snowres, &
    !$acc spec_grnd_cond(:,:), &
    !$acc epsilon_t(:,:,:), &
    !$acc epsilon_t_old(:,:,:), &
    !$acc source(:,:,:), &
    !$acc source_old(:,:,:), &
    !$acc liqfrac(:,:), &
    !$acc snowdiff, &
    !$acc pondres)

    end associate

  end subroutine ch4_tran

  !-----------------------------------------------------------------------
  subroutine get_jwt (bounds, num_methc, filter_methc, jwt, &
       soilstate_vars )
    !
    ! !DESCRIPTION:
    ! Finds the first unsaturated layer going up. Also allows a perched water table over ice.
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)  :: bounds
    integer                , intent(in)  :: num_methc           ! number of column soil points in column filter
    integer                , intent(in)  :: filter_methc(:)     ! column filter for soil points
    integer                , intent(out) :: jwt( bounds%begc: ) ! index of the soil layer right above the water table (-) [col]
    type(soilstate_type)   , intent(in)  :: soilstate_vars
    !
    ! !LOCAL VARIABLES:
    real(r8) :: f_sat    ! volumetric soil water defining top of water table or where production is allowed
    integer  :: c,j,perch! indices
    integer  :: fc       ! filter column index
    !-----------------------------------------------------------------------


    associate(                                          &
         watsat     => soilstate_vars%watsat_col      , & ! Input:  [real(r8) (:,:)  ] volumetric soil water at saturation (porosity)
         h2osoi_vol => col_ws%h2osoi_vol , & ! Input:  [real(r8) (:,:)  ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
         t_soisno   => col_es%t_soisno    & ! Input:  [real(r8) (: ,:) ]  soil temperature (Kelvin)  (-nlevsno+1:nlevsoi)
         )

      f_sat = CH4ParamsInst%f_sat

      ! The layer index of the first unsaturated layer, i.e., the layer right above
      ! the water table.
      ! ZS: Loop is currently not vectorized.
      !$acc parallel loop independent gang vector default(present)
      do fc = 1, num_methc
         c = filter_methc(fc)

         ! Check to see if any soil layers are frozen and saturated.  If so, start looking at the first layer above the top
         ! such layer.  This is potentially important for perched water tables in the Tundra.

         perch = nlevsoi
         !$acc loop seq
         do j = nlevsoi, 1, -1
            if (t_soisno(c,j) < tfrz .and. h2osoi_vol(c,j) > f_sat * watsat(c,j)) then
               ! strictly less than freezing because it could be permeable otherwise
               perch = j-1
            end if
         end do
         jwt(c) = perch

         !$acc loop seq
         do j = perch, 2, -1
            if(h2osoi_vol(c,j) > f_sat * watsat(c,j) .and. h2osoi_vol(c,j-1) < f_sat * watsat(c,j-1)) then
               jwt(c) = j-1
               exit
            end if
         enddo
         if (jwt(c) == perch .and. h2osoi_vol(c,1) > f_sat * watsat(c,1)) then ! missed that the top layer is saturated
            jwt(c) = 0
         endif
      end do

    end associate

  end subroutine get_jwt

  !-----------------------------------------------------------------------
  subroutine ch4_annualupdate(bounds, &
       num_methc, filter_methc, &
       num_methp, filter_methp, ch4_vars, dt, dayspyr)
    !
    ! !DESCRIPTION: Annual mean fields.
    !
    ! !USES:
    use elm_varcon      , only: secspday
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_methc         ! number of soil columns in filter
    integer                , intent(in)    :: filter_methc(:)   ! filter for soil columns
    integer                , intent(in)    :: num_methp         ! number of soil points in pft filter
    integer                , intent(in)    :: filter_methp(:)   ! patch filter for soil points
    type(ch4_type)         , intent(inout) :: ch4_vars
    real(r8),  intent(in)   :: dt        ! time step (seconds)
    real(r8), intent(in) :: dayspyr


    !
    ! !LOCAL VARIABLES:
    integer :: c,p       ! indices
    integer :: fc        ! soil column filter indices
    integer :: fp        ! soil pft filter indices
    real(r8):: secsperyear
    logical :: newrun
    !-----------------------------------------------------------------------

    associate(                                                      &
         agnpp          =>    veg_cf%agnpp         , & ! Input:  [real(r8) (:) ]  (gC/m2/s) aboveground NPP
         bgnpp          =>    veg_cf%bgnpp         , & ! Input:  [real(r8) (:) ]  (gC/m2/s) belowground NPP
         somhr          =>    col_cf%somhr           , & ! Input:  [real(r8) (:) ]  (gC/m2/s) soil organic matter heterotrophic respiration
         tempavg_agnpp  =>    veg_cf%tempavg_agnpp , & ! Output: [real(r8) (:) ]  temporary average above-ground NPP (gC/m2/s)
         annavg_agnpp   =>    veg_cf%annavg_agnpp  , & ! Output: [real(r8) (:) ]  annual average above-ground NPP (gC/m2/s)
         tempavg_bgnpp  =>    veg_cf%tempavg_bgnpp , & ! Output: [real(r8) (:) ]  temporary average below-ground NPP (gC/m2/s)
         annavg_bgnpp   =>    veg_cf%annavg_bgnpp  , & ! Output: [real(r8) (:) ]  annual average below-ground NPP (gC/m2/s)

         finundated     =>    ch4_vars%finundated_col             , & ! Input:  [real(r8) (:) ]  fractional inundated area in soil column
         annsum_counter =>    ch4_vars%annsum_counter_col         , & ! Output: [real(r8) (:) ]  seconds since last annual accumulator turnover
         tempavg_somhr  =>    ch4_vars%tempavg_somhr_col          , & ! Output: [real(r8) (:) ]  temporary average SOM heterotrophic resp. (gC/m2/s)
         annavg_somhr   =>    ch4_vars%annavg_somhr_col           , & ! Output: [real(r8) (:) ]  annual average SOM heterotrophic resp. (gC/m2/s)
         tempavg_finrw  =>    ch4_vars%tempavg_finrw_col          , & ! Output: [real(r8) (:) ]  respiration-weighted annual average of finundated
         annavg_finrw   =>    ch4_vars%annavg_finrw_col             & ! Output: [real(r8) (:) ]  respiration-weighted annual average of finundated
         )

      ! set time steps
      secsperyear = real( dayspyr * secspday, r8)

      newrun = .false.

      ! column loop
      !$acc parallel loop independent gang vector default(present)
      do fc = 1,num_methc
         c = filter_methc(fc)

         if (annsum_counter(c) == spval) then
            ! These variables are now in restart files for completeness, but might not be in inicFile and are not.
            ! set for arbinit.
            newrun = .true.
            annsum_counter(c)    = 0._r8
            tempavg_somhr(c)     = 0._r8
            tempavg_finrw(c)     = 0._r8
         end if

         annsum_counter(c) = annsum_counter(c) + dt
      end do

      ! patch loop
      !$acc parallel loop independent gang vector default(present)
      do fp = 1,num_methp
         p = filter_methp(fp)

         if (newrun .or. tempavg_agnpp(p) == spval) then ! Extra check needed because for back-compatibility
            tempavg_agnpp(p) = 0._r8
            tempavg_bgnpp(p) = 0._r8
         end if
      end do

      !$acc parallel loop independent gang vector default(present)
      do fc = 1,num_methc
         c = filter_methc(fc)
         if (annsum_counter(c) >= secsperyear) then

            ! update annual average somhr
            annavg_somhr(c)      =  tempavg_somhr(c)
            tempavg_somhr(c)     = 0._r8

            ! update annual average finrw
            if (annavg_somhr(c) > 0._r8) then
               annavg_finrw(c)      =  tempavg_finrw(c) / annavg_somhr(c)
            else
               annavg_finrw(c)      = 0._r8
            end if
            tempavg_finrw(c)     = 0._r8
         else
            tempavg_somhr(c)     = tempavg_somhr(c) + dt/secsperyear * somhr(c)
            tempavg_finrw(c)     = tempavg_finrw(c) + dt/secsperyear * finundated(c) * somhr(c)
         end if
      end do

      !$acc parallel loop independent gang vector default(present)
      do fp = 1,num_methp
         p = filter_methp(fp)
         c = veg_pp%column(p)
         if (annsum_counter(c) >= secsperyear) then

            annavg_agnpp(p) = tempavg_agnpp(p)
            tempavg_agnpp(p) = 0._r8

            annavg_bgnpp(p) = tempavg_bgnpp(p)
            tempavg_bgnpp(p) = 0._r8

         else
            tempavg_agnpp(p) = tempavg_agnpp(p) + dt/secsperyear * agnpp(p)
            tempavg_bgnpp(p) = tempavg_bgnpp(p) + dt/secsperyear * bgnpp(p)
         end if
      end do

      ! column loop
      !$acc parallel loop independent gang vector default(present)
      do fc = 1,num_methc
         c = filter_methc(fc)
         if (annsum_counter(c) >= secsperyear) annsum_counter(c) = 0._r8
      end do

    end associate

  end subroutine ch4_annualupdate

end module CH4Mod

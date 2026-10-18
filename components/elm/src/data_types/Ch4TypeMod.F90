module Ch4TypeMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Methane (CH4) state and flux data type, moved out of CH4Mod.
  !
  ! !USES:
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_infnan_mod    , only : nan => shr_infnan_nan
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use elm_varpar        , only : nlevsoi, ngases, nlevdecomp
  use elm_varcon        , only : spval, grlnd
  use landunit_varcon   , only : istdlak
  use elm_varctl        , only : iulog, use_lch4
  use abortutils        , only : endrun
  use decompMod         , only : bounds_type
  use GridcellType      , only : grc_pp
  use ColumnType        , only : col_pp
  use topounit_varcon   , only : max_topounits
  !
  implicit none
  save
  private

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
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init( this, bounds, cellorg_col )

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
    use shr_infnan_mod, only: nan  => shr_infnan_nan, assignment(=)
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
    integer               :: j ,g,l,c,p,t,ti,topi ! indices
    type(file_desc_t)     :: ncid        ! netcdf id
    real(r8)     ,pointer :: zwt0_in (:,:) ! read in - zwt0 
    real(r8)     ,pointer :: f0_in (:,:)   ! read in - f0 
    real(r8)     ,pointer :: p3_in (:,:)   ! read in - p3 
    real(r8)     ,pointer :: pH_in (:,:)   ! read in - pH 
    logical               :: readvar 
    character(len=256)    :: locfn
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(cellorg_col) == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))

    !----------------------------------------
    ! Initialize time constant variables
    !----------------------------------------

    allocate(zwt0_in(bounds%begg:bounds%endg,1:max_topounits))
    allocate(f0_in(bounds%begg:bounds%endg,1:max_topounits))
    allocate(p3_in(bounds%begg:bounds%endg,1:max_topounits))
    if (usephfact) allocate(ph_in(bounds%begg:bounds%endg,1:max_topounits))

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
       t = col_pp%topounit(c)
       topi = grc_pp%topi(g)
       ti = t - topi + 1

       if (.not. fin_use_fsat) then
          this%zwt0_col(c) = zwt0_in(g,ti)
          this%f0_col(c)   = f0_in(g,ti)
          this%p3_col(c)   = p3_in(g,ti)
       end if
       if (usephfact) this%pH_col(c) = pH_in(g,ti)
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
       if (col_pp%is_soil(c) .or. col_pp%is_crop(c)) then

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

       else if (col_pp%is_lake(c)) then

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
       this%ch4_ebul_depth_sat_col     (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%ch4_ebul_depth_unsat_col   (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%conc_ch4_lake_col          (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%conc_o2_lake_col           (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%ch4stress_unsat_col        (c,nlevsoi+1:nlevgrnd) = 0._r8
       this%ch4stress_sat_col          (c,nlevsoi+1:nlevgrnd) = 0._r8

       if (col_pp%is_soil(c) .or. col_pp%is_crop(c)) then

          this%conc_ch4_lake_col       (c,:) = spval
          this%conc_o2_lake_col        (c,:) = spval
          this%ch4_surf_diff_lake_col  (c)   = spval
          this%ch4_surf_ebul_lake_col  (c)   = spval
          this%ch4_prod_depth_lake_col (c,:) = spval
          this%ch4_oxid_depth_lake_col (c,:) = spval

       else if (col_pp%is_lake(c) .and. allowlakeprod) then

          this%ch4_prod_depth_unsat_col   (c,:) = spval
          this%ch4_oxid_depth_unsat_col   (c,:) = spval
          this%o2_oxid_depth_unsat_col    (c,:) = spval
          this%o2_decomp_depth_unsat_col  (c,:) = spval
          this%o2_aere_depth_unsat_col    (c,:) = spval
          this%co2_decomp_depth_unsat_col (c,:) = spval
          this%co2_oxid_depth_unsat_col   (c,:) = spval
          this%ch4_aere_depth_unsat_col   (c,:) = spval
          this%ch4_tran_depth_unsat_col   (c,:) = spval
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

end module Ch4TypeMod

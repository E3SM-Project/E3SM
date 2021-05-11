module ColumnDataType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Column data type allocation and initialization
  ! -------------------------------------------------------- 
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use clm_varpar      , only : ndecomp_cascade_transitions, ndecomp_pools, nlevcan
  use clm_varpar     , only : nlevsno, nlevgrnd, nlevlak
  use clm_varpar     , only : nlevdecomp_full, crop_prog, nlevdecomp
  use clm_varcon     , only : spval, ispval

  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private

  !
  ! NOTE(bandre, 2013-10) according to Charlie Koven, nfix_timeconst
  ! is currently used as a flag and rate constant. Rate constant: time
  ! over which to exponentially relax the npp flux for N fixation term
  ! flag: (if  <=  0. or  >=  365; use old annual method). Default value is
  ! junk that should always be overwritten by the namelist or init function!
  !
  ! (days) time over which to exponentially relax the npp flux for N fixation term
  real(r8), public :: nfix_timeconst = -1.2345_r8 
  !
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds energy state information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_energy_state
    ! temperature variables
    real(r8), pointer :: t_soisno      (:,:) => null() ! soil temperature (K)  (-nlevsno+1:nlevgrnd) 
    real(r8), pointer :: t_ssbef       (:,:) => null() ! soil/snow temperature before update (K) (-nlevsno+1:nlevgrnd) 
    real(r8), pointer :: t_h2osfc      (:)   => null() ! surface water temperature (K)
    real(r8), pointer :: t_h2osfc_bef  (:)   => null() ! surface water temperature at start of time step (K)
    real(r8), pointer :: t_soi10cm     (:)   => null() ! soil temperature in top 10cm of soil (K)
    real(r8), pointer :: t_soi17cm     (:)   => null() ! soil temperature in top 17cm of soil (K)
    real(r8), pointer :: t_grnd        (:)   => null() ! ground temperature (K)
    real(r8), pointer :: t_lake        (:,:) => null() ! lake temperature (K)  (1:nlevlak)          
    real(r8), pointer :: t_grnd_r      (:)   => null() ! rural ground temperature (K)
    real(r8), pointer :: t_grnd_u      (:)   => null() ! urban ground temperature (K) (needed by Hydrology2Mod)
    real(r8), pointer :: snot_top      (:)   => null() ! temperature of top snow layer (K)
    real(r8), pointer :: dTdz_top      (:)   => null() ! temperature gradient in top layer  (K m-1)
    real(r8), pointer :: thv           (:)   => null() ! virtual potential temperature (K)
    ! heat content variables (diagnostic)
    real(r8), pointer :: hc_soi        (:)   => null() ! soil heat content (MJ/m2)
    real(r8), pointer :: hc_soisno     (:)   => null() ! soil plus snow heat content (MJ/m2)
    ! other quantities related to energy state
    real(r8), pointer :: emg           (:)   => null() ! ground emissivity (unitless)
    real(r8), pointer :: fact          (:,:) => null() ! factors used in computing tridiagonal matrix
    real(r8), pointer :: c_h2osfc      (:)   => null() ! heat capacity of surface water (J/K)
    ! For coupling with pflotran TH
    real(r8), pointer :: t_nearsurf    (:)   => null() ! near-surface air temperature averaged over bare-veg (K)
    
  contains
    procedure, public :: Init    => col_es_init
    procedure, public :: Clean   => col_es_clean
  end type column_energy_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds water state information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_water_state
    ! Note: All units given as kg in this data type imply kg of H2O
    ! Primary water and ice state variables for soil/snow column
    real(r8), pointer :: h2osoi_liq         (:,:) => null() ! liquid water (-nlevsno+1:nlevgrnd) (kg/m2)     
    real(r8), pointer :: h2osoi_ice         (:,:) => null() ! ice lens (-nlevsno+1:nlevgrnd) (kg/m2)    
    real(r8), pointer :: h2osoi_vol         (:,:) => null() ! volumetric soil water (0<=h2osoi_vol<=watsat) (1:nlevgrnd) (m3/m3)  
    real(r8), pointer :: h2osfc             (:)   => null() ! surface water (kg/m2)
    real(r8), pointer :: h2ocan             (:)   => null() ! canopy water integrated to column (kg/m2)
    real(r8), pointer :: total_plant_stored_h2o(:)=> null() ! total water in plants (used??)
    ! Derived water and ice state variables for soil/snow column, depth varying
    real(r8), pointer :: h2osoi_liqvol      (:,:) => null() ! volumetric liquid water content (-nlevsno+1:nlevgrnd) (m3/m3)
    real(r8), pointer :: h2osoi_icevol      (:,:) => null() ! volumetric ice content (-nlevsno+1:nlevgrnd) (m3/m3)     
    real(r8), pointer :: h2osoi_liq_old     (:,:) => null() ! values from previous time step
    real(r8), pointer :: h2osoi_ice_old     (:,:) => null() ! values from previous time step
    real(r8), pointer :: bw                 (:,:) => null() ! partial density of water in the snow pack (ice + liquid) [kg/m3] 
    real(r8), pointer :: smp_l              (:,:) => null() ! liquid phase soil matric potential (-nlevsno+1:nlevgrnd) (mm h2o)
    real(r8), pointer :: soilp              (:,:) => null() ! soil pressure (1:nlevgrnd) (Pa)
    real(r8), pointer :: swe_old            (:,:) => null() ! initial snow water content (-nlevsno+1:0) (kg/m2)
    real(r8), pointer :: snw_rds            (:,:) => null() ! col snow grain radius (-nlevsno+1:0) (m^-6, or microns)
    real(r8), pointer :: air_vol            (:,:) => null() ! air filled porosity (m3/m3)    
    ! Derived water, ice, and snow variables, column aggregate
    real(r8), pointer :: qg_snow            (:)   => null() ! specific humidity over snow (kg H2O/kg moist air)
    real(r8), pointer :: qg_soil            (:)   => null() ! specific humidity over soil (kg H2O/kg moist air)
    real(r8), pointer :: qg_h2osfc          (:)   => null() ! specific humidity over surface water (kg H2O/kg moist air)
    real(r8), pointer :: qg                 (:)   => null() ! average surface specific humidity (kg H2O/kg moist air)
    real(r8), pointer :: dqgdT              (:)   => null() ! rate of change in specific humidity with temperature (kg H2O/kg moist air/K)
    real(r8), pointer :: h2osoi_liqice_10cm (:)   => null() ! liquid water + ice in top 10cm of soil (kg/m2)
    real(r8), pointer :: h2osno             (:)   => null() ! snow mass (kg/m2)
    real(r8), pointer :: h2osno_old         (:)   => null() ! snow mass for previous time step (kg/m2)
    real(r8), pointer :: h2osno_top         (:)   => null() ! top-layer mass of snow  (kg/m2)
    real(r8), pointer :: sno_liq_top        (:)   => null() ! snow liquid water fraction, by mass, top layer (kg/kg)
    real(r8), pointer :: snowice            (:)   => null() ! total snow ice (kg/m2)
    real(r8), pointer :: snowliq            (:)   => null() ! total snow liquid water (kg/m2)
    real(r8), pointer :: int_snow           (:)   => null() ! integrated snowfall (kg/m2)
    real(r8), pointer :: snow_depth         (:)   => null() ! snow height of snow covered area (m)
    real(r8), pointer :: snowdp             (:)   => null() ! snow height averaged for area with and without snow cover(m)
    real(r8), pointer :: snow_persistence   (:)   => null() ! length of time that ground has had non-zero snow thickness (sec)
    real(r8), pointer :: snw_rds_top        (:)   => null() ! snow grain radius (top layer)  (m^-6, microns)
    logical , pointer :: do_capsnow         (:)   => null() ! true => do snow capping
    ! Area fractions
    real(r8), pointer :: frac_sno           (:)   => null() ! fraction of ground covered by snow (0 to 1)
    real(r8), pointer :: frac_sno_eff       (:)   => null() ! fraction of ground covered by snow (0 to 1)
    real(r8), pointer :: frac_iceold        (:,:) => null() ! fraction of ice relative to the tot water (-nlevsno+1:nlevgrnd) 
    real(r8), pointer :: frac_h2osfc        (:)   => null() ! fractional area with surface water greater than zero
    real(r8), pointer :: wf                 (:)   => null() ! soil water as frac. of whc for top 0.05 m (0-1) 
    real(r8), pointer :: wf2                (:)   => null() ! soil water as frac. of whc for top 0.17 m (0-1) 
    real(r8), pointer :: finundated         (:)   => null() ! fraction of column inundated, for bgc caclulation (0-1)
    ! Balance checks
    real(r8), pointer :: begwb              (:)   => null() ! water mass begining of the time step (kg/m2)
    real(r8), pointer :: endwb              (:)   => null() ! water mass end of the time step (kg/m2)
    real(r8), pointer :: errh2o             (:)   => null() ! water conservation error (kg/m2)
    real(r8), pointer :: errh2osno          (:)   => null() ! snow water conservation error(kg/m2)
    real(r8), pointer :: h2osoi_liq_depth_intg(:) => null() ! grid-level depth integrated liquid soil water
    real(r8), pointer :: h2osoi_ice_depth_intg(:) => null() ! grid-level depth integrated ice soil water
    ! VSFM
    real(r8), pointer :: vsfm_fliq_col_1d   (:)   => null() ! fraction of liquid saturation for VSFM [-]
    real(r8), pointer :: vsfm_sat_col_1d    (:)   => null() ! liquid saturation from VSFM [-]
    real(r8), pointer :: vsfm_mass_col_1d   (:)   => null() ! liquid mass per unit area from VSFM [kg H2O/m^2]
    real(r8), pointer :: vsfm_smpl_col_1d   (:)   => null() ! 1D soil matrix potential liquid from VSFM [m]
    real(r8), pointer :: vsfm_soilp_col_1d  (:)   => null() ! 1D soil liquid pressure from VSFM [Pa]
   
  contains
    procedure, public :: Init    => col_ws_init
    procedure, public :: Clean   => col_ws_clean
  end type column_water_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds carbon state information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_carbon_state
    integer :: species  ! c12, c13, c14
    real(r8), pointer :: rootc                (:)    => null() ! (gC/m2) root carbon at column level (fire)
    real(r8), pointer :: totvegc              (:)    => null() ! (gC/m2) column-level totvegc (fire)
    real(r8), pointer :: leafc                (:)    => null() ! (gC/m2) column-level leafc (fire)
    real(r8), pointer :: deadstemc            (:)    => null() ! (gC/m2) column-level deadstemc (fire)
    real(r8), pointer :: fuelc                (:)    => null() ! fuel avalability factor for Reg.C (0-1)
    real(r8), pointer :: fuelc_crop           (:)    => null() ! fuel avalability factor for Reg.A (0-1)
    real(r8), pointer :: decomp_cpools_vr     (:,:,:)=> null() ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
    real(r8), pointer :: ctrunc_vr            (:,:)  => null() ! (gC/m3) vertically-resolved column-level sink for C truncation
    real(r8), pointer :: frootc               (:)    => null() ! (gC/m2) column-level C pool for fine root
    real(r8), pointer :: seedc                (:)    => null() ! (gC/m2) column-level pool for seeding new Patches
    real(r8), pointer :: prod1c               (:)    => null() ! (gC/m2) crop product C pool, 1-year lifespan
    real(r8), pointer :: prod10c              (:)    => null() ! (gC/m2) wood product C pool, 10-year lifespan
    real(r8), pointer :: prod100c             (:)    => null() ! (gC/m2) wood product C pool, 100-year lifespan
    real(r8), pointer :: totprodc             (:)    => null() ! (gC/m2) total wood product C
    real(r8), pointer :: dyn_cbal_adjustments (:)    => null() ! (gC/m2) adjustments to each column made in this timestep via dynamic column area adjustments
    real(r8), pointer :: totpftc              (:)    => null() ! (gC/m2) total patch-level carbon, including cpool averaged to column (p2c)
    real(r8), pointer :: decomp_cpools_1m     (:,:)  => null() ! (gC/m2)  Diagnostic: decomposing (litter, cwd, soil) c pools to 1 meter
    real(r8), pointer :: decomp_cpools        (:,:)  => null() ! (gC/m2)  decomposing (litter, cwd, soil) c pools
    real(r8), pointer :: cwdc                 (:)    => null() ! (gC/m2) Diagnostic: coarse woody debris C
    real(r8), pointer :: ctrunc               (:)    => null() ! (gC/m2) column-level sink for C truncation
    real(r8), pointer :: totlitc              (:)    => null() ! (gC/m2) total litter carbon
    real(r8), pointer :: totsomc              (:)    => null() ! (gC/m2) total soil organic matter carbon
    real(r8), pointer :: totlitc_1m           (:)    => null() ! (gC/m2) total litter carbon to 1 meter
    real(r8), pointer :: totsomc_1m           (:)    => null() ! (gC/m2) total soil organic matter carbon to 1 meter
    real(r8), pointer :: totecosysc           (:)    => null() ! (gC/m2) total ecosystem carbon, incl veg but excl cpool
    real(r8), pointer :: totcolc              (:)    => null() ! (gC/m2) total column carbon, incl veg and cpool
    real(r8), pointer :: totabgc              (:)    => null() ! (gC/m2) total column above ground carbon, excluding som 
    real(r8), pointer :: totblgc              (:)    => null() ! (gc/m2) total column non veg carbon
    real(r8), pointer :: totvegc_abg          (:)    => null() ! (gC/m2) total above vegetation carbon, excluding cpool averaged to column (p2c)
    real(r8), pointer :: begcb                (:)    => null() ! (gC/m2) carbon mass, beginning of time step 
    real(r8), pointer :: endcb                (:)    => null() ! (gc/m2) carbon mass, end of time step 
    real(r8), pointer :: errcb                (:)    => null() ! (gC/m2) carbon balance error for the timestep 
    real(r8), pointer :: totpftc_beg          (:)    => null() 
    real(r8), pointer :: cwdc_beg             (:)    => null() 
    real(r8), pointer :: totlitc_beg          (:)    => null() 
    real(r8), pointer :: totsomc_beg          (:)    => null() 
    real(r8), pointer :: totpftc_end          (:)    => null() 
    real(r8), pointer :: cwdc_end             (:)    => null() 
    real(r8), pointer :: totlitc_end          (:)    => null() 
    real(r8), pointer :: totsomc_end          (:)    => null() 
    real(r8), pointer :: decomp_som2c_vr      (:,:)  => null()
    real(r8), pointer :: cropseedc_deficit    (:)    => null()    

  contains
    procedure, public :: Init    => col_cs_init
    procedure, public :: Clean   => col_cs_clean
  end type column_carbon_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds nitrogen state information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_nitrogen_state
    real(r8), pointer :: decomp_npools_vr         (:,:,:) => null() ! (gN/m3) vertically-resolved decomposing (litter, cwd, soil) N pools
    real(r8), pointer :: decomp_npools            (:,:)   => null() ! (gN/m2)  decomposing (litter, cwd, soil) N pools
    real(r8), pointer :: decomp_npools_1m         (:,:)   => null() ! (gN/m2)  diagnostic: decomposing (litter, cwd, soil) N pools to 1 meter
    real(r8), pointer :: ntrunc_vr                (:,:)   => null() ! (gN/m3) vertically-resolved column-level sink for N truncation
    real(r8), pointer :: sminn_vr                 (:,:)   => null() ! (gN/m3) vertically-resolved soil mineral N
    real(r8), pointer :: smin_no3_vr              (:,:)   => null() ! (gN/m3) vertically-resolved soil mineral NO3
    real(r8), pointer :: smin_nh4_vr              (:,:)   => null() ! (gN/m3) vertically-resolved soil mineral NH4
    real(r8), pointer :: smin_nh4sorb_vr          (:,:)   => null() ! (gN/m3) vertically-resolved soil mineral NH4 absorbed
    real(r8), pointer :: smin_no3                 (:)     => null() ! (gN/m2) soil mineral NO3 pool
    real(r8), pointer :: smin_nh4                 (:)     => null() ! (gN/m2) soil mineral NH4 pool
    real(r8), pointer :: smin_nh4sorb             (:)     => null() ! (gN/m2) soil mineral NH4 pool absorbed
    real(r8), pointer :: sminn                    (:)     => null() ! (gN/m2) soil mineral N
    real(r8), pointer :: ntrunc                   (:)     => null() ! (gN/m2) column-level sink for N truncation
    real(r8), pointer :: cwdn                     (:)     => null() ! (gN/m2) Diagnostic: coarse woody debris N
    real(r8), pointer :: totlitn                  (:)     => null() ! (gN/m2) total litter nitrogen
    real(r8), pointer :: totsomn                  (:)     => null() ! (gN/m2) total soil organic matter nitrogen
    real(r8), pointer :: totlitn_1m               (:)     => null() ! (gN/m2) total litter nitrogen to 1 meter
    real(r8), pointer :: totsomn_1m               (:)     => null() ! (gN/m2) total soil organic matter nitrogen to 1 meter
    real(r8), pointer :: totecosysn               (:)     => null() ! (gN/m2) total ecosystem nitrogen, incl veg 
    real(r8), pointer :: totcoln                  (:)     => null() ! (gN/m2) total column nitrogen, incl veg
    real(r8), pointer :: totabgn                  (:)     => null() ! (gN/m2)
    real(r8), pointer :: totblgn                  (:)     => null() ! (gN/m2) total below ground nitrogen
    real(r8), pointer :: totvegn                  (:)     => null() ! (gN/m2) total vegetation nitrogen (p2c)
    real(r8), pointer :: totpftn                  (:)     => null() ! (gN/m2) total pft-level nitrogen (p2c)
    real(r8), pointer :: plant_n_buffer           (:)     => null() ! (gN/m2) col-level abstract N storage
    real(r8), pointer :: plant_nbuffer            (:)     => null() ! (gN/m2) plant nitrogen buffer, (gN/m2), used to exchange info with betr 
    real(r8), pointer :: seedn                    (:)     => null() ! (gN/m2) column-level pool for seeding new Patches
    real(r8), pointer :: cropseedn_deficit        (:)     => null() ! (gN/m2) column-level pool for seed N deficit (negative pool)
    real(r8), pointer :: prod1n                   (:)     => null() ! (gN/m2) crop product N pool, 1-year lifespan
    real(r8), pointer :: prod10n                  (:)     => null() ! (gN/m2) wood product N pool, 10-year lifespan
    real(r8), pointer :: prod100n                 (:)     => null() ! (gN/m2) wood product N pool, 100-year lifespan
    real(r8), pointer :: totprodn                 (:)     => null() ! (gN/m2) total wood product N
    real(r8), pointer :: dyn_nbal_adjustments     (:)     => null() ! (gN/m2) adjustments to each column made in this timestep via dynamic column area adjustments
    real(r8), pointer :: totpftn_beg              (:)     => null() 
    real(r8), pointer :: totpftn_end              (:)     => null() 
    real(r8), pointer :: cwdn_beg                 (:)     => null() 
    real(r8), pointer :: cwdn_end                 (:)     => null() 
    real(r8), pointer :: totlitn_beg              (:)     => null() 
    real(r8), pointer :: totlitn_end              (:)     => null() 
    real(r8), pointer :: totsomn_beg              (:)     => null() 
    real(r8), pointer :: totsomn_end              (:)     => null() 
    real(r8), pointer :: sminn_beg                (:)     => null() 
    real(r8), pointer :: sminn_end                (:)     => null() 
    real(r8), pointer :: smin_no3_beg             (:)     => null() 
    real(r8), pointer :: smin_no3_end             (:)     => null() 
    real(r8), pointer :: smin_nh4_beg             (:)     => null() 
    real(r8), pointer :: smin_nh4_end             (:)     => null() 
    real(r8), pointer :: totprodn_beg             (:)     => null() 
    real(r8), pointer :: totprodn_end             (:)     => null() 
    real(r8), pointer :: seedn_beg                (:)     => null() 
    real(r8), pointer :: seedn_end                (:)     => null() 
    real(r8), pointer :: ntrunc_beg               (:)     => null() 
    real(r8), pointer :: ntrunc_end               (:)     => null() 
    real(r8), pointer :: begnb                    (:)     => null() ! col nitrogen mass, beginning of time step (gN/m**2)
    real(r8), pointer :: endnb                    (:)     => null() ! col nitrogen mass, end of time step (gN/m**2)
    real(r8), pointer :: errnb                    (:)     => null() ! colnitrogen balance error for the timestep (gN/m**2)
  contains
    procedure, public :: Init       => col_ns_init
    procedure, public :: Clean      => col_ns_clean
  end type column_nitrogen_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds phosphorus state information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_phosphorus_state
    real(r8), pointer :: decomp_ppools_vr         (:,:,:)  => null() ! (gP/m3) vertically-resolved decomposing (litter, cwd, soil) P pools 
    real(r8), pointer :: solutionp_vr             (:,:)    => null() ! (gP/m3) vertically-resolved soil solution P
    real(r8), pointer :: labilep_vr               (:,:)    => null() ! (gP/m3) vertically-resolved soil labile mineral P
    real(r8), pointer :: secondp_vr               (:,:)    => null() ! (gP/m3) vertically-resolved soil secondary mineralP
    real(r8), pointer :: occlp_vr                 (:,:)    => null() ! (gP/m3) vertically-resolved soil occluded mineral P
    real(r8), pointer :: primp_vr                 (:,:)    => null() ! (gP/m3) vertically-resolved soil parimary mineral P
    real(r8), pointer :: sminp_vr                 (:,:)    => null() ! (gP/m3) vertically-resolved soil parimary mineral P
    real(r8), pointer :: ptrunc_vr                (:,:)    => null() ! (gP/m3) vertically-resolved column-level sink for P truncation
    real(r8), pointer :: seedp                    (:)      => null() ! (gP/m2) column-level pool for seeding new Patches
    real(r8), pointer :: prod1p                   (:)      => null() ! (gN/m2) crop product N pool, 1-year lifespan
    real(r8), pointer :: prod10p                  (:)      => null() ! (gP/m2) wood product P pool, 10-year lifespan
    real(r8), pointer :: prod100p                 (:)      => null() ! (gP/m2) wood product P pool, 100-year lifespan
    real(r8), pointer :: totprodp                 (:)      => null() ! (gP/m2) total wood product P
    real(r8), pointer :: dyn_pbal_adjustments     (:)      => null() ! (gP/m2) adjustments to each column made in this timestep via dynamic column area adjustments
    real(r8), pointer :: decomp_ppools            (:,:)    => null() ! (gP/m2)  decomposing (litter, cwd, soil) P pools
    real(r8), pointer :: decomp_ppools_1m         (:,:)    => null() ! (gP/m2)  diagnostic: decomposing (litter, cwd, soil) P pools to 1 meter
    real(r8), pointer :: sminp                    (:)      => null() ! (gP/m2) soil mineral P
    real(r8), pointer :: solutionp                (:)      => null() ! (gP/m2) soil solution P
    real(r8), pointer :: labilep                  (:)      => null() ! (gP/m2) soil labile mineral P
    real(r8), pointer :: secondp                  (:)      => null() ! (gP/m2) soil secondary mineralP
    real(r8), pointer :: occlp                    (:)      => null() ! (gP/m2) soil occluded mineral P
    real(r8), pointer :: primp                    (:)      => null() ! (gP/m2) soil parimary mineral P
    real(r8), pointer :: ptrunc                   (:)      => null() ! (gP/m2) column-level sink for P truncation
    real(r8), pointer :: cwdp                     (:)      => null() ! (gP/m2) Diagnostic: coarse woody debris P
    real(r8), pointer :: totlitp                  (:)      => null() ! (gP/m2) total litter phosphorus
    real(r8), pointer :: totsomp                  (:)      => null() ! (gP/m2) total soil organic matter phosphorus
    real(r8), pointer :: totlitp_1m               (:)      => null() ! (gP/m2) total litter phosphorus to 1 meter
    real(r8), pointer :: totsomp_1m               (:)      => null() ! (gP/m2) total soil organic matter phosphorus to 1 meter
    real(r8), pointer :: totecosysp               (:)      => null() ! (gP/m2) total ecosystem phosphorus, incl veg 
    real(r8), pointer :: totcolp                  (:)      => null() ! (gP/m2) total column phosphorus, incl veg
    real(r8), pointer :: totvegp                  (:)      => null() ! (gP/m2) total vegetation phosphorus (p2c)
    real(r8), pointer :: totpftp                  (:)      => null() ! (gP/m2) total pft-level phosphorus (p2c)
    real(r8), pointer :: begpb                    (:)      => null() ! phosphorus mass, beginning of time step (gP/m**2)
    real(r8), pointer :: endpb                    (:)      => null() ! phosphorus mass, end of time step (gP/m**2)
    real(r8), pointer :: errpb                    (:)      => null() ! phosphorus balance error for the timestep (gP/m**2)
    real(r8), pointer :: solutionp_vr_cur         (:,:)    => null() 
    real(r8), pointer :: solutionp_vr_prev        (:,:)    => null() 
    real(r8), pointer :: labilep_vr_cur           (:,:)    => null() 
    real(r8), pointer :: labilep_vr_prev          (:,:)    => null() 
    real(r8), pointer :: secondp_vr_cur           (:,:)    => null() 
    real(r8), pointer :: secondp_vr_prev          (:,:)    => null() 
    real(r8), pointer :: occlp_vr_cur             (:,:)    => null() 
    real(r8), pointer :: occlp_vr_prev            (:,:)    => null() 
    real(r8), pointer :: primp_vr_cur             (:,:)    => null() 
    real(r8), pointer :: primp_vr_prev            (:,:)    => null() 
    real(r8), pointer :: totpftp_beg              (:)      => null() 
    real(r8), pointer :: solutionp_beg            (:)      => null() 
    real(r8), pointer :: labilep_beg              (:)      => null() 
    real(r8), pointer :: secondp_beg              (:)      => null() 
    real(r8), pointer :: totlitp_beg              (:)      => null() 
    real(r8), pointer :: cwdp_beg                 (:)      => null() 
    real(r8), pointer :: totsomp_beg              (:)      => null() 
    real(r8), pointer :: totlitp_end              (:)      => null() 
    real(r8), pointer :: totpftp_end              (:)      => null() 
    real(r8), pointer :: solutionp_end            (:)      => null() 
    real(r8), pointer :: labilep_end              (:)      => null() 
    real(r8), pointer :: secondp_end              (:)      => null() 
    real(r8), pointer :: cwdp_end                 (:)      => null() 
    real(r8), pointer :: totsomp_end              (:)      => null()
    real(r8), pointer :: cropseedp_deficit        (:)      => null() ! (gP/m2) negative pool tracking seed P for DWT    
  contains
    procedure, public :: Init      => col_ps_init
    procedure, public :: Clean     => col_ps_clean
  end type column_phosphorus_state

  !-----------------------------------------------------------------------
  ! Define the data structure that holds energy flux information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_energy_flux
    real(r8), pointer :: eflx_h2osfc_to_snow     (:)   => null() ! snow melt to h2osfc heat flux (W/m2)
    real(r8), pointer :: eflx_snomelt            (:)   => null() ! snow melt heat flux (W/m**2)
    real(r8), pointer :: eflx_snomelt_r          (:)   => null() ! rural snow melt heat flux (W/m2)
    real(r8), pointer :: eflx_snomelt_u          (:)   => null() ! urban snow melt heat flux (W/m2)
    real(r8), pointer :: eflx_bot                (:)   => null() ! heat flux from beneath the soil or ice column (W/m2)
    real(r8), pointer :: eflx_fgr12              (:)   => null() ! ground heat flux between soil layers 1 and 2 (W/m2)
    real(r8), pointer :: eflx_fgr                (:,:) => null() ! (rural) soil downward heat flux (W/m2) (1:nlevgrnd)  (pos upward; usually eflx_bot >= 0)
    real(r8), pointer :: eflx_building_heat      (:)   => null() ! heat flux from urban building interior to urban walls, roof (W/m2)
    real(r8), pointer :: eflx_urban_ac           (:)   => null() ! urban air conditioning flux (W/m2)
    real(r8), pointer :: eflx_urban_heat         (:)   => null() ! urban heating flux (W/m**2)
    real(r8), pointer :: eflx_hs_h2osfc          (:)   => null() ! heat flux on standing water (W/m2)
    real(r8), pointer :: eflx_hs_top_snow        (:)   => null() ! heat flux on top snow layer (W/m2)
    real(r8), pointer :: eflx_hs_soil            (:)   => null() ! heat flux on soil [W/m2
    real(r8), pointer :: eflx_sabg_lyr           (:,:) => null() ! absorbed solar radiation (col,lyr) (W/m2)
    ! Derivatives of energy fluxes                      
    real(r8), pointer :: eflx_dhsdT              (:)   => null() ! deriv. of energy flux into surface layer wrt temp (W/m2/K)
    ! Latent heat terms
    real(r8), pointer :: htvp                    (:)   => null() ! latent heat of vapor of water (or sublimation) [j/kg]
    real(r8), pointer :: xmf                     (:)   => null() ! total latent heat of phase change of ground water
    real(r8), pointer :: xmf_h2osfc              (:)   => null() ! latent heat of phase change of surface water
    integer , pointer :: imelt                   (:,:) ! flag for melting (=1), freezing (=2), Not=0 (-nlevsno+1:nlevgrnd) 
    ! for couplig with pflotran                         
    real(r8), pointer :: eflx_soil_grnd          (:)   => null() ! integrated soil ground heat flux (W/m2)  [+ = into ground]
    real(r8), pointer :: eflx_rnet_soil          (:)   => null() ! soil net (sw+lw) radiation flux (W/m2) [+ = into soil]
    real(r8), pointer :: eflx_fgr0_soil          (:)   => null() ! soil-air heat flux (W/m2) [+ = into soil]
    real(r8), pointer :: eflx_fgr0_snow          (:)   => null() ! soil-snow heat flux (W/m2) [+ = into soil]
    real(r8), pointer :: eflx_fgr0_h2osfc        (:)   => null() ! soil-surfacewater heat flux (W/m2) [+ = into soil]
    ! Balance Checks                                    
    real(r8), pointer :: errsoi                  (:)   => null() ! soil/lake energy conservation error   (W/m2)
    real(r8), pointer :: errseb                  (:)   => null() ! surface energy conservation error     (W/m2)
    real(r8), pointer :: errsol                  (:)   => null() ! solar radiation conservation error    (W/m2)
    real(r8), pointer :: errlon                  (:)   => null() ! longwave radiation conservation error (W/m2)

  contains
    procedure, public :: Init    => col_ef_init
    procedure, public :: Clean   => col_ef_clean
  end type column_energy_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds water flux information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_water_flux
    real(r8), pointer :: qflx_prec_grnd       (:)   => null() ! water onto ground including canopy runoff [kg/(m2 s)]
    real(r8), pointer :: qflx_rain_grnd       (:)   => null() ! rain on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: qflx_snow_grnd       (:)   => null() ! snow on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: qflx_sub_snow        (:)   => null() ! sublimation rate from snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_sub_snow_vol    (:)   => null() !
    real(r8), pointer :: qflx_evap_soi        (:)   => null() ! soil evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_veg        (:)   => null() ! vegetation evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_can        (:)   => null() ! evaporation from leaves and stems (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_tot        (:)   => null() ! col_qflx_evap_soi + col_qflx_evap_veg + qflx_tran_veg
    real(r8), pointer :: qflx_evap_grnd       (:)   => null() ! ground surface evaporation rate (mm H2O/s) [+]
    real(r8), pointer :: qflx_snwcp_liq       (:)   => null() ! excess rainfall due to snow capping (mm H2O /s)
    real(r8), pointer :: qflx_snwcp_ice       (:)   => null() ! excess snowfall due to snow capping (mm H2O /s)
    real(r8), pointer :: qflx_tran_veg        (:)   => null() ! vegetation transpiration (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_dew_snow        (:)   => null() ! surface dew added to snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_dew_grnd        (:)   => null() ! ground surface dew formation (mm H2O /s) [+] (+ = to atm); usually eflx_bot >= 0)
    real(r8), pointer :: qflx_prec_intr       (:)   => null() ! interception of precipitation [mm/s]
    real(r8), pointer :: qflx_ev_snow         (:)   => null() ! evaporation heat flux from snow         (W/m**2) [+ to atm] ! NOTE: unit shall be mm H2O/s for water NOT heat
    real(r8), pointer :: qflx_ev_soil         (:)   => null() ! evaporation heat flux from soil         (W/m**2) [+ to atm] ! NOTE: unit shall be mm H2O/s for water NOT heat
    real(r8), pointer :: qflx_ev_h2osfc       (:)   => null() ! evaporation heat flux from soil         (W/m**2) [+ to atm] ! NOTE: unit shall be mm H2O/s for water NOT heat
    real(r8), pointer :: qflx_gross_evap_soil (:)   => null() ! gross infiltration from soil, this satisfies the relationship qflx_infl = qflx_gross_infl_soil-qflx_gross_evap_soil
    real(r8), pointer :: qflx_gross_infl_soil (:)   => null() ! gross infiltration, before considering the evaporation
    real(r8), pointer :: qflx_adv             (:,:) => null() ! advective flux across different soil layer interfaces [mm H2O/s] [+ downward]
    real(r8), pointer :: qflx_rootsoi         (:,:) => null() ! root and soil water exchange [mm H2O/s] [+ into root]     
    real(r8), pointer :: dwb                  (:)   => null() !  water mass change [+ increase] [mm H2O/s] 
    real(r8), pointer :: qflx_infl            (:)   => null() ! infiltration (mm H2O /s)
    real(r8), pointer :: qflx_surf            (:)   => null() ! surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_drain           (:)   => null() ! sub-surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_totdrain        (:)   => null() 
    real(r8), pointer :: qflx_top_soil        (:)   => null() ! net water input into soil from top (mm/s)
    real(r8), pointer :: qflx_h2osfc_to_ice   (:)   => null() ! conversion of h2osfc to ice
    real(r8), pointer :: qflx_h2osfc_surf     (:)   => null() ! surface water runoff
    real(r8), pointer :: qflx_snow_h2osfc     (:)   => null() ! snow falling on surface water
    real(r8), pointer :: qflx_drain_perched   (:)   => null() ! sub-surface runoff from perched wt (mm H2O /s)
    real(r8), pointer :: qflx_deficit         (:)   => null() ! water deficit to keep non-negative liquid water content (mm H2O)   
    real(r8), pointer :: qflx_floodc          (:)   => null() ! flood water flux at column level
    real(r8), pointer :: qflx_sl_top_soil     (:)   => null() ! liquid water + ice from layer above soil to top soil layer or sent to qflx_qrgwl (mm H2O/s)
    real(r8), pointer :: qflx_snomelt         (:)   => null() ! snow melt (mm H2O /s)
    real(r8), pointer :: qflx_snow_melt       (:)   => null() ! snow melt (net)
    real(r8), pointer :: qflx_qrgwl           (:)   => null() ! qflx_surf at glaciers, wetlands, lakes
    real(r8), pointer :: qflx_runoff          (:)   => null() ! total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
    real(r8), pointer :: qflx_runoff_r        (:)   => null() ! Rural total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
    real(r8), pointer :: qflx_runoff_u        (:)   => null() ! urban total runoff (qflx_drain+qflx_surf) (mm H2O /s) 
    real(r8), pointer :: qflx_rsub_sat        (:)   => null() ! soil saturation excess [mm/s]
    real(r8), pointer :: qflx_snofrz_lyr      (:,:) => null() ! snow freezing rate (positive definite) (col,lyr) [kg m-2 s-1]
    real(r8), pointer :: qflx_snofrz          (:)   => null() ! column-integrated snow freezing rate (positive definite) (col) [kg m-2 s-1]
    real(r8), pointer :: qflx_glcice          (:)   => null() ! net flux of new glacial ice (growth - melt) (mm H2O/s), passed to GLC
    real(r8), pointer :: qflx_glcice_frz      (:)   => null() ! ice growth (positive definite) (mm H2O/s)
    real(r8), pointer :: qflx_glcice_melt     (:)   => null() ! ice melt (positive definite) (mm H2O/s)
    real(r8), pointer :: qflx_drain_vr        (:,:) => null() ! liquid water lost as drainage (m /time step)
    real(r8), pointer :: qflx_h2osfc2topsoi   (:)   => null() ! liquid water coming from surface standing water top soil (mm H2O/s)
    real(r8), pointer :: qflx_snow2topsoi     (:)   => null() ! liquid water coming from residual snow to topsoil (mm H2O/s)
    real(r8), pointer :: qflx_lateral         (:)   => null() ! lateral subsurface flux (mm H2O /s)
    real(r8), pointer :: snow_sources         (:)   => null() ! snow sources (mm H2O/s)
    real(r8), pointer :: snow_sinks           (:)   => null() ! snow sinks (mm H2O/s)
    real(r8), pointer :: qflx_irrig           (:)   => null() ! irrigation flux (mm H2O/s)
    real(r8), pointer :: mflx_infl_1d         (:)   => null() ! infiltration source in top soil control volume (kg H2O /s)
    real(r8), pointer :: mflx_dew_1d          (:)   => null() ! liquid+snow dew source in top soil control volume (kg H2O /s)
    real(r8), pointer :: mflx_et_1d           (:)   => null() ! evapotranspiration sink from all soil coontrol volumes (kg H2O /s)
    real(r8), pointer :: mflx_drain_1d        (:)   => null() ! drainage from groundwater table (kg H2O /s)
    real(r8), pointer :: mflx_drain_perched_1d(:)   => null() ! drainage from perched water table (kg H2O /s)
    real(r8), pointer :: mflx_snowlyr_1d      (:)   => null() ! mass flux to top soil layer due to disappearance of snow (kg H2O /s)
    real(r8), pointer :: mflx_sub_snow_1d     (:)   => null() ! mass flux from top soil layer due to sublimation of snow (kg H2O /s)
    real(r8), pointer :: mflx_neg_snow_1d     (:)   => null() ! mass flux from top soil layer due to negative water content in snow layers (kg H2O /s)
    real(r8), pointer :: mflx_snowlyr         (:)   => null() ! mass flux to top soil layer due to disappearance of snow (kg H2O /s). This is for restart
    real(r8), pointer :: mflx_infl            (:)   => null() ! infiltration source in top soil control volume (kg H2O /s)
    real(r8), pointer :: mflx_dew             (:)   => null() ! liquid+snow dew source in top soil control volume (kg H2O /s)
    real(r8), pointer :: mflx_snowlyr_disp    (:)   => null() ! mass flux to top soil layer due to disappearance of snow (kg H2O /s)
    real(r8), pointer :: mflx_sub_snow        (:)   => null() ! mass flux from top soil layer due to sublimation of snow (kg H2O /s)
    real(r8), pointer :: mflx_et              (:,:) => null() ! evapotranspiration sink from all soil coontrol volumes (kg H2O /s)
    real(r8), pointer :: mflx_drain           (:,:) => null() ! drainage from groundwater table (kg H2O /s)
    real(r8), pointer :: mflx_recharge        (:)   => null() ! recharge from soil column to unconfined aquifer (kg H2O /s)

  contains
    procedure, public :: Init    => col_wf_init
    procedure, public :: Clean   => col_wf_clean
  end type column_water_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds carbon flux information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_carbon_flux
    ! primary decomposition and vertical movement fluxes
    real(r8), pointer :: decomp_cpools_sourcesink              (:,:,:) => null() ! change in decomposing c pools. Used to update concentrations concurrently with vertical transport (gC/m3/timestep)  
    real(r8), pointer :: decomp_cascade_hr_vr                  (:,:,:) => null() ! vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
    real(r8), pointer :: decomp_cascade_ctransfer_vr           (:,:,:) => null() ! vertically-resolved C transferred along deomposition cascade (gC/m3/s)
    real(r8), pointer :: decomp_k                              (:,:,:) => null() ! rate constant for decomposition (1./sec)
    real(r8), pointer :: decomp_cpools_transport_tendency      (:,:,:) => null() ! C tendency due to vertical transport in decomposing C pools (gC/m^3/s)
    real(r8), pointer :: decomp_cascade_hr                     (:,:)   => null() ! vertically-integrated (diagnostic) het. resp. from decomposing C pools (gC/m2/s)
    real(r8), pointer :: decomp_cascade_ctransfer              (:,:)   => null() ! vertically-integrated (diagnostic) C transferred along deomposition cascade (gC/m2/s)
    real(r8), pointer :: hr_vr                                 (:,:)   => null() ! total vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
    real(r8), pointer :: o_scalar                              (:,:)   => null() ! fraction by which decomposition is limited by anoxia
    real(r8), pointer :: w_scalar                              (:,:)   => null() ! fraction by which decomposition is limited by moisture availability
    real(r8), pointer :: t_scalar                              (:,:)   => null() ! fraction by which decomposition is limited by temperature
    real(r8), pointer :: decomp_cpools_leached                 (:,:)   => null() ! C loss from vertical transport from each decomposing C pool (gC/m^2/s)
    real(r8), pointer :: phr_vr                                (:,:)   => null() ! potential hr (not N-limited) (gC/m3/s)
    real(r8), pointer :: fphr                                  (:,:)   => null() ! fraction of potential heterotrophic respiration
    real(r8), pointer :: som_c_leached                         (:)     => null() ! total SOM C loss from vertical transport (gC/m^2/s)
    ! phenology: litterfall and crop fluxes
    real(r8), pointer :: phenology_c_to_litr_met_c             (:,:)   => null() ! C fluxes associated with phenology (litterfall and crop) to litter metabolic pool (gC/m3/s)
    real(r8), pointer :: phenology_c_to_litr_cel_c             (:,:)   => null() ! C fluxes associated with phenology (litterfall and crop) to litter cellulose pool (gC/m3/s)
    real(r8), pointer :: phenology_c_to_litr_lig_c             (:,:)   => null() ! C fluxes associated with phenology (litterfall and crop) to litter lignin pool (gC/m3/s)
    ! gap mortality
    real(r8), pointer :: gap_mortality_c_to_litr_met_c         (:,:)   => null() ! C fluxes associated with gap mortality to litter metabolic pool (gC/m3/s)
    real(r8), pointer :: gap_mortality_c_to_litr_cel_c         (:,:)   => null() ! C fluxes associated with gap mortality to litter cellulose pool (gC/m3/s)
    real(r8), pointer :: gap_mortality_c_to_litr_lig_c         (:,:)   => null() ! C fluxes associated with gap mortality to litter lignin pool (gC/m3/s)
    real(r8), pointer :: gap_mortality_c_to_cwdc               (:,:)   => null() ! C fluxes associated with gap mortality to CWD pool (gC/m3/s)
    ! column-level fire fluxes
    real(r8), pointer :: m_decomp_cpools_to_fire_vr            (:,:,:) => null() ! vertically-resolved decomposing C fire loss (gC/m3/s)
    real(r8), pointer :: m_decomp_cpools_to_fire               (:,:)   => null() ! vertically-integrated (diagnostic) decomposing C fire loss (gC/m2/s)
    real(r8), pointer :: m_c_to_litr_met_fire                  (:,:)   => null() ! C from leaf, froot, xfer and storage C to litter labile C by fire (gC/m3/s) 
    real(r8), pointer :: m_c_to_litr_cel_fire                  (:,:)   => null() ! C from leaf, froot, xfer and storage C to litter cellulose C by fire (gC/m3/s) 
    real(r8), pointer :: m_c_to_litr_lig_fire                  (:,:)   => null() ! C from leaf, froot, xfer and storage C to litter lignin C by fire (gC/m3/s) 
    real(r8), pointer :: fire_mortality_c_to_cwdc              (:,:)   => null() ! C fluxes associated with fire mortality to CWD pool (gC/m3/s)
    real(r8), pointer :: somc_fire                             (:)     => null() ! (gC/m2/s) carbon emissions due to peat burning
    ! dynamic LULCC fluxes: harvest
    real(r8), pointer :: harvest_c_to_litr_met_c               (:,:)   => null() ! C fluxes associated with harvest to litter metabolic pool (gC/m3/s)
    real(r8), pointer :: harvest_c_to_litr_cel_c               (:,:)   => null() ! C fluxes associated with harvest to litter cellulose pool (gC/m3/s)
    real(r8), pointer :: harvest_c_to_litr_lig_c               (:,:)   => null() ! C fluxes associated with harvest to litter lignin pool (gC/m3/s)
    real(r8), pointer :: harvest_c_to_cwdc                     (:,:)   => null() ! C fluxes associated with harvest to CWD pool (gC/m3/s)
    real(r8), pointer :: hrv_deadstemc_to_prod10c              (:)     => null() ! dead stem C harvest mortality to 10-year product pool (gC/m2/s)        
    real(r8), pointer :: hrv_deadstemc_to_prod100c             (:)     => null() ! dead stem C harvest mortality to 100-year product pool (gC/m2/s)        
    real(r8), pointer :: hrv_cropc_to_prod1c                   (:)     => null() ! crop C harvest mortality to 1-year product pool (gC/m2/s)
    ! dynamic LULCC fluxes: land cover change
    real(r8), pointer :: dwt_frootc_to_litr_met_c              (:,:)   => null() ! (gC/m3/s) fine root to litter due to landcover change
    real(r8), pointer :: dwt_frootc_to_litr_cel_c              (:,:)   => null() ! (gC/m3/s) fine root to litter due to landcover change
    real(r8), pointer :: dwt_frootc_to_litr_lig_c              (:,:)   => null() ! (gC/m3/s) fine root to litter due to landcover change
    real(r8), pointer :: dwt_livecrootc_to_cwdc                (:,:)   => null() ! (gC/m3/s) live coarse root to CWD due to landcover change
    real(r8), pointer :: dwt_deadcrootc_to_cwdc                (:,:)   => null() ! (gC/m3/s) dead coarse root to CWD due to landcover change
    real(r8), pointer :: dwt_slash_cflux                       (:)     => null() ! (gC/m2/s) conversion slash flux due to landcover change
    real(r8), pointer :: dwt_conv_cflux                        (:)     => null() ! (gC/m2/s) conversion C flux (immediate loss to atm)
    real(r8), pointer :: dwt_prod10c_gain                      (:)     => null() ! (gC/m2/s) addition to 10-yr wood product pool
    real(r8), pointer :: dwt_prod100c_gain                     (:)     => null() ! (gC/m2/s) addition to 100-yr wood product pool
    real(r8), pointer :: dwt_closs                             (:)     => null() ! (gC/m2/s) total carbon loss from product pools and conversion
    real(r8), pointer :: landuseflux                           (:)     => null() ! (gC/m2/s) dwt_closs+product_closs
    real(r8), pointer :: landuptake                            (:)     => null() ! (gC/m2/s) nee-landuseflux
    ! wood product pool loss fluxes
    real(r8), pointer :: prod1c_loss                           (:)     => null() ! (gC/m2/s) decomposition loss from 1-year product pool
    real(r8), pointer :: prod10c_loss                          (:)     => null() ! (gC/m2/s) decomposition loss from 10-yr wood product pool
    real(r8), pointer :: prod100c_loss                         (:)     => null() ! (gC/m2/s) decomposition loss from 100-yr wood product pool
    real(r8), pointer :: product_closs                         (:)     => null() ! (gC/m2/s) total wood product carbon loss
    ! summary (diagnostic) flux variables, not involved in mass balance
    real(r8), pointer :: lithr                                 (:)     => null() ! (gC/m2/s) litter heterotrophic respiration 
    real(r8), pointer :: somhr                                 (:)     => null() ! (gC/m2/s) soil organic matter heterotrophic respiration
    real(r8), pointer :: hr                                    (:)     => null() ! (gC/m2/s) total heterotrophic respiration
    real(r8), pointer :: sr                                    (:)     => null() ! (gC/m2/s) total soil respiration (HR + root resp)
    real(r8), pointer :: er                                    (:)     => null() ! (gC/m2/s) total ecosystem respiration, autotrophic + heterotrophic
    real(r8), pointer :: litfire                               (:)     => null() ! (gC/m2/s) litter fire losses
    real(r8), pointer :: somfire                               (:)     => null() ! (gC/m2/s) soil organic matter fire losses
    real(r8), pointer :: totfire                               (:)     => null() ! (gC/m2/s) total ecosystem fire losses
    real(r8), pointer :: nep                                   (:)     => null() ! (gC/m2/s) net ecosystem production, excludes fire, landuse, and harvest flux, positive for sink
    real(r8), pointer :: nbp                                   (:)     => null() ! (gC/m2/s) net biome production, includes fire, landuse, and harvest flux, positive for sink
    real(r8), pointer :: nee                                   (:)     => null() ! (gC/m2/s) net ecosystem exchange of carbon, includes fire, landuse, harvest, and hrv_xsmrpool flux, positive for source
    ! CLAMP summary (diagnostic) flux variables, not involved in mass balance
    real(r8), pointer :: bgc_cpool_ext_inputs_vr               (:,:,:) => null() ! col-level extneral organic carbon input gC/m3 /time step
    real(r8), pointer :: bgc_cpool_ext_loss_vr                 (:,:,:) => null() ! col-level extneral organic carbon loss gC/m3 /time step
    real(r8), pointer :: cwdc_hr                               (:)     => null() ! (gC/m2/s) col-level coarse woody debris C heterotrophic respiration
    real(r8), pointer :: cwdc_loss                             (:)     => null() ! (gC/m2/s) col-level coarse woody debris C loss
    real(r8), pointer :: litterc_loss                          (:)     => null() ! (gC/m2/s) col-level litter C loss
    ! patch averaged to column variables - to remove need for pcf_a instance
    real(r8), pointer :: rr                                    (:)     => null() ! column (gC/m2/s) root respiration (fine root MR + total root GR) (p2c)
    real(r8), pointer :: ar                                    (:)     => null() ! column (gC/m2/s) autotrophic respiration (MR + GR) (p2c)      
    real(r8), pointer :: gpp                                   (:)     => null() ! column (gC/m2/s) GPP flux before downregulation  (p2c)         
    real(r8), pointer :: npp                                   (:)     => null() ! column (gC/m2/s) net primary production (p2c)                  
    real(r8), pointer :: fire_closs_p2c                        (:)     => null() ! column (gC/m2/s) patch2col averaged column-level fire C loss (p2c)
    real(r8), pointer :: fire_closs                            (:)     => null() ! column (gC/m2/s) total patch-level fire C loss 
    real(r8), pointer :: fire_decomp_closs                     (:)     => null() ! column (gC/m2/s) carbon loss to fire for decomposable pools
    real(r8), pointer :: litfall                               (:)     => null() ! column (gC/m2/s) total patch-level litterfall C loss (p2c)       
    real(r8), pointer :: vegfire                               (:)     => null() ! column (gC/m2/s) patch-level fire loss (obsolete, mark for removal) (p2c)
    real(r8), pointer :: wood_harvestc                         (:)     => null() ! column (p2c)                                                  
    real(r8), pointer :: hrv_xsmrpool_to_atm                   (:)     => null() ! column excess MR pool harvest mortality (gC/m2/s) (p2c)
    real(r8), pointer :: plant_to_litter_cflux                 (:)     => null() ! for the purpose of mass balance check
    real(r8), pointer :: plant_to_cwd_cflux		               (:)     => null() ! for the purpose of mass balance check
    ! Temporary and annual sums
    real(r8), pointer :: annsum_npp                            (:)     => null() ! col annual sum of NPP, averaged from pft-level (gC/m2/yr)
    real(r8), pointer :: lag_npp                               (:)     => null() ! col lagged net primary production (gC/m2/s)
    ! Variables for elm_interface_funcsMod & pflotran
    real(r8), pointer :: externalc_to_decomp_cpools            (:,:,:) => null() ! col (gC/m3/s) net C fluxes associated with litter/som-adding/removal to decomp pools
    real(r8), pointer :: externalc_to_decomp_delta             (:)     => null() ! col (gC/m2) summarized net change of whole column C i/o to decomposing pool bwtn time-step
    real(r8), pointer :: f_co2_soil_vr                         (:,:)   => null() ! total vertically-resolved soil-atm. CO2 exchange (gC/m3/s)
    real(r8), pointer :: f_co2_soil                            (:)     => null() ! total soil-atm. CO2 exchange (gC/m2/s)

  contains
    procedure, public :: Init       => col_cf_init
    procedure, public :: Clean      => col_cf_clean
  end type column_carbon_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds nitrogen flux information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_nitrogen_flux
    ! harvest fluxes
    real(r8), pointer :: hrv_deadstemn_to_prod10n              (:)     => null() ! dead stem N harvest mortality to 10-year product pool (gN/m2/s)
    real(r8), pointer :: hrv_deadstemn_to_prod100n             (:)     => null() ! dead stem N harvest mortality to 100-year product pool (gN/m2/s)
    real(r8), pointer :: m_n_to_litr_met_fire                  (:,:)   => null() ! N from leaf, froot, xfer and storage N to litter labile N by fire (gN/m3/s) 
    real(r8), pointer :: m_n_to_litr_cel_fire                  (:,:)   => null() ! N from leaf, froot, xfer and storage N to litter cellulose N by fire (gN/m3/s) 
    real(r8), pointer :: m_n_to_litr_lig_fire                  (:,:)   => null() ! N from leaf, froot, xfer and storage N to litter lignin N by fire (gN/m3/s) 
    real(r8), pointer :: harvest_n_to_litr_met_n               (:,:)   => null() ! N fluxes associated with harvest to litter metabolic pool (gN/m3/s)
    real(r8), pointer :: harvest_n_to_litr_cel_n               (:,:)   => null() ! N fluxes associated with harvest to litter cellulose pool (gN/m3/s)
    real(r8), pointer :: harvest_n_to_litr_lig_n               (:,:)   => null() ! N fluxes associated with harvest to litter lignin pool (gN/m3/s)
    real(r8), pointer :: harvest_n_to_cwdn                     (:,:)   => null() ! N fluxes associated with harvest to CWD pool (gN/m3/s)
    real(r8), pointer :: hrv_cropn_to_prod1n                   (:)     => null() ! crop N harvest mortality to 1-yr product pool (gN/m2/s)
    ! fire N fluxes 
    real(r8), pointer :: m_decomp_npools_to_fire_vr            (:,:,:) => null() ! vertically-resolved decomposing N fire loss (gN/m3/s)
    real(r8), pointer :: m_decomp_npools_to_fire               (:,:)   => null() ! vertically-integrated (diagnostic) decomposing N fire loss (gN/m2/s)
    real(r8), pointer :: fire_nloss                            (:)     => null() ! total column-level fire N loss (gN/m2/s)
    real(r8), pointer :: fire_decomp_nloss                     (:)     => null() ! fire N loss from decomposable pools (gN/m2/s)
    real(r8), pointer :: fire_nloss_p2c                        (:)     => null() ! patch2col column-level fire N loss (gN/m2/s) (p2c)
    real(r8), pointer :: fire_mortality_n_to_cwdn              (:,:)   => null() ! N fluxes associated with fire mortality to CWD pool (gN/m3/s)
    ! summary (diagnostic) flux variables, not involved in mass balance
    real(r8), pointer :: wood_harvestn                         (:)     => null() ! total N losses to wood product pools (gN/m2/s) (p2c)
    ! deposition fluxes
    real(r8), pointer :: ndep_to_sminn                         (:)     => null() ! atmospheric N deposition to soil mineral N (gN/m2/s)
    real(r8), pointer :: nfix_to_sminn                         (:)     => null() ! symbiotic/asymbiotic N fixation to soil mineral N (gN/m2/s) 
    real(r8), pointer :: nfix_to_ecosysn                       (:)     => null() ! total nitrogen fixation
    real(r8), pointer :: fert_to_sminn                         (:)     => null() ! fertilizer N to soil mineral N (gN/m2/s)
    real(r8), pointer :: soyfixn_to_sminn                      (:)     => null() ! soybean fixation to soil mineral N (gN/m2/s)
    ! phenology: litterfall and crop fluxes
    real(r8), pointer :: phenology_n_to_litr_met_n             (:,:)   => null() ! N fluxes associated with phenology (litterfall and crop) to litter metabolic pool (gN/m3/s)
    real(r8), pointer :: phenology_n_to_litr_cel_n             (:,:)   => null() ! N fluxes associated with phenology (litterfall and crop) to litter cellulose pool (gN/m3/s)
    real(r8), pointer :: phenology_n_to_litr_lig_n             (:,:)   => null() ! N fluxes associated with phenology (litterfall and crop) to litter lignin pool (gN/m3/s)
    ! gap mortality fluxes
    real(r8), pointer :: gap_mortality_n_to_litr_met_n         (:,:)   => null() ! N fluxes associated with gap mortality to litter metabolic pool (gN/m3/s)
    real(r8), pointer :: gap_mortality_n_to_litr_cel_n         (:,:)   => null() ! N fluxes associated with gap mortality to litter cellulose pool (gN/m3/s)
    real(r8), pointer :: gap_mortality_n_to_litr_lig_n         (:,:)   => null() ! N fluxes associated with gap mortality to litter lignin pool (gN/m3/s)
    real(r8), pointer :: gap_mortality_n_to_cwdn               (:,:)   => null() ! N fluxes associated with gap mortality to CWD pool (gN/m3/s)
    ! decomposition fluxes
    real(r8), pointer :: decomp_cascade_ntransfer_vr           (:,:,:) => null() ! vert-res transfer of N from donor to receiver pool along decomp. cascade (gN/m3/s)
    real(r8), pointer :: decomp_cascade_ntransfer              (:,:)   => null() ! vert-int (diagnostic) transfer of N from donor to receiver pool along decomp. cascade (gN/m2/s)
    real(r8), pointer :: decomp_cascade_sminn_flux_vr          (:,:,:) => null() ! vert-res mineral N flux for transition along decomposition cascade (gN/m3/s)
    real(r8), pointer :: decomp_cascade_sminn_flux             (:,:)   => null() ! vert-int (diagnostic) mineral N flux for transition along decomposition cascade (gN/m2/s)
    ! vertically-resolved immobilization fluxes
    real(r8), pointer :: potential_immob_vr                    (:,:)   => null() ! vertically-resolved potential N immobilization (gN/m3/s) at each level
    real(r8), pointer :: potential_immob                       (:)     => null() ! vert-int (diagnostic) potential N immobilization (gN/m2/s)
    real(r8), pointer :: actual_immob_vr                       (:,:)   => null() ! vertically-resolved actual N immobilization (gN/m3/s) at each level
    real(r8), pointer :: actual_immob                          (:)     => null() ! vert-int (diagnostic) actual N immobilization (gN/m2/s)
    real(r8), pointer :: sminn_to_plant_vr                     (:,:)   => null() ! vertically-resolved plant uptake of soil mineral N (gN/m3/s)
    real(r8), pointer :: sminn_to_plant                        (:)     => null() ! vert-int (diagnostic) plant uptake of soil mineral N (gN/m2/s)
    real(r8), pointer :: supplement_to_sminn_vr                (:,:)   => null() ! vertically-resolved supplemental N supply (gN/m3/s)
    real(r8), pointer :: supplement_to_sminn                   (:)     => null() ! vert-int (diagnostic) supplemental N supply (gN/m2/s)
    real(r8), pointer :: gross_nmin_vr                         (:,:)   => null() ! vertically-resolved gross rate of N mineralization (gN/m3/s)
    real(r8), pointer :: gross_nmin                            (:)     => null() ! vert-int (diagnostic) gross rate of N mineralization (gN/m2/s)
    real(r8), pointer :: net_nmin_vr                           (:,:)   => null() ! vertically-resolved net rate of N mineralization (gN/m3/s)
    real(r8), pointer :: net_nmin                              (:)     => null() ! vert-int (diagnostic) net rate of N mineralization (gN/m2/s)
    real(r8), pointer :: sminn_no3_input_vr                    (:,:)   => null() ! no3 input, gN/m3/time step
    real(r8), pointer :: sminn_nh4_input_vr                    (:,:)   => null() ! nh4 input, gN/m3/time step
    real(r8), pointer :: sminn_no3_input                       (:)     => null() ! no3 input, gN/m2
    real(r8), pointer :: sminn_nh4_input                       (:)     => null() ! nh4 input, gN/m2
    real(r8), pointer :: sminn_input                           (:)     => null() ! minn input, gN/m2
    real(r8), pointer :: bgc_npool_ext_inputs_vr               (:,:,:) => null() ! organic nitrogen input, gN/m3/time step
    real(r8), pointer :: bgc_npool_ext_loss_vr                 (:,:,:) => null() ! extneral organic nitrogen loss, gN/m3/time step
    real(r8), pointer :: bgc_npool_inputs                      (:,:)   => null() ! organic N input, gN/m2/time step
    ! nitrification / denitrification fluxes
    real(r8), pointer :: f_nit_vr                              (:,:)   => null() ! (gN/m3/s) soil nitrification flux
    real(r8), pointer :: f_denit_vr                            (:,:)   => null() ! (gN/m3/s) soil denitrification flux
    real(r8), pointer :: f_nit                                 (:)     => null() ! (gN/m2/s) soil nitrification flux
    real(r8), pointer :: f_denit                               (:)     => null() ! (gN/m2/s) soil denitrification flux
    real(r8), pointer :: pot_f_nit_vr                          (:,:)   => null() ! (gN/m3/s) potential soil nitrification flux
    real(r8), pointer :: pot_f_denit_vr                        (:,:)   => null() ! (gN/m3/s) potential soil denitrification flux
    real(r8), pointer :: pot_f_nit                             (:)     => null() ! (gN/m2/s) potential soil nitrification flux
    real(r8), pointer :: pot_f_denit                           (:)     => null() ! (gN/m2/s) potential soil denitrification flux
    real(r8), pointer :: n2_n2o_ratio_denit_vr                 (:,:)   => null() ! ratio of N2 to N2O production by denitrification [gN/gN]
    real(r8), pointer :: f_n2o_denit_vr                        (:,:)   => null() ! flux of N2o from denitrification [gN/m^3/s]
    real(r8), pointer :: f_n2o_denit                           (:)     => null() ! flux of N2o from denitrification [gN/m^2/s]
    real(r8), pointer :: f_n2o_nit_vr                          (:,:)   => null() ! flux of N2o from nitrification [gN/m^3/s]
    real(r8), pointer :: f_n2o_nit                             (:)     => null() ! flux of N2o from nitrification [gN/m^2/s]
    ! immobilization / uptake fluxes
    real(r8), pointer :: actual_immob_no3_vr                   (:,:)   => null() ! vertically-resolved actual immobilization of NO3 (gN/m3/s)
    real(r8), pointer :: actual_immob_nh4_vr                   (:,:)   => null() ! vertically-resolved actual immobilization of NH4 (gN/m3/s)
    real(r8), pointer :: smin_no3_to_plant_vr                  (:,:)   => null() ! vertically-resolved plant uptake of soil NO3 (gN/m3/s)
    real(r8), pointer :: smin_nh4_to_plant_vr                  (:,:)   => null() ! vertically-resolved plant uptake of soil NH4 (gN/m3/s)
    real(r8), pointer :: actual_immob_no3                      (:)     => null() ! actual immobilization of NO3 (gN/m2/s)
    real(r8), pointer :: actual_immob_nh4                      (:)     => null() ! actual immobilization of NH4 (gN/m2/s)
    real(r8), pointer :: smin_no3_to_plant                     (:)     => null() ! plant uptake of soil NO3 (gN/m2/s)
    real(r8), pointer :: smin_nh4_to_plant                     (:)     => null() ! plant uptake of soil Nh4 (gN/m2/s)
    ! leaching fluxes
    real(r8), pointer :: smin_no3_leached_vr                   (:,:)   => null() ! vertically-resolved soil mineral NO3 loss to leaching (gN/m3/s)
    real(r8), pointer :: smin_no3_leached                      (:)     => null() ! soil mineral NO3 pool loss to leaching (gN/m2/s)
    real(r8), pointer :: smin_no3_runoff_vr                    (:,:)   => null() ! vertically-resolved rate of mineral NO3 loss with runoff (gN/m3/s)
    real(r8), pointer :: smin_no3_runoff                       (:)     => null() ! soil mineral NO3 pool loss to runoff (gN/m2/s)
    ! nitrification /denitrification diagnostic quantities
    real(r8), pointer :: smin_no3_massdens_vr                  (:,:)   => null() ! (ugN / g soil) soil nitrate concentration
    real(r8), pointer :: soil_bulkdensity                      (:,:)   => null() ! (kg soil / m3) bulk density of soil
    real(r8), pointer :: k_nitr_t_vr                           (:,:)   => null() 
    real(r8), pointer :: k_nitr_ph_vr                          (:,:)   => null() 
    real(r8), pointer :: k_nitr_h2o_vr                         (:,:)   => null() 
    real(r8), pointer :: k_nitr_vr                             (:,:)   => null() 
    real(r8), pointer :: wfps_vr                               (:,:)   => null() 
    real(r8), pointer :: fmax_denit_carbonsubstrate_vr         (:,:)   => null() 
    real(r8), pointer :: fmax_denit_nitrate_vr                 (:,:)   => null() 
    real(r8), pointer :: f_denit_base_vr                       (:,:)   => null() ! nitrification and denitrification fluxes
    real(r8), pointer :: diffus                                (:,:)   => null() ! diffusivity (m2/s)
    real(r8), pointer :: ratio_k1                              (:,:)   => null() 
    real(r8), pointer :: ratio_no3_co2                         (:,:)   => null() 
    real(r8), pointer :: soil_co2_prod                         (:,:)   => null() 
    real(r8), pointer :: fr_WFPS                               (:,:)   => null() 
    real(r8), pointer :: r_psi                                 (:,:)   => null() 
    real(r8), pointer :: anaerobic_frac                        (:,:)   => null() 
    ! (old version) denitrification fluxes
    real(r8), pointer :: sminn_to_denit_decomp_cascade_vr      (:,:,:) => null() ! vertically-resolved denitrification along decomp cascade (gN/m3/s) 
    real(r8), pointer :: sminn_to_denit_decomp_cascade         (:,:)   => null() ! vertically-integrated (diagnostic) denitrification along decomp cascade (gN/m2/s) 
    real(r8), pointer :: sminn_to_denit_excess_vr              (:,:)   => null() ! vertically-resolved denitrification from excess mineral N pool (gN/m3/s)
    real(r8), pointer :: sminn_to_denit_excess                 (:)     => null() ! vertically-integrated (diagnostic) denitrification from excess mineral N pool (gN/m2/s)
    ! (old version) leaching fluxes
    real(r8), pointer :: sminn_leached_vr                      (:,:)   => null() ! vertically-resolved soil mineral N pool loss to leaching (gN/m3/s)
    real(r8), pointer :: sminn_leached                         (:)     => null() ! soil mineral N pool loss to leaching (gN/m2/s)
    ! dynamic landcover fluxes
    real(r8), pointer :: dwt_slash_nflux                       (:)     => null() ! (gN/m2/s) conversion slash flux due to landcover change
    real(r8), pointer :: dwt_conv_nflux                        (:)     => null() ! (gN/m2/s) conversion N flux (immediate loss to atm)
    real(r8), pointer :: dwt_prod10n_gain                      (:)     => null() ! (gN/m2/s) addition to 10-yr wood product pool
    real(r8), pointer :: dwt_prod100n_gain                     (:)     => null() ! (gN/m2/s) addition to 100-yr wood product pool
    real(r8), pointer :: dwt_frootn_to_litr_met_n              (:,:)   => null() ! (gN/m3/s) fine root to litter due to landcover change
    real(r8), pointer :: dwt_frootn_to_litr_cel_n              (:,:)   => null() ! (gN/m3/s) fine root to litter due to landcover change
    real(r8), pointer :: dwt_frootn_to_litr_lig_n              (:,:)   => null() ! (gN/m3/s) fine root to litter due to landcover change
    real(r8), pointer :: dwt_livecrootn_to_cwdn                (:,:)   => null() ! (gN/m3/s) live coarse root to CWD due to landcover change
    real(r8), pointer :: dwt_deadcrootn_to_cwdn                (:,:)   => null() ! (gN/m3/s) dead coarse root to CWD due to landcover change
    real(r8), pointer :: dwt_nloss                             (:)     => null() ! (gN/m2/s) total nitrogen loss from product pools and conversion
    ! wood product pool loss fluxes
    real(r8), pointer :: prod1n_loss                           (:)     => null() ! (gN/m2/s) decomposition loss from 1-yr crop product pool
    real(r8), pointer :: prod10n_loss                          (:)     => null() ! (gN/m2/s) decomposition loss from 10-yr wood product pool
    real(r8), pointer :: prod100n_loss                         (:)     => null() ! (gN/m2/s) decomposition loss from 100-yr wood product pool
    real(r8), pointer :: product_nloss                         (:)     => null() ! (gN/m2/s) total wood product nitrogen loss
    ! summary (diagnostic) flux variables, not involved in mass balance
    real(r8), pointer :: denit                                 (:)     => null() ! total rate of denitrification (gN/m2/s)
    real(r8), pointer :: ninputs                               (:)     => null() ! column-level N inputs (gN/m2/s)
    real(r8), pointer :: noutputs                              (:)     => null() ! column-level N outputs (gN/m2/s)
    real(r8), pointer :: som_n_leached                         (:)     => null() ! total SOM N loss from vertical transport (gN/m^2/s)
    real(r8), pointer :: decomp_npools_leached                 (:,:)   => null() ! N loss from vertical transport from each decomposing N pool (gN/m^2/s)
    real(r8), pointer :: decomp_npools_transport_tendency      (:,:,:) => null() ! N tendency due to vertical transport in decomposing N pools (gN/m^3/s)
    ! all n pools involved in decomposition
    real(r8), pointer :: decomp_npools_sourcesink              (:,:,:) => null() ! (gN/m3) change in decomposing n pools 
    ! bgc interface/pflotran
    real(r8), pointer :: plant_ndemand                         (:)     => null() ! N flux required to support initial GPP (gN/m2/s)
    real(r8), pointer :: plant_ndemand_vr                      (:,:)   => null() ! vertically-resolved N flux required to support initial GPP (gN/m3/s)
    real(r8), pointer :: f_ngas_decomp_vr                      (:,:)   => null() ! vertically-resolved N emission from excess mineral N pool due to mineralization (gN/m3/s)
    real(r8), pointer :: f_ngas_decomp                         (:)     => null() ! N emission from excess mineral N pool due to mineralization (gN/m2/s)
    real(r8), pointer :: f_ngas_nitri_vr                       (:,:)   => null() ! vertically-resolved N emission from nitrification (gN/m3/s)
    real(r8), pointer :: f_ngas_nitri                          (:)     => null() ! vertically-resolved N emission from nitrification (gN/m2/s)
    real(r8), pointer :: f_ngas_denit_vr                       (:,:)   => null() ! vertically-resolved N emission from denitrification (gN/m3/s)
    real(r8), pointer :: f_ngas_denit                          (:)     => null() ! vertically-resolved N emission from denitrification (gN/m2/s)
    real(r8), pointer :: f_n2o_soil_vr                         (:,:)   => null() ! flux of N2o from soil-N processes [gN/m^3/s]
    real(r8), pointer :: f_n2o_soil                            (:)     => null() ! flux of N2o from soil-N processes [gN/m^2/s]
    real(r8), pointer :: f_n2_soil_vr                          (:,:)   => null() ! flux of N2 from soil-N processes [gN/m^3/s]
    real(r8), pointer :: f_n2_soil                             (:)     => null() ! flux of N2 from soil-N processes [gN/m^2/s]
    real(r8), pointer :: externaln_to_decomp_npools            (:,:,:) => null() ! net N fluxes associated with litter/som-adding/removal to decomp pools (gN/m3/s)
    real(r8), pointer :: externaln_to_decomp_delta             (:)     => null() ! summarized net N i/o changes associated with litter/som-adding/removal to decomp pools  btw time-step (gN/m2)
    real(r8), pointer :: no3_net_transport_vr                  (:,:)   => null() ! net NO3 transport associated with runoff/leaching/diffusion (gN/m3/s)
    real(r8), pointer :: nh4_net_transport_vr                  (:,:)   => null() ! net NH4 transport associated with runoff/leaching/diffusion (gN/m3/s)
    real(r8), pointer :: col_plant_ndemand_vr                  (:,:)   => null() ! plant N demand
    real(r8), pointer :: col_plant_nh4demand_vr                (:,:)   => null() ! plant NH4 demand
    real(r8), pointer :: col_plant_no3demand_vr                (:,:)   => null() ! plant NO3 demand
    real(r8), pointer :: col_plant_pdemand_vr                  (:,:)   => null() ! plant P demand
    real(r8), pointer :: pmnf_decomp_cascade                   (:,:,:) => null() ! potential mineral N flux, from one pool to another
    real(r8), pointer :: plant_n_uptake_flux                   (:)     => null() ! for the purpose of mass balance check  
    real(r8), pointer :: soil_n_immob_flux                     (:)     => null() ! for the purpose of mass balance check
    real(r8), pointer :: soil_n_immob_flux_vr                  (:,:)   => null() ! for the purpose of mass balance check
    real(r8), pointer :: soil_n_grossmin_flux                  (:)     => null() ! for the purpose of mass balance check
    real(r8), pointer :: plant_to_litter_nflux                 (:)     => null() ! for the purpose of mass balance check
    real(r8), pointer :: plant_to_cwd_nflux                    (:)     => null() ! for the purpose of mass balance check
  contains
    procedure, public :: Init       => col_nf_init
    procedure, public :: Clean      => col_nf_clean
  end type column_nitrogen_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds phosphorus flux information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_phosphorus_flux
    real(r8), pointer :: hrv_deadstemp_to_prod10p              (:)     ! dead stem P harvest mortality to 10-year product pool (gP/m2/s)
    real(r8), pointer :: hrv_deadstemp_to_prod100p             (:)     ! dead stem P harvest mortality to 100-year product pool (gP/m2/s)
    real(r8), pointer :: m_p_to_litr_met_fire                  (:,:)   ! P from leaf, froot, xfer and storage P to litter labile P by fire (gP/m3/s) 
    real(r8), pointer :: m_p_to_litr_cel_fire                  (:,:)   ! P from leaf, froot, xfer and storage P to litter cellulose P by fire (gP/m3/s) 
    real(r8), pointer :: m_p_to_litr_lig_fire                  (:,:)   ! P from leaf, froot, xfer and storage P to litter lignin P by fire (gP/m3/s) 
    real(r8), pointer :: harvest_p_to_litr_met_p               (:,:)   ! P fluxes associated with harvest to litter metabolic pool (gP/m3/s)
    real(r8), pointer :: harvest_p_to_litr_cel_p               (:,:)   ! P fluxes associated with harvest to litter cellulose pool (gP/m3/s)
    real(r8), pointer :: harvest_p_to_litr_lig_p               (:,:)   ! P fluxes associated with harvest to litter lignin pool (gP/m3/s)
    real(r8), pointer :: harvest_p_to_cwdp                     (:,:)   ! P fluxes associated with harvest to CWD pool (gP/m3/s)
    real(r8), pointer :: hrv_cropp_to_prod1p                   (:)     ! crop P harvest mortality to 1-yr product pool (gP/m2/s)
    real(r8), pointer :: m_decomp_ppools_to_fire_vr            (:,:,:) ! vertically-resolved decomposing P fire loss (gP/m3/s)
    real(r8), pointer :: m_decomp_ppools_to_fire               (:,:)   ! vertically-integrated (diagnostic) decomposing P fire loss (gP/m2/s)
    real(r8), pointer :: fire_ploss                            (:)     ! total column-level fire P loss (gP/m2/s)
    real(r8), pointer :: fire_decomp_ploss                     (:)     ! fire p loss from decomposable pools (gP/m2/s)
    real(r8), pointer :: fire_ploss_p2c                        (:)     ! patch2col column-level fire P loss (gP/m2/s) (p2c)
    real(r8), pointer :: fire_mortality_p_to_cwdp              (:,:)   ! P fluxes associated with fire mortality to CWD pool (gP/m3/s)
    real(r8), pointer :: wood_harvestp                         (:)     ! total P losses to wood product pools (gP/m2/s) (p2c)
    real(r8), pointer :: pdep_to_sminp                         (:)     ! atmospheric P deposition to soil mineral P (gP/m2/s)
    real(r8), pointer :: fert_p_to_sminp                       (:)     ! fertilizer P to soil mineral P (gP/m2/s)
    real(r8), pointer :: phenology_p_to_litr_met_p             (:,:)   ! P fluxes associated with phenology (litterfall and crop) to litter metabolic pool (gP/m3/s)
    real(r8), pointer :: phenology_p_to_litr_cel_p             (:,:)   ! P fluxes associated with phenology (litterfall and crop) to litter cellulose pool (gP/m3/s)
    real(r8), pointer :: phenology_p_to_litr_lig_p             (:,:)   ! P fluxes associated with phenology (litterfall and crop) to litter lignin pool (gP/m3/s)
    real(r8), pointer :: gap_mortality_p_to_litr_met_p         (:,:)   ! P fluxes associated with gap mortality to litter metabolic pool (gP/m3/s)
    real(r8), pointer :: gap_mortality_p_to_litr_cel_p         (:,:)   ! P fluxes associated with gap mortality to litter cellulose pool (gP/m3/s)
    real(r8), pointer :: gap_mortality_p_to_litr_lig_p         (:,:)   ! P fluxes associated with gap mortality to litter lignin pool (gP/m3/s)
    real(r8), pointer :: gap_mortality_p_to_cwdp               (:,:)   ! P fluxes associated with gap mortality to CWD pool (gP/m3/s)
    real(r8), pointer :: decomp_cascade_ptransfer_vr           (:,:,:) ! vert-res transfer of P from donor to receiver pool along decomp. cascade (gP/m3/s)
    real(r8), pointer :: decomp_cascade_ptransfer              (:,:)   ! vert-int (diagnostic) transfer of P from donor to receiver pool along decomp. cascade (gP/m2/s)
    real(r8), pointer :: decomp_cascade_sminp_flux_vr          (:,:,:) ! vert-res mineral P flux for transition along decomposition cascade (gP/m3/s)
    real(r8), pointer :: decomp_cascade_sminp_flux             (:,:)   ! vert-int (diagnostic) mineral P flux for transition along decomposition cascade (gP/m2/s)
    real(r8), pointer :: potential_immob_p_vr                  (:,:)   ! vertically-resolved potential P immobilization (gP/m3/s) at each level
    real(r8), pointer :: potential_immob_p                     (:)     ! vert-int (diagnostic) potential P immobilization (gP/m2/s)
    real(r8), pointer :: actual_immob_p_vr                     (:,:)   ! vertically-resolved actual P immobilization (gP/m3/s) at each level
    real(r8), pointer :: actual_immob_p                        (:)     ! vert-int (diagnostic) actual P immobilization (gP/m2/s)
    real(r8), pointer :: sminp_to_plant_vr                     (:,:)   ! vertically-resolved plant uptake of soil mineral P (gP/m3/s)
    real(r8), pointer :: sminp_to_plant                        (:)     ! vert-int (diagnostic) plant uptake of soil mineral P (gP/m2/s)
    real(r8), pointer :: supplement_to_sminp_vr                (:,:)   ! vertically-resolved supplemental P supply (gP/m3/s)
    real(r8), pointer :: supplement_to_sminp                   (:)     ! vert-int (diagnostic) supplemental P supply (gP/m2/s)
    real(r8), pointer :: gross_pmin_vr                         (:,:)   ! vertically-resolved gross rate of P mineralization (gP/m3/s)
    real(r8), pointer :: gross_pmin                            (:)     ! vert-int (diagnostic) gross rate of P mineralization (gP/m2/s)
    real(r8), pointer :: net_pmin_vr                           (:,:)   ! vertically-resolved net rate of P mineralization (gP/m3/s)
    real(r8), pointer :: net_pmin                              (:)     ! vert-int (diagnostic) net rate of P mineralization (gP/m2/s)
    real(r8), pointer :: biochem_pmin_ppools_vr                (:,:,:) ! vertically-resolved biochemical P mineralization for each soi pool (gP/m3/s)
    real(r8), pointer :: biochem_pmin_vr                       (:,:)   ! vertically-resolved total biochemical P mineralization (gP/m3/s)
    real(r8), pointer :: biochem_pmin_to_ecosysp_vr            (:,:)   ! biochemical P mineralization directly goes to soil (gP/m3/s)
    real(r8), pointer :: biochem_pmin                          (:)     ! vert-int (diagnostic) total biochemical P mineralization (gP/m3/s)
    real(r8), pointer :: primp_to_labilep_vr                   (:,:)   ! (gP/m3/s) flux of P from primary mineral to labile 
    real(r8), pointer :: primp_to_labilep                      (:)     ! (gP/m3/s) flux of P from primary mineral to labile 
    real(r8), pointer :: labilep_to_secondp_vr                 (:,:)   ! (gP/m3/s) flux of labile P to secondary mineral P 
    real(r8), pointer :: labilep_to_secondp                    (:)     ! (gP/m3/s) flux of labile P to secondary mineral P 
    real(r8), pointer :: secondp_to_labilep_vr                 (:,:)   ! (gP/m3/s) flux of the desorption of secondary mineral P to labile P
    real(r8), pointer :: secondp_to_labilep                    (:)     ! (gP/m3/s) flux of the desorption of secondary mineral P to labile P
    real(r8), pointer :: secondp_to_occlp_vr                   (:,:)   ! (gP/m3/s) flux of the occlusion of secondary P to occluded P
    real(r8), pointer :: secondp_to_occlp                      (:)     ! (gP/m3/s) flux of the occlusion of secondary P to occluded P
    real(r8), pointer :: sminp_leached_vr                      (:,:)   ! vertically-resolved soil mineral P pool loss to leaching (gP/m3/s)
    real(r8), pointer :: sminp_leached                         (:)     ! soil mineral P pool loss to leaching (gP/m2/s)
    real(r8), pointer :: dwt_slash_pflux                       (:)     ! (gP/m2/s) conversion slash flux due to landcover change
    real(r8), pointer :: dwt_conv_pflux                        (:)     ! (gP/m2/s) conversion P flux (immediate loss to atm)
    real(r8), pointer :: dwt_prod10p_gain                      (:)     ! (gP/m2/s) addition to 10-yr wood product pool
    real(r8), pointer :: dwt_prod100p_gain                     (:)     ! (gP/m2/s) addition to 100-yr wood product pool
    real(r8), pointer :: dwt_frootp_to_litr_met_p              (:,:)   ! (gP/m3/s) fine root to litter due to landcover change
    real(r8), pointer :: dwt_frootp_to_litr_cel_p              (:,:)   ! (gP/m3/s) fine root to litter due to landcover change
    real(r8), pointer :: dwt_frootp_to_litr_lig_p              (:,:)   ! (gP/m3/s) fine root to litter due to landcover change
    real(r8), pointer :: dwt_livecrootp_to_cwdp                (:,:)   ! (gP/m3/s) live coarse root to CWD due to landcover change
    real(r8), pointer :: dwt_deadcrootp_to_cwdp                (:,:)   ! (gP/m3/s) dead coarse root to CWD due to landcover change
    real(r8), pointer :: dwt_ploss                             (:)     ! (gP/m2/s) total phosphorus loss from product pools and conversion
    real(r8), pointer :: prod1p_loss                           (:)     ! (gP/m2/s) decomposition loss from 1-yr crop product pool
    real(r8), pointer :: prod10p_loss                          (:)     ! (gP/m2/s) decomposition loss from 10-yr wood product pool
    real(r8), pointer :: prod100p_loss                         (:)     ! (gP/m2/s) decomposition loss from 100-yr wood product pool
    real(r8), pointer :: product_ploss                         (:)     ! (gP/m2/s) total wood product phosphorus loss
    real(r8), pointer :: pinputs                               (:)     ! column-level P inputs (gP/m2/s)
    real(r8), pointer :: poutputs                              (:)     ! column-level P outputs (gP/m2/s)
    real(r8), pointer :: som_p_leached                         (:)     ! total SOM P loss from vertical transport (gP/m^2/s)
    real(r8), pointer :: decomp_ppools_leached                 (:,:)   ! P loss from vertical transport from each decomposing P pool (gP/m^2/s)
    real(r8), pointer :: decomp_ppools_transport_tendency      (:,:,:) ! P tendency due to vertical transport in decomposing P pools (gP/m^3/s)
    real(r8), pointer :: decomp_ppools_sourcesink              (:,:,:) ! (gP/m3) change in decomposing P pools
    real(r8), pointer :: plant_pdemand                         (:)     ! P flux required to support initial GPP (gN/m2/s)
    real(r8), pointer :: plant_pdemand_vr                      (:,:)   ! vertically-resolved P flux required to support initial GPP (gP/m3/s)
    real(r8), pointer :: externalp_to_decomp_ppools            (:,:,:) ! net N fluxes associated with litter/som-adding/removal to decomp pools (gP/m3/s)
    real(r8), pointer :: externalp_to_decomp_delta             (:)     ! summarized net N i/o changes associated with litter/som-adding/removal to decomp pools  btw time-step (gP/m2)
    real(r8), pointer :: sminp_net_transport_vr                (:,:)   ! net sminp transport associated with runoff/leaching (gP/m3/s)
    real(r8), pointer :: sminp_net_transport_delta             (:)     ! summarized net change of column-level sminp leaching bwtn time-step (for balance checking) (gP/m2)
    real(r8), pointer :: adsorb_to_labilep_vr                  (:,:)
    real(r8), pointer :: desorb_to_solutionp_vr                (:,:)
    real(r8), pointer :: adsorb_to_labilep                     (:)
    real(r8), pointer :: desorb_to_solutionp                   (:)
    real(r8), pointer :: pmpf_decomp_cascade                   (:,:,:)
    real(r8), pointer :: plant_p_uptake_flux                   (:)     ! for the purpose of mass balance check  
    real(r8), pointer :: soil_p_immob_flux                     (:)     ! for the purpose of mass balance check
    real(r8), pointer :: soil_p_immob_flux_vr                  (:,:)   ! for the purpose of mass balance check
    real(r8), pointer :: soil_p_grossmin_flux                  (:)     ! for the purpose of mass balance check
    real(r8), pointer :: smin_p_to_plant                       (:)     ! for the purpose of mass balance check
    real(r8), pointer :: plant_to_litter_pflux                 (:)     ! for the purpose of mass balance check
    real(r8), pointer :: plant_to_cwd_pflux                    (:)     ! for the purpose of mass balance check
  contains
    procedure, public :: Init       => col_pf_init
    procedure, public :: Clean      => col_pf_clean
  end type column_phosphorus_flux
   
  !-----------------------------------------------------------------------
  ! declare the public instances of column-level data types
  !-----------------------------------------------------------------------
  ! State types 
  type(column_energy_state)          , public, target :: col_es     ! column energy state
  type(column_water_state)           , public, target :: col_ws     ! column water state
  type(column_carbon_state)          , public, target :: col_cs     ! column carbon state
  type(column_carbon_state)          , public, target :: c13_col_cs ! column carbon state (C13)
  type(column_carbon_state)          , public, target :: c14_col_cs ! column carbon state (C14)
  type(column_nitrogen_state)        , public, target :: col_ns     ! column nitrogen state
  type(column_phosphorus_state)      , public, target :: col_ps     ! column phosphorus state
  ! Flux types 
  type(column_energy_flux)           , public, target :: col_ef     ! column energy flux
  type(column_water_flux)            , public, target :: col_wf     ! column water flux
  type(column_carbon_flux)           , public, target :: col_cf     ! column carbon flux
  type(column_carbon_flux)           , public, target :: c13_col_cf ! column carbon flux
  type(column_carbon_flux)           , public, target :: c14_col_cf ! column carbon flux
  type(column_nitrogen_flux)         , public, target :: col_nf     ! column nitrogen flux
  type(column_phosphorus_flux)       , public, target :: col_pf     ! column phosphorus flux

  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column energy state data structure
  !------------------------------------------------------------------------
  subroutine col_es_init(this, begc, endc)
    !
    ! !USES:
    use landunit_varcon, only : istice, istwet, istsoil, istdlak, istice_mec
    use clm_varctl     , only : iulog, use_cn, use_vancouver, use_mexicocity
    use column_varcon  , only : icol_road_perv, icol_road_imperv, icol_roof, icol_sunwall, icol_shadewall
    !
    ! !ARGUMENTS:
    class(column_energy_state) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    !
    ! !LOCAL VARIABLES:
    integer           :: c,l,j                        ! indices
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays

    !-----------------------------------------------------------------------
    ! allocate for each member of col_es
    !-----------------------------------------------------------------------
    allocate(this%t_soisno         (begc:endc,-nlevsno+1:nlevgrnd)) ; this%t_soisno           (:,:) = nan
    allocate(this%t_ssbef          (begc:endc,-nlevsno+1:nlevgrnd)) ; this%t_ssbef            (:,:) = nan
    allocate(this%t_h2osfc         (begc:endc))                     ; this%t_h2osfc           (:)   = nan
    allocate(this%t_h2osfc_bef     (begc:endc))                     ; this%t_h2osfc_bef       (:)   = nan
    allocate(this%t_soi10cm        (begc:endc))                     ; this%t_soi10cm          (:)   = nan
    allocate(this%t_soi17cm        (begc:endc))                     ; this%t_soi17cm          (:)   = spval
    allocate(this%t_grnd           (begc:endc))                     ; this%t_grnd             (:)   = nan
    allocate(this%t_lake           (begc:endc,1:nlevlak))           ; this%t_lake             (:,:) = nan
    allocate(this%t_grnd_r         (begc:endc))                     ; this%t_grnd_r           (:)   = nan
    allocate(this%t_grnd_u         (begc:endc))                     ; this%t_grnd_u           (:)   = nan
    allocate(this%snot_top         (begc:endc))                     ; this%snot_top           (:)   = nan
    allocate(this%dTdz_top         (begc:endc))                     ; this%dTdz_top           (:)   = nan
    allocate(this%thv              (begc:endc))                     ; this%thv                (:)   = nan
    allocate(this%hc_soi           (begc:endc))                     ; this%hc_soi             (:)   = nan
    allocate(this%hc_soisno        (begc:endc))                     ; this%hc_soisno          (:)   = nan
    allocate(this%emg              (begc:endc))                     ; this%emg                (:)   = nan
    allocate(this%fact             (begc:endc, -nlevsno+1:nlevgrnd)); this%fact               (:,:) = nan
    allocate(this%c_h2osfc         (begc:endc))                     ; this%c_h2osfc           (:)   = nan
    allocate(this%t_nearsurf       (begc:endc))                     ; this%t_nearsurf         (:)   = nan

  end subroutine col_es_init

  !------------------------------------------------------------------------
  subroutine col_es_clean(this)
    !
    ! !ARGUMENTS:
    class(column_energy_state) :: this
    !------------------------------------------------------------------------
    deallocate(this%t_h2osfc)
  end subroutine col_es_clean
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column water state data structure
  !------------------------------------------------------------------------
  subroutine col_ws_init(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_water_state) :: this
    integer , intent(in)      :: begc,endc
    !
    ! !LOCAL VARIABLES:
    real(r8), pointer  :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    real(r8)           :: snowbd      ! temporary calculation of snow bulk density (kg/m3)
    real(r8)           :: fmelt       ! snowbd/100
    integer            :: c,l,j,nlevs,nlevbed, ncells
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of col_ws
    !-----------------------------------------------------------------------
    allocate(this%h2osoi_liq         (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_liq         (:,:) = nan
    allocate(this%h2osoi_ice         (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_ice         (:,:) = nan
    allocate(this%h2osoi_vol         (begc:endc, 1:nlevgrnd))         ; this%h2osoi_vol         (:,:) = nan
    allocate(this%h2osfc             (begc:endc))                     ; this%h2osfc             (:)   = nan   
    allocate(this%h2ocan             (begc:endc))                     ; this%h2ocan             (:)   = nan  
    allocate(this%total_plant_stored_h2o(begc:endc))                  ; this%total_plant_stored_h2o(:)= nan  
    allocate(this%h2osoi_liqvol      (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_liqvol      (:,:) = nan
    allocate(this%h2osoi_icevol      (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_icevol      (:,:) = nan    
    allocate(this%h2osoi_liq_old     (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_liq_old     (:,:) = nan
    allocate(this%h2osoi_ice_old     (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_ice_old     (:,:) = nan 
    allocate(this%bw                 (begc:endc,-nlevsno+1:0))        ; this%bw                 (:,:) = nan   
    allocate(this%smp_l              (begc:endc,-nlevsno+1:nlevgrnd)) ; this%smp_l              (:,:) = nan
    allocate(this%soilp              (begc:endc,1:nlevgrnd))          ; this%soilp              (:,:) = 0._r8
    allocate(this%swe_old            (begc:endc,-nlevsno+1:0))        ; this%swe_old            (:,:) = nan   
    allocate(this%snw_rds            (begc:endc,-nlevsno+1:0))        ; this%snw_rds            (:,:) = nan
    allocate(this%air_vol            (begc:endc, 1:nlevgrnd))         ; this%air_vol            (:,:) = nan
    allocate(this%qg_snow            (begc:endc))                     ; this%qg_snow            (:)   = nan   
    allocate(this%qg_soil            (begc:endc))                     ; this%qg_soil            (:)   = nan   
    allocate(this%qg_h2osfc          (begc:endc))                     ; this%qg_h2osfc          (:)   = nan   
    allocate(this%qg                 (begc:endc))                     ; this%qg                 (:)   = nan   
    allocate(this%dqgdT              (begc:endc))                     ; this%dqgdT              (:)   = nan   
    allocate(this%h2osoi_liqice_10cm (begc:endc))                     ; this%h2osoi_liqice_10cm (:)   = nan
    allocate(this%h2osno             (begc:endc))                     ; this%h2osno             (:)   = nan   
    allocate(this%h2osno_old         (begc:endc))                     ; this%h2osno_old         (:)   = nan   
    allocate(this%h2osno_top         (begc:endc))                     ; this%h2osno_top         (:)   = nan
    allocate(this%sno_liq_top        (begc:endc))                     ; this%sno_liq_top        (:)   = nan
    allocate(this%snowice            (begc:endc))                     ; this%snowice            (:)   = nan   
    allocate(this%snowliq            (begc:endc))                     ; this%snowliq            (:)   = nan   
    allocate(this%int_snow           (begc:endc))                     ; this%int_snow           (:)   = nan   
    allocate(this%snow_depth         (begc:endc))                     ; this%snow_depth         (:)   = nan
    allocate(this%snowdp             (begc:endc))                     ; this%snowdp             (:)   = nan
    allocate(this%snow_persistence   (begc:endc))                     ; this%snow_persistence   (:)   = nan
    allocate(this%snw_rds_top        (begc:endc))                     ; this%snw_rds_top        (:)   = nan
    allocate(this%do_capsnow         (begc:endc))                   
    allocate(this%frac_sno           (begc:endc))                     ; this%frac_sno           (:)   = nan
    allocate(this%frac_sno_eff       (begc:endc))                     ; this%frac_sno_eff       (:)   = nan
    allocate(this%frac_iceold        (begc:endc,-nlevsno+1:nlevgrnd)) ; this%frac_iceold        (:,:) = nan
    allocate(this%frac_h2osfc        (begc:endc))                     ; this%frac_h2osfc        (:)   = nan 
    allocate(this%wf                 (begc:endc))                     ; this%wf                 (:)   = nan
    allocate(this%wf2                (begc:endc))                     ; this%wf2                (:)   = nan
    allocate(this%finundated         (begc:endc))                     ; this%finundated         (:)   = nan
    allocate(this%begwb              (begc:endc))                     ; this%begwb              (:)   = nan
    allocate(this%endwb              (begc:endc))                     ; this%endwb              (:)   = nan
    allocate(this%errh2o             (begc:endc))                     ; this%errh2o             (:)   = nan
    allocate(this%errh2osno          (begc:endc))                     ; this%errh2osno          (:)   = nan
    allocate(this%h2osoi_liq_depth_intg(begc:endc))                   ; this%h2osoi_liq_depth_intg(:) = nan
    allocate(this%h2osoi_ice_depth_intg(begc:endc))                   ; this%h2osoi_ice_depth_intg(:) = nan
    ncells = (endc - begc + 1)*nlevgrnd
    allocate(this%vsfm_fliq_col_1d   (ncells))                        ; this%vsfm_fliq_col_1d   (:)   = nan
    allocate(this%vsfm_sat_col_1d    (ncells))                        ; this%vsfm_sat_col_1d    (:)   = nan
    allocate(this%vsfm_mass_col_1d   (ncells))                        ; this%vsfm_mass_col_1d   (:)   = nan
    allocate(this%vsfm_smpl_col_1d   (ncells))                        ; this%vsfm_smpl_col_1d   (:)   = nan
    allocate(this%vsfm_soilp_col_1d  (ncells))                        ; this%vsfm_soilp_col_1d  (:)   = nan
    
  end subroutine col_ws_init

  !------------------------------------------------------------------------
  subroutine col_ws_clean(this)
    !
    ! !ARGUMENTS:
    class(column_water_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine col_ws_clean
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column carbon state data structure
  !------------------------------------------------------------------------
  subroutine col_cs_init(this, begc, endc, carbon_type, ratio, c12_carbonstate_vars)
    !
    ! !ARGUMENTS:
    class(column_carbon_state)    :: this
    integer          , intent(in) :: begc,endc
    character(len=3) , intent(in) :: carbon_type ! one of ['c12', c13','c14']
    real(r8)         , intent(in) :: ratio
    type(column_carbon_state), optional, intent(in) :: c12_carbonstate_vars
    !
    ! !LOCAL VARIABLES:
    integer           :: c,l,j,k
    character(24)     :: fieldname
    character(100)    :: longname
    real(r8), pointer :: data1dptr(:)   ! temp. pointer for slicing larger arrays
    real(r8), pointer :: data2dptr(:,:) ! temp. pointer for slicing larger arrays
    integer :: fc                                        ! filter index
    integer :: num_special_col                           ! number of good values in special_col filter
    integer :: special_col(endc-begc+1)    ! special landunit filter - columns
    !------------------------------------------------------------------------
    
    !-----------------------------------------------------------------------
    ! allocate for each member of col_cs
    !-----------------------------------------------------------------------
    allocate(this%rootc                (begc:endc))     ; this%rootc                (:)     = nan
    allocate(this%totvegc              (begc:endc))     ; this%totvegc              (:)     = nan
    allocate(this%leafc                (begc:endc))     ; this%leafc                (:)     = nan
    allocate(this%deadstemc            (begc:endc))     ; this%deadstemc            (:)     = nan
    allocate(this%fuelc                (begc:endc))     ; this%fuelc                (:)     = nan
    allocate(this%fuelc_crop           (begc:endc))     ; this%fuelc_crop           (:)     = nan
    allocate(this%frootc               (begc:endc))     ; this%frootc               (:)     = nan
    allocate(this%seedc                (begc:endc))     ; this%seedc                (:)     = nan
    allocate(this%prod1c               (begc:endc))     ; this%prod1c               (:)     = nan
    allocate(this%prod10c              (begc:endc))     ; this%prod10c              (:)     = nan
    allocate(this%prod100c             (begc:endc))     ; this%prod100c             (:)     = nan
    allocate(this%totprodc             (begc:endc))     ; this%totprodc             (:)     = nan
    allocate(this%dyn_cbal_adjustments (begc:endc))     ; this%dyn_cbal_adjustments (:)     = nan
    allocate(this%totpftc              (begc:endc))     ; this%totpftc              (:)     = nan
    allocate(this%cwdc                 (begc:endc))     ; this%cwdc                 (:)     = nan
    allocate(this%ctrunc               (begc:endc))     ; this%ctrunc               (:)     = nan
    allocate(this%totlitc              (begc:endc))     ; this%totlitc              (:)     = nan
    allocate(this%totsomc              (begc:endc))     ; this%totsomc              (:)     = nan
    allocate(this%totlitc_1m           (begc:endc))     ; this%totlitc_1m           (:)     = nan
    allocate(this%totsomc_1m           (begc:endc))     ; this%totsomc_1m           (:)     = nan
    allocate(this%totecosysc           (begc:endc))     ; this%totecosysc           (:)     = nan
    allocate(this%totcolc              (begc:endc))     ; this%totcolc              (:)     = nan
    allocate(this%totabgc              (begc:endc))     ; this%totabgc              (:)     = nan
    allocate(this%totblgc              (begc:endc))     ; this%totblgc              (:)     = nan
    allocate(this%totvegc_abg          (begc:endc))     ; this%totvegc_abg          (:)     = nan
    allocate(this%begcb                (begc:endc))     ; this%begcb                (:)     = nan
    allocate(this%endcb                (begc:endc))     ; this%endcb                (:)     = nan
    allocate(this%errcb                (begc:endc))     ; this%errcb                (:)     = nan
    allocate(this%totpftc_beg          (begc:endc))     ; this%totpftc_beg          (:)     = nan
    allocate(this%cwdc_beg             (begc:endc))     ; this%cwdc_beg             (:)     = nan
    allocate(this%totlitc_beg          (begc:endc))     ; this%totlitc_beg          (:)     = nan
    allocate(this%totsomc_beg          (begc:endc))     ; this%totsomc_beg          (:)     = nan
    allocate(this%totpftc_end          (begc:endc))     ; this%totpftc_end          (:)     = nan
    allocate(this%cwdc_end             (begc:endc))     ; this%cwdc_end             (:)     = nan
    allocate(this%totlitc_end          (begc:endc))     ; this%totlitc_end          (:)     = nan
    allocate(this%totsomc_end          (begc:endc))     ; this%totsomc_end          (:)     = nan
    allocate(this%cropseedc_deficit    (begc:endc))     ; this%cropseedc_deficit    (:)     = nan
    allocate(this%decomp_cpools_vr (begc:endc,1:nlevdecomp_full,1:ndecomp_pools)) ; this%decomp_cpools_vr (:,:,:) = nan
    allocate(this%ctrunc_vr        (begc:endc,1:nlevdecomp_full))                 ; this%ctrunc_vr        (:,:)   = nan
    allocate(this%decomp_som2c_vr  (begc:endc,1:nlevdecomp_full))                 ; this%decomp_som2c_vr  (:,:)   = nan
    allocate(this%decomp_cpools_1m (begc:endc,1:ndecomp_pools))                   ; this%decomp_cpools_1m (:,:)   = nan
    allocate(this%decomp_cpools    (begc:endc,1:ndecomp_pools))                   ; this%decomp_cpools    (:,:)   = nan

          
  end subroutine col_cs_init
    

  !------------------------------------------------------------------------
  subroutine col_cs_clean(this)
    !
    ! !ARGUMENTS:
    class(column_carbon_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine col_cs_clean
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column nitrogen state data structure
  !------------------------------------------------------------------------
  subroutine col_ns_init(this, begc, endc, col_cs)
    !
    ! !ARGUMENTS:
    class(column_nitrogen_state)          :: this
    integer, intent(in)                   :: begc,endc
    type(column_carbon_state), intent(in) :: col_cs
    !
    ! !LOCAL VARIABLES:
    integer           :: c,l,j,k,fc
    integer           :: num_special_col           ! number of good values in special_col filter
    integer           :: special_col (endc-begc+1) ! special landunit filter - columns
    character(24)     :: fieldname
    character(100)    :: longname
    character(8)      :: vr_suffix
    real(r8), pointer :: data1dptr(:)   ! temp. pointer for slicing larger arrays
    real(r8), pointer :: data2dptr(:,:) ! temp. pointer for slicing larger arrays
    !------------------------------------------------------------------------
    
    !-----------------------------------------------------------------------
    ! allocate for each member of col_ns
    !-----------------------------------------------------------------------

    allocate(this%decomp_npools_vr      (begc:endc,1:nlevdecomp_full,1:ndecomp_pools)) ; this%decomp_npools_vr(:,:,:) = nan
    allocate(this%ntrunc_vr             (begc:endc,1:nlevdecomp_full))   ; this%ntrunc_vr             (:,:) = nan
    allocate(this%sminn_vr              (begc:endc,1:nlevdecomp_full))   ; this%sminn_vr              (:,:) = nan
    allocate(this%smin_no3_vr           (begc:endc,1:nlevdecomp_full))   ; this%smin_no3_vr           (:,:) = nan
    allocate(this%smin_nh4_vr           (begc:endc,1:nlevdecomp_full))   ; this%smin_nh4_vr           (:,:) = nan
    allocate(this%smin_nh4sorb_vr       (begc:endc,1:nlevdecomp_full))   ; this%smin_nh4sorb_vr       (:,:) = nan
    allocate(this%decomp_npools         (begc:endc,1:ndecomp_pools))     ; this%decomp_npools         (:,:) = nan
    allocate(this%decomp_npools_1m      (begc:endc,1:ndecomp_pools))     ; this%decomp_npools_1m      (:,:) = nan
    allocate(this%smin_no3              (begc:endc))                     ; this%smin_no3              (:)   = nan
    allocate(this%smin_nh4              (begc:endc))                     ; this%smin_nh4              (:)   = nan
    allocate(this%smin_nh4sorb          (begc:endc))                     ; this%smin_nh4sorb          (:)   = nan
    allocate(this%sminn                 (begc:endc))                     ; this%sminn                 (:)   = nan
    allocate(this%ntrunc                (begc:endc))                     ; this%ntrunc                (:)   = nan
    allocate(this%cwdn                  (begc:endc))                     ; this%cwdn                  (:)   = nan
    allocate(this%totlitn               (begc:endc))                     ; this%totlitn               (:)   = nan
    allocate(this%totsomn               (begc:endc))                     ; this%totsomn               (:)   = nan
    allocate(this%totlitn_1m            (begc:endc))                     ; this%totlitn_1m            (:)   = nan
    allocate(this%totsomn_1m            (begc:endc))                     ; this%totsomn_1m            (:)   = nan
    allocate(this%totecosysn            (begc:endc))                     ; this%totecosysn            (:)   = nan
    allocate(this%totcoln               (begc:endc))                     ; this%totcoln               (:)   = nan
    allocate(this%totabgn               (begc:endc))                     ; this%totabgn               (:)   = nan
    allocate(this%totblgn               (begc:endc))                     ; this%totblgn               (:)   = nan
    allocate(this%totvegn               (begc:endc))                     ; this%totvegn               (:)   = nan
    allocate(this%totpftn               (begc:endc))                     ; this%totpftn               (:)   = nan
    allocate(this%plant_n_buffer        (begc:endc))                     ; this%plant_n_buffer        (:)   = nan
    allocate(this%plant_nbuffer         (begc:endc))                     ; this%plant_nbuffer         (:)   = nan
    allocate(this%seedn                 (begc:endc))                     ; this%seedn                 (:)   = nan
    allocate(this%cropseedn_deficit     (begc:endc))                     ; this%cropseedn_deficit     (:)   = nan
    allocate(this%prod1n                (begc:endc))                     ; this%prod1n                (:)   = nan
    allocate(this%prod10n               (begc:endc))                     ; this%prod10n               (:)   = nan
    allocate(this%prod100n              (begc:endc))                     ; this%prod100n              (:)   = nan
    allocate(this%totprodn              (begc:endc))                     ; this%totprodn              (:)   = nan
    allocate(this%dyn_nbal_adjustments  (begc:endc))                     ; this%dyn_nbal_adjustments  (:)   = nan
    allocate(this%totpftn_beg           (begc:endc))                     ; this%totpftn_beg           (:)   = nan
    allocate(this%totpftn_end           (begc:endc))                     ; this%totpftn_end           (:)   = nan
    allocate(this%cwdn_beg              (begc:endc))                     ; this%cwdn_beg              (:)   = nan
    allocate(this%cwdn_end              (begc:endc))                     ; this%cwdn_end              (:)   = nan
    allocate(this%totlitn_beg           (begc:endc))                     ; this%totlitn_beg           (:)   = nan
    allocate(this%totlitn_end           (begc:endc))                     ; this%totlitn_end           (:)   = nan
    allocate(this%totsomn_beg           (begc:endc))                     ; this%totsomn_beg           (:)   = nan
    allocate(this%totsomn_end           (begc:endc))                     ; this%totsomn_end           (:)   = nan
    allocate(this%sminn_beg             (begc:endc))                     ; this%sminn_beg             (:)   = nan
    allocate(this%sminn_end             (begc:endc))                     ; this%sminn_end             (:)   = nan
    allocate(this%smin_no3_beg          (begc:endc))                     ; this%smin_no3_beg          (:)   = nan
    allocate(this%smin_no3_end          (begc:endc))                     ; this%smin_no3_end          (:)   = nan
    allocate(this%smin_nh4_beg          (begc:endc))                     ; this%smin_nh4_beg          (:)   = nan
    allocate(this%smin_nh4_end          (begc:endc))                     ; this%smin_nh4_end          (:)   = nan
    allocate(this%totprodn_beg          (begc:endc))                     ; this%totprodn_beg          (:)   = nan
    allocate(this%totprodn_end          (begc:endc))                     ; this%totprodn_end          (:)   = nan
    allocate(this%seedn_beg             (begc:endc))                     ; this%seedn_beg             (:)   = nan
    allocate(this%seedn_end             (begc:endc))                     ; this%seedn_end             (:)   = nan
    allocate(this%ntrunc_beg            (begc:endc))                     ; this%ntrunc_beg            (:)   = nan
    allocate(this%ntrunc_end            (begc:endc))                     ; this%ntrunc_end            (:)   = nan
    allocate(this%begnb                 (begc:endc))                     ; this%begnb                 (:)   = nan
    allocate(this%endnb                 (begc:endc))                     ; this%endnb                 (:)   = nan
    allocate(this%errnb                 (begc:endc))                     ; this%errnb                 (:)   = nan

  
  end subroutine col_ns_init
    

  !------------------------------------------------------------------------
  subroutine col_ns_clean(this)
    !
    ! !ARGUMENTS:
    class(column_nitrogen_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine col_ns_clean
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column phosphorus state data structure
  !------------------------------------------------------------------------
  subroutine col_ps_init(this, begc, endc, col_cs)
    !
    ! !ARGUMENTS:
    class(column_phosphorus_state)        :: this
    integer, intent(in)                   :: begc,endc
    type(column_carbon_state), intent(in) :: col_cs
    !
    ! !LOCAL VARIABLES:
    integer           :: fc,l,c,j,k               ! indices
    integer           :: num_special_col          ! number of good values in special_col filter
    integer           :: special_col(endc-begc+1) ! special landunit filter - columns
    character(8)      :: vr_suffix
    character(24)     :: fieldname
    character(100)    :: longname
    real(r8), pointer :: data1dptr(:)             ! temp. pointer for slicing larger arrays
    real(r8), pointer :: data2dptr(:,:)           ! temp. pointer for slicing larger arrays
    !------------------------------------------------------------------------
    
    !-----------------------------------------------------------------------
    ! allocate for each member of col_ps
    !-----------------------------------------------------------------------
    allocate(this%ptrunc_vr            (begc:endc,1:nlevdecomp_full)) ; this%ptrunc_vr            (:,:) = nan
    allocate(this%solutionp_vr         (begc:endc,1:nlevdecomp_full)) ; this%solutionp_vr         (:,:) = nan  
    allocate(this%labilep_vr           (begc:endc,1:nlevdecomp_full)) ; this%labilep_vr           (:,:) = nan  
    allocate(this%secondp_vr           (begc:endc,1:nlevdecomp_full)) ; this%secondp_vr           (:,:) = nan  
    allocate(this%occlp_vr             (begc:endc,1:nlevdecomp_full)) ; this%occlp_vr             (:,:) = nan  
    allocate(this%primp_vr             (begc:endc,1:nlevdecomp_full)) ; this%primp_vr             (:,:) = nan  
    allocate(this%sminp_vr             (begc:endc,1:nlevdecomp_full)) ; this%sminp_vr             (:,:) = nan  
    allocate(this%solutionp            (begc:endc))                   ; this%solutionp            (:)   = nan
    allocate(this%labilep              (begc:endc))                   ; this%labilep              (:)   = nan
    allocate(this%secondp              (begc:endc))                   ; this%secondp              (:)   = nan
    allocate(this%occlp                (begc:endc))                   ; this%occlp                (:)   = nan
    allocate(this%primp                (begc:endc))                   ; this%primp                (:)   = nan
    allocate(this%cwdp                 (begc:endc))                   ; this%cwdp                 (:)   = nan
    allocate(this%sminp                (begc:endc))                   ; this%sminp                (:)   = nan
    allocate(this%ptrunc               (begc:endc))                   ; this%ptrunc               (:)   = nan
    allocate(this%seedp                (begc:endc))                   ; this%seedp                (:)   = nan
    allocate(this%prod1p               (begc:endc))                   ; this%prod1p               (:)   = nan
    allocate(this%prod10p              (begc:endc))                   ; this%prod10p              (:)   = nan
    allocate(this%prod100p             (begc:endc))                   ; this%prod100p             (:)   = nan
    allocate(this%totprodp             (begc:endc))                   ; this%totprodp             (:)   = nan
    allocate(this%dyn_pbal_adjustments (begc:endc))                   ; this%dyn_pbal_adjustments (:)   = nan
    allocate(this%totlitp              (begc:endc))                   ; this%totlitp              (:)   = nan
    allocate(this%totsomp              (begc:endc))                   ; this%totsomp              (:)   = nan
    allocate(this%totlitp_1m           (begc:endc))                   ; this%totlitp_1m           (:)   = nan
    allocate(this%totsomp_1m           (begc:endc))                   ; this%totsomp_1m           (:)   = nan
    allocate(this%totecosysp           (begc:endc))                   ; this%totecosysp           (:)   = nan
    allocate(this%totcolp              (begc:endc))                   ; this%totcolp              (:)   = nan
    allocate(this%decomp_ppools        (begc:endc,1:ndecomp_pools))   ; this%decomp_ppools        (:,:) = nan
    allocate(this%decomp_ppools_1m     (begc:endc,1:ndecomp_pools))   ; this%decomp_ppools_1m     (:,:) = nan
    allocate(this%totpftp              (begc:endc))                   ; this%totpftp              (:)   = nan
    allocate(this%totvegp              (begc:endc))                   ; this%totvegp              (:)   = nan
    allocate(this%decomp_ppools_vr     (begc:endc,1:nlevdecomp_full,1:ndecomp_pools)); this%decomp_ppools_vr(:,:,:)= nan
    allocate(this%begpb                (begc:endc))                   ; this%begpb                (:)   = nan
    allocate(this%endpb                (begc:endc))                   ; this%endpb                (:)   = nan
    allocate(this%errpb                (begc:endc))                   ; this%errpb                (:)   = nan 
    allocate(this%solutionp_vr_cur     (begc:endc,1:nlevdecomp_full)) ; this%solutionp_vr_cur     (:,:) = nan
    allocate(this%solutionp_vr_prev    (begc:endc,1:nlevdecomp_full)) ; this%solutionp_vr_prev    (:,:) = nan
    allocate(this%labilep_vr_cur       (begc:endc,1:nlevdecomp_full)) ; this%labilep_vr_cur       (:,:) = nan
    allocate(this%labilep_vr_prev      (begc:endc,1:nlevdecomp_full)) ; this%labilep_vr_prev      (:,:) = nan
    allocate(this%secondp_vr_cur       (begc:endc,1:nlevdecomp_full)) ; this%secondp_vr_cur       (:,:) = nan
    allocate(this%secondp_vr_prev      (begc:endc,1:nlevdecomp_full)) ; this%secondp_vr_prev      (:,:) = nan
    allocate(this%occlp_vr_cur         (begc:endc,1:nlevdecomp_full)) ; this%occlp_vr_cur         (:,:) = nan
    allocate(this%occlp_vr_prev        (begc:endc,1:nlevdecomp_full)) ; this%occlp_vr_prev        (:,:) = nan
    allocate(this%primp_vr_cur         (begc:endc,1:nlevdecomp_full)) ; this%primp_vr_cur         (:,:) = nan
    allocate(this%primp_vr_prev        (begc:endc,1:nlevdecomp_full)) ; this%primp_vr_prev        (:,:) = nan
    allocate(this%totpftp_beg          (begc:endc))                   ; this%totpftp_beg          (:)   = nan
    allocate(this%solutionp_beg        (begc:endc))                   ; this%solutionp_beg        (:)   = nan
    allocate(this%labilep_beg          (begc:endc))                   ; this%labilep_beg          (:)   = nan
    allocate(this%secondp_beg          (begc:endc))                   ; this%secondp_beg          (:)   = nan
    allocate(this%totlitp_beg          (begc:endc))                   ; this%totlitp_beg          (:)   = nan
    allocate(this%cwdp_beg             (begc:endc))                   ; this%cwdp_beg             (:)   = nan
    allocate(this%totsomp_beg          (begc:endc))                   ; this%totsomp_beg          (:)   = nan
    allocate(this%totlitp_end          (begc:endc))                   ; this%totlitp_end          (:)   = nan
    allocate(this%totpftp_end          (begc:endc))                   ; this%totpftp_end          (:)   = nan
    allocate(this%labilep_end          (begc:endc))                   ; this%labilep_end          (:)   = nan
    allocate(this%secondp_end          (begc:endc))                   ; this%secondp_end          (:)   = nan
    allocate(this%solutionp_end        (begc:endc))                   ; this%solutionp_end        (:)   = nan
    allocate(this%cwdp_end             (begc:endc))                   ; this%cwdp_end             (:)   = nan
    allocate(this%totsomp_end          (begc:endc))                   ; this%totsomp_end          (:)   = nan
    allocate(this%cropseedp_deficit    (begc:endc))                   ; this%cropseedp_deficit    (:)   = nan

  end subroutine col_ps_init
    
  
  !------------------------------------------------------------------------
  subroutine col_ps_clean(this)
    !
    ! !ARGUMENTS:
    class(column_phosphorus_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine col_ps_clean

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column energy flux data structure
  !------------------------------------------------------------------------
  subroutine col_ef_init(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_energy_flux) :: this
    integer, intent(in) :: begc,endc
    ! !LOCAL VARIABLES:
    integer  :: l,c
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of col_ef
    !-----------------------------------------------------------------------
    allocate(this%eflx_h2osfc_to_snow  (begc:endc))              ; this%eflx_h2osfc_to_snow  (:)   = nan
    allocate(this%eflx_snomelt         (begc:endc))              ; this%eflx_snomelt         (:)   = nan
    allocate(this%eflx_snomelt_r       (begc:endc))              ; this%eflx_snomelt_r       (:)   = nan
    allocate(this%eflx_snomelt_u       (begc:endc))              ; this%eflx_snomelt_u       (:)   = nan
    allocate(this%eflx_bot             (begc:endc))              ; this%eflx_bot             (:)   = nan
    allocate(this%eflx_fgr12           (begc:endc))              ; this%eflx_fgr12           (:)   = nan
    allocate(this%eflx_fgr             (begc:endc, 1:nlevgrnd))  ; this%eflx_fgr             (:,:) = nan
    allocate(this%eflx_building_heat   (begc:endc))              ; this%eflx_building_heat   (:)   = nan
    allocate(this%eflx_urban_ac        (begc:endc))              ; this%eflx_urban_ac        (:)   = nan
    allocate(this%eflx_urban_heat      (begc:endc))              ; this%eflx_urban_heat      (:)   = nan
    allocate(this%eflx_hs_h2osfc       (begc:endc))              ; this%eflx_hs_h2osfc       (:)   = nan
    allocate(this%eflx_hs_top_snow     (begc:endc))              ; this%eflx_hs_top_snow     (:)   = nan
    allocate(this%eflx_hs_soil         (begc:endc))              ; this%eflx_hs_soil         (:)   = nan
    allocate(this%eflx_sabg_lyr        (begc:endc, -nlevsno+1:1)); this%eflx_sabg_lyr        (:,:) = nan
    allocate(this%eflx_dhsdT           (begc:endc))              ; this%eflx_dhsdT           (:)   = nan
    allocate(this%htvp                 (begc:endc))              ; this%htvp                 (:)   = nan
    allocate(this%xmf                  (begc:endc))              ; this%xmf                  (:)   = nan
    allocate(this%xmf_h2osfc           (begc:endc))              ; this%xmf_h2osfc           (:)   = nan
    allocate(this%imelt                (begc:endc,-nlevsno+1:nlevgrnd))  ; this%imelt        (:,:) = huge(1)
    allocate(this%eflx_soil_grnd       (begc:endc))              ; this%eflx_soil_grnd       (:)   = nan
    allocate(this%eflx_rnet_soil       (begc:endc))              ; this%eflx_rnet_soil       (:)   = nan
    allocate(this%eflx_fgr0_soil       (begc:endc))              ; this%eflx_fgr0_soil       (:)   = nan
    allocate(this%eflx_fgr0_snow       (begc:endc))              ; this%eflx_fgr0_snow       (:)   = nan
    allocate(this%eflx_fgr0_h2osfc     (begc:endc))              ; this%eflx_fgr0_h2osfc     (:)   = nan
    allocate(this%errsoi               (begc:endc))              ; this%errsoi               (:)   = nan
    allocate(this%errseb               (begc:endc))              ; this%errseb               (:)   = nan
    allocate(this%errsol               (begc:endc))              ; this%errsol               (:)   = nan
    allocate(this%errlon               (begc:endc))              ; this%errlon               (:)   = nan

  end subroutine col_ef_init

  !------------------------------------------------------------------------
  subroutine col_ef_clean(this)
    !
    ! !ARGUMENTS:
    class(column_energy_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine col_ef_clean
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column water flux data structure
  !------------------------------------------------------------------------
  subroutine col_wf_init(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_water_flux) :: this
    integer, intent(in) :: begc,endc
    ! !LOCAL VARIABLES:
    integer  :: l,c
    integer  :: ncells
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of col_wf
    !-----------------------------------------------------------------------
    allocate(this%qflx_prec_grnd         (begc:endc))             ; this%qflx_prec_grnd       (:)   = nan
    allocate(this%qflx_rain_grnd         (begc:endc))             ; this%qflx_rain_grnd       (:)   = nan
    allocate(this%qflx_snow_grnd         (begc:endc))             ; this%qflx_snow_grnd       (:)   = nan
    allocate(this%qflx_sub_snow          (begc:endc))             ; this%qflx_sub_snow        (:)   = nan
    allocate(this%qflx_sub_snow_vol      (begc:endc))             ; this%qflx_sub_snow_vol    (:)   = nan
    allocate(this%qflx_evap_soi          (begc:endc))             ; this%qflx_evap_soi        (:)   = nan
    allocate(this%qflx_evap_veg          (begc:endc))             ; this%qflx_evap_veg        (:)   = nan
    allocate(this%qflx_evap_can          (begc:endc))             ; this%qflx_evap_can        (:)   = nan
    allocate(this%qflx_evap_tot          (begc:endc))             ; this%qflx_evap_tot        (:)   = nan
    allocate(this%qflx_evap_grnd         (begc:endc))             ; this%qflx_evap_grnd       (:)   = nan
    allocate(this%qflx_snwcp_liq         (begc:endc))             ; this%qflx_snwcp_liq       (:)   = nan
    allocate(this%qflx_snwcp_ice         (begc:endc))             ; this%qflx_snwcp_ice       (:)   = nan
    allocate(this%qflx_tran_veg          (begc:endc))             ; this%qflx_tran_veg        (:)   = nan
    allocate(this%qflx_dew_snow          (begc:endc))             ; this%qflx_dew_snow        (:)   = nan
    allocate(this%qflx_dew_grnd          (begc:endc))             ; this%qflx_dew_grnd        (:)   = nan
    allocate(this%qflx_prec_intr         (begc:endc))             ; this%qflx_prec_intr       (:)   = nan
    allocate(this%qflx_ev_snow           (begc:endc))             ; this%qflx_ev_snow         (:)   = nan
    allocate(this%qflx_ev_soil           (begc:endc))             ; this%qflx_ev_soil         (:)   = nan
    allocate(this%qflx_ev_h2osfc         (begc:endc))             ; this%qflx_ev_h2osfc       (:)   = nan
    allocate(this%qflx_gross_evap_soil   (begc:endc))             ; this%qflx_gross_evap_soil (:)   = nan
    allocate(this%qflx_gross_infl_soil   (begc:endc))             ; this%qflx_gross_infl_soil (:)   = nan
    allocate(this%qflx_adv               (begc:endc,0:nlevgrnd))  ; this%qflx_adv             (:,:) = nan
    allocate(this%qflx_rootsoi           (begc:endc,1:nlevgrnd))  ; this%qflx_rootsoi         (:,:) = nan
    allocate(this%dwb                    (begc:endc))             ; this%dwb                  (:)   = nan
    allocate(this%qflx_infl              (begc:endc))             ; this%qflx_infl            (:)   = nan
    allocate(this%qflx_surf              (begc:endc))             ; this%qflx_surf            (:)   = nan
    allocate(this%qflx_drain             (begc:endc))             ; this%qflx_drain           (:)   = nan
    allocate(this%qflx_totdrain          (begc:endc))             ; this%qflx_totdrain        (:)   = nan
    allocate(this%qflx_top_soil          (begc:endc))             ; this%qflx_top_soil        (:)   = nan
    allocate(this%qflx_h2osfc_to_ice     (begc:endc))             ; this%qflx_h2osfc_to_ice   (:)   = nan
    allocate(this%qflx_h2osfc_surf       (begc:endc))             ; this%qflx_h2osfc_surf     (:)   = nan
    allocate(this%qflx_snow_h2osfc       (begc:endc))             ; this%qflx_snow_h2osfc     (:)   = nan
    allocate(this%qflx_drain_perched     (begc:endc))             ; this%qflx_drain_perched   (:)   = nan
    allocate(this%qflx_deficit           (begc:endc))             ; this%qflx_deficit         (:)   = nan
    allocate(this%qflx_floodc            (begc:endc))             ; this%qflx_floodc          (:)   = nan
    allocate(this%qflx_sl_top_soil       (begc:endc))             ; this%qflx_sl_top_soil     (:)   = nan
    allocate(this%qflx_snomelt           (begc:endc))             ; this%qflx_snomelt         (:)   = nan
    allocate(this%qflx_snow_melt         (begc:endc))             ; this%qflx_snow_melt       (:)   = nan
    allocate(this%qflx_qrgwl             (begc:endc))             ; this%qflx_qrgwl           (:)   = nan
    allocate(this%qflx_runoff            (begc:endc))             ; this%qflx_runoff          (:)   = nan
    allocate(this%qflx_runoff_r          (begc:endc))             ; this%qflx_runoff_r        (:)   = nan
    allocate(this%qflx_runoff_u          (begc:endc))             ; this%qflx_runoff_u        (:)   = nan
    allocate(this%qflx_rsub_sat          (begc:endc))             ; this%qflx_rsub_sat        (:)   = nan
    allocate(this%qflx_snofrz_lyr        (begc:endc,-nlevsno+1:0)); this%qflx_snofrz_lyr      (:,:) = nan
    allocate(this%qflx_snofrz            (begc:endc))             ; this%qflx_snofrz          (:)   = nan
    allocate(this%qflx_glcice            (begc:endc))             ; this%qflx_glcice          (:)   = nan
    allocate(this%qflx_glcice_frz        (begc:endc))             ; this%qflx_glcice_frz      (:)   = nan
    allocate(this%qflx_glcice_melt       (begc:endc))             ; this%qflx_glcice_melt     (:)   = nan
    allocate(this%qflx_drain_vr          (begc:endc,1:nlevgrnd))  ; this%qflx_drain_vr        (:,:) = nan
    allocate(this%qflx_h2osfc2topsoi     (begc:endc))             ; this%qflx_h2osfc2topsoi   (:)   = nan
    allocate(this%qflx_snow2topsoi       (begc:endc))             ; this%qflx_snow2topsoi     (:)   = nan
    allocate(this%qflx_lateral           (begc:endc))             ; this%qflx_lateral         (:)   = 0._r8
    allocate(this%snow_sources           (begc:endc))             ; this%snow_sources         (:)   = nan
    allocate(this%snow_sinks             (begc:endc))             ; this%snow_sinks           (:)   = nan
    allocate(this%qflx_irrig             (begc:endc))             ; this%qflx_irrig           (:)   = nan
    
    !VSFM variables
    ncells = endc - begc + 1
    allocate(this%mflx_infl_1d           (ncells))                ; this%mflx_infl_1d                    (:)   = nan
    allocate(this%mflx_dew_1d            (ncells))                ; this%mflx_dew_1d                     (:)   = nan
    allocate(this%mflx_snowlyr_1d        (ncells))                ; this%mflx_snowlyr_1d                 (:)   = nan
    allocate(this%mflx_sub_snow_1d       (ncells))                ; this%mflx_sub_snow_1d                (:)   = nan
    allocate(this%mflx_neg_snow_1d       (ncells))                ; this%mflx_neg_snow_1d                (:)   = nan
    ncells = (endc - begc + 1)*nlevgrnd
    allocate(this%mflx_et_1d             (ncells))                ; this%mflx_et_1d                      (:)   = nan
    allocate(this%mflx_drain_1d          (ncells))                ; this%mflx_drain_1d                   (:)   = nan
    allocate(this%mflx_drain_perched_1d  (ncells))                ; this%mflx_drain_perched_1d           (:)   = nan

    allocate(this%mflx_snowlyr           (begc:endc))             ; this%mflx_snowlyr                    (:)   = 0._r8
    allocate(this%mflx_infl              (begc:endc))             ; this%mflx_infl                       (:)   = nan
    allocate(this%mflx_dew               (begc:endc))             ; this%mflx_dew                        (:)   = nan
    allocate(this%mflx_snowlyr_disp      (begc:endc))             ; this%mflx_snowlyr_disp               (:)   = nan
    allocate(this%mflx_sub_snow          (begc:endc))             ; this%mflx_sub_snow                   (:)   = nan
    allocate(this%mflx_et                (begc:endc,1:nlevgrnd))  ; this%mflx_et                         (:,:) = nan
    allocate(this%mflx_drain             (begc:endc,1:nlevgrnd))  ; this%mflx_drain                      (:,:) = nan
    allocate(this%mflx_recharge          (begc:endc))             ; this%mflx_recharge                   (:)   = nan

  end subroutine col_wf_init

  !------------------------------------------------------------------------
  subroutine col_wf_clean(this)
    !
    ! !ARGUMENTS:
    class(column_water_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine col_wf_clean
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column carbon flux data structure
  !------------------------------------------------------------------------
  subroutine col_cf_init(this, begc, endc, carbon_type)
    !
    ! !ARGUMENTS:
    class(column_carbon_flux) :: this
    integer, intent(in) :: begc,endc
    character(len=3) , intent(in) :: carbon_type ! one of ['c12', c13','c14']
    !
    ! !LOCAL VARIABLES:
    integer           :: c,j                        ! indices
    integer           :: k,l,ii,jj 
    character(8)      :: vr_suffix
    character(10)     :: active
    character(24)     :: fieldname
    character(100)    :: longname
    real(r8), pointer :: data1dptr(:)   ! temp. pointer for slicing larger arrays
    real(r8), pointer :: data2dptr(:,:) ! temp. pointer for slicing larger arrays
    character(len=3)  :: ctag
    integer           :: fc                                        ! filter index
    integer           :: num_special_col                           ! number of good values in special_col filter
    integer           :: special_col(endc-begc+1)    ! special landunit filter - columns
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of col_cf
    !-----------------------------------------------------------------------
    allocate(this%decomp_cpools_sourcesink          (begc:endc,1:nlevdecomp_full,1:ndecomp_pools))               ; this%decomp_cpools_sourcesink   (:,:,:) = nan
    allocate(this%decomp_cascade_hr_vr              (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)) ; this%decomp_cascade_hr_vr       (:,:,:) = spval
    allocate(this%decomp_cascade_ctransfer_vr       (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)) ; this%decomp_cascade_ctransfer_vr(:,:,:) = nan
    allocate(this%decomp_k                          (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)) ; this%decomp_k                   (:,:,:) = spval
    allocate(this%decomp_cpools_transport_tendency  (begc:endc,1:nlevdecomp_full,1:ndecomp_pools)) ; this%decomp_cpools_transport_tendency(:,:,:) = nan
    allocate(this%decomp_cascade_hr                 (begc:endc,1:ndecomp_cascade_transitions))     ; this%decomp_cascade_hr               (:,:)   = nan  
    allocate(this%decomp_cascade_ctransfer          (begc:endc,1:ndecomp_cascade_transitions))     ; this%decomp_cascade_ctransfer        (:,:)   = nan  
    allocate(this%o_scalar                          (begc:endc,1:nlevdecomp_full)); this%o_scalar                     (:,:) = spval  
    allocate(this%w_scalar                          (begc:endc,1:nlevdecomp_full)); this%w_scalar                     (:,:) = spval  
    allocate(this%t_scalar                          (begc:endc,1:nlevdecomp_full)); this%t_scalar                     (:,:) = spval  
    allocate(this%decomp_cpools_leached             (begc:endc,1:ndecomp_pools))  ; this%decomp_cpools_leached        (:,:) = nan  
    allocate(this%phr_vr                            (begc:endc,1:nlevdecomp_full)); this%phr_vr                       (:,:) = nan  
    allocate(this%fphr                              (begc:endc,1:nlevgrnd))       ; this%fphr                         (:,:) = nan  
    allocate(this%som_c_leached                     (begc:endc))                  ; this%som_c_leached                (:)   = nan    
    allocate(this%phenology_c_to_litr_met_c         (begc:endc,1:nlevdecomp_full)); this%phenology_c_to_litr_met_c    (:,:) = nan  
    allocate(this%phenology_c_to_litr_cel_c         (begc:endc,1:nlevdecomp_full)); this%phenology_c_to_litr_cel_c    (:,:) = nan  
    allocate(this%phenology_c_to_litr_lig_c         (begc:endc,1:nlevdecomp_full)); this%phenology_c_to_litr_lig_c    (:,:) = nan  
    allocate(this%gap_mortality_c_to_litr_met_c     (begc:endc,1:nlevdecomp_full)); this%gap_mortality_c_to_litr_met_c(:,:) = nan  
    allocate(this%gap_mortality_c_to_litr_cel_c     (begc:endc,1:nlevdecomp_full)); this%gap_mortality_c_to_litr_cel_c(:,:) = nan  
    allocate(this%gap_mortality_c_to_litr_lig_c     (begc:endc,1:nlevdecomp_full)); this%gap_mortality_c_to_litr_lig_c(:,:) = nan  
    allocate(this%gap_mortality_c_to_cwdc           (begc:endc,1:nlevdecomp_full)); this%gap_mortality_c_to_cwdc      (:,:) = nan  
    allocate(this%m_decomp_cpools_to_fire_vr        (begc:endc,1:nlevdecomp_full,1:ndecomp_pools)) ; this%m_decomp_cpools_to_fire_vr(:,:,:)= nan
    allocate(this%m_decomp_cpools_to_fire           (begc:endc,1:ndecomp_pools))  ; this%m_decomp_cpools_to_fire      (:,:) = nan  
    allocate(this%m_c_to_litr_met_fire              (begc:endc,1:nlevdecomp_full)); this%m_c_to_litr_met_fire         (:,:) = nan  
    allocate(this%m_c_to_litr_cel_fire              (begc:endc,1:nlevdecomp_full)); this%m_c_to_litr_cel_fire         (:,:) = nan  
    allocate(this%m_c_to_litr_lig_fire              (begc:endc,1:nlevdecomp_full)); this%m_c_to_litr_lig_fire         (:,:) = nan  
    allocate(this%fire_mortality_c_to_cwdc          (begc:endc,1:nlevdecomp_full)); this%fire_mortality_c_to_cwdc     (:,:) = nan  
    allocate(this%somc_fire                         (begc:endc))                  ; this%somc_fire                    (:)   = nan    
    allocate(this%harvest_c_to_litr_met_c           (begc:endc,1:nlevdecomp_full)); this%harvest_c_to_litr_met_c      (:,:) = nan  
    allocate(this%harvest_c_to_litr_cel_c           (begc:endc,1:nlevdecomp_full)); this%harvest_c_to_litr_cel_c      (:,:) = nan  
    allocate(this%harvest_c_to_litr_lig_c           (begc:endc,1:nlevdecomp_full)); this%harvest_c_to_litr_lig_c      (:,:) = nan  
    allocate(this%harvest_c_to_cwdc                 (begc:endc,1:nlevdecomp_full)); this%harvest_c_to_cwdc            (:,:) = nan  
    allocate(this%hrv_deadstemc_to_prod10c          (begc:endc))                  ; this%hrv_deadstemc_to_prod10c     (:)   = nan    
    allocate(this%hrv_deadstemc_to_prod100c         (begc:endc))                  ; this%hrv_deadstemc_to_prod100c    (:)   = nan   
    allocate(this%hrv_cropc_to_prod1c               (begc:endc))                  ; this%hrv_cropc_to_prod1c          (:)   = nan   
    allocate(this%dwt_frootc_to_litr_met_c          (begc:endc,1:nlevdecomp_full)); this%dwt_frootc_to_litr_met_c     (:,:) = nan  
    allocate(this%dwt_frootc_to_litr_cel_c          (begc:endc,1:nlevdecomp_full)); this%dwt_frootc_to_litr_cel_c     (:,:) = nan  
    allocate(this%dwt_frootc_to_litr_lig_c          (begc:endc,1:nlevdecomp_full)); this%dwt_frootc_to_litr_lig_c     (:,:) = nan  
    allocate(this%dwt_livecrootc_to_cwdc            (begc:endc,1:nlevdecomp_full)); this%dwt_livecrootc_to_cwdc       (:,:) = nan  
    allocate(this%dwt_deadcrootc_to_cwdc            (begc:endc,1:nlevdecomp_full)); this%dwt_deadcrootc_to_cwdc       (:,:) = nan  
    allocate(this%dwt_slash_cflux                   (begc:endc))                  ; this%dwt_slash_cflux              (:)   = nan    
    allocate(this%dwt_conv_cflux                    (begc:endc))                  ; this%dwt_conv_cflux               (:)   = nan    
    allocate(this%dwt_prod10c_gain                  (begc:endc))                  ; this%dwt_prod10c_gain             (:)   = nan    
    allocate(this%dwt_prod100c_gain                 (begc:endc))                  ; this%dwt_prod100c_gain            (:)   = nan    
    allocate(this%dwt_closs                         (begc:endc))                  ; this%dwt_closs                    (:)   = nan    
    allocate(this%landuseflux                       (begc:endc))                  ; this%landuseflux                  (:)   = nan    
    allocate(this%landuptake                        (begc:endc))                  ; this%landuptake                   (:)   = nan    
    allocate(this%prod1c_loss                       (begc:endc))                  ; this%prod1c_loss                  (:)   = nan    
    allocate(this%prod10c_loss                      (begc:endc))                  ; this%prod10c_loss                 (:)   = nan    
    allocate(this%prod100c_loss                     (begc:endc))                  ; this%prod100c_loss                (:)   = nan    
    allocate(this%product_closs                     (begc:endc))                  ; this%product_closs                (:)   = nan    
    allocate(this%hr_vr                             (begc:endc,1:nlevdecomp_full)); this%hr_vr                        (:,:) = nan  
    allocate(this%lithr                             (begc:endc))                  ; this%lithr                        (:)   = nan    
    allocate(this%somhr                             (begc:endc))                  ; this%somhr                        (:)   = nan    
    allocate(this%hr                                (begc:endc))                  ; this%hr                           (:)   = nan    
    allocate(this%sr                                (begc:endc))                  ; this%sr                           (:)   = nan    
    allocate(this%er                                (begc:endc))                  ; this%er                           (:)   = nan    
    allocate(this%litfire                           (begc:endc))                  ; this%litfire                      (:)   = nan    
    allocate(this%somfire                           (begc:endc))                  ; this%somfire                      (:)   = nan    
    allocate(this%totfire                           (begc:endc))                  ; this%totfire                      (:)   = nan    
    allocate(this%nep                               (begc:endc))                  ; this%nep                          (:)   = nan    
    allocate(this%nbp                               (begc:endc))                  ; this%nbp                          (:)   = nan    
    allocate(this%nee                               (begc:endc))                  ; this%nee                          (:)   = nan    
    allocate(this%bgc_cpool_ext_inputs_vr           (begc:endc, 1:nlevdecomp_full,ndecomp_pools)) ; this%bgc_cpool_ext_inputs_vr(:,:,:) = nan
    allocate(this%bgc_cpool_ext_loss_vr             (begc:endc, 1:nlevdecomp_full,ndecomp_pools)) ; this%bgc_cpool_ext_loss_vr  (:,:,:) = nan
    allocate(this%cwdc_hr                           (begc:endc))                  ; this%cwdc_hr                      (:)   = nan    
    allocate(this%cwdc_loss                         (begc:endc))                  ; this%cwdc_loss                    (:)   = nan    
    allocate(this%litterc_loss                      (begc:endc))                  ; this%litterc_loss                 (:)   = nan    
    allocate(this%rr                                (begc:endc))                  ; this%rr                           (:)   = nan    
    allocate(this%ar                                (begc:endc))                  ; this%ar                           (:)   = nan    
    allocate(this%gpp                               (begc:endc))                  ; this%gpp                          (:)   = nan    
    allocate(this%npp                               (begc:endc))                  ; this%npp                          (:)   = nan    
    allocate(this%fire_closs_p2c                    (begc:endc))                  ; this%fire_closs_p2c               (:)   = nan    
    allocate(this%fire_closs                        (begc:endc))                  ; this%fire_closs                   (:)   = nan    
    allocate(this%fire_decomp_closs                 (begc:endc))                  ; this%fire_decomp_closs            (:)   = nan    
    allocate(this%litfall                           (begc:endc))                  ; this%litfall                      (:)   = nan    
    allocate(this%vegfire                           (begc:endc))                  ; this%vegfire                      (:)   = nan    
    allocate(this%wood_harvestc                     (begc:endc))                  ; this%wood_harvestc                (:)   = nan    
    allocate(this%hrv_xsmrpool_to_atm               (begc:endc))                  ; this%hrv_xsmrpool_to_atm          (:)   = nan    
    allocate(this%plant_to_litter_cflux             (begc:endc))                  ; this%plant_to_litter_cflux        (:)   = nan              
    allocate(this%plant_to_cwd_cflux	             (begc:endc))                  ; this%plant_to_cwd_cflux		       (:)   = nan 
    allocate(this%annsum_npp                        (begc:endc))                  ; this%annsum_npp                   (:)   = nan 
    allocate(this%lag_npp                           (begc:endc))                  ; this%lag_npp                      (:)   = spval 
    allocate(this%externalc_to_decomp_cpools        (begc:endc,1:nlevdecomp_full,1:ndecomp_pools)) ; this%externalc_to_decomp_cpools(:,:,:) = spval
    allocate(this%externalc_to_decomp_delta         (begc:endc))                  ; this%externalc_to_decomp_delta    (:)   = spval    
    allocate(this%f_co2_soil_vr                     (begc:endc,1:nlevdecomp_full)); this%f_co2_soil_vr                (:,:) = nan  
    allocate(this%f_co2_soil                        (begc:endc))                  ; this%f_co2_soil                   (:)   = nan    

end subroutine col_cf_init

  !------------------------------------------------------------------------
  subroutine col_cf_clean(this)
    !
    ! !ARGUMENTS:
    class(column_carbon_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine col_cf_clean
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column nitrogen flux data structure
  !------------------------------------------------------------------------
  subroutine col_nf_init(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_nitrogen_flux) :: this
    integer, intent(in) :: begc,endc
    !
    ! !LOCAL VARIABLES:
    integer        :: k,l
    character(10)  :: active
    character(24)  :: fieldname
    character(100) :: longname
    character(8)   :: vr_suffix
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    integer :: c
    integer :: fc                          ! filter indices
    integer :: num_special_col             ! number of good values in special_col filter
    integer :: special_col(endc-begc+1)    ! special landunit filter - columns
    !------------------------------------------------------------------------
    
    !-----------------------------------------------------------------------
    ! allocate for each member of col_nf
    !-----------------------------------------------------------------------
    allocate(this%ndep_to_sminn                   (begc:endc))                   ; this%ndep_to_sminn	                 (:)   = nan
    allocate(this%nfix_to_sminn                   (begc:endc))                   ; this%nfix_to_sminn	                 (:)   = nan
    allocate(this%nfix_to_ecosysn                 (begc:endc))                   ; this%nfix_to_ecosysn                (:)   = nan
    allocate(this%fert_to_sminn                   (begc:endc))                   ; this%fert_to_sminn	                 (:)   = nan
    allocate(this%soyfixn_to_sminn                (begc:endc))                   ; this%soyfixn_to_sminn               (:)   = nan
    allocate(this%hrv_deadstemn_to_prod10n        (begc:endc))                   ; this%hrv_deadstemn_to_prod10n       (:)   = nan
    allocate(this%hrv_deadstemn_to_prod100n       (begc:endc))                   ; this%hrv_deadstemn_to_prod100n      (:)   = nan
    allocate(this%hrv_cropn_to_prod1n             (begc:endc))                   ; this%hrv_cropn_to_prod1n            (:)   = nan
    allocate(this%sminn_to_plant                  (begc:endc))                   ; this%sminn_to_plant	              (:)   = nan
    allocate(this%potential_immob                 (begc:endc))                   ; this%potential_immob                (:)   = nan
    allocate(this%actual_immob                    (begc:endc))                   ; this%actual_immob                   (:)   = nan
    allocate(this%gross_nmin                      (begc:endc))                   ; this%gross_nmin                     (:)   = nan
    allocate(this%net_nmin                        (begc:endc))                   ; this%net_nmin                       (:)   = nan
    allocate(this%denit                           (begc:endc))                   ; this%denit		                    (:)   = nan
    allocate(this%supplement_to_sminn             (begc:endc))                   ; this%supplement_to_sminn            (:)   = nan
    allocate(this%prod1n_loss                     (begc:endc))                   ; this%prod1n_loss                    (:)   = nan
    allocate(this%prod10n_loss                    (begc:endc))                   ; this%prod10n_loss                   (:)   = nan
    allocate(this%prod100n_loss                   (begc:endc))                   ; this%prod100n_loss	                 (:)   = nan
    allocate(this%product_nloss                   (begc:endc))                   ; this%product_nloss	                 (:)   = nan
    allocate(this%ninputs                         (begc:endc))                   ; this%ninputs                        (:)   = nan
    allocate(this%noutputs                        (begc:endc))                   ; this%noutputs                       (:)   = nan
    allocate(this%fire_nloss                      (begc:endc))                   ; this%fire_nloss                     (:)   = nan
    allocate(this%fire_decomp_nloss               (begc:endc))                   ; this%fire_decomp_nloss              (:)   = nan
    allocate(this%fire_nloss_p2c                  (begc:endc))                   ; this%fire_nloss_p2c                 (:)   = nan
    allocate(this%som_n_leached                   (begc:endc))                   ; this%som_n_leached	                 (:)   = nan
    allocate(this%m_n_to_litr_met_fire            (begc:endc,1:nlevdecomp_full)) ; this%m_n_to_litr_met_fire           (:,:) = nan
    allocate(this%m_n_to_litr_cel_fire            (begc:endc,1:nlevdecomp_full)) ; this%m_n_to_litr_cel_fire           (:,:) = nan
    allocate(this%m_n_to_litr_lig_fire            (begc:endc,1:nlevdecomp_full)) ; this%m_n_to_litr_lig_fire           (:,:) = nan
    allocate(this%r_psi                           (begc:endc,1:nlevdecomp_full)) ; this%r_psi                          (:,:) = spval
    allocate(this%anaerobic_frac                  (begc:endc,1:nlevdecomp_full)) ; this%anaerobic_frac                 (:,:) = spval
    allocate(this%potential_immob_vr              (begc:endc,1:nlevdecomp_full)) ; this%potential_immob_vr             (:,:) = nan
    allocate(this%actual_immob_vr                 (begc:endc,1:nlevdecomp_full)) ; this%actual_immob_vr                (:,:) = nan
    allocate(this%sminn_to_plant_vr               (begc:endc,1:nlevdecomp_full)) ; this%sminn_to_plant_vr              (:,:) = nan
    allocate(this%supplement_to_sminn_vr          (begc:endc,1:nlevdecomp_full)) ; this%supplement_to_sminn_vr         (:,:) = nan
    allocate(this%gross_nmin_vr                   (begc:endc,1:nlevdecomp_full)) ; this%gross_nmin_vr                  (:,:) = nan
    allocate(this%net_nmin_vr                     (begc:endc,1:nlevdecomp_full)) ; this%net_nmin_vr                    (:,:) = nan
    allocate(this%dwt_slash_nflux                 (begc:endc))                   ; this%dwt_slash_nflux                (:)   = nan
    allocate(this%dwt_conv_nflux                  (begc:endc))                   ; this%dwt_conv_nflux                 (:)   = nan
    allocate(this%dwt_prod10n_gain                (begc:endc))                   ; this%dwt_prod10n_gain               (:)   = nan
    allocate(this%dwt_prod100n_gain               (begc:endc))                   ; this%dwt_prod100n_gain              (:)   = nan
    allocate(this%dwt_nloss                       (begc:endc))                   ; this%dwt_nloss                      (:)   = nan
    allocate(this%wood_harvestn                   (begc:endc))                   ; this%wood_harvestn                  (:)   = nan
    allocate(this%dwt_frootn_to_litr_met_n        (begc:endc,1:nlevdecomp_full)) ; this%dwt_frootn_to_litr_met_n       (:,:) = nan
    allocate(this%dwt_frootn_to_litr_cel_n        (begc:endc,1:nlevdecomp_full)) ; this%dwt_frootn_to_litr_cel_n       (:,:) = nan
    allocate(this%dwt_frootn_to_litr_lig_n        (begc:endc,1:nlevdecomp_full)) ; this%dwt_frootn_to_litr_lig_n       (:,:) = nan
    allocate(this%dwt_livecrootn_to_cwdn          (begc:endc,1:nlevdecomp_full)) ; this%dwt_livecrootn_to_cwdn         (:,:) = nan
    allocate(this%dwt_deadcrootn_to_cwdn          (begc:endc,1:nlevdecomp_full)) ; this%dwt_deadcrootn_to_cwdn         (:,:) = nan
    allocate(this%f_nit_vr                        (begc:endc,1:nlevdecomp_full)) ; this%f_nit_vr                       (:,:) = nan
    allocate(this%f_denit_vr                      (begc:endc,1:nlevdecomp_full)) ; this%f_denit_vr                     (:,:) = nan
    allocate(this%smin_no3_leached_vr             (begc:endc,1:nlevdecomp_full)) ; this%smin_no3_leached_vr            (:,:) = nan
    allocate(this%smin_no3_leached                (begc:endc))                   ; this%smin_no3_leached               (:)   = nan
    allocate(this%smin_no3_runoff_vr              (begc:endc,1:nlevdecomp_full)) ; this%smin_no3_runoff_vr             (:,:) = nan
    allocate(this%smin_no3_runoff                 (begc:endc))                   ; this%smin_no3_runoff                (:)   = nan
    allocate(this%pot_f_nit_vr                    (begc:endc,1:nlevdecomp_full)) ; this%pot_f_nit_vr                   (:,:) = nan
    allocate(this%pot_f_nit                       (begc:endc))                   ; this%pot_f_nit                      (:)   = nan
    allocate(this%pot_f_denit_vr                  (begc:endc,1:nlevdecomp_full)) ; this%pot_f_denit_vr                 (:,:) = nan
    allocate(this%pot_f_denit                     (begc:endc))                   ; this%pot_f_denit                    (:)   = nan
    allocate(this%actual_immob_no3_vr             (begc:endc,1:nlevdecomp_full)) ; this%actual_immob_no3_vr            (:,:) = nan
    allocate(this%actual_immob_nh4_vr             (begc:endc,1:nlevdecomp_full)) ; this%actual_immob_nh4_vr            (:,:) = nan
    allocate(this%smin_no3_to_plant_vr            (begc:endc,1:nlevdecomp_full)) ; this%smin_no3_to_plant_vr           (:,:) = nan
    allocate(this%smin_nh4_to_plant_vr            (begc:endc,1:nlevdecomp_full)) ; this%smin_nh4_to_plant_vr           (:,:) = nan
    allocate(this%smin_no3_to_plant               (begc:endc))                   ; this%smin_no3_to_plant              (:)   = nan
    allocate(this%smin_nh4_to_plant               (begc:endc))                   ; this%smin_nh4_to_plant              (:)   = nan
    allocate(this%f_nit                           (begc:endc))                   ; this%f_nit                          (:)   = nan
    allocate(this%f_denit                         (begc:endc))                   ; this%f_denit                        (:)   = nan
    allocate(this%n2_n2o_ratio_denit_vr           (begc:endc,1:nlevdecomp_full)) ; this%n2_n2o_ratio_denit_vr          (:,:) = nan
    allocate(this%f_n2o_denit                     (begc:endc))                   ; this%f_n2o_denit                    (:)   = nan
    allocate(this%f_n2o_denit_vr                  (begc:endc,1:nlevdecomp_full)) ; this%f_n2o_denit_vr                 (:,:) = nan
    allocate(this%f_n2o_nit                       (begc:endc))                   ; this%f_n2o_nit                      (:)   = nan
    allocate(this%f_n2o_nit_vr                    (begc:endc,1:nlevdecomp_full)) ; this%f_n2o_nit_vr                   (:,:) = nan
    allocate(this%sminn_no3_input_vr              (begc:endc,1:nlevdecomp_full)) ; this%sminn_no3_input_vr             (:,:) = nan
    allocate(this%sminn_nh4_input_vr              (begc:endc,1:nlevdecomp_full)) ; this%sminn_nh4_input_vr             (:,:) = nan
    allocate(this%sminn_nh4_input                 (begc:endc))                   ; this%sminn_nh4_input                (:)   = nan
    allocate(this%sminn_no3_input                 (begc:endc))                   ; this%sminn_no3_input                (:)   = nan
    allocate(this%sminn_input                     (begc:endc))                   ; this%sminn_input                    (:)   = nan
    allocate(this%bgc_npool_inputs                (begc:endc,ndecomp_pools))     ; this%bgc_npool_inputs               (:,:) = nan
    allocate(this%smin_no3_massdens_vr            (begc:endc,1:nlevdecomp_full)) ; this%smin_no3_massdens_vr           (:,:) = nan
    allocate(this%soil_bulkdensity                (begc:endc,1:nlevdecomp_full)) ; this%soil_bulkdensity               (:,:) = nan
    allocate(this%k_nitr_t_vr                     (begc:endc,1:nlevdecomp_full)) ; this%k_nitr_t_vr                    (:,:) = nan
    allocate(this%k_nitr_ph_vr                    (begc:endc,1:nlevdecomp_full)) ; this%k_nitr_ph_vr                   (:,:) = nan
    allocate(this%k_nitr_h2o_vr                   (begc:endc,1:nlevdecomp_full)) ; this%k_nitr_h2o_vr                  (:,:) = nan
    allocate(this%k_nitr_vr                       (begc:endc,1:nlevdecomp_full)) ; this%k_nitr_vr                      (:,:) = nan
    allocate(this%wfps_vr                         (begc:endc,1:nlevdecomp_full)) ; this%wfps_vr                        (:,:) = nan
    allocate(this%f_denit_base_vr                 (begc:endc,1:nlevdecomp_full)) ; this%f_denit_base_vr                (:,:) = nan
    allocate(this%diffus                          (begc:endc,1:nlevdecomp_full)) ; this%diffus                         (:,:) = spval
    allocate(this%ratio_k1                        (begc:endc,1:nlevdecomp_full)) ; this%ratio_k1                       (:,:) = nan
    allocate(this%ratio_no3_co2                   (begc:endc,1:nlevdecomp_full)) ; this%ratio_no3_co2                  (:,:) = spval
    allocate(this%soil_co2_prod                   (begc:endc,1:nlevdecomp_full)) ; this%soil_co2_prod                  (:,:) = nan
    allocate(this%fr_WFPS                         (begc:endc,1:nlevdecomp_full)) ; this%fr_WFPS                        (:,:) = spval
    allocate(this%fmax_denit_carbonsubstrate_vr   (begc:endc,1:nlevdecomp_full)) ; this%fmax_denit_carbonsubstrate_vr  (:,:) = nan
    allocate(this%fmax_denit_nitrate_vr           (begc:endc,1:nlevdecomp_full)) ; this%fmax_denit_nitrate_vr          (:,:) = nan
    allocate(this%phenology_n_to_litr_met_n       (begc:endc, 1:nlevdecomp_full)) ; this%phenology_n_to_litr_met_n     (:,:) = nan
    allocate(this%phenology_n_to_litr_cel_n       (begc:endc, 1:nlevdecomp_full)) ; this%phenology_n_to_litr_cel_n     (:,:) = nan
    allocate(this%phenology_n_to_litr_lig_n       (begc:endc, 1:nlevdecomp_full)) ; this%phenology_n_to_litr_lig_n     (:,:) = nan
    allocate(this%gap_mortality_n_to_litr_met_n   (begc:endc, 1:nlevdecomp_full)) ; this%gap_mortality_n_to_litr_met_n (:,:) = nan
    allocate(this%gap_mortality_n_to_litr_cel_n   (begc:endc, 1:nlevdecomp_full)) ; this%gap_mortality_n_to_litr_cel_n (:,:) = nan
    allocate(this%gap_mortality_n_to_litr_lig_n   (begc:endc, 1:nlevdecomp_full)) ; this%gap_mortality_n_to_litr_lig_n (:,:) = nan
    allocate(this%gap_mortality_n_to_cwdn         (begc:endc, 1:nlevdecomp_full)) ; this%gap_mortality_n_to_cwdn       (:,:) = nan
    allocate(this%fire_mortality_n_to_cwdn        (begc:endc, 1:nlevdecomp_full)) ; this%fire_mortality_n_to_cwdn      (:,:) = nan
    allocate(this%harvest_n_to_litr_met_n         (begc:endc, 1:nlevdecomp_full)) ; this%harvest_n_to_litr_met_n       (:,:) = nan
    allocate(this%harvest_n_to_litr_cel_n         (begc:endc, 1:nlevdecomp_full)) ; this%harvest_n_to_litr_cel_n       (:,:) = nan
    allocate(this%harvest_n_to_litr_lig_n         (begc:endc, 1:nlevdecomp_full)) ; this%harvest_n_to_litr_lig_n       (:,:) = nan
    allocate(this%harvest_n_to_cwdn               (begc:endc, 1:nlevdecomp_full)) ; this%harvest_n_to_cwdn             (:,:) = nan
    allocate(this%plant_ndemand                   (begc:endc))                    ; this%plant_ndemand                 (:)   = nan
    allocate(this%plant_ndemand_vr                (begc:endc,1:nlevdecomp_full))  ; this%plant_ndemand_vr              (:,:) = nan
    allocate(this%f_ngas_decomp_vr                (begc:endc,1:nlevdecomp_full))  ; this%f_ngas_decomp_vr              (:,:) = nan
    allocate(this%f_ngas_nitri_vr                 (begc:endc,1:nlevdecomp_full))  ; this%f_ngas_nitri_vr               (:,:) = nan
    allocate(this%f_ngas_denit_vr                 (begc:endc,1:nlevdecomp_full))  ; this%f_ngas_denit_vr               (:,:) = nan
    allocate(this%f_n2o_soil_vr                   (begc:endc,1:nlevdecomp_full))  ; this%f_n2o_soil_vr                 (:,:) = nan
    allocate(this%f_n2_soil_vr                    (begc:endc,1:nlevdecomp_full))  ; this%f_n2_soil_vr                  (:,:) = nan
    allocate(this%f_ngas_decomp                   (begc:endc))                    ; this%f_ngas_decomp                 (:)   = nan
    allocate(this%f_ngas_nitri                    (begc:endc))                    ; this%f_ngas_nitri                  (:)   = nan
    allocate(this%f_ngas_denit                    (begc:endc))                    ; this%f_ngas_denit                  (:)   = nan
    allocate(this%f_n2o_soil                      (begc:endc))                    ; this%f_n2o_soil                    (:)   = nan
    allocate(this%f_n2_soil                       (begc:endc))                    ; this%f_n2_soil                     (:)   = nan
    allocate(this%externaln_to_decomp_delta       (begc:endc))                    ; this%externaln_to_decomp_delta     (:)   = spval
    allocate(this%no3_net_transport_vr            (begc:endc,1:nlevdecomp_full))  ; this%no3_net_transport_vr          (:,:) = spval
    allocate(this%nh4_net_transport_vr            (begc:endc,1:nlevdecomp_full))  ; this%nh4_net_transport_vr          (:,:) = spval
    allocate(this%col_plant_ndemand_vr            (begc:endc,1:nlevdecomp))       ; this%col_plant_ndemand_vr          (:,:) = nan
    allocate(this%col_plant_nh4demand_vr          (begc:endc,1:nlevdecomp))       ; this%col_plant_nh4demand_vr        (:,:) = nan
    allocate(this%col_plant_no3demand_vr          (begc:endc,1:nlevdecomp))       ; this%col_plant_no3demand_vr        (:,:) = nan
    allocate(this%col_plant_pdemand_vr            (begc:endc,1:nlevdecomp))       ; this%col_plant_pdemand_vr          (:,:) = nan
    allocate(this%plant_n_uptake_flux             (begc:endc))                    ; this%plant_n_uptake_flux           (:)   = nan
    allocate(this%soil_n_immob_flux               (begc:endc))                    ; this%soil_n_immob_flux	           (:)   = nan
    allocate(this%soil_n_immob_flux_vr            (begc:endc,1:nlevdecomp))       ; this%soil_n_immob_flux_vr          (:,:) = nan
    allocate(this%soil_n_grossmin_flux            (begc:endc))                    ; this%soil_n_grossmin_flux          (:)   = nan
    allocate(this%actual_immob_no3                (begc:endc))                    ; this%actual_immob_no3              (:)   = nan
    allocate(this%actual_immob_nh4                (begc:endc))                    ; this%actual_immob_nh4              (:)   = nan
    allocate(this%smin_no3_to_plant               (begc:endc))                    ; this%smin_no3_to_plant             (:)   = nan
    allocate(this%smin_nh4_to_plant               (begc:endc))                    ; this%smin_nh4_to_plant             (:)   = nan 
    allocate(this%plant_to_litter_nflux           (begc:endc))                    ; this%plant_to_litter_nflux         (:)   = nan
    allocate(this%plant_to_cwd_nflux              (begc:endc))                    ; this%plant_to_cwd_nflux            (:)   = nan
    allocate(this%bgc_npool_ext_inputs_vr         (begc:endc,1:nlevdecomp_full,ndecomp_pools                )) ; this%bgc_npool_ext_inputs_vr          (:,:,:) = nan
    allocate(this%bgc_npool_ext_loss_vr           (begc:endc,1:nlevdecomp_full,ndecomp_pools                )) ; this%bgc_npool_ext_loss_vr            (:,:,:) = nan
    allocate(this%decomp_cascade_ntransfer_vr     (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)) ; this%decomp_cascade_ntransfer_vr      (:,:,:) = nan
    allocate(this%decomp_cascade_sminn_flux_vr    (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)) ; this%decomp_cascade_sminn_flux_vr     (:,:,:) = nan
    allocate(this%m_decomp_npools_to_fire_vr      (begc:endc,1:nlevdecomp_full,1:ndecomp_pools              )) ; this%m_decomp_npools_to_fire_vr       (:,:,:) = nan
    allocate(this%decomp_cascade_ntransfer        (begc:endc,1:ndecomp_cascade_transitions                  )) ; this%decomp_cascade_ntransfer         (:,:)   = nan
    allocate(this%decomp_cascade_sminn_flux       (begc:endc,1:ndecomp_cascade_transitions                  )) ; this%decomp_cascade_sminn_flux        (:,:)   = nan
    allocate(this%m_decomp_npools_to_fire         (begc:endc,1:ndecomp_pools                                )) ; this%m_decomp_npools_to_fire          (:,:)   = nan
    allocate(this%sminn_to_denit_decomp_cascade_vr(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)) ; this%sminn_to_denit_decomp_cascade_vr (:,:,:) = nan
    allocate(this%sminn_to_denit_decomp_cascade   (begc:endc,1:ndecomp_cascade_transitions                  )) ; this%sminn_to_denit_decomp_cascade    (:,:)   = nan
    allocate(this%sminn_to_denit_excess_vr        (begc:endc,1:nlevdecomp_full                              )) ; this%sminn_to_denit_excess_vr         (:,:)   = nan
    allocate(this%sminn_to_denit_excess           (begc:endc                                                )) ; this%sminn_to_denit_excess            (:)     = nan
    allocate(this%sminn_leached_vr                (begc:endc,1:nlevdecomp_full                              )) ; this%sminn_leached_vr                 (:,:)   = nan
    allocate(this%sminn_leached                   (begc:endc                                                )) ; this%sminn_leached                    (:)     = nan
    allocate(this%decomp_npools_leached           (begc:endc,1:ndecomp_pools                                )) ; this%decomp_npools_leached            (:,:)   = nan
    allocate(this%decomp_npools_transport_tendency(begc:endc,1:nlevdecomp_full,1:ndecomp_pools              )) ; this%decomp_npools_transport_tendency (:,:,:) = nan
    allocate(this%decomp_npools_sourcesink        (begc:endc,1:nlevdecomp_full,1:ndecomp_pools              )) ; this%decomp_npools_sourcesink         (:,:,:) = nan
    allocate(this%externaln_to_decomp_npools      (begc:endc,1:nlevdecomp_full, 1:ndecomp_pools             )) ; this%externaln_to_decomp_npools       (:,:,:) = spval
    allocate(this%pmnf_decomp_cascade             (begc:endc,1:nlevdecomp,1:ndecomp_cascade_transitions     )) ; this%pmnf_decomp_cascade              (:,:,:) = nan

  end subroutine col_nf_init

  !------------------------------------------------------------------------
  subroutine col_nf_clean(this)
    !
    ! !ARGUMENTS:
    class(column_nitrogen_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine col_nf_clean
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column phosphorus flux data structure
  !------------------------------------------------------------------------
  subroutine col_pf_init(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_phosphorus_flux) :: this
    integer, intent(in) :: begc,endc
    !
    ! !LOCAL VARIABLES:
    integer           :: k,l
    character(24)     :: fieldname
    character(100)    :: longname
    character(8)      :: vr_suffix
    character(1)      :: aa 
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    integer           :: c
    integer           :: fc                           ! filter indices
    integer           :: num_special_col              ! number of good values in special_col filter
    integer           :: special_col(endc-begc+1)     ! special landunit filter - columns
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of col_pf
    !-----------------------------------------------------------------------
    allocate(this%pdep_to_sminp                    (begc:endc))                   ; this%pdep_to_sminp                 (:)   = nan
    allocate(this%fert_p_to_sminp                  (begc:endc))                   ; this%fert_p_to_sminp               (:)   = nan
    allocate(this%hrv_deadstemp_to_prod10p         (begc:endc))                   ; this%hrv_deadstemp_to_prod10p      (:)   = nan
    allocate(this%hrv_deadstemp_to_prod100p        (begc:endc))                   ; this%hrv_deadstemp_to_prod100p     (:)   = nan
    allocate(this%hrv_cropp_to_prod1p              (begc:endc))                   ; this%hrv_cropp_to_prod1p           (:)   = nan
    allocate(this%sminp_to_plant                   (begc:endc))                   ; this%sminp_to_plant                (:)   = nan
    allocate(this%potential_immob_p                (begc:endc))                   ; this%potential_immob_p             (:)   = nan
    allocate(this%actual_immob_p                   (begc:endc))                   ; this%actual_immob_p                (:)   = nan
    allocate(this%gross_pmin                       (begc:endc))                   ; this%gross_pmin                    (:)   = nan
    allocate(this%net_pmin                         (begc:endc))                   ; this%net_pmin                      (:)   = nan
    allocate(this%supplement_to_sminp              (begc:endc))                   ; this%supplement_to_sminp           (:)   = nan
    allocate(this%prod1p_loss                      (begc:endc))                   ; this%prod1p_loss                   (:)   = nan
    allocate(this%prod10p_loss                     (begc:endc))                   ; this%prod10p_loss                  (:)   = nan
    allocate(this%prod100p_loss                    (begc:endc))                   ; this%prod100p_loss                 (:)   = nan
    allocate(this%product_ploss                    (begc:endc))                   ; this%product_ploss                 (:)   = nan
    allocate(this%pinputs                          (begc:endc))                   ; this%pinputs                       (:)   = nan
    allocate(this%poutputs                         (begc:endc))                   ; this%poutputs                      (:)   = nan
    allocate(this%fire_ploss                       (begc:endc))                   ; this%fire_ploss                    (:)   = nan
    allocate(this%fire_decomp_ploss                (begc:endc))                   ; this%fire_decomp_ploss             (:)   = nan
    allocate(this%fire_ploss_p2c                   (begc:endc))                   ; this%fire_ploss_p2c                (:)   = nan
    allocate(this%som_p_leached                    (begc:endc))                   ; this%som_p_leached                 (:)   = nan
    allocate(this%m_p_to_litr_met_fire             (begc:endc,1:nlevdecomp_full)) ; this%m_p_to_litr_met_fire          (:,:) = nan
    allocate(this%m_p_to_litr_cel_fire             (begc:endc,1:nlevdecomp_full)) ; this%m_p_to_litr_cel_fire          (:,:) = nan
    allocate(this%m_p_to_litr_lig_fire             (begc:endc,1:nlevdecomp_full)) ; this%m_p_to_litr_lig_fire          (:,:) = nan
    allocate(this%potential_immob_p_vr             (begc:endc,1:nlevdecomp_full)) ; this%potential_immob_p_vr          (:,:) = nan
    allocate(this%actual_immob_p_vr                (begc:endc,1:nlevdecomp_full)) ; this%actual_immob_p_vr             (:,:) = nan
    allocate(this%sminp_to_plant_vr                (begc:endc,1:nlevdecomp_full)) ; this%sminp_to_plant_vr             (:,:) = nan
    allocate(this%supplement_to_sminp_vr           (begc:endc,1:nlevdecomp_full)) ; this%supplement_to_sminp_vr        (:,:) = nan
    allocate(this%gross_pmin_vr                    (begc:endc,1:nlevdecomp_full)) ; this%gross_pmin_vr                 (:,:) = nan
    allocate(this%net_pmin_vr                      (begc:endc,1:nlevdecomp_full)) ; this%net_pmin_vr                   (:,:) = nan
    allocate(this%biochem_pmin_to_ecosysp_vr       (begc:endc,1:nlevdecomp_full)) ; this%biochem_pmin_to_ecosysp_vr    (:,:) = nan
    allocate(this%biochem_pmin_ppools_vr           (begc:endc,1:nlevdecomp_full,1:ndecomp_pools)) ; this%biochem_pmin_ppools_vr      (:,:,:) = nan
    allocate(this%biochem_pmin_vr                  (begc:endc,1:nlevdecomp_full)) ; this%biochem_pmin_vr               (:,:) = nan
    allocate(this%biochem_pmin                     (begc:endc))                   ; this%biochem_pmin                  (:)   = nan
    allocate(this%dwt_slash_pflux                  (begc:endc))                   ; this%dwt_slash_pflux               (:)   = nan
    allocate(this%dwt_conv_pflux                   (begc:endc))                   ; this%dwt_conv_pflux                (:)   = nan
    allocate(this%dwt_prod10p_gain                 (begc:endc))                   ; this%dwt_prod10p_gain              (:)   = nan
    allocate(this%dwt_prod100p_gain                (begc:endc))                   ; this%dwt_prod100p_gain             (:)   = nan
    allocate(this%dwt_ploss                        (begc:endc))                   ; this%dwt_ploss                     (:)   = nan
    allocate(this%wood_harvestp                    (begc:endc))                   ; this%wood_harvestp                 (:)   = nan
    allocate(this%dwt_frootp_to_litr_met_p         (begc:endc,1:nlevdecomp_full)) ; this%dwt_frootp_to_litr_met_p      (:,:) = nan
    allocate(this%dwt_frootp_to_litr_cel_p         (begc:endc,1:nlevdecomp_full)) ; this%dwt_frootp_to_litr_cel_p      (:,:) = nan
    allocate(this%dwt_frootp_to_litr_lig_p         (begc:endc,1:nlevdecomp_full)) ; this%dwt_frootp_to_litr_lig_p      (:,:) = nan
    allocate(this%dwt_livecrootp_to_cwdp           (begc:endc,1:nlevdecomp_full)) ; this%dwt_livecrootp_to_cwdp        (:,:) = nan
    allocate(this%dwt_deadcrootp_to_cwdp           (begc:endc,1:nlevdecomp_full)) ; this%dwt_deadcrootp_to_cwdp        (:,:) = nan
    allocate(this%decomp_cascade_ptransfer_vr      (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions )) ; this%decomp_cascade_ptransfer_vr  (:,:,:) = nan
    allocate(this%decomp_cascade_sminp_flux_vr     (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions )) ; this%decomp_cascade_sminp_flux_vr (:,:,:) = nan
    allocate(this%m_decomp_ppools_to_fire_vr       (begc:endc,1:nlevdecomp_full,1:ndecomp_pools               )) ; this%m_decomp_ppools_to_fire_vr   (:,:,:) = nan
    allocate(this%decomp_cascade_ptransfer         (begc:endc,1:ndecomp_cascade_transitions                   )) ; this%decomp_cascade_ptransfer     (:,:)   = nan
    allocate(this%decomp_cascade_sminp_flux        (begc:endc,1:ndecomp_cascade_transitions                   )) ; this%decomp_cascade_sminp_flux    (:,:)   = nan
    allocate(this%m_decomp_ppools_to_fire          (begc:endc,1:ndecomp_pools ))  ; this%m_decomp_ppools_to_fire       (:,:) = nan
    allocate(this%phenology_p_to_litr_met_p        (begc:endc,1:nlevdecomp_full)) ; this%phenology_p_to_litr_met_p     (:,:) = nan
    allocate(this%phenology_p_to_litr_cel_p        (begc:endc,1:nlevdecomp_full)) ; this%phenology_p_to_litr_cel_p     (:,:) = nan
    allocate(this%phenology_p_to_litr_lig_p        (begc:endc,1:nlevdecomp_full)) ; this%phenology_p_to_litr_lig_p     (:,:) = nan
    allocate(this%gap_mortality_p_to_litr_met_p    (begc:endc,1:nlevdecomp_full)) ; this%gap_mortality_p_to_litr_met_p (:,:) = nan
    allocate(this%gap_mortality_p_to_litr_cel_p    (begc:endc,1:nlevdecomp_full)) ; this%gap_mortality_p_to_litr_cel_p (:,:) = nan
    allocate(this%gap_mortality_p_to_litr_lig_p    (begc:endc,1:nlevdecomp_full)) ; this%gap_mortality_p_to_litr_lig_p (:,:) = nan
    allocate(this%gap_mortality_p_to_cwdp          (begc:endc,1:nlevdecomp_full)) ; this%gap_mortality_p_to_cwdp       (:,:) = nan
    allocate(this%fire_mortality_p_to_cwdp         (begc:endc,1:nlevdecomp_full)) ; this%fire_mortality_p_to_cwdp      (:,:) = nan
    allocate(this%harvest_p_to_litr_met_p          (begc:endc,1:nlevdecomp_full)) ; this%harvest_p_to_litr_met_p       (:,:) = nan
    allocate(this%harvest_p_to_litr_cel_p          (begc:endc,1:nlevdecomp_full)) ; this%harvest_p_to_litr_cel_p       (:,:) = nan
    allocate(this%harvest_p_to_litr_lig_p          (begc:endc,1:nlevdecomp_full)) ; this%harvest_p_to_litr_lig_p       (:,:) = nan
    allocate(this%harvest_p_to_cwdp                (begc:endc,1:nlevdecomp_full)) ; this%harvest_p_to_cwdp             (:,:) = nan
    allocate(this%primp_to_labilep_vr              (begc:endc,1:nlevdecomp_full)) ; this%primp_to_labilep_vr           (:,:) = nan
    allocate(this%primp_to_labilep                 (begc:endc))                   ; this%primp_to_labilep              (:)   = nan
    allocate(this%labilep_to_secondp_vr            (begc:endc,1:nlevdecomp_full)) ; this%labilep_to_secondp_vr         (:,:) = nan
    allocate(this%labilep_to_secondp               (begc:endc))                   ; this%labilep_to_secondp            (:)   = nan
    allocate(this%secondp_to_labilep_vr            (begc:endc,1:nlevdecomp_full)) ; this%secondp_to_labilep_vr         (:,:) = nan
    allocate(this%secondp_to_labilep               (begc:endc))                   ; this%secondp_to_labilep            (:)   = nan
    allocate(this%secondp_to_occlp_vr              (begc:endc,1:nlevdecomp_full)) ; this%secondp_to_occlp_vr           (:,:) = nan
    allocate(this%secondp_to_occlp                 (begc:endc))                   ; this%secondp_to_occlp              (:)   = nan
    allocate(this%sminp_leached_vr                 (begc:endc,1:nlevdecomp_full)) ; this%sminp_leached_vr              (:,:) = nan
    allocate(this%sminp_leached                    (begc:endc))                   ; this%sminp_leached                 (:)   = nan
    allocate(this%decomp_ppools_leached            (begc:endc,1:ndecomp_pools  )) ; this%decomp_ppools_leached         (:,:) = nan
    allocate(this%decomp_ppools_transport_tendency (begc:endc,1:nlevdecomp_full,1:ndecomp_pools           )) ; this%decomp_ppools_transport_tendency (:,:,:) = nan 
    allocate(this%decomp_ppools_sourcesink         (begc:endc,1:nlevdecomp_full,1:ndecomp_pools           )) ; this%decomp_ppools_sourcesink         (:,:,:) = nan
    allocate(this%adsorb_to_labilep_vr             (begc:endc,1:nlevdecomp_full)) ; this%adsorb_to_labilep_vr          (:,:) = nan
    allocate(this%desorb_to_solutionp_vr           (begc:endc,1:nlevdecomp_full)) ; this%desorb_to_solutionp_vr        (:,:) = nan
    allocate(this%adsorb_to_labilep                (begc:endc))                   ; this%adsorb_to_labilep             (:)   = nan
    allocate(this%desorb_to_solutionp              (begc:endc))                   ; this%desorb_to_solutionp           (:)   = nan
    allocate(this%pmpf_decomp_cascade              (begc:endc,1:nlevdecomp,1:ndecomp_cascade_transitions  )) ; this%pmpf_decomp_cascade              (:,:,:) = nan
    allocate(this%plant_p_uptake_flux              (begc:endc))                   ; this%plant_p_uptake_flux           (:)   = nan
    allocate(this%soil_p_immob_flux                (begc:endc))                   ; this%soil_p_immob_flux             (:)   = nan
    allocate(this%soil_p_immob_flux_vr             (begc:endc,1:nlevdecomp_full)) ; this%soil_p_immob_flux_vr          (:,:) = nan
    allocate(this%soil_p_grossmin_flux             (begc:endc))                   ; this%soil_p_grossmin_flux          (:)   = nan
    allocate(this%smin_p_to_plant                  (begc:endc))                   ; this%smin_p_to_plant               (:)   = nan
    allocate(this%plant_to_litter_pflux            (begc:endc))                   ; this%plant_to_litter_pflux         (:)   = nan
    allocate(this%plant_to_cwd_pflux               (begc:endc))                   ; this%plant_to_cwd_pflux            (:)   = nan
    allocate(this%plant_pdemand                    (begc:endc))                   ; this%plant_pdemand                 (:)   = nan
    allocate(this%plant_pdemand_vr                 (begc:endc,1:nlevdecomp_full)) ; this%plant_pdemand_vr              (:,:) = nan
    allocate(this%externalp_to_decomp_ppools       (begc:endc,1:nlevdecomp_full, 1:ndecomp_pools          )) ;    this%externalp_to_decomp_ppools    (:,:,:) = spval
    allocate(this%externalp_to_decomp_delta        (begc:endc))                   ; this%externalp_to_decomp_delta     (:)   = spval
    allocate(this%sminp_net_transport_vr           (begc:endc,1:nlevdecomp_full)) ; this%sminp_net_transport_vr        (:,:) = spval
    allocate(this%sminp_net_transport_delta        (begc:endc))                   ; this%sminp_net_transport_delta     (:)   = spval

  end subroutine col_pf_init
    
  !------------------------------------------------------------------------
  subroutine col_pf_clean(this)
    !
    ! !ARGUMENTS:
    class(column_phosphorus_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine col_pf_clean
  
    !------------------------------------------------------------------------
    
end module ColumnDataType


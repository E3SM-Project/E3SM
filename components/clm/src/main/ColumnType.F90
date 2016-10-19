module ColumnType
  
  ! -------------------------------------------------------- 
  ! ALM sub-grid hierarchy:
  ! Define Soil column unit data types, with Init and Clean for each
  ! -------------------------------------------------------- 
  ! 10 Oct 2016, DW  

  ! Moved from ColumnType.F90
  ! -------------------------------------------------------- 
  ! column types can have values of
  ! -------------------------------------------------------- 
  !   1  => (istsoil)          soil (vegetated or bare soil)
  !   2  => (istcrop)          crop (only for crop configuration)
  !   3  => (istice)           land ice
  !   4  => (istice_mec)       land ice (multiple elevation classes)   
  !   5  => (istdlak)          deep lake
  !   6  => (istwet)           wetland
  !   71 => (icol_roof)        urban roof
  !   72 => (icol_sunwall)     urban sunwall
  !   73 => (icol_shadewall)   urban shadewall
  !   74 => (icol_road_imperv) urban impervious road
  !   75 => (icol_road_perv)   urban pervious road
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use clm_varcon     , only : ispval, spval
  use clm_varpar     , only : nlevsno, nlevgrnd, nlevlak


!  moved from CNCarbonstate.F90
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use clm_varpar             , only : ndecomp_cascade_transitions, ndecomp_pools, nlevcan
  use clm_varpar             , only : nlevdecomp_full, crop_prog, nlevdecomp
  use clm_varcon             , only : spval, ispval, dzsoi_decomp, zisoi, zsoi
  use landunit_varcon        , only : istcrop 
  use clm_varctl             , only : iulog, use_vertsoilc, use_cndv, spinup_state 
  use decompMod              , only : bounds_type
!  use CNStateType            , only : cnstate_type
!  use pftvarcon              , only : npcropmin
  use CNDecompCascadeConType , only : decomp_cascade_con
  use EcophysConType         , only : ecophyscon
  use abortutils             , only : endrun
  use spmdMod                , only : masterproc 
!  use subgridAveMod          , only : p2c
  use LandunitType           , only : lun                
! use ColumnType             , only : col                
  use PatchType              , only : pft
  use clm_varctl             , only : nu_com

! Moved from WaterStateType.F90  
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use decompMod      , only : bounds_type
  use clm_varctl     , only : use_vancouver, use_mexicocity, use_cn, iulog
  use clm_varpar     , only : nlevgrnd, nlevurb, nlevsno   
  use clm_varcon     , only : spval
  use LandunitType   , only : lun                
!  use ColumnType     , only : col   


  implicit none
  save
  private

  ! sub-grid geospatial and physical properties defined at the soil_column level
  ! migrate variables list from ColumnType.F90

  type, public :: soilcol_properties
     ! g/l/c/p hierarchy, local g/l/c/p cells only
     integer , pointer :: landunit             (:)   ! index into landunit level quantities
     real(r8), pointer :: wtlunit              (:)   ! weight (relative to landunit)
     integer , pointer :: gridcell             (:)   ! index into gridcell level quantities
     real(r8), pointer :: wtgcell              (:)   ! weight (relative to gridcell)
     integer , pointer :: pfti                 (:)   ! beginning pft index for each column
     integer , pointer :: pftf                 (:)   ! ending pft index for each column
     integer , pointer :: npfts                (:)   ! number of patches for each column

     ! topological mapping functionality
     integer , pointer :: itype                (:)   ! column type
     logical , pointer :: active               (:)   ! true=>do computations on this column 

     ! topography
     real(r8), pointer :: glc_topo             (:)   ! surface elevation (m)
     real(r8), pointer :: micro_sigma          (:)   ! microtopography pdf sigma (m)
     real(r8), pointer :: n_melt               (:)   ! SCA shape parameter
     real(r8), pointer :: topo_slope           (:)   ! gridcell topographic slope
     real(r8), pointer :: topo_std             (:)   ! gridcell elevation standard deviation

     ! vertical levels
     integer , pointer :: snl                  (:)   ! number of snow layers
     real(r8), pointer :: dz                   (:,:) ! layer thickness (m)  (-nlevsno+1:nlevgrnd) 
     real(r8), pointer :: z                    (:,:) ! layer depth (m) (-nlevsno+1:nlevgrnd) 
     real(r8), pointer :: zi                   (:,:) ! interface level below a "z" level (m) (-nlevsno+0:nlevgrnd) 
     real(r8), pointer :: zii                  (:)   ! convective boundary height [m]
     real(r8), pointer :: dz_lake              (:,:) ! lake layer thickness (m)  (1:nlevlak)
     real(r8), pointer :: z_lake               (:,:) ! layer depth for lake (m)
     real(r8), pointer :: lakedepth            (:)   ! variable lake depth (m)                             

   contains

     procedure, public :: Init => init_col_pp
     procedure, public :: Clean => clean_col_pp

  end type soilcol_properties


  type, public :: soilcol_energy_state
 (***DW***: No energy state Type?)  
  contains
     procedure, public :: Init => init_col_es
     procedure, public :: Clean => clean_col_es
  end type soilcol_energy_state


  type, public :: soilcol_water_state
     logical , pointer :: do_capsnow_col         (:)   ! col true => do snow capping
     real(r8), pointer :: snow_depth_col         (:)   ! col snow height of snow covered area (m)
     real(r8), pointer :: snow_persistence_col   (:)   ! col length of time that ground has had non-zero snow thickness (sec)
     real(r8), pointer :: snowdp_col             (:)   ! col gridcell averaged snow height (m)
     real(r8), pointer :: snowice_col            (:)   ! col average snow ice lens
     real(r8), pointer :: snowliq_col            (:)   ! col average snow liquid water
     real(r8), pointer :: int_snow_col           (:)   ! col integrated snowfall (mm H2O)
     real(r8), pointer :: snow_layer_unity_col   (:,:) ! value 1 for each snow layer, used for history diagnostics
     real(r8), pointer :: bw_col                 (:,:) ! col partial density of water in the snow pack (ice + liquid) [kg/m3] 
     real(r8), pointer :: finundated_col         (:)   ! fraction of column that is inundated, this is for bgc caclulation in betr

     real(r8), pointer :: smp_l_col              (:,:) ! col liquid phase soil matric potential, mm
     real(r8), pointer :: h2osno_col             (:)   ! col snow water (mm H2O)
     real(r8), pointer :: h2osno_old_col         (:)   ! col snow mass for previous time step (kg/m2) (new)
     real(r8), pointer :: h2osoi_liq_col         (:,:) ! col liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)    
     real(r8), pointer :: h2osoi_ice_col         (:,:) ! col ice lens (kg/m2) (new) (-nlevsno+1:nlevgrnd)    
     real(r8), pointer :: h2osoi_liqice_10cm_col (:)   ! col liquid water + ice lens in top 10cm of soil (kg/m2)
     real(r8), pointer :: h2osoi_vol_col         (:,:) ! col volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevgrnd)
     real(r8), pointer :: air_vol_col            (:,:) ! col air filled porosity     
     real(r8), pointer :: h2osoi_liqvol_col      (:,:) ! col volumetric liquid water content (v/v)
     real(r8), pointer :: h2osoi_icevol_col      (:,:) ! col volumetric ice content (v/v)     
     real(r8), pointer :: h2ocan_patch           (:)   ! patch canopy water (mm H2O)
     real(r8), pointer :: h2ocan_col             (:)   ! col canopy water (mm H2O)
     real(r8), pointer :: h2osfc_col             (:)   ! col surface water (mm H2O)
     real(r8), pointer :: swe_old_col            (:,:) ! col initial snow water
     real(r8), pointer :: liq1_grc               (:)   ! grc initial gridcell total h2o liq content
     real(r8), pointer :: liq2_grc               (:)   ! grc post land cover change total liq content
     real(r8), pointer :: ice1_grc               (:)   ! grc initial gridcell total h2o ice content
     real(r8), pointer :: ice2_grc               (:)   ! grc post land cover change total ice content
     real(r8), pointer :: tws_grc                (:)   ! grc total water storage (mm H2O)

     real(r8), pointer :: snw_rds_col            (:,:) ! col snow grain radius (col,lyr)    [m^-6, microns]
     real(r8), pointer :: snw_rds_top_col        (:)   ! col snow grain radius (top layer)  [m^-6, microns]
     real(r8), pointer :: h2osno_top_col         (:)   ! col top-layer mass of snow  [kg]
     real(r8), pointer :: sno_liq_top_col        (:)   ! col snow liquid water fraction (mass), top layer  [fraction]

     real(r8), pointer :: q_ref2m_patch          (:)   ! patch 2 m height surface specific humidity (kg/kg)
     real(r8), pointer :: rh_ref2m_patch         (:)   ! patch 2 m height surface relative humidity (%)
     real(r8), pointer :: rh_ref2m_r_patch       (:)   ! patch 2 m height surface relative humidity - rural (%)
     real(r8), pointer :: rh_ref2m_u_patch       (:)   ! patch 2 m height surface relative humidity - urban (%)
     real(r8), pointer :: rh_af_patch            (:)   ! patch fractional humidity of canopy air (dimensionless) ! private
     real(r8), pointer :: qg_snow_col            (:)   ! col ground specific humidity [kg/kg]
     real(r8), pointer :: qg_soil_col            (:)   ! col ground specific humidity [kg/kg]
     real(r8), pointer :: qg_h2osfc_col          (:)   ! col ground specific humidity [kg/kg]
     real(r8), pointer :: qg_col                 (:)   ! col ground specific humidity [kg/kg]
     real(r8), pointer :: dqgdT_col              (:)   ! col d(qg)/dT
     real(r8), pointer :: qaf_lun                (:)   ! lun urban canopy air specific humidity (kg/kg)

     ! Fractions
     real(r8), pointer :: frac_sno_col           (:)   ! col fraction of ground covered by snow (0 to 1)
     real(r8), pointer :: frac_sno_eff_col       (:)   ! col fraction of ground covered by snow (0 to 1)
     real(r8), pointer :: frac_iceold_col        (:,:) ! col fraction of ice relative to the tot water (new) (-nlevsno+1:nlevgrnd) 
     real(r8), pointer :: frac_h2osfc_col        (:)   ! col fractional area with surface water greater than zero
     real(r8), pointer :: wf_col                 (:)   ! col soil water as frac. of whc for top 0.05 m (0-1) 
     real(r8), pointer :: wf2_col                (:)   ! col soil water as frac. of whc for top 0.17 m (0-1) 
     real(r8), pointer :: fwet_patch             (:)   ! patch canopy fraction that is wet (0 to 1)
     real(r8), pointer :: fdry_patch             (:)   ! patch canopy fraction of foliage that is green and dry [-] (new)

     ! Balance Checks

     real(r8), pointer :: begwb_patch            (:)   ! water mass begining of the time step
     real(r8), pointer :: begwb_col              (:)   ! water mass begining of the time step
     real(r8), pointer :: endwb_patch            (:)   ! water mass end of the time step
     real(r8), pointer :: endwb_col              (:)   ! water mass end of the time step
     real(r8), pointer :: errh2o_patch           (:)   ! water conservation error (mm H2O)
     real(r8), pointer :: errh2o_col             (:)   ! water conservation error (mm H2O)
     real(r8), pointer :: errh2osno_col          (:)   ! snow water conservation error(mm H2O)

     ! For VSFM
     real(r8), pointer :: vsfm_fliq_col_1d       (:)   ! fraction of liquid saturation for VSFM [-]
     real(r8), pointer :: vsfm_sat_col_1d        (:)   ! liquid saturation from VSFM [-]
     real(r8), pointer :: vsfm_mass_col_1d       (:)   ! liquid mass per unit area from VSFM [kg H2O/m^2]
     real(r8), pointer :: vsfm_smpl_col_1d       (:)   ! 1D soil matrix potential liquid from VSFM [m]
     real(r8), pointer :: vsfm_soilp_col_1d      (:)   ! 1D soil liquid pressure from VSFM [Pa]
     real(r8), pointer :: soilp_col              (:,:) ! col soil pressure (Pa)

  contains
     procedure, public :: Init => init_col_ws
     procedure, public :: Clean => clean_col_ws
!     procedure , public  :: Restart => restart_col_ws
     procedure , public  :: Reset => reset_col_ws  
     procedure , private :: InitAllocate => initallocate_col_ws
!     procedure , private :: InitHistory => inithistory_col_ws
!     procedure , private :: InitCold => initcold_col_ws

  end type soilcol_water_state


  type, public :: soilcol_carbon_state
!   Moved from CNCarbonStateType.F90
     real(r8), pointer :: rootc_col                (:)     ! col (gC/m2) root carbon at column level (fire)
     real(r8), pointer :: totvegc_col              (:)     ! col (gC/m2) column-level totvegc (fire)
     real(r8), pointer :: leafc_col                (:)     ! col (gC/m2) column-level leafc (fire)
     real(r8), pointer :: deadstemc_col            (:)     ! col (gC/m2) column-level deadstemc (fire)
     real(r8), pointer :: fuelc_col                (:)     ! col fuel avalability factor for Reg.C (0-1)
     real(r8), pointer :: fuelc_crop_col           (:)     ! col fuel avalability factor for Reg.A (0-1)

     ! all c pools involved in decomposition
     real(r8), pointer :: decomp_cpools_vr_col    (:,:,:)  ! col (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
     real(r8), pointer :: ctrunc_vr_col           (:,:)    ! col (gC/m3) vertically-resolved column-level sink for C truncation

     ! pools for dynamic landcover
     real(r8), pointer :: frootc_col               (:)     ! col (gC/m2) column-level C pool for fine root
     real(r8), pointer :: seedc_col                (:)     ! col (gC/m2) column-level pool for seeding new Patches
     real(r8), pointer :: prod10c_col              (:)     ! col (gC/m2) wood product C pool, 10-year lifespan
     real(r8), pointer :: prod100c_col             (:)     ! col (gC/m2) wood product C pool, 100-year lifespan
     real(r8), pointer :: totprodc_col             (:)     ! col (gC/m2) total wood product C
     ! pools for crop harvest
     real(r8), pointer :: prod1c_col               (:)     ! col (gC/m2) crop product C pool, 1-year lifespan

     ! summary (diagnostic) state variables, not involved in mass balance
     real(r8), pointer :: totvegcc_col             (:)     ! (gC/m2) total vegetation carbon, excluding cpool averaged to column (p2c)
     real(r8), pointer :: totpftc_col              (:)     ! (gC/m2) total patch-level carbon, including cpool averaged to column (p2c)
     real(r8), pointer :: decomp_cpools_1m_col     (:,:)   ! col (gC/m2)  Diagnostic: decomposing (litter, cwd, soil) c pools to 1 meter
     real(r8), pointer :: decomp_cpools_col        (:,:)   ! col (gC/m2)  decomposing (litter, cwd, soil) c pools
     real(r8), pointer :: cwdc_col                 (:)     ! col (gC/m2) Diagnostic: coarse woody debris C
     real(r8), pointer :: ctrunc_col               (:)     ! col (gC/m2) column-level sink for C truncation
     real(r8), pointer :: totlitc_col              (:)     ! col (gC/m2) total litter carbon
     real(r8), pointer :: totsomc_col              (:)     ! col (gC/m2) total soil organic matter carbon
     real(r8), pointer :: totlitc_1m_col           (:)     ! col (gC/m2) total litter carbon to 1 meter
     real(r8), pointer :: totsomc_1m_col           (:)     ! col (gC/m2) total soil organic matter carbon to 1 meter
     real(r8), pointer :: totecosysc_col           (:)     ! col (gC/m2) total ecosystem carbon, incl veg but excl cpool
     real(r8), pointer :: totcolc_col              (:)     ! col (gC/m2) total column carbon, incl veg and cpool
     real(r8), pointer :: totabgc_col              (:)     ! col (gC/m2) total column above ground carbon, excluding som 

     ! Balance checks
     real(r8), pointer :: begcb_col                (:)     ! patch carbon mass, beginning of time step (gC/m**2)
     real(r8), pointer :: endcb_col                (:)     ! patch carbon mass, end of time step (gC/m**2)
     real(r8), pointer :: errcb_col                (:)     ! patch carbon balance error for the timestep (gC/m**2)
     
     real(r8), pointer :: totpftc_beg_col(:)
     real(r8), pointer :: cwdc_beg_col(:)
     real(r8), pointer :: totlitc_beg_col(:)
     real(r8), pointer :: totsomc_beg_col(:)
     
     real(r8), pointer :: totpftc_end_col(:)
     real(r8), pointer :: cwdc_end_col(:)
     real(r8), pointer :: totlitc_end_col(:)
     real(r8), pointer :: totsomc_end_col(:)

  contains
     procedure, public :: Init => init_col_cs
     procedure, public :: Clean => clean_col_cs
     procedure , public  :: SetValues => setvalues_col_cs
!     procedure , public  :: ZeroDWT => zerodwt_col_cs  !zerodwt only at patch level
!     procedure , public  :: Restart => restart_col_cs
     procedure , public  :: Summary => summary_col_cs  !moved summary at patch level
     procedure , private :: InitAllocate => initallocate_col_cs
!     procedure , private :: InitHistory => inithistory_col_cs
!     procedure , private :: InitCold => initcold_col_cs
  end type soilcol_carbon_state

  type, public :: soilcol_nitrogen_state

  contains
     procedure, public :: Init => init_col_ns
     procedure, public :: Clean => clean_col_ns
  end type soilcol_nitrogen_state

  type, public :: soilcol_phosphorus_state

  contains
     procedure, public :: Init => init_col_ps
     procedure, public :: Clean => clean_col_ps
  end type soilcol_phosphorus_state


  ! declare the public instances of soilcolumn types
  type(soilcol_properties) , public, target :: col_pp
  type(soilcol_energy_state), public, target :: col_es
  type(soilcol_water_state), public, target :: col_ws
  type(soilcol_carbon_state) , public, target :: col_cp
  type(soilcol_nitrogen_state), public, target :: col_ns
  type(soilcol_phosphorus_state), public, target :: col_ps

  type(column_type), public, target :: col !column data structure (soil/snow/canopy columns)
  !------------------------------------------------------------------------

  contains 

  subroutine init_col_pp(this, begc, endc)
    class(soilcol_properties) :: this
    integer, intent(in) :: begc   ! beginning soil column index
    integer, intent(in) :: endc   ! ending soil column index

    ! The following is set in initGridCellsMod  (***DW***: Need to considere Topounit?)
    allocate(this%gridcell    (begc:endc))                     ; this%gridcell    (:)   = ispval
    allocate(this%wtgcell     (begc:endc))                     ; this%wtgcell     (:)   = nan
    allocate(this%landunit    (begc:endc))                     ; this%landunit    (:)   = ispval
    allocate(this%wtlunit     (begc:endc))                     ; this%wtlunit     (:)   = nan
    allocate(this%pfti        (begc:endc))                     ; this%pfti        (:)   = ispval
    allocate(this%pftf        (begc:endc))                     ; this%pftf        (:)   = ispval
    allocate(this%npfts       (begc:endc))                     ; this%npfts       (:)   = ispval
    allocate(this%itype       (begc:endc))                     ; this%itype       (:)   = ispval
    allocate(this%active      (begc:endc))                     ; this%active      (:)   = .false.

    ! The following is set in initVerticalMod
    allocate(this%snl         (begc:endc))                     ; this%snl         (:)   = ispval  !* cannot be averaged up
    allocate(this%dz          (begc:endc,-nlevsno+1:nlevgrnd)) ; this%dz          (:,:) = nan
    allocate(this%z           (begc:endc,-nlevsno+1:nlevgrnd)) ; this%z           (:,:) = nan
    allocate(this%zi          (begc:endc,-nlevsno+0:nlevgrnd)) ; this%zi          (:,:) = nan
    allocate(this%zii         (begc:endc))                     ; this%zii         (:)   = nan
    allocate(this%lakedepth   (begc:endc))                     ; this%lakedepth   (:)   = spval  
    allocate(this%dz_lake     (begc:endc,nlevlak))             ; this%dz_lake     (:,:) = nan
    allocate(this%z_lake      (begc:endc,nlevlak))             ; this%z_lake      (:,:) = nan

    allocate(this%glc_topo    (begc:endc))                     ; this%glc_topo    (:)   = nan
    allocate(this%micro_sigma (begc:endc))                     ; this%micro_sigma (:)   = nan
    allocate(this%n_melt      (begc:endc))                     ; this%n_melt      (:)   = nan 
    allocate(this%topo_slope  (begc:endc))                     ; this%topo_slope  (:)   = nan
    allocate(this%topo_std    (begc:endc))                     ; this%topo_std    (:)   = nan
  
  subroutine clean_col_pp(this)
    class(soilcol_properties) :: this
  
! !ARGUMENTS:
    class(column_type) :: this
    !------------------------------------------------------------------------

    deallocate(this%gridcell   )
    deallocate(this%wtgcell    )
    deallocate(this%landunit   )
    deallocate(this%wtlunit    )
    deallocate(this%pfti       )
    deallocate(this%pftf       )
    deallocate(this%npfts      )
    deallocate(this%itype      )
    deallocate(this%active     )
    deallocate(this%snl        )
    deallocate(this%dz         )
    deallocate(this%z          )
    deallocate(this%zi         )
    deallocate(this%zii        )
    deallocate(this%lakedepth  )
    deallocate(this%dz_lake    )
    deallocate(this%z_lake     )
    deallocate(this%glc_topo   )
    deallocate(this%micro_sigma)
    deallocate(this%n_melt     )
    deallocate(this%topo_slope )
    deallocate(this%topo_std   )

  end subroutine clean_col_pp


  !------------------------------------------------------------------------
  !!!!!  subroutine Init(this, bounds, carbon_type, ratio, c12_carbonstate_vars)
  subroutine init_col_ns(this, bounds, ratio, c12_carbonstate_vars)
    class(soilcol_carbonstate_type)               :: this
    type(bounds_type)      , intent(in)           :: bounds  
  !!!!  character(len=3)       , intent(in)           :: carbon_type
    real(r8)               , intent(in)           :: ratio
    type(carbonstate_type) , intent(in), optional :: c12_carbonstate_vars

    call this%InitAllocate ( bounds)
  !!!!!    call this%InitHistory ( bounds, carbon_type)
    call this%InitHistory( bounds)
    if (present(c12_carbonstate_vars)) then
       call this%InitCold  ( bounds, ratio, c12_carbonstate_vars)
    else
       call this%InitCold  ( bounds, ratio)
    end if

  end subroutine init_col_ns


  !------------------------------------------------------------------------
  subroutine initallocate_col_ns(this, bounds)
    !
    ! !ARGUMENTS:
    class (carbonstate_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
!    integer           :: begp,endp
    integer           :: begc,endc
    !------------------------------------------------------------------------

!    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc

    allocate(this%cwdc_col                 (begc :endc))                   ;     this%cwdc_col                 (:)   = nan
    allocate(this%ctrunc_col               (begc :endc))                   ;     this%ctrunc_col               (:)   = nan
    allocate(this%ctrunc_vr_col            (begc :endc,1:nlevdecomp_full)) ;     this%ctrunc_vr_col            (:,:) = nan
    allocate(this%seedc_col                (begc :endc))                   ;     this%seedc_col                (:)   = nan
    allocate(this%prod10c_col              (begc :endc))                   ;     this%prod10c_col              (:)   = nan
    allocate(this%prod100c_col             (begc :endc))                   ;     this%prod100c_col             (:)   = nan
    allocate(this%prod1c_col               (begc :endc))                   ;     this%prod1c_col               (:)   = nan
    allocate(this%totprodc_col             (begc :endc))                   ;     this%totprodc_col             (:)   = nan
    allocate(this%totlitc_col              (begc :endc))                   ;     this%totlitc_col              (:)   = nan
    allocate(this%totsomc_col              (begc :endc))                   ;     this%totsomc_col              (:)   = nan
    allocate(this%totlitc_1m_col           (begc :endc))                   ;     this%totlitc_1m_col           (:)   = nan
    allocate(this%totsomc_1m_col           (begc :endc))                   ;     this%totsomc_1m_col           (:)   = nan
    allocate(this%totecosysc_col           (begc :endc))                   ;     this%totecosysc_col           (:)   = nan
    allocate(this%totcolc_col              (begc :endc))                   ;     this%totcolc_col              (:)   = nan
    allocate(this%rootc_col                (begc :endc))                   ;     this%rootc_col                (:)   = nan
    allocate(this%totvegc_col              (begc :endc))                   ;     this%totvegc_col              (:)   = nan
    allocate(this%leafc_col                (begc :endc))                   ;     this%leafc_col                (:)   = nan
    allocate(this%deadstemc_col            (begc :endc))                   ;     this%deadstemc_col            (:)   = nan
    allocate(this%fuelc_col                (begc :endc))                   ;     this%fuelc_col                (:)   = nan
    allocate(this%fuelc_crop_col           (begc :endc))                   ;     this%fuelc_crop_col           (:)   = nan
    allocate(this%decomp_cpools_col        (begc :endc,1:ndecomp_pools))   ;     this%decomp_cpools_col        (:,:) = nan
    allocate(this%decomp_cpools_1m_col     (begc :endc,1:ndecomp_pools))   ;     this%decomp_cpools_1m_col     (:,:) = nan
    allocate(this%totpftc_col              (begc :endc))                   ;     this%totpftc_col              (:)   = nan
    allocate(this%totvegc_col              (begc :endc))                   ;     this%totvegc_col              (:)   = nan

    allocate(this%totabgc_col              (begc :endc))                   ;     this%totabgc_col              (:)   = nan
    allocate(this%decomp_cpools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))  
    this%decomp_cpools_vr_col(:,:,:)= nan

    allocate(this%begcb_col   (begc:endc));     this%begcb_col   (:) = nan
    allocate(this%endcb_col   (begc:endc));     this%endcb_col   (:) = nan
    allocate(this%errcb_col   (begc:endc));     this%errcb_col   (:) = nan

    allocate(this%totpftc_beg_col(begc:endc));  this%totpftc_beg_col (:) = nan
    allocate(this%cwdc_beg_col   (begc:endc));  this%cwdc_beg_col    (:) = nan
    allocate(this%totlitc_beg_col(begc:endc));  this%totlitc_beg_col (:) = nan
    allocate(this%totsomc_beg_col(begc:endc));  this%totsomc_beg_col (:) = nan
    
    allocate(this%totpftc_end_col(begc:endc));  this%totpftc_end_col (:) = nan
    allocate(this%cwdc_end_col   (begc:endc));  this%cwdc_end_col    (:) = nan
    allocate(this%totlitc_end_col(begc:endc));  this%totlitc_end_col (:) = nan
    allocate(this%totsomc_end_col(begc:endc));  this%totsomc_end_col (:) = nan

  end subroutine initallocate_col_ns

!!!!!!  subroutine initcold_col_ns (this, bounds, ratio, c12_carbonstate_vars)
!!!!!!  end subroutine initcold_col_ns

!!!!!!  subroutine inithistory_col_ns (this, bounds, ratio)
!!!!!!  end subroutine inithistory_col_ns

!!!!!!  subroutine restart_col_ns (this, bounds, ncid, flag, c12_carbonstate_vars, cnstate_var)
!!!!!!  end subroutine restart_col_ns


  !-----------------------------------------------------------------------
!   subroutine SetValues ( this, &
!       num_patch, filter_patch, value_patch, &
!       num_column, filter_column, value_column)

   subroutine setvalues_col_ns ( this, &
       num_column, filter_column, value_column)

    !
    ! !DESCRIPTION:
    ! Set carbon state variables at soil column level
    !
    ! !ARGUMENTS:
    class (soilcol_carbonstate_type) :: this
!    integer , intent(in) :: num_patch
!    integer , intent(in) :: filter_patch(:)
!    real(r8), intent(in) :: value_patch
    integer , intent(in) :: num_column
    integer , intent(in) :: filter_column(:)
    real(r8), intent(in) :: value_column
    !
    ! !LOCAL VARIABLES:
    integer :: fi,i,j,k,l     ! loop index
    !------------------------------------------------------------------------

    do fi = 1,num_column
       i = filter_column(fi)
       this%cwdc_col(i)       = value_column
       this%ctrunc_col(i)     = value_column
       this%totlitc_col(i)    = value_column
       this%totsomc_col(i)    = value_column
       this%totecosysc_col(i) = value_column
       this%totcolc_col(i)    = value_column
       this%rootc_col(i)      = value_column
       this%totvegc_col(i)    = value_column
       this%leafc_col(i)      = value_column
       this%deadstemc_col(i)  = value_column
       this%fuelc_col(i)      = value_column
       this%fuelc_crop_col(i) = value_column
       this%totlitc_1m_col(i) = value_column
       this%totsomc_1m_col(i) = value_column
    end do

    do j = 1,nlevdecomp_full
       do fi = 1,num_column
          i = filter_column(fi)
          this%ctrunc_vr_col(i,j) = value_column
       end do
    end do

    do k = 1, ndecomp_pools
       do fi = 1,num_column
          i = filter_column(fi)
          this%decomp_cpools_col(i,k) = value_column
          this%decomp_cpools_1m_col(i,k) = value_column
       end do
    end do

    do j = 1,nlevdecomp_full
       do k = 1, ndecomp_pools
          do fi = 1,num_column
             i = filter_column(fi)
             this%decomp_cpools_vr_col(i,j,k) = value_column
          end do
       end do
    end do

  end subroutine setvalues_col_ns

  !-----------------------------------------------------------------------
 !  subroutine summary_col_ns(this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
  subroutine summary_col_ns(this, bounds, num_soilc, filter_soilc)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, perform column-level carbon summary calculations
    !
    ! !USES:
    use clm_varctl       , only: iulog
    use clm_time_manager , only: get_step_size
    use clm_varcon       , only: secspday
    use clm_varpar       , only: nlevdecomp, ndecomp_pools 
    !
    ! !ARGUMENTS:
    class(carbonstate_type) :: this
    type(bounds_type)      , intent(in)    :: bounds          
    integer                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
!    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
!    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    !
    ! !LOCAL VARIABLES:
    real(r8) :: nfixlags, dtime ! temp variables for making lagged npp
    integer  :: c,p,j,k,l       ! indices
    integer  :: fp,fc           ! lake filter indices
    real(r8) :: maxdepth        ! depth to integrate soil variables
    !-----------------------------------------------------------------------

    ! column level summary

      ! vertically integrate each of the decomposing C pools
      do l = 1, ndecomp_pools
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%decomp_cpools_col(c,l) = 0._r8
       end do
      end do
      do l = 1, ndecomp_pools
       do j = 1, nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%decomp_cpools_col(c,l) = &
                  this%decomp_cpools_col(c,l) + &
                  this%decomp_cpools_vr_col(c,j,l) * dzsoi_decomp(j)
          end do
       end do
      end do

      if ( nlevdecomp > 1) then

       ! vertically integrate each of the decomposing C pools to 1 meter
       maxdepth = 1._r8
       do l = 1, ndecomp_pools
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%decomp_cpools_1m_col(c,l) = 0._r8
          end do
       end do
       do l = 1, ndecomp_pools
          do j = 1, nlevdecomp
             if ( zisoi(j) <= maxdepth ) then
                do fc = 1,num_soilc
                   c = filter_soilc(fc)
                   this%decomp_cpools_1m_col(c,l) = &
                        this%decomp_cpools_1m_col(c,l) + &
                        this%decomp_cpools_vr_col(c,j,l) * dzsoi_decomp(j)
                end do
             elseif ( zisoi(j-1) < maxdepth ) then
                do fc = 1,num_soilc
                   c = filter_soilc(fc)
                   this%decomp_cpools_1m_col(c,l) = &
                        this%decomp_cpools_1m_col(c,l) + &
                        this%decomp_cpools_vr_col(c,j,l) * (maxdepth - zisoi(j-1))
                end do
             endif
          end do
       end do

       ! total litter carbon in the top meter (TOTLITC_1m)
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%totlitc_1m_col(c)         = 0._r8
       end do
       do l = 1, ndecomp_pools
          if ( decomp_cascade_con%is_litter(l) ) then
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                this%totlitc_1m_col(c) = &
                     this%totlitc_1m_col(c) + &
                     this%decomp_cpools_1m_col(c,l)
             end do
          endif
       end do

       ! total soil organic matter carbon in the top meter (TOTSOMC_1m)
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%totsomc_1m_col(c) = 0._r8
       end do
       do l = 1, ndecomp_pools
          if ( decomp_cascade_con%is_soil(l) ) then
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                this%totsomc_1m_col(c) = &
                     this%totsomc_1m_col(c) + &
                     this%decomp_cpools_1m_col(c,l)
             end do
          end if
       end do

      endif
    
      ! total litter carbon (TOTLITC)
      do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%totlitc_col(c) = 0._r8
      end do
      do l = 1, ndecomp_pools
       if ( decomp_cascade_con%is_litter(l) ) then
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%totlitc_col(c) = &
                  this%totlitc_col(c) + &
                  this%decomp_cpools_col(c,l)
          end do
       endif
      end do

      ! total soil organic matter carbon (TOTSOMC)
      do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%totsomc_col(c) = 0._r8
      end do
      do l = 1, ndecomp_pools
       if ( decomp_cascade_con%is_soil(l) ) then
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%totsomc_col(c) = &
                  this%totsomc_col(c) + &
                  this%decomp_cpools_col(c,l)
          end do
       end if
      end do

      ! coarse woody debris carbon
      do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%cwdc_col(c) = 0._r8
      end do
      do l = 1, ndecomp_pools
       if ( decomp_cascade_con%is_cwd(l) ) then
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%cwdc_col(c) = &
                  this%cwdc_col(c) + &
                  this%decomp_cpools_col(c,l)
          end do
       end if
      end do

    ! truncation carbon
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%ctrunc_col(c) = 0._r8
    end do
    do j = 1, nlevdecomp
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%ctrunc_col(c) = &
               this%ctrunc_col(c) + &
               this%ctrunc_vr_col(c,j) * dzsoi_decomp(j)
       end do
    end do

    do fc = 1,num_soilc
       c = filter_soilc(fc)

       ! total product carbon
       this%totprodc_col(c) = &
            this%prod10c_col(c)  + &
            this%prod100c_col(c) + &
            this%prod1c_col(c) 

       ! total ecosystem carbon, including veg but excluding cpool (TOTECOSYSC)
       this%totecosysc_col(c) = &
            this%cwdc_col(c)     + &
            this%totlitc_col(c)  + &
            this%totsomc_col(c)  + &
            this%totprodc_col(c) + &
            this%totvegc_col(c)

       ! total column carbon, including veg and cpool (TOTCOLC)
       ! adding col_ctrunc, seedc
       this%totcolc_col(c) = &
            this%totpftc_col(c)  + &
            this%cwdc_col(c)     + &
            this%totlitc_col(c)  + &
            this%totsomc_col(c)  + &
            this%totprodc_col(c) + &
            this%seedc_col(c)    + &
            this%ctrunc_col(c)
            
       this%totabgc_col(c) = &
            this%totpftc_col(c)  + &
            this%totprodc_col(c) + &
            this%seedc_col(c)    + &
            this%ctrunc_col(c)    
    end do

  end subroutine summary_col_ns

end module ColumnType


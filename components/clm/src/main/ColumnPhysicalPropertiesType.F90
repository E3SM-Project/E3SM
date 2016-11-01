module ColumnPhysicalPropertiesType

  ! -------------------------------------------------------- 
  ! ALM sub-grid hierarchy:
  ! Define column physical properties data type 
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


! moved from ChemStateType
  !------------------------------------------------------------------------------
  ! !USES:
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use decompMod       , only : bounds_type
  use abortutils      , only : endrun
  use clm_varpar      , only : nlevsoi

!  moved from CNDecompCascadeConType.F90

  use decompMod      , only : bounds_type
  use clm_varpar     , only : ndecomp_cascade_transitions, ndecomp_pools


  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Lake data types and associated procesures
  !
  ! !USES:
!  use shr_kind_mod , only : r8 => shr_kind_r8
!  use shr_log_mod  , only : errMsg => shr_log_errMsg
  use clm_varcon   , only : spval, grlnd
!  use decompMod    , only : bounds_type
  use spmdMod      , only : masterproc
!  use abortUtils   , only : endrun
!****** NEED to change to NEW LANDUNITTYPE.
!  use LandunitType , only : lun                
!******  use ColumnType   , only : col  



  !------------------------------------------------------------------------------
  ! !USES:    SoilStateType.F90
  use spmdMod         , only : mpicom, MPI_INTEGER, masterproc
  use ncdio_pio       , only : file_desc_t, ncd_defvar, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
  use ncdio_pio       , only : ncd_pio_openfile, ncd_inqfdims, ncd_pio_closefile, ncd_inqdid, ncd_inqdlen
  use clm_varpar      , only : more_vertlayers, numpft, numrad 
  use clm_varpar      , only : nlevsoi, nlevgrnd, nlevlak, nlevsoifl, nlayer, nlayert, nlevurb, nlevsno
  use landunit_varcon , only : istice, istdlak, istwet, istsoil, istcrop, istice_mec
  use column_varcon   , only : icol_roof, icol_sunwall, icol_shadewall, icol_road_perv, icol_road_imperv 
  use clm_varcon      , only : zsoi, dzsoi, zisoi, spval
  use clm_varcon      , only : secspday, pc, mu, denh2o, denice, grlnd
  use clm_varctl      , only : use_cn, use_lch4,use_dynroot
  use clm_varctl      , only : iulog, fsurdat, hist_wrtch4diag
  use ch4varcon       , only : allowlakeprod
!****  use LandunitType    , only : lun                
!  use ColumnType      , only : col                
!****  use PatchType       , only : pft 

!   From SoilHydrologyType.F90
  use spmdMod               , only : masterproc, mpicom
  use abortutils            , only : endrun
  use clm_varpar            , only : nlevgrnd, nlayer, nlayert, nlevsoi 
  use clm_varpar            , only : more_vertlayers, nlevsoifl, toplev_equalspace 
  use clm_varcon            , only : zsoi, dzsoi, zisoi, spval
  use clm_varctl            , only : iulog 
  use CNSharedParamsMod     , only : CNParamsShareInst
 !!!!!!!    Need to use New LandUnitType here
 ! use LandunitType          , only : lun                
 ! use ColumnType            , only : col   

  implicit none
  save
  private

  ! !PRIVATE MEMBER FUNCTIONS:   From SoilHydrologyType.F90
  private :: initSoilParVIC    ! Convert default CLM soil properties to VIC parameters
  private :: initCLMVICMap     ! Initialize map from VIC to CLM layers
  private :: linear_interp     ! function for linear interperation 


  ! sub-grid geospatial and physical properties defined at the soil_column level
  ! migrate variables list from ColumnType.F90

  type, public :: column_phyisical_properties
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
!  moved from chemStateType
     real(r8), pointer :: soil_pH              (:,:) ! soil pH (-nlevsno+1:nlevgrnd)                   

!  moved from CNDecompCascadeConType.F90
     !-- properties of each pathway along decomposition cascade 
     character(len=8)  , pointer :: cascade_step_name(:)               ! name of transition
     integer           , pointer :: cascade_donor_pool(:)              ! which pool is C taken from for a given decomposition step
     integer           , pointer :: cascade_receiver_pool(:)           ! which pool is C added to for a given decomposition step

     !-- properties of each decomposing pool
     logical           , pointer  :: floating_cn_ratio_decomp_pools(:) ! TRUE => pool has fixed C:N ratio
     logical           , pointer  :: floating_cp_ratio_decomp_pools(:) ! TRUE => pool has fixed C:N ratio
     character(len=8)  , pointer  :: decomp_pool_name_restart(:)       ! name of pool for restart files
     character(len=8)  , pointer  :: decomp_pool_name_history(:)       ! name of pool for history files
     character(len=20) , pointer  :: decomp_pool_name_long(:)          ! name of pool for netcdf long names
     character(len=8)  , pointer  :: decomp_pool_name_short(:)         ! name of pool for netcdf short names
     logical           , pointer  :: is_litter(:)                      ! TRUE => pool is a litter pool
     logical           , pointer  :: is_soil(:)                        ! TRUE => pool is a soil pool
     logical           , pointer  :: is_cwd(:)                         ! TRUE => pool is a cwd pool
     real(r8)          , pointer  :: initial_cn_ratio(:)               ! c:n ratio for initialization of pools
     real(r8)          , pointer  :: initial_cp_ratio(:)               ! c:n ratio for initialization of pools
     real(r8)          , pointer  :: initial_stock(:)                  ! initial concentration for seeding at spinup
     logical           , pointer  :: is_metabolic(:)                   ! TRUE => pool is metabolic material
     logical           , pointer  :: is_cellulose(:)                   ! TRUE => pool is cellulose
     logical           , pointer  :: is_lignin(:)                      ! TRUE => pool is lignin
     real(r8)          , pointer  :: spinup_factor(:)                  ! factor by which to scale AD and relevant processes by

!   moved from LakeStateType.F90
     ! Time constant variables
     real(r8), pointer :: lakefetch_col     (:)   ! col lake fetch from surface data (m)                    
     real(r8), pointer :: etal_col          (:)   ! col lake extinction coefficient from surface data (1/m) 

     ! Time varying variables
     real(r8), pointer :: lake_raw_col      (:)   ! col aerodynamic resistance for moisture (s/m)
     real(r8), pointer :: ks_col            (:)   ! col coefficient for calculation of decay of eddy diffusivity with depth
     real(r8), pointer :: ws_col            (:)   ! col surface friction velocity (m/s)
     real(r8), pointer :: ust_lake_col      (:)   ! col friction velocity (m/s)          
     real(r8), pointer :: betaprime_col     (:)   ! col effective beta: sabg_lyr(p,jtop) for snow layers, beta otherwise
     real(r8), pointer :: savedtke1_col     (:)   ! col top level eddy conductivity from previous timestep (W/mK)
     real(r8), pointer :: lake_icefrac_col  (:,:) ! col mass fraction of lake layer that is frozen
     real(r8), pointer :: lake_icethick_col (:)   ! col ice thickness (m) (integrated if lakepuddling)
     real(r8), pointer :: lakeresist_col    (:)   ! col [s/m] (Needed for calc. of grnd_ch4_cond)
!*******     real(r8), pointer :: ram1_lake_patch   (:)   ! patch aerodynamical resistance (s/m)


!   moved from SoilStateType.F90
     ! sand/ clay/ organic matter   ! deleted patch level variables
     real(r8), pointer :: mss_frc_cly_vld_col  (:)   ! col mass fraction clay limited to 0.20
     real(r8), pointer :: cellorg_col          (:,:) ! col organic matter for gridcell containing column (1:nlevsoi)
     real(r8), pointer :: cellsand_col         (:,:) ! sand value for gridcell containing column (1:nlevsoi)
     real(r8), pointer :: cellclay_col         (:,:) ! clay value for gridcell containing column (1:nlevsoi)
     real(r8), pointer :: bd_col               (:,:) ! col bulk density of dry soil material [kg/m^3] (CN)

     ! hydraulic properties
     real(r8), pointer :: hksat_col            (:,:) ! col hydraulic conductivity at saturation (mm H2O /s) 
     real(r8), pointer :: hksat_min_col        (:,:) ! col mineral hydraulic conductivity at saturation (hksat) (mm/s)
     real(r8), pointer :: hk_l_col             (:,:) ! col hydraulic conductivity (mm/s)
     real(r8), pointer :: smp_l_col            (:,:) ! col soil matric potential (mm)
     real(r8), pointer :: smpmin_col           (:)   ! col restriction for min of soil potential (mm) 
     real(r8), pointer :: bsw_col              (:,:) ! col Clapp and Hornberger "b" (nlevgrnd)  
     real(r8), pointer :: watsat_col           (:,:) ! col volumetric soil water at saturation (porosity) 
     real(r8), pointer :: watdry_col           (:,:) ! col btran parameter for btran = 0
     real(r8), pointer :: watopt_col           (:,:) ! col btran parameter for btran = 1
     real(r8), pointer :: watfc_col            (:,:) ! col volumetric soil water at field capacity (nlevsoi)
     real(r8), pointer :: sucsat_col           (:,:) ! col minimum soil suction (mm) (nlevgrnd) 
     real(r8), pointer :: soilbeta_col         (:)   ! col factor that reduces ground evaporation L&P1992(-)
     real(r8), pointer :: soilalpha_col        (:)   ! col factor that reduces ground saturated specific humidity (-)
     real(r8), pointer :: soilalpha_u_col      (:)   ! col urban factor that reduces ground saturated specific humidity (-) 
     real(r8), pointer :: soilpsi_col          (:,:) ! col soil water potential in each soil layer (MPa) (CN)
     real(r8), pointer :: wtfact_col           (:)   ! col maximum saturated fraction for a gridcell
     real(r8), pointer :: porosity_col         (:,:) ! col soil porisity (1-bulk_density/soil_density) (VIC)
     real(r8), pointer :: eff_porosity_col     (:,:) ! col effective porosity = porosity - vol_ice (nlevgrnd) 
     real(r8), pointer :: gwc_thr_col          (:)   ! col threshold soil moisture based on clay content

     ! thermal conductivity / heat capacity
     real(r8), pointer :: thk_col              (:,:) ! col thermal conductivity of each layer [W/m-K] 
     real(r8), pointer :: tkmg_col             (:,:) ! col thermal conductivity, soil minerals  [W/m-K] (new) (nlevgrnd) 
     real(r8), pointer :: tkdry_col            (:,:) ! col thermal conductivity, dry soil (W/m/Kelvin) (nlevgrnd) 
     real(r8), pointer :: tksatu_col           (:,:) ! col thermal conductivity, saturated soil [W/m-K] (new) (nlevgrnd) 
     real(r8), pointer :: csol_col             (:,:) ! col heat capacity, soil solids (J/m**3/Kelvin) (nlevgrnd) 

     ! roots   deleted patch level variables
     real(r8), pointer :: rootr_col            (:,:) ! col effective fraction of roots in each soil layer (nlevgrnd)  
     real(r8), pointer :: rootfr_col           (:,:) ! col fraction of roots in each soil layer (nlevgrnd) 
     real(r8), pointer :: rootr_road_perv_col  (:,:) ! col effective fraction of roots in each soil layer of urban pervious road
     real(r8), pointer :: rootfr_road_perv_col (:,:) ! col effective fraction of roots in each soil layer of urban pervious road

!    moved from SoilorderConType.F90
     real(r8), allocatable :: smax(:)
     real(r8), allocatable :: ks_sorption(:)
     real(r8), allocatable :: r_weather(:)
     real(r8), allocatable :: r_adsorp(:)
     real(r8), allocatable :: r_desorp(:)
     real(r8), allocatable :: r_occlude(:)
     real(r8), allocatable :: k_s1_biochem(:)
     real(r8), allocatable :: k_s2_biochem(:)
     real(r8), allocatable :: k_s3_biochem(:)
     real(r8), allocatable :: k_s4_biochem(:)


!    moved from SoilHydrologyType.F90

     integer :: h2osfcflag              ! true => surface water is active (namelist)       
     integer :: origflag                ! used to control soil hydrology properties (namelist)

     ! NON-VIC
     real(r8), pointer :: frost_table_col   (:)     ! col frost table depth                    
     real(r8), pointer :: zwt_col           (:)     ! col water table depth
     real(r8), pointer :: zwts_col          (:)     ! col water table depth, the shallower of the two water depths     
     real(r8), pointer :: zwt_perched_col   (:)     ! col perched water table depth
     real(r8), pointer :: wa_col            (:)     ! col water in the unconfined aquifer (mm)
     real(r8), pointer :: qcharge_col       (:)     ! col aquifer recharge rate (mm/s) 
     real(r8), pointer :: fracice_col       (:,:)   ! col fractional impermeability (-)
     real(r8), pointer :: icefrac_col       (:,:)   ! col fraction of ice       
     real(r8), pointer :: fcov_col          (:)     ! col fractional impermeable area
     real(r8), pointer :: fsat_col          (:)     ! col fractional area with water table at surface
     real(r8), pointer :: h2osfc_thresh_col (:)     ! col level at which h2osfc "percolates"   (time constant)

     ! VIC 
     real(r8), pointer :: hkdepth_col       (:)     ! col VIC decay factor (m) (time constant)                    
     real(r8), pointer :: b_infil_col       (:)     ! col VIC b infiltration parameter (time constant)                    
     real(r8), pointer :: ds_col            (:)     ! col VIC fracton of Dsmax where non-linear baseflow begins (time constant)                    
     real(r8), pointer :: dsmax_col         (:)     ! col VIC max. velocity of baseflow (mm/day) (time constant)
     real(r8), pointer :: Wsvic_col         (:)     ! col VIC fraction of maximum soil moisutre where non-liear base flow occurs (time constant)
     real(r8), pointer :: porosity_col      (:,:)   ! col VIC porosity (1-bulk_density/soil_density)
     real(r8), pointer :: vic_clm_fract_col (:,:,:) ! col VIC fraction of VIC layers in CLM layers 
     real(r8), pointer :: depth_col         (:,:)   ! col VIC layer depth of upper layer  
     real(r8), pointer :: c_param_col       (:)     ! col VIC baseflow exponent (Qb) 
     real(r8), pointer :: expt_col          (:,:)   ! col VIC pore-size distribution related paramter(Q12) 
     real(r8), pointer :: ksat_col          (:,:)   ! col VIC Saturated hydrologic conductivity 
     real(r8), pointer :: phi_s_col         (:,:)   ! col VIC soil moisture dissusion parameter 
     real(r8), pointer :: moist_col         (:,:)   ! col VIC soil moisture (kg/m2) for VIC soil layers 
     real(r8), pointer :: moist_vol_col     (:,:)   ! col VIC volumetric soil moisture for VIC soil layers 
     real(r8), pointer :: max_moist_col     (:,:)   ! col VIC max layer moist + ice (mm) 
     real(r8), pointer :: max_infil_col     (:)     ! col VIC maximum infiltration rate calculated in VIC
     real(r8), pointer :: i_0_col           (:)     ! col VIC average saturation in top soil layers 
     real(r8), pointer :: ice_col           (:,:)   ! col VIC soil ice (kg/m2) for VIC soil layers

   contains

     procedure, public :: Init => init_col_pp
     procedure, public :: Clean => clean_col_pp

!    moved from LakeStateType.F90
!     procedure, public  :: Init        
     procedure, public  :: Restart  => restart_col_pp
     procedure, private :: InitAllocate => initallocation_col_pp
     procedure, private :: InitHistory  => inithistory_col_pp
     procedure, private :: InitCold  => initcold_col_pp

!!!!!****   moved from SoilHydrologyType.F90
     procedure, private :: ReadSoilHydrologyNL

  end type column_physical_properties

  ! declare the public instances of soilcolumn types
  type(comlumn_physical_properties) , public, target :: col_pp

  contains 
!!!!***********DW***********change to (this, bounds) based on ChemStateType.F90
!  subroutine init_col_pp(this, begc, endc)
!    class(column_phyiscal_properties) :: this
!    integer, intent(in) :: begc   ! beginning soil column index
!    integer, intent(in) :: endc   ! ending soil column index


  !------------------------------------------------------------------------
  subroutine init_col_pp(this, bounds)

    class(column_physical_properties) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate ( bounds )
    call this%InitHistory ( bounds )
    call this%InitCold ( bounds ) 

  end subroutine init_col_pps

!  
   subroutine initallocation_col_pp(this, bounds)
!   moved from ChemStateType.f90

!    use clm_varpar            , only : nlevsoi

    ! !USES:  moved from SoilorderConType.F90
    use clm_varpar, only : nsoilorder
    use soilorder_varcon, only : smax
    use soilorder_varcon, only : ks_sorption
    use soilorder_varcon, only : r_weather,r_adsorp,r_desorp,r_occlude
    use soilorder_varcon, only : k_s1_biochem,k_s2_biochem,k_s3_biochem,k_s4_biochem

    class(column_physical_properties)         :: this
    type(bounds_type), intent(in) :: bounds 

    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg
    integer :: lbj,  ubj

        ! !LOCAL VARIABLES:  SoilorderConType.F90
    integer :: m, ib

    begc = bounds%begc;
    endc = bounds%endc
    lbj  = 1;
    ubj  = nlevsoi

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
!   moved from ChemStateType.F90
    allocate(this%soil_pH(begc:endc, lbj:ubj))

!  moved from CNDecompCascadeConType.F90
    ! !DESCRIPTION:
    ! Initialize decomposition cascade state
    !------------------------------------------------------------------------

    !-- properties of each pathway along decomposition cascade 
    allocate(this%cascade_step_name(1:ndecomp_cascade_transitions))
    allocate(this%cascade_donor_pool(1:ndecomp_cascade_transitions))
    allocate(this%cascade_receiver_pool(1:ndecomp_cascade_transitions))

    !-- properties of each decomposing pool
    allocate(this%floating_cn_ratio_decomp_pools(0:ndecomp_pools))
    allocate(this%floating_cp_ratio_decomp_pools(0:ndecomp_pools))
    allocate(this%decomp_pool_name_restart(0:ndecomp_pools))
    allocate(this%decomp_pool_name_history(0:ndecomp_pools))
    allocate(this%decomp_pool_name_long(0:ndecomp_pools))
    allocate(this%decomp_pool_name_short(0:ndecomp_pools))
    allocate(this%is_litter(0:ndecomp_pools))
    allocate(this%is_soil(0:ndecomp_pools))
    allocate(this%is_cwd(0:ndecomp_pools))
    allocate(this%initial_cn_ratio(0:ndecomp_pools))
    allocate(this%initial_cp_ratio(0:ndecomp_pools))
    allocate(this%initial_stock(0:ndecomp_pools))
    allocate(this%is_metabolic(0:ndecomp_pools))
    allocate(this%is_cellulose(0:ndecomp_pools))
    allocate(this%is_lignin(0:ndecomp_pools))
    allocate(this%spinup_factor(0:ndecomp_pools))

    !-- properties of each pathway along decomposition cascade 
    this%cascade_step_name(1:ndecomp_cascade_transitions) = ''
    this%cascade_donor_pool(1:ndecomp_cascade_transitions) = 0
    this%cascade_receiver_pool(1:ndecomp_cascade_transitions) = 0

    !-- first initialization of properties of each decomposing pool
    this%floating_cn_ratio_decomp_pools(0:ndecomp_pools) = .false.
    this%floating_cp_ratio_decomp_pools(0:ndecomp_pools) = .false.
    this%decomp_pool_name_history(0:ndecomp_pools)       = ''
    this%decomp_pool_name_restart(0:ndecomp_pools)       = ''
    this%decomp_pool_name_long(0:ndecomp_pools)          = ''
    this%decomp_pool_name_short(0:ndecomp_pools)         = ''
    this%is_litter(0:ndecomp_pools)                      = .false.
    this%is_soil(0:ndecomp_pools)                        = .false.
    this%is_cwd(0:ndecomp_pools)                         = .false.
    this%initial_cn_ratio(0:ndecomp_pools)               = nan
    this%initial_cp_ratio(0:ndecomp_pools)               = nan
    this%initial_stock(0:ndecomp_pools)                  = nan
    this%is_metabolic(0:ndecomp_pools)                   = .false.
    this%is_cellulose(0:ndecomp_pools)                   = .false.
    this%is_lignin(0:ndecomp_pools)                      = .false.
    this%spinup_factor(0:ndecomp_pools)                  = nan
  

 ! moved from LakeStateType.F90

    allocate(this%etal_col           (begc:endc))           ; this%etal_col           (:)   = nan
    allocate(this%lakefetch_col      (begc:endc))           ; this%lakefetch_col      (:)   = nan
    allocate(this%lakeresist_col     (begc:endc))           ; this%lakeresist_col     (:)   = nan
    allocate(this%savedtke1_col      (begc:endc))           ; this%savedtke1_col      (:)   = spval  
    allocate(this%lake_icefrac_col   (begc:endc,1:nlevlak)) ; this%lake_icefrac_col   (:,:) = nan
    allocate(this%lake_icethick_col  (begc:endc))           ; this%lake_icethick_col  (:)   = nan
    allocate(this%ust_lake_col       (begc:endc))           ; this%ust_lake_col       (:)   = spval   
    allocate(this%ram1_lake_patch      (begp:endp))           ; this%ram1_lake_patch      (:)   = nan
    allocate(this%lake_raw_col       (begc:endc))           ; this%lake_raw_col       (:)   = nan
    allocate(this%ks_col             (begc:endc))           ; this%ks_col             (:)   = nan
    allocate(this%ws_col             (begc:endc))           ; this%ws_col             (:)   = nan
    allocate(this%betaprime_col      (begc:endc))           ; this%betaprime_col      (:)   = nan


!  moved from SoilStateType.F90

    allocate(this%mss_frc_cly_vld_col  (begc:endc))                     ; this%mss_frc_cly_vld_col  (:)   = nan
 !   allocate(this%sandfrac_patch       (begp:endp))                     ; this%sandfrac_patch       (:)   = nan
 !   allocate(this%clayfrac_patch       (begp:endp))                     ; this%clayfrac_patch       (:)   = nan
    allocate(this%cellorg_col          (begc:endc,nlevsoi))             ; this%cellorg_col          (:,:) = nan 
    allocate(this%cellsand_col         (begc:endc,nlevsoi))             ; this%cellsand_col         (:,:) = nan 
    allocate(this%cellclay_col         (begc:endc,nlevsoi))             ; this%cellclay_col         (:,:) = nan 
    allocate(this%bd_col               (begc:endc,nlevgrnd))            ; this%bd_col               (:,:) = nan

    allocate(this%hksat_col            (begc:endc,nlevgrnd))            ; this%hksat_col            (:,:) = spval
    allocate(this%hksat_min_col        (begc:endc,nlevgrnd))            ; this%hksat_min_col        (:,:) = spval
    allocate(this%hk_l_col             (begc:endc,nlevgrnd))            ; this%hk_l_col             (:,:) = nan   
    allocate(this%smp_l_col            (begc:endc,nlevgrnd))            ; this%smp_l_col            (:,:) = nan   
    allocate(this%smpmin_col           (begc:endc))                     ; this%smpmin_col           (:)   = nan

    allocate(this%bsw_col              (begc:endc,nlevgrnd))            ; this%bsw_col              (:,:) = nan
    allocate(this%watsat_col           (begc:endc,nlevgrnd))            ; this%watsat_col           (:,:) = nan
    allocate(this%watdry_col           (begc:endc,nlevgrnd))            ; this%watdry_col           (:,:) = spval
    allocate(this%watopt_col           (begc:endc,nlevgrnd))            ; this%watopt_col           (:,:) = spval
    allocate(this%watfc_col            (begc:endc,nlevgrnd))            ; this%watfc_col            (:,:) = nan
    allocate(this%sucsat_col           (begc:endc,nlevgrnd))            ; this%sucsat_col           (:,:) = spval
    allocate(this%soilbeta_col         (begc:endc))                     ; this%soilbeta_col         (:)   = nan   
    allocate(this%soilalpha_col        (begc:endc))                     ; this%soilalpha_col        (:)   = nan
    allocate(this%soilalpha_u_col      (begc:endc))                     ; this%soilalpha_u_col      (:)   = nan
    allocate(this%soilpsi_col          (begc:endc,nlevgrnd))            ; this%soilpsi_col          (:,:) = nan
    allocate(this%wtfact_col           (begc:endc))                     ; this%wtfact_col           (:)   = nan
    allocate(this%porosity_col         (begc:endc,nlayer))              ; this%porosity_col         (:,:) = spval
    allocate(this%eff_porosity_col     (begc:endc,nlevgrnd))            ; this%eff_porosity_col     (:,:) = spval
    allocate(this%gwc_thr_col          (begc:endc))                     ; this%gwc_thr_col          (:)   = nan

    allocate(this%thk_col              (begc:endc,-nlevsno+1:nlevgrnd)) ; this%thk_col              (:,:) = nan
    allocate(this%tkmg_col             (begc:endc,nlevgrnd))            ; this%tkmg_col             (:,:) = nan
    allocate(this%tkdry_col            (begc:endc,nlevgrnd))            ; this%tkdry_col            (:,:) = nan
    allocate(this%tksatu_col           (begc:endc,nlevgrnd))            ; this%tksatu_col           (:,:) = nan
    allocate(this%csol_col             (begc:endc,nlevgrnd))            ; this%csol_col             (:,:) = nan

  !  allocate(this%rootr_patch          (begp:endp,1:nlevgrnd))          ; this%rootr_patch          (:,:) = nan
    allocate(this%rootr_col            (begc:endc,nlevgrnd))            ; this%rootr_col            (:,:) = nan
    allocate(this%rootr_road_perv_col  (begc:endc,1:nlevgrnd))          ; this%rootr_road_perv_col  (:,:) = nan
  !  allocate(this%rootfr_patch         (begp:endp,1:nlevgrnd))          ; this%rootfr_patch         (:,:) = nan
    allocate(this%rootfr_col           (begc:endc,1:nlevgrnd))          ; this%rootfr_col           (:,:) = nan 
    allocate(this%rootfr_road_perv_col (begc:endc,1:nlevgrnd))          ; this%rootfr_road_perv_col (:,:) = nan
  !  allocate(this%root_depth_patch     (begp:endp))                     ; this%root_depth_patch     (:)   = spval

!   moved from SoilorderConType.F90
    allocate(this%smax           (0:nsoilorder))        ; this%smax(:)        =nan
    allocate(this%ks_sorption    (0:nsoilorder))        ; this%ks_sorption(:) =nan
    allocate(this%r_weather      (0:nsoilorder))        ; this%r_weather(:)   =nan
    allocate(this%r_adsorp       (0:nsoilorder))        ; this%r_adsorp(:)    =nan
    allocate(this%r_desorp       (0:nsoilorder))        ; this%r_desorp(:)    =nan
    allocate(this%r_occlude      (0:nsoilorder))        ; this%r_occlude(:)    =nan

    allocate(this%k_s1_biochem      (0:nsoilorder))        ; this%k_s1_biochem(:)    =nan
    allocate(this%k_s2_biochem      (0:nsoilorder))        ; this%k_s2_biochem(:)    =nan
    allocate(this%k_s3_biochem      (0:nsoilorder))        ; this%k_s3_biochem(:)    =nan
    allocate(this%k_s4_biochem      (0:nsoilorder))        ; this%k_s4_biochem(:)    =nan


    do m = 0,nsoilorder


       this%smax(m)         = smax(m)
       this%ks_sorption(m)         = ks_sorption(m)
       this%r_weather(m)         = r_weather(m)
       this%r_adsorp(m)         = r_adsorp(m)
       this%r_desorp(m)         = r_desorp(m)
       this%r_occlude(m)         = r_occlude(m)
       this%k_s1_biochem(m)         = k_s1_biochem(m)
       this%k_s2_biochem(m)         = k_s2_biochem(m)
       this%k_s3_biochem(m)         = k_s3_biochem(m)
       this%k_s4_biochem(m)         = k_s4_biochem(m)


    end do


  end subroutine initallocation_col_pp


    !-----------------------------------------------------------------------
  !! Moved from LakeStateType.F90

  subroutine inithistory_col_pp(this, bounds)
    !
    ! History fields initialization
    !
    ! !USES:
    use shr_infnan_mod, only: nan => shr_infnan_nan, assignment(=)
    use histFileMod   , only: hist_addfld1d, hist_addfld2d

!   mvoed from SoilStateType.F90
    use histFileMod   , only: hist_addfld1d, hist_addfld2d, no_snow_normal
    !
    ! !ARGUMENTS:
    class(column_physical_properties) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !---------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    this%lake_icefrac_col(begc:endc,:) = spval
    call hist_addfld2d (fname='LAKEICEFRAC',  units='unitless', type2d='levlak', &
         avgflag='A', long_name='lake layer ice mass fraction', &
         ptr_col=this%lake_icefrac_col)
    
    this%lake_icethick_col(begc:endc) = spval ! This will be more useful than LAKEICEFRAC for many users.
    call hist_addfld1d (fname='LAKEICETHICK', units='m', &
         avgflag='A', long_name='thickness of lake ice (including physical expansion on freezing)', &
         ptr_col=this%lake_icethick_col, set_nolake=spval)

    this%savedtke1_col(begc:endc) = spval
    call hist_addfld1d (fname='TKE1',  units='W/(mK)', &
         avgflag='A', long_name='top lake level eddy thermal conductivity', &
         ptr_col=this%savedtke1_col)

!**********need to moved to patch level
!    this%ram1_lake_patch(begp:endp) = spval
!   call hist_addfld1d (fname='RAM_LAKE', units='s/m', &
!         avgflag='A', long_name='aerodynamic resistance for momentum (lakes only)', &
!         ptr_patch=this%ram1_lake_patch, set_nolake=spval, default='inactive')

    this%ust_lake_col(begc:endc) = spval
    call hist_addfld1d (fname='UST_LAKE', units='m/s', &
         avgflag='A', long_name='friction velocity (lakes only)', &
         ptr_col=this%ust_lake_col, set_nolake=spval, default='inactive')


!   moved from SoilStateType.F90

      if (use_lch4) then
       if (hist_wrtch4diag) then
          active = "active"
       else
          active = "inactive"
       end if
    else
       active = "inactive"
    end if
    call hist_addfld2d (fname='SMP',  units='mm', type2d='levgrnd',  &
         avgflag='A', long_name='soil matric potential (vegetated landunits only)', &
         ptr_col=this%smp_l_col, set_spec=spval, l2g_scale_type='veg', default=active)

    if (use_cn) then
       this%bsw_col(begc:endc,:) = spval 
       call hist_addfld2d (fname='bsw', units='unitless', type2d='levgrnd', &
            avgflag='A', long_name='clap and hornberger B', &
            ptr_col=this%bsw_col, default='inactive')
    end if

  !  if (use_cn) then
  !     this%rootfr_patch(begp:endp,:) = spval
  1     call hist_addfld2d (fname='ROOTFR', units='proportion', type2d='levgrnd', &
  !          avgflag='A', long_name='fraction of roots in each soil layer', &
  !          ptr_patch=this%rootfr_patch, default='inactive')
    end if

  !  if (use_cn) then
  !     this%rootr_patch(begp:endp,:) = spval
  !     call hist_addfld2d (fname='ROOTR', units='proportion', type2d='levgrnd', &
  !          avgflag='A', long_name='effective fraction of roots in each soil layer', &
  !          ptr_patch=this%rootr_patch, default='inactive')
  !  end if

    if (use_cn) then
       this%rootr_col(begc:endc,:) = spval
       call hist_addfld2d (fname='ROOTR_COLUMN', units='proportion', type2d='levgrnd', &
            avgflag='A', long_name='effective fraction of roots in each soil layer', &
            ptr_col=this%rootr_col, default='inactive')
       
    end if

 !   if (use_dynroot) then
 !      this%root_depth_patch(begp:endp) = spval
 !      call hist_addfld1d (fname='ROOT_DEPTH', units="m", &
 !           avgflag='A', long_name='rooting depth', &
 !           ptr_patch=this%root_depth_patch, default='inactive' )
 !   end if

    if (use_cn) then
       this%soilpsi_col(begc:endc,:) = spval
       call hist_addfld2d (fname='SOILPSI', units='MPa', type2d='levgrnd', &
            avgflag='A', long_name='soil water potential in each soil layer', &
            ptr_col=this%soilpsi_col)
    end if

    this%thk_col(begc:endc,-nlevsno+1:0) = spval
    data2dptr => this%thk_col(:,-nlevsno+1:0)
    call hist_addfld2d (fname='SNO_TK', units='W/m-K', type2d='levsno', &
         avgflag='A', long_name='Thermal conductivity', &
         ptr_col=data2dptr, no_snow_behavior=no_snow_normal, default='inactive')

    this%hk_l_col(begc:endc,:) = spval
    call hist_addfld2d (fname='HK',  units='mm/s', type2d='levgrnd',  &
         avgflag='A', long_name='hydraulic conductivity (vegetated landunits only)', &
         ptr_col=this%hk_l_col, set_spec=spval, l2g_scale_type='veg', default='inactive')

    this%soilalpha_col(begc:endc) = spval
    call hist_addfld1d (fname='SoilAlpha',  units='unitless',  &
         avgflag='A', long_name='factor limiting ground evap', &
         ptr_col=this%soilalpha_col, set_urb=spval)

    this%soilalpha_u_col(begc:endc) = spval
    call hist_addfld1d (fname='SoilAlpha_U',  units='unitless',  &
         avgflag='A', long_name='urban factor limiting ground evap', &
         ptr_col=this%soilalpha_u_col, set_nourb=spval)

    if (use_cn) then
       this%watsat_col(begc:endc,:) = spval 
       call hist_addfld2d (fname='watsat', units='m^3/m^3', type2d='levgrnd', &
            avgflag='A', long_name='water saturated', &
            ptr_col=this%watsat_col, default='inactive')
    end if

    if (use_cn) then
       this%eff_porosity_col(begc:endc,:) = spval
       call hist_addfld2d (fname='EFF_POROSITY', units='proportion', type2d='levgrnd', &
            avgflag='A', long_name='effective porosity = porosity - vol_ice', &
            ptr_col=this%eff_porosity_col, default='inactive')
    end if

    if (use_cn) then
       this%watfc_col(begc:endc,:) = spval 
       call hist_addfld2d (fname='watfc', units='m^3/m^3', type2d='levgrnd', &
            avgflag='A', long_name='water field capacity', &
            ptr_col=this%watfc_col, default='inactive')
    end if


  end subroutine inithistory_col_pp


  !-----------------------------------------------------------------------
  !*******************DW************* need to deal with landunit and gridcell variable

  subroutine initcold_col_pp(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize time constant and time varying module variables
    !
    ! !USES:
    use clm_varctl , only : fsurdat
    use clm_varctl , only : iulog
    use clm_varpar , only : nlevlak
    use clm_varcon , only : tkwat 
    use fileutils  , only : getfil
    use ncdio_pio  , only : file_desc_t, ncd_defvar, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
    use ncdio_pio  , only : ncd_pio_openfile, ncd_inqfdims, ncd_pio_closefile, ncd_inqdid, ncd_inqdlen

    ! moved from SoilStateType.F90
    ! Initialize module surface albedos to reasonable values
    !
    ! !USES:
 !*****   use pftvarcon           , only : noveg, roota_par, rootb_par
    use fileutils           , only : getfil
    use organicFileMod      , only : organicrd 
    use CNSharedParamsMod   , only : CNParamsShareInst
    use FuncPedotransferMod , only : pedotransf, get_ipedof
    use RootBiophysMod      , only : init_vegrootfr

    !
    ! !ARGUMENTS:
    class(column_physical_properties) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer             :: c,g,i,j,l,lev
    logical             :: readvar 
    type(file_desc_t)   :: ncid                   ! netcdf id
    character(len=256)  :: locfn                  ! local filename
    real(r8)            :: depthratio             ! ratio of lake depth to standard deep lake depth
    real(r8) ,pointer   :: lakefetch_in (:)       ! read in - lakefetch 
    real(r8) ,pointer   :: etal_in (:)            ! read in - etal 
    !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  
    !  moved from SoilStateType.F90
                                                        ! !LOCAL VARIABLES:
    integer            :: p, lev, c, l, g, j            ! indices
    real(r8)           :: om_frac                       ! organic matter fraction
    real(r8)           :: om_tkm         = 0.25_r8      ! thermal conductivity of organic soil (Farouki, 1986) [W/m/K]
    real(r8)           :: om_watsat_lake = 0.9_r8       ! porosity of organic soil
    real(r8)           :: om_hksat_lake  = 0.1_r8       ! saturated hydraulic conductivity of organic soil [mm/s]
    real(r8)           :: om_sucsat_lake = 10.3_r8      ! saturated suction for organic matter (Letts, 2000)
    real(r8)           :: om_b_lake      = 2.7_r8       ! Clapp Hornberger paramater for oragnic soil (Letts, 2000) (lake)
    real(r8)           :: om_watsat                     ! porosity of organic soil
    real(r8)           :: om_hksat                      ! saturated hydraulic conductivity of organic soil [mm/s]
    real(r8)           :: om_sucsat                     ! saturated suction for organic matter (mm)(Letts, 2000)
    real(r8)           :: om_csol        = 2.5_r8       ! heat capacity of peat soil *10^6 (J/K m3) (Farouki, 1986)
    real(r8)           :: om_tkd         = 0.05_r8      ! thermal conductivity of dry organic soil (Farouki, 1981)
    real(r8)           :: om_b                          ! Clapp Hornberger paramater for oragnic soil (Letts, 2000)
    real(r8)           :: zsapric        = 0.5_r8       ! depth (m) that organic matter takes on characteristics of sapric peat
    real(r8)           :: csol_bedrock   = 2.0e6_r8     ! vol. heat capacity of granite/sandstone  J/(m3 K)(Shabbir, 2000)
    real(r8)           :: pcalpha        = 0.5_r8       ! percolation threshold
    real(r8)           :: pcbeta         = 0.139_r8     ! percolation exponent
    real(r8)           :: pc_lake        = 0.5_r8       ! percolation threshold
    real(r8)           :: perc_frac                     ! "percolating" fraction of organic soil
    real(r8)           :: perc_norm                     ! normalize to 1 when 100% organic soil
    real(r8)           :: uncon_hksat                   ! series conductivity of mineral/organic soil
    real(r8)           :: uncon_frac                    ! fraction of "unconnected" soil
    real(r8)           :: bd                            ! bulk density of dry soil material [kg/m^3]
    real(r8)           :: tkm                           ! mineral conductivity
    real(r8)           :: xksat                         ! maximum hydraulic conductivity of soil [mm/s]
    real(r8)           :: clay,sand                     ! temporaries
    real(r8)           :: organic_max                   ! organic matter (kg/m3) where soil is assumed to act like peat
    integer            :: dimid                         ! dimension id
    logical            :: readvar 
    type(file_desc_t)  :: ncid                          ! netcdf id
    real(r8) ,pointer  :: zsoifl (:)                    ! Output: [real(r8) (:)]  original soil midpoint 
    real(r8) ,pointer  :: zisoifl (:)                   ! Output: [real(r8) (:)]  original soil interface depth 
    real(r8) ,pointer  :: dzsoifl (:)                   ! Output: [real(r8) (:)]  original soil thickness 
    real(r8) ,pointer  :: gti (:)                       ! read in - fmax 
    real(r8) ,pointer  :: sand3d (:,:)                  ! read in - soil texture: percent sand (needs to be a pointer for use in ncdio)
    real(r8) ,pointer  :: clay3d (:,:)                  ! read in - soil texture: percent clay (needs to be a pointer for use in ncdio)
    real(r8) ,pointer  :: organic3d (:,:)               ! read in - organic matter: kg/m3 (needs to be a pointer for use in ncdio)
    character(len=256) :: locfn                         ! local filename
    integer            :: ipedof  
    integer            :: begc, endc
    integer            :: begg, endg
    !-----------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg




    !-------------------------------------------------
    ! Initialize time constant variables
    !-------------------------------------------------

    call getfil (fsurdat, locfn, 0)
    call ncd_pio_openfile (ncid, locfn, 0)

    ! Read lake eta
    allocate(etal_in(bounds%begg:bounds%endg))
    call ncd_io(ncid=ncid, varname='ETALAKE', flag='read', data=etal_in, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       if (masterproc) then
          write(iulog,*) 'WARNING:: ETALAKE not found on surface data set. All lake columns will have eta', &
               ' set equal to default value as a function of depth.'
       end if
       etal_in(:) = -1._r8
    end if
    do c = bounds%begc, bounds%endc
       g = col%gridcell(c)
       this%etal_col(c) = etal_in(g)
    end do
    deallocate(etal_in)

    ! Read lake fetch
    allocate(lakefetch_in(bounds%begg:bounds%endg))
    call ncd_io(ncid=ncid, varname='LAKEFETCH', flag='read', data=lakefetch_in, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       if (masterproc) then
          write(iulog,*) 'WARNING:: LAKEFETCH not found on surface data set. All lake columns will have fetch', &
               ' set equal to default value as a function of depth.'
       end if
       lakefetch_in(:) = -1._r8
    end if
    do c = bounds%begc, bounds%endc
       g = col%gridcell(c)
       this%lakefetch_col(c) = lakefetch_in(g)
    end do
    deallocate(lakefetch_in)

    call ncd_pio_closefile(ncid)

    !-------------------------------------------------
    ! Initialize time varying variables
    !-------------------------------------------------
         
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%lakpoi(l)) then

          ! Set lake ice fraction and top eddy conductivity from previous timestep
          ! Always initialize with no ice to prevent excessive ice sheets from forming when
          ! starting with old lake model that has unrealistically cold lake conseratures.
          ! Keep lake temperature as is, and the energy deficit below freezing (which is no smaller
          ! than it would have been with prognostic ice, as the temperature would then have been higher
          ! and more heat would have flowed out of the lake) will be converted to ice in the first timestep.
          this%lake_icefrac_col(c,1:nlevlak) = 0._r8

          ! Set lake top eddy conductivity from previous timestep
          this%savedtke1_col(c) = tkwat

          ! Set column friction vlocity 
          this%ust_lake_col(c)  = 0.1_r8
       end if
    end do

!   moved from SoilStateType.F90


    do c = bounds%begc, bounds%endc
       this%smpmin_col(c) = -1.e8_r8
    end do

    ! --------------------------------------------------------------------
    ! Initialize root fraction (computing from surface, d is depth in meter):
    ! --------------------------------------------------------------------

    ! Currently pervious road has same properties as soil
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)

!!!!!!!   ****** lun and col type here

       if (lun%urbpoi(l) .and. col%itype(c) == icol_road_perv) then 
          do lev = 1, nlevgrnd
             this%rootfr_road_perv_col(c,lev) = 0._r8
          enddo
          do lev = 1,nlevsoi
             this%rootfr_road_perv_col(c,lev) = 0.1_r8  ! uniform profile
          end do
       end if
    end do

    do c = bounds%begc,bounds%endc
       this%rootfr_col (c,nlevsoi+1:nlevgrnd) = 0._r8
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          this%rootfr_col (c,nlevsoi+1:nlevgrnd) = 0._r8
       else if (lun%itype(l) == istdlak .and. allowlakeprod) then
          this%rootfr_col (c,:) = spval
       else  ! Inactive CH4 columns
          this%rootfr_col (c,:) = spval
       end if
    end do

   ! Initialize root fraction 
   
!!***   call init_vegrootfr(bounds, nlevsoi, nlevgrnd, &
!!***        this%rootfr_patch(bounds%begp:bounds%endp,1:nlevgrnd))

    ! --------------------------------------------------------------------
    ! dynamic memory allocation
    ! --------------------------------------------------------------------

    allocate(sand3d(begg:endg,nlevsoifl))
    allocate(clay3d(begg:endg,nlevsoifl))

    ! --------------------------------------------------------------------
    ! Read surface dataset
    ! --------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Attempting to read soil color, sand and clay boundary data .....'
    end if

    call getfil (fsurdat, locfn, 0)
    call ncd_pio_openfile (ncid, locfn, 0)

    call ncd_inqdlen(ncid,dimid,nlevsoifl,name='nlevsoi')
    if ( .not. more_vertlayers )then
       if ( nlevsoifl /= nlevsoi )then
          call endrun(msg=' ERROR: Number of soil layers on file does NOT match the number being used'//&
               errMsg(__FILE__, __LINE__))
       end if
    else
       ! read in layers, interpolate to high resolution grid later
    end if

    ! Read in organic matter dataset 

    organic_max = CNParamsShareInst%organic_max

    allocate(organic3d(bounds%begg:bounds%endg,nlevsoifl))
    call organicrd(organic3d)

    ! Read in sand and clay data

    call ncd_io(ncid=ncid, varname='PCT_SAND', flag='read', data=sand3d, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=' ERROR: PCT_SAND NOT on surfdata file'//errMsg(__FILE__, __LINE__)) 
    end if

    call ncd_io(ncid=ncid, varname='PCT_CLAY', flag='read', data=clay3d, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=' ERROR: PCT_CLAY NOT on surfdata file'//errMsg(__FILE__, __LINE__)) 
    end if


!!!!!!  need to move to patch level 
!***    do p = bounds%begp,bounds%endp
!       g = pft%gridcell(p)
!       if ( sand3d(g,1)+clay3d(g,1) == 0.0_r8 )then
!          if ( any( sand3d(g,:)+clay3d(g,:) /= 0.0_r8 ) )then
!             call endrun(msg='found depth points that do NOT sum to zero when surface does'//&
!                  errMsg(__FILE__, __LINE__)) 
!          end if
!          sand3d(g,:) = 1.0_r8
!          clay3d(g,:) = 1.0_r8
!       end if
!       if ( any( sand3d(g,:)+clay3d(g,:) == 0.0_r8 ) )then
!          call endrun(msg='after setting, found points sum to zero'//errMsg(__FILE__, __LINE__)) 
!       end if
!
!       this%sandfrac_patch(p) = sand3d(g,1)/100.0_r8
!       this%clayfrac_patch(p) = clay3d(g,1)/100.0_r8
!    end do

    ! Read fmax

    allocate(gti(bounds%begg:bounds%endg))
    call ncd_io(ncid=ncid, varname='FMAX', flag='read', data=gti, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=' ERROR: FMAX NOT on surfdata file'//errMsg(__FILE__, __LINE__)) 
    end if
    do c = bounds%begc, bounds%endc
       g = col%gridcell(c)
       this%wtfact_col(c) = gti(g)
    end do
    deallocate(gti)

    ! Close file

    call ncd_pio_closefile(ncid)

    ! --------------------------------------------------------------------
    ! get original soil depths to be used in interpolation of sand and clay
    ! --------------------------------------------------------------------

    allocate(zsoifl(1:nlevsoifl), zisoifl(0:nlevsoifl), dzsoifl(1:nlevsoifl))
    do j = 1, nlevsoifl
       zsoifl(j) = 0.025*(exp(0.5_r8*(j-0.5_r8))-1._r8)    !node depths
    enddo

    dzsoifl(1) = 0.5_r8*(zsoifl(1)+zsoifl(2))             !thickness b/n two interfaces
    do j = 2,nlevsoifl-1
       dzsoifl(j)= 0.5_r8*(zsoifl(j+1)-zsoifl(j-1))
    enddo
    dzsoifl(nlevsoifl) = zsoifl(nlevsoifl)-zsoifl(nlevsoifl-1)

    zisoifl(0) = 0._r8
    do j = 1, nlevsoifl-1
       zisoifl(j) = 0.5_r8*(zsoifl(j)+zsoifl(j+1))         !interface depths
    enddo
    zisoifl(nlevsoifl) = zsoifl(nlevsoifl) + 0.5_r8*dzsoifl(nlevsoifl)

    ! --------------------------------------------------------------------
    ! Set soil hydraulic and thermal properties: non-lake
    ! --------------------------------------------------------------------

    !   urban roof, sunwall and shadewall thermal properties used to 
    !   derive thermal conductivity and heat capacity are set to special 
    !   value because thermal conductivity and heat capacity for urban 
    !   roof, sunwall and shadewall are prescribed in SoilThermProp.F90 
    !   in SoilPhysicsMod.F90


    do c = bounds%begc, bounds%endc
       g = col%gridcell(c)
       l = col%landunit(c)

       if (lun%itype(l)==istwet .or. lun%itype(l)==istice .or. lun%itype(l)==istice_mec) then

          do lev = 1,nlevgrnd
             this%bsw_col(c,lev)    = spval
             this%watsat_col(c,lev) = spval
             this%watfc_col(c,lev)  = spval
             this%hksat_col(c,lev)  = spval
             this%sucsat_col(c,lev) = spval
             this%watdry_col(c,lev) = spval 
             this%watopt_col(c,lev) = spval 
             this%bd_col(c,lev)     = spval 
             if (lev <= nlevsoi) then
                this%cellsand_col(c,lev) = spval
                this%cellclay_col(c,lev) = spval
                this%cellorg_col(c,lev)  = spval
             end if
          end do

          do lev = 1,nlevgrnd
             this%tkmg_col(c,lev)   = spval
             this%tksatu_col(c,lev) = spval
             this%tkdry_col(c,lev)  = spval
             if (lun%itype(l)==istwet .and. lev > nlevsoi) then
                this%csol_col(c,lev) = csol_bedrock
             else
                this%csol_col(c,lev)= spval
             endif
          end do

       else if (lun%urbpoi(l) .and. (col%itype(c) /= icol_road_perv) .and. (col%itype(c) /= icol_road_imperv) )then

          ! Urban Roof, sunwall, shadewall properties set to special value
          do lev = 1,nlevgrnd
             this%watsat_col(c,lev) = spval
             this%watfc_col(c,lev)  = spval
             this%bsw_col(c,lev)    = spval
             this%hksat_col(c,lev)  = spval
             this%sucsat_col(c,lev) = spval
             this%watdry_col(c,lev) = spval 
             this%watopt_col(c,lev) = spval 
             this%bd_col(c,lev) = spval 
             if (lev <= nlevsoi) then
                this%cellsand_col(c,lev) = spval
                this%cellclay_col(c,lev) = spval
                this%cellorg_col(c,lev)  = spval
             end if
          end do

          do lev = 1,nlevgrnd
             this%tkmg_col(c,lev)   = spval
             this%tksatu_col(c,lev) = spval
             this%tkdry_col(c,lev)  = spval
             this%csol_col(c,lev)   = spval
          end do

       else

          do lev = 1,nlevgrnd

             if ( more_vertlayers )then ! duplicate clay and sand values from last soil layer

                if (lev .eq. 1) then
                   clay = clay3d(g,1)
                   sand = sand3d(g,1)
                   om_frac = organic3d(g,1)/organic_max 
                else if (lev <= nlevsoi) then
                   do j = 1,nlevsoifl-1
                      if (zisoi(lev) >= zisoifl(j) .AND. zisoi(lev) < zisoifl(j+1)) then
                         clay = clay3d(g,j+1)
                         sand = sand3d(g,j+1)
                         om_frac = organic3d(g,j+1)/organic_max    
                      endif
                   end do
                else
                   clay = clay3d(g,nlevsoifl)
                   sand = sand3d(g,nlevsoifl)
                   om_frac = 0._r8
                endif
             else
                if (lev <= nlevsoi) then ! duplicate clay and sand values from 10th soil layer
                   clay = clay3d(g,lev)
                   sand = sand3d(g,lev)
                   om_frac = (organic3d(g,lev)/organic_max)**2._r8
                else
                   clay = clay3d(g,nlevsoi)
                   sand = sand3d(g,nlevsoi)
                   om_frac = 0._r8
                endif
             end if

             if (lun%itype(l) == istdlak) then

                if (lev <= nlevsoi) then
                   this%cellsand_col(c,lev) = sand
                   this%cellclay_col(c,lev) = clay
                   this%cellorg_col(c,lev)  = om_frac*organic_max
                end if

             else if (lun%itype(l) /= istdlak) then  ! soil columns of both urban and non-urban types

                if (lun%urbpoi(l)) then
                   om_frac = 0._r8 ! No organic matter for urban
                end if

                if (lev <= nlevsoi) then
                   this%cellsand_col(c,lev) = sand
                   this%cellclay_col(c,lev) = clay
                   this%cellorg_col(c,lev)  = om_frac*organic_max
                end if

                ! Note that the following properties are overwritten for urban impervious road 
                ! layers that are not soil in SoilThermProp.F90 within SoilTemperatureMod.F90

                !determine the type of pedotransfer function to be used based on soil order
                !I will use the following implementation to further explore the ET problem, now
                !I set soil order to 0 for all soils. Jinyun Tang, Mar 20, 2014

                ipedof=get_ipedof(0)
                call pedotransf(ipedof, sand, clay, &
                     this%watsat_col(c,lev), this%bsw_col(c,lev), this%sucsat_col(c,lev), xksat)

                om_watsat         = max(0.93_r8 - 0.1_r8   *(zsoi(lev)/zsapric), 0.83_r8)
                om_b              = min(2.7_r8  + 9.3_r8   *(zsoi(lev)/zsapric), 12.0_r8)
                om_sucsat         = min(10.3_r8 - 0.2_r8   *(zsoi(lev)/zsapric), 10.1_r8)
                om_hksat          = max(0.28_r8 - 0.2799_r8*(zsoi(lev)/zsapric), 0.0001_r8)

                this%bd_col(c,lev)        = (1._r8 - this%watsat_col(c,lev))*2.7e3_r8 
                this%watsat_col(c,lev)    = (1._r8 - om_frac) * this%watsat_col(c,lev) + om_watsat*om_frac
                tkm                       = (1._r8-om_frac) * (8.80_r8*sand+2.92_r8*clay)/(sand+clay)+om_tkm*om_frac ! W/(m K)
                this%bsw_col(c,lev)       = (1._r8-om_frac) * (2.91_r8 + 0.159_r8*clay) + om_frac*om_b   
                this%sucsat_col(c,lev)    = (1._r8-om_frac) * this%sucsat_col(c,lev) + om_sucsat*om_frac  
                this%hksat_min_col(c,lev) = xksat

                ! perc_frac is zero unless perf_frac greater than percolation threshold
                if (om_frac > pcalpha) then
                   perc_norm=(1._r8 - pcalpha)**(-pcbeta)
                   perc_frac=perc_norm*(om_frac - pcalpha)**pcbeta
                else
                   perc_frac=0._r8
                endif

                ! uncon_frac is fraction of mineral soil plus fraction of "nonpercolating" organic soil
                uncon_frac=(1._r8-om_frac)+(1._r8-perc_frac)*om_frac

                ! uncon_hksat is series addition of mineral/organic conductivites
                if (om_frac < 1._r8) then
                   uncon_hksat=uncon_frac/((1._r8-om_frac)/xksat &
                        +((1._r8-perc_frac)*om_frac)/om_hksat)
                else
                   uncon_hksat = 0._r8
                end if
                this%hksat_col(c,lev)  = uncon_frac*uncon_hksat + (perc_frac*om_frac)*om_hksat

                this%tkmg_col(c,lev)   = tkm ** (1._r8- this%watsat_col(c,lev))           

                this%tksatu_col(c,lev) = this%tkmg_col(c,lev)*0.57_r8**this%watsat_col(c,lev)

                this%tkdry_col(c,lev)  = ((0.135_r8*this%bd_col(c,lev) + 64.7_r8) / &
                     (2.7e3_r8 - 0.947_r8*this%bd_col(c,lev)))*(1._r8-om_frac) + om_tkd*om_frac  

                this%csol_col(c,lev)   = ((1._r8-om_frac)*(2.128_r8*sand+2.385_r8*clay) / (sand+clay) + &
                     om_csol*om_frac)*1.e6_r8  ! J/(m3 K)

                if (lev > nlevsoi) then
                   this%csol_col(c,lev) = csol_bedrock
                endif

                this%watdry_col(c,lev) = this%watsat_col(c,lev) * &
                     (316230._r8/this%sucsat_col(c,lev)) ** (-1._r8/this%bsw_col(c,lev)) 
                this%watopt_col(c,lev) = this%watsat_col(c,lev) * &
                     (158490._r8/this%sucsat_col(c,lev)) ** (-1._r8/this%bsw_col(c,lev)) 

                !! added by K.Sakaguchi for beta from Lee and Pielke, 1992
                ! water content at field capacity, defined as hk = 0.1 mm/day
                ! used eqn (7.70) in CLM3 technote with k = 0.1 (mm/day) / secspday (day/sec)
                this%watfc_col(c,lev) = this%watsat_col(c,lev) * &
                     (0.1_r8 / (this%hksat_col(c,lev)*secspday))**(1._r8/(2._r8*this%bsw_col(c,lev)+3._r8))
             end if
          end do

          ! Urban pervious and impervious road
          if (col%itype(c) == icol_road_imperv) then
             ! Impervious road layers -- same as above except set watdry and watopt as missing
             do lev = 1,nlevgrnd
                this%watdry_col(c,lev) = spval 
                this%watopt_col(c,lev) = spval 
             end do
          else if (col%itype(c) == icol_road_perv) then 
             ! pervious road layers  - set in UrbanInitTimeConst
          end if

       end if
    end do

    ! --------------------------------------------------------------------
    ! Set soil hydraulic and thermal properties: lake
    ! --------------------------------------------------------------------

    do c = bounds%begc, bounds%endc
       g = col%gridcell(c)
       l = col%landunit(c)

       if (lun%itype(l)==istdlak) then

          do lev = 1,nlevgrnd
             if ( lev <= nlevsoi )then
                clay    =  this%cellclay_col(c,lev)
                sand    =  this%cellsand_col(c,lev)
                om_frac = (this%cellorg_col(c,lev)/organic_max)**2._r8
             else
                clay    = this%cellclay_col(c,nlevsoi)
                sand    = this%cellsand_col(c,nlevsoi)
                om_frac = 0.0_r8
             end if

             this%watsat_col(c,lev) = 0.489_r8 - 0.00126_r8*sand
             this%bsw_col(c,lev)    = 2.91 + 0.159*clay
             this%sucsat_col(c,lev) = 10._r8 * ( 10._r8**(1.88_r8-0.0131_r8*sand) )
             bd                     = (1._r8-this%watsat_col(c,lev))*2.7e3_r8
             this%watsat_col(c,lev) = (1._r8 - om_frac)*this%watsat_col(c,lev) + om_watsat_lake * om_frac
             tkm                    = (1._r8-om_frac)*(8.80_r8*sand+2.92_r8*clay)/(sand+clay) + om_tkm * om_frac ! W/(m K)
             this%bsw_col(c,lev)    = (1._r8-om_frac)*(2.91_r8 + 0.159_r8*clay) + om_frac * om_b_lake
             this%sucsat_col(c,lev) = (1._r8-om_frac)*this%sucsat_col(c,lev) + om_sucsat_lake * om_frac
             xksat                  = 0.0070556 *( 10.**(-0.884+0.0153*sand) ) ! mm/s

             ! perc_frac is zero unless perf_frac greater than percolation threshold
             if (om_frac > pc_lake) then
                perc_norm = (1._r8 - pc_lake)**(-pcbeta)
                perc_frac = perc_norm*(om_frac - pc_lake)**pcbeta
             else
                perc_frac = 0._r8
             endif

             ! uncon_frac is fraction of mineral soil plus fraction of "nonpercolating" organic soil
             uncon_frac = (1._r8-om_frac) + (1._r8-perc_frac)*om_frac

             ! uncon_hksat is series addition of mineral/organic conductivites
             if (om_frac < 1._r8) then
                xksat = 0.0070556 *( 10.**(-0.884+0.0153*sand) ) ! mm/s
                uncon_hksat = uncon_frac/((1._r8-om_frac)/xksat + ((1._r8-perc_frac)*om_frac)/om_hksat_lake)
             else
                uncon_hksat = 0._r8
             end if

             this%hksat_col(c,lev)  = uncon_frac*uncon_hksat + (perc_frac*om_frac)*om_hksat_lake
             this%tkmg_col(c,lev)   = tkm ** (1._r8- this%watsat_col(c,lev))
             this%tksatu_col(c,lev) = this%tkmg_col(c,lev)*0.57_r8**this%watsat_col(c,lev)
             this%tkdry_col(c,lev)  = ((0.135_r8*bd + 64.7_r8) / (2.7e3_r8 - 0.947_r8*bd))*(1._r8-om_frac) + &
                                       om_tkd * om_frac
             this%csol_col(c,lev)   = ((1._r8-om_frac)*(2.128_r8*sand+2.385_r8*clay) / (sand+clay) +   &
                                       om_csol * om_frac)*1.e6_r8  ! J/(m3 K)
             if (lev > nlevsoi) then
                this%csol_col(c,lev) = csol_bedrock
             endif

             this%watdry_col(c,lev) = this%watsat_col(c,lev) * (316230._r8/this%sucsat_col(c,lev)) ** (-1._r8/this%bsw_col(c,lev))
             this%watopt_col(c,lev) = this%watsat_col(c,lev) * (158490._r8/this%sucsat_col(c,lev)) ** (-1._r8/this%bsw_col(c,lev))

             !! added by K.Sakaguchi for beta from Lee and Pielke, 1992
             ! water content at field capacity, defined as hk = 0.1 mm/day
             ! used eqn (7.70) in CLM3 technote with k = 0.1 (mm/day) / (# seconds/day)
             this%watfc_col(c,lev) = this%watsat_col(c,lev) * (0.1_r8 / &
                               (this%hksat_col(c,lev)*secspday))**(1._r8/(2._r8*this%bsw_col(c,lev)+3._r8))
          end do
       endif

    end do

    ! --------------------------------------------------------------------
    ! Initialize threshold soil moisture and mass fracion of clay limited to 0.20
    ! --------------------------------------------------------------------

    do c = bounds%begc, bounds%endc
       g = col%gridcell(c)

       this%gwc_thr_col(c) = 0.17_r8 + 0.14_r8 * clay3d(g,1) * 0.01_r8
       this%mss_frc_cly_vld_col(c) = min(clay3d(g,1) * 0.01_r8, 0.20_r8)
    end do

    ! --------------------------------------------------------------------
    ! Deallocate memory
    ! --------------------------------------------------------------------

    deallocate(sand3d, clay3d, organic3d)
    deallocate(zisoifl, zsoifl, dzsoifl)



  end subroutine initcold_col_pp


  !------------------------------------------------------------------------
  subroutine restart_col_pp(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use ncdio_pio  , only : file_desc_t, ncd_defvar, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
    use restUtilMod
    
!   moved from SoilStateType.F90

    ! !USES:
    use shr_log_mod, only : errMsg => shr_log_errMsg
    use spmdMod    , only : masterproc
    use abortutils , only : endrun
    use restUtilMod
    use ncdio_pio
    use clm_varctl,  only : use_dynroot
    use RootBiophysMod      , only : init_vegrootfr
    !


    ! !ARGUMENTS:
    class(column_physical_properties) :: this
    type(bounds_type), intent(in)    :: bounds  
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer :: j,c ! indices
    logical :: readvar      ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='LAKE_ICEFRAC', xtype=ncd_double,  &
         dim1name='column', dim2name='levlak', switchdim=.true., &
         long_name='lake layer ice fraction', units='kg/kg', &
         interpinic_flag='interp', readvar=readvar, data=this%lake_icefrac_col)

    call restartvar(ncid=ncid, flag=flag, varname='SAVEDTKE1', xtype=ncd_double,  &
         dim1name='column', &
         long_name='top lake layer eddy conductivity', units='W/(m K)', &
         interpinic_flag='interp', readvar=readvar, data=this%savedtke1_col)

    call restartvar(ncid=ncid, flag=flag, varname='USTLAKE', xtype=ncd_double,  &
         dim1name='column', &
         long_name='friction velocity for lakes', units='m/s', &
         interpinic_flag='interp', readvar=readvar, data=this%ust_lake_col)

!!!   *********DW********** moved (use_dynroot) into patch level datastructure
!!!    if(use_dynroot) then
!!!    end if

  end subroutine restart_col_pp


  subroutine clean_col_pp(this)
    class(column_physical_properties) :: this
  
! !ARGUMENTS:
!    class(column_physical_properties) :: this
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
    deallocate(this%soil_pH   )

    deallocate(this%cascade_step_name    )
    deallocate(this%cascade_donor_pool   )
    deallocate(this%cascade_receiver_pool)

    !-- properties of each decomposing pool
    deallocate(this%floating_cn_ratio_decomp_pools          )
    deallocate(this%floating_cp_ratio_decomp_pools          )
    deallocate(this%decomp_pool_name_restart(0:ndecomp_pools))
    deallocate(this%decomp_pool_name_history(0:ndecomp_pools))
    deallocate(this%decomp_pool_name_long(0:ndecomp_pools))
    deallocate(this%decomp_pool_name_short(0:ndecomp_pools))
    deallocate(this%is_litter(0:ndecomp_pools))
    deallocate(this%is_soil(0:ndecomp_pools))
    deallocate(this%is_cwd(0:ndecomp_pools))
    deallocate(this%initial_cn_ratio(0:ndecomp_pools))
    deallocate(this%initial_cp_ratio(0:ndecomp_pools))
    deallocate(this%initial_stock(0:ndecomp_pools))
    deallocate(this%is_metabolic(0:ndecomp_pools))
    deallocate(this%is_cellulose(0:ndecomp_pools))
    deallocate(this%is_lignin(0:ndecomp_pools))
    deallocate(this%spinup_factor(0:ndecomp_pools))

    deallocate(this%etal_col           (begc:endc))          
    deallocate(this%lakefetch_col      (begc:endc))           
    deallocate(this%lakeresist_col     (begc:endc))           
    deallocate(this%savedtke1_col      (begc:endc))           
    deallocate(this%lake_icefrac_col   (begc:endc,1:nlevlak)) 
    deallocate(this%lake_icethick_col  (begc:endc))           
    deallocate(this%ust_lake_col       (begc:endc))           
    deallocate(this%ram1_lake_patch      (begp:endp))         
    deallocate(this%lake_raw_col       (begc:endc))           
    deallocate(this%ks_col             (begc:endc))          
    deallocate(this%ws_col             (begc:endc))          
    deallocate(this%betaprime_col      (begc:endc))  


!  moved from SoilStateType.F90

    deallocate(this%mss_frc_cly_vld_col)                                   
    deallocate(this%cellorg_col        )          
    allocate(this%cellsand_col         )            
    allocate(this%cellclay_col         )
    allocate(this%bd_col               )

    allocate(this%hksat_col            )
    allocate(this%hksat_min_col        )
    allocate(this%hk_l_col             )  
    allocate(this%smp_l_col            )  
    allocate(this%smpmin_col           )

    allocate(this%bsw_col              )
    allocate(this%watsat_col           )
    allocate(this%watdry_col           )
    allocate(this%watopt_col           )
    allocate(this%watfc_col            )
    allocate(this%sucsat_col           )
    allocate(this%soilbeta_col         )  
    allocate(this%soilalpha_col        )
    allocate(this%soilalpha_u_col      )
    allocate(this%soilpsi_col          )
    allocate(this%wtfact_col           )
    allocate(this%porosity_col         )
    allocate(this%eff_porosity_col     )
    allocate(this%gwc_thr_col          )

    allocate(this%thk_col              )
    allocate(this%tkmg_col             )
    allocate(this%tkdry_col            )
    allocate(this%tksatu_col           )
    allocate(this%csol_col             )

  !  allocate(this%rootr_patch         )
    allocate(this%rootr_col            )
    allocate(this%rootr_road_perv_col  )
  !  allocate(this%rootfr_patch        )
    allocate(this%rootfr_col           )
    allocate(this%rootfr_road_perv_col )
  !  allocate(this%root_depth_patch    )

    deallocate(this%smax           )
    deallocate(this%ks_sorption    )
    deallocate(this%r_weather      )
    deallocate(this%r_adsorp       )
    deallocate(this%r_desorp       )
    deallocate(this%r_occlude      )

    deallocate(this%k_s1_biochem      )
    deallocate(this%k_s2_biochem      )
    deallocate(this%k_s3_biochem      )
    deallocate(this%k_s4_biochem      )

  end subroutine clean_col_pp

end module ColumnPhysicalPropertiesType


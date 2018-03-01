module FatesInterfaceMod

   ! ------------------------------------------------------------------------------------
   ! This is the FATES public API
   ! A host land model has defined and allocated a structure "fates" as
   ! defined by fates_interface_type
   !
   ! It is also likely/possible that this type is defined as a vector
   ! which is allocated by thread
   ! ------------------------------------------------------------------------------------

   use EDTypesMod          , only : ed_site_type
   use EDTypesMod          , only : maxPatchesPerSite
   use EDTypesMod          , only : maxCohortsPerPatch
   use EDTypesMod          , only : maxSWb
   use EDTypesMod          , only : ivis
   use EDTypesMod          , only : inir
   use EDTypesMod          , only : nclmax
   use EDTypesMod          , only : nlevleaf
   use EDTypesMod          , only : maxpft
   use FatesConstantsMod   , only : r8 => fates_r8
   use FatesConstantsMod   , only : itrue,ifalse
   use FatesGlobals        , only : fates_global_verbose
   use FatesGlobals        , only : fates_log
   use FatesGlobals        , only : endrun => fates_endrun
   use EDPftvarcon         , only : FatesReportPFTParams
   use EDPftvarcon         , only : EDPftvarcon_inst
   use EDParamsMod         , only : FatesReportParams


   ! CIME Globals
   use shr_log_mod         , only : errMsg => shr_log_errMsg
   use shr_infnan_mod      , only : nan => shr_infnan_nan, assignment(=)

   implicit none

   public :: FatesInterfaceInit
   public :: set_fates_ctrlparms
   public :: SetFatesTime
   public :: set_fates_global_elements
   public :: FatesReportParameters

   character(len=*), parameter, private :: sourcefile = &
         __FILE__
   
   ! -------------------------------------------------------------------------------------
   ! Parameters that are dictated by the Host Land Model
   ! THESE ARE NOT DYNAMIC. SHOULD BE SET ONCE DURING INTIALIZATION.
   ! -------------------------------------------------------------------------------------

  
   integer, protected :: hlm_numSWb  ! Number of broad-bands in the short-wave radiation
                                     ! specturm to track 
                                     ! (typically 2 as a default, VIS/NIR, in ED variants <2016)

   integer, protected :: hlm_ivis    ! The HLMs assumption of the array index associated with the 
                                     ! visible portion of the spectrum in short-wave radiation arrays

   integer, protected :: hlm_inir    ! The HLMs assumption of the array index associated with the 
                                     ! NIR portion of the spectrum in short-wave radiation arrays


   integer, protected :: hlm_numlevgrnd   ! Number of ground layers
   integer, protected :: hlm_numlevsoil   ! Number of soil layers

   
   integer, protected :: hlm_numlevdecomp_full ! Number of GROUND layers for the purposes
                                               ! of biogeochemistry; can be either 1 
                                               ! or the total number of soil layers
                                               ! (includes bedrock)
   
   
   integer, protected :: hlm_numlevdecomp ! Number of SOIL layers for the purposes of 
                                          ! biogeochemistry; can be either 1 or the total
                                          ! number of soil layers

   integer, protected :: hlm_is_restart   ! Is the HLM signalling that this is a restart
                                          ! type simulation?
                                          ! 1=TRUE, 0=FALSE
   
   character(len=16), protected :: hlm_name ! This character string passed by the HLM
                                            ! is used during the processing of IO data, 
                                            ! so that FATES knows which IO variables it 
                                            ! should prepare.  For instance
                                            ! ATS, ALM and CLM will only want variables 
                                            ! specficially packaged for them.
                                            ! This string sets which filter is enacted.
   
  
   real(r8), protected :: hlm_hio_ignore_val  ! This value can be flushed to history 
                                              ! diagnostics, such that the
                                              ! HLM will interpret that the value should not 
                                              ! be included in the average.
   
   integer, protected :: hlm_masterproc  ! Is this the master processor, typically useful
                                         ! for knowing if the current machine should be 
                                         ! printing out messages to the logs or terminals
                                         ! 1 = TRUE (is master) 0 = FALSE (is not master)

   integer, protected :: hlm_ipedof      ! The HLM pedotransfer index
                                         ! this is only used by the plant hydraulics
                                         ! submodule to check and/or enable consistency
                                         ! between the pedotransfer functions of the HLM
                                         ! and how it moves and stores water in its
                                         ! rhizosphere shells
   
   integer, protected :: hlm_max_patch_per_site ! The HLM needs to exchange some patch
                                                ! level quantities with FATES
                                                ! FATES does not dictate those allocations
                                                ! since it happens pretty early in
                                                ! the model initialization sequence.
                                                ! So we want to at least query it,
                                                ! compare it to our maxpatchpersite,
                                                ! and gracefully halt if we are over-allocating

   integer, protected :: hlm_use_vertsoilc ! This flag signals whether or not the 
                                           ! host model is using vertically discretized
                                           ! soil carbon
                                           ! 1 = TRUE,  0 = FALSE
   
   integer, protected :: hlm_use_spitfire  ! This flag signals whether or not to use SPITFIRE
                                           ! 1 = TRUE, 0 = FALSE


   integer, protected :: hlm_use_logging       ! This flag signals whether or not to use
                                               ! the logging module

   integer, protected :: hlm_use_planthydro    ! This flag signals whether or not to use
                                               ! plant hydraulics (bchristo/xu methods)
                                               ! 1 = TRUE, 0 = FALSE
                                               ! THIS IS CURRENTLY NOT SUPPORTED 

   integer, protected :: hlm_use_ed_st3        ! This flag signals whether or not to use
                                               ! (ST)atic (ST)and (ST)ructure mode (ST3)
                                               ! Essentially, this gives us the ability
                                               ! to turn off "dynamics", ie growth, disturbance
                                               ! recruitment and mortality.
                                               ! (EXPERIMENTAL!!!!! - RGK 07-2017)
                                               ! 1 = TRUE, 0 = FALSE
                                               ! default should be FALSE (dynamics on)
                                               ! cannot be true with prescribed_phys

   integer, protected :: hlm_use_ed_prescribed_phys ! This flag signals whether or not to use
                                                    ! prescribed physiology, somewhat the opposite
                                                    ! to ST3, in this case can turn off
                                                    ! fast processes like photosynthesis and respiration
                                                    ! and prescribe NPP
                                                    ! (NOT CURRENTLY IMPLEMENTED - PLACEHOLDER)
                                                    ! 1 = TRUE, 0 = FALSE
                                                    ! default should be FALSE (biophysics on)
                                                    ! cannot be true with st3 mode

   integer, protected :: hlm_use_inventory_init     ! Initialize this simulation from
                                                    ! an inventory file. If this is toggled on
                                                    ! an inventory control file must be specified
                                                    ! as well.
                                                    ! 1 = TRUE, 0 = FALSE
   
   character(len=256), protected :: hlm_inventory_ctrl_file ! This is the full path to the
                                                            ! inventory control file that
                                                            ! specifieds the availabel inventory datasets
                                                            ! there locations and their formats
                                                            ! This need only be defined when
                                                            ! hlm_use_inventory_init = 1

   ! -------------------------------------------------------------------------------------
   ! Parameters that are dictated by FATES and known to be required knowledge
   !  needed by the HLMs
   ! -------------------------------------------------------------------------------------

   ! Variables mostly used for dimensioning host land model (HLM) array spaces
   
   integer, protected :: fates_maxElementsPerPatch ! maxElementsPerPatch is the value that is ultimately
                                                   ! used to set the size of the largest arrays necessary
                                                   ! in things like restart files (probably hosted by the 
                                                   ! HLM). The size of these arrays are not a parameter
                                                   ! because it is simply the maximum of several different
                                                   ! dimensions. It is possible that this would be the
                                                   ! maximum number of cohorts per patch, but
                                                   ! but it could be other things.

   integer, protected :: fates_maxElementsPerSite  ! This is the max number of individual items one can store per 
                                                   ! each grid cell and effects the striding in the ED restart 
                                                   ! data as some fields are arrays where each array is
                                                   ! associated with one cohort

   ! -------------------------------------------------------------------------------------
   ! These vectors are used for history output mapping
   ! CLM/ALM have limited support for multi-dimensional history output arrays.
   ! FATES structure and composition is multi-dimensional, so we end up "multi-plexing"
   ! multiple dimensions into one dimension.  These new dimensions need definitions,
   ! mapping to component dimensions, and definitions for those component dimensions as
   ! well.
   ! -------------------------------------------------------------------------------------
   
   real(r8), allocatable :: fates_hdim_levsclass(:)        ! plant size class lower bound dimension
   integer , allocatable :: fates_hdim_pfmap_levscpf(:)    ! map of pfts into size-class x pft dimension
   integer , allocatable :: fates_hdim_scmap_levscpf(:)    ! map of size-class into size-class x pft dimension
   real(r8), allocatable :: fates_hdim_levage(:)           ! patch age lower bound dimension
   integer , allocatable :: fates_hdim_levpft(:)           ! plant pft dimension
   integer , allocatable :: fates_hdim_levfuel(:)          ! fire fuel class dimension
   integer , allocatable :: fates_hdim_levcwdsc(:)         ! cwd class dimension
   integer , allocatable :: fates_hdim_levcan(:)           ! canopy-layer dimension 
   integer , allocatable :: fates_hdim_canmap_levcnlf(:)   ! canopy-layer map into the canopy-layer x leaf-layer dim
   integer , allocatable :: fates_hdim_lfmap_levcnlf(:)    ! leaf-layer map into the can-layer x leaf-layer dimension
   integer , allocatable :: fates_hdim_canmap_levcnlfpf(:) ! can-layer map into the can-layer x pft x leaf-layer dim
   integer , allocatable :: fates_hdim_lfmap_levcnlfpf(:)  ! leaf-layer map into the can-layer x pft x leaf-layer dim
   integer , allocatable :: fates_hdim_pftmap_levcnlfpf(:) ! pft map into the canopy-layer x pft x leaf-layer dim
   integer , allocatable :: fates_hdim_scmap_levscag(:)    ! map of size-class into size-class x patch age dimension
   integer , allocatable :: fates_hdim_agmap_levscag(:)    ! map of patch-age into size-class x patch age dimension

   ! ------------------------------------------------------------------------------------
   !                              DYNAMIC BOUNDARY CONDITIONS
   ! ------------------------------------------------------------------------------------


   ! -------------------------------------------------------------------------------------
   ! Scalar Timing Variables
   ! It is assumed that all of the sites on a given machine will be synchronous.
   ! It is also assumed that the HLM will control time.
   ! -------------------------------------------------------------------------------------
   integer, protected  :: hlm_current_year    ! Current year
   integer, protected  :: hlm_current_month   ! month of year
   integer, protected  :: hlm_current_day     ! day of month
   integer, protected  :: hlm_current_tod     ! time of day (seconds past 0Z)
   integer, protected  :: hlm_current_date    ! time of day (seconds past 0Z)
   integer, protected  :: hlm_reference_date  ! YYYYMMDD
   real(r8), protected :: hlm_model_day       ! elapsed days between current date and ref
   integer, protected  :: hlm_day_of_year     ! The integer day of the year
   integer, protected  :: hlm_days_per_year   ! The HLM controls time, some HLMs may 
                                              ! include a leap
   real(r8), protected :: hlm_freq_day        ! fraction of year for daily time-step 
                                              ! (1/days_per_year_, this is a frequency
   

   ! -------------------------------------------------------------------------------------
   !
   ! Constant parameters that are dictated by the fates parameter file
   !
   ! -------------------------------------------------------------------------------------

   integer, protected :: numpft          ! The total number of PFTs defined in the simulation
   integer, protected :: nlevsclass   ! The total number of cohort size class bins output to history
   integer, protected :: nlevage      ! The total number of patch age bins output to history
   

   ! -------------------------------------------------------------------------------------
   ! Structured Boundary Conditions (SITE/PATCH SCALE)
   ! For floating point arrays, it is sometimes the convention to define the arrays as
   ! POINTER instead of ALLOCATABLE.  This usually achieves the same result with subtle
   ! differences.  POINTER arrays can point to scalar values, discontinuous array slices
   ! or alias other variables, ALLOCATABLES cannnot.  According to S. Lionel 
   ! (Intel-Forum Post), ALLOCATABLES are better perfomance wise as long as they point 
   ! to contiguous memory spaces and do not alias other variables, the case here.
   ! Naming conventions:   _gl  means ground layer dimensions
   !                       _si  means site dimensions (scalar in that case)
   !                       _pa  means patch dimensions
   !                       _rb  means radiation band
   ! ------------------------------------------------------------------------------------

   type, public :: bc_in_type

      ! The actual number of FATES' ED patches
      integer :: npatches


      ! Soil layer structure
      real(r8),allocatable :: zi_sisl(:)         ! interface level below a "z" level (m)
                                                 ! this contains a zero index for surface.
      real(r8),allocatable :: dz_sisl(:)         ! layer thickness (m)
      real(r8),allocatable :: z_sisl(:)          ! layer depth (m) (1:hlm_nlevsoil) 

      ! Decomposition Layer Structure
      real(r8), allocatable :: dz_decomp_sisl(:)

      ! Vegetation Dynamics
      ! ---------------------------------------------------------------------------------

      ! The site level 24 hour vegetation temperature is used for various purposes during vegetation 
      ! dynamics.  However, we are currently using the bare ground patch's value [K]
      ! TO-DO: Get some consensus on the correct vegetation temperature used for phenology.
      ! It is possible that the bare-ground value is where the average is being stored.
      ! (RGK-01-2017)
      real(r8)             :: t_veg24_si

      ! Patch 24 hour vegetation temperature [K]
      real(r8),allocatable :: t_veg24_pa(:)  
      
      ! Fire Model

      ! Average precipitation over the last 24 hours [mm/s]
      real(r8), allocatable :: precip24_pa(:)

      ! Average relative humidity over past 24 hours [-]
      real(r8), allocatable :: relhumid24_pa(:)

      ! Patch 24-hour running mean of wind (m/s ?)
      real(r8), allocatable :: wind24_pa(:)


      ! Radiation variables for calculating sun/shade fractions
      ! ---------------------------------------------------------------------------------

      ! Downwelling direct beam radiation (patch,radiation-band) [W/m2]
      real(r8), allocatable :: solad_parb(:,:)  

      ! Downwelling diffuse (I-ndirect) radiation (patch,radiation-band) [W/m2]
      real(r8), allocatable :: solai_parb(:,:)



      ! Photosynthesis variables
      ! ---------------------------------------------------------------------------------

      ! Patch level filter flag for photosynthesis calculations
      ! has a short memory, flags:
      ! 1 = patch has not been called
      ! 2 = patch is currently marked for photosynthesis
      ! 3 = patch has been called for photosynthesis at least once
      integer, allocatable  :: filter_photo_pa(:)

      ! atmospheric pressure (Pa)
      real(r8)              :: forc_pbot             

      ! daylength scaling factor (0-1)
      real(r8), allocatable :: dayl_factor_pa(:)
      
      ! saturation vapor pressure at t_veg (Pa)
      real(r8), allocatable :: esat_tv_pa(:)

      ! vapor pressure of canopy air (Pa)
      real(r8), allocatable :: eair_pa(:)

      ! Atmospheric O2 partial pressure (Pa)
      real(r8), allocatable :: oair_pa(:)

      ! Atmospheric CO2 partial pressure (Pa)
      real(r8), allocatable :: cair_pa(:)

      ! boundary layer resistance (s/m)
      real(r8), allocatable :: rb_pa(:)

      ! vegetation temperature (Kelvin)
      real(r8), allocatable :: t_veg_pa(:)
             
      ! air temperature at agcm reference height (kelvin)
      real(r8), allocatable :: tgcm_pa(:)

      ! soil temperature (Kelvin)
      real(r8), allocatable :: t_soisno_gl(:)

      ! Canopy Radiation Boundaries
      ! ---------------------------------------------------------------------------------
      
      ! Filter for vegetation patches with a positive zenith angle (daylight)
      logical, allocatable :: filter_vegzen_pa(:)

      ! Cosine of the zenith angle (0-1), by patch
      ! Note RGK: It does not seem like the code would currently generate
      !           different zenith angles for different patches (nor should it)
      !           I am leaving it at this scale for simplicity.  Patches should
      !           have no spacially variable information
      real(r8), allocatable :: coszen_pa(:)
      
      ! Abledo of the ground for direct radiation, by site broadband (0-1)
      real(r8), allocatable :: albgr_dir_rb(:)

      ! Albedo of the ground for diffuse radiation, by site broadband (0-1)
      real(r8), allocatable :: albgr_dif_rb(:)
      
      ! LitterFlux Boundaries
      ! the index of the deepest model soil level where roots may be
      ! due to permafrost or bedrock constraints
      integer  :: max_rooting_depth_index_col

      ! BGC Accounting

      real(r8) :: tot_het_resp  ! total heterotrophic respiration  (gC/m2/s)
      real(r8) :: tot_somc      ! total soil organic matter carbon (gc/m2)
      real(r8) :: tot_litc      ! total litter carbon tracked in the HLM (gc/m2)

      ! Canopy Structure

      real(r8) :: snow_depth_si    ! Depth of snow in snowy areas of site (m)
      real(r8) :: frac_sno_eff_si  ! Fraction of ground covered by snow (0-1)

      ! Hydrology variables for BTRAN
      ! ---------------------------------------------------------------------------------

      ! Soil suction potential of layers in each site, negative, [mm]
      real(r8), allocatable :: smp_gl(:)

      ! Effective porosity = porosity - vol_ic, of layers in each site [-]
      real(r8), allocatable :: eff_porosity_gl(:)

      ! volumetric soil water at saturation (porosity)
      real(r8), allocatable :: watsat_gl(:)

      ! Temperature of ground layers [K]
      real(r8), allocatable :: tempk_gl(:)

      ! Liquid volume in ground layer (m3/m3)
      real(r8), allocatable :: h2o_liqvol_gl(:)

      ! Site level filter for uptake response functions
      logical               :: filter_btran

      ! Plant-Hydro
      ! ---------------------------------------------------------------------------------

      
      real(r8),allocatable :: qflx_transp_pa(:)    ! Transpiration flux as dictated by the HLM's
                                                   ! canopy solver. [mm H2O/s] [+ into root]
      real(r8),allocatable :: swrad_net_pa(:)      ! Net absorbed shortwave radiation (W/m2)
      real(r8),allocatable :: lwrad_net_pa(:)      ! Net absorbed longwave radiation (W/m2)
      real(r8),allocatable :: watsat_sisl(:)       ! volumetric soil water at saturation (porosity)
      real(r8),allocatable :: watres_sisl(:)       ! volumetric residual soil water
      real(r8),allocatable :: sucsat_sisl(:)       ! minimum soil suction (mm) (hlm_nlevsoil) 
      real(r8),allocatable :: bsw_sisl(:)          ! Clapp and Hornberger "b" (hlm_nlevsoil)
      real(r8),allocatable :: hksat_sisl(:)        ! hydraulic conductivity at saturation (mm H2O /s)
      real(r8),allocatable :: h2o_liq_sisl(:)      ! Liquid water mass in each layer (kg/m2)
      real(r8) :: smpmin_si                        ! restriction for min of soil potential (mm)
      
   end type bc_in_type


   type, public :: bc_out_type

      ! Sunlit fraction of the canopy for this patch [0-1]
      real(r8),allocatable :: fsun_pa(:)

      ! Sunlit canopy LAI
      real(r8),allocatable :: laisun_pa(:)
      
      ! Shaded canopy LAI
      real(r8),allocatable :: laisha_pa(:)
      
      ! Logical stating whether a ground layer can have water uptake by plants
      ! The only condition right now is that liquid water exists
      ! The name (suction) is used to indicate that soil suction should be calculated
      logical, allocatable :: active_suction_gl(:)

      ! Effective fraction of roots in each soil layer 
      real(r8), allocatable :: rootr_pagl(:,:)

      ! Integrated (vertically) transpiration wetness factor (0 to 1) 
      ! (diagnostic, should not be used by HLM)
      real(r8), allocatable :: btran_pa(:)

      ! Sunlit canopy resistance [s/m]
      real(r8), allocatable :: rssun_pa(:)

      ! Shaded canopy resistance [s/m]
      real(r8), allocatable :: rssha_pa(:)

      ! leaf photosynthesis (umol CO2 /m**2/ s)
      ! (NOT CURRENTLY USED, PLACE-HOLDER)
      !real(r8), allocatable :: psncanopy_pa(:)

      ! leaf maintenance respiration rate (umol CO2/m**2/s) 
      ! (NOT CURRENTLY USED, PLACE-HOLDER)
      !real(r8), allocatable :: lmrcanopy_pa(:)

      ! Canopy Radiation Boundaries
      ! ---------------------------------------------------------------------------------
      
      ! Surface albedo (direct) (HLMs use this for atm coupling and balance checks)
      real(r8), allocatable :: albd_parb(:,:)
      
      ! Surface albedo (diffuse) (HLMs use this for atm coupling and balance checks)
      real(r8), allocatable :: albi_parb(:,:)                 
      
      ! Flux absorbed by canopy per unit direct flux (HLMs use this for balance checks)
      real(r8), allocatable :: fabd_parb(:,:) 
      
      ! Flux absorbed by canopy per unit diffuse flux (HLMs use this for balance checks)
      real(r8), allocatable :: fabi_parb(:,:)

      ! Down direct flux below canopy per unit direct flx (HLMs use this for balance checks)
      real(r8), allocatable :: ftdd_parb(:,:)

      ! Down diffuse flux below canopy per unit direct flx (HLMs use this for balance checks)
      real(r8), allocatable :: ftid_parb(:,:)
      
      ! Down diffuse flux below canopy per unit diffuse flx (HLMs use this for balance checks)
      real(r8), allocatable :: ftii_parb(:,:)


      ! litterfall fluxes of C from FATES patches to BGC columns

      ! total labile    litter coming from ED. gC/m3/s
      real(r8), allocatable :: FATES_c_to_litr_lab_c_col(:)      

      !total cellulose litter coming from ED. gC/m3/s
      real(r8), allocatable :: FATES_c_to_litr_cel_c_col(:)      
      
      !total lignin    litter coming from ED. gC/m3/s
      real(r8), allocatable :: FATES_c_to_litr_lig_c_col(:)      

      

      ! Canopy Structure

      real(r8), allocatable :: elai_pa(:)  ! exposed leaf area index
      real(r8), allocatable :: esai_pa(:)  ! exposed stem area index
      real(r8), allocatable :: tlai_pa(:)  ! total leaf area index
      real(r8), allocatable :: tsai_pa(:)  ! total stem area index
      real(r8), allocatable :: htop_pa(:)  ! top of the canopy [m]
      real(r8), allocatable :: hbot_pa(:)  ! bottom of canopy? [m]

      real(r8), allocatable :: z0m_pa(:)   ! roughness length [m]
      real(r8), allocatable :: displa_pa(:) ! displacement height [m]
      real(r8), allocatable :: dleaf_pa(:)  ! leaf characteristic dimension/width/diameter [m]

      real(r8), allocatable :: canopy_fraction_pa(:) ! Area fraction of each patch in the site
                                                     ! Use most likely for weighting
                                                     ! This is currently the projected canopy
                                                     ! area of each patch [0-1]

      real(r8), allocatable :: frac_veg_nosno_alb_pa(:) ! This is not really a fraction
                                                        ! this is actually binary based on if any
                                                        ! vegetation in the patch is exposed.
                                                        ! [0,1]

      ! FATES Hydraulics

      real(r8) :: plant_stored_h2o_si             ! stored water in vegetation (kg/m2 H2O)
                                                  ! Assuming density of 1Mg/m3 ~= mm/m2 H2O
                                                  ! This must be set and transfered prior to clm_drv()
                                                  ! following the calls to ed_update_site()
                                                  ! ed_update_site() is called during both the restart
                                                  ! and coldstart process
      
      real(r8),allocatable :: qflx_soil2root_sisl(:)   ! Water flux from soil into root by site and soil layer
                                                       ! [mm H2O/s] [+ into root]

   end type bc_out_type


   type, public :: fates_interface_type
      
      ! This is the root of the ED/FATES hierarchy of instantaneous state variables
      ! ie the root of the linked lists. Each path list is currently associated with a 
      ! grid-cell, this is intended to be migrated to columns 

      integer                         :: nsites

      type(ed_site_type), pointer :: sites(:)

      ! These are boundary conditions that the FATES models are required to be filled.  
      ! These values are filled by the driver or HLM.  Once filled, these have an 
      ! intent(in) status.  Each site has a derived type structure, which may include 
      ! a scalar for site level data, a patch vector, potentially cohort vectors (but 
      ! not yet atm) and other dimensions such as soil-depth or pft.  These vectors 
      ! are initialized by maximums, and the allocations are static in time to avoid
      ! having to allocate/de-allocate memory

      type(bc_in_type), allocatable   :: bc_in(:)

      ! These are the boundary conditions that the FATES model returns to its HLM or 
      ! driver. It has the same allocation strategy and similar vector types.
      
      type(bc_out_type), allocatable  :: bc_out(:)

   contains
      
      procedure, public :: zero_bcs

   end type fates_interface_type

  


contains

   ! ====================================================================================
  subroutine FatesInterfaceInit(log_unit,global_verbose)

    use FatesGlobals, only : FatesGlobalsInit

    implicit none

    integer, intent(in) :: log_unit
    logical, intent(in) :: global_verbose

    call FatesGlobalsInit(log_unit,global_verbose)

  end subroutine FatesInterfaceInit

   ! ====================================================================================

   ! INTERF-TODO: THIS IS A PLACE-HOLDER ROUTINE, NOT CALLED YET...
   subroutine fates_clean(this)
      
      implicit none
      
      ! Input Arguments
      class(fates_interface_type), intent(inout) :: this
      
      ! Incrementally walk through linked list and deallocate
      
      
      
      ! Deallocate the site list
!      deallocate (this%sites)
      
      return
   end subroutine fates_clean


   ! ====================================================================================
   

   subroutine allocate_bcin(bc_in)
      
      ! ---------------------------------------------------------------------------------
      ! Allocate and Initialze the FATES boundary condition vectors
      ! ---------------------------------------------------------------------------------
      
      implicit none
      type(bc_in_type), intent(inout) :: bc_in
      
      ! Allocate input boundaries
      allocate(bc_in%zi_sisl(0:hlm_numlevsoil))
      allocate(bc_in%dz_sisl(hlm_numlevsoil))
      allocate(bc_in%z_sisl(hlm_numlevsoil))

      allocate(bc_in%dz_decomp_sisl(hlm_numlevdecomp_full))

      ! Vegetation Dynamics
      allocate(bc_in%t_veg24_pa(maxPatchesPerSite))

      allocate(bc_in%wind24_pa(maxPatchesPerSite))
      allocate(bc_in%relhumid24_pa(maxPatchesPerSite))
      allocate(bc_in%precip24_pa(maxPatchesPerSite))
      
      ! Radiation
      allocate(bc_in%solad_parb(maxPatchesPerSite,hlm_numSWb))
      allocate(bc_in%solai_parb(maxPatchesPerSite,hlm_numSWb))
      
      ! Hydrology
      allocate(bc_in%smp_gl(hlm_numlevgrnd))
      allocate(bc_in%eff_porosity_gl(hlm_numlevgrnd))
      allocate(bc_in%watsat_gl(hlm_numlevgrnd))
      allocate(bc_in%tempk_gl(hlm_numlevgrnd))
      allocate(bc_in%h2o_liqvol_gl(hlm_numlevgrnd))

      ! Photosynthesis
      allocate(bc_in%filter_photo_pa(maxPatchesPerSite))
      allocate(bc_in%dayl_factor_pa(maxPatchesPerSite))
      allocate(bc_in%esat_tv_pa(maxPatchesPerSite))
      allocate(bc_in%eair_pa(maxPatchesPerSite))
      allocate(bc_in%oair_pa(maxPatchesPerSite))
      allocate(bc_in%cair_pa(maxPatchesPerSite))
      allocate(bc_in%rb_pa(maxPatchesPerSite))
      allocate(bc_in%t_veg_pa(maxPatchesPerSite))
      allocate(bc_in%tgcm_pa(maxPatchesPerSite))
      allocate(bc_in%t_soisno_gl(hlm_numlevgrnd))

      ! Canopy Radiation
      allocate(bc_in%filter_vegzen_pa(maxPatchesPerSite))
      allocate(bc_in%coszen_pa(maxPatchesPerSite))
      allocate(bc_in%albgr_dir_rb(hlm_numSWb))
      allocate(bc_in%albgr_dif_rb(hlm_numSWb))

      ! Plant-Hydro BC's
      if (hlm_use_planthydro.eq.itrue) then
      
         allocate(bc_in%qflx_transp_pa(maxPatchesPerSite))
         allocate(bc_in%swrad_net_pa(maxPatchesPerSite))
         allocate(bc_in%lwrad_net_pa(maxPatchesPerSite))
         allocate(bc_in%watsat_sisl(hlm_numlevsoil))
         allocate(bc_in%watres_sisl(hlm_numlevsoil))
         allocate(bc_in%sucsat_sisl(hlm_numlevsoil))
         allocate(bc_in%bsw_sisl(hlm_numlevsoil))
         allocate(bc_in%hksat_sisl(hlm_numlevsoil))
         allocate(bc_in%h2o_liq_sisl(hlm_numlevsoil)); bc_in%h2o_liq_sisl = nan
      end if

      return
   end subroutine allocate_bcin
   
   subroutine allocate_bcout(bc_out)

      ! ---------------------------------------------------------------------------------
      ! Allocate and Initialze the FATES boundary condition vectors
      ! ---------------------------------------------------------------------------------
      
      implicit none
      type(bc_out_type), intent(inout) :: bc_out
      
      
      ! Radiation
      allocate(bc_out%fsun_pa(maxPatchesPerSite))
      allocate(bc_out%laisun_pa(maxPatchesPerSite))
      allocate(bc_out%laisha_pa(maxPatchesPerSite))
      
      ! Hydrology
      allocate(bc_out%active_suction_gl(hlm_numlevgrnd))
      allocate(bc_out%rootr_pagl(maxPatchesPerSite,hlm_numlevgrnd))
      allocate(bc_out%btran_pa(maxPatchesPerSite))
      
      ! Photosynthesis

      allocate(bc_out%rssun_pa(maxPatchesPerSite))
      allocate(bc_out%rssha_pa(maxPatchesPerSite))
      
      ! Canopy Radiation
      allocate(bc_out%albd_parb(maxPatchesPerSite,hlm_numSWb))
      allocate(bc_out%albi_parb(maxPatchesPerSite,hlm_numSWb))
      allocate(bc_out%fabd_parb(maxPatchesPerSite,hlm_numSWb))
      allocate(bc_out%fabi_parb(maxPatchesPerSite,hlm_numSWb))
      allocate(bc_out%ftdd_parb(maxPatchesPerSite,hlm_numSWb))
      allocate(bc_out%ftid_parb(maxPatchesPerSite,hlm_numSWb))
      allocate(bc_out%ftii_parb(maxPatchesPerSite,hlm_numSWb))

      ! biogeochemistry
      allocate(bc_out%FATES_c_to_litr_lab_c_col(hlm_numlevdecomp_full))        
      allocate(bc_out%FATES_c_to_litr_cel_c_col(hlm_numlevdecomp_full))
      allocate(bc_out%FATES_c_to_litr_lig_c_col(hlm_numlevdecomp_full))

      ! Canopy Structure
      allocate(bc_out%elai_pa(maxPatchesPerSite))
      allocate(bc_out%esai_pa(maxPatchesPerSite))
      allocate(bc_out%tlai_pa(maxPatchesPerSite))
      allocate(bc_out%tsai_pa(maxPatchesPerSite))
      allocate(bc_out%htop_pa(maxPatchesPerSite))
      allocate(bc_out%hbot_pa(maxPatchesPerSite))
      allocate(bc_out%dleaf_pa(maxPatchesPerSite))

      allocate(bc_out%displa_pa(maxPatchesPerSite))
      allocate(bc_out%z0m_pa(maxPatchesPerSite))

      allocate(bc_out%canopy_fraction_pa(maxPatchesPerSite))
      allocate(bc_out%frac_veg_nosno_alb_pa(maxPatchesPerSite))

      ! Plant-Hydro BC's
      if (hlm_use_planthydro.eq.itrue) then
         allocate(bc_out%qflx_soil2root_sisl(hlm_numlevsoil))
      end if

      return
   end subroutine allocate_bcout

   ! ====================================================================================

   subroutine zero_bcs(this,s)

      implicit none
      class(fates_interface_type), intent(inout) :: this
      integer, intent(in) :: s

      ! Input boundaries
      this%bc_in(s)%zi_sisl(:)     = 0.0_r8
      this%bc_in(s)%dz_sisl(:)     = 0.0_r8
      this%bc_in(s)%z_sisl(:)      = 0.0_r8
      this%bc_in(s)%dz_decomp_sisl = 0.0_r8
      
      this%bc_in(s)%t_veg24_si     = 0.0_r8
      this%bc_in(s)%t_veg24_pa(:)  = 0.0_r8
      this%bc_in(s)%precip24_pa(:) = 0.0_r8
      this%bc_in(s)%relhumid24_pa(:) = 0.0_r8
      this%bc_in(s)%wind24_pa(:)     = 0.0_r8

      this%bc_in(s)%solad_parb(:,:)     = 0.0_r8
      this%bc_in(s)%solai_parb(:,:)     = 0.0_r8
      this%bc_in(s)%smp_gl(:)           = 0.0_r8
      this%bc_in(s)%eff_porosity_gl(:)  = 0.0_r8
      this%bc_in(s)%watsat_gl(:)        = 0.0_r8
      this%bc_in(s)%tempk_gl(:)         = 0.0_r8
      this%bc_in(s)%h2o_liqvol_gl(:)    = 0.0_r8
      this%bc_in(s)%filter_vegzen_pa(:) = .false.
      this%bc_in(s)%coszen_pa(:)        = 0.0_r8
      this%bc_in(s)%albgr_dir_rb(:)     = 0.0_r8
      this%bc_in(s)%albgr_dif_rb(:)     = 0.0_r8
      this%bc_in(s)%max_rooting_depth_index_col = 0
      this%bc_in(s)%tot_het_resp        = 0.0_r8
      this%bc_in(s)%tot_somc            = 0.0_r8 
      this%bc_in(s)%tot_litc            = 0.0_r8
      this%bc_in(s)%snow_depth_si       = 0.0_r8
      this%bc_in(s)%frac_sno_eff_si     = 0.0_r8

      if (hlm_use_planthydro.eq.itrue) then
  
         this%bc_in(s)%qflx_transp_pa(:) = 0.0_r8
         this%bc_in(s)%swrad_net_pa(:) = 0.0_r8
         this%bc_in(s)%lwrad_net_pa(:) = 0.0_r8
         this%bc_in(s)%watsat_sisl(:) = 0.0_r8
         this%bc_in(s)%watres_sisl(:) = 0.0_r8
         this%bc_in(s)%sucsat_sisl(:) = 0.0_r8
         this%bc_in(s)%bsw_sisl(:) = 0.0_r8
         this%bc_in(s)%hksat_sisl(:) = 0.0_r8
      end if


      ! Output boundaries
      this%bc_out(s)%active_suction_gl(:) = .false.
      this%bc_out(s)%fsun_pa(:)      = 0.0_r8
      this%bc_out(s)%laisun_pa(:)    = 0.0_r8
      this%bc_out(s)%laisha_pa(:)    = 0.0_r8
      this%bc_out(s)%rootr_pagl(:,:) = 0.0_r8
      this%bc_out(s)%btran_pa(:)     = 0.0_r8

      this%bc_out(s)%FATES_c_to_litr_lab_c_col(:) = 0.0_r8
      this%bc_out(s)%FATES_c_to_litr_cel_c_col(:) = 0.0_r8
      this%bc_out(s)%FATES_c_to_litr_lig_c_col(:) = 0.0_r8

      this%bc_out(s)%rssun_pa(:)     = 0.0_r8
      this%bc_out(s)%rssha_pa(:)     = 0.0_r8

      this%bc_out(s)%albd_parb(:,:) = 0.0_r8
      this%bc_out(s)%albi_parb(:,:) = 0.0_r8
      this%bc_out(s)%fabd_parb(:,:) = 0.0_r8
      this%bc_out(s)%fabi_parb(:,:) = 0.0_r8
      this%bc_out(s)%ftdd_parb(:,:) = 0.0_r8
      this%bc_out(s)%ftid_parb(:,:) = 0.0_r8
      this%bc_out(s)%ftii_parb(:,:) = 0.0_r8

      this%bc_out(s)%elai_pa(:)   = 0.0_r8
      this%bc_out(s)%esai_pa(:)   = 0.0_r8
      this%bc_out(s)%tlai_pa(:)   = 0.0_r8
      this%bc_out(s)%tsai_pa(:)   = 0.0_r8
      this%bc_out(s)%htop_pa(:)   = 0.0_r8
      this%bc_out(s)%hbot_pa(:)   = 0.0_r8
      this%bc_out(s)%displa_pa(:) = 0.0_r8
      this%bc_out(s)%z0m_pa(:)    = 0.0_r8
      this%bc_out(s)%dleaf_pa(:)   = 0.0_r8

      this%bc_out(s)%canopy_fraction_pa(:) = 0.0_r8
      this%bc_out(s)%frac_veg_nosno_alb_pa(:) = 0.0_r8

      if (hlm_use_planthydro.eq.itrue) then
         this%bc_out(s)%qflx_soil2root_sisl(:) = 0.0_r8
      end if
      this%bc_out(s)%plant_stored_h2o_si = 0.0_r8

      return
   end subroutine zero_bcs


    ! ===================================================================================
    
    subroutine set_fates_global_elements(use_fates)

       ! --------------------------------------------------------------------------------
       !
       ! This subroutine is called directly from the HLM, and is the first FATES routine
       ! that is called.
       !
       ! This subroutine MUST BE CALLED AFTER the FATES PFT parameter file has been read in,
       ! and the EDPftvarcon_inst structure has been made.
       ! This subroutine must ALSO BE CALLED BEFORE the history file dimensions
       ! are set.
       ! 
       ! This routine requires no information from the HLM. This routine is responsible
       ! for generating the globals that are required by the HLM that are entirely
       ! FATES derived.
       !
       ! --------------------------------------------------------------------------------

      use EDParamsMod, only : ED_val_history_sizeclass_bin_edges, ED_val_history_ageclass_bin_edges
      use CLMFatesParamInterfaceMod         , only : FatesReadParameters
      implicit none
      
      logical,intent(in) :: use_fates    ! Is fates turned on?
      
      integer :: i
      
      if (use_fates) then

         ! first read the non-PFT parameters
         call FatesReadParameters()

         ! Identify the number of PFTs by evaluating a pft array
         ! Using wood density as that is not expected to be deprecated any time soon

         if(lbound(EDPftvarcon_inst%wood_density(:),dim=1) .eq. 0 ) then
            numpft = size(EDPftvarcon_inst%wood_density,dim=1)-1
         elseif(lbound(EDPftvarcon_inst%wood_density(:),dim=1) .eq. 1 ) then
            numpft = size(EDPftvarcon_inst%wood_density,dim=1)
         else
            write(fates_log(), *) 'While assessing the number of FATES PFTs,'
            write(fates_log(), *) 'it was found that the lower bound was neither 0 or 1?'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(numpft>maxpft) then
            write(fates_log(), *) 'The number of PFTs dictated by the FATES parameter file'
            write(fates_log(), *) 'is larger than the maximum allowed. Increase the FATES parameter constant'
            write(fates_log(), *) 'FatesInterfaceMod.F90:maxpft accordingly'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         

         ! These values are used to define the restart file allocations and general structure
         ! of memory for the cohort arrays
         
         fates_maxElementsPerPatch = max(maxCohortsPerPatch, &
               numpft * nclmax * nlevleaf)
         
         fates_maxElementsPerSite = maxPatchesPerSite * fates_maxElementsPerPatch

         ! Identify number of size and age class bins for history output
         ! assume these arrays are 1-indexed
         nlevsclass = size(ED_val_history_sizeclass_bin_edges,dim=1)
         nlevage = size(ED_val_history_ageclass_bin_edges,dim=1)

         ! do some checks on the size and age bin arrays to make sure they make sense:
         ! make sure that both start at zero, and that both are monotonically increasing
         if ( ED_val_history_sizeclass_bin_edges(1) .ne. 0._r8 ) then
            write(fates_log(), *) 'size class bins specified in parameter file must start at zero'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         endif
         if ( ED_val_history_ageclass_bin_edges(1) .ne. 0._r8 ) then
            write(fates_log(), *) 'age class bins specified in parameter file must start at zero'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         endif
         do i = 2,nlevsclass
            if ( (ED_val_history_sizeclass_bin_edges(i) - ED_val_history_sizeclass_bin_edges(i-1)) .le. 0._r8) then
               write(fates_log(), *) 'size class bins specified in parameter file must be monotonically increasing'
               call endrun(msg=errMsg(sourcefile, __LINE__))
            end if
         end do
         do i = 2,nlevage
            if ( (ED_val_history_ageclass_bin_edges(i) - ED_val_history_ageclass_bin_edges(i-1)) .le. 0._r8) then
               write(fates_log(), *) 'age class bins specified in parameter file must be monotonically increasing'
               call endrun(msg=errMsg(sourcefile, __LINE__))
            end if
         end do

         ! Set Various Mapping Arrays used in history output as well
         ! These will not be used if use_ed or use_fates is false
         call fates_history_maps()


      else
         ! If we are not using FATES, the cohort dimension is still
         ! going to be initialized, lets set it to the smallest value
         ! possible so that the dimensioning info takes up little space

         fates_maxElementsPerPatch = 1
      
         fates_maxElementsPerSite = 1
         

      end if


    end subroutine set_fates_global_elements

    !==============================================================================================
    
    subroutine fates_history_maps
       
       use EDTypesMod, only : NFSC
       use EDTypesMod, only : NCWD
       use EDTypesMod, only : nclmax
       use EDTypesMod, only : nlevleaf
       use EDParamsMod, only : ED_val_history_sizeclass_bin_edges
       use EDParamsMod, only : ED_val_history_ageclass_bin_edges

       ! ------------------------------------------------------------------------------------------
       ! This subroutine allocates and populates the variables
       ! that define the mapping of variables in history files in multiplexed dimensions liked
       ! the "scpf" format
       ! back to
       ! their respective single component dimensions, like size-class "sc" and pft "pf"
       ! ------------------------------------------------------------------------------------------

       integer :: i
       integer :: isc
       integer :: ipft
       integer :: icwd
       integer :: ifuel
       integer :: ican
       integer :: ileaf
       integer :: iage

       allocate( fates_hdim_levsclass(1:nlevsclass   ))
       allocate( fates_hdim_pfmap_levscpf(1:nlevsclass*numpft))
       allocate( fates_hdim_scmap_levscpf(1:nlevsclass*numpft))
       allocate( fates_hdim_levpft(1:numpft   ))
       allocate( fates_hdim_levfuel(1:NFSC   ))
       allocate( fates_hdim_levcwdsc(1:NCWD   ))
       allocate( fates_hdim_levage(1:nlevage   ))

       allocate( fates_hdim_levcan(nclmax))
       allocate( fates_hdim_canmap_levcnlf(nlevleaf*nclmax))
       allocate( fates_hdim_lfmap_levcnlf(nlevleaf*nclmax))
       allocate( fates_hdim_canmap_levcnlfpf(nlevleaf*nclmax*numpft))
       allocate( fates_hdim_lfmap_levcnlfpf(nlevleaf*nclmax*numpft))
       allocate( fates_hdim_pftmap_levcnlfpf(nlevleaf*nclmax*numpft))
       allocate( fates_hdim_scmap_levscag(nlevsclass * nlevage ))
       allocate( fates_hdim_agmap_levscag(nlevsclass * nlevage ))

       ! Fill the IO array of plant size classes
       fates_hdim_levsclass(:) = ED_val_history_sizeclass_bin_edges(:)
       fates_hdim_levage(:) = ED_val_history_ageclass_bin_edges(:)

       ! make pft array
       do ipft=1,numpft
          fates_hdim_levpft(ipft) = ipft
       end do

       ! make fuel array
       do ifuel=1,NFSC
          fates_hdim_levfuel(ifuel) = ifuel
       end do

       ! make cwd array
       do icwd=1,NCWD
          fates_hdim_levcwdsc(icwd) = icwd
       end do

       ! make canopy array
       do ican = 1,nclmax
          fates_hdim_levcan(ican) = ican
       end do

       ! Fill the IO arrays that match pft and size class to their combined array
       i=0
       do ipft=1,numpft
          do isc=1,nlevsclass
             i=i+1
             fates_hdim_pfmap_levscpf(i) = ipft
             fates_hdim_scmap_levscpf(i) = isc
          end do
       end do

       i=0
       do ican=1,nclmax
          do ileaf=1,nlevleaf
             i=i+1
             fates_hdim_canmap_levcnlf(i) = ican
             fates_hdim_lfmap_levcnlf(i) = ileaf
          end do
       end do

       i=0
       do iage=1,nlevage
          do isc=1,nlevsclass
             i=i+1
             fates_hdim_scmap_levscag(i) = isc
             fates_hdim_agmap_levscag(i) = iage
          end do
       end do

       i=0
       do ipft=1,numpft
          do ican=1,nclmax
             do ileaf=1,nlevleaf
                i=i+1
                fates_hdim_canmap_levcnlfpf(i) = ican
                fates_hdim_lfmap_levcnlfpf(i) = ileaf
                fates_hdim_pftmap_levcnlfpf(i) = ipft
             end do
          end do
       end do

    end subroutine fates_history_maps

    ! ===================================================================================

    subroutine SetFatesTime(current_year_in, current_month_in, &
                          current_day_in, current_tod_in, &
                          current_date_in, reference_date_in, &
                          model_day_in, day_of_year_in, &
                          days_per_year_in, freq_day_in)

     ! This subroutine should be called directly from the HLM
     
     integer,  intent(in) :: current_year_in
     integer,  intent(in) :: current_month_in
     integer,  intent(in) :: current_day_in
     integer,  intent(in) :: current_tod_in
     integer,  intent(in) :: current_date_in
     integer,  intent(in) :: reference_date_in
     real(r8), intent(in) :: model_day_in
     integer,  intent(in) :: day_of_year_in
     integer,  intent(in) :: days_per_year_in
     real(r8), intent(in) :: freq_day_in

     hlm_current_year   = current_year_in
     hlm_current_month  = current_month_in
     hlm_current_day    = current_day_in
     hlm_current_tod    = current_tod_in
     hlm_current_date   = current_date_in
     hlm_reference_date = reference_date_in
     hlm_model_day      = model_day_in
     hlm_day_of_year    = day_of_year_in
     hlm_days_per_year  = days_per_year_in
     hlm_freq_day       = freq_day_in

  end subroutine SetFatesTime

  ! ==================================================================================== 

  subroutine set_fates_ctrlparms(tag,ival,rval,cval)
      
      ! ---------------------------------------------------------------------------------
      ! Certain model control parameters and dimensions used by FATES are dictated by 
      ! the the driver or the host mode. To see which parameters should be filled here
      ! please also look at the ctrl_parms_type in FATESTYpeMod, in the section listing
      ! components dictated by the host model.
      !
      ! Some important points:
      ! 1. Calls to this function are likely from the clm_fates module in the HLM.
      ! 2. The calls should be preceeded by a flush function.
      ! 3. All values in ctrl_parm (FATESTypesMod.F90) that are classified as 
      !    'dictated by the HLM' must be listed in this subroutine
      ! 4. Should look like this:
      ! 
      ! call set_fates_ctrlparms('flush_to_unset')
      ! call set_fates_ctrlparms('num_sw_bbands',numrad)  ! or other variable
      ! ...
      ! call set_fates_ctrlparms('num_lev_ground',nlevgrnd)   ! or other variable
      ! call set_fates_ctrlparms('check_allset') 
      !
      ! RGK-2016
      ! ---------------------------------------------------------------------------------

      ! Arguments
      integer, optional, intent(in)         :: ival
      real(r8), optional, intent(in)        :: rval
      character(len=*),optional, intent(in) :: cval
      character(len=*),intent(in)           :: tag
      
      ! local variables
      logical              :: all_set
      integer,  parameter  :: unset_int = -999
      real(r8), parameter  :: unset_double = -999.9
      
      
      select case (trim(tag))
      case('flush_to_unset')
         if (fates_global_verbose()) then
            write(fates_log(), *) 'Flushing FATES control parameters prior to transfer from host'
         end if

         hlm_numSWb     = unset_int
         hlm_inir       = unset_int
         hlm_ivis       = unset_int
         hlm_is_restart = unset_int
         hlm_numlevgrnd = unset_int
         hlm_numlevsoil = unset_int
         hlm_numlevdecomp_full = unset_int
         hlm_numlevdecomp = unset_int
         hlm_name         = 'unset'
         hlm_hio_ignore_val   = unset_double
         hlm_masterproc   = unset_int
         hlm_ipedof       = unset_int
         hlm_max_patch_per_site = unset_int
         hlm_use_vertsoilc = unset_int
         hlm_use_spitfire  = unset_int
         hlm_use_planthydro = unset_int
         hlm_use_logging   = unset_int
         hlm_use_ed_st3    = unset_int
         hlm_use_ed_prescribed_phys = unset_int
         hlm_use_inventory_init = unset_int
         hlm_inventory_ctrl_file = 'unset'

      case('check_allset')
         
         if(hlm_numSWb .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES dimension/parameter unset: num_sw_rad_bbands'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_masterproc .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES parameter unset: hlm_masterproc'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_numSWb > maxSWb) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES sets a maximum number of shortwave bands'
               write(fates_log(), *) 'for some scratch-space, maxSWb'
               write(fates_log(), *) 'it defaults to 2, but can be increased as needed'
               write(fates_log(), *) 'your driver or host model is intending to drive'
               write(fates_log(), *) 'FATES with:',hlm_numSWb,' bands.'
               write(fates_log(), *) 'please increase maxSWb in EDTypes to match'
               write(fates_log(), *) 'or exceed this value'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         
         if (  .not.((hlm_use_planthydro.eq.1).or.(hlm_use_planthydro.eq.0))    ) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'The FATES namelist planthydro flag must be 0 or 1, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if ( .not.((hlm_use_logging .eq.1).or.(hlm_use_logging.eq.0))    ) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'The FATES namelist use_logging flag must be 0 or 1, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if (  .not.((hlm_use_ed_st3.eq.1).or.(hlm_use_ed_st3.eq.0))    ) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'The FATES namelist stand structure flag must be 0 or 1, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if (  .not.((hlm_use_ed_prescribed_phys.eq.1).or.(hlm_use_ed_prescribed_phys.eq.0))    ) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'The FATES namelist prescribed physiology flag must be 0 or 1, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if ( hlm_use_ed_prescribed_phys.eq.1 .and. hlm_use_ed_st3.eq.1 ) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES ST3 and prescribed physiology cannot both be turned on.'
               write(fates_log(), *) 'Review the namelist entries, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         
         if (  .not.((hlm_use_inventory_init.eq.1).or.(hlm_use_inventory_init.eq.0))    ) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'The FATES NL inventory flag must be 0 or 1, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         
         if(trim(hlm_inventory_ctrl_file) .eq. 'unset') then
            if (fates_global_verbose()) then
               write(fates_log(),*) 'namelist entry for fates inventory control file is unset, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_ivis .ne. ivis) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES assumption about the index of visible shortwave'
               write(fates_log(), *) 'radiation is different from the HLM, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         
         if(hlm_inir .ne. inir) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES assumption about the index of NIR shortwave'
               write(fates_log(), *) 'radiation is different from the HLM, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_is_restart .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES parameter unset: hlm_is_restart, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_numlevgrnd .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES dimension/parameter unset: numlevground, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_numlevsoil .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES dimension/parameter unset: numlevground, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_numlevdecomp_full .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES dimension/parameter unset: numlevdecomp_full, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_numlevdecomp .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES dimension/parameter unset: numlevdecomp, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(trim(hlm_name) .eq. 'unset') then
            if (fates_global_verbose()) then
               write(fates_log(),*) 'FATES dimension/parameter unset: hlm_name, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if( abs(hlm_hio_ignore_val-unset_double)<1e-10 ) then
            if (fates_global_verbose()) then
               write(fates_log(),*) 'FATES dimension/parameter unset: hio_ignore'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_ipedof .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'index for the HLMs pedotransfer function unset: hlm_ipedof, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_max_patch_per_site .eq. unset_int ) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'the number of patch-space per site unset: hlm_max_patch_per_site, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         elseif(hlm_max_patch_per_site < maxPatchesPerSite ) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'FATES is trying to allocate space for more patches per site, than the HLM has space for.'
               write(fates_log(), *) 'hlm_max_patch_per_site (HLM side): ', hlm_max_patch_per_site
               write(fates_log(), *) 'maxPatchesPerSite (FATES side): ', maxPatchesPerSite
               write(fates_log(), *)
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_use_vertsoilc .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'switch for the HLMs soil carbon discretization unset: hlm_use_vertsoilc, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if(hlm_use_spitfire .eq. unset_int) then
            if (fates_global_verbose()) then
               write(fates_log(), *) 'switch for SPITFIRE unset: hlm_use_spitfire, exiting'
            end if
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         if (fates_global_verbose()) then
            write(fates_log(), *) 'Checked. All control parameters sent to FATES.'
         end if

         
      case default

         if(present(ival))then
            select case (trim(tag))

            case('masterproc')
               hlm_masterproc = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering masterproc = ',ival,' to FATES'
               end if

            case('num_sw_bbands')
               hlm_numSwb = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering num_sw_bbands = ',ival,' to FATES'
               end if
               
            case('vis_sw_index')
               hlm_ivis = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering index associated with visible SW rad = ',ival,' to FATES'
               end if
            
            case('nir_sw_index')
               hlm_inir = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering index associated with NIR SW rad = ',ival,' to FATES'
               end if

            case('is_restart')
               hlm_is_restart = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering flag signaling restart / not-restart = ',ival,' to FATES'
               end if

            case('num_lev_ground')
               hlm_numlevgrnd = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering num_lev_ground = ',ival,' to FATES'
               end if

            case('num_lev_soil')
               hlm_numlevsoil = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering num_lev_ground = ',ival,' to FATES'
               end if

            case('num_levdecomp_full')
               hlm_numlevdecomp_full = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering num_levdecomp_full = ',ival,' to FATES'
               end if
            
            case('num_levdecomp')
               hlm_numlevdecomp = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering num_levdecomp = ',ival,' to FATES'
               end if

            case('soilwater_ipedof')
               hlm_ipedof = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_ipedof = ',ival,' to FATES'
               end if

            case('max_patch_per_site')
               hlm_max_patch_per_site = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_max_patch_per_site = ',ival,' to FATES'
               end if

            case('use_vertsoilc')
               hlm_use_vertsoilc = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_use_vertsoilc= ',ival,' to FATES'
               end if

            case('use_spitfire')
               hlm_use_spitfire = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_use_spitfire= ',ival,' to FATES'
               end if
               
            case('use_planthydro')
               hlm_use_planthydro = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_use_planthydro= ',ival,' to FATES'
               end if

            case('use_logging')
               hlm_use_logging = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_use_logging= ',ival,' to FATES'
               end if

            case('use_ed_st3')
               hlm_use_ed_st3 = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_use_ed_st3= ',ival,' to FATES'
               end if

            case('use_ed_prescribed_phys')
               hlm_use_ed_prescribed_phys = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_use_ed_prescribed_phys= ',ival,' to FATES'
               end if

            case('use_inventory_init')
               hlm_use_inventory_init = ival
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hlm_use_inventory_init= ',ival,' to FATES'
               end if

            case default
               if (fates_global_verbose()) then
                  write(fates_log(), *) 'tag not recognized:',trim(tag)
               end if
               ! end_run
            end select
            
         end if
         
         if(present(rval))then
            select case (trim(tag))
            case ('hio_ignore_val')
               hlm_hio_ignore_val = rval
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering hio_ignore_val = ',rval,' to FATES'
               end if
            case default
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'tag not recognized:',trim(tag)
               end if
               ! end_run
            end select
         end if

         if(present(cval))then
            select case (trim(tag))
               
            case('hlm_name')
               hlm_name = trim(cval)
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering the HLM name = ',trim(cval)
               end if

            case('inventory_ctrl_file')
               hlm_inventory_ctrl_file = trim(cval)
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'Transfering the name of the inventory control file = ',trim(cval)
               end if
               
            case default
               if (fates_global_verbose()) then
                  write(fates_log(),*) 'tag not recognized:',trim(tag)
               end if
               ! end_run
            end select
         end if

      end select
            
      return
   end subroutine set_fates_ctrlparms
   
   ! ====================================================================================

   subroutine FatesReportParameters(masterproc)
      
      ! -----------------------------------------------------
      ! Simple parameter reporting functions
      ! A debug like print flag is contained in each routine
      ! -----------------------------------------------------

      logical,intent(in) :: masterproc

      call FatesReportPFTParams(masterproc)
      call FatesReportParams(masterproc)
      
      
      return
   end subroutine FatesReportParameters

end module FatesInterfaceMod

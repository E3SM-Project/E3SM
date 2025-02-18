module pftvarcon

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing vegetation constants and method to
  ! read and initialize vegetation (PFT) constants.
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_log_mod , only : errMsg => shr_log_errMsg
  use abortutils  , only : endrun
  use elm_varpar  , only : mxpft, numrad, ivis, inir
  use elm_varpar  , only:  mxpft_nc
  use elm_varctl  , only : iulog, use_vertsoilc
  use elm_varpar  , only : nlevdecomp_full, nsoilorder
  use elm_varctl  , only : nu_com
  !-------------------------------------------------------------------------------------------
  use elm_varpar  , only : crop_prog
  use elm_varpar  , only : natpft_size, natpft_lb, natpft_ub
  use elm_varpar  , only : cft_size, cft_lb, cft_ub
  use elm_varpar  , only : surfpft_size, surfpft_lb, surfpft_ub
  use elm_varpar  , only : numpft, numcft, maxpatch_pft, max_patch_per_col, maxpatch_urb
  use elm_varctl  , only : create_crop_landunit
  !-------------------------------------------------------------------------------------------
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! Vegetation type constants
  !
  integer :: noveg                  !value for not vegetated
  integer :: ndllf_evr_tmp_tree     !value for Needleleaf evergreen temperate tree
  integer :: ndllf_evr_brl_tree     !value for Needleleaf evergreen boreal tree
  integer :: ndllf_dcd_brl_tree     !value for Needleleaf deciduous boreal tree
  integer :: nbrdlf_evr_trp_tree    !value for Broadleaf evergreen tropical tree
  integer :: nbrdlf_evr_tmp_tree    !value for Broadleaf evergreen temperate tree
  integer :: nbrdlf_dcd_trp_tree    !value for Broadleaf deciduous tropical tree
  integer :: nbrdlf_dcd_tmp_tree    !value for Broadleaf deciduous temperate tree
  integer :: nbrdlf_dcd_brl_tree    !value for Broadleaf deciduous boreal tree
  integer :: ntree                  !value for last type of tree
  integer :: nbrdlf_evr_shrub       !value for Broadleaf evergreen shrub
  integer :: nbrdlf_dcd_tmp_shrub   !value for Broadleaf deciduous temperate shrub
  integer :: nbrdlf_dcd_brl_shrub   !value for Broadleaf deciduous boreal shrub
  integer :: nc3_arctic_grass       !value for C3 arctic grass
  integer :: nc3_nonarctic_grass    !value for C3 non-arctic grass
  integer :: nc4_grass              !value for C4 grass
  integer :: npcropmin              !value for first crop
  integer :: nppercropmin           !value for first perennial crop
  integer :: ncorn                  !value for corn, rain fed (rf)
  integer :: ncornirrig             !value for corn, irrigated (ir)
  integer :: nscereal               !value for spring temperate cereal (rf)
  integer :: nscerealirrig          !value for spring temperate cereal (ir)
  integer :: nwcereal               !value for winter temperate cereal (rf)
  integer :: nwcerealirrig          !value for winter temperate cereal (ir)
  integer :: nsoybean               !value for soybean (rf)
  integer :: nsoybeanirrig          !value for soybean (ir)
  integer :: npcropmax              !value for last prognostic crop in list
  integer :: nppercropmax           !value for last prognostic perennial crop in list
  integer :: nc3crop                !value for generic crop (rf)
  integer :: nc3irrig               !value for irrigated generic crop (ir)
  integer :: ncassava               !value for cassava, rain fed (rf)
  integer :: ncassavairrig          !value for cassava, irrigated (ir)
  integer :: ncotton                !value for cotton, rain fed (rf)
  integer :: ncottonirrig           !value for cotton, irrigated (ir)
  integer :: nfoddergrass           !value for foddergrass, rain fed (rf)
  integer :: nfoddergrassirrig      !value for foddergrass, irrigated (ir)
  integer :: noilpalm               !value for oilpalm, rain fed (rf)
  integer :: noilpalmirrig          !value for oilpalm, irrigated (ir)
  integer :: nograins               !value for other grains, rain fed (rf)
  integer :: nograinsirrig          !value for other grains, irrigated (ir)
  integer :: nrapeseed              !value for rapeseed, rain fed (rf)
  integer :: nrapeseedirrig         !value for rapeseed, irrigated (ir)
  integer :: nrice                  !value for rice, rain fed (rf)
  integer :: nriceirrig             !value for rice, irrigated (ir)
  integer :: nrtubers               !value for root tubers, rain fed (rf)
  integer :: nrtubersirrig          !value for root tubers, irrigated (ir)
  integer :: nsugarcane             !value for sugarcane, rain fed (rf)
  integer :: nsugarcaneirrig        !value for sugarcane, irrigated (ir)
  integer :: nmiscanthus            !value for miscanthus, rain fed (rf)
  integer :: nmiscanthusirrig       !value for miscanthus, irrigated (ir)
  integer :: nswitchgrass           !value for switchgrass, rain fed (rf)
  integer :: nswitchgrassirrig      !value for switchgrass, irrigated (ir)
  integer :: npoplar                !value for poplar, rain fed (rf)
  integer :: npoplarirrig           !value for poplar, irrigated (ir)
  integer :: nwillow                !value for willow, rain fed (rf)
  integer :: nwillowirrig           !value for willow, irrigated (ir)

  ! Number of crop functional types actually used in the model. This includes each CFT for
  ! which is_pft_known_to_model is true. Note that this includes irrigated crops even if
  ! irrigation is turned off in this run: it just excludes crop types that aren't handled
  ! at all, as given by the mergetoelmpft list.
  integer :: num_cfts_known_to_model

  real(r8), allocatable :: dleaf(:)       !characteristic leaf dimension (m)
  real(r8), allocatable :: c3psn(:)       !photosynthetic pathway: 0. = c4, 1. = c3
  real(r8), allocatable :: xl(:)          !leaf/stem orientation index
  real(r8), allocatable :: rhol(:,:)      !leaf reflectance: 1=vis, 2=nir
  real(r8), allocatable :: rhos(:,:)      !stem reflectance: 1=vis, 2=nir
  real(r8), allocatable :: taul(:,:)      !leaf transmittance: 1=vis, 2=nir
  real(r8), allocatable :: taus(:,:)      !stem transmittance: 1=vis, 2=nir
  real(r8), allocatable :: z0mr(:)        !ratio of momentum roughness length to canopy top height (-)
  real(r8), allocatable :: displar(:)     !ratio of displacement height to canopy top height (-)
  real(r8), allocatable :: roota_par(:)   !CLM rooting distribution parameter [1/m]
  real(r8), allocatable :: rootb_par(:)   !CLM rooting distribution parameter [1/m]
  real(r8), allocatable :: crop(:)        !crop pft: 0. = not crop, 1. = crop pft
  real(r8), allocatable :: percrop(:)     !perennial crop pft: 0. = not annual crop, 1. = perennial crop pft
  real(r8), allocatable :: irrigated(:)   !irrigated pft: 0. = not, 1. = irrigated
  real(r8), allocatable :: smpso(:)       !soil water potential at full stomatal opening (mm)
  real(r8), allocatable :: smpsc(:)       !soil water potential at full stomatal closure (mm)
  real(r8), allocatable :: fnitr(:)       !foliage nitrogen limitation factor (-)

  ! begin new pft parameters for CN code
  real(r8), allocatable :: slatop(:)      !SLA at top of canopy [m^2/gC]
  real(r8), allocatable :: dsladlai(:)    !dSLA/dLAI [m^2/gC]
  real(r8), allocatable :: leafcn(:)      !leaf C:N [gC/gN]
  real(r8), allocatable :: flnr(:)        !fraction of leaf N in Rubisco [no units]
  real(r8), allocatable :: woody(:)       !woody lifeform flag (0 = non-woody, 1 = tree, 2 = shrub)
  real(r8), allocatable :: lflitcn(:)     !leaf litter C:N (gC/gN)
  real(r8), allocatable :: frootcn(:)     !fine root C:N (gC/gN)
  real(r8), allocatable :: livewdcn(:)    !live wood (phloem and ray parenchyma) C:N (gC/gN)
  real(r8), allocatable :: deadwdcn(:)    !dead wood (xylem and heartwood) C:N (gC/gN)
  real(r8), allocatable :: grperc(:)      !growth respiration parameter
  real(r8), allocatable :: grpnow(:)      !growth respiration parameter
  real(r8), allocatable :: rootprof_beta(:) !CLM rooting distribution parameter for C and N inputs [unitless]


  ! add pft dependent parameters for phosphorus -X.YANG
  real(r8), allocatable :: leafcp(:)      !leaf C:P [gC/gP]
  real(r8), allocatable :: lflitcp(:)     !leaf litter C:P (gC/gP)
  real(r8), allocatable :: frootcp(:)     !fine root C:P (gC/gP)
  real(r8), allocatable :: livewdcp(:)    !live wood (phloem and ray parenchyma) C:P (gC/gP)
  real(r8), allocatable :: deadwdcp(:)    !dead wood (xylem and heartwood) C:P (gC/gP)

  ! for crop

  ! These arrays give information about the merge of unused crop types to the types CLM
  ! knows about. mergetoelmpft(m) gives the crop type that CLM uses to simulate input
  ! type m (and mergetoelmpft(m) == m implies that CLM simulates crop type m
  ! directly). is_pft_known_to_model(m) is true if CLM simulates crop type m, and false
  ! otherwise. Note that these do NOT relate to whether irrigation is on or off in a
  ! given simulation - that is handled separately.
  integer , allocatable :: mergetoelmpft         (:)
  logical , allocatable :: is_pft_known_to_model (:)

  real(r8), allocatable :: graincn(:)      !grain C:N (gC/gN)
  real(r8), allocatable :: graincp(:)      !grain C:N (gC/gN)
  real(r8), allocatable :: mxtmp(:)        !parameter used in accFlds
  real(r8), allocatable :: baset(:)        !parameter used in accFlds
  real(r8), allocatable :: declfact(:)     !parameter used in CNAllocation
  real(r8), allocatable :: bfact(:)        !parameter used in CNAllocation
  real(r8), allocatable :: aleaff(:)       !parameter used in CNAllocation
  real(r8), allocatable :: arootf(:)       !parameter used in CNAllocation
  real(r8), allocatable :: astemf(:)       !parameter used in CNAllocation
  real(r8), allocatable :: arooti(:)       !parameter used in CNAllocation
  real(r8), allocatable :: fleafi(:)       !parameter used in CNAllocation
  real(r8), allocatable :: allconsl(:)     !parameter used in CNAllocation
  real(r8), allocatable :: allconss(:)     !parameter used in CNAllocation
  real(r8), allocatable :: ztopmx(:)       !parameter used in VegStructUpdate
  real(r8), allocatable :: laimx(:)        !parameter used in VegStructUpdate
  real(r8), allocatable :: gddmin(:)       !parameter used in Phenology
  real(r8), allocatable :: hybgdd(:)       !parameter used in Phenology
  real(r8), allocatable :: lfemerg(:)      !parameter used in Phenology
  real(r8), allocatable :: grnfill(:)      !parameter used in Phenology
  integer , allocatable :: mxmat(:)        !parameter used in Phenology
  integer , allocatable :: mnNHplantdate(:)!minimum planting date for NorthHemisphere (YYYYMMDD)
  integer , allocatable :: mxNHplantdate(:)!maximum planting date for NorthHemisphere (YYYYMMDD)
  integer , allocatable :: mnSHplantdate(:)!minimum planting date for SouthHemisphere (YYYYMMDD)
  integer , allocatable :: mxSHplantdate(:)!maximum planting date for SouthHemisphere (YYYYMMDD)
  real(r8), allocatable :: planttemp(:)    !planting temperature used in Phenology (K)
  real(r8), allocatable :: minplanttemp(:) !mininum planting temperature used in Phenology (K)
  real(r8), allocatable :: senestemp(:)    !senescence temperature for perennial crops used in Phenology (K)
  real(r8), allocatable :: min_days_senes(:)   !minimum leaf age to allow for leaf senescence
  real(r8), allocatable :: froot_leaf(:)   !allocation parameter: new fine root C per new leaf C (gC/gC)
  real(r8), allocatable :: stem_leaf(:)    !allocation parameter: new stem c per new leaf C (gC/gC)
  real(r8), allocatable :: croot_stem(:)   !allocation parameter: new coarse root C per new stem C (gC/gC)
  real(r8), allocatable :: flivewd(:)      !allocation parameter: fraction of new wood that is live (phloem and ray parenchyma) (no units)
  real(r8), allocatable :: fcur(:)         !allocation parameter: fraction of allocation that goes to currently displayed growth, remainder to storage
  real(r8), allocatable :: lf_flab(:)      !leaf litter labile fraction
  real(r8), allocatable :: lf_fcel(:)      !leaf litter cellulose fraction
  real(r8), allocatable :: lf_flig(:)      !leaf litter lignin fraction
  real(r8), allocatable :: fr_flab(:)      !fine root litter labile fraction
  real(r8), allocatable :: fr_fcel(:)      !fine root litter cellulose fraction
  real(r8), allocatable :: fr_flig(:)      !fine root litter lignin fraction
  real(r8), allocatable :: leaf_long(:)    !leaf longevity (yrs)
  real(r8), allocatable :: froot_long(:)   !fine root longevity(yrs)
  real(r8), allocatable :: evergreen(:)    !binary flag for evergreen leaf habit (0 or 1)
  real(r8), allocatable :: stress_decid(:) !binary flag for stress-deciduous leaf habit (0 or 1)
  real(r8), allocatable :: season_decid(:) !binary flag for seasonal-deciduous leaf habit (0 or 1)
  real(r8), allocatable :: pconv(:)        !proportion of deadstem to conversion flux
  real(r8), allocatable :: pprod10(:)      !proportion of deadstem to 10-yr product pool
  real(r8), allocatable :: pprod100(:)     !proportion of deadstem to 100-yr product pool
  real(r8), allocatable :: pprodharv10(:)  !harvest mortality proportion of deadstem to 10-yr pool
  ! pft paraemeters for fire code
  real(r8), allocatable :: cc_leaf(:)
  real(r8), allocatable :: cc_lstem(:)
  real(r8), allocatable :: cc_dstem(:)
  real(r8), allocatable :: cc_other(:)
  real(r8), allocatable :: fm_leaf(:)
  real(r8), allocatable :: fm_lstem(:)
  real(r8), allocatable :: fm_dstem(:)
  real(r8), allocatable :: fm_other(:)
  real(r8), allocatable :: fm_root(:)
  real(r8), allocatable :: fm_lroot(:)
  real(r8), allocatable :: fm_droot(:)
  real(r8), allocatable :: fsr_pft(:)
  real(r8), allocatable :: fd_pft(:)
  ! pft parameters for crop code
  real(r8), allocatable :: manunitro(:)    !fertilizer
  real(r8), allocatable :: fleafcn(:)      !C:N during grain fill; leaf
  real(r8), allocatable :: ffrootcn(:)     !C:N during grain fill; fine root
  real(r8), allocatable :: fstemcn(:)      !C:N during grain fill; stem
  real(r8), allocatable :: presharv(:)     !porportion of residue harvested
  real(r8), allocatable :: convfact(:)     !conversion factor to bu/acre
  real(r8), allocatable :: fyield(:)       !fraction of grain that is actually harvested
  real(r8), allocatable :: root_dmx(:)     !maximum root depth

  integer, parameter :: pftname_len = 40    ! max length of pftname
  character(len=:), allocatable :: pftname(:) !PFT description

  real(r8), parameter :: reinickerp = 1.6_r8 !parameter in allometric equation
  real(r8), parameter :: dwood  = 2.5e5_r8   !cn wood density (gC/m3); lpj:2.0e5
  real(r8), parameter :: allom1 = 100.0_r8   !parameters in
  real(r8), parameter :: allom2 =  40.0_r8   !...allometric
  real(r8), parameter :: allom3 =   0.5_r8   !...equations
  real(r8), parameter :: allom1s = 250.0_r8  !modified for shrubs by
  real(r8), parameter :: allom2s =   8.0_r8  !X.D.Z

  ! Q. Zhu add pft dependent parameters for phosphorus for nutrient competition
  real(r8), allocatable :: VMAX_PLANT_NH4(:)   ! VMAX for plant NH4 uptake
  real(r8), allocatable :: VMAX_PLANT_NO3(:)   ! VMAX for plant NO3 uptake
  real(r8), allocatable :: VMAX_PLANT_P(:)     ! VMAX for plant P uptake
  real(r8), allocatable :: VMAX_MINSURF_P_vr(:,:)! VMAX for P adsorption -> move to soilorder_varcon
  real(r8), allocatable :: KM_PLANT_NH4(:)     ! KM for plant NH4 uptake
  real(r8), allocatable :: KM_PLANT_NO3(:)     ! KM for plant NO3 uptake
  real(r8), allocatable :: KM_PLANT_P(:)       ! KM for plant P uptake
  real(r8), allocatable :: KM_MINSURF_P_vr(:,:)! KM for P adsorption -> move to soilorder_varcon
  real(r8)              :: KM_DECOMP_NH4       ! KM for microbial decomposer NH4 uptake
  real(r8)              :: KM_DECOMP_NO3       ! KM for microbial decomposer NO3 uptake
  real(r8)              :: KM_DECOMP_P         ! KM for microbial decomposer P uptake
  real(r8)              :: KM_NIT              ! KM for nitrifier NH4 uptake
  real(r8)              :: KM_DEN              ! KM for denitrifier NO3 uptake
  real(r8), allocatable :: decompmicc_patch_vr(:,:) ! microbial decomposer biomass gC/m3
  real(r8), allocatable :: alpha_nfix(:)            ! fraction of fixed N goes directly to plant
  real(r8), allocatable :: alpha_ptase(:)           ! fraction of phosphatase produced P goes directly to plant
  real(r8), allocatable :: ccost_nfix(:)            ! plant C cost per unit N produced by N2 fixation
  real(r8), allocatable :: pcost_nfix(:)            ! plant P cost per unit N produced by N2 fixation
  real(r8), allocatable :: ccost_ptase(:)           ! plant C cost per unit P produced by phosphatase
  real(r8), allocatable :: ncost_ptase(:)           ! plant N cost per unit P produced by phosphatase
  real(r8), allocatable :: VMAX_NFIX(:)        ! VMAX of symbiotic N2 fixation
  real(r8), allocatable :: KM_NFIX(:)          ! KM of symbiotic N2 fixation
  real(r8), allocatable :: VMAX_PTASE(:)       ! VMAX of biochemical P production
  real(r8)              :: KM_PTASE            ! KM of biochemical P production
  real(r8)              :: lamda_ptase         ! critical value that incur biochemical production
  real(r8), allocatable :: i_vc(:)             ! intercept of photosynthesis vcmax ~ leaf N content regression model
  real(r8), allocatable :: s_vc(:)             ! slope of photosynthesis vcmax ~ leaf N content regression model
  real(r8), allocatable :: nsc_rtime(:)        ! non-structural carbon residence time
  real(r8), allocatable :: pinit_beta1(:)      ! shaping parameter for P initialization
  real(r8), allocatable :: pinit_beta2(:)      ! shaping parameter for P initialization
  ! new stoichiometry
  real(r8), allocatable :: leafcn_obs(:)       !leaf C:N [gC/gN]
  real(r8), allocatable :: frootcn_obs(:)      !fine root C:N (gC/gN)
  real(r8), allocatable :: livewdcn_obs(:)     !live wood (phloem and ray parenchyma) C:N (gC/gN)
  real(r8), allocatable :: deadwdcn_obs(:)     !dead wood (xylem and heartwood) C:N (gC/gN)
  real(r8), allocatable :: leafcp_obs(:)       !leaf C:P [gC/gP]
  real(r8), allocatable :: frootcp_obs(:)      !fine root C:P (gC/gP)
  real(r8), allocatable :: livewdcp_obs(:)     !live wood (phloem and ray parenchyma) C:P (gC/gP)
  real(r8), allocatable :: deadwdcp_obs(:)     !dead wood (xylem and heartwood) C:P (gC/gP)
  real(r8), allocatable :: leafcn_obs_flex(:,:)       !upper and lower range of leaf C:N [gC/gN]
  real(r8), allocatable :: frootcn_obs_flex(:,:)      !upper and lower range of fine root C:N (gC/gN)
  real(r8), allocatable :: livewdcn_obs_flex(:,:)     !upper and lower range of live wood (phloem and ray parenchyma) C:N (gC/gN)
  real(r8), allocatable :: deadwdcn_obs_flex(:,:)     !upper and lower range of dead wood (xylem and heartwood) C:N (gC/gN)
  real(r8), allocatable :: leafcp_obs_flex(:,:)       !upper and lower range of leaf C:P [gC/gP]
  real(r8), allocatable :: frootcp_obs_flex(:,:)      !upper and lower range of fine root C:P (gC/gP)
  real(r8), allocatable :: livewdcp_obs_flex(:,:)     !upper and lower range of live wood (phloem and ray parenchyma) C:P (gC/gP)
  real(r8), allocatable :: deadwdcp_obs_flex(:,:)     !upper and lower range of dead wood (xylem and heartwood) C:P (gC/gP)
  ! Photosynthesis parameters
  real(r8), allocatable :: fnr(:)              !fraction of nitrogen in RuBisCO
  real(r8), allocatable :: act25(:)
  real(r8), allocatable :: kcha(:)             !Activation energy for kc
  real(r8), allocatable :: koha(:)             !Activation energy for ko
  real(r8), allocatable :: cpha(:)             !Activation energy for cp
  real(r8), allocatable :: vcmaxha(:)          !Activation energy for vcmax
  real(r8), allocatable :: jmaxha(:)           !Activation energy for jmax
  real(r8), allocatable :: tpuha(:)            !Activation energy for tpu
  real(r8), allocatable :: lmrha(:)            !Acitivation energy for lmr
  real(r8), allocatable :: vcmaxhd(:)          !Deactivation energy for vcmax
  real(r8), allocatable :: jmaxhd(:)           !Deactivation energy for jmax
  real(r8), allocatable :: tpuhd(:)            !Deactivation energy for tpu
  real(r8), allocatable :: lmrhd(:)            !Deacitivation energy for lmr
  real(r8), allocatable :: lmrse(:)            !SE for lmr
  real(r8), allocatable :: qe(:)               !Quantum efficiency
  real(r8), allocatable :: theta_cj(:)         !
  real(r8), allocatable :: bbbopt(:)           !Ball-Berry stomatal conductance intercept
  real(r8), allocatable :: mbbopt(:)           !Ball-Berry stomatal conductance slope
  real(r8), allocatable :: nstor(:)            !Nitrogen storage pool timescale
  real(r8), allocatable :: br_xr(:)            !Base rate for excess respiration
  real(r8)              :: tc_stress           !Critial temperature for moisture stress
  real(r8), allocatable :: vcmax_np1(:)        !vcmax~np relationship coefficient
  real(r8), allocatable :: vcmax_np2(:)        !vcmax~np relationship coefficient
  real(r8), allocatable :: vcmax_np3(:)        !vcmax~np relationship coefficient
  real(r8), allocatable :: vcmax_np4(:)        !vcmax~np relationship coefficient
  real(r8)              :: jmax_np1            !jmax~np relationship coefficient
  real(r8)              :: jmax_np2            !jmax~np relationship coefficient
  real(r8)              :: jmax_np3            !jmax~np relationship coefficient
  real(r8)              :: laimax
  ! Hydrology
  real(r8)              :: rsub_top_globalmax
  ! Soil erosion ground cover
  real(r8), allocatable :: gcbc_p(:)           !effectiveness of surface cover in reducing rainfall-driven erosion
  real(r8), allocatable :: gcbc_q(:)           !effectiveness of surface cover in reducing runoff-driven erosion
  real(r8), allocatable :: gcbr_p(:)           !effectiveness of roots in reducing rainfall-driven erosion
  real(r8), allocatable :: gcbr_q(:)           !effectiveness of roots in reducing runoff-driven erosion

  ! NGEE Arctic snow-vegetation interactions
  real(r8), allocatable :: bendresist(:)       ! vegetation resistance to bending under snow loading, 0 to 1 (e.g., Liston and Hiemstra 2011; Sturm et al. 2005)
  real(r8), allocatable :: vegshape(:)         ! shape parameter to modify shrub burial by snow (1 = parabolic, 2 = hemispheric)
  real(r8), allocatable :: stocking(:)         ! stocking density for pft (stems / hectare)
  real(r8), allocatable :: taper(:)            ! ratio of height:radius_breast_height (woody vegetation allometry)
  logical               :: taper_defaults      ! set flag to use taper defaults if not on params file (necessary as import and set values are in different places)

  ! new pft properties, together with woody, crop, percrop, evergreen, stress_decid, season_decid, defined above,
  ! are introduced to define vegetation properties. This will be well defineing a pft so that no indices needed for codes.
  real(r8), allocatable :: climatezone(:)      ! distributed climate zone  (0 = any zone, 1 = tropical, 2 = temperate, 3 = boreal, 4 = arctic)
  real(r8), allocatable :: nonvascular(:)      ! nonvascular lifeform flag (0 = vascular, 1 = moss, 2 = lichen)
  real(r8), allocatable :: needleleaf(:)       ! needleleaf lifeform flag  (0 = broadleaf, 1 = needleleaf)
  real(r8), allocatable :: graminoid(:)        ! graminoid lifeform flag   (0 = nonvascular+woody+crop+percrop, 1 = graminoid)
  logical , allocatable :: iscft(:)            ! crop function type flag   (.false. = non_crop or generic crop, i.e. crop when use_crop=false, .true. = prognostic crop with cft created)
  real(r8), allocatable :: temp_iscft(:)       ! for file read, translated to logical afterwards
  real(r8), allocatable :: nfixer(:)           ! nitrogen fixer flag  (0 = inable, 1 = able to nitrogen fixation from atm. N2)


  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: pftconrd ! Read and initialize vegetation (PFT) constants
  !
  ! !REVISION HISTORY:
  ! Created by Sam Levis (put into module form by Mariana Vertenstein)
  ! 10/21/03, Peter Thornton: Added new variables for CN code
  ! 06/24/09, Erik Kluzek: Add indices for all pft types, and add expected_pftnames array and comparision
  ! 09/17/10, David Lawrence: Modified code to read in netCDF pft physiology file
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine pftconrd
    !
    ! !DESCRIPTION:
    ! Read and initialize vegetation (PFT) constants
    !
    ! !USES:
    use fileutils ,  only : getfil
    use ncdio_pio ,  only : ncd_io, ncd_pio_closefile, ncd_pio_openfile, file_desc_t, &
                            ncd_inqdid, ncd_inqdlen
    use elm_varctl,  only : paramfile, use_fates
    use elm_varctl,  only : use_crop, use_dynroot
    use elm_varcon,  only : tfrz
    use spmdMod   ,  only : masterproc

    !
    ! !ARGUMENTS:
    implicit none
    !
    ! !REVISION HISTORY:
    ! Created by Gordon Bonan
    ! F. Li and S. Levis (11/06/12)
    !
    ! !LOCAL VARIABLES:
    character(len=256) :: locfn ! local file name
    integer :: i,n              ! loop indices
    integer :: ier              ! error code
    type(file_desc_t) :: ncid   ! pio netCDF file id
    integer :: dimid            ! netCDF dimension id
    integer :: npft             ! number of pfts on pft-physiology file
    logical :: readv            ! read variable in or not
    logical :: PFT_DEFAULT      ! pft names are default, i.e. NOT user-defined
    integer :: ncft0, ncft      ! crop pft index of first/last when 'create_crop_landunit' is true
    integer :: noncropmax       ! max non-crop pft index (to check when 'create_crop_landunit' is true)
    real(r8) :: local_iscft      ! a local transfer of iscft from logical to integer for use in error checks
    character(len=32) :: subname = 'pftconrd'              ! subroutine name
    !
    ! Expected PFT names: The names expected on the paramfile file and the order they are expected to be in.
    ! NOTE: similar types are assumed to be together, first trees (ending with broadleaf_deciduous_boreal_tree
    !       then shrubs, ending with broadleaf_deciduous_boreal_shrub, then grasses starting with c3_arctic_grass
    !       and finally crops, ending with soybean
    ! DO NOT CHANGE THE ORDER -- WITHOUT MODIFYING OTHER PARTS OF THE CODE WHERE THE ORDER MATTERS!
    !
    character(len=pftname_len) :: expected_pftnames(0:mxpft)
!-----------------------------------------------------------------------

    expected_pftnames( 0) = 'not_vegetated                      '
    expected_pftnames( 1) = 'needleleaf_evergreen_temperate_tree'
    expected_pftnames( 2) = 'needleleaf_evergreen_boreal_tree   '
    expected_pftnames( 3) = 'needleleaf_deciduous_boreal_tree   '
    expected_pftnames( 4) = 'broadleaf_evergreen_tropical_tree  '
    expected_pftnames( 5) = 'broadleaf_evergreen_temperate_tree '
    expected_pftnames( 6) = 'broadleaf_deciduous_tropical_tree  '
    expected_pftnames( 7) = 'broadleaf_deciduous_temperate_tree '
    expected_pftnames( 8) = 'broadleaf_deciduous_boreal_tree    '
    expected_pftnames( 9) = 'broadleaf_evergreen_shrub          '
    expected_pftnames(10) = 'broadleaf_deciduous_temperate_shrub'
    expected_pftnames(11) = 'broadleaf_deciduous_boreal_shrub   '
    expected_pftnames(12) = 'c3_arctic_grass                    '
    expected_pftnames(13) = 'c3_non-arctic_grass                '
    expected_pftnames(14) = 'c4_grass                           '
    expected_pftnames(15) = 'c3_crop                            '
    expected_pftnames(16) = 'c3_irrigated                       '
    expected_pftnames(17) = 'corn                               '
    expected_pftnames(18) = 'irrigated_corn                     '
    expected_pftnames(19) = 'spring_temperate_cereal            '
    expected_pftnames(20) = 'irrigated_spring_temperate_cereal  '
    expected_pftnames(21) = 'winter_temperate_cereal            '
    expected_pftnames(22) = 'irrigated_winter_temperate_cereal  '
    expected_pftnames(23) = 'soybean                            '
    expected_pftnames(24) = 'irrigated_soybean                  '
    expected_pftnames(25) = 'cassava                            '
    expected_pftnames(26) = 'irrigated_cassava                  '
    expected_pftnames(27) = 'cotton                             '
    expected_pftnames(28) = 'irrigated_cotton                   '
    expected_pftnames(29) = 'foddergrass                        '
    expected_pftnames(30) = 'irrigated_foddergrass              '
    expected_pftnames(31) = 'oilpalm                            '
    expected_pftnames(32) = 'irrigated_oilpalm                  '
    expected_pftnames(33) = 'other_grains                       '
    expected_pftnames(34) = 'irrigated_other_grains             '
    expected_pftnames(35) = 'rapeseed                           '
    expected_pftnames(36) = 'irrigated_rapeseed                 '
    expected_pftnames(37) = 'rice                               '
    expected_pftnames(38) = 'irrigated_rice                     '
    expected_pftnames(39) = 'root_tubers                        '
    expected_pftnames(40) = 'irrigated_root_tubers              '
    expected_pftnames(41) = 'sugarcane                          '
    expected_pftnames(42) = 'irrigated_sugarcane                '
    expected_pftnames(43) = 'miscanthus                         '
    expected_pftnames(44) = 'irrigated_miscanthus               '
    expected_pftnames(45) = 'switchgrass                        '
    expected_pftnames(46) = 'irrigated_switchgrass              '
    expected_pftnames(47) = 'poplar                             '
    expected_pftnames(48) = 'irrigated_poplar                   '
    expected_pftnames(49) = 'willow                             '
    expected_pftnames(50) = 'irrigated_willow                   '

    ! read actual 'npft' from parameter file
    if (masterproc) then
       write(iulog,*) 'Attempting to read PFT physiological data .....'
    end if
    call getfil (paramfile, locfn, 0)
    call ncd_pio_openfile (ncid, trim(locfn), 0)
    call ncd_inqdid(ncid,'pft',dimid)
    call ncd_inqdlen(ncid,dimid,npft)

    ! now 'mxpft' in 'elm_varpar' updated by npft here
    mxpft = npft - 1
    mxpft_nc = npft - 1     ! when .not.use_crop, so here temporarily set and may be changed after read-through PFT-physiology

    ! NOTES:
    !    In parameter file, 'npft' (mxpft) [number of pfts in parameter file] can be greater than 'maxpatch_pft' [number of pfts in surfdata file] which also set by namelist 'maxpft',
    !    but cannot be less, i.e. there cannot be more pfts in the surface file than there are set of pft parameters, if running the non-crop version of model.
    if (npft<maxpatch_pft .and. .not.use_crop) call endrun(msg=' ERROR: pft number less than maxpatch_pft: '//errMsg(__FILE__, __LINE__))

    allocate( dleaf         (0:mxpft) )
    allocate( c3psn         (0:mxpft) )
    allocate( xl            (0:mxpft) )
    allocate( rhol          (0:mxpft,numrad) )
    allocate( rhos          (0:mxpft,numrad) )
    allocate( taul          (0:mxpft,numrad) )
    allocate( taus          (0:mxpft,numrad) )
    allocate( z0mr          (0:mxpft) )
    allocate( displar       (0:mxpft) )
    allocate( roota_par     (0:mxpft) )
    allocate( rootb_par     (0:mxpft) )
    allocate( crop          (0:mxpft) )
    allocate( percrop       (0:mxpft) )
    allocate( irrigated     (0:mxpft) )
    allocate( smpso         (0:mxpft) )
    allocate( smpsc         (0:mxpft) )
    allocate( fnitr         (0:mxpft) )
    allocate( slatop        (0:mxpft) )
    allocate( dsladlai      (0:mxpft) )
    allocate( leafcn        (0:mxpft) )
    allocate( flnr          (0:mxpft) )
    allocate( woody         (0:mxpft) )
    allocate( lflitcn       (0:mxpft) )
    allocate( frootcn       (0:mxpft) )
    allocate( livewdcn      (0:mxpft) )
    allocate( deadwdcn      (0:mxpft) )

    ! add phosphorus
    allocate( leafcp        (0:mxpft) )
    allocate( lflitcp       (0:mxpft) )
    allocate( frootcp       (0:mxpft) )
    allocate( livewdcp      (0:mxpft) )
    allocate( deadwdcp      (0:mxpft) )

    allocate( grperc        (0:mxpft) )
    allocate( grpnow        (0:mxpft) )
    allocate( rootprof_beta (0:mxpft) )

    allocate( mergetoelmpft (0:mxpft) )
    allocate( is_pft_known_to_model  (0:mxpft) )

    allocate( graincn       (0:mxpft) )
    allocate( graincp       (0:mxpft) )
    allocate( mxtmp         (0:mxpft) )
    allocate( baset         (0:mxpft) )
    allocate( declfact      (0:mxpft) )
    allocate( bfact         (0:mxpft) )
    allocate( aleaff        (0:mxpft) )
    allocate( arootf        (0:mxpft) )
    allocate( astemf        (0:mxpft) )
    allocate( arooti        (0:mxpft) )
    allocate( fleafi        (0:mxpft) )
    allocate( allconsl      (0:mxpft) )
    allocate( allconss      (0:mxpft) )
    allocate( ztopmx        (0:mxpft) )
    allocate( laimx         (0:mxpft) )
    allocate( gddmin        (0:mxpft) )
    allocate( hybgdd        (0:mxpft) )
    allocate( lfemerg       (0:mxpft) )
    allocate( grnfill       (0:mxpft) )
    allocate( mxmat         (0:mxpft) )
    allocate( mnNHplantdate (0:mxpft) )
    allocate( mxNHplantdate (0:mxpft) )
    allocate( mnSHplantdate (0:mxpft) )
    allocate( mxSHplantdate (0:mxpft) )
    allocate( planttemp     (0:mxpft) )
    allocate( minplanttemp  (0:mxpft) )
    allocate( senestemp     (0:mxpft) )
    allocate( min_days_senes (0:mxpft) )
    allocate( froot_leaf    (0:mxpft) )
    allocate( stem_leaf     (0:mxpft) )
    allocate( croot_stem    (0:mxpft) )
    allocate( flivewd       (0:mxpft) )
    allocate( fcur          (0:mxpft) )
    allocate( lf_flab       (0:mxpft) )
    allocate( lf_fcel       (0:mxpft) )
    allocate( lf_flig       (0:mxpft) )
    allocate( fr_flab       (0:mxpft) )
    allocate( fr_fcel       (0:mxpft) )
    allocate( fr_flig       (0:mxpft) )
    allocate( leaf_long     (0:mxpft) )
    allocate( froot_long    (0:mxpft) )
    allocate( evergreen     (0:mxpft) )
    allocate( stress_decid  (0:mxpft) )
    allocate( season_decid  (0:mxpft) )
    allocate( pconv         (0:mxpft) )
    allocate( pprod10       (0:mxpft) )
    allocate( pprod100      (0:mxpft) )
    allocate( pprodharv10   (0:mxpft) )
    allocate( cc_leaf       (0:mxpft) )
    allocate( cc_lstem      (0:mxpft) )
    allocate( cc_dstem      (0:mxpft) )
    allocate( cc_other      (0:mxpft) )
    allocate( fm_leaf       (0:mxpft) )
    allocate( fm_lstem      (0:mxpft) )
    allocate( fm_dstem      (0:mxpft) )
    allocate( fm_other      (0:mxpft) )
    allocate( fm_root       (0:mxpft) )
    allocate( fm_lroot      (0:mxpft) )
    allocate( fm_droot      (0:mxpft) )
    allocate( fsr_pft       (0:mxpft) )
    allocate( fd_pft        (0:mxpft) )
    allocate( manunitro     (0:mxpft) )
    allocate( fleafcn       (0:mxpft) )
    allocate( ffrootcn      (0:mxpft) )
    allocate( fstemcn       (0:mxpft) )
    allocate( presharv      (0:mxpft) )
    allocate( convfact      (0:mxpft) )
    allocate( fyield        (0:mxpft) )
    allocate( root_dmx      (0:mxpft) )

    if (use_crop) then
       allocate(character(pftname_len) :: pftname(0:mxpft))
    else
       allocate(character(pftname_len) :: pftname(0:mxpft_nc))
    end if

    allocate( VMAX_PLANT_NH4(0:mxpft) )
    allocate( VMAX_PLANT_NO3(0:mxpft) )
    allocate( VMAX_PLANT_P(0:mxpft) )
    allocate( VMAX_MINSURF_P_vr(1:nlevdecomp_full,0:nsoilorder))
    allocate( KM_PLANT_NH4(0:mxpft) )
    allocate( KM_PLANT_NO3(0:mxpft) )
    allocate( KM_PLANT_P(0:mxpft) )
    allocate( KM_MINSURF_P_vr(1:nlevdecomp_full,0:nsoilorder))
    allocate( decompmicc_patch_vr (1:nlevdecomp_full,0:mxpft))
    allocate( VMAX_PTASE(0:mxpft))
    allocate( i_vc               (0:mxpft) )
    allocate( s_vc               (0:mxpft) )
    allocate( nsc_rtime          (0:mxpft) )
    allocate( pinit_beta1        (0:nsoilorder))
    allocate( pinit_beta2        (0:nsoilorder))
    allocate( alpha_nfix         (0:mxpft) )
    allocate( alpha_ptase        (0:mxpft) )
    allocate( ccost_nfix         (0:mxpft) )
    allocate( pcost_nfix         (0:mxpft) )
    allocate( ccost_ptase        (0:mxpft) )
    allocate( ncost_ptase        (0:mxpft) )
    allocate( VMAX_NFIX          (0:mxpft) )
    allocate( KM_NFIX            (0:mxpft) )
    ! new stoichiometry
    allocate( leafcn_obs         (0:mxpft) )
    allocate( frootcn_obs        (0:mxpft) )
    allocate( livewdcn_obs       (0:mxpft) )
    allocate( deadwdcn_obs       (0:mxpft) )
    allocate( leafcp_obs         (0:mxpft) )
    allocate( frootcp_obs        (0:mxpft) )
    allocate( livewdcp_obs       (0:mxpft) )
    allocate( deadwdcp_obs       (0:mxpft) )
    allocate( leafcn_obs_flex         (0:mxpft,1:2) )
    allocate( frootcn_obs_flex        (0:mxpft,1:2) )
    allocate( livewdcn_obs_flex       (0:mxpft,1:2) )
    allocate( deadwdcn_obs_flex       (0:mxpft,1:2) )
    allocate( leafcp_obs_flex         (0:mxpft,1:2) )
    allocate( frootcp_obs_flex        (0:mxpft,1:2) )
    allocate( livewdcp_obs_flex       (0:mxpft,1:2) )
    allocate( deadwdcp_obs_flex       (0:mxpft,1:2) )
    allocate( vcmax_np1          (0:mxpft) )
    allocate( vcmax_np2          (0:mxpft) )
    allocate( vcmax_np3          (0:mxpft) )
    allocate( vcmax_np4          (0:mxpft) )
    ! Photosynthesis
    allocate( fnr                (0:mxpft) )
    allocate( act25              (0:mxpft) )
    allocate( kcha               (0:mxpft) )
    allocate( koha               (0:mxpft) )
    allocate( cpha               (0:mxpft) )
    allocate( vcmaxha            (0:mxpft) )
    allocate( jmaxha             (0:mxpft) )
    allocate( tpuha              (0:mxpft) )
    allocate( lmrha              (0:mxpft) )
    allocate( vcmaxhd            (0:mxpft) )
    allocate( jmaxhd             (0:mxpft) )
    allocate( tpuhd              (0:mxpft) )
    allocate( lmrhd              (0:mxpft) )
    allocate( lmrse              (0:mxpft) )
    allocate( qe                 (0:mxpft) )
    allocate( theta_cj           (0:mxpft) )
    allocate( bbbopt             (0:mxpft) )
    allocate( mbbopt             (0:mxpft) )
    allocate( nstor              (0:mxpft) )
    allocate( br_xr              (0:mxpft) )
    ! Ground cover for soil erosion
    allocate( gcbc_p             (0:mxpft) )
    allocate( gcbc_q             (0:mxpft) )
    allocate( gcbr_p             (0:mxpft) )
    allocate( gcbr_q             (0:mxpft) )
    ! new pft properties
    allocate( climatezone        (0:mxpft) )
    allocate( nonvascular        (0:mxpft) )
    allocate( graminoid          (0:mxpft) )
    allocate( iscft              (0:mxpft) )
    allocate( temp_iscft         (0:mxpft) )
    allocate( needleleaf         (0:mxpft) )
    allocate( nfixer             (0:mxpft) )

   ! NGEE arctic snow-vegetation interactions
    allocate( bendresist         (0:mxpft) )
    allocate( vegshape           (0:mxpft) )
    allocate( stocking           (0:mxpft) )
    allocate( taper              (0:mxpft) )

    ! Set specific vegetation type values


    call ncd_io('pftname',pftname, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    ! checking if using standard pft-names and default properties
    PFT_DEFAULT = .TRUE.
    do i = 0, min(npft,size(expected_pftnames))-1
       if ( trim(adjustl(pftname(i))) /= trim(expected_pftnames(i)) )then
          PFT_DEFAULT = .FALSE.
       end if
    end do

    call ncd_io('z0mr',z0mr, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('displar',displar, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('dleaf',dleaf, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('c3psn',c3psn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('rholvis',rhol(:,ivis), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('rholnir',rhol(:,inir), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('rhosvis',rhos(:,ivis), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('rhosnir', rhos(:,inir), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('taulvis',taul(:,ivis), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('taulnir',taul(:,inir), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('tausvis',taus(:,ivis), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('tausnir',taus(:,inir), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('xl',xl, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('roota_par',roota_par, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('rootb_par',rootb_par, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('slatop',slatop, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('dsladlai',dsladlai, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('leafcn',leafcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('flnr',flnr, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('smpso',smpso, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('smpsc',smpsc, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('fnitr',fnitr, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('woody',woody, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('lflitcn',lflitcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('frootcn',frootcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('livewdcn',livewdcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('deadwdcn',deadwdcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))


    call ncd_io('leafcp',leafcp, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('lflitcp',lflitcp, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('frootcp',frootcp, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('livewdcp',livewdcp, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('deadwdcp',deadwdcp, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))


    call ncd_io('grperc',grperc, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('grpnow',grpnow, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('froot_leaf',froot_leaf, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('stem_leaf',stem_leaf, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('croot_stem',croot_stem, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('flivewd',flivewd, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('fcur',fcur, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('lf_flab',lf_flab, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('lf_fcel',lf_fcel, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('lf_flig',lf_flig, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('fr_flab',fr_flab, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('fr_fcel',fr_fcel, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('fr_flig',fr_flig, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('leaf_long',leaf_long, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('froot_long',froot_long, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if (.not. readv) froot_long = leaf_long
    !if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('evergreen',evergreen, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('stress_decid',stress_decid, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('season_decid',season_decid, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('fertnitro',manunitro, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('fleafcn',fleafcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('ffrootcn',ffrootcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('fstemcn',fstemcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    if (use_crop) then
       call ncd_io('presharv',presharv, 'read', ncid, readvar=readv, posNOTonfile=.true.)
       if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
       call ncd_io('convfact',convfact, 'read', ncid, readvar=readv, posNOTonfile=.true.)
       if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
       call ncd_io('fyield',fyield, 'read', ncid, readvar=readv, posNOTonfile=.true.)
       if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
       call ncd_io('percrop', percrop, 'read', ncid, readvar=readv)
       if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
       call ncd_io('senescence_temp', senestemp, 'read', ncid, readvar=readv)
       if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
       call ncd_io('min_days_senescence', min_days_senes, 'read', ncid, readvar=readv)
       if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    else
       call ncd_io('percrop', percrop, 'read', ncid, readvar=readv)
       if ( .not. readv ) percrop(:) = 0._r8
    end if
    if(use_dynroot)then
       call ncd_io('root_dmx',root_dmx, 'read', ncid, readvar=readv, posNOTonfile=.true.)
       if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    endif

    if (use_vertsoilc) then
       call ncd_io('rootprof_beta',rootprof_beta, 'read', ncid, readvar=readv, posNOTonfile=.true.)
       if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    end if
    call ncd_io('pconv',pconv, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('pprod10',pprod10, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('pprodharv10',pprodharv10, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('pprod100',pprod100, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('graincn',graincn, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('graincp',graincp, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('mxtmp',mxtmp, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('baset',baset, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('declfact',declfact, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('bfact',bfact, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('aleaff',aleaff, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('arootf',arootf, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('astemf',astemf, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('arooti',arooti, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('fleafi',fleafi, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('allconsl',allconsl, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('allconss',allconss, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('crop',crop, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('irrigated',irrigated, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('ztopmx',ztopmx, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('laimx',laimx, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('gddmin',gddmin, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('hybgdd',hybgdd, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('lfemerg',lfemerg, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('grnfill',grnfill, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('mxmat',mxmat, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('cc_leaf', cc_leaf, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('cc_lstem',cc_lstem, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('cc_dstem',cc_dstem, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('cc_other',cc_other, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('fm_leaf', fm_leaf, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('fm_lstem',fm_lstem, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('fm_dstem',fm_dstem, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('fm_other',fm_other, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('fm_root', fm_root, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('fm_lroot',fm_lroot, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('fm_droot',fm_droot, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('fsr_pft', fsr_pft, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('fd_pft',  fd_pft, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('planting_temp',planttemp, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('min_planting_temp',minplanttemp, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('min_NH_planting_date',mnNHplantdate, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('min_SH_planting_date',mnSHplantdate, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('max_NH_planting_date',mxNHplantdate, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    call ncd_io('max_SH_planting_date',mxSHplantdate, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    if (nu_com .ne. 'RD' ) then

        ! These are soil parameters and used for both FATES and big leaf ELM
        call ncd_io('VMAX_MINSURF_P_vr',VMAX_MINSURF_P_vr, 'read', ncid, readvar=readv)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in soil order VMAX_MINSURF_P_vr'//errMsg(__FILE__, __LINE__))
        call ncd_io('KM_MINSURF_P_vr',KM_MINSURF_P_vr, 'read', ncid, readvar=readv)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in soil order KM_MINSURF_P_vr'//errMsg(__FILE__, __LINE__))
        call ncd_io('KM_DECOMP_NH4',KM_DECOMP_NH4, 'read', ncid, readvar=readv)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in KM_DECOMP_NH4'//errMsg(__FILE__, __LINE__))
        call ncd_io('KM_DECOMP_NO3',KM_DECOMP_NO3, 'read', ncid, readvar=readv)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in KM_DECOMP_NO3'//errMsg(__FILE__, __LINE__))
        call ncd_io('KM_DECOMP_P',KM_DECOMP_P, 'read', ncid, readvar=readv)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in KM_DECOMP_P'//errMsg(__FILE__, __LINE__))
        call ncd_io('KM_NIT',KM_NIT, 'read', ncid, readvar=readv)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in KM_NIT'//errMsg(__FILE__, __LINE__))
        call ncd_io('KM_DEN',KM_DEN, 'read', ncid, readvar=readv)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in KM_DEN'//errMsg(__FILE__, __LINE__))
        call ncd_io('pinit_beta1',pinit_beta1, 'read', ncid, readvar=readv)
        if ( .not. readv ) pinit_beta1(:) = 0.5_r8
        call ncd_io('pinit_beta2',pinit_beta2, 'read', ncid, readvar=readv)
        if ( .not. readv ) pinit_beta2(:) = 0.1_r8

        if(.not.use_fates) then

        call ncd_io('VMAX_PLANT_NH4',VMAX_PLANT_NH4, 'read', ncid, readvar=readv)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft VMAX_PLANT_NH4'//errMsg(__FILE__, __LINE__))
        call ncd_io('VMAX_PLANT_NO3',VMAX_PLANT_NO3, 'read', ncid, readvar=readv)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft VMAX_PLANT_NO3'//errMsg(__FILE__, __LINE__))
        call ncd_io('VMAX_PLANT_P',VMAX_PLANT_P, 'read', ncid, readvar=readv)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft VMAX_PLANT_P'//errMsg(__FILE__, __LINE__))

        call ncd_io('KM_PLANT_NH4',KM_PLANT_NH4, 'read', ncid, readvar=readv)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft KM_PLANT_NH4'//errMsg(__FILE__, __LINE__))
        call ncd_io('KM_PLANT_NO3',KM_PLANT_NO3, 'read', ncid, readvar=readv)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft KM_PLANT_NO3'//errMsg(__FILE__, __LINE__))
        call ncd_io('KM_PLANT_P',KM_PLANT_P, 'read', ncid, readvar=readv)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft KM_PLANT_P'//errMsg(__FILE__, __LINE__))

        call ncd_io('decompmicc_patch_vr',decompmicc_patch_vr, 'read', ncid, readvar=readv)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft decompmicc_patch_vr'//errMsg(__FILE__, __LINE__))
        call ncd_io('alpha_nfix',alpha_nfix, 'read', ncid, readvar=readv)
        if ( .not. readv ) alpha_nfix(:)=0._r8
        call ncd_io('alpha_ptase',alpha_ptase, 'read', ncid, readvar=readv)
        if ( .not. readv ) alpha_ptase(:)=0._r8
        call ncd_io('ccost_nfix',ccost_nfix, 'read', ncid, readvar=readv)
        if ( .not. readv ) ccost_nfix(:)=0._r8
        call ncd_io('pcost_nfix',pcost_nfix, 'read', ncid, readvar=readv)
        if ( .not. readv ) pcost_nfix(:)=0._r8
        call ncd_io('ccost_ptase',ccost_ptase, 'read', ncid, readvar=readv)
        if ( .not. readv ) ccost_ptase(:)=0._r8
        call ncd_io('ncost_ptase',ncost_ptase, 'read', ncid, readvar=readv)
        if ( .not. readv ) ncost_ptase(:)=0._r8
        call ncd_io('VMAX_NFIX',VMAX_NFIX, 'read', ncid, readvar=readv)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in VMAX_NFIX'//errMsg(__FILE__, __LINE__))
        call ncd_io('KM_NFIX',KM_NFIX, 'read', ncid, readvar=readv)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in KM_NFIX'//errMsg(__FILE__, __LINE__))
        call ncd_io('VMAX_PTASE',VMAX_PTASE, 'read', ncid, readvar=readv)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in VMAX_PTASE'//errMsg(__FILE__, __LINE__))
        call ncd_io('KM_PTASE',KM_PTASE, 'read', ncid, readvar=readv)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in KM_PTASE'//errMsg(__FILE__, __LINE__))
        call ncd_io('lamda_ptase',lamda_ptase, 'read', ncid, readvar=readv)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in lamda_ptase'//errMsg(__FILE__, __LINE__))
        call ncd_io('i_vc',i_vc, 'read', ncid, readvar=readv)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in i_vc'//errMsg(__FILE__, __LINE__))
        call ncd_io('s_vc',s_vc, 'read', ncid, readvar=readv)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in s_vc'//errMsg(__FILE__, __LINE__))
        call ncd_io('nsc_rtime',nsc_rtime, 'read', ncid, readvar=readv)
        if ( .not. readv ) nsc_rtime(:) = 1.0_r8
        ! new stoichiometry
        call ncd_io('leafcn_obs',leafcn_obs, 'read', ncid, readvar=readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
        call ncd_io('frootcn_obs',frootcn_obs, 'read', ncid, readvar=readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
        call ncd_io('livewdcn_obs',livewdcn_obs, 'read', ncid, readvar=readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
        call ncd_io('deadwdcn_obs',deadwdcn_obs, 'read', ncid, readvar=readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
        call ncd_io('leafcp_obs',leafcp_obs, 'read', ncid, readvar=readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
        call ncd_io('frootcp_obs',frootcp_obs, 'read', ncid, readvar=readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
        call ncd_io('livewdcp_obs',livewdcp_obs, 'read', ncid, readvar=readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
        call ncd_io('deadwdcp_obs',deadwdcp_obs, 'read', ncid, readvar=readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
         call ncd_io('leafcn_obs_flex',leafcn_obs_flex, 'read', ncid, readvar=readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
        call ncd_io('frootcn_obs_flex',frootcn_obs_flex, 'read', ncid, readvar=readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
        call ncd_io('livewdcn_obs_flex',livewdcn_obs_flex, 'read', ncid, readvar=readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
        call ncd_io('deadwdcn_obs_flex',deadwdcn_obs_flex, 'read', ncid, readvar=readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
        call ncd_io('leafcp_obs_flex',leafcp_obs_flex, 'read', ncid, readvar=readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
        call ncd_io('frootcp_obs_flex',frootcp_obs_flex, 'read', ncid, readvar=readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
        call ncd_io('livewdcp_obs_flex',livewdcp_obs_flex, 'read', ncid, readvar=readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
        call ncd_io('deadwdcp_obs_flex',deadwdcp_obs_flex, 'read', ncid, readvar=readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

        call ncd_io('vcmax_np1',vcmax_np1, 'read', ncid, readvar=readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in vcmax_np data'//errMsg(__FILE__, __LINE__))
        call ncd_io('vcmax_np2',vcmax_np2, 'read', ncid, readvar=readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in vcmax_np data'//errMsg(__FILE__, __LINE__))
        call ncd_io('vcmax_np3',vcmax_np3, 'read', ncid, readvar=readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in vcmax_np data'//errMsg(__FILE__, __LINE__))
        call ncd_io('vcmax_np4',vcmax_np4, 'read', ncid, readvar=readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in vcmax_np data'//errMsg(__FILE__, __LINE__))
        call ncd_io('jmax_np1',jmax_np1, 'read', ncid, readvar=readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in vcmax_np data'//errMsg(__FILE__, __LINE__))
        call ncd_io('jmax_np2',jmax_np2, 'read', ncid, readvar=readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in vcmax_np data'//errMsg(__FILE__, __LINE__))
        call ncd_io('jmax_np3',jmax_np3, 'read', ncid, readvar=readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in vcmax_np data'//errMsg(__FILE__, __LINE__))
        call ncd_io('laimax',laimax, 'read', ncid, readvar=readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=' ERROR: error in reading in laimax data'//errMsg(__FILE__, __LINE__))
    end if
    end if
    call ncd_io('rsub_top_globalmax', rsub_top_globalmax, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if (.not. readv) rsub_top_globalmax = 10._r8
    !if ( .not. readv) call endrun(msg='ERROR:  error in reading in pft data'//errMsg(__FILE__,__LINE__))
    call ncd_io('fnr', fnr, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv) call endrun(msg='ERROR:  error in reading in pft data'//errMsg(__FILE__,__LINE__))
    call ncd_io('act25', act25, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv) call endrun(msg='ERROR:  error in reading in pft data'//errMsg(__FILE__,__LINE__))
    call ncd_io('kcha', kcha, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv) call endrun(msg='ERROR:  error in reading in pft data'//errMsg(__FILE__,__LINE__))
    call ncd_io('koha', koha, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv) call endrun(msg='ERROR:  error in reading in pft data'//errMsg(__FILE__,__LINE__))
    call ncd_io('cpha', cpha, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv) call endrun(msg='ERROR:  error in reading in pft data'//errMsg(__FILE__,__LINE__))
    call ncd_io('vcmaxha', vcmaxha, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv) call endrun(msg='ERROR:  error in reading in pft data'//errMsg(__FILE__,__LINE__))
    call ncd_io('jmaxha', jmaxha, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv) call endrun(msg='ERROR:  error in reading in pft data'//errMsg(__FILE__,__LINE__))
    call ncd_io('tpuha', tpuha, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv) call endrun(msg='ERROR:  error in reading in pft data'//errMsg(__FILE__,__LINE__))
    call ncd_io('lmrha', lmrha, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv) call endrun(msg='ERROR:  error in reading in pft data'//errMsg(__FILE__,__LINE__))
    call ncd_io('vcmaxhd', vcmaxhd, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv) call endrun(msg='ERROR:  error in reading in pft data'//errMsg(__FILE__,__LINE__))
    call ncd_io('jmaxhd', jmaxhd, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv) call endrun(msg='ERROR:  error in reading in pft data'//errMsg(__FILE__,__LINE__))
    call ncd_io('tpuhd', tpuhd, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv) call endrun(msg='ERROR:  error in reading in pftdata'//errMsg(__FILE__,__LINE__))
    call ncd_io('lmrhd', lmrhd, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv) call endrun(msg='ERROR:  error in reading in pft data'//errMsg(__FILE__,__LINE__))
    call ncd_io('lmrse', lmrse, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv) call endrun(msg='ERROR:  error in reading in pft data'//errMsg(__FILE__,__LINE__))
    call ncd_io('qe', qe, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv) call endrun(msg='ERROR:  error in reading in pft data'//errMsg(__FILE__,__LINE__))
    call ncd_io('theta_cj', theta_cj, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv) call endrun(msg='ERROR:  error in reading in pft data'//errMsg(__FILE__,__LINE__))
    call ncd_io('bbbopt', bbbopt, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv) call endrun(msg='ERROR:  error in reading in pft data'//errMsg(__FILE__,__LINE__))
    call ncd_io('mbbopt', mbbopt, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv) call endrun(msg='ERROR:  error in reading in pft data'//errMsg(__FILE__,__LINE__))
    call ncd_io('nstor', nstor, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv) call endrun(msg='ERROR:  error in reading in pft data'//errMsg(__FILE__,__LINE__))
    call ncd_io('br_xr', br_xr, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    !if ( .not. readv) call endrun(msg='ERROR:  error in reading in pft data'//errMsg(__FILE__,__LINE__))
    if (.not. readv) br_xr(:) = 0._r8
    call ncd_io('tc_stress', tc_stress, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv) call endrun(msg='ERROR:  error in reading in pft data'//errMsg(__FILE__,__LINE__))
    call ncd_io('gcbc_p',gcbc_p, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) gcbc_p(:) = 0._r8
    call ncd_io('gcbc_q',gcbc_q, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) gcbc_q(:) = 0._r8
    call ncd_io('gcbr_p',gcbr_p, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) gcbr_p(:) = 0._r8
    call ncd_io('gcbr_q',gcbr_q, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) gcbr_q(:) = 0._r8

    call ncd_io('mergetoclmpft', mergetoelmpft, 'read', ncid, readvar=readv)
    ! in case parameter file is using 'mergetoelmpft'
    if ( .not. readv ) call ncd_io('mergetoelmpft', mergetoelmpft, 'read', ncid, readvar=readv)
    if ( .not. readv ) then
       do i = 0, mxpft
          mergetoelmpft(i) = i
       end do
    end if

    ! NGEE-Arctic snow-vegetation parameters
    call ncd_io('bendresist', bendresist, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if (.not. readv ) bendresist(:) = 1._r8
    call ncd_io('vegshape', vegshape, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if (.not. readv ) vegshape(:) = 1._r8
   ! check validity
    do i = 0, npft-1
      if (bendresist(i) .gt. 1.0_r8 .or. bendresist(i) .le. 0._r8) then
         call endrun(msg="Non-physical selection of bendresist parameter, set between 0 and 1"//errMsg(__FILE__, __LINE__))
      end if
      if (vegshape(i) .gt. 2.0_r8 .or. vegshape(i) .le. 0._r8) then
         call endrun(msg="Non-physical selection of vegshape parameter, set between 0 and 2"//errMsg(__FILE__, __LINE__))
      end if
   end do
    call ncd_io('stocking', stocking, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if (.not. readv ) stocking(:) = 0.1_r8 ! convert previous default of 1000 stems/ha to stems/m2 as had been done in VegStructUpdateMod.F90
   taper_defaults = .false.
    call ncd_io('taper', taper, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if (.not. readv ) then
      taper(:) = 200._r8 ! pftnames not set to integers yet, so reassign further down.
      taper_defaults = .true.
    end if

    ! NOTE: the following 5 PFT flags/options are addtions to 'woody', 'stress_decid', 'season_decid',
    !     'evergreen', and 'crop', 'percop'. For default ELM, will be hard-coded; while for user-defined
    !      it must be included in physiology file.

    ! if 'needleleaf' flag is defined for each PFT
    call ncd_io('needleleaf', needleleaf, 'read', ncid, readvar=readv)
    if ( .not. readv ) then
       if (PFT_DEFAULT) then
          needleleaf(:) = 0   ! will assign a value below
       else
          call endrun(msg='ERROR:  error in reading in user-defined pft data'//errMsg(__FILE__,__LINE__))
       end if
    end if

    ! if 'climatezone' flag is defined for each PFT
    call ncd_io('climatezone', climatezone, 'read', ncid, readvar=readv)
    if ( .not. readv ) then
       if (PFT_DEFAULT) then
          climatezone(:) = 0   ! will assign a value below
       else
          call endrun(msg='ERROR:  error in reading in user-defined pft data'//errMsg(__FILE__,__LINE__))
       end if
    end if

    ! if 'nfixer' flag is defined for each PFT
    call ncd_io('nfixer', nfixer, 'read', ncid, readvar=readv)
    if ( .not. readv ) then
       if (PFT_DEFAULT) then
          nfixer(:) = 0   ! will assign a value below
       else
          call endrun(msg='ERROR:  error in reading in user-defined pft data'//errMsg(__FILE__,__LINE__))
       end if
    end if

    ! if 'nonvascular' flag is defined for each PFT (0: vascular, 1: moss, 2: lichen)
    call ncd_io('nonvascular', nonvascular, 'read', ncid, readvar=readv)
    if ( .not. readv ) then
       if (PFT_DEFAULT) then
          nonvascular(:) = 0   ! will assign a value below
       else
          call endrun(msg='ERROR:  error in reading in user-defined pft data'//errMsg(__FILE__,__LINE__))
       end if
    end if

    call ncd_io('graminoid', graminoid, 'read', ncid, readvar=readv)
    if ( .not. readv ) then
       if (PFT_DEFAULT) then
          graminoid(:) = 0   ! will assign a value below
       else
          call endrun(msg='ERROR:  error in reading in user-defined pft data'//errMsg(__FILE__,__LINE__))
       end if
    end if

    call ncd_io('iscft', temp_iscft, 'read', ncid, readvar=readv)   ! read-in is 'CFT' or not for crop type
    if ( .not. readv ) then
       if (PFT_DEFAULT) then
          temp_iscft(:) = 0._r8   ! will assign a value below
       else
          call endrun(msg='ERROR:  error in reading in user-defined pft data'//errMsg(__FILE__,__LINE__))
       end if
    end if

    call ncd_pio_closefile(ncid)

    ! transfer the temporary real to logical
    do i=0, mxpft
       if (temp_iscft(i) == 1._r8) iscft(i) = .true.
       if (temp_iscft(i) == 0._r8) iscft(i) = .false.
    end do

   if ( PFT_DEFAULT ) then
   ! if still reading in default PFT physiology file,
   ! pft indexing will be as old way

    do i = 0, mxpft

       if(.not. use_crop .and. i > mxpft_nc) EXIT ! exit the do loop

       ! (FATES-INTERF) Later, depending on how the team plans to structure the crop model
       ! or other modules that co-exist while FATES is on, we may want to preserve these pft definitions
       ! on non-fates columns.  For now, they are incompatible, and this check is warranted (rgk 04-2017)

       ! avd - this should be independent of FATES because it fails for non-crop config otherwise
       if(.not. use_fates)then
          if ( trim(adjustl(pftname(i))) /= trim(expected_pftnames(i)) )then
             write(iulog,*)'pftconrd: pftname is NOT what is expected, name = ', &
                  trim(pftname(i)), ', expected name = ', trim(expected_pftnames(i))
             call endrun(msg='pftconrd: bad name for pft on paramfile dataset'//errMsg(__FILE__, __LINE__))
          end if
       end if

       if ( trim(pftname(i)) == 'not_vegetated'                       ) noveg                = i
       if ( trim(pftname(i)) == 'needleleaf_evergreen_temperate_tree' ) ndllf_evr_tmp_tree   = i
       if ( trim(pftname(i)) == 'needleleaf_evergreen_boreal_tree'    ) ndllf_evr_brl_tree   = i
       if ( trim(pftname(i)) == 'needleleaf_deciduous_boreal_tree'    ) ndllf_dcd_brl_tree   = i
       if ( trim(pftname(i)) == 'broadleaf_evergreen_tropical_tree'   ) nbrdlf_evr_trp_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_evergreen_temperate_tree'  ) nbrdlf_evr_tmp_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_deciduous_tropical_tree'   ) nbrdlf_dcd_trp_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_deciduous_temperate_tree'  ) nbrdlf_dcd_tmp_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_deciduous_boreal_tree'     ) nbrdlf_dcd_brl_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_evergreen_shrub'           ) nbrdlf_evr_shrub     = i
       if ( trim(pftname(i)) == 'broadleaf_deciduous_temperate_shrub' ) nbrdlf_dcd_tmp_shrub = i
       if ( trim(pftname(i)) == 'broadleaf_deciduous_boreal_shrub'    ) nbrdlf_dcd_brl_shrub = i
       if ( trim(pftname(i)) == 'c3_arctic_grass'                     ) nc3_arctic_grass     = i
       if ( trim(pftname(i)) == 'c3_non-arctic_grass'                 ) nc3_nonarctic_grass  = i
       if ( trim(pftname(i)) == 'c4_grass'                            ) nc4_grass            = i
       if ( trim(pftname(i)) == 'c3_crop'                             ) nc3crop              = i
       if ( trim(pftname(i)) == 'c3_irrigated'                        ) nc3irrig             = i
       if ( trim(pftname(i)) == 'corn'                                ) ncorn                = i
       if ( trim(pftname(i)) == 'irrigated_corn'                      ) ncornirrig           = i
       if ( trim(pftname(i)) == 'spring_temperate_cereal'             ) nscereal             = i
       if ( trim(pftname(i)) == 'irrigated_spring_temperate_cereal'   ) nscerealirrig        = i
       if ( trim(pftname(i)) == 'winter_temperate_cereal'             ) nwcereal             = i
       if ( trim(pftname(i)) == 'irrigated_winter_temperate_cereal'   ) nwcerealirrig        = i
       if ( trim(pftname(i)) == 'soybean'                             ) nsoybean             = i
       if ( trim(pftname(i)) == 'irrigated_soybean'                   ) nsoybeanirrig        = i
       if ( trim(pftname(i)) == 'cassava'                             ) ncassava             = i
       if ( trim(pftname(i)) == 'irrigated_cassava'                   ) ncassavairrig        = i
       if ( trim(pftname(i)) == 'cotton'                              ) ncotton              = i
       if ( trim(pftname(i)) == 'irrigated_cotton'                    ) ncottonirrig         = i
       if ( trim(pftname(i)) == 'foddergrass'                         ) nfoddergrass         = i
       if ( trim(pftname(i)) == 'irrigated_foddergrass'               ) nfoddergrassirrig    = i
       if ( trim(pftname(i)) == 'oilpalm'                             ) noilpalm             = i
       if ( trim(pftname(i)) == 'irrigated_oilpalm'                   ) noilpalmirrig        = i
       if ( trim(pftname(i)) == 'other_grains'                        ) nograins             = i
       if ( trim(pftname(i)) == 'irrigated_other_grains'              ) nograinsirrig        = i
       if ( trim(pftname(i)) == 'rapeseed'                            ) nrapeseed            = i
       if ( trim(pftname(i)) == 'irrigated_rapeseed'                  ) nrapeseedirrig       = i
       if ( trim(pftname(i)) == 'rice'                                ) nrice                = i
       if ( trim(pftname(i)) == 'irrigated_rice'                      ) nriceirrig           = i
       if ( trim(pftname(i)) == 'root_tubers'                         ) nrtubers             = i
       if ( trim(pftname(i)) == 'irrigated_root_tubers'               ) nrtubersirrig        = i
       if ( trim(pftname(i)) == 'sugarcane'                           ) nsugarcane           = i
       if ( trim(pftname(i)) == 'irrigated_sugarcane'                 ) nsugarcaneirrig      = i
       if ( trim(pftname(i)) == 'miscanthus'                          ) nmiscanthus          = i
       if ( trim(pftname(i)) == 'irrigated_miscanthus'                ) nmiscanthusirrig     = i
       if ( trim(pftname(i)) == 'switchgrass'                         ) nswitchgrass         = i
       if ( trim(pftname(i)) == 'irrigated_switchgrass'               ) nswitchgrassirrig    = i
       if ( trim(pftname(i)) == 'poplar'                              ) npoplar              = i
       if ( trim(pftname(i)) == 'irrigated_poplar'                    ) npoplarirrig         = i
       if ( trim(pftname(i)) == 'willow'                              ) nwillow              = i
       if ( trim(pftname(i)) == 'irrigated_willow'                    ) nwillowirrig         = i
    end do

    ntree                = nbrdlf_dcd_brl_tree  ! value for last type of tree
    npcropmin            = ncorn                ! first prognostic crop
    if( .not. use_crop) then
       npcropmax            = nsoybeanirrig        ! last prognostic crop in list
       nppercropmax         = 0                    ! set value for iscft test below
    else
       npcropmax            = nsugarcaneirrig      ! last prognostic crop in list
       nppercropmin         = nmiscanthus          ! first prognostic perennial crop
       nppercropmax         = nwillowirrig         ! last prognostic perennial crop in list
    end if

   !-------------------------------------------------------------------------------------------
    ! default: for tree and shrub always 1, now hard-coded as following
    woody(ndllf_evr_tmp_tree:nbrdlf_dcd_brl_tree) = 1
    woody(nbrdlf_evr_shrub:nbrdlf_dcd_brl_shrub)  = 2

    ! the following is initialized as 0 above for all PFT.
    ! here the hard-coded values (or flags) for default ELM PFT physiology will be working as original
    ! when not using those indexing of PFT orders anymore in other codes than here.
    needleleaf(noveg+1:ndllf_dcd_brl_tree) = 1
    graminoid(nc3_arctic_grass:nc4_grass) = 1
    iscft(npcropmin:max(npcropmax,nppercropmax)) = .true.
    nfixer(nsoybean)       = 1
    nfixer(nsoybeanirrig)  = 1

    climatezone(ndllf_evr_tmp_tree)  = 2
    climatezone(ndllf_evr_brl_tree)  = 3
    climatezone(ndllf_dcd_brl_tree)  = 3
    climatezone(nbrdlf_evr_trp_tree) = 1
    climatezone(nbrdlf_evr_tmp_tree) = 2
    climatezone(nbrdlf_dcd_trp_tree) = 1
    climatezone(nbrdlf_dcd_tmp_tree) = 2
    climatezone(nbrdlf_dcd_brl_tree) = 3
    climatezone(nbrdlf_dcd_tmp_shrub)= 2
    climatezone(nbrdlf_dcd_brl_shrub)= 3
    climatezone(nc3_arctic_grass)    = 4
    !-------------------------------------------------------------------------------------------


   ! NOT default PFT file
   else

       ! not vegetated checking
       noveg = -1

       ! when user-defined PFTs, if crop included, it must be the default way:
       ! cropts must be in one block and after nat-pft if 'create_crop_landunit' or 'use_crop' is true
       npcropmin            = -1                   ! first prognostic crop
       npcropmax            = -1                   ! last prognostic crop in list
       nppercropmin         = -1                   ! first prognostic perennial crop
       nppercropmax         = -1                   ! last prognostic perennial crop in list

       numcft = 0
       ncft   = -1
       ncft0  = -1
       noncropmax = 0
       do i = 0, npft-1
          if (crop(i)>=1 .or. percrop(i)>=1 .or. iscft(i)) then
              numcft = numcft + 1  ! includes generic_crop, while cft_size NOT (???? todo checking)

              if(use_crop) then

                 ! the following assumes that all crop pfts are in a block

                 ! if 'generic crop' (crop=1) specifically flagged by iscft=.false.
                 ! 'crop' will not be counted into prognostic
                 if (crop(i)>=1 .and. iscft(i)) then
                    npcropmax = i
                    if(npcropmin<=0) npcropmin = i
                 end if

                 ! NOTE: there is a misunderstanding in pft physiology parameter file: 'perennial crop' IS NOT 'crop',
                 !       but still counted into 'numcft'
                 if(percrop(i)==1) then
                    nppercropmax = i
                    if(nppercropmin<=0) nppercropmin = i
                 end if

              else
                 if(crop(i)>=1 .or. percrop(i)>=1) then
                    ! in case either 'crop' or 'generic crop' or both defined, it must be generic, when not use_crop=.true.
                    iscft(i) = .false.
                 end if
              end if

              !
              if (create_crop_landunit) then
                 ! make sure all crop-pft are in one block and following nat-pft SO THAT creating crop landunit is possible
                 ncft = ncft + 1
                 if (ncft0<0) ncft0=i
              end if

          !
          else if (create_crop_landunit) then
              noncropmax = i

          end if

          ! need to check 'noveg'
          if ( woody(i)<=0 .and. graminoid(i)<=0 .and. nonvascular(i)<=0 .and. &
               .not. iscft(i) .and. crop(i)<=0 .and. percrop(i)<=0) then
              if (noveg>=0) then
                 ! not yet support multiple non-vegetated PFT
                 ! this also will catch error of no actual PFT if npft>1
                 call endrun(msg=' ERROR: more than 1 not vegetated in physiology parameter nc file.'//errMsg(__FILE__, __LINE__))
              else
                 noveg = i
              end if
          end if

          !
       end do

       ! make sure non-generic crop indices always beyond natural-pft, even if not available - used in filterMod.F90)
       if (npcropmin < 0 .and. npcropmax < 0) then
          npcropmin = npft
          npcropmax = npft
       end if

       ! MUST re-do some constants which already set in 'elm_varpar.F90:elm_varpar_init()'
       mxpft_nc     = min(maxpatch_pft,npft) - 1      ! user-defined is what max.
       numpft       = min(maxpatch_pft,npft) - 1      ! actual # of patches (without bare)

       if (create_crop_landunit) then
          if (ncft0 /= noncropmax+1) then
             call endrun(msg=' ERROR: when create_crop_landunit is true, crop must be following non-crop PFT .'//errMsg(__FILE__, __LINE__))
          end if
          if (ncft /= npft-1) then
             call endrun(msg=' ERROR: when create_crop_landunit is true, last crop must be the last one of all PFTs .'//errMsg(__FILE__, __LINE__))
          end if

          natpft_size = (numpft + 1) - numcft    ! note that numpft doesn't include bare ground -- thus we add 1
          cft_size    = numcft
       else
          natpft_size = numpft + 1               ! note that numpft doesn't include bare ground -- thus we add 1
          cft_size    = 0
       end if
       natpft_lb = 0
       natpft_ub = natpft_lb + natpft_size - 1
       cft_lb = natpft_ub + 1
       cft_ub = max(cft_lb, cft_lb + cft_size - 1)            ! NOTE: if cft_size is ZERO, could be issue (but so far so good)
       surfpft_lb  = natpft_lb
       surfpft_ub  = natpft_ub
       surfpft_size = natpft_size
       max_patch_per_col= max(numpft+1, numcft, maxpatch_urb)


   end if  ! end if 'PFT_DEFAULT'

     ! checking of pft flags' conflict
     if ( .not. use_fates ) then
        do i = 0, mxpft
          if (iscft(i)) local_iscft = 1._r8
          if (.not. iscft(i)) local_iscft = 0._r8
          if (i == noveg) then
             if ( (nonvascular(i)+woody(i)+graminoid(i)+max(local_iscft,crop(i)+percrop(i))) >= 1 .or. &
                  (needleleaf(i)+evergreen(i)+stress_decid(i)+season_decid(i)+nfixer(i)) >= 1 ) then
                print *, 'ERROR: Incorrect not-vegetated PFT flags: ', i, ' ', trim(pftname(i))
                call endrun(msg=' ERROR: not_vegetated has at least one positive PFT flag '//errMsg(__FILE__, __LINE__))
             end if

          else if ( (nonvascular(i)+woody(i)+graminoid(i)+max(local_iscft,crop(i)+percrop(i))) >= 1) then
             if (nonvascular(i) >= 1 .and. (woody(i)+graminoid(i)+max(local_iscft,crop(i)+percrop(i))) >= 1) then
                print *, 'ERROR: Incorrect nonvasculr PFT flags: ', i, ' ', trim(pftname(i))
                call endrun(msg=' ERROR: nonvascular PFT cannot be any of woody/graminoid/crop type '//errMsg(__FILE__, __LINE__))
             else if (woody(i) >= 1 .and. (nonvascular(i)+graminoid(i)+max(local_iscft,crop(i)+percrop(i))) >= 1) then
                print *, 'ERROR: Incorrect woody PFT flags: ', i, ' ', trim(pftname(i))
                call endrun(msg=' ERROR: woody PFT cannot be any of nonvascular/graminoid/crop type - '//errMsg(__FILE__, __LINE__))
             else if (graminoid(i) >= 1 .and. (nonvascular(i)+woody(i)+max(local_iscft,crop(i)+percrop(i))) >=1 ) then
                print *, 'ERROR: Incorrect graminoid PFT flags: ', i, ' ', trim(pftname(i))
                call endrun(msg=' ERROR: graminoid PFT cannot be any of nonvascular/woody/crop type - '//errMsg(__FILE__, __LINE__))
             else if ( (max(local_iscft,crop(i)+percrop(i))) >= 1 .and. (nonvascular(i)+woody(i)+graminoid(i)) >= 1) then
                print *, 'ERROR: Incorrect crop PFT flags: ', i, ' ', trim(pftname(i))
                call endrun(msg=' ERROR: crop PFT cannot be any of nonvascular/woody/graminoid type - '//errMsg(__FILE__, __LINE__))
             end if

             if( (stress_decid(i)*season_decid(i)) >= 1 ) then
                print *, 'ERROR: Incorrect stress_decid or season_decid flags: ', i, ' ', trim(pftname(i))
                call endrun(msg=' ERROR: stress_decid AND season_decid cannot be both 1 - '//errMsg(__FILE__, __LINE__))
             elseif( (evergreen(i)*(stress_decid(i)+season_decid(i)) ) >= 1 ) then
                print *, 'ERROR: Incorrect evergreen AND season_/stress_decid flags: ', i, ' ', trim(pftname(i))
                call endrun(msg=' ERROR: evergreen AND (stress_decid OR season_decid) cannot be both 1 - '//errMsg(__FILE__, __LINE__))
             end if

          else

             call endrun(msg=' ERROR: not_vegetated AND none of vegetation type for PFT - '//errMsg(__FILE__, __LINE__))

          end if

        end do
     end if

     ! information
     if (masterproc) then
           write(iulog,*)
           write(iulog,*) 'Using PFT physiological parameters from: ', paramfile
           write(iulog,*) '        -- index -- name                                 -- climate zone --    -- woody --     -- needleleaf --    -- evergreen --    -- stress_decid --    -- season_decid --    -- graminoid--    -- iscft --   -- crop --    -- perennial crop --    -- nfixer --'
           do i = 0, npft-1
               write(iulog,*) i, pftname(i), int(climatezone(i)), int(woody(i)), int(needleleaf(i)), &
                int(evergreen(i)), int(stress_decid(i)), int(season_decid(i)), &
                int(graminoid(i)), int(temp_iscft(i)), int(crop(i)), int(percrop(i)), int(nfixer(i))
           end do
           write(iulog,*)
     end if


    !-------------------------------------------------------------------------------------------

    call set_is_pft_known_to_model()
    if (cft_size>0) call set_num_cfts_known_to_model()

    if( .not. use_fates ) then
       if( .not. use_crop) then
          if ( npcropmax /= mxpft_nc .and. crop_prog)then
             call endrun(msg=' ERROR: npcropmax is NOT the last value'//errMsg(__FILE__, __LINE__))
          end if
       else
          if ( nppercropmax /= mxpft )then
             call endrun(msg=' ERROR: nppercropmax is NOT the last value'//errMsg(__FILE__, __LINE__))
          end if
       end if
       do i = 0, mxpft
          if(.not. use_crop .and. i > mxpft_nc) EXIT ! exit the do loop
          if(.not.PFT_DEFAULT) EXIT  ! no checking of indexing PFTs for user-defined
          if( .not. use_crop) then
             if ( irrigated(i) == 1.0_r8  .and. (i == nc3irrig .or. &
                                                 i == ncornirrig .or. &
                                                 i == nscerealirrig .or. &
                                                 i == nwcerealirrig .or. &
                                                 i == nsoybeanirrig ) ) then
                ! correct
             else if ( irrigated(i) == 0.0_r8 )then
                ! correct
             else
                call endrun(msg=' ERROR: irrigated has wrong values'//errMsg(__FILE__, __LINE__))
             end if
          else
             if ( irrigated(i) == 1.0_r8  .and. (i == nc3irrig .or. &
                                                 i == ncornirrig .or. &
                                                 i == nscerealirrig .or. &
                                                 i == nwcerealirrig .or. &
                                                 i == nsoybeanirrig .or. &
                                                 i == ncassavairrig .or. &
                                                 i == ncottonirrig .or. &
                                                 i == nfoddergrassirrig .or. &
                                                 i == noilpalmirrig .or. &
                                                 i == nograinsirrig .or. &
                                                 i == nrapeseedirrig .or. &
                                                 i == nriceirrig .or. &
                                                 i == nrtubersirrig .or. &
                                                 i == nsugarcaneirrig .or. &
                                                 i == nmiscanthusirrig .or. &
                                                 i == nswitchgrassirrig .or. &
                                                 i == npoplarirrig .or. &
                                                 i == nwillowirrig ) )then
                ! correct
             else if ( irrigated(i) == 0.0_r8 )then
                ! correct
             else
                call endrun(msg=' ERROR: irrigated has wrong values'//errMsg(__FILE__, __LINE__))
             end if
             if (      percrop(i) == 1.0_r8 .and. (i >= nmiscanthus .and. i <= nppercropmax))then
                ! correct
             else if ( percrop(i) == 0.0_r8 )then
                ! correct
             else
                call endrun(msg=' ERROR: perennial crop has wrong values'//errMsg(__FILE__, __LINE__))
             end if
          end if
          if (      crop(i) == 1.0_r8 .and. (i >= nc3crop .and. i <= npcropmax) )then
             ! correct
          else if ( crop(i) == 0.0_r8 )then
             ! correct
          else
             call endrun(msg=' ERROR: crop has wrong values'//errMsg(__FILE__, __LINE__))
          end if
          if ( (i /= noveg) .and. (i < npcropmin) .and. &
               abs(pconv(i)+pprod10(i)+pprod100(i) - 1.0_r8) > 1.e-7_r8 )then
             call endrun(msg=' ERROR: pconv+pprod10+pprod100 do NOT sum to one.'//errMsg(__FILE__, __LINE__))
          end if
          if ( pprodharv10(i) > 1.0_r8 .or. pprodharv10(i) < 0.0_r8 )then
             call endrun(msg=' ERROR: pprodharv10 outside of range.'//errMsg(__FILE__, __LINE__))
          end if
       end do
    end if

    ! reassign taper values for shrubs - RPF
    if (taper_defaults) then
      do i = 0, npft-1
         if (woody(i) == 2._r8) then
            taper(i) = 10._r8 ! shrubs
         else if (woody(i) == 1._r8) then
            taper(i) = 200._r8
         end if
      end do
   end if

    if (masterproc) then
       write(iulog,*) 'Successfully read PFT physiological data'
       write(iulog,*)
    end if

  end subroutine pftconrd

  !-----------------------------------------------------------------------
  subroutine set_is_pft_known_to_model()
    !
    ! !DESCRIPTION:
    ! Set is_pft_known_to_model based on mergetoelmpft
    !
    ! !USES:
    !
    ! !LOCAL VARIABLES:
    integer :: m, merge_type

    character(len=*), parameter :: subname = 'set_is_pft_known_to_model'
    !-----------------------------------------------------------------------

    is_pft_known_to_model(:) = .false.

    ! NOTE(wjs, 2015-10-04) Currently, type 0 has mergetoelmpft = _FillValue in the file,
    ! so we can't handle it in the general loop below. But CLM always uses type 0, so
    ! handle it specially here.
    is_pft_known_to_model(0) = .true.

    ! NOTE(wjs, 2015-10-04) Currently, mergetoelmpft is only used for crop types.
    ! However, we handle it more generally here (treating ALL pft types), in case its use
    ! is ever extended to work with non-crop types as well.
    do m = 1, mxpft
       merge_type                        = mergetoelmpft(m)
       is_pft_known_to_model(merge_type) = .true.
    end do

  end subroutine set_is_pft_known_to_model

  !-----------------------------------------------------------------------
  subroutine set_num_cfts_known_to_model()
    !
    ! !DESCRIPTION:
    ! Set the module-level variable, num_cfts_known_to_model
    !
    ! !USES:
    !
    ! !LOCAL VARIABLES:
    integer :: m

    character(len=*), parameter :: subname = 'set_num_cfts_known_to_model'
    !-----------------------------------------------------------------------

    num_cfts_known_to_model = 0
    do m = cft_lb, cft_ub
       if (is_pft_known_to_model(m)) then
          num_cfts_known_to_model = num_cfts_known_to_model + 1
       end if
    end do

  end subroutine set_num_cfts_known_to_model

end module pftvarcon


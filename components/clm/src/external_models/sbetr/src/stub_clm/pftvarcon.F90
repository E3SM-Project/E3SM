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
  use clm_varpar  , only : mxpft, numrad, ivis, inir
  use clm_varctl  , only : iulog, use_cndv, use_vertsoilc
  use clm_varpar  , only : nlevdecomp, nsoilorder
  use clm_varctl  , only : nu_com
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
  integer :: ncorn                  !value for corn, rain fed (rf)
  integer :: ncornirrig             !value for corn, irrigated (ir)
  integer :: nscereal               !value for spring temperate cereal (rf)
  integer :: nscerealirrig          !value for spring temperate cereal (ir)
  integer :: nwcereal               !value for winter temperate cereal (rf)
  integer :: nwcerealirrig          !value for winter temperate cereal (ir)
  integer :: nsoybean               !value for soybean (rf)
  integer :: nsoybeanirrig          !value for soybean (ir)
  integer :: npcropmax              !value for last prognostic crop in list
  integer :: nc3crop                !value for generic crop (rf)
  integer :: nc3irrig               !value for irrigated generic crop (ir)

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
  real(r8), pointer :: crop(:)        => null()    !crop pft: 0. = not crop, 1. = crop pft
  real(r8), allocatable :: irrigated(:)   !irrigated pft: 0. = not, 1. = irrigated
  real(r8), allocatable :: smpso(:)       !soil water potential at full stomatal opening (mm)
  real(r8), allocatable :: smpsc(:)       !soil water potential at full stomatal closure (mm)
  real(r8), allocatable :: fnitr(:)       !foliage nitrogen limitation factor (-)

  ! begin new pft parameters for CN code
  real(r8), allocatable :: slatop(:)      !SLA at top of canopy [m^2/gC]
  real(r8), allocatable :: dsladlai(:)    !dSLA/dLAI [m^2/gC]
  real(r8), allocatable :: leafcn(:)      !leaf C:N [gC/gN]
  real(r8), allocatable :: flnr(:)        !fraction of leaf N in Rubisco [no units]
  real(r8), allocatable :: woody(:)       !woody lifeform flag (0 or 1)
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
  real(r8), allocatable :: ztopmx(:)       !parameter used in CNVegStructUpdate
  real(r8), allocatable :: laimx(:)        !parameter used in CNVegStructUpdate
  real(r8), allocatable :: gddmin(:)       !parameter used in CNPhenology
  real(r8), allocatable :: hybgdd(:)       !parameter used in CNPhenology
  real(r8), allocatable :: lfemerg(:)      !parameter used in CNPhenology
  real(r8), allocatable :: grnfill(:)      !parameter used in CNPhenology
  integer , allocatable :: mxmat(:)        !parameter used in CNPhenology
  integer , allocatable :: mnNHplantdate(:)!minimum planting date for NorthHemisphere (YYYYMMDD)
  integer , allocatable :: mxNHplantdate(:)!maximum planting date for NorthHemisphere (YYYYMMDD)
  integer , allocatable :: mnSHplantdate(:)!minimum planting date for SouthHemisphere (YYYYMMDD)
  integer , allocatable :: mxSHplantdate(:)!maximum planting date for SouthHemisphere (YYYYMMDD)
  real(r8), allocatable :: planttemp(:)    !planting temperature used in CNPhenology (K)
  real(r8), allocatable :: minplanttemp(:) !mininum planting temperature used in CNPhenology (K)
  real(r8), allocatable :: froot_leaf(:)   !allocation parameter: new fine root C per new leaf C (gC/gC)
  real(r8), allocatable :: stem_leaf(:)    !allocation parameter: new stem c per new leaf C (gC/gC)
  real(r8), allocatable :: croot_stem(:)   !allocation parameter: new coarse root C per new stem C (gC/gC)
  real(r8), allocatable :: flivewd(:)      !allocation parameter: fraction of new wood that is live (phloem and ray parenchyma) (no units)
  real(r8), allocatable :: fcur(:)         !allocation parameter: fraction of allocation that goes to currently displayed growth, remainder to storage
  real(r8), allocatable :: fcurdv(:)       !alternate fcur for use with cndv
  real(r8), allocatable :: lf_flab(:)      !leaf litter labile fraction
  real(r8), allocatable :: lf_fcel(:)      !leaf litter cellulose fraction
  real(r8), allocatable :: lf_flig(:)      !leaf litter lignin fraction
  real(r8), allocatable :: fr_flab(:)      !fine root litter labile fraction
  real(r8), allocatable :: fr_fcel(:)      !fine root litter cellulose fraction
  real(r8), allocatable :: fr_flig(:)      !fine root litter lignin fraction
  real(r8), allocatable :: leaf_long(:)    !leaf longevity (yrs)
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
  real(r8), allocatable :: fertnitro(:)    !fertilizer
  real(r8), allocatable :: fleafcn(:)      !C:N during grain fill; leaf
  real(r8), allocatable :: ffrootcn(:)     !C:N during grain fill; fine root
  real(r8), allocatable :: fstemcn(:)      !C:N during grain fill; stem
  real(r8), allocatable :: presharv(:)     !porportion of residue harvested
  real(r8), allocatable :: convfact(:)     !conversion factor to bu/acre
  real(r8), allocatable :: fyield(:)       !fraction of grain that is actually harvested
  real(r8), allocatable :: root_dmx(:)     !maximum root depth

  ! pft parameters for CNDV code
  ! from LPJ subroutine pftparameters
  real(r8), allocatable :: pftpar20(:)       !tree maximum crown area (m2)
  real(r8), allocatable :: pftpar28(:)       !min coldest monthly mean temperature
  real(r8), allocatable :: pftpar29(:)       !max coldest monthly mean temperature
  real(r8), allocatable :: pftpar30(:)       !min growing degree days (>= 5 deg C)
  real(r8), allocatable :: pftpar31(:)       !upper limit of temperature of the warmest month (twmax)

  integer, parameter :: pftname_len = 40    ! max length of pftname
  character(len=pftname_len) :: pftname(0:mxpft) !PFT description

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
  real(r8)              :: VMAX_NFIX           ! VMAX of symbiotic N2 fixation
  real(r8)              :: KM_NFIX             ! KM of symbiotic N2 fixation
  real(r8), allocatable :: VMAX_PTASE_vr(:)    ! VMAX of biochemical P production
  real(r8)              :: KM_PTASE            ! KM of biochemical P production
  real(r8)              :: lamda_ptase         ! critical value that incur biochemical production
  real(r8), allocatable :: i_vc(:)             ! intercept of photosynthesis vcmax ~ leaf N content regression model
  real(r8), allocatable :: s_vc(:)             ! slope of photosynthesis vcmax ~ leaf N content regression model
  ! new stoichiometry
  real(r8), allocatable :: leafcn_obs(:)       !leaf C:N [gC/gN]
  real(r8), allocatable :: frootcn_obs(:)      !fine root C:N (gC/gN)
  real(r8), allocatable :: livewdcn_obs(:)     !live wood (phloem and ray parenchyma) C:N (gC/gN)
  real(r8), allocatable :: deadwdcn_obs(:)     !dead wood (xylem and heartwood) C:N (gC/gN)
  real(r8), allocatable :: leafcp_obs(:)       !leaf C:P [gC/gP]
  real(r8), allocatable :: frootcp_obs(:)      !fine root C:P (gC/gP)
  real(r8), allocatable :: livewdcp_obs(:)     !live wood (phloem and ray parenchyma) C:P (gC/gP)
  real(r8), allocatable :: deadwdcp_obs(:)     !dead wood (xylem and heartwood) C:P (gC/gP)

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
    use ncdio_pio ,  only : ncd_io, ncd_pio_closefile, ncd_pio_openfile, file_desc_t, &
                            ncd_inqdid, ncd_inqdlen
    use clm_varctl,  only : paramfile, use_ed
    use clm_varctl,  only : use_crop, use_dynroot
    use clm_varcon,  only : tfrz
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
    allocate( crop          (0:mxpft) ); crop(:) = 0
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
    allocate( froot_leaf    (0:mxpft) )
    allocate( stem_leaf     (0:mxpft) )
    allocate( croot_stem    (0:mxpft) )
    allocate( flivewd       (0:mxpft) )
    allocate( fcur          (0:mxpft) )
    allocate( fcurdv        (0:mxpft) )
    allocate( lf_flab       (0:mxpft) )
    allocate( lf_fcel       (0:mxpft) )
    allocate( lf_flig       (0:mxpft) )
    allocate( fr_flab       (0:mxpft) )
    allocate( fr_fcel       (0:mxpft) )
    allocate( fr_flig       (0:mxpft) )
    allocate( leaf_long     (0:mxpft) )
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
    allocate( fertnitro     (0:mxpft) )
    allocate( fleafcn       (0:mxpft) )
    allocate( ffrootcn      (0:mxpft) )
    allocate( fstemcn       (0:mxpft) )
    allocate( presharv      (0:mxpft) )
    allocate( convfact      (0:mxpft) )
    allocate( fyield        (0:mxpft) )
    allocate( root_dmx      (0:mxpft) )
    allocate( pftpar20      (0:mxpft) )
    allocate( pftpar28      (0:mxpft) )
    allocate( pftpar29      (0:mxpft) )
    allocate( pftpar30      (0:mxpft) )
    allocate( pftpar31      (0:mxpft) )

    allocate( VMAX_PLANT_NH4(0:mxpft) )
    allocate( VMAX_PLANT_NO3(0:mxpft) )
    allocate( VMAX_PLANT_P(0:mxpft) )
    allocate( VMAX_MINSURF_P_vr(1:nlevdecomp,0:nsoilorder))
    allocate( KM_PLANT_NH4(0:mxpft) )
    allocate( KM_PLANT_NO3(0:mxpft) )
    allocate( KM_PLANT_P(0:mxpft) )
    allocate( KM_MINSURF_P_vr(1:nlevdecomp,0:nsoilorder))
    allocate( decompmicc_patch_vr (1:nlevdecomp,0:mxpft))
    allocate( VMAX_PTASE_vr(1:nlevdecomp))
    allocate( i_vc               (0:mxpft) )
    allocate( s_vc               (0:mxpft) )
    ! new stoichiometry
    allocate( leafcn_obs         (0:mxpft) )
    allocate( frootcn_obs        (0:mxpft) )
    allocate( livewdcn_obs       (0:mxpft) )
    allocate( deadwdcn_obs       (0:mxpft) )
    allocate( leafcp_obs         (0:mxpft) )
    allocate( frootcp_obs        (0:mxpft) )
    allocate( livewdcp_obs       (0:mxpft) )
    allocate( deadwdcp_obs       (0:mxpft) )


    noveg  = 1

  end subroutine pftconrd

end module pftvarcon

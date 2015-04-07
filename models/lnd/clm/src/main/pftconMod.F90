module pftconMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing vegetation constants and method to
  ! read and initialize vegetation (PFT) constants.
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils  , only : endrun
  use clm_varpar  , only : mxpft, numrad, ivis, inir
  use clm_varctl  , only : iulog, use_cndv, use_vertsoilc
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! Vegetation type constants
  !
  integer :: noveg                  ! value for not vegetated 
  integer :: ndllf_evr_tmp_tree     ! value for Needleleaf evergreen temperate tree
  integer :: ndllf_evr_brl_tree     ! value for Needleleaf evergreen boreal tree
  integer :: ndllf_dcd_brl_tree     ! value for Needleleaf deciduous boreal tree
  integer :: nbrdlf_evr_trp_tree    ! value for Broadleaf evergreen tropical tree
  integer :: nbrdlf_evr_tmp_tree    ! value for Broadleaf evergreen temperate tree
  integer :: nbrdlf_dcd_trp_tree    ! value for Broadleaf deciduous tropical tree
  integer :: nbrdlf_dcd_tmp_tree    ! value for Broadleaf deciduous temperate tree
  integer :: nbrdlf_dcd_brl_tree    ! value for Broadleaf deciduous boreal tree
  integer :: ntree                  ! value for last type of tree
  integer :: nbrdlf_evr_shrub       ! value for Broadleaf evergreen shrub
  integer :: nbrdlf_dcd_tmp_shrub   ! value for Broadleaf deciduous temperate shrub
  integer :: nbrdlf_dcd_brl_shrub   ! value for Broadleaf deciduous boreal shrub
  integer :: nc3_arctic_grass       ! value for C3 arctic grass
  integer :: nc3_nonarctic_grass    ! value for C3 non-arctic grass
  integer :: nc4_grass              ! value for C4 grass
  integer :: npcropmin              ! value for first crop
  integer :: ncorn                  ! value for corn, rain fed (rf)
  integer :: ncornirrig             ! value for corn, irrigated (ir)
  integer :: nscereal               ! value for spring temperate cereal (rf)
  integer :: nscerealirrig          ! value for spring temperate cereal (ir)
  integer :: nwcereal               ! value for winter temperate cereal (rf)
  integer :: nwcerealirrig          ! value for winter temperate cereal (ir)
  integer :: nsoybean               ! value for soybean (rf)
  integer :: nsoybeanirrig          ! value for soybean (ir)
  integer :: npcropmax              ! value for last prognostic crop in list
  integer :: nc3crop                ! value for generic crop (rf)
  integer :: nc3irrig               ! value for irrigated generic crop (ir)

  ! !PUBLIC TYPES:
  type, public :: pftcon_type

     integer , allocatable :: noveg         (:)   ! value for not vegetated
     integer , allocatable :: tree          (:)   ! tree or not?

     real(r8), allocatable :: dleaf         (:)   ! characteristic leaf dimension (m)
     real(r8), allocatable :: c3psn         (:)   ! photosynthetic pathway: 0. = c4, 1. = c3
     real(r8), allocatable :: xl            (:)   ! leaf/stem orientation index
     real(r8), allocatable :: rhol          (:,:) ! leaf reflectance: 1=vis, 2=nir
     real(r8), allocatable :: rhos          (:,:) ! stem reflectance: 1=vis, 2=nir
     real(r8), allocatable :: taul          (:,:) ! leaf transmittance: 1=vis, 2=nir
     real(r8), allocatable :: taus          (:,:) ! stem transmittance: 1=vis, 2=nir
     real(r8), allocatable :: z0mr          (:)   ! ratio of momentum roughness length to canopy top height (-)
     real(r8), allocatable :: displar       (:)   ! ratio of displacement height to canopy top height (-)
     real(r8), allocatable :: roota_par     (:)   ! CLM rooting distribution parameter [1/m]
     real(r8), allocatable :: rootb_par     (:)   ! CLM rooting distribution parameter [1/m]
     real(r8), allocatable :: crop          (:)   ! crop pft: 0. = not crop, 1. = crop pft
     real(r8), allocatable :: irrigated     (:)   ! irrigated pft: 0. = not, 1. = irrigated
     real(r8), allocatable :: smpso         (:)   ! soil water potential at full stomatal opening (mm)
     real(r8), allocatable :: smpsc         (:)   ! soil water potential at full stomatal closure (mm)
     real(r8), allocatable :: fnitr         (:)   ! foliage nitrogen limitation factor (-)

     !  CN code
     real(r8), allocatable :: dwood         (:)   ! wood density (gC/m3)
     real(r8), allocatable :: slatop        (:)   ! SLA at top of canopy [m^2/gC]
     real(r8), allocatable :: dsladlai      (:)   ! dSLA/dLAI [m^2/gC]
     real(r8), allocatable :: leafcn        (:)   ! leaf C:N [gC/gN]
     real(r8), allocatable :: flnr          (:)   ! fraction of leaf N in Rubisco [no units]
     real(r8), allocatable :: woody         (:)   ! woody lifeform flag (0 or 1)
     real(r8), allocatable :: lflitcn       (:)   ! leaf litter C:N (gC/gN)
     real(r8), allocatable :: frootcn       (:)   ! fine root C:N (gC/gN)
     real(r8), allocatable :: livewdcn      (:)   ! live wood (phloem and ray parenchyma) C:N (gC/gN)
     real(r8), allocatable :: deadwdcn      (:)   ! dead wood (xylem and heartwood) C:N (gC/gN)
     real(r8), allocatable :: grperc        (:)   ! growth respiration parameter
     real(r8), allocatable :: grpnow        (:)   ! growth respiration parameter
     real(r8), allocatable :: rootprof_beta (:)   ! CLM rooting distribution parameter for C and N inputs [unitless]

     !  crop
     real(r8), allocatable :: graincn       (:)   ! grain C:N (gC/gN)
     real(r8), allocatable :: mxtmp         (:)   ! parameter used in accFlds
     real(r8), allocatable :: baset         (:)   ! parameter used in accFlds
     real(r8), allocatable :: declfact      (:)   ! parameter used in CNAllocation
     real(r8), allocatable :: bfact         (:)   ! parameter used in CNAllocation
     real(r8), allocatable :: aleaff        (:)   ! parameter used in CNAllocation
     real(r8), allocatable :: arootf        (:)   ! parameter used in CNAllocation
     real(r8), allocatable :: astemf        (:)   ! parameter used in CNAllocation
     real(r8), allocatable :: arooti        (:)   ! parameter used in CNAllocation
     real(r8), allocatable :: fleafi        (:)   ! parameter used in CNAllocation
     real(r8), allocatable :: allconsl      (:)   ! parameter used in CNAllocation
     real(r8), allocatable :: allconss      (:)   ! parameter used in CNAllocation
     real(r8), allocatable :: ztopmx        (:)   ! parameter used in CNVegStructUpdate
     real(r8), allocatable :: laimx         (:)   ! parameter used in CNVegStructUpdate
     real(r8), allocatable :: gddmin        (:)   ! parameter used in CNPhenology
     real(r8), allocatable :: hybgdd        (:)   ! parameter used in CNPhenology
     real(r8), allocatable :: lfemerg       (:)   ! parameter used in CNPhenology
     real(r8), allocatable :: grnfill       (:)   ! parameter used in CNPhenology
     integer , allocatable :: mxmat         (:)   ! parameter used in CNPhenology
     integer , allocatable :: mnNHplantdate (:)   ! minimum planting date for NorthHemisphere (YYYYMMDD)
     integer , allocatable :: mxNHplantdate (:)   ! maximum planting date for NorthHemisphere (YYYYMMDD)
     integer , allocatable :: mnSHplantdate (:)   ! minimum planting date for SouthHemisphere (YYYYMMDD)
     integer , allocatable :: mxSHplantdate (:)   ! maximum planting date for SouthHemisphere (YYYYMMDD)
     real(r8), allocatable :: planttemp     (:)   ! planting temperature used in CNPhenology (K)
     real(r8), allocatable :: minplanttemp  (:)   ! mininum planting temperature used in CNPhenology (K)
     real(r8), allocatable :: froot_leaf    (:)   ! allocation parameter: new fine root C per new leaf C (gC/gC) 
     real(r8), allocatable :: stem_leaf     (:)   ! allocation parameter: new stem c per new leaf C (gC/gC)
     real(r8), allocatable :: croot_stem    (:)   ! allocation parameter: new coarse root C per new stem C (gC/gC)
     real(r8), allocatable :: flivewd       (:)   ! allocation parameter: fraction of new wood that is live (phloem and ray parenchyma) (no units)
     real(r8), allocatable :: fcur          (:)   ! allocation parameter: fraction of allocation that goes to currently displayed growth, remainder to storage
     real(r8), allocatable :: fcurdv        (:)   ! alternate fcur for use with cndv
     real(r8), allocatable :: lf_flab       (:)   ! leaf litter labile fraction
     real(r8), allocatable :: lf_fcel       (:)   ! leaf litter cellulose fraction
     real(r8), allocatable :: lf_flig       (:)   ! leaf litter lignin fraction
     real(r8), allocatable :: fr_flab       (:)   ! fine root litter labile fraction
     real(r8), allocatable :: fr_fcel       (:)   ! fine root litter cellulose fraction
     real(r8), allocatable :: fr_flig       (:)   ! fine root litter lignin fraction
     real(r8), allocatable :: leaf_long     (:)   ! leaf longevity (yrs)
     real(r8), allocatable :: evergreen     (:)   ! binary flag for evergreen leaf habit (0 or 1)
     real(r8), allocatable :: stress_decid  (:)   ! binary flag for stress-deciduous leaf habit (0 or 1)
     real(r8), allocatable :: season_decid  (:)   ! binary flag for seasonal-deciduous leaf habit (0 or 1)
     real(r8), allocatable :: pconv         (:)   ! proportion of deadstem to conversion flux
     real(r8), allocatable :: pprod10       (:)   ! proportion of deadstem to 10-yr product pool
     real(r8), allocatable :: pprod100      (:)   ! proportion of deadstem to 100-yr product pool
     real(r8), allocatable :: pprodharv10   (:)   ! harvest mortality proportion of deadstem to 10-yr pool

     ! pft paraemeters for fire code
     real(r8), allocatable :: cc_leaf       (:)
     real(r8), allocatable :: cc_lstem      (:)
     real(r8), allocatable :: cc_dstem      (:)
     real(r8), allocatable :: cc_other      (:)
     real(r8), allocatable :: fm_leaf       (:)
     real(r8), allocatable :: fm_lstem      (:)
     real(r8), allocatable :: fm_dstem      (:)
     real(r8), allocatable :: fm_other      (:)
     real(r8), allocatable :: fm_root       (:)
     real(r8), allocatable :: fm_lroot      (:)
     real(r8), allocatable :: fm_droot      (:)
     real(r8), allocatable :: fsr_pft       (:)
     real(r8), allocatable :: fd_pft        (:)

     ! pft parameters for crop code
     real(r8), allocatable :: fertnitro     (:)   ! fertilizer
     real(r8), allocatable :: fleafcn       (:)   ! C:N during grain fill; leaf
     real(r8), allocatable :: ffrootcn      (:)   ! C:N during grain fill; fine root
     real(r8), allocatable :: fstemcn       (:)   ! C:N during grain fill; stem

     ! pft parameters for CNDV code (from LPJ subroutine pftparameters)
     real(r8), allocatable :: pftpar20      (:)   ! tree maximum crown area (m2)
     real(r8), allocatable :: pftpar28      (:)   ! min coldest monthly mean temperature
     real(r8), allocatable :: pftpar29      (:)   ! max coldest monthly mean temperature
     real(r8), allocatable :: pftpar30      (:)   ! min growing degree days (>= 5 deg C)
     real(r8), allocatable :: pftpar31      (:)   ! upper limit of temperature of the warmest month (twmax)

   contains

     procedure, public  :: Init
     procedure, private :: InitAllocate   
     procedure, private :: InitRead

  end type pftcon_type

  type(pftcon_type), public :: pftcon ! pft type constants structure

  integer, parameter :: pftname_len = 40         ! max length of pftname       
  character(len=pftname_len) :: pftname(0:mxpft) ! PFT description

  real(r8), parameter :: reinickerp = 1.6_r8     ! parameter in allometric equation
  real(r8), parameter :: dwood  = 2.5e5_r8       ! cn wood density (gC/m3); lpj:2.0e5
  real(r8), parameter :: allom1 = 100.0_r8       ! parameters in
  real(r8), parameter :: allom2 =  40.0_r8       ! ...allometric
  real(r8), parameter :: allom3 =   0.5_r8       ! ...equations
  real(r8), parameter :: allom1s = 250.0_r8      ! modified for shrubs by
  real(r8), parameter :: allom2s =   8.0_r8      ! X.D.Z
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this)

    class(pftcon_type) :: this

    call this%InitAllocate()
    call this%InitRead()

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate (this)
    !
    ! !DESCRIPTION:
    ! Read and initialize vegetation (PFT) constants
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(pftcon_type) :: this
    !-----------------------------------------------------------------------

    allocate( this%noveg         (0:mxpft)); this%noveg (:)   =huge(1)
    allocate( this%tree          (0:mxpft)); this%tree  (:)   =huge(1)

    allocate( this%dleaf         (0:mxpft) )       
    allocate( this%c3psn         (0:mxpft) )       
    allocate( this%xl            (0:mxpft) )          
    allocate( this%rhol          (0:mxpft,numrad) ) 
    allocate( this%rhos          (0:mxpft,numrad) ) 
    allocate( this%taul          (0:mxpft,numrad) ) 
    allocate( this%taus          (0:mxpft,numrad) ) 
    allocate( this%z0mr          (0:mxpft) )        
    allocate( this%displar       (0:mxpft) )     
    allocate( this%roota_par     (0:mxpft) )   
    allocate( this%rootb_par     (0:mxpft) )   
    allocate( this%crop          (0:mxpft) )        
    allocate( this%irrigated     (0:mxpft) )   
    allocate( this%smpso         (0:mxpft) )       
    allocate( this%smpsc         (0:mxpft) )       
    allocate( this%fnitr         (0:mxpft) )       
    allocate( this%slatop        (0:mxpft) )      
    allocate( this%dsladlai      (0:mxpft) )    
    allocate( this%leafcn        (0:mxpft) )      
    allocate( this%flnr          (0:mxpft) )        
    allocate( this%woody         (0:mxpft) )       
    allocate( this%lflitcn       (0:mxpft) )      
    allocate( this%frootcn       (0:mxpft) )      
    allocate( this%livewdcn      (0:mxpft) )     
    allocate( this%deadwdcn      (0:mxpft) )     
    allocate( this%grperc        (0:mxpft) )       
    allocate( this%grpnow        (0:mxpft) )       
    allocate( this%rootprof_beta (0:mxpft) )
    allocate( this%graincn       (0:mxpft) )      
    allocate( this%mxtmp         (0:mxpft) )        
    allocate( this%baset         (0:mxpft) )        
    allocate( this%declfact      (0:mxpft) )     
    allocate( this%bfact         (0:mxpft) )        
    allocate( this%aleaff        (0:mxpft) )       
    allocate( this%arootf        (0:mxpft) )       
    allocate( this%astemf        (0:mxpft) )       
    allocate( this%arooti        (0:mxpft) )       
    allocate( this%fleafi        (0:mxpft) )       
    allocate( this%allconsl      (0:mxpft) )     
    allocate( this%allconss      (0:mxpft) )     
    allocate( this%ztopmx        (0:mxpft) )       
    allocate( this%laimx         (0:mxpft) )        
    allocate( this%gddmin        (0:mxpft) )       
    allocate( this%hybgdd        (0:mxpft) )       
    allocate( this%lfemerg       (0:mxpft) )      
    allocate( this%grnfill       (0:mxpft) )      
    allocate( this%mxmat         (0:mxpft) )        
    allocate( this%mnNHplantdate (0:mxpft) )
    allocate( this%mxNHplantdate (0:mxpft) )
    allocate( this%mnSHplantdate (0:mxpft) )
    allocate( this%mxSHplantdate (0:mxpft) )
    allocate( this%planttemp     (0:mxpft) )    
    allocate( this%minplanttemp  (0:mxpft) ) 
    allocate( this%froot_leaf    (0:mxpft) )   
    allocate( this%stem_leaf     (0:mxpft) )    
    allocate( this%croot_stem    (0:mxpft) )   
    allocate( this%flivewd       (0:mxpft) )      
    allocate( this%fcur          (0:mxpft) )         
    allocate( this%fcurdv        (0:mxpft) )       
    allocate( this%lf_flab       (0:mxpft) )      
    allocate( this%lf_fcel       (0:mxpft) )      
    allocate( this%lf_flig       (0:mxpft) )      
    allocate( this%fr_flab       (0:mxpft) )      
    allocate( this%fr_fcel       (0:mxpft) )      
    allocate( this%fr_flig       (0:mxpft) )      
    allocate( this%leaf_long     (0:mxpft) )   
    allocate( this%evergreen     (0:mxpft) )    
    allocate( this%stress_decid  (0:mxpft) ) 
    allocate( this%season_decid  (0:mxpft) ) 
    allocate( this%dwood         (0:mxpft) )
    allocate( this%pconv         (0:mxpft) )        
    allocate( this%pprod10       (0:mxpft) )      
    allocate( this%pprod100      (0:mxpft) )     
    allocate( this%pprodharv10   (0:mxpft) )  
    allocate( this%cc_leaf       (0:mxpft) )
    allocate( this%cc_lstem      (0:mxpft) )
    allocate( this%cc_dstem      (0:mxpft) )
    allocate( this%cc_other      (0:mxpft) )
    allocate( this%fm_leaf       (0:mxpft) )
    allocate( this%fm_lstem      (0:mxpft) )
    allocate( this%fm_dstem      (0:mxpft) )
    allocate( this%fm_other      (0:mxpft) )
    allocate( this%fm_root       (0:mxpft) )
    allocate( this%fm_lroot      (0:mxpft) )
    allocate( this%fm_droot      (0:mxpft) )
    allocate( this%fsr_pft       (0:mxpft) )
    allocate( this%fd_pft        (0:mxpft) )
    allocate( this%fertnitro     (0:mxpft) )
    allocate( this%fleafcn       (0:mxpft) )  
    allocate( this%ffrootcn      (0:mxpft) ) 
    allocate( this%fstemcn       (0:mxpft) )  
    allocate( this%pftpar20      (0:mxpft) )   
    allocate( this%pftpar28      (0:mxpft) )   
    allocate( this%pftpar29      (0:mxpft) )   
    allocate( this%pftpar30      (0:mxpft) )   
    allocate( this%pftpar31      (0:mxpft) )   

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitRead(this)
    !
    ! !DESCRIPTION:
    ! Read and initialize vegetation (PFT) constants
    !
    ! !USES:
    use shr_log_mod , only : errMsg => shr_log_errMsg
    use fileutils   , only : getfil
    use ncdio_pio   , only : ncd_io, ncd_pio_closefile, ncd_pio_openfile, file_desc_t
    use ncdio_pio   , only : ncd_inqdid, ncd_inqdlen
    use clm_varctl  , only : paramfile, use_ed
    use spmdMod     , only : masterproc
    use EDPftvarcon , only : EDpftconrd
    !
    ! !ARGUMENTS:
    class(pftcon_type) :: this
    !
    ! !LOCAL VARIABLES:
    character(len=256) :: locfn                ! local file name
    integer            :: i,n,m                ! loop indices
    integer            :: ier                  ! error code
    type(file_desc_t)  :: ncid                 ! pio netCDF file id
    integer            :: dimid                ! netCDF dimension id
    integer            :: npft                 ! number of pfts on pft-physiology file
    logical            :: readv                ! read variable in or not
    character(len=32)  :: subname = 'InitRead' ! subroutine name
    character(len=pftname_len) :: expected_pftnames(0:mxpft) 
    !-----------------------------------------------------------------------
    !
    ! Expected PFT names: The names expected on the paramfile file and the order they are expected to be in.
    ! NOTE: similar types are assumed to be together, first trees (ending with broadleaf_deciduous_boreal_tree
    !       then shrubs, ending with broadleaf_deciduous_boreal_shrub, then grasses starting with c3_arctic_grass
    !       and finally crops, ending with soybean
    ! DO NOT CHANGE THE ORDER -- WITHOUT MODIFYING OTHER PARTS OF THE CODE WHERE THE ORDER MATTERS!

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

    ! Set specific vegetation type values

    if (masterproc) then
       write(iulog,*) 'Attempting to read PFT physiological data .....'
    end if
    call getfil (paramfile, locfn, 0)
    call ncd_pio_openfile (ncid, trim(locfn), 0)
    call ncd_inqdid(ncid, 'pft', dimid)
    call ncd_inqdlen(ncid, dimid, npft)

    call ncd_io('pftname',pftname, 'read', ncid, readvar=readv, posNOTonfile=.true.) 
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('z0mr', this%z0mr, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('displar', this%displar, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('dleaf', this%dleaf, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('c3psn', this%c3psn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('rholvis', this%rhol(:,ivis), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('rholnir', this%rhol(:,inir), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('rhosvis', this%rhos(:,ivis), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('rhosnir', this% rhos(:,inir), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('taulvis', this%taul(:,ivis), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('taulnir', this%taul(:,inir), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('tausvis', this%taus(:,ivis), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('tausnir', this%taus(:,inir), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('xl', this%xl, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('roota_par', this%roota_par, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('rootb_par', this%rootb_par, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('slatop', this%slatop, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('dsladlai', this%dsladlai, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('leafcn', this%leafcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('flnr', this%flnr, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('smpso', this%smpso, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('smpsc', this%smpsc, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('fnitr', this%fnitr, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('woody', this%woody, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('lflitcn', this%lflitcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('frootcn', this%frootcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('livewdcn', this%livewdcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('deadwdcn', this%deadwdcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('grperc', this%grperc, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('grpnow', this%grpnow, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('froot_leaf', this%froot_leaf, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('stem_leaf', this%stem_leaf, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('croot_stem', this%croot_stem, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('flivewd', this%flivewd, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('fcur', this%fcur, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('fcurdv', this%fcurdv, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('lf_flab', this%lf_flab, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('lf_fcel', this%lf_fcel, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('lf_flig', this%lf_flig, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('fr_flab', this%fr_flab, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('fr_fcel', this%fr_fcel, 'read', ncid, readvar=readv, posNOTonfile=.true.)    
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('fr_flig', this%fr_flig, 'read', ncid, readvar=readv, posNOTonfile=.true.)    
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('leaf_long', this%leaf_long, 'read', ncid, readvar=readv, posNOTonfile=.true.)    
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('evergreen', this%evergreen, 'read', ncid, readvar=readv, posNOTonfile=.true.)    
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('stress_decid', this%stress_decid, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('season_decid', this%season_decid, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('pftpar20', this%pftpar20, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('pftpar28', this%pftpar28, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('pftpar29', this%pftpar29, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('pftpar30', this%pftpar30, 'read', ncid, readvar=readv, posNOTonfile=.true.)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('pftpar31', this%pftpar31, 'read', ncid, readvar=readv, posNOTonfile=.true.)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('fertnitro', this%fertnitro, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('fleafcn', this%fleafcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('ffrootcn', this%ffrootcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('fstemcn', this%fstemcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    if (use_vertsoilc) then
       call ncd_io('rootprof_beta', this%rootprof_beta, 'read', ncid, readvar=readv, posNOTonfile=.true.)
       if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))
    end if

    call ncd_io('pconv', this%pconv, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('pprod10', this%pprod10, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('pprodharv10', this%pprodharv10, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('pprod100', this%pprod100, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('graincn', this%graincn, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('mxtmp', this%mxtmp, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('baset', this%baset, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('declfact', this%declfact, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('bfact', this%bfact, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('aleaff', this%aleaff, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('arootf', this%arootf, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('astemf', this%astemf, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('arooti', this%arooti, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('fleafi', this%fleafi, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('allconsl', this%allconsl, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('allconss', this%allconss, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('crop', this%crop, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('irrigated', this%irrigated, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('ztopmx', this%ztopmx, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('laimx', this%laimx, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('gddmin', this%gddmin, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('hybgdd', this%hybgdd, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('lfemerg', this%lfemerg, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('grnfill', this%grnfill, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('mxmat', this%mxmat, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('cc_leaf', this% cc_leaf, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('cc_lstem', this%cc_lstem, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('cc_dstem', this%cc_dstem, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('cc_other', this%cc_other, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('fm_leaf', this% fm_leaf, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('fm_lstem', this%fm_lstem, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('fm_dstem', this%fm_dstem, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('fm_other', this%fm_other, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('fm_root', this% fm_root, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('fm_lroot', this%fm_lroot, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('fm_droot', this%fm_droot, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('fsr_pft', this% fsr_pft, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('fd_pft', this%  fd_pft, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('planting_temp', this%planttemp, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('min_planting_temp', this%minplanttemp, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('min_NH_planting_date', this%mnNHplantdate, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('min_SH_planting_date', this%mnSHplantdate, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('max_NH_planting_date', this%mxNHplantdate, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    call ncd_io('max_SH_planting_date', this%mxSHplantdate, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(__FILE__, __LINE__))

    !
    ! Constants
    !
    !MV (10-08-14) TODO is this right - used to be numpft - is it okay to set it to mxpft?
    do m = 0,mxpft 
       this%dwood(m) = dwood
       if (m <= ntree) then
          this%tree(m) = 1
       else
          this%tree(m) = 0
       end if
    end do
    !
    ! ED variables
    !
    if ( use_ed ) then
       ! The following sets the module variable EDpftcon_inst in EDPftcon
       call EDpftconrd ( ncid )
    endif
       
    call ncd_pio_closefile(ncid)

    do i = 0, mxpft
       if (.not. use_ed)then
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
    end do

    ntree                = nbrdlf_dcd_brl_tree  ! value for last type of tree
    npcropmin            = ncorn                ! first prognostic crop
    npcropmax            = nsoybeanirrig        ! last prognostic crop in list

    if (use_cndv) then
       this%fcur(:) = this%fcurdv(:)
    end if
    !
    ! Do some error checking, but not if ED is on.
    !
    ! FIX(SPM,032414) double check if some of these should be on...

    if( .not. use_ed ) then
       if ( npcropmax /= mxpft )then
          call endrun(msg=' ERROR: npcropmax is NOT the last value'//errMsg(__FILE__, __LINE__))
       end if
       do i = 0, mxpft
          if ( this%irrigated(i) == 1.0_r8  .and. &
               (i == nc3irrig      .or. &
                i == ncornirrig    .or. &
                i == nscerealirrig .or. &
                i == nwcerealirrig .or. &
                i == nsoybeanirrig) )then
             ! correct
          else if ( this%irrigated(i) == 0.0_r8 )then
             ! correct
          else
             call endrun(msg=' ERROR: irrigated has wrong values'//errMsg(__FILE__, __LINE__))
          end if
          if (      this%crop(i) == 1.0_r8 .and. (i >= nc3crop .and. i <= npcropmax) )then
             ! correct
          else if ( this%crop(i) == 0.0_r8 )then
             ! correct
          else
             call endrun(msg=' ERROR: crop has wrong values'//errMsg(__FILE__, __LINE__))
          end if
          if ( (i /= noveg) .and. (i < npcropmin) .and. &
               abs(this%pconv(i) + this%pprod10(i) + this%pprod100(i) - 1.0_r8) > 1.e-7_r8 )then
             call endrun(msg=' ERROR: pconv+pprod10+pprod100 do NOT sum to one.'//errMsg(__FILE__, __LINE__))
          end if
          if ( this%pprodharv10(i) > 1.0_r8 .or. this%pprodharv10(i) < 0.0_r8 )then
             call endrun(msg=' ERROR: pprodharv10 outside of range.'//errMsg(__FILE__, __LINE__))
          end if
       end do
    end if

    if (masterproc) then
       write(iulog,*) 'Successfully read PFT physiological data'
       write(iulog,*)
    end if

  end subroutine InitRead

end module pftconMod


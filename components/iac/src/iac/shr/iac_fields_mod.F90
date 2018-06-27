
Module iac_fields_mod
  
!---------------------------------------------------------------------------
!BOP
!
! !MODULE: iac_fields_mod
!
!  iac coupling fields definition
!
! !DESCRIPTION:
!
! !USES:

  implicit none
  SAVE
  private                              ! By default make data private

! !PUBLIC TYPES:

   public iac_cdata_type

! !PUBLIC MEMBER FUNCTIONS:
  public :: iac_fields_init

! !PUBLIC DATA MEMBERS:
  !--- this should be namelist or something like it ---

  real*8,  parameter, public :: iac_spval = -999.0
  integer, parameter, public :: iac_gcam_nreg =  14
  integer, parameter, public :: iac_gcam_nsector =  12
  integer, parameter, public :: iac_gcamoemis_nemis=  1
  integer, parameter, public :: iac_gcam_naez =  18
  integer, parameter, public :: iac_gcam_ncrops =  27
  integer, parameter, public :: iac_gcam_timestep =  5
  integer, parameter, public :: iac_gcam_ioyears =  15
  integer, parameter, public :: iac_gcamo_ntime =  2
  integer, parameter, public :: iac_glm_nx  = 720
  integer, parameter, public :: iac_glm_ny  = 360
  integer, parameter, public :: iac_iaco_npfts =  17
  integer, parameter, public :: iac_iac_npfts  = 16
  character(len=*), parameter, public :: iac_gcam2glm_map = 'some_map_file'

! ! cdata datatype
  type iac_cdata_type
    character(len=640),pointer :: c(:)
    real*8 , pointer :: r(:)
    integer, pointer :: i(:)
    logical, pointer :: l(:)
  end type
  
!-----------------------------------------------------------------

  integer, parameter, public :: iac_eclock_size  =  8
  integer, parameter, public :: iac_eclock_ymd   =  1
  integer, parameter, public :: iac_eclock_tod   =  2
  integer, parameter, public :: iac_eclock_dt    =  3
  integer, parameter, public :: iac_eclock_AclmC =  4
  integer, parameter, public :: iac_eclock_Agcam =  5
  integer, parameter, public :: iac_eclock_Aglm  =  6
  integer, parameter, public :: iac_eclock_Agcamsetden =  7

  integer, parameter, public :: iac_cdata_size            = 32
  !--- characters ---
  integer, parameter, public :: iac_cdatac_casename        =  1
  integer, parameter, public :: iac_cdatac_clm2gcam        =  2
  integer, parameter, public :: iac_cdatac_ibclmfile       =  3
  integer, parameter, public :: iac_cdatac_clmcbfndir      =  4
  integer, parameter, public :: iac_cdatac_gcam2glm_basecrop= 5
  integer, parameter, public :: iac_cdatac_gcam2glm_basepast= 6
  integer, parameter, public :: iac_cdatac_gcam2glm_baseothr= 7
  integer, parameter, public :: iac_cdatac_gcam2glm_aezmap =  8
  integer, parameter, public :: iac_cdatac_gcam2glm_basebiomass= 9
  integer, parameter, public :: iac_cdatac_gcam2emisfile_co2base2000 = 10
  integer, parameter, public :: iac_cdatac_gcam2emisfile_grid720x360 = 11
  integer, parameter, public :: iac_cdatac_gcam2emisfile_grid288x192 = 12
  integer, parameter, public :: iac_cdatac_gcam2emisfile_co2shipbase2000 = 13
  integer, parameter, public :: iac_cdatac_gcam2emisfile_lut720x360map = 14
  integer, parameter, public :: iac_cdatac_gcam2emisfile_downscaleinfo = 15
  integer, parameter, public :: iac_cdatac_gcam2emisfile_rcp45allsteps = 16


  !--- iofields ---
  integer, parameter, public :: iac_cdataio_gcami_crops     =  1
  integer, parameter, public :: iac_cdataio_gcami_flds      =  2
  integer, parameter, public :: iac_cdataio_gcamo_flds      =  3
  integer, parameter, public :: iac_cdataio_gcamo_regs      =  4
  integer, parameter, public :: iac_cdataio_glmi_flds       =  5
  integer, parameter, public :: iac_cdataio_glmo_flds       =  6

  !--- reals ---
  !--- integers ---
  integer, parameter, public :: iac_cdatai_logunit         =  1
  integer, parameter, public :: iac_cdatai_iac_nx          =  2
  integer, parameter, public :: iac_cdatai_iac_ny          =  3
  integer, parameter, public :: iac_cdatai_iac_size        =  4
  integer, parameter, public :: iac_cdatai_glm_nx          =  5
  integer, parameter, public :: iac_cdatai_glm_ny          =  6
  integer, parameter, public :: iac_cdatai_glm_size        =  7
  integer, parameter, public :: iac_cdatai_gcam_yr1        =  8
  integer, parameter, public :: iac_cdatai_gcam_yr2        =  9
  integer, parameter, public :: iac_cdatai_gcam_nreg       =  10
  integer, parameter, public :: iac_cdatai_gcam_naez       =  11
  integer, parameter, public :: iac_cdatai_gcam_timestep   =  12
  integer, parameter, public :: iac_cdatai_gcamo_ntime     =  13
  integer, parameter, public :: iac_cdatai_gcamo_nflds     =  14
  integer, parameter, public :: iac_cdatai_gcamo_size      =  15
  integer, parameter, public :: iac_cdatai_gcami_nflds     =  16
  integer, parameter, public :: iac_cdatai_gcami_size      =  17
  integer, parameter, public :: iac_cdatai_gcami_ncrops    =  18
  integer, parameter, public :: iac_cdatai_gcamoemis_size  =  19
  !--- logicals ---
  integer, parameter, public :: iac_cdatal_rest            =  1
  integer, parameter, public :: iac_cdatal_iac_present     =  2
  integer, parameter, public :: iac_cdatal_iac_prognostic  =  3
  integer, parameter, public :: iac_cdatal_glm_present     =  4
  integer, parameter, public :: iac_cdatal_glm_prognostic  =  5
  integer, parameter, public :: iac_cdatal_gcam_present    =  6
  integer, parameter, public :: iac_cdatal_gcam_prognostic =  7
  integer, parameter, public :: iac_cdatal_fastiac         =  8
  integer, parameter, public :: iac_cdatal_npphr           =  9
  integer, parameter, public :: iac_cdatal_initrun         =  10
  integer, parameter, public :: iac_cdatal_sneakermode     =  11
  integer, parameter, public :: iac_cdatal_nocarbonscale   =  12
  integer, parameter, public :: iac_cdatal_co2flux_coupling=  13
  integer, parameter, public :: iac_cdatal_write_rest      =  14

  integer           , public :: iac_iaci_nflds
  integer           , public :: iac_iaci_area
  integer, pointer  , public :: iac_iaci_pft(:)
  integer, pointer  , public :: iac_iaci_AGB(:)
  integer, pointer  , public :: iac_iaci_BGB(:)

  integer           , public :: iac_iaco_nflds
  integer, pointer  , public :: iac_iaco_pft(:)

  character(len=80) , pointer  , public :: iac_gcami_crop_names(:)
  character(len=80) , pointer  , public :: iac_gcami_fld_names(:)
  character(len=80) , pointer  , public :: iac_gcamo_reg_names(:)
  character(len=80) , pointer  , public :: iac_gcamo_fld_names(:)
  character(len=80) , pointer  , public :: iac_glmi_fld_names(:)
  character(len=80) , pointer  , public :: iac_glmo_fld_names(:)

  integer           , public :: iac_gcami_nflds
  integer           , public :: iac_gcami_biomass
  integer           , public :: iac_gcami_eucalyptus
  integer           , public :: iac_gcami_miscanthus
  integer           , public :: iac_gcami_willow
  integer           , public :: iac_gcami_Corn
  integer           , public :: iac_gcami_FiberCrop
  integer           , public :: iac_gcami_FodderGrass
  integer           , public :: iac_gcami_FodderHerb
  integer           , public :: iac_gcami_Forest
  integer           , public :: iac_gcami_Grassland
  integer           , public :: iac_gcami_Jatropha
  integer           , public :: iac_gcami_MiscCrop
  integer           , public :: iac_gcami_OilCrop
  integer           , public :: iac_gcami_OtherArableLand
  integer           , public :: iac_gcami_OtherGrain
  integer           , public :: iac_gcami_PalmFruit
  integer           , public :: iac_gcami_Pasture
  integer           , public :: iac_gcami_Rice
  integer           , public :: iac_gcami_RockIceDesert
  integer           , public :: iac_gcami_Root_Tuber
  integer           , public :: iac_gcami_Shrubland
  integer           , public :: iac_gcami_SugarCrop
  integer           , public :: iac_gcami_Tundra
  integer           , public :: iac_gcami_UnmanagedForest
  integer           , public :: iac_gcami_UnmanagedPasture
  integer           , public :: iac_gcami_UrbanLand
  integer           , public :: iac_gcami_Wheat
  integer,save      , public :: iac_gcami_above_ground_carbon=1
  integer,save      , public :: iac_gcami_below_ground_carbon=2

  integer, pointer  , public :: iac_gcami_cdens(:)

  integer           , public :: iac_gcamo_nflds
  integer           , public :: iac_gcamo_buildup
  integer           , public :: iac_gcamo_crop
  integer           , public :: iac_gcamo_pasture
  integer           , public :: iac_gcamo_woodharv
  integer           , public :: iac_gcamo_forest
  integer           , public :: iac_gcamo_other

  integer           , public :: iac_gcamo_reg_usa
  integer           , public :: iac_gcamo_reg_canada
  integer           , public :: iac_gcamo_reg_western_europe
  integer           , public :: iac_gcamo_reg_japan
  integer           , public :: iac_gcamo_reg_australia_nz
  integer           , public :: iac_gcamo_reg_former_soviet_union
  integer           , public :: iac_gcamo_reg_china
  integer           , public :: iac_gcamo_reg_middle_east
  integer           , public :: iac_gcamo_reg_africa
  integer           , public :: iac_gcamo_reg_latin_america
  integer           , public :: iac_gcamo_reg_southeast_asia
  integer           , public :: iac_gcamo_reg_eastern_europe
  integer           , public :: iac_gcamo_reg_korea
  integer           , public :: iac_gcamo_reg_india

  integer           , public :: iac_glmi_nflds
  integer           , public :: iac_glmi_natveg
  integer           , public :: iac_glmi_cropland
  integer           , public :: iac_glmi_pasture
  integer           , public :: iac_glmi_woodharv

  integer           , public :: iac_glmo_nflds
  integer           , public :: iac_glmo_gcrop
  integer           , public :: iac_glmo_gpast
  integer           , public :: iac_glmo_gothr
  integer           , public :: iac_glmo_gsecd
  integer           , public :: iac_glmo_gfvh1
  integer           , public :: iac_glmo_gfvh2
  integer           , public :: iac_glmo_gfsh1
  integer           , public :: iac_glmo_gfsh2
  integer           , public :: iac_glmo_gfsh3

! !REVISION HISTORY:
! Author: T Craig


! !PRIVATE DATA MEMBERS:

!EOP
!===============================================================
contains
!===============================================================

!---------------------------------------------------------------------------
!BOP

! !IROUTINE: iac_fields_init()

! !INTERFACE:
  subroutine iac_fields_init()

! !DESCRIPTION:
! Initialize coupling field indices

! !USES:
    implicit none

! !ARGUMENTS:

    logical :: npp_on=.true.

! !LOCAL VARIABLES:

    integer :: n,i
    character(len=*),parameter :: subname='(iac_fields_init)'

! !REVISION HISTORY:
! Author: T Craig

!EOP
!-----------------------------------------------------------------------

    iac_iaci_nflds = 3*iac_iac_npfts + 1
    iac_iaci_area = 1
    allocate(iac_iaci_pft(iac_iac_npfts))
    allocate(iac_iaci_AGB(iac_iac_npfts))
    allocate(iac_iaci_BGB(iac_iac_npfts))
    do n = 1,iac_iac_npfts
       iac_iaci_pft(n) = 0*iac_iac_npfts + n + 1
       iac_iaci_AGB(n) = 1*iac_iac_npfts + n + 1
       iac_iaci_BGB(n) = 2*iac_iac_npfts + n + 1
    enddo

    iac_iaco_nflds = iac_iac_npfts
    allocate(iac_iaco_pft(iac_iac_npfts))
    do n = 1,iac_iac_npfts
       iac_iaco_pft(n) = n
    enddo

    allocate(iac_gcami_cdens(iac_iaco_npfts))
    do n = 1,iac_iaco_npfts
       iac_gcami_cdens(n) = n
    enddo

    iac_gcamo_nflds    = 6
    allocate(iac_gcamo_fld_names(iac_gcamo_nflds))
    iac_gcamo_buildup  = 1
    iac_gcamo_fld_names(1) = "buildup"
    iac_gcamo_crop     = 2
    iac_gcamo_fld_names(2) = "crop"
    iac_gcamo_pasture  = 3
    iac_gcamo_fld_names(3) = "pasture"
    iac_gcamo_woodharv = 4
    iac_gcamo_fld_names(4) = "woodharv"
    iac_gcamo_forest   = 5
    iac_gcamo_fld_names(5) = "forest"
    iac_gcamo_other    = 6
    iac_gcamo_fld_names(6) = "other"

    allocate(iac_gcamo_reg_names(iac_gcam_nreg))
    iac_gcamo_reg_usa                 = 1
    iac_gcamo_reg_names(1)            = "USA" 
    iac_gcamo_reg_canada              = 2
    iac_gcamo_reg_names(2)            = "Canada"
    iac_gcamo_reg_western_europe      = 3
    iac_gcamo_reg_names(3)            = "Western Europe"
    iac_gcamo_reg_japan               = 4
    iac_gcamo_reg_names(4)            = "Japan"
    iac_gcamo_reg_australia_nz        = 5
    iac_gcamo_reg_names(5)            = "Australia NZ"
    iac_gcamo_reg_former_soviet_union = 6
    iac_gcamo_reg_names(6)            = "Former Soviet Union = 6"
    iac_gcamo_reg_china               = 7
    iac_gcamo_reg_names(7)            = "China"
    iac_gcamo_reg_middle_east         = 8
    iac_gcamo_reg_names(8)            = "Middle East"
    iac_gcamo_reg_africa              = 9
    iac_gcamo_reg_names(9)            = "Africa"
    iac_gcamo_reg_latin_america       = 10
    iac_gcamo_reg_names(10)            = "Latin America"
    iac_gcamo_reg_southeast_asia      = 11
    iac_gcamo_reg_names(11)            = "Southeast Asia"
    iac_gcamo_reg_eastern_europe      = 12
    iac_gcamo_reg_names(12)            = "Eastern Europe"
    iac_gcamo_reg_korea               = 13
    iac_gcamo_reg_names(13)            = "Korea"
    iac_gcamo_reg_india               = 14
    iac_gcamo_reg_names(14)            = "India"

    iac_gcami_nflds                   = 2
    allocate(iac_gcami_fld_names(iac_gcami_nflds))
    if(npp_on)then
       iac_gcami_fld_names(1) = "NPP"
       iac_gcami_fld_names(2) = "HR"
    else
       iac_gcami_fld_names(1) = "AboveGroundCarbon"
       iac_gcami_fld_names(2) = "BelowGroundCarbon"
    end if

    i=0
    allocate(iac_gcami_crop_names(iac_gcam_ncrops))
    i=i+1
    iac_gcami_biomass       = i
    iac_gcami_crop_names(i) = "biomass"
    i=i+1
    iac_gcami_corn          = i
    iac_gcami_crop_names(i) = "Corn"
    i=i+1
    iac_gcami_eucalyptus    = i
    iac_gcami_crop_names(i) = "eucalyptus"
    i=i+1
    iac_gcami_fibercrop     = i
    iac_gcami_crop_names(i) = "FiberCrop"
    i=i+1
    iac_gcami_FodderGrass   = i
    iac_gcami_crop_names(i) = "FodderGrass"
    i=i+1
    iac_gcami_FodderHerb    = i
    iac_gcami_crop_names(i) = "FodderHerb"
    i=i+1
    iac_gcami_Forest        = i
    iac_gcami_crop_names(i) = "Forest"
    i=i+1
    iac_gcami_Grassland     = i
    iac_gcami_crop_names(i) = "Grassland"
    i=i+1
    iac_gcami_Jatropha      = i
    iac_gcami_crop_names(i) = "Jatropha"
    i=i+1
    iac_gcami_miscanthus    = i
    iac_gcami_crop_names(i) = "miscanthus"
    i=i+1
    iac_gcami_MiscCrop      = i
    iac_gcami_crop_names(i) = "MiscCrop"
    i=i+1
    iac_gcami_OilCrop       = i
    iac_gcami_crop_names(i) = "OilCrop"
    i=i+1
    iac_gcami_OtherArableLand = i
    iac_gcami_crop_names(i) = "OtherArableLand"
    i=i+1
    iac_gcami_OtherGrain    = i
    iac_gcami_crop_names(i) = "OtherGrain"
    i=i+1
    iac_gcami_PalmFruit     = i
    iac_gcami_crop_names(i) = "PalmFruit"
    i=i+1
    iac_gcami_Pasture       = i
    iac_gcami_crop_names(i) = "Pasture"
    i=i+1
    iac_gcami_Rice          = i
    iac_gcami_crop_names(i) = "Rice"
    i=i+1
    iac_gcami_RockIceDesert = i
    iac_gcami_crop_names(i) = "RockIceDesert"
    i=i+1
    iac_gcami_Root_Tuber    = i
    iac_gcami_crop_names(i) = "Root_Tuber"
    i=i+1
    iac_gcami_Shrubland     = i
    iac_gcami_crop_names(i) = "Shrubland"
    i=i+1
    iac_gcami_SugarCrop     = i
    iac_gcami_crop_names(i) = "SugarCrop"
    i=i+1
    iac_gcami_Tundra        = i
    iac_gcami_crop_names(i) = "Tundra"
    i=i+1
    iac_gcami_UnmanagedForest = i
    iac_gcami_crop_names(i) = "UnmanagedForest"
    i=i+1
    iac_gcami_UnmanagedPasture = i
    iac_gcami_crop_names(i) = "UnmanagedPasture"
    i=i+1
    iac_gcami_UrbanLand     = i
    iac_gcami_crop_names(i) = "UrbanLand"
    i=i+1
    iac_gcami_Wheat         = i
    iac_gcami_crop_names(i) = "Wheat"
    i=i+1
    iac_gcami_willow        = i
    iac_gcami_crop_names(i) = "willow"

    iac_glmi_nflds     = 3
    allocate(iac_glmi_fld_names(iac_glmi_nflds))
    iac_glmi_fld_names(1) = "natveg"
    iac_glmi_natveg    = 1
    iac_glmi_fld_names(2) = "cropland"
    iac_glmi_cropland  = 2
    iac_glmi_fld_names(3) = "pasture"
    iac_glmi_pasture   = 3

    iac_glmo_nflds     = 9
    allocate(iac_glmo_fld_names(iac_glmo_nflds))
    iac_glmo_fld_names(1) = "gcrop"
    iac_glmo_gcrop     = 1
    iac_glmo_fld_names(2) = "gpast"
    iac_glmo_gpast     = 2
    iac_glmo_fld_names(3) = "gothr"
    iac_glmo_gothr     = 3
    iac_glmo_fld_names(4) = "gsecd"
    iac_glmo_gsecd     = 4
    iac_glmo_fld_names(5) = "gfvh1"
    iac_glmo_gfvh1     = 5
    iac_glmo_fld_names(6) = "gfvh2"
    iac_glmo_gfvh2     = 6
    iac_glmo_fld_names(7) = "gfsh1"
    iac_glmo_gfsh1     = 7
    iac_glmo_fld_names(8) = "gfsh2"
    iac_glmo_gfsh2     = 8
    iac_glmo_fld_names(9) = "gfsh3"
    iac_glmo_gfsh3     = 9

  end subroutine iac_fields_init

!---------------------------------------------------------------------------


end module iac_fields_mod

module SoilStateType

  !------------------------------------------------------------------------------
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
  use decompMod       , only : bounds_type
  use abortutils      , only : endrun
  use spmdMod         , only : mpicom, MPI_INTEGER, masterproc
  use ncdio_pio       , only : file_desc_t, ncd_defvar, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
  use ncdio_pio       , only : ncd_pio_openfile, ncd_inqfdims, ncd_pio_closefile, ncd_inqdid, ncd_inqdlen
  use clm_varpar      , only : more_vertlayers, numpft, numrad 
  use clm_varpar      , only : nlevsoi, nlevgrnd, nlevlak, nlevsoifl, nlayer, nlayert, nlevurb, nlevsno
  use landunit_varcon , only : istice, istdlak, istwet, istsoil, istcrop, istice_mec
  use column_varcon   , only : icol_roof, icol_sunwall, icol_shadewall, icol_road_perv, icol_road_imperv 
  use clm_varcon      , only : zsoi, dzsoi, zisoi, spval
  use clm_varcon      , only : secspday, pc, mu, denh2o, denice, grlnd
  use clm_varctl      , only : use_cn, use_lch4,use_dynroot, use_fates
  use clm_varctl      , only : use_var_soil_thick
  use clm_varctl      , only : iulog, fsurdat, hist_wrtch4diag
  use CH4varcon       , only : allowlakeprod
  use LandunitType    , only : lun_pp                
  use ColumnType      , only : col_pp                
  use VegetationType       , only : veg_pp                
  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  type, public :: soilstate_type

     ! sand/ clay/ organic matter
     real(r8), pointer :: sandfrac_patch       (:)   ! patch sand fraction
     real(r8), pointer :: clayfrac_patch       (:)   ! patch clay fraction
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
     real(r8), pointer :: watmin_col           (:,:) ! col minimum volumetric soil water (nlevsoi)
     real(r8), pointer :: sucsat_col           (:,:) ! col minimum soil suction (mm) (nlevgrnd)
     real(r8), pointer :: sucmin_col           (:,:) ! col minimum allowable soil liquid suction pressure (mm) [Note: sucmin_col is a negative value, while sucsat_col is a positive quantity]
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

     ! roots
     real(r8), pointer :: rootr_patch          (:,:) ! patch effective fraction of roots in each soil layer (nlevgrnd)
     real(r8), pointer :: rootr_col            (:,:) ! col effective fraction of roots in each soil layer (nlevgrnd)  
     real(r8), pointer :: rootfr_col           (:,:) ! col fraction of roots in each soil layer (nlevgrnd) 
     real(r8), pointer :: rootfr_patch         (:,:) ! patch fraction of roots in each soil layer (nlevgrnd)
     real(r8), pointer :: rootr_road_perv_col  (:,:) ! col effective fraction of roots in each soil layer of urban pervious road
     real(r8), pointer :: rootfr_road_perv_col (:,:) ! col effective fraction of roots in each soil layer of urban pervious road
     real(r8), pointer :: root_depth_patch     (:)   ! rooting depth of each PFT (m)

   contains

     procedure, public  :: Init
     procedure, private :: InitAllocate
     procedure, private :: InitHistory  
     procedure, private :: InitCold    
     procedure, public  :: Restart
     procedure, public  :: InitColdGhost

  end type soilstate_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(soilstate_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    call this%InitCold(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !ARGUMENTS:
    class(soilstate_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg
    integer :: begc_all, endc_all
    !------------------------------------------------------------------------

    begp     = bounds%begp    ; endp     = bounds%endp
    begc     = bounds%begc    ; endc     = bounds%endc
    begg     = bounds%begg    ; endg     = bounds%endg
    begc_all = bounds%begc_all; endc_all = bounds%endc_all

    allocate(this%mss_frc_cly_vld_col  (begc:endc))                     ; this%mss_frc_cly_vld_col  (:)   = nan
    allocate(this%sandfrac_patch       (begp:endp))                     ; this%sandfrac_patch       (:)   = nan
    allocate(this%clayfrac_patch       (begp:endp))                     ; this%clayfrac_patch       (:)   = nan
    allocate(this%cellorg_col          (begc:endc,nlevgrnd))            ; this%cellorg_col          (:,:) = nan 
    allocate(this%cellsand_col         (begc:endc,nlevgrnd))            ; this%cellsand_col         (:,:) = nan 
    allocate(this%cellclay_col         (begc:endc,nlevgrnd))            ; this%cellclay_col         (:,:) = nan 
    allocate(this%bd_col               (begc:endc,nlevgrnd))            ; this%bd_col               (:,:) = nan

    allocate(this%hksat_col            (begc_all:endc_all,nlevgrnd))    ; this%hksat_col            (:,:) = spval
    allocate(this%hksat_min_col        (begc:endc,nlevgrnd))            ; this%hksat_min_col        (:,:) = spval
    allocate(this%hk_l_col             (begc:endc,nlevgrnd))            ; this%hk_l_col             (:,:) = nan   
    allocate(this%smp_l_col            (begc:endc,nlevgrnd))            ; this%smp_l_col            (:,:) = nan   
    allocate(this%smpmin_col           (begc:endc))                     ; this%smpmin_col           (:)   = nan

    allocate(this%bsw_col              (begc_all:endc_all,nlevgrnd))    ; this%bsw_col              (:,:) = nan
    allocate(this%watsat_col           (begc_all:endc_all,nlevgrnd))    ; this%watsat_col           (:,:) = nan
    allocate(this%watdry_col           (begc:endc,nlevgrnd))            ; this%watdry_col           (:,:) = spval
    allocate(this%watopt_col           (begc:endc,nlevgrnd))            ; this%watopt_col           (:,:) = spval
    allocate(this%watfc_col            (begc:endc,nlevgrnd))            ; this%watfc_col            (:,:) = nan
    allocate(this%watmin_col           (begc:endc,nlevgrnd))            ; this%watmin_col           (:,:) = nan
    allocate(this%sucsat_col           (begc:endc,nlevgrnd))            ; this%sucsat_col           (:,:) = spval
    allocate(this%sucmin_col           (begc:endc,nlevgrnd))            ; this%sucmin_col           (:,:) = spval
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

    allocate(this%rootr_patch          (begp:endp,1:nlevgrnd))          ; this%rootr_patch          (:,:) = nan
    allocate(this%rootr_col            (begc:endc,nlevgrnd))            ; this%rootr_col            (:,:) = nan
    allocate(this%rootr_road_perv_col  (begc:endc,1:nlevgrnd))          ; this%rootr_road_perv_col  (:,:) = nan
    allocate(this%rootfr_patch         (begp:endp,1:nlevgrnd))          ; this%rootfr_patch         (:,:) = nan
    allocate(this%rootfr_col           (begc:endc,1:nlevgrnd))          ; this%rootfr_col           (:,:) = nan 
    allocate(this%rootfr_road_perv_col (begc:endc,1:nlevgrnd))          ; this%rootfr_road_perv_col (:,:) = nan
    allocate(this%root_depth_patch     (begp:endp))                     ; this%root_depth_patch     (:)   = spval

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! History fields initialization
    !
    ! !USES:
    use histFileMod   , only: hist_addfld1d, hist_addfld2d, no_snow_normal
    !
    ! !ARGUMENTS:
    class(soilstate_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc
    integer :: begp, endp
    character(10)     :: active
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    !---------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc


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

    if (use_cn) then
       this%rootfr_patch(begp:endp,:) = spval
       call hist_addfld2d (fname='ROOTFR', units='proportion', type2d='levgrnd', &
            avgflag='A', long_name='fraction of roots in each soil layer', &
            ptr_patch=this%rootfr_patch, default='inactive')
    end if

    if (use_cn) then
       this%rootr_patch(begp:endp,:) = spval
       call hist_addfld2d (fname='ROOTR', units='proportion', type2d='levgrnd', &
            avgflag='A', long_name='effective fraction of roots in each soil layer', &
            ptr_patch=this%rootr_patch, default='inactive')
    end if

    if (use_cn) then
       this%rootr_col(begc:endc,:) = spval
       call hist_addfld2d (fname='ROOTR_COLUMN', units='proportion', type2d='levgrnd', &
            avgflag='A', long_name='effective fraction of roots in each soil layer', &
            ptr_col=this%rootr_col, default='inactive')
       
    end if

    if (use_dynroot) then
       this%root_depth_patch(begp:endp) = spval
       call hist_addfld1d (fname='ROOT_DEPTH', units="m", &
            avgflag='A', long_name='rooting depth', &
            ptr_patch=this%root_depth_patch, default='inactive' )
    end if

    if (use_cn .or. use_fates) then
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

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! Initialize module surface albedos to reasonable values
    !
    ! !USES:
    use pftvarcon           , only : noveg, roota_par, rootb_par
    use fileutils           , only : getfil
    use organicFileMod      , only : organicrd 
    use SharedParamsMod   , only : ParamsShareInst
    use FuncPedotransferMod , only : pedotransf, get_ipedof
    use RootBiophysMod      , only : init_vegrootfr
    !
    ! !ARGUMENTS:
    class(soilstate_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
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
    integer            :: nlevbed                       ! # of layers above bedrock
    integer            :: ipedof  
    integer            :: begc, endc
    integer            :: begg, endg
    real(r8), parameter :: min_liquid_pressure = -10132500._r8 ! Minimum soil liquid water pressure [mm]
    !-----------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    do c = bounds%begc, bounds%endc
       this%smpmin_col(c) = -1.e8_r8
    end do

    ! --------------------------------------------------------------------
    ! Initialize root fraction (computing from surface, d is depth in meter):
    ! --------------------------------------------------------------------

    ! Currently pervious road has same properties as soil
    do c = bounds%begc, bounds%endc
       l = col_pp%landunit(c)

       if (lun_pp%urbpoi(l) .and. col_pp%itype(c) == icol_road_perv) then 
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
       if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then
          this%rootfr_col (c,nlevsoi+1:nlevgrnd) = 0._r8
       else if (lun_pp%itype(l) == istdlak .and. allowlakeprod) then
          this%rootfr_col (c,:) = spval
       else  ! Inactive CH4 columns
          this%rootfr_col (c,:) = spval
       end if
    end do

   ! Initialize root fraction 
   
   call init_vegrootfr(bounds, nlevsoi, nlevgrnd, &
        col_pp%nlevbed(bounds%begc:bounds%endc)    , &
        this%rootfr_patch(bounds%begp:bounds%endp,1:nlevgrnd))

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

    organic_max = ParamsShareInst%organic_max

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

    do p = bounds%begp,bounds%endp
       g = veg_pp%gridcell(p)
       if ( sand3d(g,1)+clay3d(g,1) == 0.0_r8 )then
          if ( any( sand3d(g,:)+clay3d(g,:) /= 0.0_r8 ) )then
             call endrun(msg='found depth points that do NOT sum to zero when surface does'//&
                  errMsg(__FILE__, __LINE__)) 
          end if
          sand3d(g,:) = 1.0_r8
          clay3d(g,:) = 1.0_r8
       end if
       if ( any( sand3d(g,:)+clay3d(g,:) == 0.0_r8 ) )then
          call endrun(msg='after setting, found points sum to zero'//errMsg(__FILE__, __LINE__)) 
       end if

       this%sandfrac_patch(p) = sand3d(g,1)/100.0_r8
       this%clayfrac_patch(p) = clay3d(g,1)/100.0_r8
    end do

    ! Read fmax

    allocate(gti(bounds%begg:bounds%endg))
    call ncd_io(ncid=ncid, varname='FMAX', flag='read', data=gti, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=' ERROR: FMAX NOT on surfdata file'//errMsg(__FILE__, __LINE__)) 
    end if
    do c = bounds%begc, bounds%endc
       g = col_pp%gridcell(c)
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
       g = col_pp%gridcell(c)
       l = col_pp%landunit(c)

       if (lun_pp%itype(l)==istwet .or. lun_pp%itype(l)==istice .or. lun_pp%itype(l)==istice_mec) then

          do lev = 1,nlevgrnd
             this%bsw_col(c,lev)    = spval
             this%watsat_col(c,lev) = spval
             this%watfc_col(c,lev)  = spval
             this%watmin_col(c,lev) = spval
             this%hksat_col(c,lev)  = spval
             this%sucsat_col(c,lev) = spval
             this%sucmin_col(c,lev) = spval
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
             if (lun_pp%itype(l)==istwet .and. lev > nlevsoi) then
                this%csol_col(c,lev) = csol_bedrock
             else
                this%csol_col(c,lev)= spval
             endif
          end do

       else if (lun_pp%urbpoi(l) .and. (col_pp%itype(c) /= icol_road_perv) .and. (col_pp%itype(c) /= icol_road_imperv) )then

          ! Urban Roof, sunwall, shadewall properties set to special value
          do lev = 1,nlevgrnd
             this%watsat_col(c,lev) = spval
             this%watfc_col(c,lev)  = spval
             this%watmin_col(c,lev) = spval
             this%bsw_col(c,lev)    = spval
             this%hksat_col(c,lev)  = spval
             this%sucsat_col(c,lev) = spval
             this%sucmin_col(c,lev) = spval
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
             ! Number of soil layers in hydrologically active columns = NLEV2BED
	     nlevbed = col_pp%nlevbed(c)
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

             if (lun_pp%itype(l) == istdlak) then

                if (lev <= nlevsoi) then
                   this%cellsand_col(c,lev) = sand
                   this%cellclay_col(c,lev) = clay
                   this%cellorg_col(c,lev)  = om_frac*organic_max
                end if

             else if (lun_pp%itype(l) /= istdlak) then  ! soil columns of both urban and non-urban types

                if (lun_pp%urbpoi(l)) then
                   om_frac = 0._r8 ! No organic matter for urban
                end if

                if (lev <= nlevbed) then
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

                if (lev > nlevbed) then
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

                this%sucmin_col(c,lev) = min_liquid_pressure

                this%watmin_col(c,lev) = &
                     this%watsat_col(c,lev)*(-min_liquid_pressure/this%sucsat_col(c,lev))**(-1._r8/this%bsw_col(c,lev))

             end if
          end do

          ! Urban pervious and impervious road
          if (col_pp%itype(c) == icol_road_imperv) then
             ! Impervious road layers -- same as above except set watdry and watopt as missing
             do lev = 1,nlevgrnd
                this%watdry_col(c,lev) = spval 
                this%watopt_col(c,lev) = spval 
             end do
          else if (col_pp%itype(c) == icol_road_perv) then 
             ! pervious road layers  - set in UrbanInitTimeConst
          end if

       end if
    end do

    ! --------------------------------------------------------------------
    ! Set soil hydraulic and thermal properties: lake
    ! --------------------------------------------------------------------

    do c = bounds%begc, bounds%endc
       g = col_pp%gridcell(c)
       l = col_pp%landunit(c)

       if (lun_pp%itype(l)==istdlak) then

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
       g = col_pp%gridcell(c)

       this%gwc_thr_col(c) = 0.17_r8 + 0.14_r8 * clay3d(g,1) * 0.01_r8
       this%mss_frc_cly_vld_col(c) = min(clay3d(g,1) * 0.01_r8, 0.20_r8)
    end do

    ! --------------------------------------------------------------------
    ! Deallocate memory
    ! --------------------------------------------------------------------

    deallocate(sand3d, clay3d, organic3d)
    deallocate(zisoifl, zsoifl, dzsoifl)

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    !
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
    class(soilstate_type) :: this
    type(bounds_type), intent(in)    :: bounds
    type(file_desc_t), intent(inout) :: ncid
    character(len=*) , intent(in)    :: flag
    !
    ! !LOCAL VARIABLES:
    logical          :: readvar   ! determine if variable is on initial file
    logical          :: readrootfr = .false.
    !-----------------------------------------------------------------------

    if(use_dynroot) then
       call restartvar(ncid=ncid, flag=flag, varname='root_depth', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='root depth', units='m', &
            interpinic_flag='interp', readvar=readvar, data=this%root_depth_patch)

       call restartvar(ncid=ncid, flag=flag, varname='rootfr', xtype=ncd_double,  &
            dim1name='pft', dim2name='levgrnd', switchdim=.true., &
            long_name='root fraction', units='', &
            interpinic_flag='interp', readvar=readrootfr, data=this%rootfr_patch)
    else
       readrootfr = .false.
    end if
    if (flag=='read' .and. .not. readrootfr) then
       if (masterproc) then
          write(iulog,*) "can't find rootfr in restart (or initial) file..."
          write(iulog,*) "Initialize rootfr to default"
       end if
       call init_vegrootfr(bounds, nlevsoi, nlevgrnd, &
            col_pp%nlevbed(bounds%begc:bounds%endc), &
            this%rootfr_patch(bounds%begp:bounds%endp,1:nlevgrnd))
    end if
  end subroutine Restart


  !------------------------------------------------------------------------
#ifdef USE_PETSC_LIB
  subroutine InitColdGhost(this, bounds_proc)
    !
    ! !DESCRIPTION:
    ! Assign soil properties for ghost/halo columns
    !
    ! !USES:
    use domainLateralMod       , only : ExchangeColumnLevelGhostData
    use shr_infnan_mod         , only : shr_infnan_isnan
    use shr_infnan_mod         , only : isnan => shr_infnan_isnan
    use landunit_varcon        , only : max_lunit
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(soilstate_type)            :: this
    type(bounds_type), intent(in)    :: bounds_proc
    !
    integer             :: c,j                     ! indices
    integer             :: nvals_col               ! number of values per subgrid category
    integer             :: beg_idx, end_idx        ! begin/end index for accessing values in data_send/data_recv
    real(r8) , parameter:: FILL_VALUE = -999999.d0 ! temporary
    real(r8) , pointer  :: data_send_col(:)        ! data sent by local mpi rank
    real(r8) , pointer  :: data_recv_col(:)        ! data received by local mpi rank

    ! Number of values per soil column
    nvals_col = 4*nlevgrnd ! (watsat + hksat + bsw + sucsat) * nlevgrnd

    ! Allocate value
    allocate(data_send_col((bounds_proc%endc     - bounds_proc%begc     + 1)*nvals_col))
    allocate(data_recv_col((bounds_proc%endc_all - bounds_proc%begc_all + 1)*nvals_col))

    ! Assemble the data to send
    do c = bounds_proc%begc, bounds_proc%endc

       beg_idx = (c - bounds_proc%begc)*nvals_col

       do j = 1, nlevgrnd

          beg_idx = beg_idx + 1
          if (.not. isnan(this%watsat_col(c,j)) .and. this%watsat_col(c,j) /= spval) then
             data_send_col(beg_idx) = this%watsat_col(c,j)
          else
             data_send_col(beg_idx) = FILL_VALUE
          endif

          beg_idx = beg_idx + 1
          if (.not. isnan(this%hksat_col(c,j)) .and. this%hksat_col(c,j) /= spval) then
             data_send_col(beg_idx) = this%hksat_col(c,j)
          else
             data_send_col(beg_idx) = FILL_VALUE
          endif

          beg_idx = beg_idx + 1
          if (.not. isnan(this%bsw_col(c,j)) .and. this%bsw_col(c,j) /= spval) then
             data_send_col(beg_idx) = this%bsw_col(c,j)
          else
             data_send_col(beg_idx) = FILL_VALUE
          endif

          beg_idx = beg_idx + 1
          if (.not. isnan(this%sucsat_col(c,j)) .and. this%sucsat_col(c,j) /= spval) then
             data_send_col(beg_idx) = this%sucsat_col(c,j)
          else
             data_send_col(beg_idx) = FILL_VALUE
          endif
       enddo
    enddo

    ! Send the data
    call ExchangeColumnLevelGhostData(bounds_proc, nvals_col, data_send_col, data_recv_col)

    ! Assign data corresponding to ghost/halo soil columns
    do c = bounds_proc%endc + 1, bounds_proc%endc_all
       beg_idx = (c - bounds_proc%begc)*nvals_col
       do j = 1, nlevgrnd
          beg_idx = beg_idx + 1
          this%watsat_col(c,j) = data_recv_col(beg_idx)

          beg_idx = beg_idx + 1
          this%hksat_col(c,j) = data_recv_col(beg_idx)

          beg_idx = beg_idx + 1
          this%bsw_col(c,j) = data_recv_col(beg_idx)

          beg_idx = beg_idx + 1
          this%sucsat_col(c,j) = data_recv_col(beg_idx)
       enddo
    enddo

    ! Free up memory
    deallocate(data_send_col)
    deallocate(data_recv_col)

  end subroutine InitColdGhost

#else

  !------------------------------------------------------------------------
  subroutine InitColdGhost(this, bounds_proc)
    !
    ! !DESCRIPTION:
    ! Assign soil properties for ghost/halo columns
    !
    ! !USES:
    implicit none
    !
    ! !ARGUMENTS:
    class(soilstate_type)            :: this
    type(bounds_type), intent(in)    :: bounds_proc

    character(len=*), parameter :: subname = 'InitColdGhost'

    call endrun(msg='ERROR ' // trim(subname) //': Requires '//&
         'PETSc, but the code was compiled without -DUSE_PETSC_LIB')

  end subroutine InitColdGhost

#endif
  !------------------------------------------------------------------------

end module SoilStateType

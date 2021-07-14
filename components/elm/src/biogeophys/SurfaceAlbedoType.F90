module SurfaceAlbedoType

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use decompMod      , only : bounds_type
  use elm_varpar     , only : numrad, nlevcan, nlevsno
  use abortutils     , only : endrun
  use elm_varctl     , only : fsurdat, iulog
  use elm_varcon     , only : grlnd
  use ColumnType     , only : col_pp
   

  !
  ! !PUBLIC TYPES:
  implicit none
  save

  ! !PUBLIC DATA MEMBERS:
  ! The ELM default albice values are too high.
  ! Full-spectral albedo for land ice is ~0.5 (Paterson, Physics of Glaciers, 1994, p. 59)
  ! This is the value used in CAM3 by Pritchard et al., GRL, 35, 2008.

  ! albedo land ice by waveband (1=vis, 2=nir)
  real(r8), public  :: albice(numrad) = (/ 0.80_r8, 0.55_r8 /)
  !$acc declare copyin(albice)

  ! namelist default setting for inputting alblakwi
  real(r8), public  :: lake_melt_icealb(numrad) = (/ 0.10_r8, 0.10_r8/)
  !$acc declare copyin(lake_melt_icealb)

  ! albedo frozen lakes by waveband (1=vis, 2=nir)
  ! unclear what the reference is for this
  real(r8), public :: alblak(numrad) = (/0.60_r8, 0.40_r8/)
  !$acc declare copyin(alblak)

  ! albedo of melting lakes due to puddling, open water, or white ice
  ! From D. Mironov (2010) Boreal Env. Research
  ! To revert albedo of melting lakes to the cold snow-free value, set
  ! lake_melt_icealb namelist to 0.60, 0.40 like alblak above.
  real(r8), public :: alblakwi(numrad)
  !$acc declare create(alblakwi)

  !
  ! DATA FUNCTIONS:
  real(r8), allocatable :: albsat (:,:) ! wet soil albedo by color class and waveband (1=vis,2=nir)
  real(r8), allocatable :: albdry (:,:) ! dry soil albedo by color class and waveband (1=vis,2=nir)
  integer , allocatable :: isoicol(:)  ! column soil color class
  !$acc declare create(albsat )
  !$acc declare create(albdry )
  !$acc declare create(isoicol)


  public :: SurfaceAlbedoInitTimeConst
  !
  ! !PUBLIC DATA MEMBERS:
  type, public :: surfalb_type

     real(r8), pointer :: coszen_col           (:)   => null() ! col cosine of solar zenith angle
     real(r8), pointer :: albd_patch           (:,:) => null() ! patch surface albedo (direct)   (numrad)
     real(r8), pointer :: albi_patch           (:,:) => null() ! patch surface albedo (diffuse)  (numrad)
     real(r8), pointer :: albgrd_pur_col       (:,:) => null() ! col pure snow ground direct albedo     (numrad)
     real(r8), pointer :: albgri_pur_col       (:,:) => null() ! col pure snow ground diffuse albedo    (numrad)
     real(r8), pointer :: albgrd_bc_col        (:,:) => null() ! col ground direct  albedo without BC   (numrad)
     real(r8), pointer :: albgri_bc_col        (:,:) => null() ! col ground diffuse albedo without BC   (numrad)
     real(r8), pointer :: albgrd_oc_col        (:,:) => null() ! col ground direct  albedo without OC   (numrad)
     real(r8), pointer :: albgri_oc_col        (:,:) => null() ! col ground diffuse albedo without OC   (numrad)
     real(r8), pointer :: albgrd_dst_col       (:,:) => null() ! col ground direct  albedo without dust (numrad)
     real(r8), pointer :: albgri_dst_col       (:,:) => null() ! col ground diffuse albedo without dust (numrad)
     real(r8), pointer :: albgrd_col           (:,:) => null() ! col ground albedo (direct)  (numrad)
     real(r8), pointer :: albgri_col           (:,:) => null() ! col ground albedo (diffuse) (numrad)
     real(r8), pointer :: albsod_col           (:,:) => null() ! col soil albedo: direct  (col,bnd) [frc]
     real(r8), pointer :: albsoi_col           (:,:) => null() ! col soil albedo: diffuse (col,bnd) [frc]
     real(r8), pointer :: albsnd_hst_col       (:,:) => null() ! col snow albedo, direct , for history files (col,bnd) [frc]
     real(r8), pointer :: albsni_hst_col       (:,:) => null() ! col snow albedo, diffuse, for history files (col,bnd) [frc]

     real(r8), pointer :: ftdd_patch           (:,:) => null() ! patch down direct flux below canopy per unit direct flx    (numrad)
     real(r8), pointer :: ftid_patch           (:,:) => null() ! patch down diffuse flux below canopy per unit direct flx   (numrad)
     real(r8), pointer :: ftii_patch           (:,:) => null() ! patch down diffuse flux below canopy per unit diffuse flx  (numrad)
     real(r8), pointer :: fabd_patch           (:,:) => null() ! patch flux absorbed by canopy per unit direct flux         (numrad)
     real(r8), pointer :: fabd_sun_patch       (:,:) => null() ! patch flux absorbed by sunlit canopy per unit direct flux  (numrad)
     real(r8), pointer :: fabd_sha_patch       (:,:) => null() ! patch flux absorbed by shaded canopy per unit direct flux  (numrad)
     real(r8), pointer :: fabi_patch           (:,:) => null() ! patch flux absorbed by canopy per unit diffuse flux        (numrad)
     real(r8), pointer :: fabi_sun_patch       (:,:) => null() ! patch flux absorbed by sunlit canopy per unit diffuse flux (numrad)
     real(r8), pointer :: fabi_sha_patch       (:,:) => null() ! patch flux absorbed by shaded canopy per unit diffuse flux (numrad)
     real(r8), pointer :: fabd_sun_z_patch     (:,:) => null() ! patch absorbed sunlit leaf direct  PAR (per unit lai+sai) for each canopy layer
     real(r8), pointer :: fabd_sha_z_patch     (:,:) => null() ! patch absorbed shaded leaf direct  PAR (per unit lai+sai) for each canopy layer
     real(r8), pointer :: fabi_sun_z_patch     (:,:) => null() ! patch absorbed sunlit leaf diffuse PAR (per unit lai+sai) for each canopy layer
     real(r8), pointer :: fabi_sha_z_patch     (:,:) => null() ! patch absorbed shaded leaf diffuse PAR (per unit lai+sai) for each canopy layer
     real(r8), pointer :: flx_absdv_col        (:,:) => null() ! col absorbed flux per unit incident direct flux:  VIS (col,lyr) [frc]
     real(r8), pointer :: flx_absdn_col        (:,:) => null() ! col absorbed flux per unit incident direct flux:  NIR (col,lyr) [frc]
     real(r8), pointer :: flx_absiv_col        (:,:) => null() ! col absorbed flux per unit incident diffuse flux: VIS (col,lyr) [frc]
     real(r8), pointer :: flx_absin_col        (:,:) => null() ! col absorbed flux per unit incident diffuse flux: NIR (col,lyr) [frc]

     real(r8) , pointer :: fsun_z_patch        (:,:) => null() ! patch patch sunlit fraction of canopy layer
     real(r8) , pointer :: tlai_z_patch        (:,:) => null() ! patch tlai increment for canopy layer
     real(r8) , pointer :: tsai_z_patch        (:,:) => null() ! patch tsai increment for canopy layer
     integer  , pointer :: ncan_patch          (:)   => null() ! patch number of canopy layers
     integer  , pointer :: nrad_patch          (:)   => null() ! patch number of canopy layers, above snow for radiative transfer
     real(r8) , pointer :: vcmaxcintsun_patch  (:)   => null() ! patch leaf to canopy scaling coefficient, sunlit leaf vcmax
     real(r8) , pointer :: vcmaxcintsha_patch  (:)   => null() ! patch leaf to canopy scaling coefficient, shaded leaf vcmax

   contains

     procedure, public  :: Init
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold
     procedure, public  :: Restart


  end type surfalb_type
  !-----------------------------------------------------------------------



contains

  !-----------------------------------------------------------------------
  subroutine SurfaceAlbedoInitTimeConst(bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module time constant variables
    !
    ! !USES:
    use shr_log_mod, only : errMsg => shr_log_errMsg
    use fileutils  , only : getfil
    use abortutils , only : endrun
    use ncdio_pio  , only : file_desc_t, ncd_defvar, ncd_io, ncd_pio_openfile, ncd_pio_closefile
    use spmdMod    , only : masterproc
    use topounit_varcon, only : max_topounits
    use GridcellType , only : grc_pp
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer            :: c,g,t,ti,topi          ! indices
    integer            :: mxsoil_color ! maximum number of soil color classes
    type(file_desc_t)  :: ncid         ! netcdf id
    character(len=256) :: locfn        ! local filename
    integer            :: ier          ! error status
    logical            :: readvar
    integer  ,pointer  :: soic2d (:,:)   ! read in - soil color
    !---------------------------------------------------------------------

    ! Allocate module variable for soil color

    allocate(isoicol(bounds%begc:bounds%endc))

    ! Determine soil color and number of soil color classes
    ! if number of soil color classes is not on input dataset set it to 8

    call getfil (fsurdat, locfn, 0)
    call ncd_pio_openfile (ncid, locfn, 0)

    call ncd_io(ncid=ncid, varname='mxsoil_color', flag='read', data=mxsoil_color, readvar=readvar)
    if ( .not. readvar ) mxsoil_color = 8

    allocate(soic2d(bounds%begg:bounds%endg,max_topounits))
    call ncd_io(ncid=ncid, varname='SOIL_COLOR', flag='read', data=soic2d, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=' ERROR: SOIL_COLOR NOT on surfdata file'//errMsg(__FILE__, __LINE__))
    end if
    do c = bounds%begc, bounds%endc
       g = col_pp%gridcell(c)
       t = col_pp%topounit(c)
       topi = grc_pp%topi(g)
       ti = t - topi + 1
       isoicol(c) = soic2d(g,ti)
       
    end do
    deallocate(soic2d)

    call ncd_pio_closefile(ncid)

    ! Determine saturated and dry soil albedos for n color classes and
    ! numrad wavebands (1=vis, 2=nir)

    allocate(albsat(mxsoil_color,numrad), albdry(mxsoil_color,numrad), stat=ier)
    if (ier /= 0) then
       write(iulog,*)'allocation error for albsat, albdry'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    if (masterproc) then
       write(iulog,*) 'Attempting to read soil colo data .....'
    end if

    if (mxsoil_color == 8) then
       albsat(1:8,1) = (/0.12_r8,0.11_r8,0.10_r8,0.09_r8,0.08_r8,0.07_r8,0.06_r8,0.05_r8/)
       albsat(1:8,2) = (/0.24_r8,0.22_r8,0.20_r8,0.18_r8,0.16_r8,0.14_r8,0.12_r8,0.10_r8/)
       albdry(1:8,1) = (/0.24_r8,0.22_r8,0.20_r8,0.18_r8,0.16_r8,0.14_r8,0.12_r8,0.10_r8/)
       albdry(1:8,2) = (/0.48_r8,0.44_r8,0.40_r8,0.36_r8,0.32_r8,0.28_r8,0.24_r8,0.20_r8/)
    else if (mxsoil_color == 20) then
       albsat(1:20,1) = (/0.25_r8,0.23_r8,0.21_r8,0.20_r8,0.19_r8,0.18_r8,0.17_r8,0.16_r8,&
            0.15_r8,0.14_r8,0.13_r8,0.12_r8,0.11_r8,0.10_r8,0.09_r8,0.08_r8,0.07_r8,0.06_r8,0.05_r8,0.04_r8/)
       albsat(1:20,2) = (/0.50_r8,0.46_r8,0.42_r8,0.40_r8,0.38_r8,0.36_r8,0.34_r8,0.32_r8,&
            0.30_r8,0.28_r8,0.26_r8,0.24_r8,0.22_r8,0.20_r8,0.18_r8,0.16_r8,0.14_r8,0.12_r8,0.10_r8,0.08_r8/)
       albdry(1:20,1) = (/0.36_r8,0.34_r8,0.32_r8,0.31_r8,0.30_r8,0.29_r8,0.28_r8,0.27_r8,&
            0.26_r8,0.25_r8,0.24_r8,0.23_r8,0.22_r8,0.20_r8,0.18_r8,0.16_r8,0.14_r8,0.12_r8,0.10_r8,0.08_r8/)
       albdry(1:20,2) = (/0.61_r8,0.57_r8,0.53_r8,0.51_r8,0.49_r8,0.48_r8,0.45_r8,0.43_r8,&
            0.41_r8,0.39_r8,0.37_r8,0.35_r8,0.33_r8,0.31_r8,0.29_r8,0.27_r8,0.25_r8,0.23_r8,0.21_r8,0.16_r8/)
    else
      write(iulog,*)'maximum color class = ',mxsoil_color,' is not supported'
      call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ! Set alblakwi
    alblakwi(:) = lake_melt_icealb(:)

  !$acc enter data copyin(albsat,albdry,isoicol,alblakwi )
  end subroutine SurfaceAlbedoInitTimeConst

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(surfalb_type) :: this
    type(bounds_type), intent(in) :: bounds

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    call this%InitCold(bounds)

  end subroutine Init


  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! Allocate module variables and data structures
    !
    ! !USES:
    use elm_varcon    , only: spval, ispval
    !
    ! !ARGUMENTS:
    class(surfalb_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc

    allocate(this%coszen_col         (begc:endc))              ; this%coszen_col         (:)   = spval
    allocate(this%albgrd_col         (begc:endc,numrad))       ; this%albgrd_col         (:,:) =spval
    allocate(this%albgri_col         (begc:endc,numrad))       ; this%albgri_col         (:,:) =spval
    allocate(this%albsnd_hst_col     (begc:endc,numrad))       ; this%albsnd_hst_col     (:,:) = spval
    allocate(this%albsni_hst_col     (begc:endc,numrad))       ; this%albsni_hst_col     (:,:) = spval
    allocate(this%albsod_col         (begc:endc,numrad))       ; this%albsod_col         (:,:) = spval
    allocate(this%albsoi_col         (begc:endc,numrad))       ; this%albsoi_col         (:,:) = spval
    allocate(this%albgrd_pur_col     (begc:endc,numrad))       ; this%albgrd_pur_col     (:,:) =spval
    allocate(this%albgri_pur_col     (begc:endc,numrad))       ; this%albgri_pur_col     (:,:) =spval
    allocate(this%albgrd_bc_col      (begc:endc,numrad))       ; this%albgrd_bc_col      (:,:) =spval
    allocate(this%albgri_bc_col      (begc:endc,numrad))       ; this%albgri_bc_col      (:,:) =spval
    allocate(this%albgrd_oc_col      (begc:endc,numrad))       ; this%albgrd_oc_col      (:,:) =spval
    allocate(this%albgri_oc_col      (begc:endc,numrad))       ; this%albgri_oc_col      (:,:) =spval
    allocate(this%albgrd_dst_col     (begc:endc,numrad))       ; this%albgrd_dst_col     (:,:) =spval
    allocate(this%albgri_dst_col     (begc:endc,numrad))       ; this%albgri_dst_col     (:,:) =spval
    allocate(this%albd_patch         (begp:endp,numrad))       ; this%albd_patch         (:,:) =spval
    allocate(this%albi_patch         (begp:endp,numrad))       ; this%albi_patch         (:,:) =spval

    allocate(this%ftdd_patch         (begp:endp,numrad))       ; this%ftdd_patch         (:,:) =spval
    allocate(this%ftid_patch         (begp:endp,numrad))       ; this%ftid_patch         (:,:) =spval
    allocate(this%ftii_patch         (begp:endp,numrad))       ; this%ftii_patch         (:,:) =spval
    allocate(this%fabd_patch         (begp:endp,numrad))       ; this%fabd_patch         (:,:) =spval
    allocate(this%fabd_sun_patch     (begp:endp,numrad))       ; this%fabd_sun_patch     (:,:) =spval
    allocate(this%fabd_sha_patch     (begp:endp,numrad))       ; this%fabd_sha_patch     (:,:) =spval
    allocate(this%fabi_patch         (begp:endp,numrad))       ; this%fabi_patch         (:,:) =spval
    allocate(this%fabi_sun_patch     (begp:endp,numrad))       ; this%fabi_sun_patch     (:,:) =spval
    allocate(this%fabi_sha_patch     (begp:endp,numrad))       ; this%fabi_sha_patch     (:,:) =spval
    allocate(this%fabd_sun_z_patch   (begp:endp,nlevcan))      ; this%fabd_sun_z_patch   (:,:) = 0._r8
    allocate(this%fabd_sha_z_patch   (begp:endp,nlevcan))      ; this%fabd_sha_z_patch   (:,:) = 0._r8
    allocate(this%fabi_sun_z_patch   (begp:endp,nlevcan))      ; this%fabi_sun_z_patch   (:,:) = 0._r8
    allocate(this%fabi_sha_z_patch   (begp:endp,nlevcan))      ; this%fabi_sha_z_patch   (:,:) = 0._r8
    allocate(this%flx_absdv_col      (begc:endc,-nlevsno+1:1)) ; this%flx_absdv_col      (:,:) = spval
    allocate(this%flx_absdn_col      (begc:endc,-nlevsno+1:1)) ; this%flx_absdn_col      (:,:) = spval
    allocate(this%flx_absiv_col      (begc:endc,-nlevsno+1:1)) ; this%flx_absiv_col      (:,:) = spval
    allocate(this%flx_absin_col      (begc:endc,-nlevsno+1:1)) ; this%flx_absin_col      (:,:) = spval

    allocate(this%fsun_z_patch       (begp:endp,nlevcan))      ; this%fsun_z_patch       (:,:) = 0._r8
    allocate(this%tlai_z_patch       (begp:endp,nlevcan))      ; this%tlai_z_patch       (:,:) = 0._r8
    allocate(this%tsai_z_patch       (begp:endp,nlevcan))      ; this%tsai_z_patch       (:,:) = 0._r8
    allocate(this%ncan_patch         (begp:endp))              ; this%ncan_patch         (:)   = 0
    allocate(this%nrad_patch         (begp:endp))              ; this%nrad_patch         (:)   = 0
    allocate(this%vcmaxcintsun_patch (begp:endp))              ; this%vcmaxcintsun_patch (:)   =spval
    allocate(this%vcmaxcintsha_patch (begp:endp))              ; this%vcmaxcintsha_patch (:)   =spval

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! History fields initialization
    !
    ! !USES:
    use shr_infnan_mod, only: nan => shr_infnan_nan, assignment(=)
    use elm_varcon    , only: spval
    use histFileMod   , only: hist_addfld1d, hist_addfld2d
    !
    ! !ARGUMENTS:
    class(surfalb_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc

    this%coszen_col(begc:endc) = spval
    call hist_addfld1d (fname='COSZEN', units='1', &
         avgflag='A', long_name='cosine of solar zenith angle', &
         ptr_col=this%coszen_col, default='inactive')

    this%albgri_col(begc:endc,:) = spval
    call hist_addfld2d (fname='ALBGRD', units='proportion', type2d='numrad', &
         avgflag='A', long_name='ground albedo (direct)', &
         ptr_col=this%albgrd_col, default='inactive')

    this%albgri_col(begc:endc,:) = spval
    call hist_addfld2d (fname='ALBGRI', units='proportion', type2d='numrad', &
         avgflag='A', long_name='ground albedo (indirect)', &
         ptr_col=this%albgri_col, default='inactive')

    this%albd_patch(begp:endp,:) = spval
    call hist_addfld2d (fname='ALBD', units='proportion', type2d='numrad', &
         avgflag='A', long_name='surface albedo (direct)', &
         ptr_patch=this%albd_patch, default='inactive', c2l_scale_type='urbanf')

    this%albi_patch(begp:endp,:) = spval
    call hist_addfld2d (fname='ALBI', units='proportion', type2d='numrad', &
         avgflag='A', long_name='surface albedo (indirect)', &
         ptr_patch=this%albi_patch, default='inactive', c2l_scale_type='urbanf')

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! Initialize module surface albedos to reasonable values
    !
    ! !ARGUMENTS:
    class(surfalb_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc
    integer :: begp, endp
    !-----------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    this%albgrd_col     (begc:endc, :) = 0.2_r8
    this%albgri_col     (begc:endc, :) = 0.2_r8
    this%albsod_col     (begc:endc, :) = 0.2_r8
    this%albsoi_col     (begc:endc, :) = 0.2_r8
    this%albsnd_hst_col (begc:endc, :) = 0.6_r8
    this%albsni_hst_col (begc:endc, :) = 0.6_r8
    this%albd_patch     (begp:endp, :) = 0.2_r8
    this%albi_patch     (begp:endp, :) = 0.2_r8

    this%albgrd_pur_col (begc:endc, :) = 0.2_r8
    this%albgri_pur_col (begc:endc, :) = 0.2_r8
    this%albgrd_bc_col  (begc:endc, :) = 0.2_r8
    this%albgri_bc_col  (begc:endc, :) = 0.2_r8
    this%albgrd_oc_col  (begc:endc, :) = 0.2_r8
    this%albgri_oc_col  (begc:endc, :) = 0.2_r8
    this%albgrd_dst_col (begc:endc, :) = 0.2_r8
    this%albgri_dst_col (begc:endc, :) = 0.2_r8

    this%fabi_patch     (begp:endp, :) = 0.0_r8
    this%fabd_patch     (begp:endp, :) = 0.0_r8
    this%fabi_sun_patch (begp:endp, :) = 0.0_r8
    this%fabd_sun_patch (begp:endp, :) = 0.0_r8
    this%fabd_sha_patch (begp:endp, :) = 0.0_r8
    this%fabi_sha_patch (begp:endp, :) = 0.0_r8
    this%ftdd_patch     (begp:endp, :) = 1.0_r8
    this%ftid_patch     (begp:endp, :) = 0.0_r8
    this%ftii_patch     (begp:endp, :) = 1.0_r8

  end subroutine InitCold

  !---------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag, &
       tlai_patch, tsai_patch)
    !
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use elm_varctl , only : use_snicar_frc, iulog
    use spmdMod    , only : masterproc
    use decompMod  , only : bounds_type
    use abortutils , only : endrun
    use ncdio_pio  , only : file_desc_t, ncd_defvar, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(surfalb_type)               :: this
    type(bounds_type) , intent(in)    :: bounds
    type(file_desc_t) , intent(inout) :: ncid ! netcdf id
    character(len=*)  , intent(in)    :: flag ! 'read' or 'write'
    real(r8)          , intent(in)    :: tlai_patch(bounds%begp:)
    real(r8)          , intent(in)    :: tsai_patch(bounds%begp:)
    !
    ! !LOCAL VARIABLES:
    logical :: readvar      ! determine if variable is on initial file
    integer :: iv
    integer :: begp, endp
    integer :: begc, endc
    !---------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(tlai_patch)  == (/bounds%endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(tsai_patch)  == (/bounds%endp/)), errMsg(__FILE__, __LINE__))

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc

    call restartvar(ncid=ncid, flag=flag, varname='coszen', xtype=ncd_double,  &
         dim1name='column', &
         long_name='cosine of solar zenith angle', units='1', &
         interpinic_flag='interp', readvar=readvar, data=this%coszen_col)

    call restartvar(ncid=ncid, flag=flag, varname='albd', xtype=ncd_double,  &
         dim1name='pft', dim2name='numrad', switchdim=.true., &
         long_name='surface albedo (direct) (0 to 1)', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%albd_patch)

    call restartvar(ncid=ncid, flag=flag, varname='albi', xtype=ncd_double,  &
         dim1name='pft', dim2name='numrad', switchdim=.true., &
         long_name='surface albedo (diffuse) (0 to 1)', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%albi_patch)

    call restartvar(ncid=ncid, flag=flag, varname='albgrd', xtype=ncd_double,  &
         dim1name='column', dim2name='numrad', switchdim=.true., &
         long_name='ground albedo (direct) (0 to 1)', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%albgrd_col)

    call restartvar(ncid=ncid, flag=flag, varname='albgri', xtype=ncd_double,  &
         dim1name='column', dim2name='numrad', switchdim=.true., &
         long_name='ground albedo (indirect) (0 to 1)', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%albgri_col)

    call restartvar(ncid=ncid, flag=flag, varname='albsod', xtype=ncd_double,  &
         dim1name='column', dim2name='numrad', switchdim=.true., &
         long_name='soil albedo (direct) (0 to 1)', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%albsod_col)

    call restartvar(ncid=ncid, flag=flag, varname='albsoi', xtype=ncd_double,  &
         dim1name='column', dim2name='numrad', switchdim=.true., &
         long_name='soil albedo (indirect) (0 to 1)', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%albsoi_col)

    call restartvar(ncid=ncid, flag=flag, varname='albsnd_hst', xtype=ncd_double,  &
         dim1name='column', dim2name='numrad', switchdim=.true., &
         long_name='snow albedo (direct) (0 to 1)', units='proportion', &
         interpinic_flag='interp', readvar=readvar, data=this%albsnd_hst_col)

    call restartvar(ncid=ncid, flag=flag, varname='albsni_hst', xtype=ncd_double,  &
         dim1name='column', dim2name='numrad', switchdim=.true., &
         long_name='snow albedo (diffuse) (0 to 1)', units='proportion', &
         interpinic_flag='interp', readvar=readvar, data=this%albsni_hst_col)

    call restartvar(ncid=ncid, flag=flag, varname='tlai_z', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='tlai increment for canopy layer', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%tlai_z_patch)
    if (flag=='read' .and. .not. readvar) then
       if (masterproc) then
          write(iulog,*) "can't find tlai_z in restart (or initial) file..."
          write(iulog,*) "Initialize tlai_z to tlai/nlevcan"
       end if
       do iv=1,nlevcan
          this%tlai_z_patch(begp:endp,iv) =  tlai_patch(begp:endp) / nlevcan
       end do
    end if

    call restartvar(ncid=ncid, flag=flag, varname='tsai_z', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='tsai increment for canopy layer', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%tsai_z_patch)
    if (flag=='read' .and. .not. readvar) then
       if (masterproc) then
          write(iulog,*) "can't find tsai_z in restart (or initial) file..."
          write(iulog,*) "Initialize tsai_z to tsai/nlevcan"
       end if
       do iv=1,nlevcan
          this%tsai_z_patch(begp:endp,iv) = tsai_patch(begp:endp) / nlevcan
       end do
    end if

    call restartvar(ncid=ncid, flag=flag, varname='ncan', xtype=ncd_int,  &
         dim1name='pft', long_name='number of canopy layers', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%ncan_patch)
    if (flag=='read' .and. .not. readvar) then
       if (masterproc) write(iulog,*) "can't find ncan in restart (or initial) file..."
       if (masterproc) write(iulog,*) "Initialize ncan to nlevcan"
       this%ncan_patch(begp:endp) = nlevcan
    end if

    call restartvar(ncid=ncid, flag=flag, varname='nrad', xtype=ncd_int,  &
         dim1name='pft', long_name='number of canopy layers, above snow for radiative transfer', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%nrad_patch)
    if (flag=='read' .and. .not. readvar) then
       if (masterproc) write(iulog,*) "can't find nrad in restart (or initial) file..."
       if (masterproc) write(iulog,*) "Initialize nrad to nlevcan"
       this%nrad_patch(begp:endp) = nlevcan
    end if

    call restartvar(ncid=ncid, flag=flag, varname='fsun_z', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='sunlit fraction for canopy layer', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%fsun_z_patch)
    if (flag=='read' .and. .not. readvar) then
       if (masterproc) write(iulog,*) "can't find fsun_z in restart (or initial) file..."
       if (masterproc) write(iulog,*) "Initialize fsun_z to 0"
       do iv=1,nlevcan
          this%fsun_z_patch(begp:endp,iv) = 0._r8
       end do
    end if

    call restartvar(ncid=ncid, flag=flag, varname='vcmaxcintsun', xtype=ncd_double, &
         dim1name='pft', long_name='sunlit canopy scaling coefficient', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%vcmaxcintsun_patch)
    if (flag=='read' .and. .not. readvar) then
       if (masterproc) write(iulog,*) "can't find vcmaxcintsun in restart (or initial) file..."
       if (masterproc) write(iulog,*) "Initialize vcmaxcintsun to 1"
       this%vcmaxcintsun_patch(begp:endp) = 1._r8
    end if

    call restartvar(ncid=ncid, flag=flag, varname='vcmaxcintsha', xtype=ncd_double,  &
         dim1name='pft', long_name='shaded canopy scaling coefficient', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%vcmaxcintsha_patch)
    if (flag=='read' .and. .not. readvar) then
       if (masterproc) write(iulog,*) "can't find vcmaxcintsha in restart (or initial) file..."
       if (masterproc) write(iulog,*) "Initialize vcmaxcintsha to 1"
       this%vcmaxcintsha_patch(begp:endp) = 1._r8
    end if

    if (use_snicar_frc) then

       call restartvar(ncid=ncid, flag=flag, varname='albgrd_bc', xtype=ncd_double,  &
            dim1name='column', dim2name='numrad', switchdim=.true., &
            long_name='ground albedo without BC (direct) (0 to 1)', units='', &
            interpinic_flag='interp',readvar=readvar, data=this%albgrd_bc_col)
       if (flag=='read' .and. .not. readvar) then
          if (masterproc) write(iulog,*) "SNICAR: can't find albgrd_bc in initial file..."
          if (masterproc) write(iulog,*) "Initialize albgrd_bc to albgrd"
          this%albgrd_bc_col(begc:endc,:) = this%albgrd_col(begc:endc,:)
       end if

       call restartvar(ncid=ncid, flag=flag, varname='albgri_bc', xtype=ncd_double,  &
            dim1name='column', dim2name='numrad', switchdim=.true., &
            long_name='ground albedo without BC (diffuse) (0 to 1)', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%albgri_bc_col)
       if (flag=='read' .and. .not. readvar) then
          if (masterproc) write(iulog,*) "SNICAR: can't find albgri_bc in initial file..."
          if (masterproc) write(iulog,*) "Initialize albgri_bc to albgri"
          this%albgri_bc_col(begc:endc,:) = this%albgri_col(begc:endc,:)
       end if

       call restartvar(ncid=ncid, flag=flag, varname='albgrd_pur', xtype=ncd_double,  &
            dim1name='column', dim2name='numrad', switchdim=.true., &
            long_name='pure snow ground albedo (direct) (0 to 1)', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%albgrd_pur_col)
       if (flag=='read' .and. .not. readvar) then
          if (masterproc) write(iulog,*) "SNICAR: can't find albgrd_pur in initial file..."
          if (masterproc) write(iulog,*) "Initialize albgrd_pur to albgrd"
          this%albgrd_pur_col(begc:endc,:) = this%albgrd_col(begc:endc,:)
       end if

       call restartvar(ncid=ncid, flag=flag, varname='albgri_pur', xtype=ncd_double,  &
            dim1name='column', dim2name='numrad', switchdim=.true., &
            long_name='pure snow ground albedo (diffuse) (0 to 1)', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%albgri_pur_col)
       if (flag=='read' .and. .not. readvar) then
          if (masterproc) write(iulog,*) "SNICAR: can't find albgri_pur in initial file..."
          if (masterproc) write(iulog,*) "Initialize albgri_pur to albgri"
          this%albgri_pur_col(begc:endc,:) = this%albgri_col(begc:endc,:)
       end if

       call restartvar(ncid=ncid, flag=flag, varname='albgrd_oc', xtype=ncd_double,  &
            dim1name='column', dim2name='numrad', switchdim=.true., &
            long_name='ground albedo without OC (direct) (0 to 1)', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%albgrd_oc_col)
       if (flag=='read' .and. .not. readvar) then
          if (masterproc) write(iulog,*) "SNICAR: can't find albgrd_oc in initial file..."
          if (masterproc) write(iulog,*) "Initialize albgrd_oc to albgrd"
          this%albgrd_oc_col(begc:endc,:) = this%albgrd_col(begc:endc,:)
       end if

       call restartvar(ncid=ncid, flag=flag, varname='albgri_oc', xtype=ncd_double,  &
            dim1name='column', dim2name='numrad', switchdim=.true., &
            long_name='ground albedo without OC (diffuse) (0 to 1)', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%albgri_oc_col)
       if (flag=='read' .and. .not. readvar) then
          if (masterproc) write(iulog,*) "SNICAR: can't find albgri_oc in restart (or initial) file..."
          if (masterproc) write(iulog,*) "Initialize albgri_oc to albgri"
          this%albgri_oc_col(begc:endc,:) = this%albgri_col(begc:endc,:)
       end if

       call restartvar(ncid=ncid, flag=flag, varname='albgrd_dst', xtype=ncd_double,  &
            dim1name='column', dim2name='numrad', switchdim=.true., &
            long_name='ground albedo without dust (direct) (0 to 1)', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%albgrd_dst_col)
       if (flag=='read' .and. .not. readvar) then
          if (masterproc) write(iulog,*) "SNICAR: can't find albgrd_dst in initial file..."
          if (masterproc) write(iulog,*) "Initialize albgrd_dst to albgrd"
          this%albgrd_dst_col(begc:endc,:) = this%albgrd_col(begc:endc,:)
       end if

       call restartvar(ncid=ncid, flag=flag, varname='albgri_dst', xtype=ncd_double,  &
            dim1name='column', dim2name='numrad', switchdim=.true., &
            long_name='ground albedo without dust (diffuse) (0 to 1)', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%albgri_dst_col)
       if (flag=='read' .and. .not. readvar) then
          if (masterproc) write(iulog,*) "SNICAR: can't find albgri_dst in initial file..."
          if (masterproc) write(iulog,*) "Initialize albgri_dst to albgri"
          this%albgri_dst_col(begc:endc,:) = this%albgri_col(begc:endc,:)
       end if

    end if  ! end of if-use_snicar_frc

    ! patch type physical state variable - fabd
    call restartvar(ncid=ncid, flag=flag, varname='fabd', xtype=ncd_double,  &
         dim1name='pft', dim2name='numrad', switchdim=.true., &
         long_name='flux absorbed by veg per unit direct flux', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%fabd_patch)

    call restartvar(ncid=ncid, flag=flag, varname='fabi', xtype=ncd_double,  &
         dim1name='pft', dim2name='numrad', switchdim=.true., &
         long_name='flux absorbed by veg per unit diffuse flux', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%fabi_patch)

    call restartvar(ncid=ncid, flag=flag, varname='fabd_sun', xtype=ncd_double,  &
         dim1name='pft', dim2name='numrad', switchdim=.true., &
         long_name='flux absorbed by sunlit leaf per unit direct flux', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%fabd_sun_patch)
    if (flag=='read' .and. .not. readvar) then
       if (masterproc) write(iulog,*) "can't find fabd_sun in restart (or initial) file..."
       if (masterproc) write(iulog,*) "Initialize fabd_sun to fabd/2"
       this%fabd_sun_patch(begp:endp,:) = this%fabd_patch(begp:endp,:)/2._r8
    end if

    call restartvar(ncid=ncid, flag=flag, varname='fabd_sha', xtype=ncd_double,  &
         dim1name='pft', dim2name='numrad', switchdim=.true., &
         long_name='flux absorbed by shaded leaf per unit direct flux', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%fabd_sha_patch)
    if (flag=='read' .and. .not. readvar) then
       if (masterproc) write(iulog,*) "can't find fabd_sha in restart (or initial) file..."
       if (masterproc) write(iulog,*) "Initialize fabd_sha to fabd/2"
       this%fabd_sha_patch(begp:endp,:) = this%fabd_patch(begp:endp,:)/2._r8
    end if

    call restartvar(ncid=ncid, flag=flag, varname='fabi_sun', xtype=ncd_double,  &
         dim1name='pft', dim2name='numrad', switchdim=.true., &
         long_name='flux absorbed by sunlit leaf per unit diffuse flux', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%fabi_sun_patch)
    if (flag=='read' .and. .not. readvar) then
       if (masterproc) write(iulog,*) "can't find fabi_sun in restart (or initial) file..."
       if (masterproc) write(iulog,*) "Initialize fabi_sun to fabi/2"
       this%fabi_sun_patch(begp:endp,:) = this%fabi_patch(begp:endp,:)/2._r8
    end if

    call restartvar(ncid=ncid, flag=flag, varname='fabi_sha', xtype=ncd_double,  &
         dim1name='pft', dim2name='numrad', switchdim=.true., &
         long_name='flux absorbed by shaded leaf per unit diffuse flux', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%fabi_sha_patch)
    if (flag=='read' .and. .not. readvar) then
       if (masterproc) write(iulog,*) "can't find fabi_sha in restart (or initial) file..."
       if (masterproc) write(iulog,*) "Initialize fabi_sha to fabi/2"
       this%fabi_sha_patch(begp:endp,:) = this%fabi_patch(begp:endp,:)/2._r8
    end if

    call restartvar(ncid=ncid, flag=flag, varname='fabd_sun_z', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='absorbed sunlit leaf direct PAR (per unit lai+sai) for canopy layer', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%fabd_sun_z_patch)
    if (flag=='read' .and. .not. readvar) then
       if (masterproc) write(iulog,*) "can't find fabd_sun_z in restart (or initial) file..."
       if (masterproc) write(iulog,*) "Initialize fabd_sun_z to (fabd/2)/nlevcan"
       do iv=1,nlevcan
          this%fabd_sun_z_patch(begp:endp,iv) = (this%fabd_patch(begp:endp,1)/2._r8)/nlevcan
       end do
    end if

    call restartvar(ncid=ncid, flag=flag, varname='fabd_sha_z', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='absorbed shaded leaf direct PAR (per unit lai+sai) for canopy layer', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%fabd_sha_z_patch)
    if (flag=='read' .and. .not. readvar) then
       if (masterproc) write(iulog,*) "can't find fabd_sha_z in restart (or initial) file..."
       if (masterproc) write(iulog,*) "Initialize fabd_sha_z to (fabd/2)/nlevcan"
       do iv=1,nlevcan
          this%fabd_sha_z_patch(begp:endp,iv) = (this%fabd_patch(begp:endp,1)/2._r8)/nlevcan
       end do
    end if

    call restartvar(ncid=ncid, flag=flag, varname='fabi_sun_z', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='absorbed sunlit leaf diffuse PAR (per unit lai+sai) for canopy layer', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%fabi_sun_z_patch)
    if (flag=='read' .and. .not. readvar) then
       if (masterproc) write(iulog,*) "can't find fabi_sun_z in restart (or initial) file..."
       if (masterproc) write(iulog,*) "Initialize fabi_sun_z to (fabi/2)/nlevcan"
       do iv=1,nlevcan
          this%fabi_sun_z_patch(begp:endp,iv) = (this%fabi_patch(begp:endp,1)/2._r8)/nlevcan
       end do
    end if

    call restartvar(ncid=ncid, flag=flag, varname='fabi_sha_z', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='absorbed shaded leaf diffuse PAR (per unit lai+sai) for canopy layer', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%fabi_sha_z_patch)
    if (flag=='read' .and. .not. readvar) then
       if (masterproc) write(iulog,*) "can't find fabi_sha_z in restart (or initial) file..."
       if (masterproc) write(iulog,*) "Initialize fabi_sha_z to (fabi/2)/nlevcan"
       do iv=1,nlevcan
          this%fabi_sha_z_patch(begp:endp,iv) = &
               (this%fabi_patch(begp:endp,1)/2._r8)/nlevcan
       end do
    end if

    call restartvar(ncid=ncid, flag=flag, varname='ftdd', xtype=ncd_double,  &
         dim1name='pft', dim2name='numrad', switchdim=.true., &
         long_name='down direct flux below veg per unit direct flux', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%ftdd_patch)

    call restartvar(ncid=ncid, flag=flag, varname='ftid', xtype=ncd_double,  &
         dim1name='pft', dim2name='numrad', switchdim=.true., &
         long_name='down diffuse flux below veg per unit direct flux', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%ftid_patch)

    call restartvar(ncid=ncid, flag=flag, varname='ftii', xtype=ncd_double,  &
         dim1name='pft', dim2name='numrad', switchdim=.true., &
         long_name='down diffuse flux below veg per unit diffuse flux', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%ftii_patch)

    !--------------------------------
    ! variables needed for SNICAR
    !--------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='flx_absdv', xtype=ncd_double,  &
         dim1name='column', dim2name='levsno1', switchdim=.true., lowerb2=-nlevsno+1, upperb2=1, &
         long_name='snow layer flux absorption factors (direct, VIS)', units='fraction', &
         interpinic_flag='interp', readvar=readvar, data=this%flx_absdv_col)

    call restartvar(ncid=ncid, flag=flag, varname='flx_absdn', xtype=ncd_double,  &
         dim1name='column', dim2name='levsno1', switchdim=.true., lowerb2=-nlevsno+1, upperb2=1, &
         long_name='snow layer flux absorption factors (direct, NIR)', units='fraction', &
         interpinic_flag='interp', readvar=readvar, data=this%flx_absdn_col)

    call restartvar(ncid=ncid, flag=flag, varname='flx_absiv', xtype=ncd_double,  &
         dim1name='column', dim2name='levsno1', switchdim=.true., lowerb2=-nlevsno+1, upperb2=1, &
         long_name='snow layer flux absorption factors (diffuse, VIS)', units='fraction', &
         interpinic_flag='interp', readvar=readvar, data=this%flx_absiv_col)

    call restartvar(ncid=ncid, flag=flag, varname='flx_absin', xtype=ncd_double,  &
         dim1name='column', dim2name='levsno1', switchdim=.true., lowerb2=-nlevsno+1, upperb2=1, &
         long_name='snow layer flux absorption factors (diffuse, NIR)', units='fraction', &
         interpinic_flag='interp', readvar=readvar, data=this%flx_absin_col)

  end subroutine Restart


end module SurfaceAlbedoType

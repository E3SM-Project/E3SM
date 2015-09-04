module FrictionVelocityType

  !------------------------------------------------------------------------------
  !
  ! !USES:
  use shr_sys_mod    , only : shr_sys_flush
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use clm_varctl     , only : use_cn
  use clm_varpar     , only : nlevcan, nlevsno, nlevgrnd, nlevsoi  
  use clm_varcon     , only : spval
  use decompMod      , only : bounds_type
  use abortutils     , only : endrun
  use spmdMod        , only : masterproc
  use LandunitType   , only : lun                
  use ColumnType     , only : col                
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  type, public :: frictionvel_type

     ! Roughness length/resistance for friction velocity calculation

     real(r8), pointer :: forc_hgt_u_patch (:)   ! patch wind forcing height (10m+z0m+d) (m)
     real(r8), pointer :: forc_hgt_t_patch (:)   ! patch temperature forcing height (10m+z0m+d) (m)
     real(r8), pointer :: forc_hgt_q_patch (:)   ! patch specific humidity forcing height (10m+z0m+d) (m)
     real(r8), pointer :: u10_patch        (:)   ! patch 10-m wind (m/s) (for dust model)
     real(r8), pointer :: u10_clm_patch    (:)   ! patch 10-m wind (m/s) (for clm_map2gcell)
     real(r8), pointer :: va_patch         (:)   ! patch atmospheric wind speed plus convective velocity (m/s)
     real(r8), pointer :: vds_patch        (:)   ! patch deposition velocity term (m/s) (for dry dep SO4, NH4NO3)
     real(r8), pointer :: fv_patch         (:)   ! patch friction velocity (m/s) (for dust model)
     real(r8), pointer :: rb1_patch        (:)   ! patch aerodynamical resistance (s/m) (for dry deposition of chemical tracers)
     real(r8), pointer :: ram1_patch       (:)   ! patch aerodynamical resistance (s/m)
     real(r8), pointer :: z0m_patch        (:)   ! patch momentum roughness length (m)
     real(r8), pointer :: z0mv_patch       (:)   ! patch roughness length over vegetation, momentum [m]
     real(r8), pointer :: z0hv_patch       (:)   ! patch roughness length over vegetation, sensible heat [m]
     real(r8), pointer :: z0qv_patch       (:)   ! patch roughness length over vegetation, latent heat [m]
     real(r8), pointer :: z0mg_col         (:)   ! col roughness length over ground, momentum  [m] 
     real(r8), pointer :: z0hg_col         (:)   ! col roughness length over ground, sensible heat [m]
     real(r8), pointer :: z0qg_col         (:)   ! col roughness length over ground, latent heat [m]

   contains

     procedure, public  :: Init         
     procedure, public  :: Restart      
     procedure, private :: InitAllocate 
     procedure, private :: InitHistory  
     procedure, private :: InitCold     

  end type frictionvel_type
  !------------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(frictionvel_type) :: this
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
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(frictionvel_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    allocate(this%forc_hgt_u_patch (begp:endp)) ; this%forc_hgt_u_patch (:)   = nan
    allocate(this%forc_hgt_t_patch (begp:endp)) ; this%forc_hgt_t_patch (:)   = nan
    allocate(this%forc_hgt_q_patch (begp:endp)) ; this%forc_hgt_q_patch (:)   = nan
    allocate(this%u10_patch        (begp:endp)) ; this%u10_patch        (:)   = nan
    allocate(this%u10_clm_patch    (begp:endp)) ; this%u10_clm_patch    (:)   = nan
    allocate(this%va_patch         (begp:endp)) ; this%va_patch         (:)   = nan
    allocate(this%vds_patch        (begp:endp)) ; this%vds_patch        (:)   = nan
    allocate(this%fv_patch         (begp:endp)) ; this%fv_patch         (:)   = nan
    allocate(this%rb1_patch        (begp:endp)) ; this%rb1_patch        (:)   = nan
    allocate(this%ram1_patch       (begp:endp)) ; this%ram1_patch       (:)   = nan
    allocate(this%z0m_patch        (begp:endp)) ; this%z0m_patch        (:)   = nan
    allocate(this%z0mv_patch       (begp:endp)) ; this%z0mv_patch       (:)   = nan
    allocate(this%z0hv_patch       (begp:endp)) ; this%z0hv_patch       (:)   = nan
    allocate(this%z0qv_patch       (begp:endp)) ; this%z0qv_patch       (:)   = nan
    allocate(this%z0mg_col         (begc:endc)) ; this%z0mg_col         (:)   = nan
    allocate(this%z0qg_col         (begc:endc)) ; this%z0qg_col         (:)   = nan
    allocate(this%z0hg_col         (begc:endc)) ; this%z0hg_col         (:)   = nan

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! History fields initialization
    !
    ! !USES:
    use shr_infnan_mod, only: nan => shr_infnan_nan, assignment(=)
    use histFileMod   , only: hist_addfld1d, hist_addfld2d
    !
    ! !ARGUMENTS:
    class(frictionvel_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc
    integer :: begp, endp
    !---------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    this%z0mg_col(begc:endc) = spval
    call hist_addfld1d (fname='Z0MG', units='m', &
         avgflag='A', long_name='roughness length over ground, momentum', &
         ptr_col=this%z0mg_col, default='inactive')

    this%z0hg_col(begc:endc) = spval
    call hist_addfld1d (fname='Z0HG', units='m', &
         avgflag='A', long_name='roughness length over ground, sensible heat', &
         ptr_col=this%z0hg_col, default='inactive')

    this%z0qg_col(begc:endc) = spval
    call hist_addfld1d (fname='Z0QG', units='m', &
         avgflag='A', long_name='roughness length over ground, latent heat', &
         ptr_col=this%z0qg_col, default='inactive')

    this%va_patch(begp:endp) = spval
    call hist_addfld1d (fname='VA', units='m/s', &
         avgflag='A', long_name='atmospheric wind speed plus convective velocity', &
         ptr_patch=this%va_patch, default='inactive')

    this%u10_clm_patch(begp:endp) = spval
    call hist_addfld1d (fname='U10', units='m/s', &
         avgflag='A', long_name='10-m wind', &
         ptr_patch=this%u10_clm_patch)

    if (use_cn) then
       this%u10_patch(begp:endp) = spval
       call hist_addfld1d (fname='U10_DUST', units='m/s', &
            avgflag='A', long_name='10-m wind for dust model', &
            ptr_patch=this%u10_patch, default='inactive')
    end if

    if (use_cn) then
       this%ram1_patch(begp:endp) = spval
       call hist_addfld1d (fname='RAM1', units='s/m', &
            avgflag='A', long_name='aerodynamical resistance ', &
            ptr_patch=this%ram1_patch, default='inactive')
    end if

    if (use_cn) then
       this%fv_patch(begp:endp) = spval
       call hist_addfld1d (fname='FV', units='m/s', &
            avgflag='A', long_name='friction velocity for dust model', &
            ptr_patch=this%fv_patch, default='inactive')
    end if

    if (use_cn) then
       this%z0hv_patch(begp:endp) = spval
       call hist_addfld1d (fname='Z0HV', units='m', &
            avgflag='A', long_name='roughness length over vegetation, sensible heat', &
            ptr_patch=this%z0hv_patch, default='inactive')
    end if

    if (use_cn) then
       this%z0m_patch(begp:endp) = spval
       call hist_addfld1d (fname='Z0M', units='m', &
            avgflag='A', long_name='momentum roughness length', &
            ptr_patch=this%z0m_patch, default='inactive')
    end if

    if (use_cn) then
       this%z0mv_patch(begp:endp) = spval
       call hist_addfld1d (fname='Z0MV', units='m', &
            avgflag='A', long_name='roughness length over vegetation, momentum', &
            ptr_patch=this%z0mv_patch, default='inactive')
    end if

    if (use_cn) then
       this%z0qv_patch(begp:endp) = spval
       call hist_addfld1d (fname='Z0QV', units='m', &
            avgflag='A', long_name='roughness length over vegetation, latent heat', &
            ptr_patch=this%z0qv_patch, default='inactive')
    end if

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! Initialize module surface albedos to reasonable values
    !
    ! !ARGUMENTS:
    class(frictionvel_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: p, c, l                         ! indices
    !-----------------------------------------------------------------------

    ! Added 5/4/04, PET: initialize forc_hgt_u (gridcell-level),
    ! since this is not initialized before first call to CNVegStructUpdate,
    ! and it is required to set the upper bound for canopy top height.
    ! Changed 3/21/08, KO: still needed but don't have sufficient information 
    ! to set this properly (e.g., patch-level displacement height and roughness 
    ! length). So leave at 30m.

    if (use_cn) then
       do p = bounds%begp, bounds%endp
          this%forc_hgt_u_patch(p) = 30._r8
       end do
    end if

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%lakpoi(l)) then !lake
          this%z0mg_col(c) = 0.0004_r8
       end if
    end do
   
  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use spmdMod    , only : masterproc
    use ncdio_pio  , only : file_desc_t, ncd_defvar, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(frictionvel_type) :: this
    type(bounds_type) , intent(in)    :: bounds 
    type(file_desc_t) , intent(inout) :: ncid   ! netcdf id
    character(len=*)  , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer :: j,c ! indices
    logical :: readvar      ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='Z0MG', xtype=ncd_double,  &
         dim1name='column', &
         long_name='ground momentum roughness length', units='m', &
         interpinic_flag='interp', readvar=readvar, data=this%z0mg_col)

  end subroutine Restart


end module FrictionVelocityType


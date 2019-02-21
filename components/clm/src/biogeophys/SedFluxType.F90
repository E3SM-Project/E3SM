module SedFluxType

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Hold sediment, POC, PON and POP dynamic fluxes induced by soil erosion 
  !
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use clm_varcon        , only : spval
  use clm_varctl        , only : use_erosion 
  use decompMod         , only : bounds_type
  use abortutils        , only : endrun
  use ColumnType        , only : col_pp
  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  type, public :: sedflux_type

     ! Soil erosion model parameters
     real(r8), pointer :: prnsf_col(:)              ! col rainfall-driven erosion scaling factor
     real(r8), pointer :: qrnsf_col(:)              ! col runoff-driven erosion scaling factor
     real(r8), pointer :: tcsf_col(:)               ! col transport capacity scaling factor
     
     ! Soil erosion fluxes
     real(r8), pointer :: flx_p_ero_col(:)          ! col sed detach driven by rainfall (kg/m2/s)
     real(r8), pointer :: flx_q_ero_col(:)          ! col sed detach driven by runoff (kg/m2/s)
     real(r8), pointer :: flx_lndsld_ero_col(:)     ! col sed detach driven by landslides (kg/m2/s)
     real(r8), pointer :: flx_sed_ero_col(:)        ! col total sed detach (kg/m2/s)
     real(r8), pointer :: flx_lndsld_yld_col(:)     ! col sed yield driven by landslides (kg/m2/s)
     real(r8), pointer :: flx_sed_yld_col(:)        ! col total sed yield (kg/m2/s) 

  contains

     procedure, public  :: Init
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold
     procedure, public  :: Restart

  end type sedflux_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)
    
    class(sedflux_type)            :: this
    type(bounds_type) , intent(in) :: bounds

    call this%InitAllocate ( bounds )
    if (use_erosion) then
       call this%InitHistory ( bounds )
       call this%InitCold ( bounds )
    end if

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
    class(sedflux_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc
    !------------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc

    allocate( this%flx_p_ero_col       (begc:endc))      ; this%flx_p_ero_col          (:) = nan
    allocate( this%flx_q_ero_col       (begc:endc))      ; this%flx_q_ero_col          (:) = nan
    allocate( this%flx_lndsld_ero_col  (begc:endc))      ; this%flx_lndsld_ero_col     (:) = nan
    allocate( this%flx_sed_ero_col     (begc:endc))      ; this%flx_sed_ero_col        (:) = nan
    allocate( this%flx_lndsld_yld_col  (begc:endc))      ; this%flx_lndsld_yld_col     (:) = nan
    allocate( this%flx_sed_yld_col     (begc:endc))      ; this%flx_sed_yld_col        (:) = nan

    allocate( this%prnsf_col           (begc:endc))      ; this%prnsf_col              (:) = nan
    allocate( this%qrnsf_col           (begc:endc))      ; this%qrnsf_col              (:) = nan
    allocate( this%tcsf_col            (begc:endc))      ; this%tcsf_col               (:) = nan

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use histFileMod    , only : hist_addfld1d
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(sedflux_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer           :: begc, endc
    character(24)     :: fieldname
    character(100)    :: longname
    real(r8), pointer :: data1dptr(:)   ! temp. pointer for slicing larger arrays
    !------------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc

    this%flx_p_ero_col(begc:endc) = spval
    call hist_addfld1d (fname='P_ERO',  units='kg/m^2/s',  &
         avgflag='A', long_name='hillslope rainfall-driven erosion', &
         ptr_col=this%flx_p_ero_col, l2g_scale_type='veg', &
         default='inactive')

    this%flx_q_ero_col(begc:endc) = spval
    call hist_addfld1d (fname='Q_ERO',  units='kg/m^2/s',  &
         avgflag='A', long_name='hillslope runoff-driven erosion', &
         ptr_col=this%flx_q_ero_col, l2g_scale_type= 'veg', &
         default='inactive')

    this%flx_lndsld_ero_col(begc:endc) = spval
    call hist_addfld1d (fname='LNDSLD_ERO',  units='kg/m^2/s',  &
         avgflag='A', long_name='hillslope landslide-driven erosion', &
         ptr_col=this%flx_lndsld_ero_col, l2g_scale_type= 'veg', &
         default='inactive')

    this%flx_sed_ero_col(begc:endc) = spval
    call hist_addfld1d (fname='SED_ERO',  units='kg/m^2/s',  &
         avgflag='A', long_name='hillslope total erosion', &
         ptr_col=this%flx_sed_ero_col, l2g_scale_type= 'veg', &
         default='inactive')

    this%flx_lndsld_yld_col(begc:endc) = spval
    call hist_addfld1d (fname='LNDSLD_YLD',  units='kg/m^2/s',  &
         avgflag='A', long_name='hillslope landslide sediment yield', &
         ptr_col=this%flx_lndsld_yld_col, l2g_scale_type= 'veg', &
         default='inactive')

    this%flx_sed_yld_col(begc:endc) = spval
    call hist_addfld1d (fname='SED_YLD',  units='kg/m^2/s',  &
         avgflag='A', long_name='hillslope total sediment yield', &
         ptr_col=this%flx_sed_yld_col, l2g_scale_type= 'veg', &
         default='inactive')

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize cold start conditions for module variables
    !
    ! !USES:
    use shr_kind_mod   , only : r8 => shr_kind_r8
    use clm_varcon     , only : grlnd
    use clm_varctl     , only : fsurdat
    use fileutils      , only : getfil
    use ncdio_pio      , only : file_desc_t, ncd_io, ncd_pio_openfile, ncd_pio_closefile
    !
    ! !ARGUMENTS:
    class(sedflux_type)         :: this
    type(bounds_type) , intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer            :: c, g
    logical            :: readvar
    type(file_desc_t)  :: ncid
    character(len=256) :: locfn
    real(r8) ,pointer  :: prnsf2d         (:)      ! read in - prnsf
    real(r8) ,pointer  :: qrnsf2d         (:)      ! read in - qrnsf
    real(r8) ,pointer  :: tcsf2d          (:)      ! read in - tcsf
    !-----------------------------------------------------------------------

    allocate(prnsf2d(bounds%begg:bounds%endg))
    allocate(qrnsf2d(bounds%begg:bounds%endg))
    allocate(tcsf2d(bounds%begg:bounds%endg))

    call getfil (fsurdat, locfn, 0)
    call ncd_pio_openfile (ncid, locfn, 0)
    call ncd_io(ncid=ncid, varname='parEro_c1', flag='read', data=prnsf2d, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=' ERROR: parEro_c1 NOT on surfdata file'//errMsg(__FILE__, __LINE__)) 
    end if
    call ncd_io(ncid=ncid, varname='parEro_c2', flag='read', data=qrnsf2d, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=' ERROR: parEro_c2 NOT on surfdata file'//errMsg(__FILE__, __LINE__)) 
    end if
    call ncd_io(ncid=ncid, varname='parEro_c3', flag='read', data=tcsf2d, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=' ERROR: parEro_c3 NOT on surfdata file'//errMsg(__FILE__, __LINE__)) 
    end if
    call ncd_pio_closefile(ncid)

    do c = bounds%begc, bounds%endc
       this%flx_p_ero_col(c)                 = 0._r8
       this%flx_q_ero_col(c)                 = 0._r8
       this%flx_lndsld_ero_col(c)            = 0._r8
       this%flx_sed_ero_col(c)               = 0._r8
       this%flx_lndsld_yld_col(c)            = 0._r8
       this%flx_sed_yld_col(c)               = 0._r8

       g = col_pp%gridcell(c)
       this%prnsf_col(c)                     = prnsf2d(g)
       this%qrnsf_col(c)                     = qrnsf2d(g)
       this%tcsf_col(c)                      = tcsf2d(g)
    end do

    deallocate(prnsf2d, qrnsf2d, tcsf2d)

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use ncdio_pio       , only : file_desc_t
    !
    ! !ARGUMENTS:
    class(sedflux_type) :: this
    type(bounds_type), intent(in)    :: bounds
    type(file_desc_t), intent(inout) :: ncid
    character(len=*) , intent(in)    :: flag
    !
    ! !LOCAL VARIABLES:
    !-----------------------------------------------------------------------

    ! do nothing

  end subroutine Restart

end module SedFluxType

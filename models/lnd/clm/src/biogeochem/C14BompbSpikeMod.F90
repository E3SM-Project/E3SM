module C14BombSpikeMod

  !-----------------------------------------------------------------------
  ! Module for transient pulse simulation 
  !
  ! !USES:
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use clm_time_manager , only : get_curr_date,get_days_per_year
  use clm_varcon       , only : c14ratio, secspday
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: C14BombSpike
  public:: C14_init_BombSpike

  ! !PUBLIC TYPES:
  logical            , public :: use_c14_bombspike = .false. ! do we use time-varying atmospheric C14?
  character(len=256) , public :: atm_c14_filename = ' '      ! file name of C14 input data

  ! !PRIVATE TYPES:
  real(r8), allocatable, private :: atm_c14file_time(:)
  real(r8), allocatable, private :: atm_delta_c14(:)
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine C14BombSpike( rc14_atm )
    !
    ! !DESCRIPTION:
    ! for transient pulse simulation, impose a simplified bomb spike
    !
    ! !ARGUMENTS:
    real(r8), intent(out) :: rc14_atm
    !
    ! !LOCAL VARIABLES:
    integer  :: yr, mon, day, tod, offset
    real(r8) :: dateyear
    real(r8) :: delc14o2_atm 
    real(r8) :: days_per_year ! days per year
    integer  :: fp, p, nt
    integer  :: ind_below
    integer  :: ntim_atm_ts
    real(r8) :: twt_1, twt_2  ! weighting fractions for interpolating
    !-----------------------------------------------------------------------

    ! get current date
    call get_curr_date(yr, mon, day, tod, offset)
    days_per_year = get_days_per_year()
    dateyear = real(yr) + real(mon)/12._r8 + real(day)/days_per_year + real(tod)/(secspday*days_per_year)

    ! find points in atm timeseries to interpolate between
    ntim_atm_ts = size(atm_c14file_time)
    ind_below = 0
    do nt = 1, ntim_atm_ts
       if (dateyear >= atm_c14file_time(nt) ) then
          ind_below = ind_below+1
       endif
    end do

    ! interpolate between nearest two points in atm c14 timeseries
    if (ind_below .eq. 0 ) then 
       delc14o2_atm = atm_delta_c14(1)
    elseif (ind_below .eq. ntim_atm_ts ) then
       delc14o2_atm = atm_delta_c14(ntim_atm_ts)
    else
       twt_2 = min(1._r8, max(0._r8,(dateyear-atm_c14file_time(ind_below)) &
            / (atm_c14file_time(ind_below+1)-atm_c14file_time(ind_below))))
       twt_1 = 1._r8 - twt_2
       delc14o2_atm = atm_delta_c14(ind_below) * twt_1 +  atm_delta_c14(ind_below+1) * twt_2
    endif

    ! change delta units to ratio, put on patch loop

    rc14_atm = (delc14o2_atm * 1.e-3_r8 + 1._r8) * c14ratio

  end subroutine C14BombSpike

  !-----------------------------------------------------------------------
  subroutine C14_init_BombSpike()
    !
    ! !DESCRIPTION:
    ! read netcdf file containing a timeseries of atmospheric delta C14 values; save in module-level array 
    !
    ! !USES:
    use ncdio_pio
    use fileutils   , only : getfil
    use abortutils  , only : endrun
    use clm_varctl  , only : iulog
    use spmdMod     , only : masterproc
    use shr_log_mod , only : errMsg => shr_log_errMsg
    !
    ! !LOCAL VARIABLES:
    character(len=256) :: locfn           ! local file name
    type(file_desc_t)  :: ncid            ! netcdf id
    integer :: dimid,varid                ! input netCDF id's
    integer :: ntim                       ! number of input data time samples
    integer :: t
    !-----------------------------------------------------------------------

    if ( masterproc ) then
       write(iulog, *) 'C14_init_BombSpike: preparing to open file:'
       write(iulog, *) trim(locfn)
    endif

    call getfil(atm_c14_filename, locfn, 0)

    call ncd_pio_openfile (ncid, trim(locfn), 0)

    call ncd_inqdlen(ncid,dimid,ntim,'time')

    !! allocate arrays based on size of netcdf timeseries
    allocate(atm_c14file_time(ntim))
    allocate(atm_delta_c14(ntim))

    call ncd_io(ncid=ncid, varname='time', flag='read', data=atm_c14file_time)

    call ncd_io(ncid=ncid, varname='atm_delta_c14', flag='read', data=atm_delta_c14)

    call ncd_pio_closefile(ncid)

    ! check to make sure that time dimension is well behaved
    do t = 2, ntim
       if ( atm_c14file_time(t) - atm_c14file_time(t-1) <= 0._r8 ) then
          write(iulog, *) 'C14_init_BombSpike: error.  time axis must be monotonically increasing'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif
    end do

  end subroutine C14_init_BombSpike

end module C14BombSpikeMod

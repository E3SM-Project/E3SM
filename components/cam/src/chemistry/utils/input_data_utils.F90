!--------------------------------------------------------------------------------
! module where input data utility classes reside
!
! classes:
!     time_coordinate -- manages the time coordinate of input data sets
!--------------------------------------------------------------------------------
module input_data_utils
  use shr_kind_mod,   only : r8 => shr_kind_r8, cs => shr_kind_cs, cl=> shr_kind_cl
  use cam_abortutils, only : endrun
  use cam_logfile,    only : iulog
  use pio,            only : file_desc_t, pio_inq_dimid, pio_inq_dimlen, pio_get_att
  use pio,            only : pio_seterrorhandling, pio_get_var, pio_inq_varid
  use pio,            only : PIO_NOWRITE, PIO_BCAST_ERROR, PIO_INTERNAL_ERROR, PIO_NOERR
  use time_manager,   only : timemgr_get_calendar_cf, set_time_float_from_date, get_curr_date
  use spmd_utils,     only : masterproc
  
  implicit none

  private
  public :: time_coordinate

  type :: time_coordinate
     integer :: ntimes
     real(r8) :: wghts(2)
     integer :: indxs(2)
     real(r8), allocatable :: times(:)
     real(r8), allocatable :: time_bnds(:,:)
     logical :: time_interp = .true.
     logical :: fixed = .false.
     integer :: fixed_ymd, fixed_tod
     character(len=cl) :: filename
     real(r8) :: dtime ! time shift in interpolation point (days)
   contains
     procedure :: initialize
     procedure :: advance
     procedure :: read_more
     procedure :: copy
     procedure :: destroy
  end type time_coordinate

contains

  !-----------------------------------------------------------------------------
  ! initializer
  !-----------------------------------------------------------------------------
  subroutine initialize( this, filepath, fixed, fixed_ymd, fixed_tod, force_time_interp, set_weights, try_dates, delta_days )
    use ioFileMod,      only : getfil
    use cam_pio_utils,  only : cam_pio_openfile, cam_pio_closefile
    use string_utils,   only : to_upper

    class(time_coordinate), intent(inout) :: this
    character(len=*), intent(in) :: filepath
    logical, optional,intent(in) :: fixed
    integer, optional,intent(in) :: fixed_ymd
    integer, optional,intent(in) :: fixed_tod
    logical, optional,intent(in) :: force_time_interp
    logical, optional,intent(in) :: set_weights
    logical, optional,intent(in) :: try_dates
    real(r8), optional,intent(in) :: delta_days !  time shift in interpolation point (days) -- for previous day set this to -1.

    character(len=cl) :: filen
    character(len=cl) :: time_units, err_str
    character(len=cs) :: time_calendar, model_calendar
    character(len=4) :: yr_str
    character(len=2) :: mon_str, day_str, hr_str, min_str, sec_str
    integer :: ref_yr, ref_mon, ref_day, ref_hr, ref_min, ref_sec, tod
    integer :: varid, ierr
    real(r8) :: ref_time

    integer,  allocatable :: dates(:)
    integer,  allocatable :: datesecs(:)
    real(r8), allocatable :: times_file(:)
    real(r8), allocatable :: times_modl(:)
    real(r8), allocatable :: time_bnds_file(:,:)
    type(file_desc_t) :: fileid
    logical :: force_interp
    logical :: set_wghts
    logical :: use_time, adj_times, use_time_bnds
    integer :: i, ind

    if (present(fixed)) this%fixed = fixed
    if (present(fixed_ymd)) this%fixed_ymd = fixed_ymd
    if (present(fixed_tod)) this%fixed_tod = fixed_tod

    this%dtime = 0._r8
    if (present(delta_days)) this%dtime = delta_days

    if (present(force_time_interp)) then
       force_interp = force_time_interp
    else
       force_interp = .false.
    endif

    if (present(set_weights)) then
       set_wghts = set_weights
    else
       set_wghts = .true.
    endif

    this%filename = trim(filepath)

    call getfil( filepath, filen, 0 )
    call cam_pio_openfile( fileid, filen, PIO_NOWRITE )

    call pio_seterrorhandling( fileid, PIO_BCAST_ERROR)

    call get_dimension( fileid, 'time', this%ntimes )
    allocate ( times_file( this%ntimes ) )
    allocate ( times_modl( this%ntimes ) )
    allocate ( this%times( this%ntimes ) )

    ierr =  pio_inq_varid( fileid, 'time', varid )
    use_time = ierr.eq.PIO_NOERR
    ierr = pio_get_att( fileid, varid, 'calendar', time_calendar)
    use_time = ierr.eq.PIO_NOERR .and. use_time
    ierr = pio_get_att( fileid, varid, 'units', time_units)
    use_time = ierr.eq.PIO_NOERR .and. use_time
    if (use_time) then
       use_time = time_units(1:10).eq.'days since'
    endif

    if (present(try_dates)) then
       if (try_dates) then
          ierr = pio_inq_varid( fileid, 'date', varid  )
          use_time = ierr .ne. PIO_NOERR
       endif
    endif

    adj_times = .false.
    use_time_bnds =  .false.

    time_var_use: if (use_time) then

       ! check the calendar attribute - must match model calendar
       model_calendar = timemgr_get_calendar_cf()

       if (this%ntimes>2) then 
          ! if only 2 time records then it is assumed that the input has 2 identical time records
          !  -- climatological or solar-cycle avaraged

          adj_times = (to_upper(time_calendar(1:6)) .ne. to_upper(model_calendar(1:6)))

          if (adj_times .and. masterproc) then
             write(iulog,*) 'time_coordinate%initialize: model calendar '//trim(model_calendar)// &
                            ' does not match input data calendar '//trim(time_calendar)
             write(iulog,*) ' -- will try to use date and datesec in the input file to adjust the time coordinate.'
          end if
       end if

       ! parse out ref date and time
       !  time:units = "days since YYYY-MM-DD hh:mm:ss" ;

       yr_str  = time_units(12:15)
       mon_str = time_units(17:18)
       day_str = time_units(20:21)
       hr_str  = time_units(23:24)
       min_str = time_units(26:27)
       
       ind = index(hr_str,':')
       if( ind == 1 ) hr_str = hr_str(2:2)
       if( ind == 2 ) hr_str = hr_str(1:1)

       ind = index(min_str,':')
       if( ind == 1 ) min_str = min_str(2:2)
       if( ind == 2 ) min_str = min_str(1:1)
       
       read( yr_str,  * ) ref_yr
       read( mon_str, * ) ref_mon
       read( day_str, * ) ref_day
       read( hr_str,  * ) ref_hr
       read( min_str, * ) ref_min
       if (len_trim(time_units).ge.30) then
          sec_str = time_units(29:30)
          read( sec_str, * ) ref_sec
       else
          ref_sec = 0
       endif

       tod = ref_hr*3600 + ref_min*60 + ref_sec
       call set_time_float_from_date( ref_time, ref_yr, ref_mon, ref_day, tod )

       ierr = pio_get_var( fileid, varid, times_file )
       if (ierr.ne.PIO_NOERR) then
          call endrun('time_coordinate%initialize: not able to read times')
       endif

       times_file = times_file + ref_time

       ierr =  pio_inq_varid( fileid, 'time_bnds', varid )

       use_time_bnds = (ierr==PIO_NOERR .and. .not.force_interp)
       
       if (use_time_bnds) then
          allocate ( this%time_bnds( 2, this%ntimes ) )
          allocate ( time_bnds_file( 2, this%ntimes ) )
          ierr =  pio_get_var( fileid, varid, time_bnds_file )
          time_bnds_file = time_bnds_file + ref_time
          this%time_interp = .false.
          do i = 1,this%ntimes
            if (.not. (time_bnds_file(1,i)<times_file(i) &
                 .and. time_bnds_file(2,i)>times_file(i)) ) then
                write(err_str,*) 'incorrect time_bnds -- time index: ',i,' file: '//trim(filepath)
                call endrun(err_str)
            endif
          enddo
       else
          this%time_interp = .true.
       endif
    else
       this%time_interp = .true.
    endif time_var_use

    read_dates: if (adj_times .or. .not.use_time) then

       ! try using date and datesec
       allocate(dates(this%ntimes), stat=ierr )
       if( ierr /= 0 ) then
          write(iulog,*) 'time_coordinate%initialize: failed to allocate dates; error = ',ierr
          call endrun('time_coordinate%initialize: failed to allocate dates')
       end if

       allocate(datesecs(this%ntimes), stat=ierr )
       if( ierr /= 0 ) then
          write(iulog,*) 'time_coordinate%initialize: failed to allocate datesecs; error = ',ierr
          call endrun('time_coordinate%initialize: failed to allocate datesecs')
       end if

       ierr = pio_inq_varid( fileid, 'date', varid  )
       if (ierr/=PIO_NOERR) then
          call endrun('time_coordinate%initialize: input file must contain time or date variable '//trim(filepath))
       endif
       ierr = pio_get_var( fileid, varid, dates )
       ierr = pio_inq_varid( fileid, 'datesec', varid )
       if (ierr==PIO_NOERR) then
          ierr = pio_get_var( fileid, varid, datesecs )
       else
          datesecs(:) = 0
       endif

       !datesecs(:) = 0
       call convert_dates( dates, datesecs, times_modl )

       deallocate( dates, datesecs )

    endif read_dates

    if (adj_times) then
       ! time_bnds_modl - time_bnds_file = times_modl - times_file
       ! time_bnds_modl = time_bnds_file + times_modl - times_file
       this%times(:) = times_modl(:)
       if (use_time_bnds) then 
          this%time_bnds(1,:) = time_bnds_file(1,:) + times_modl(:) - times_file(:)
          this%time_bnds(2,:) = time_bnds_file(2,:) + times_modl(:) - times_file(:)
       endif
    else if (use_time) then
       this%times(:) = times_file(:)
       if (use_time_bnds) then 
          this%time_bnds(1,:) = time_bnds_file(1,:)
          this%time_bnds(2,:) = time_bnds_file(2,:)
       endif
    else
       this%times(:) = times_modl(:)
    endif

    deallocate( times_modl, times_file )
    if (use_time_bnds) deallocate(time_bnds_file)

    call pio_seterrorhandling(fileid, PIO_INTERNAL_ERROR)

    call cam_pio_closefile(fileid)

    this%indxs(1)=1
    if (set_wghts) call set_wghts_indices(this)

  end subroutine initialize

  !-----------------------------------------------------------------------------
  ! advance the time coordinate
  !-----------------------------------------------------------------------------
  subroutine advance( this )
    class(time_coordinate) :: this

    if (.not.this%fixed) call set_wghts_indices(this)

  end subroutine advance

  !-----------------------------------------------------------------------------
  ! determine if need to read more data from input data set
  !-----------------------------------------------------------------------------
  function read_more(this) result(check)
    class(time_coordinate), intent(in) :: this
    logical :: check

    real(r8) :: model_time

    model_time = get_model_time() + this%dtime

    if (.not.this%fixed) then
       if (allocated(this%time_bnds)) then
          check = model_time > this%time_bnds(2,this%indxs(1))
       else
          check = model_time > this%times(this%indxs(2))
       endif
    else
       check = .false.
    endif

  end function read_more

  !-----------------------------------------------------------------------------
  ! destroy method -- deallocate memory and revert to default settings
  !-----------------------------------------------------------------------------
  subroutine destroy( this )
        class(time_coordinate), intent(inout) :: this
    
    if (allocated(this%times)) deallocate(this%times)
    if (allocated(this%time_bnds)) deallocate(this%time_bnds)
    this%ntimes = 0
    this%filename='NONE'

  end subroutine destroy
    
  !-----------------------------------------------------------------------------
  ! produce a duplicate time coordinate object
  !-----------------------------------------------------------------------------
  subroutine copy( this, obj )
    class(time_coordinate), intent(inout) :: this
    class(time_coordinate), intent(in) :: obj

    call this%destroy()

    this%ntimes = obj%ntimes
    this%fixed  = obj%fixed
    this%fixed_ymd = obj%fixed_ymd
    this%fixed_tod = obj%fixed_tod

    allocate ( this%times( this%ntimes ) )
    this%times = obj%times

    if (allocated( obj%time_bnds )) then
       allocate ( this%time_bnds( 2, this%ntimes ) )
       this%time_bnds = obj%time_bnds
    endif
    this%filename = obj%filename

  end subroutine copy

! private methods

  !-----------------------------------------------------------------------
  ! set time interpolation weights
  !-----------------------------------------------------------------------
  subroutine set_wghts_indices(obj)

    class(time_coordinate), intent(inout) :: obj

    real(r8) :: model_time
    real(r8) :: datatm, datatp
    integer :: yr, mon, day
    integer :: index, i
    character(len=cl) :: errmsg

    !if(masterproc)write(102,*)'input_data-set_weight'
    ! set time indices and time-interpolation weights 
    fixed_time: if (obj%fixed) then
       yr = obj%fixed_ymd/10000
       mon = (obj%fixed_ymd-yr*10000) / 100
       day = obj%fixed_ymd-yr*10000-mon*100
       call set_time_float_from_date( model_time, yr, mon, day, obj%fixed_tod )
       model_time = model_time + obj%dtime
    else
       model_time = get_model_time() + obj%dtime
    endif fixed_time

    index = -1

    findtimes: do i = obj%indxs(1), obj%ntimes
       if (allocated(obj%time_bnds)) then
          datatm = obj%time_bnds(1,i)
          datatp = obj%time_bnds(2,i)
       else
          if (i .ge. obj%ntimes) then
             errmsg = '1:input_data_utils::set_wghts_indices cannot not find model time in: '&
                    // trim(obj%filename)
             write(iulog,*) trim(errmsg)
             call endrun(trim(errmsg))
          endif
          datatm = obj%times(i)
          datatp = obj%times(i+1)
       endif
       if ( model_time .lt. datatm ) then
          errmsg = '2:input_data_utils::set_wghts_indices cannot not find model time in: '&
                 // trim(obj%filename)
          write(iulog,*) trim(errmsg)
          call endrun(trim(errmsg))
       endif
       if ( model_time .ge. datatm .and. model_time .le. datatp ) then
          index = i
          obj%indxs(1) = i
          obj%indxs(2) = i+1
          exit findtimes
       endif
    enddo findtimes

    if ((allocated(obj%time_bnds)) .and. (i<obj%ntimes)) then
       if (.not.(obj%time_bnds(1,i+1) > obj%time_bnds(1,i))) then
          obj%indxs = obj%indxs+1  ! skip 29 Feb when calendar is noleap
       endif
    endif

    if (.not.(index>0.and.index<obj%ntimes)) then
       errmsg = 'input_data_utils::set_wghts_indices cannot not find time indices for input file: '&
            // trim(obj%filename)
       write(iulog,*) trim(errmsg)
       call endrun(trim(errmsg))
    endif

    if (obj%time_interp) then
       obj%wghts(2) = ( model_time - obj%times(index) ) / ( obj%times(index+1) - obj%times(index) )
       obj%wghts(1) = 1._r8 - obj%wghts(2)
    else
       obj%wghts(1) = 1._r8
       obj%wghts(2) = 0._r8       
    endif

  end subroutine set_wghts_indices

  !-----------------------------------------------------------------------
  ! returns dimension size
  !-----------------------------------------------------------------------
  subroutine get_dimension( fid, dname, dsize )
    type(file_desc_t), intent(in) :: fid
    character(*), intent(in) :: dname
    integer, intent(out) :: dsize

    integer :: dimid, ierr

    ierr = pio_inq_dimid( fid, dname, dimid )
    ierr = pio_inq_dimlen( fid, dimid, dsize )

  end subroutine get_dimension

  !-----------------------------------------------------------------------
  ! returns a real which represents the current model time
  !-----------------------------------------------------------------------
  function get_model_time() result(time)

    real(r8) :: time

    integer yr, mon, day, ncsec  ! components of a date

    call get_curr_date(yr, mon, day, ncsec)

    call set_time_float_from_date( time, yr, mon, day, ncsec )

  end function get_model_time

  !---------------------------------------------------------------------------
  ! convert a collection of dates and times to reals
  !---------------------------------------------------------------------------
  subroutine convert_dates( dates, secs, times )

    use time_manager, only: set_time_float_from_date

    integer,  intent(in)  :: dates(:)
    integer,  intent(in)  :: secs(:)

    real(r8), intent(out) :: times(:)

    integer :: year, month, day, sec,n ,i

    n = size( dates ) 

    do i=1,n
       year = dates(i)/10000
       month = (dates(i)-year*10000)/100
       day = dates(i)-year*10000-month*100
       sec = secs(i)
       call set_time_float_from_date( times(i), year, month, day, sec )
    enddo

  end subroutine convert_dates

end module input_data_utils

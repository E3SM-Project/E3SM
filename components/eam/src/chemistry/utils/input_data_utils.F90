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
   use shr_cal_mod,    only : shr_cal_useNoleapInput, shr_cal_calendarName, shr_cal_noleap, shr_cal_gregorian
   use shr_cal_mod,    only : shr_cal_advDate, shr_cal_date2ymd, shr_cal_validDate
   use shr_cal_mod,    only : shr_cal_inputTimeMethod
   use shr_cal_mod,    only : shr_cal_useInputPolicy, shr_cal_useGregorianInput, shr_cal_gregorianElapsedDays
   use shr_cal_mod,    only : shr_cal_resolveInputFile
   use shr_input_calendar_mod, only : input_calendar_check_times, input_calendar_check_date
   use shr_input_calendar_mod, only : input_calendar_flavor, input_calendar_check_range
   use shr_string_mod, only : shr_string_parseCFtunit
   use netcdf, only : nf90_enotatt
   use pio,           only : pio_inquire, pio_inquire_variable

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
   logical :: input_compatibility = .false.
   logical :: input_gregorian_view = .false.
   logical :: input_leapday_active = .false.
   integer :: input_leapday_count = 0
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
   character(64) :: source_flavor, source_units
   character(len=cs) :: declared_calendar
   logical :: source_bounds_declared
   integer :: source_epoch
   real(r8) :: source_tod
    character(len=4) :: yr_str
    character(len=2) :: mon_str, day_str, hr_str, min_str, sec_str
    integer :: ref_yr, ref_mon, ref_day, ref_hr, ref_min, ref_sec, tod
   integer :: varid, ierr, field_id, field_count
   character(len=cl) :: field_name
   character(len=2048) :: cell_methods
   character(len=32) :: time_method
   character(len=cl) :: time_bounds_name
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
   this%input_compatibility = .false.
   this%input_gregorian_view = .false.
   this%input_leapday_active = .false.
   this%input_leapday_count = 0

    call getfil( filepath, filen, 0 )
    call cam_pio_openfile( fileid, filen, PIO_NOWRITE )

    call pio_seterrorhandling( fileid, PIO_BCAST_ERROR)

    call get_dimension( fileid, 'time', this%ntimes )
      if (shr_cal_useInputPolicy() .and. this%ntimes == 0) &
          call endrun('NOLEAP_INPUT: empty time coordinate; file='//trim(filepath))
    allocate ( times_file( this%ntimes ) )
    allocate ( times_modl( this%ntimes ) )
    allocate ( this%times( this%ntimes ) )

    ierr =  pio_inq_varid( fileid, 'time', varid )
    use_time = ierr.eq.PIO_NOERR
      if (shr_cal_useInputPolicy() .and. .not. use_time) &
          call endrun('NOLEAP_INPUT: time coordinate metadata is required; file='//trim(filepath))
    time_calendar = ''
    ierr = pio_get_att( fileid, varid, 'calendar', time_calendar)
   declared_calendar = time_calendar
   source_flavor = input_calendar_flavor(time_calendar)
   if (shr_cal_useInputPolicy()) then
        if (ierr /= PIO_NOERR .and. ierr /= nf90_enotatt) &
           call endrun('NOLEAP_INPUT: cannot read calendar metadata; file='//trim(filepath))
      call shr_cal_resolveInputFile(trim(filepath),'',time_calendar,ierr == PIO_NOERR,model_calendar,source_flavor)
        time_calendar = model_calendar
        ierr = PIO_NOERR
       model_calendar = shr_cal_calendarName(timemgr_get_calendar_cf())
      if (shr_cal_useNoleapInput() .and. model_calendar /= shr_cal_gregorian) &
            call endrun('NOLEAP_INPUT: fudge_noleap_for_gregorian requires the Gregorian model calendar')
        if (shr_cal_useGregorianInput() .and. model_calendar /= shr_cal_noleap) &
           call endrun('GREGORIAN_INPUT: fudge_gregorian_to_noleap requires the no-leap model calendar')
        this%input_compatibility = shr_cal_useNoleapInput() .and. &
           shr_cal_calendarName(time_calendar) == shr_cal_noleap
        this%input_gregorian_view = shr_cal_useGregorianInput() .and. time_calendar == shr_cal_gregorian
        if (this%input_gregorian_view .and. masterproc) then
         write(iulog,*) 'WARNING GREGORIAN_INPUT: EAM chemistry file=',trim(filepath)
         write(iulog,*) 'WARNING GREGORIAN_INPUT: source=GREGORIAN model=NO_LEAP; retain February 29 brackets'
         write(iulog,*) 'WARNING GREGORIAN_INPUT: model time skips February 29; annual totals are NOT conserved'
        endif
       if (this%input_compatibility .and. masterproc) then
          write(iulog,*) 'WARNING NOLEAP_INPUT: EAM chemistry file=',trim(filepath)
          write(iulog,*) 'WARNING NOLEAP_INPUT: source=NO_LEAP model=GREGORIAN; preserve record date labels'
          write(iulog,*) 'WARNING NOLEAP_INPUT: time bounds are mapped; annual integrals are NOT conserved'
       endif
    endif
    use_time = ierr.eq.PIO_NOERR .and. use_time
    ierr = pio_get_att( fileid, varid, 'units', time_units)
    use_time = ierr.eq.PIO_NOERR .and. use_time
    if (use_time) then
       use_time = time_units(1:10).eq.'days since'
    endif
    time_bounds_name = 'time_bnds'
   source_bounds_declared = .false.
   if (this%input_compatibility .or. this%input_gregorian_view .or. source_flavor == 'STANDARD') then
       ierr = pio_get_att(fileid,varid,'bounds',time_bounds_name)
      source_bounds_declared = ierr == PIO_NOERR
      if (source_flavor == 'STANDARD' .and. ierr /= PIO_NOERR .and. ierr /= nf90_enotatt) &
         call endrun('INPUT_CALENDAR: cannot read bounds metadata; file='//trim(filepath))
       if (ierr /= PIO_NOERR) time_bounds_name = 'time_bnds'
    endif

    if (source_flavor == 'STANDARD') then
       if (.not. use_time) call endrun('INPUT_CALENDAR: mixed Gregorian input requires days since units; file='//trim(filepath))
       call shr_string_parseCFtunit(time_units,source_units,source_epoch,source_tod)
       ierr = pio_get_var(fileid,varid,times_file)
       if (ierr /= PIO_NOERR) call endrun('INPUT_CALENDAR: cannot read time; file='//trim(filepath))
       call input_calendar_check_range(source_flavor,source_epoch,source_tod,source_units,times_file,err_str)
        if (err_str /= '') call endrun(trim(err_str)//'; declared='//trim(declared_calendar)//'; file='//trim(filepath))
       ierr = pio_inq_varid(fileid,trim(time_bounds_name),field_id)
        if (ierr /= PIO_NOERR .and. source_bounds_declared) &
           call endrun('INPUT_CALENDAR: missing declared time bounds; file='//trim(filepath))
       if (ierr == PIO_NOERR) then
          allocate(time_bnds_file(2,this%ntimes))
          ierr = pio_get_var(fileid,field_id,time_bnds_file)
          if (ierr /= PIO_NOERR) call endrun('INPUT_CALENDAR: cannot read bounds; file='//trim(filepath))
          call input_calendar_check_range(source_flavor,source_epoch,source_tod,source_units, &
               reshape(time_bnds_file,[size(time_bnds_file)]),err_str)
          if (err_str /= '') call endrun(trim(err_str)//'; bounds; file='//trim(filepath))
          deallocate(time_bnds_file)
       end if
    end if

    if (present(try_dates)) then
       if (try_dates) then
          ierr = pio_inq_varid( fileid, 'date', field_id )
          if (this%input_compatibility .or. this%input_gregorian_view) then
             use_time = use_time .and. ierr /= PIO_NOERR
          else
             use_time = ierr .ne. PIO_NOERR
          endif
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

         adj_times = (to_upper(time_calendar(1:6)) .ne. to_upper(model_calendar(1:6))) &
            .and. .not. (this%input_compatibility .or. this%input_gregorian_view)

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
      if (shr_cal_useInputPolicy()) then
          call input_calendar_check_date(ref_yr*10000+ref_mon*100+ref_day,real(tod,r8),time_calendar,err_str)
          if (err_str /= '') call endrun(trim(err_str)//'; file='//trim(filepath))
       endif
       if (this%input_gregorian_view) then
          ref_time = shr_cal_gregorianElapsedDays(ref_yr*10000+ref_mon*100+ref_day,real(tod,r8))
       else
          call set_time_float_from_date( ref_time, ref_yr, ref_mon, ref_day, tod )
       endif

       ierr = pio_get_var( fileid, varid, times_file )
       if (ierr.ne.PIO_NOERR) then
          call endrun('time_coordinate%initialize: not able to read times')
       endif
      if (shr_cal_useInputPolicy()) then
          call input_calendar_check_times(times_file,err_str)
          if (err_str /= '') call endrun(trim(err_str)//'; file='//trim(filepath))
       endif

       if (this%input_compatibility) then
          call map_noleap_offsets(times_file,ref_yr*10000+ref_mon*100+ref_day,tod)
          adj_times = .false.
       else
          times_file = times_file + ref_time
       endif

      ierr =  pio_inq_varid( fileid, trim(time_bounds_name), varid )

       use_time_bnds = (ierr==PIO_NOERR .and. .not.force_interp)

       if (use_time_bnds) then
          allocate ( this%time_bnds( 2, this%ntimes ) )
          allocate ( time_bnds_file( 2, this%ntimes ) )
          ierr =  pio_get_var( fileid, varid, time_bnds_file )
         if (shr_cal_useInputPolicy() .and. ierr /= PIO_NOERR) &
            call endrun('NOLEAP_INPUT: cannot read chemistry time bounds; file='//trim(filepath))
          if (this%input_compatibility) then
             call map_noleap_offsets(time_bnds_file(1,:),ref_yr*10000+ref_mon*100+ref_day,tod)
             call map_noleap_offsets(time_bnds_file(2,:),ref_yr*10000+ref_mon*100+ref_day,tod)
          else
             time_bnds_file = time_bnds_file + ref_time
          endif
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
      if (shr_cal_useInputPolicy() .and. ierr /= PIO_NOERR) &
           call endrun('NOLEAP_INPUT: cannot read source dates; file='//trim(filepath))
       if (this%input_compatibility) then
          if (ierr /= PIO_NOERR) call endrun('NOLEAP_INPUT: cannot read source dates; file='//trim(filepath))
          do i = 1,this%ntimes
             if (.not. shr_cal_validDate(dates(i),shr_cal_noleap)) &
                  call endrun('NOLEAP_INPUT: invalid no-leap source date; file='//trim(filepath))
          enddo
       endif
       ierr = pio_inq_varid( fileid, 'datesec', varid )
       if (ierr==PIO_NOERR) then
          ierr = pio_get_var( fileid, varid, datesecs )
           if (shr_cal_useInputPolicy() .and. ierr /= PIO_NOERR) &
                call endrun('NOLEAP_INPUT: cannot read source datesecs; file='//trim(filepath))
       else
          datesecs(:) = 0
       endif

      if (shr_cal_useInputPolicy()) then
          do i = 1,this%ntimes
             call input_calendar_check_date(dates(i),real(datesecs(i),r8),time_calendar,err_str)
             if (err_str /= '') call endrun(trim(err_str)//'; file='//trim(filepath))
          enddo
       endif
       if (source_flavor == 'STANDARD') then
          do i = 1,this%ntimes
             call input_calendar_check_range(source_flavor,dates(i),real(datesecs(i),r8),'seconds',[0._r8],err_str)
             if (err_str /= '') call endrun(trim(err_str)//'; explicit dates; file='//trim(filepath))
          end do
       end if

       !datesecs(:) = 0
      call convert_dates( dates, datesecs, times_modl, this%input_gregorian_view )

       deallocate( dates, datesecs )

    endif read_dates

   if (shr_cal_useInputPolicy() .and. (adj_times .or. .not. use_time)) then
       call input_calendar_check_times(times_modl,err_str)
       if (err_str /= '') call endrun(trim(err_str)//'; file='//trim(filepath))
    endif

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

   if (this%input_compatibility .or. this%input_gregorian_view) then
       ierr = pio_inquire(fileid,nVariables=field_count)
       if (ierr /= PIO_NOERR) call endrun('NOLEAP_INPUT: cannot inquire chemistry input variables')
       do field_id = 1,field_count
          ierr = pio_inquire_variable(fileid,field_id,name=field_name)
          if (ierr /= PIO_NOERR) call endrun('NOLEAP_INPUT: cannot inquire chemistry input field')
          cell_methods = ''
          ierr = pio_get_att(fileid,field_id,'cell_methods',cell_methods)
          if (ierr /= PIO_NOERR) cycle
         time_method = shr_cal_inputTimeMethod(cell_methods)
         if (time_method /= '' .and. time_method /= 'point' .and. time_method /= 'mean') &
            call endrun('NOLEAP_INPUT: unsupported chemistry time method='//trim(time_method)//' field='// &
               trim(field_name)//' file='//trim(filepath))
         if (time_method == 'mean' .and. .not. use_time_bnds) &
               call endrun('NOLEAP_INPUT: chemistry time means require mapped time_bnds; field='//trim(field_name))
       enddo
    endif

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

   model_time = get_model_time(this%input_gregorian_view) + this%dtime

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

   if (this%input_compatibility .and. masterproc) write(iulog,*) &
       'WARNING NOLEAP_INPUT: chemistry summary; file=',trim(this%filename), &
       ' leap days in this segment=',this%input_leapday_count
    if (allocated(this%times)) deallocate(this%times)
    if (allocated(this%time_bnds)) deallocate(this%time_bnds)
    this%ntimes = 0
    this%filename='NONE'
   this%input_compatibility = .false.
   this%input_gregorian_view = .false.
   this%input_leapday_active = .false.
   this%input_leapday_count = 0

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
   this%input_compatibility = obj%input_compatibility
   this%input_gregorian_view = obj%input_gregorian_view
   this%dtime = obj%dtime
   this%time_interp = obj%time_interp

  end subroutine copy

! private methods

  !-----------------------------------------------------------------------
  ! set time interpolation weights
  !-----------------------------------------------------------------------
  subroutine set_wghts_indices(obj)

    class(time_coordinate), intent(inout) :: obj

    real(r8) :: model_time
    real(r8) :: datatm, datatp
   integer :: yr, mon, day, current_tod
   logical :: leapday_active
    integer :: index, i
    character(len=cl) :: errmsg

    ! set time indices and time-interpolation weights
    fixed_time: if (obj%fixed) then
       yr = obj%fixed_ymd/10000
       mon = (obj%fixed_ymd-yr*10000) / 100
       day = obj%fixed_ymd-yr*10000-mon*100
       if (obj%input_gregorian_view) then
          model_time = shr_cal_gregorianElapsedDays(obj%fixed_ymd,real(obj%fixed_tod,r8))
       else
          call set_time_float_from_date( model_time, yr, mon, day, obj%fixed_tod )
       endif
       model_time = model_time + obj%dtime
    else
      model_time = get_model_time(obj%input_gregorian_view) + obj%dtime
    endif fixed_time

    if (obj%input_compatibility) then
       call get_curr_date(yr,mon,day,current_tod)
       leapday_active = mon == 2 .and. day == 29
       if (leapday_active .neqv. obj%input_leapday_active) then
          if (leapday_active) obj%input_leapday_count = obj%input_leapday_count + 1
          if (masterproc) then
             if (leapday_active) then
                write(iulog,*) 'WARNING NOLEAP_INPUT: chemistry entering February 29; file=',trim(obj%filename)
             else
                write(iulog,*) 'WARNING NOLEAP_INPUT: chemistry leaving February 29; file=',trim(obj%filename)
             endif
          endif
          obj%input_leapday_active = leapday_active
       endif
    endif

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

  subroutine map_noleap_offsets(values,reference_date,reference_tod)
      use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    real(r8), intent(inout) :: values(:)
    integer, intent(in) :: reference_date, reference_tod
    integer :: time_index, source_date, source_tod, year, month, day
    real(r8) :: seconds, mapped_time

    if (.not. shr_cal_validDate(reference_date,shr_cal_noleap)) &
       call endrun('NOLEAP_INPUT: invalid no-leap reference date')
    do time_index = 1,size(values)
    if (.not. ieee_is_finite(values(time_index))) &
      call endrun('NOLEAP_INPUT: non-finite no-leap source time offset')
     call shr_cal_advDate(values(time_index),'days',reference_date,real(reference_tod,r8), &
        source_date,seconds,shr_cal_noleap)
     source_tod = floor(seconds)
     call shr_cal_date2ymd(source_date,year,month,day)
     call set_time_float_from_date(mapped_time,year,month,day,source_tod)
     values(time_index) = mapped_time + (seconds-real(source_tod,r8))/86400._r8
    enddo
  end subroutine map_noleap_offsets

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
   function get_model_time(gregorian_view) result(time)

    real(r8) :: time
      logical, intent(in) :: gregorian_view

    integer yr, mon, day, ncsec  ! components of a date

    call get_curr_date(yr, mon, day, ncsec)

    if (gregorian_view) then
       time = shr_cal_gregorianElapsedDays(yr*10000+mon*100+day,real(ncsec,r8))
    else
       call set_time_float_from_date( time, yr, mon, day, ncsec )
    endif

  end function get_model_time

  !---------------------------------------------------------------------------
  ! convert a collection of dates and times to reals
  !---------------------------------------------------------------------------
   subroutine convert_dates( dates, secs, times, gregorian_view )

    use time_manager, only: set_time_float_from_date

    integer,  intent(in)  :: dates(:)
    integer,  intent(in)  :: secs(:)

    real(r8), intent(out) :: times(:)
   logical, intent(in) :: gregorian_view

    integer :: year, month, day, sec,n ,i

    n = size( dates )

    do i=1,n
       year = dates(i)/10000
       month = (dates(i)-year*10000)/100
       day = dates(i)-year*10000-month*100
       sec = secs(i)
       if (gregorian_view) then
          times(i) = shr_cal_gregorianElapsedDays(dates(i),real(sec,r8))
       else
          call set_time_float_from_date( times(i), year, month, day, sec )
       endif
    enddo

  end subroutine convert_dates

end module input_data_utils

!-----------------------------------------------------------------------
! chlorine loading data for LINOZ
!-----------------------------------------------------------------------
module chlorine_loading_data
  use shr_kind_mod, only: r8 => shr_kind_r8
  use spmd_utils,   only: masterproc
  use abortutils,   only: endrun
  use cam_logfile,  only: iulog
  use linoz_data,   only: has_linoz_data

  implicit none

  save
  private
  public :: chlorine_loading_init
  public :: chlorine_loading_advance

  public :: chlorine_loading

  real(r8) :: chlorine_loading

  integer :: ntimes

  integer, parameter :: nt = 2
  real(r8) :: iloading(nt)
  
  real(r8), allocatable :: data_times(:)

  integer :: last_index = 1

  logical :: initialized = .false.

  logical  :: fixed
  real(r8) :: offset_time

! namelist vars
  character(len=256) :: chlorine_loading_file = ''
  character(len=8)   :: chlorine_loading_type = 'SERIAL' ! "FIXED" or "SERIAL"
  integer            :: chlorine_loading_ymd = 0         ! YYYYMMDD for "FIXED" type
  integer            :: chlorine_loading_tod = 0         ! seconds of day for "FIXED" type

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine chlorine_loading_init( file, type, ymd, tod )
    use cam_pio_utils, only : cam_pio_openfile
    use pio, only : file_desc_t, pio_noerr, pio_inq_dimid, pio_inq_varid, pio_get_var, pio_inq_dimlen, &
         pio_internal_error, pio_bcast_error, pio_seterrorhandling, pio_closefile
    use ioFileMod, only : getfil

  ! inputs
    character(len=256),         intent(in) :: file
    character(len=8), optional, intent(in) :: type
    integer,          optional, intent(in) :: ymd
    integer,          optional, intent(in) :: tod

  ! local vars
    integer :: astat, dimid, vid, data_vid
    character(len=256) :: filen   
    real(r8), allocatable :: loading(:)
    integer,  allocatable :: dates(:)
    integer,  allocatable :: datesecs(:)
    type(file_desc_t)     :: file_id
    
    integer :: i, ierr

    real(r8) :: model_time, time

    chlorine_loading_file = file
  
    if (.not.has_linoz_data) return

    if ( present(type) ) then
       chlorine_loading_type = type
    endif
    if ( present(ymd) ) then
       chlorine_loading_ymd = ymd
    endif
    if ( present(tod) ) then
       chlorine_loading_tod = tod
    endif

    fixed = trim(chlorine_loading_type) == 'FIXED'
    
    if (masterproc) then
       write(iulog,*) 'chlorine_loading_init: chlorine_loading_file = ',trim(chlorine_loading_file)
       write(iulog,*) 'chlorine_loading_init: chlorine_loading_type = ',trim(chlorine_loading_type)
       write(iulog,*) 'chlorine_loading_init: chlorine_loading_ymd = ',chlorine_loading_ymd
       write(iulog,*) 'chlorine_loading_init: chlorine_loading_tod = ',chlorine_loading_tod
    endif


    call getfil( chlorine_loading_file, filen, 0 )
    call cam_pio_openfile( file_id, filen, 0 )
    if ( masterproc ) write(iulog,*)'chlorine_loading_init: data file = ',trim(filen)
    ierr = pio_inq_dimid( file_id, 'time', dimid )
    ierr = pio_inq_dimlen( file_id, dimid, ntimes )
    ierr = pio_inq_varid( file_id, 'chlorine_loading', data_vid )

    allocate(data_times(ntimes), stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'chlorine_loading_init: failed to allocate data_times; error = ',astat
       call endrun('chlorine_loading_init')
    end if

    allocate(dates(ntimes), stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'chlorine_loading_init: failed to allocate dates; error = ',astat
       call endrun('chlorine_loading_init')
    end if

    allocate(datesecs(ntimes), stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'chlorine_loading_init: failed to allocate datesecs; error = ',astat
       call endrun('chlorine_loading_init')
    end if


    ierr = pio_inq_varid( file_id, 'date', vid  )
    ierr = pio_get_var( file_id, vid, dates )
    call pio_seterrorhandling(file_id, PIO_BCAST_ERROR)
    ierr = pio_inq_varid( file_id, 'datesec', vid )
    call pio_seterrorhandling(file_id, PIO_INTERNAL_ERROR)
    if (ierr == PIO_NOERR) then
       ierr = pio_get_var( file_id, vid, datesecs )
    else
       datesecs(:) = 0
    endif
    call pio_closefile(file_id)


    offset_time = 0._r8

    if ( (chlorine_loading_ymd > 0) .or. (chlorine_loading_tod > 0) ) then
       if (fixed) then
          call get_model_time( model_time )
          call convert_date( chlorine_loading_ymd, chlorine_loading_tod, time )
          offset_time = time - model_time
       else
          call endrun('chlorine_loading_init: cannont specify chlorine_loading_fixed_ymd ' &
                   // 'or chlorine_loading_fixed_tod if chlorine_loading_type is not FIXED' )
       endif
    endif

    call convert_dates( dates, datesecs, data_times )

    data_times = data_times - offset_time

    deallocate(dates)

    ! need to force data loading when the model starts at a time =/ 00:00:00.000
    ! -- may occur in restarts also
    call chlorine_loading_advance()
    initialized = .true.

  end subroutine chlorine_loading_init

!-----------------------------------------------------------------------
! Reads in the ETF data for the current date.  
!-----------------------------------------------------------------------
  subroutine chlorine_loading_advance( )
    use cam_pio_utils, only : cam_pio_openfile
    use pio, only : file_desc_t, pio_closefile, pio_inq_varid, pio_get_var
    use physconst,    only : cday
    use ioFileMod, only : getfil

    integer  :: year, month, day, sec
    integer  :: index, i
    integer  :: offset(1), count(1)
    logical  :: do_adv, read_data
    real(r8) :: time, delt
    type(file_desc_t) :: file_id
    integer :: data_vid, ierr
    character(len=256) :: filen


    if (.not.has_linoz_data) return
    if ( fixed .and. initialized ) return

    index = -1
    call get_model_time( time, year=year, month=month, day=day, seconds=sec )

    read_data = time > data_times(last_index) .or. .not.initialized

    if ( read_data ) then

       find_ndx: do i = last_index, ntimes
          if ( data_times(i) - time > 1.e-6_r8 ) then
             index = i-1
             exit find_ndx
          endif
       enddo find_ndx

       last_index = index+1

       if ( index < 1 ) then
          write(iulog,102) year,month,day,sec
          call endrun('chlorine_loading_advance: failed to read data from '//trim(chlorine_loading_file))
       endif

       ! get the surrounding time slices
       offset = (/ index /)
       count =  (/ nt /)

       call getfil( chlorine_loading_file, filen, 0 )
       call cam_pio_openfile( file_id, filen, 0 )
       if ( masterproc ) write(iulog,*)'chlorine_loading_advance: data file = ',trim(filen)
       ierr = pio_inq_varid( file_id, 'chlorine_loading', data_vid )
       ierr = pio_get_var( file_id, data_vid, offset, count, iloading )
       call pio_closefile(file_id)
    else
       index = last_index - 1
    endif

    delt = ( time - data_times(index) ) / ( data_times(index+1) - data_times(index) )

    ! this assures that FIXED data are b4b on restarts
    if ( fixed ) then
       delt = dble(int(delt*cday+.5_r8))/dble(cday)
    endif

    chlorine_loading = iloading(1) + delt*( iloading(2) - iloading(1) )
    
    if ( masterproc ) then
       write(iulog,101) year, month, day, sec, chlorine_loading
    endif

101 FORMAT('chlorine_loading_advance: date, loading : ',i4.4,'-',i2.2,'-',i2.2,'-',i5.5,',  ',f12.6)
102 FORMAT('chlorine_loading_advance: not able to find data for : ',i4.4,'-',i2.2,'-',i2.2,'-',i5.5)

  end subroutine chlorine_loading_advance
  
  !---------------------------------------------------------------------------
  ! private methods
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

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  subroutine convert_date( date, sec, time )

    integer,  intent(in)  :: date
    integer,  intent(in)  :: sec
    real(r8), intent(out) :: time

    integer :: dates(1), secs(1)
    real(r8) :: times(1)
    dates(1) = date
    secs(1) = sec
    call convert_dates( dates, secs, times )
    time = times(1)
  end subroutine convert_date

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine get_model_time( time, year, month, day, seconds )

    use time_manager, only: get_curr_date

    real(r8), intent(out) :: time
    integer, optional, intent(out) :: year, month, day, seconds

    integer  :: yr, mn, dy, sc, date

    call get_curr_date(yr, mn, dy, sc)
    date = yr*10000 + mn*100 + dy
    call convert_date( date, sc, time )

    if (present(year))    year = yr
    if (present(month))   month = mn
    if (present(day))     day = dy
    if (present(seconds)) seconds = sc

  end subroutine get_model_time

end module chlorine_loading_data

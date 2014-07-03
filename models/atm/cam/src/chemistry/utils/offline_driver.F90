!================================================================================
! offline unit driver utility module
!
!================================================================================
module offline_driver

  use shr_kind_mod, only: r8=>SHR_KIND_R8, cl=>SHR_KIND_CL
  use abortutils,   only: endrun
  use spmd_utils,   only: masterproc
  use cam_logfile,  only: iulog

  implicit none
  private
  save

  public :: offline_driver_init
  public :: offline_driver_end_of_data
  public :: offline_driver_run
  public :: offline_driver_dorun
  public :: offline_driver_done
  public :: offline_driver_readnl

  integer :: recno = 0
  logical :: offline_driver_dorun = .false.
  logical :: offline_driver_done = .false.

  character(len=cl) :: current_file = ' '
  character(len=cl) :: next_file = ' '
  character(len=cl) :: offline_driver_fileslist = ' '

  real(r8) :: time0

contains

!================================================================================
!================================================================================
  subroutine offline_driver_run( phys_state, pbuf2d, cam_out, cam_in, dtime )

    use physics_types,    only: physics_state
    use ppgrid,           only: begchunk, endchunk
    use drv_input_data,   only: dates, secs, ntimes, in_data_dtime
    use camsrfexch,       only: cam_out_t, cam_in_t     
    use physics_buffer,   only: physics_buffer_desc
    use time_manager,     only: get_step_size, timemgr_set_date_time
    use unit_driver,      only: unit_driver_run
    use tracer_data,      only: incr_filename
    use drv_input_data,   only: drv_input_data_open, drv_input_data_close

    type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
    type(cam_out_t),     intent(inout) :: cam_out(begchunk:endchunk)
    type(cam_in_t),      intent(inout) :: cam_in(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    real(r8) ,           intent(inout) :: dtime
    
    real(r8) :: time

    if (.not.offline_driver_dorun) return

    offline_driver_done = offline_driver_end_of_data()

    if (offline_driver_done) return

    dtime = get_step_size()

    recno = recno+1
    if ( recno > ntimes ) then
       call drv_input_data_close()
       current_file = next_file
       if ( current_file/='NOT_FOUND' ) then
          call drv_input_data_open( current_file )
          recno = 1
          next_file = incr_filename( current_file, filenames_list=offline_driver_fileslist, abort=.false.)
       else
          return
       endif
    endif
    call timemgr_set_date_time(dates(recno),secs(recno))
    
    time = get_time_val( dates(recno), secs(recno) )
    
    if ( time > time0 ) then
       in_data_dtime = time - time0
       time0 = time
    endif

    call unit_driver_run(phys_state, pbuf2d, cam_out, cam_in, dtime, recno)

  end subroutine offline_driver_run

!================================================================================
!================================================================================
  subroutine offline_driver_init
    use unit_driver,    only: unit_driver_init
    use drv_input_data, only: dates, secs, ntimes, in_data_dtime

    real(r8) :: time

    if (.not.offline_driver_dorun) return

    call unit_driver_init()

    time0 = get_time_val( dates(1), secs(1) )

    if (ntimes>1) then
       time = get_time_val( dates(2), secs(2) )
       in_data_dtime = time - time0
    endif

  endsubroutine offline_driver_init

!================================================================================
!================================================================================
  function offline_driver_end_of_data()

    use drv_input_data,     only: ntimes

    logical :: offline_driver_end_of_data

    offline_driver_end_of_data = recno>=ntimes .and.  trim(next_file)=='NOT_FOUND' 

  end function offline_driver_end_of_data

!=================================================================================
!=================================================================================
  subroutine offline_driver_readnl( nlfile )

    use namelist_utils,only: find_group_name
    use units,         only: getunit, freeunit
    use abortutils,    only: endrun
    use drv_input_data,only: drv_input_data_open, do_fdh
    use tracer_data,   only: incr_filename

    ! arguments
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! local vars
    integer :: unitn, ierr

    character(len=cl) :: offline_driver_infile = ' '
    logical :: offline_driver_do_fdh = .false.

    namelist /offline_driver_nl/ offline_driver_infile, offline_driver_fileslist, offline_driver_do_fdh

    unitn = getunit()
    open( unitn, file=trim(nlfile), status='old' )
    call find_group_name(unitn, 'offline_driver_nl', status=ierr)
    if (ierr == 0) then
       read(unitn, offline_driver_nl, iostat=ierr)
       if (ierr /= 0) then
          call endrun('offline_driver_readnl: ERROR reading namelist')
       end if
    end if
    close(unitn)
    call freeunit(unitn)

    current_file = 'NOT_FOUND'
    next_file = 'NOT_FOUND'
    offline_driver_dorun = .false.
    
    if ( len_trim(offline_driver_infile) > 0 ) then
       current_file = trim(offline_driver_infile)
    elseif ( len_trim(offline_driver_fileslist) > 0 ) then
       current_file = incr_filename( offline_driver_infile, filenames_list=offline_driver_fileslist )
    else 
       offline_driver_dorun = .false.
       return
    endif

    if ( trim(current_file)/='NOT_FOUND' .and. len_trim(current_file) > 0 ) then
       call drv_input_data_open( current_file )
       offline_driver_dorun = .true.
       if ( len_trim(offline_driver_fileslist) > 0 ) then
          next_file = incr_filename( current_file, filenames_list=offline_driver_fileslist, abort=.false.)
       endif
    endif

    do_fdh = offline_driver_do_fdh

  endsubroutine offline_driver_readnl

! private methods
!================================================================================
!================================================================================ 
  function get_time_val( date, secs ) result(time)
    use time_manager,  only: set_time_float_from_date
    use shr_const_mod, only: SHR_CONST_CDAY

    integer, intent(in) :: date, secs

    real(r8) :: time
    integer  :: year, month, day, sec

    year = date/10000
    month = (date-year*10000)/100
    day = date-year*10000-month*100
    sec = secs

    call set_time_float_from_date( time, year, month, day, sec )
    time = SHR_CONST_CDAY*time ! seconds

  end function get_time_val

end module offline_driver

!================================================================================
! utility module for driver input data
!================================================================================
module drv_input_data

  use shr_kind_mod, only: r8=>SHR_KIND_R8, cl=>SHR_KIND_CL, cs=>SHR_KIND_CS
  use abortutils,   only: endrun
  use spmd_utils,   only: masterproc
  use ppgrid,       only: pcols, pver, pverp, begchunk, endchunk
  use cam_logfile,  only: iulog
  use pio,          only: file_desc_t
  use time_manager, only: dtime

  implicit none
  private
  save

  public :: drv_input_data_open
  public :: drv_input_data_read
  public :: drv_input_data_close

  integer, public :: ntimes
  integer, public :: nlons
  integer, public :: nlats
  integer, public :: nlevs
  integer, public :: nilevs

  integer,  public, allocatable :: dates(:)
  integer,  public, allocatable :: secs(:)
  real(r8), public, allocatable :: times(:)

  real(r8), public, allocatable :: hyam(:)
  real(r8), public, allocatable :: hybm(:)
  real(r8), public :: P0
  real(r8), public :: in_data_dtime = 0._r8 ! sec
  logical,  public :: do_fdh

  integer, public :: data_dtime
  integer, public :: ref_ymd
  integer, public :: ref_tod
 
  character(len=cs), public :: calendar  ! Calendar type

  type(file_desc_t) :: piofile

  interface drv_input_data_read
    module procedure drv_input_data_read_2d
    module procedure drv_input_data_read_3d
  end interface

contains

!=================================================================================
!=================================================================================
  subroutine drv_input_data_open( infile )

    use ioFileMod,     only: getfil
    use cam_pio_utils, only: cam_pio_openfile
    use pio,           only: PIO_NOCLOBBER, pio_inq_dimid, pio_inq_dimlen
    use pio,           only: pio_inq_varid, pio_get_var, pio_get_att
    use shr_const_mod, only: SHR_CONST_CDAY
    use shr_string_mod,only: shr_string_toUpper

    use pmgrid, only: plon,plat,plev !,plevp,plnlv

    implicit none

    character(len=*), intent(in) :: infile


    character(len=cl) :: filen
    integer :: id, ierr
    integer :: nvars, vid, ndims

    character(len=cs) :: vname
    character(len=cs) :: time_units
    character(len=10) :: ref_date_str
    character(len=4) :: ref_year_str
    character(len=2) :: ref_mon_str
    character(len=2) :: ref_day_str
    character(len=2) :: ref_hr_str
    character(len=2) :: ref_min_str
    character(len=2) :: ref_sec_str
    character(len=8) :: ref_time_str
    integer :: ref_year, ref_mon, ref_day, ref_hr, ref_min, ref_sec

    !
    ! open file and get fileid
    !
    call getfil( infile, filen, 0 )
    call cam_pio_openfile( piofile, filen, PIO_NOCLOBBER)

    if(masterproc) write(iulog,*)'drv_input_data_open: opened: ',trim(filen)

    ierr = pio_inq_dimid( piofile, 'time', id )
    ierr = pio_inq_dimlen( piofile, id, ntimes )

    ierr = pio_inq_dimid( piofile, 'lon', id )
    ierr = pio_inq_dimlen( piofile, id, nlons )

    ierr = pio_inq_dimid( piofile, 'lat', id )
    ierr = pio_inq_dimlen( piofile, id, nlats )

    ierr = pio_inq_dimid( piofile, 'lev', id )
    ierr = pio_inq_dimlen( piofile, id, nlevs )

    ierr = pio_inq_dimid( piofile, 'ilev', id )
    ierr = pio_inq_dimlen( piofile, id, nilevs )

    allocate( dates(ntimes), secs(ntimes) )

    ierr = pio_inq_varid( piofile, 'date',  id  )
    ierr = pio_get_var( piofile, id, dates )

    ierr = pio_inq_varid( piofile, 'datesec',  id  )
    ierr = pio_get_var( piofile, id, secs )

    ierr = pio_inq_varid( piofile, 'time',  id  )
    ierr = pio_get_att( piofile, id, 'units', time_units )
    ierr = pio_get_att( piofile, id, 'calendar', calendar )

    calendar = shr_string_toUpper(calendar)

    ref_date_str = time_units(12:21)
    ref_year_str = ref_date_str(1:4)
    ref_mon_str = ref_date_str(6:7)
    ref_day_str = ref_date_str(9:10)

    read(ref_year_str,fmt='(I4)') ref_year
    read(ref_mon_str, fmt='(I2)') ref_mon
    read(ref_day_str, fmt='(I2)') ref_day

    ref_ymd = ref_year*10000 + ref_mon*100 + ref_day

    ref_time_str = time_units(23:)
    ref_hr_str = ref_time_str(1:2)
    ref_min_str = ref_time_str(4:5)
    ref_sec_str = ref_time_str(7:8)

    read(ref_hr_str, fmt='(I2)') ref_hr
    read(ref_min_str,fmt='(I2)') ref_min
    read(ref_sec_str,fmt='(I2)') ref_sec

    ref_tod = ref_hr*3600 + ref_min*60 + ref_sec

    ierr = pio_inq_varid( piofile, 'mdt',  id  )
    ierr = pio_get_var( piofile, id, data_dtime )

    if (plon /= nlons) then
      call endrun('drv_input_data_open: plon /= nlons')
    endif
    if (plat /= nlats) then
      call endrun('drv_input_data_open: plat /= nlats')
    endif
    if (plev /= nlevs) then
      call endrun('drv_input_data_open: plev /= nlevs')
    endif

    allocate(hyam(nlevs))
    ierr = pio_inq_varid( piofile, 'hyam',  id  )
    ierr = pio_get_var( piofile, id, hyam )
    
    allocate(hybm(nlevs))
    ierr = pio_inq_varid( piofile, 'hybm',  id  )
    ierr = pio_get_var( piofile, id, hybm )

    ierr = pio_inq_varid( piofile, 'P0',  id  )
    ierr = pio_get_var( piofile, id, p0 )

    if (ntimes>1) then
       if ( .not. (data_dtime == dtime)) then
          write( iulog, * )  'drv_input_data_open: data freq does not match dtime... use dtime = ',data_dtime
          call endrun( 'drv_input_data_open: data freq does not match dtime.')
       endif
    endif

  end subroutine drv_input_data_open

!================================================================================================
!================================================================================================
  subroutine drv_input_data_close
    use pio, only: pio_closefile
    implicit none

    deallocate( dates, secs, hyam, hybm )

    call pio_closefile( piofile )

  end subroutine drv_input_data_close

  !=================================================================================
  !=================================================================================
  function drv_input_data_read_2d( fldname, recno, abort ) result(field_array)
    use ncdio_atm,        only: infld

    implicit none

    character(len=*), intent(in) :: fldname
    integer,          intent(in) :: recno
    logical, optional,intent(in) :: abort

    logical  :: found, abort_run
    real(r8) :: field_array(pcols,begchunk:endchunk)

    abort_run = .false.
    if (present(abort)) then
       abort_run = abort
    endif

    call infld( fldname, piofile, 'lon','lat', 1,pcols, begchunk,endchunk, &
                field_array, found, grid_map='PHYS',timelevel=recno)

    if (.not.found) then
       if ( abort_run ) then
          call endrun('drv_input_data_read_2d: did not find '// trim(fldname))
       else
          if (masterproc) write( iulog, * )  'drv_input_data_read_2d: ' // trim(fldname) // ' set to zero '
          field_array = 0._r8
       endif
    endif

  endfunction drv_input_data_read_2d

  !=================================================================================
  !=================================================================================
  function drv_input_data_read_3d( fldname, vertname, vertsize, recno, abort ) result(field_array)
    use ncdio_atm,        only: infld
    implicit none

    character(len=*), intent(in) :: fldname
    character(len=*), intent(in) :: vertname
    integer,          intent(in) :: vertsize
    integer,          intent(in) :: recno
    logical, optional,intent(in) :: abort

    logical  :: found, abort_run
    real(r8) :: field_array(pcols,vertsize,begchunk:endchunk)

    real(r8), allocatable :: tmp_array(:,:,:)
    
    abort_run = .false.
    if (present(abort)) then
       abort_run = abort
    endif

    call infld( fldname, piofile, 'lon',vertname,'lat', 1,pcols, 1,vertsize, begchunk,endchunk, &
                field_array, found, grid_map='PHYS',timelevel=recno)
    
    if (.not.found) then
       if ( abort_run ) then
          call endrun('drv_input_data_read_3d: did not find '// trim(fldname))
       else
          if (masterproc) write( iulog, * )  'drv_input_data_read_3d: ' // trim(fldname) // ' set to zero '
          field_array = 0._r8
       endif
    endif

  endfunction drv_input_data_read_3d

end module drv_input_data

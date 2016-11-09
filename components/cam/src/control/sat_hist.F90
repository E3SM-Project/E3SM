!-------------------------------------------------------------------------------
! Outputs history field columns as specified by a satellite track data file
!
! Created by Francis Vitt -- 17 Sep 2010
!-------------------------------------------------------------------------------
module sat_hist

  use perf_mod,      only: t_startf, t_stopf
  use shr_kind_mod,  only: r8 => shr_kind_r8
  use cam_logfile,   only: iulog
  use ppgrid,        only: pcols, pver, begchunk, endchunk
  use cam_history_support, only: fieldname_lenp2, max_string_len, ptapes
  use spmd_utils,    only: masterproc, iam
  use cam_abortutils,    only: endrun

  use pio,           only: file_desc_t, iosystem_desc_t, iosystem_desc_t, var_desc_t, io_desc_t
  use pio,           only: pio_openfile, pio_redef, pio_enddef, pio_inq_dimid, pio_inq_varid, pio_seterrorhandling, pio_def_var
  use pio,           only: pio_inq_dimlen, pio_get_att, pio_put_att, pio_get_var, pio_put_var, pio_write_darray
  use pio,           only: pio_real, pio_int, pio_double
  use pio,           only: PIO_WRITE,PIO_NOWRITE, PIO_NOERR, PIO_BCAST_ERROR, PIO_INTERNAL_ERROR, PIO_Rearr_box, PIO_GLOBAL
  use spmd_utils,    only: mpicom
#ifdef SPMD
  use mpishorthand,  only: mpichar, mpiint
#endif
   use physconst, only: pi 
  
  implicit none

  private
  save

  public :: sat_hist_readnl
  public :: sat_hist_init
  public :: sat_hist_write
  public :: sat_hist_define
  public :: is_satfile

  character(len=max_string_len)  :: sathist_track_infile
  type(file_desc_t) :: infile

  integer :: half_step
  logical :: has_sat_hist = .false.

  integer :: sathist_nclosest
  integer :: sathist_ntimestep

  real(r8), allocatable :: obs_lats(:)
  real(r8), allocatable :: obs_lons(:)

  logical  :: doy_format
  real(r8) :: first_datetime
  real(r8) :: last_datetime
  integer  :: last_start_index
  integer  :: time_ndx
  integer  :: t_buffer_size
  integer, allocatable :: date_buffer(:), time_buffer(:)
  integer :: sat_tape_num=ptapes-1

  
  ! input file
  integer :: n_profiles
  integer :: time_vid, date_vid, lat_vid, lon_vid, instr_vid, orbit_vid, prof_vid, zenith_vid

  integer :: in_julian_vid
  integer :: in_localtime_vid
  integer :: in_doy_vid
  integer :: in_occ_type_vid

  integer :: in_start_col


  ! output file
  type(var_desc_t) :: out_latid, out_lonid, out_dstid, out_instrid, out_zenithid, out_orbid, out_profid
  type(var_desc_t) :: out_instr_lat_vid, out_instr_lon_vid
  type(var_desc_t) :: out_obs_date_vid, out_obs_time_vid
  type(var_desc_t) :: out_julian_vid
  type(var_desc_t) :: out_localtime_vid
  type(var_desc_t) :: out_doy_vid
  type(var_desc_t) :: out_occ_type_vid

  logical, parameter :: debug = .false.

  real(r8), parameter :: rad2deg = 180._r8/pi            ! degrees per radian

contains
  
!-------------------------------------------------------------------------------

  logical function is_satfile (file_index)
    integer, intent(in) :: file_index ! index of file in question
    is_satfile = file_index == sat_tape_num
  end function is_satfile

!-------------------------------------------------------------------------------
  subroutine sat_hist_readnl(nlfile, hfilename_spec, mfilt, fincl, nhtfrq, avgflag_pertape)
    
    use namelist_utils,      only: find_group_name
    use units,               only: getunit, freeunit
    use cam_history_support, only: pflds
    use cam_instance,        only: inst_suffix

    implicit none

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
    character(len=*), intent(inout) :: hfilename_spec(:)
    character(len=*), intent(inout) :: fincl(:,:)
    character(len=1), intent(inout) :: avgflag_pertape(:)
    integer,          intent(inout) :: mfilt(:), nhtfrq(:)
    
    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'sat_hist_readnl'
    integer :: f, fcnt

    character(len=fieldname_lenp2) :: sathist_fincl(pflds)
    character(len=max_string_len)  :: sathist_hfilename_spec
    integer :: sathist_mfilt, sat_tape_num

    namelist /satellite_options_nl/ sathist_track_infile, sathist_hfilename_spec, sathist_fincl, &
         sathist_mfilt, sathist_nclosest, sathist_ntimestep

    ! set defaults

    sathist_track_infile = ' '
    sathist_hfilename_spec = '%c.cam' // trim(inst_suffix) // '.hs.%y-%m-%d-%s.nc'
    sathist_fincl(:) = ' '
    sathist_mfilt = 100000
    sathist_nclosest = 1
    sathist_ntimestep = 1

    !read namelist options

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'satellite_options_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, satellite_options_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    ! broadcast the options to all MPI tasks
    call mpibcast(sathist_track_infile,   len(sathist_track_infile),   mpichar, 0, mpicom)
    call mpibcast(sathist_hfilename_spec, len(sathist_hfilename_spec), mpichar, 0, mpicom)
    call mpibcast(sathist_fincl,          pflds*len(sathist_fincl(1)), mpichar, 0, mpicom)
    call mpibcast(sathist_mfilt,          1,                           mpiint,  0, mpicom)
    call mpibcast(sathist_nclosest,       1,                           mpiint,  0, mpicom)
    call mpibcast(sathist_ntimestep,      1,                           mpiint,  0, mpicom)
#endif

    has_sat_hist = len_trim(sathist_track_infile) > 0

    if (.not.has_sat_hist) return

     sat_tape_num=ptapes-1
     hfilename_spec(sat_tape_num) = sathist_hfilename_spec
     mfilt(sat_tape_num) = sathist_mfilt
     fcnt=0
     do f=1, pflds
        fincl(f,sat_tape_num) = sathist_fincl(f)
        if(len_trim(sathist_fincl(f)) > 0) then
           fcnt=fcnt+1
        end if
     enddo
     
     nhtfrq(sat_tape_num) = 1
     avgflag_pertape(sat_tape_num) = 'I'

     if(masterproc) then
        write(iulog,*) 'sathist_track_infile: ',trim(sathist_track_infile)
        write(iulog,*) 'sathist_hfilename_spec: ',trim(sathist_hfilename_spec)
        write(iulog,*) 'sathist_fincl: ',(trim(sathist_fincl(f))//' ', f=1,fcnt)
        write(iulog,*) 'max columns per file sathist_mfilt: ',sathist_mfilt
        write(iulog,*) 'sathist_nclosest: ',sathist_nclosest
        write(iulog,*) 'sathist_ntimestep: ',sathist_ntimestep
     end if

   end subroutine sat_hist_readnl

  
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine sat_hist_init
    use cam_pio_utils, only: cam_pio_openfile
    use ioFileMod,     only: getfil
    use spmd_utils,    only: npes
    use time_manager,  only: get_step_size
    use string_utils,  only: to_lower, GLC

    implicit none

    character(len=max_string_len)  :: locfn       ! Local filename
    integer :: ierr, dimid, i

    character(len=128) :: date_format

    if (.not.has_sat_hist) return

    call getfil (sathist_track_infile, locfn)
    call cam_pio_openfile(infile, locfn, PIO_NOWRITE)

    ierr = pio_inq_dimid(infile,'profs',dimid)
    ierr = pio_inq_dimlen(infile, dimid, n_profiles)

    ierr = pio_inq_varid( infile, 'time', time_vid )
    ierr = pio_inq_varid( infile, 'date', date_vid )

    ierr = pio_get_att( infile, date_vid, 'long_name', date_format)
    date_format = to_lower(trim( date_format(:GLC(date_format))))

    if ( index( date_format, 'yyyymmdd') > 0 ) then
       doy_format = .false.
    else if  ( index( date_format, 'yyyyddd') > 0 ) then
       doy_format = .true.
    else
       call endrun('sat_hist_init: date_format not recognized : '//trim(date_format))
    endif

    ierr = pio_inq_varid( infile, 'lat', lat_vid )
    ierr = pio_inq_varid( infile, 'lon', lon_vid )

    call pio_seterrorhandling(infile, PIO_BCAST_ERROR)
    ierr = pio_inq_varid( infile, 'instr_num', instr_vid )
    if(ierr/=PIO_NOERR) instr_vid=-1

    ierr = pio_inq_varid( infile, 'orbit_num', orbit_vid )
    if(ierr/=PIO_NOERR) orbit_vid=-1

    ierr = pio_inq_varid( infile, 'prof_num',  prof_vid )
    if(ierr/=PIO_NOERR) prof_vid=-1

    ierr = pio_inq_varid( infile, 'instr_sza', zenith_vid )
    if(ierr/=PIO_NOERR) zenith_vid=-1

    ierr = pio_inq_varid( infile, 'julian', in_julian_vid )
    if(ierr/=PIO_NOERR) in_julian_vid=-1

    ierr = pio_inq_varid( infile, 'local_time', in_localtime_vid )
    if(ierr/=PIO_NOERR) in_localtime_vid=-1

    ierr = pio_inq_varid( infile, 'doy', in_doy_vid )
    if(ierr/=PIO_NOERR) in_doy_vid=-1

    ierr = pio_inq_varid( infile, 'occ_type', in_occ_type_vid )
    if(ierr/=PIO_NOERR) in_occ_type_vid=-1

    call pio_seterrorhandling(infile, PIO_INTERNAL_ERROR)

    call read_datetime( first_datetime, 1 )
    call read_datetime( last_datetime, n_profiles )
    last_start_index = -1
    t_buffer_size = min(1000,n_profiles)
    allocate( date_buffer(t_buffer_size), time_buffer(t_buffer_size) )
    if (masterproc) write(iulog,*) "sathist_init:", n_profiles, first_datetime, last_datetime
    if ( last_datetime<first_datetime ) then
       call endrun('sat_hist_init: satellite track file has invalid date time info')
    endif

    time_ndx = 1
    half_step = get_step_size()*0.5_r8

  end subroutine sat_hist_init

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine read_datetime( datetime, index )

    real(r8), intent( out ) :: datetime
    integer,  intent( in )  :: index

    integer :: ierr
    integer :: cnt(1)
    integer :: start(1)
    integer :: date(1), time(1)

    cnt = (/ 1 /)
    start = (/index/)

    ierr = pio_get_var( infile, time_vid, start, cnt, time )
    ierr = pio_get_var( infile, date_vid, start, cnt, date )
    
    datetime = convert_date_time( date(1),time(1) )

  end subroutine read_datetime

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine read_buffered_datetime( datetime, index )

    real(r8), intent( out ) :: datetime
    integer,  intent( in )  :: index

    integer :: ii

    integer :: ierr
    integer :: cnt
    integer :: start
    integer :: date, time
    
    ! If the request is outside of the buffer then reload the buffer.
    if ((last_start_index == -1) .or. (index < last_start_index) &
         .or. (index >= (last_start_index + t_buffer_size))) then

       start = (index - 1) / t_buffer_size * t_buffer_size + 1
       if ( start+t_buffer_size-1 <= n_profiles ) then
          cnt = t_buffer_size 
       else
          cnt = n_profiles-start+1
       endif
       ierr = pio_get_var( infile, time_vid, (/ start /), (/ cnt /), time_buffer(1:cnt) )
       ierr = pio_get_var( infile, date_vid, (/ start /), (/ cnt /), date_buffer(1:cnt) )

       last_start_index = start
    endif

    ii = mod( index - 1, t_buffer_size ) + 1
    time = time_buffer(ii)
    date = date_buffer(ii)
    datetime = convert_date_time( date,time )

  end subroutine read_buffered_datetime

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  function convert_date_time( date,time )
    use time_manager, only: set_time_float_from_date

    integer, intent(in) :: date,time
    real(r8) :: convert_date_time

    real(r8) :: datetime
    integer :: yr, doy, mon, dom

    if ( doy_format ) then
       yr = date/1000
       doy = date - yr*1000
       call set_time_float_from_date( datetime, yr, 1, doy, time )
    else 
       yr = date/10000
       mon = (date - yr*10000)/100
       dom = date - yr*10000 - mon*100
       call set_time_float_from_date( datetime, yr, mon, dom, time )
    endif
    convert_date_time = datetime

  end function convert_date_time
!-------------------------------------------------------------------------------
  subroutine sat_hist_define(outfile)
    use pio, only : pio_inquire
    type(file_desc_t), intent(inout) :: outfile

    integer :: coldim
    integer :: ierr
    
    ierr = pio_inquire(outfile, unlimitedDimId=coldim)

    call pio_seterrorhandling(outfile, PIO_BCAST_ERROR)
    ierr = define_var( 'instr_lat', coldim, infile, lat_vid,  outfile, out_instr_lat_vid )
    ierr = define_var( 'instr_lon', coldim, infile, lon_vid,  outfile, out_instr_lon_vid )
    ierr = define_var( 'obs_time', coldim, infile, time_vid,  outfile, out_obs_time_vid )
    ierr = define_var( 'obs_date', coldim, infile, date_vid,  outfile, out_obs_date_vid )

    ierr = pio_inq_varid( outfile, 'distance', out_dstid )
    if (ierr /= PIO_NOERR) then
       ierr = pio_def_var  ( outfile, 'distance', PIO_REAL, (/coldim/), out_dstid )
       ierr = pio_put_att  ( outfile, out_dstid, "long_name", "distance from midpoint to observation")
       ierr = pio_put_att  ( outfile, out_dstid, "units", "km")
    end if

    if (orbit_vid>0) then
       ierr = define_var( 'orbit_num', coldim, infile, orbit_vid,  outfile, out_orbid )
    endif
    if (prof_vid>0) then
       ierr = define_var( 'prof_num', coldim, infile, prof_vid,  outfile, out_profid )
    endif
    if (instr_vid>0) then
       ierr = define_var( 'instr_num', coldim, infile, instr_vid,  outfile, out_instrid )
    endif
    if (zenith_vid>0) then
       ierr = define_var( 'instr_sza', coldim, infile, zenith_vid,  outfile, out_zenithid )
    endif
    if (in_occ_type_vid>0) then
       ierr = define_var( 'occ_type', coldim, infile, in_occ_type_vid,  outfile, out_occ_type_vid )
    endif
    if (in_julian_vid>0) then
       ierr = define_var( 'julian', coldim, infile, in_julian_vid,  outfile, out_julian_vid )
    endif
    if (in_localtime_vid>0) then
       ierr = define_var( 'local_time', coldim, infile, in_localtime_vid,  outfile, out_localtime_vid )
    endif
    if (in_doy_vid>0) then
       ierr = define_var( 'doy', coldim, infile, in_doy_vid,  outfile, out_doy_vid )
    endif

    call pio_seterrorhandling(outfile, PIO_INTERNAL_ERROR)
    ierr=pio_put_att (outfile, PIO_GLOBAL, 'satellite_track_file', sathist_track_infile)
  end subroutine sat_hist_define


!-------------------------------------------------------------------------------
  subroutine sat_hist_write( tape , nflds, nfils)

    use ppgrid,   only : pcols, begchunk, endchunk
    use phys_grid, only: phys_decomp
    use dyn_grid, only: dyn_decomp
    use cam_history_support, only : active_entry
    use pio, only : pio_file_is_open
    implicit none
    type(active_entry) :: tape
    integer, intent(in) :: nflds
    integer, intent(inout) :: nfils

    integer :: t, f, i, ncols, nocols    
    integer :: ierr

    integer, allocatable :: col_ndxs(:)
    integer, allocatable :: chk_ndxs(:)
    integer, allocatable :: fdyn_ndxs(:)
    integer, allocatable :: ldyn_ndxs(:)
    integer, allocatable :: phs_owners(:)
    integer, allocatable :: dyn_owners(:)
    real(r8),allocatable :: mlats(:)
    real(r8),allocatable :: mlons(:)
    real(r8),allocatable :: phs_dists(:)

    integer :: coldim

    integer :: io_type
    logical :: has_dyn_flds

    if (.not.has_sat_hist) return

    call read_next_position( ncols )

    if ( ncols < 1 ) return

    call t_startf ('sat_hist_write')

    ! The n closest columns to the observation will be output,
    ! so increase the size of the columns used for output/
    nocols = ncols * sathist_nclosest

    allocate( col_ndxs(nocols) )
    allocate( chk_ndxs(nocols) )
    allocate( fdyn_ndxs(nocols) )
    allocate( ldyn_ndxs(nocols) )
    allocate( phs_owners(nocols) )
    allocate( dyn_owners(nocols) )
    allocate( mlats(nocols) )
    allocate( mlons(nocols) )
    allocate( phs_dists(nocols) )

    has_dyn_flds = .false.
    dyn_flds_loop: do f=1,nflds
       if ( tape%hlist(f)%field%decomp_type == dyn_decomp ) then
          has_dyn_flds = .true.
          exit dyn_flds_loop
       endif
    enddo dyn_flds_loop

    call get_indices( obs_lats, obs_lons, ncols, nocols, has_dyn_flds, col_ndxs, chk_ndxs, &
         fdyn_ndxs, ldyn_ndxs, phs_owners, dyn_owners, mlats, mlons, phs_dists )

    if ( .not. pio_file_is_open(tape%File) ) then
       call endrun('sat file not open')
    endif

    ierr = pio_inq_dimid(tape%File,'ncol',coldim )
    
    ierr = pio_inq_varid(tape%File, 'lat', out_latid )
    ierr = pio_inq_varid(tape%File, 'lon', out_lonid )
    ierr = pio_inq_varid(tape%File, 'distance', out_dstid )

    call write_record_coord( tape, mlats(:), mlons(:), phs_dists(:), ncols, nfils )

    do f=1,nflds

       select case (tape%hlist(f)%field%decomp_type)
       case (phys_decomp)
          call dump_columns(tape%File, tape%hlist(f), nocols, nfils, col_ndxs(:), chk_ndxs(:), phs_owners(:) )
       case (dyn_decomp)
          call dump_columns(tape%File, tape%hlist(f), nocols, nfils, fdyn_ndxs(:), ldyn_ndxs(:), dyn_owners(:) )
       end select

    enddo

    deallocate( col_ndxs, chk_ndxs, fdyn_ndxs, ldyn_ndxs, phs_owners, dyn_owners )
    deallocate( mlons, mlats, phs_dists )
    deallocate( obs_lons, obs_lats )

    nfils = nfils + nocols

    call t_stopf ('sat_hist_write')

  end subroutine sat_hist_write

!-------------------------------------------------------------------------------
  subroutine dump_columns( File, hitem, ncols, nfils, fdims, ldims, owners  )
    use cam_history_support,  only: field_info, hentry, hist_coords, fillvalue
    use pionfwrite_mod, only: write_nf
    use pio,            only: pio_initdecomp, pio_freedecomp, pio_setframe, pio_iam_iotask, pio_setdebuglevel, pio_offset_kind

    type(File_desc_t),intent(inout)  :: File
    type(hentry),     intent(in), target     :: hitem
    integer,          intent(in)     :: ncols
    integer,          intent(in)     :: nfils
    integer,          intent(in)     :: fdims(:)
    integer,          intent(in)     :: ldims(:)
    integer,          intent(in)     :: owners(:)

    type(field_info), pointer :: field
    type(var_desc_t) :: vardesc
    type(iosystem_desc_t), pointer :: sat_iosystem
    type(io_desc_t) :: iodesc
    integer :: t, ierr, ndims
    integer, allocatable :: dimlens(:)

    real(r8), allocatable :: buf(:)
    integer,  allocatable :: dof(:)
    integer :: i,k, cnt

    call t_startf ('sat_hist::dump_columns')

    sat_iosystem => File%iosystem
    field => hitem%field
    vardesc = hitem%varid(1)


    ndims=1
    if(associated(field%mdims)) then
       ndims = size(field%mdims)+1
    else if(field%numlev>1) then
       ndims=2
    end if
    allocate(dimlens(ndims))
    dimlens(ndims)=ncols
    if(ndims>2) then
       do i=1,ndims-1
          dimlens(i)=hist_coords(field%mdims(i))%dimsize
       enddo
    else if(field%numlev>1) then
       dimlens(1) = field%numlev
    end if
   

    allocate( buf( product(dimlens) ) )
    allocate( dof( product(dimlens) ) )

    cnt = 0
    buf = fillvalue
    dof = 0

    do i = 1,ncols
       do k = 1,field%numlev
          cnt = cnt+1
          if ( iam == owners(i) ) then
             buf(cnt) = hitem%hbuf( fdims(i), k, ldims(i) )
             dof(cnt) = cnt
          endif
       enddo
    enddo

    call pio_setframe(File, vardesc, int(-1,kind=PIO_OFFSET_KIND))

    call pio_initdecomp(sat_iosystem, pio_double, dimlens, dof, iodesc )

    call pio_setframe(File, vardesc, int(nfils,kind=PIO_OFFSET_KIND))

    call pio_write_darray(File, vardesc, iodesc, buf, ierr, fillval=fillvalue)

    call pio_freedecomp(sat_iosystem, iodesc)

    deallocate( buf )
    deallocate( dof )
    deallocate( dimlens )

    call t_stopf ('sat_hist::dump_columns')

  end subroutine dump_columns

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine read_next_position( ncols )
    use time_manager, only: get_curr_date, get_prev_date
    use time_manager, only: set_time_float_from_date
   
    implicit none

    integer,  intent(out) :: ncols

    integer :: ierr
    integer :: yr, mon, day, tod
    real(r8) :: begdatetime, enddatetime
    integer :: beg_ndx, end_ndx, i

    real(r8) :: datetime

    call get_curr_date(yr, mon, day, tod)
    call set_time_float_from_date(begdatetime, yr, mon, day, tod-half_step*sathist_ntimestep)
    call set_time_float_from_date(enddatetime, yr, mon, day, tod+half_step*sathist_ntimestep)

    ncols = 0

    if ( first_datetime > enddatetime ) then
       if (masterproc) write(iulog,'(a,2f16.6)') &
            'sat_hist->read_next_position: all of the satellite date times are after the time window', first_datetime, enddatetime
       return
    endif
    if ( last_datetime < begdatetime ) then
       if (masterproc) write(iulog,'(a,2f16.6)') &
            'sat_hist->read_next_position: all of the satellite date times are before the time window', begdatetime, last_datetime
       return
    endif

    call t_startf ('sat_hist::read_next_position')

    beg_ndx = -99
    end_ndx = -99

    bnds_loop: do i = time_ndx,n_profiles

       call read_buffered_datetime( datetime, i )

       if ( datetime>begdatetime .and. beg_ndx<0 ) beg_ndx = i
       if ( datetime>enddatetime ) exit bnds_loop
       end_ndx = i

    enddo bnds_loop

    if (beg_ndx == -99 .and. end_ndx== -99) then
       if (masterproc) write(iulog,'(a)')  'sat_hist->read_next_position: must be beyond last position -- returning.'
       return
    endif

    ! Advance the search forward, but because of ntimesteps, it is possible
    ! for observations used here to be used again. However, we should not go
    ! back before the previous beginning time.
    if (beg_ndx>0) time_ndx = beg_ndx

    ncols = end_ndx-beg_ndx+1

    if (ncols > 0) then
       allocate( obs_lats(ncols), obs_lons(ncols) )
       in_start_col = beg_ndx

       ierr = pio_get_var( infile, lat_vid, (/beg_ndx/), (/ncols/), obs_lats )
       ierr = pio_get_var( infile, lon_vid, (/beg_ndx/), (/ncols/), obs_lons )

    endif

    call t_stopf ('sat_hist::read_next_position')
  end subroutine read_next_position
  
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine write_record_coord( tape, mod_lats, mod_lons, mod_dists, ncols, nfils )

    use time_manager,  only: get_nstep, get_curr_date, get_curr_time
    use cam_history_support, only : active_entry
    implicit none
    type(active_entry), intent(inout) :: tape

    integer,  intent(in) :: ncols
    real(r8), intent(in) :: mod_lats(ncols * sathist_nclosest)
    real(r8), intent(in) :: mod_lons(ncols * sathist_nclosest)
    real(r8), intent(in) :: mod_dists(ncols * sathist_nclosest)
    integer,  intent(in) ::  nfils

    integer :: t, ierr, i
    integer :: yr, mon, day      ! year, month, and day components of a date
    integer :: nstep             ! current timestep number
    integer :: ncdate            ! current date in integer format [yyyymmdd]
    integer :: ncsec             ! current time of day [seconds]
    integer :: ndcur             ! day component of current time
    integer :: nscur             ! seconds component of current time
    real(r8) :: time             ! current time
    integer, allocatable  :: itmp(:)
    real(r8), allocatable :: rtmp(:)
    real(r8), allocatable :: out_lats(:)
    real(r8), allocatable :: out_lons(:)

    call t_startf ('sat_hist::write_record_coord')

    nstep = get_nstep()
    call get_curr_date(yr, mon, day, ncsec)
    ncdate = yr*10000 + mon*100 + day
    call get_curr_time(ndcur, nscur)


    time = ndcur + nscur/86400._r8

    allocate( itmp(ncols * sathist_nclosest) )
    allocate( rtmp(ncols * sathist_nclosest) )
    
    itmp(:) = ncdate
    ierr = pio_put_var(tape%File, tape%dateid,(/nfils/), (/ncols * sathist_nclosest/),itmp)
    itmp(:) = ncsec
    ierr = pio_put_var(tape%File, tape%datesecid,(/nfils/),(/ncols * sathist_nclosest/),itmp)
    rtmp(:) = time
    ierr = pio_put_var(tape%File, tape%timeid, (/nfils/),(/ncols * sathist_nclosest/),rtmp)
    
    deallocate(itmp)
    deallocate(rtmp)

    ! output model column coordinates
    ierr = pio_put_var(tape%File, out_latid, (/nfils/),(/ncols * sathist_nclosest/), mod_lats)
    ierr = pio_put_var(tape%File, out_lonid, (/nfils/),(/ncols * sathist_nclosest/), mod_lons)
    ierr = pio_put_var(tape%File, out_dstid, (/nfils/),(/ncols * sathist_nclosest/), mod_dists / 1000._r8)
    
    ! output instrument location
    allocate( out_lats(ncols * sathist_nclosest) )
    allocate( out_lons(ncols * sathist_nclosest) )
    
    do i = 1, ncols
      out_lats(((i-1)*sathist_nclosest)+1 : (i*sathist_nclosest)) = obs_lats(i)
      out_lons(((i-1)*sathist_nclosest)+1 : (i*sathist_nclosest)) = obs_lons(i)
    enddo

    ierr = pio_put_var(tape%File, out_instr_lat_vid, (/nfils/),(/ncols * sathist_nclosest/), out_lats)
    ierr = pio_put_var(tape%File, out_instr_lon_vid, (/nfils/),(/ncols * sathist_nclosest/), out_lons)

    deallocate(out_lats)
    deallocate(out_lons)
    
    
    ierr = copy_data( infile, date_vid, tape%File, out_obs_date_vid, in_start_col, nfils, ncols )
    ierr = copy_data( infile, time_vid, tape%File, out_obs_time_vid, in_start_col, nfils, ncols )
    
    ! output observation identifiers
    if (instr_vid>0) then
       ierr = copy_data( infile, instr_vid, tape%File, out_instrid, in_start_col, nfils, ncols )
    endif
    if (orbit_vid>0) then
       ierr = copy_data( infile, orbit_vid, tape%File, out_orbid, in_start_col, nfils, ncols )
    endif
    if (prof_vid>0) then
       ierr = copy_data( infile, prof_vid, tape%File, out_profid, in_start_col, nfils, ncols )
    endif
    if (zenith_vid>0) then
       ierr = copy_data( infile, zenith_vid, tape%File, out_zenithid, in_start_col, nfils, ncols )
    endif
    if (in_julian_vid>0) then
       ierr = copy_data( infile, in_julian_vid, tape%File, out_julian_vid, in_start_col, nfils, ncols )
    endif
    if (in_occ_type_vid>0) then
       ierr = copy_data( infile, in_occ_type_vid, tape%File, out_occ_type_vid, in_start_col, nfils, ncols )
    endif
    if (in_localtime_vid>0) then
       ierr = copy_data( infile, in_localtime_vid, tape%File, out_localtime_vid, in_start_col, nfils, ncols )
    endif
    if (in_doy_vid>0) then
       ierr = copy_data( infile, in_doy_vid, tape%File, out_doy_vid, in_start_col, nfils, ncols )
    endif

    call t_stopf ('sat_hist::write_record_coord')
  end subroutine write_record_coord

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

  subroutine get_indices( lats, lons, ncols, nocols, has_dyn_flds, col_ndxs, chk_ndxs, &
       fdyn_ndxs, ldyn_ndxs, phs_owners, dyn_owners, mlats, mlons, phs_dists )

    use dyn_grid, only : dyn_grid_get_colndx
    use phys_grid, only: get_rlat_p, get_rlon_p

    integer,  intent(in)  :: ncols
    real(r8), intent(in)  :: lats(ncols)
    real(r8), intent(in)  :: lons(ncols)
    integer,  intent(in)  :: nocols
    logical,  intent(in)  :: has_dyn_flds
    integer,  intent(out) :: col_ndxs(nocols)
    integer,  intent(out) :: chk_ndxs(nocols)
    integer,  intent(out) :: fdyn_ndxs(nocols)
    integer,  intent(out) :: ldyn_ndxs(nocols)
    integer,  intent(out) :: phs_owners(nocols)
    integer,  intent(out) :: dyn_owners(nocols)
    real(r8), intent(out) :: mlats(nocols)
    real(r8), intent(out) :: mlons(nocols)
    real(r8), intent(out) :: phs_dists(nocols)

    integer :: i, j, ndx
    real(r8) :: lat, lon
    
    integer,  allocatable :: ichks(:),icols(:),idyn1s(:),idyn2s(:), iphs_owners(:), idyn_owners(:)
    real(r8), allocatable :: rlats(:), rlons(:), plats(:), plons(:), iphs_dists(:)

    integer :: gcols(sathist_nclosest)

    call t_startf ('sat_hist::get_indices')

    allocate(ichks(sathist_nclosest),icols(sathist_nclosest),idyn1s(sathist_nclosest), &
         idyn2s(sathist_nclosest),iphs_owners(sathist_nclosest),idyn_owners(sathist_nclosest))
    allocate(rlats(sathist_nclosest), rlons(sathist_nclosest), plats(sathist_nclosest), &
         plons(sathist_nclosest), iphs_dists(sathist_nclosest) )

    col_ndxs = -1
    chk_ndxs = -1
    fdyn_ndxs = -1
    ldyn_ndxs = -1
    phs_owners = -1
    dyn_owners = -1
    phs_dists = -1

    ndx = 0
    do i = 1,ncols

       lat = lats(i)
       lon = lons(i)

       if ( lon >= 360._r8) then
         lon = lon-360._r8
       endif
       if ( lon < 0._r8) then
         lon = lon+360._r8
       endif
       if (lat<-90._r8 .or. lat>90._r8) then
          write(iulog,*) 'sat_hist::get_indices lat = ',lat
          call endrun('sat_hist::get_indices : lat must be between -90 and 90 degrees (-90<=lat<=90)')
       endif
       if (lon<0._r8 .or. lon>=360._r8) then
          write(iulog,*) 'sat_hist::get_indices lon = ',lon
          call endrun('sat_hist::get_indices : lon must be between 0 and 360 degrees (0<=lon<360)')
       endif
       
       call find_cols( lat, lon, sathist_nclosest, iphs_owners, ichks, icols, &
                       gcols, iphs_dists, plats, plons )

       if (has_dyn_flds) then
          call dyn_grid_get_colndx( gcols, sathist_nclosest, idyn_owners, idyn1s, idyn2s )
       endif

       do j = 1, sathist_nclosest
          
          if (debug .and. iam==iphs_owners(j) ) then
             if ( abs(plats(j)-rlats(j))>1.e-3_r8 ) then
                write(*,'(a,3f20.12)') ' lat, plat, rlat = ', lat, plats(j), rlats(j)
                write(*,'(a,3f20.12)') ' lon, plon, rlon = ', lon, plons(j), rlons(j)
                call endrun('sat_hist::get_indices: dyn lat is different than phys lat ')
             endif
             if ( abs(plons(j)-rlons(j))>1.e-3_r8 ) then
                write(*,'(a,3f20.12)') ' lat, plat, rlat = ', lat, plats(j), rlats(j)
                write(*,'(a,3f20.12)') ' lon, plon, rlon = ', lon, plons(j), rlons(j)
                call endrun('sat_hist::get_indices: dyn lon is different than phys lon ')
             endif
          endif
          
          ndx = ndx+1
          
          chk_ndxs(ndx)   = ichks(j)
          col_ndxs(ndx)   = icols(j)
          fdyn_ndxs(ndx)  = idyn1s(j)
          ldyn_ndxs(ndx)  = idyn2s(j)
          mlats(ndx)      = plats(j)
          mlons(ndx)      = plons(j)
          phs_owners(ndx) = iphs_owners(j)
          dyn_owners(ndx) = idyn_owners(j)
          phs_dists(ndx)  = iphs_dists(j)
       enddo
    enddo

    deallocate(ichks, icols, idyn1s, idyn2s, iphs_owners, idyn_owners)
    deallocate(rlats, rlons, plats, plons, iphs_dists )

    call t_stopf ('sat_hist::get_indices')
  end subroutine get_indices

!-------------------------------------------------------------------------------
! utility function
!-------------------------------------------------------------------------------
  integer function define_var( var_name, coldim, infile, in_vid, outfile, out_id ) result(res)

    use pio, only: pio_inq_vartype

    character(len=*), intent(in) :: var_name
    integer,          intent(in) :: coldim
    type(File_desc_t),intent(inout) :: infile
    type(File_desc_t),intent(inout) :: outfile
    integer,          intent(in) :: in_vid
    type(var_desc_t), intent(out):: out_id

    integer :: type

    res = pio_inq_varid( outfile, var_name, out_id )
    if(res/=PIO_NOERR) then

       res = pio_inq_vartype( infile, in_vid, type )

       res = pio_def_var ( outfile, var_name, type, (/coldim/), out_id )

       res = copy_att( infile, in_vid, 'long_name', outfile, out_id )
       res = copy_att( infile, in_vid, 'units',     outfile, out_id )

    endif

  end function define_var

!-------------------------------------------------------------------------------
! utility function
!-------------------------------------------------------------------------------
  integer function copy_data( infile, in_vid, outfile, out_id, instart, outstart, ncols ) result(res)

    type(File_desc_t),intent(in) :: infile
    type(File_desc_t),intent(inout) :: outfile
    integer,          intent(in) :: in_vid
    type(var_desc_t), intent(in) :: out_id
    integer,          intent(in) :: instart, outstart, ncols

    real(r8), allocatable :: data(:)
    real(r8), allocatable :: outdata(:)
    integer               :: i

    allocate( data(ncols) )

    res = pio_get_var( infile,  in_vid, (/instart/),  (/ncols/), data )

    allocate( outdata(ncols * sathist_nclosest) )
    
    do i = 1, ncols
      outdata(((i-1)*sathist_nclosest)+1 : (i*sathist_nclosest)) = data(i)
    enddo
  
    res = pio_put_var( outfile, out_id, (/outstart/), (/ncols * sathist_nclosest/), outdata )

    deallocate(outdata)
    deallocate(data)

  end function copy_data

!-------------------------------------------------------------------------------
! utility function
! -- should be able to use pio_copy_att which does not seem to work
!-------------------------------------------------------------------------------
  integer function copy_att( infile, in_vid, att_name, outfile, out_id ) result(res)

    type(File_desc_t),intent(inout) :: infile
    type(File_desc_t),intent(inout) :: outfile
    character(len=*), intent(in) :: att_name
    integer,          intent(in) :: in_vid
    type(var_desc_t), intent(in) :: out_id

    character(len=1024) :: att
    

    res = pio_get_att( infile, in_vid, trim(att_name), att )
    if (res==PIO_NOERR) then
       res = pio_put_att ( outfile, out_id, trim(att_name), trim(att))
    endif


  end function copy_att
  
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------
  subroutine find_cols(lat, lon, nclosest, owner, lcid, icol, gcol, distmin, mlats, mlons)
    use physconst,  only: rearth
    use phys_grid,  only: get_rlon_all_p, get_rlat_all_p, get_gcol_p, get_ncols_p
    use spmd_utils, only: iam, npes, mpi_integer, mpi_real8, mpicom

    real(r8),intent(in)  :: lat, lon            ! requested location in degrees
    integer, intent(in)  :: nclosest            ! number of closest points to find
    integer, intent(out) :: owner(nclosest)     ! rank of chunk owner
    integer, intent(out) :: lcid(nclosest)      ! local chunk index
    integer, intent(out) :: icol(nclosest)      ! column index within the chunk
    integer, intent(out) :: gcol(nclosest)      ! global column index 
    real(r8),intent(out) :: distmin(nclosest)   ! the distance (m) of the closest column(s)
    real(r8),intent(out) :: mlats(nclosest)     ! the latitude of the closest column(s)
    real(r8),intent(out) :: mlons(nclosest)     ! the longitude of the closest column(s)

    real(r8) :: dist
    real(r8) :: rlats(pcols), rlons(pcols)
    real(r8) :: latr, lonr

    integer :: my_owner(nclosest)
    integer :: my_lcid(nclosest)
    integer :: my_icol(nclosest)
    integer :: my_gcol(nclosest)
    real(r8) :: my_distmin(nclosest)
    real(r8) :: my_mlats(nclosest)
    real(r8) :: my_mlons(nclosest)

    integer  :: c, i, j, k, ierr, ncols, mindx(1)
    real(r8) :: sendbufr(3)
    real(r8) :: recvbufr(3,npes)
    integer  :: sendbufi(4)
    integer  :: recvbufi(4,npes)

    call t_startf ('sat_hist::find_cols')

    latr = lat/rad2deg              ! to radians
    lonr = lon/rad2deg              ! to radians
    
    my_owner(:)   = -999
    my_lcid(:)    = -999
    my_icol(:)    = -999
    my_gcol(:)    = -999
    my_mlats(:)   = -999
    my_mlons(:)   = -999
    my_distmin(:) = 1.e10_r8

    chk_loop: do c=begchunk,endchunk
       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, pcols, rlats)
       call get_rlon_all_p(c, pcols, rlons)

       col_loop: do i = 1,ncols
          ! Use the Spherical Law of Cosines to find the great-circle distance.
          dist = acos(sin(latr) * sin(rlats(i)) + cos(latr) * cos(rlats(i)) * cos(rlons(i) - lonr)) * rearth       

          closest_loop: do j = nclosest, 1, -1
             if (dist < my_distmin(j)) then

                if (j < nclosest) then
                   my_distmin(j+1) = my_distmin(j)
                   my_owner(j+1)   = my_owner(j)
                   my_lcid(j+1)    = my_lcid(j)
                   my_icol(j+1)    = my_icol(j)
                   my_gcol(j+1)    = my_gcol(j)
                   my_mlats(j+1)   = my_mlats(j)
                   my_mlons(j+1)   = my_mlons(j)
                end if

                my_distmin(j) = dist
                my_owner(j)   = iam
                my_lcid(j)    = c
                my_icol(j)    = i
                my_gcol(j)    = get_gcol_p(c,i)
                my_mlats(j)   = rlats(i) * rad2deg
                my_mlons(j)   = rlons(i) * rad2deg
             else
                exit 
             end if
          enddo closest_loop

       enddo col_loop
    enddo chk_loop

    k = 1

    do j = 1, nclosest

       sendbufr(1) = my_distmin(k)
       sendbufr(2) = my_mlats(k)
       sendbufr(3) = my_mlons(k)

       call mpi_allgather( sendbufr, 3, mpi_real8, recvbufr, 3, mpi_real8, mpicom, ierr )

       mindx = minloc(recvbufr(1,:))
       distmin(j) = recvbufr(1,mindx(1))
       mlats(j)   = recvbufr(2,mindx(1))
       mlons(j)   = recvbufr(3,mindx(1))

       sendbufi(1) = my_owner(k)
       sendbufi(2) = my_lcid(k)
       sendbufi(3) = my_icol(k)
       sendbufi(4) = my_gcol(k)

       call mpi_allgather( sendbufi, 4, mpi_integer, recvbufi, 4, mpi_integer, mpicom, ierr )

       owner(j)   = recvbufi(1,mindx(1))
       lcid(j)    = recvbufi(2,mindx(1))
       icol(j)    = recvbufi(3,mindx(1))
       gcol(j)    = recvbufi(4,mindx(1))

       if ( iam == owner(j) ) then
          k = k+1
       endif

    enddo

    call t_stopf ('sat_hist::find_cols')

  end subroutine find_cols

end module sat_hist

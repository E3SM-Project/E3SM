module tracer_data
!----------------------------------------------------------------------- 
! module used to read (and interpolate) offline tracer data (sources and
! mixing ratios)
! Created by: Francis Vitt -- 2 May 2006
! Modified by : Jim Edwards -- 10 March 2009
! Modified by : Cheryl Craig and Chih-Chieh (Jack) Chen  -- February 2010
!----------------------------------------------------------------------- 

  use perf_mod,     only : t_startf, t_stopf
  use shr_kind_mod, only : r8 => shr_kind_r8,r4 => shr_kind_r4, shr_kind_cl, SHR_KIND_CS
  use time_manager, only : get_curr_date, get_step_size, get_curr_calday
  use spmd_utils,   only : masterproc
  use ppgrid,       only : pcols, pver, pverp, begchunk, endchunk
  use cam_abortutils,   only : endrun
  use cam_logfile,  only : iulog
  
  use scamMod, only : scm_observed_aero, single_column
  use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_get_index
  use time_manager, only : set_time_float_from_date, set_date_from_time_float
  use pio,          only : file_desc_t, var_desc_t, &
                           pio_seterrorhandling, pio_internal_error, pio_bcast_error, &
                           pio_setdebuglevel, &
                           pio_char, pio_noerr, &
                           pio_inq_dimid, pio_inq_varid, &
                           pio_def_dim, pio_def_var, &
                           pio_put_att, pio_put_var, &
                           pio_get_var, pio_get_att, pio_nowrite, pio_inq_dimlen, &
                           pio_inq_vardimid, pio_inq_dimlen, pio_closefile, &
                           pio_inquire_variable

  implicit none

  private  ! all unless made public
  save 

  public :: trfld, input3d, input2d, trfile
  public :: trcdata_init
  public :: advance_trcdata
  public :: get_fld_data
  public :: put_fld_data
  public :: get_fld_ndx
  public :: write_trc_restart
  public :: read_trc_restart
  public :: init_trc_restart
  public :: incr_filename


  ! !PUBLIC MEMBERS

  type input3d
     real(r8), dimension(:,:,:), pointer :: data => null()
  endtype input3d

  type input2d
     real(r8), dimension(:,:), pointer :: data => null()
  endtype input2d

  type trfld
     real(r8), dimension(:,:,:), pointer :: data => null()
     type(input3d), dimension(4) :: input
     character(len=32) :: srcnam
     character(len=32) :: fldnam
     character(len=32) :: units
     type(var_desc_t) :: var_id
     integer :: coords(4) ! LATDIM | LONDIM | LEVDIM | TIMDIM
     integer :: order(4) ! LATDIM | LONDIM | LEVDIM | TIMDIM
     logical :: srf_fld = .false.
     integer :: pbuf_ndx = -1
  endtype trfld

  type trfile
     type(input2d), dimension(4) :: ps_in
     character(len=shr_kind_cl) :: pathname = ' '
     character(len=shr_kind_cl) :: curr_filename = ' '
     character(len=shr_kind_cl) :: next_filename = ' '
     type(file_desc_t) :: curr_fileid
     type(file_desc_t) :: next_fileid

     type(var_desc_t), pointer :: currfnameid => null() ! pio restart file var id 
     type(var_desc_t), pointer :: nextfnameid => null() ! pio restart file var id

     character(len=shr_kind_cl) :: filenames_list = ''
     real(r8) :: datatimem = -1.e36_r8     ! time of prv. values read in
     real(r8) :: datatimep = -1.e36_r8     ! time of nxt. values read in
     real(r8) :: datatimes(4)
     integer :: interp_recs
     real(r8), pointer, dimension(:) :: curr_data_times => null()
     real(r8), pointer, dimension(:) :: next_data_times => null()
     logical :: remove_trc_file = .false.  ! delete file when finished with it
     real(r8) :: offset_time
     integer :: cyc_ndx_beg
     integer :: cyc_ndx_end
     integer :: cyc_yr = 0
     real(r8) :: one_yr = 0
     real(r8) :: curr_mod_time ! model time - calendar day
     real(r8) :: next_mod_time ! model time - calendar day - next time step
     integer :: nlon
     integer :: nlat
     integer :: nlev
     integer :: nilev
     integer :: ps_coords(3) ! LATDIM | LONDIM | TIMDIM
     integer :: ps_order(3) ! LATDIM | LONDIM | TIMDIM
     real(r8), pointer, dimension(:) :: lons => null()
     real(r8), pointer, dimension(:) :: lats => null()
     real(r8), pointer, dimension(:) :: levs => null()
     real(r8), pointer, dimension(:) :: ilevs => null()
     real(r8), pointer, dimension(:) :: hyam => null()
     real(r8), pointer, dimension(:) :: hybm => null()
     real(r8), pointer, dimension(:,:) :: ps => null()
     real(r8), pointer, dimension(:) :: hyai => null()
     real(r8), pointer, dimension(:) :: hybi => null()
     real(r8), pointer, dimension(:,:) :: weight_x => null(), weight_y => null()
     integer, pointer, dimension(:) :: count_x => null(), count_y => null()
     integer, pointer, dimension(:,:) :: index_x => null(), index_y => null()
     real(r8)                        :: p0
     type(var_desc_t) :: ps_id
     logical,  allocatable, dimension(:) :: in_pbuf
     logical :: has_ps = .false.
     logical :: zonal_ave = .false.
     logical :: alt_data = .false.     
     logical :: cyclical = .false.
     logical :: cyclical_list = .false.
     logical :: weight_by_lat = .false.
     logical :: conserve_column = .false.
     logical :: fill_in_months = .false.
     logical :: fixed = .false.
     logical :: initialized = .false.
     logical :: top_bndry = .false.
     logical :: stepTime = .false.  ! Do not interpolate in time, but use stepwise times
  endtype trfile

  integer, public, parameter :: MAXTRCRS = 100

  integer, parameter :: LONDIM = 1
  integer, parameter :: LATDIM = 2
  integer, parameter :: LEVDIM = 3
  integer, parameter :: TIMDIM = 4

  integer, parameter :: PS_TIMDIM = 3

  integer, parameter :: ZA_LATDIM = 1
  integer, parameter :: ZA_LEVDIM = 2
  integer, parameter :: ZA_TIMDIM = 3

  integer, parameter :: nm=1    ! array index for previous (minus) data
  integer, parameter :: np=2    ! array index for next (plus) data

  integer :: plon, plat

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
  subroutine trcdata_init( specifier, filename, filelist, datapath, flds, file, &
                           rmv_file, data_cycle_yr, data_fixed_ymd, data_fixed_tod, data_type )

    use mo_constants,    only : d2r
    use cam_control_mod, only : nsrest
    use dyn_grid,        only : get_dyn_grid_parm
    use string_utils,    only : to_upper
    use horizontal_interpolate, only : xy_interp_init
#if ( defined SPMD )
    use mpishorthand,    only: mpicom, mpir8, mpiint
#endif

    implicit none

    character(len=*),    intent(in)    :: specifier(:)
    character(len=*),    intent(in)    :: filename
    character(len=*),    intent(in)    :: filelist
    character(len=*),    intent(in)    :: datapath
    type(trfld), dimension(:), pointer :: flds
    type(trfile),        intent(inout) :: file
    logical,             intent(in)    :: rmv_file
    integer,             intent(in)    :: data_cycle_yr
    integer,             intent(in)    :: data_fixed_ymd
    integer,             intent(in)    :: data_fixed_tod
    character(len=*),    intent(in)    :: data_type

    integer :: f, mxnflds, astat
    integer :: str_yr, str_mon, str_day
    integer :: lon_dimid, lat_dimid, lev_dimid, tim_dimid, old_dimid
    integer :: dimids(4), did
    type(var_desc_t) :: varid
    integer :: idx
    integer :: ierr
    integer :: errcode
    real(r8) :: start_time, time1, time2
    integer :: i1,i2,j1,j2
    integer :: nvardims, vardimids(4)

    character(len=80) :: data_units

    call specify_fields( specifier, flds )

    file%datatimep=-1.e36_r8
    file%datatimem=-1.e36_r8

    mxnflds = 0 
    if (associated(flds)) mxnflds = size( flds )

    if (mxnflds < 1) return
    
    file%remove_trc_file = rmv_file
    file%pathname = trim(datapath)
    file%filenames_list = trim(filelist)

    file%fill_in_months = .false.
    file%cyclical = .false.  
    file%cyclical_list = .false.  

! does not work when compiled with pathf90
!    select case ( to_upper(data_type) )
    select case ( data_type )
    case( 'FIXED' )
       file%fixed = .true.
    case( 'INTERP_MISSING_MONTHS' )
       file%fill_in_months = .true.
    case( 'CYCLICAL' )
       file%cyclical = .true.
       file%cyc_yr = data_cycle_yr
    case( 'CYCLICAL_LIST' )
       file%cyclical_list = .true.
       file%cyc_yr = data_cycle_yr
    case( 'SERIAL' )
    case default 
       write(iulog,*) 'trcdata_init: invalid data type: '//trim(data_type)//' file: '//trim(filename)
       write(iulog,*) 'trcdata_init: valid data types: SERIAL | CYCLICAL | CYCLICAL_LIST | FIXED | INTERP_MISSING_MONTHS '
       call endrun('trcdata_init: invalid data type: '//trim(data_type)//' file: '//trim(filename))
    endselect

    if ( (.not.file%fixed) .and. ((data_fixed_ymd>0._r8) .or.(data_fixed_tod>0._r8))) then
       call endrun('trcdata_init: Cannot specify data_fixed_ymd or data_fixed_tod if data type is not FIXED')
    endif
    if ( (.not.file%cyclical) .and. (data_cycle_yr>0._r8) ) then
       call endrun('trcdata_init: Cannot specify data_cycle_yr if data type is not CYCLICAL')
    endif

    if (masterproc) then
       write(iulog,*) 'trcdata_init: data type: '//trim(data_type)//' file: '//trim(filename)
    endif

    ! if there is no list of files (len_trim(file%filenames_list)<1) then
    !  -> set curr_filename from namelist rather from restart data
    if ( len_trim(file%curr_filename)<1 .or. len_trim(file%filenames_list)<1 .or. file%fixed ) then ! initial run
       file%curr_filename = trim(filename)

       call get_model_time(file)

       if ( file%fixed ) then
          str_yr = data_fixed_ymd/10000
          str_mon = (data_fixed_ymd - str_yr*10000)/100
          str_day = data_fixed_ymd - str_yr*10000 - str_mon*100
          call set_time_float_from_date( start_time, str_yr, str_mon, str_day, data_fixed_tod )
          file%offset_time = start_time - file%curr_mod_time
       else
          file%offset_time = 0
       endif
    endif

    call set_time_float_from_date( time2, 2, 1, 1, 0 )
    call set_time_float_from_date( time1, 1, 1, 1, 0 )
    file%one_yr = time2-time1

    if ( file%cyclical .or. file%cyclical_list) then
       file%cyc_ndx_beg = -1
       file%cyc_ndx_end = -1
       if ( file%cyc_yr /= 0 ) then
          call set_time_float_from_date( time1, file%cyc_yr  , 1, 1, 0 )
          call set_time_float_from_date( time2, file%cyc_yr+1, 1, 1, 0 )
          file%one_yr = time2-time1
       endif

       call open_trc_datafile( file%curr_filename, file%pathname, file%curr_fileid, file%curr_data_times, &
            cyc_ndx_beg=file%cyc_ndx_beg, cyc_ndx_end=file%cyc_ndx_end, cyc_yr=file%cyc_yr )
    else
       call open_trc_datafile( file%curr_filename, file%pathname, file%curr_fileid, file%curr_data_times )
       file%curr_data_times = file%curr_data_times - file%offset_time
    endif

    call pio_seterrorhandling(File%curr_fileid, PIO_BCAST_ERROR)
    ierr = pio_inq_dimid( file%curr_fileid, 'lon', idx )
    call pio_seterrorhandling(File%curr_fileid, PIO_INTERNAL_ERROR)

    file%zonal_ave = (ierr/=PIO_NOERR)

    plon = get_dyn_grid_parm('plon')
    plat = get_dyn_grid_parm('plat')

    if ( .not. file%zonal_ave ) then

       call get_dimension( file%curr_fileid, 'lon', file%nlon, dimid=old_dimid, data=file%lons )

       file%lons =  file%lons * d2r

       lon_dimid = old_dimid

    endif

    ierr = pio_inq_dimid( file%curr_fileid, 'time', old_dimid)

    ! Hack to work with weird netCDF and old gcc or NAG bug.
    tim_dimid = old_dimid

    call get_dimension( file%curr_fileid, 'lat', file%nlat, dimid=old_dimid, data=file%lats )
    file%lats =  file%lats * d2r

    lat_dimid = old_dimid

    allocate( file%ps(file%nlon,file%nlat), stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'trcdata_init: file%ps allocation error = ',astat
       call endrun('trcdata_init: failed to allocate x array')
    end if

    call pio_seterrorhandling(File%curr_fileid, PIO_BCAST_ERROR)
    ierr = pio_inq_varid( file%curr_fileid, 'PS', file%ps_id )
    file%has_ps = (ierr==PIO_NOERR)
    ierr = pio_inq_dimid( file%curr_fileid, 'altitude', idx )
    file%alt_data = (ierr==PIO_NOERR)

    call pio_seterrorhandling(File%curr_fileid, PIO_INTERNAL_ERROR)

    if ( file%has_ps) then
       ierr = pio_inq_vardimid (file%curr_fileid, file%ps_id, dimids(1:3))
       do did = 1,3
          if      ( dimids(did) == lon_dimid ) then
             file%ps_coords(LONDIM) = did
             file%ps_order(did) = LONDIM
          else if ( dimids(did) == lat_dimid ) then
             file%ps_coords(LATDIM) = did
             file%ps_order(did) = LATDIM
          else if ( dimids(did) == tim_dimid ) then
             file%ps_coords(PS_TIMDIM) = did
             file%ps_order(did) = PS_TIMDIM
          endif
       enddo
    endif

    if (masterproc) then 
       write(iulog,*) 'trcdata_init: file%has_ps = ' , file%has_ps 
    endif ! masterproc

    if (file%alt_data) then
       call get_dimension( file%curr_fileid, 'altitude_int', file%nilev,  data=file%ilevs  )
       call get_dimension( file%curr_fileid, 'altitude',     file%nlev, dimid=old_dimid, data=file%levs  )
    else
       call get_dimension( file%curr_fileid, 'lev', file%nlev, dimid=old_dimid, data=file%levs  )
       if (old_dimid>0) then
          file%levs =  file%levs*100._r8 ! mbar->pascals
       endif
    endif

    ! For some bizarre reason, netCDF with older gcc is keeping a pointer to the dimid, and overwriting it later!
    ! Hackish workaround is to make a copy...
    lev_dimid = old_dimid

    if (file%has_ps) then

       allocate( file%hyam(file%nlev),  file%hybm(file%nlev), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'trcdata_init: file%hyam,file%hybm allocation error = ',astat
          call endrun('trcdata_init: failed to allocate file%hyam and file%hybm arrays')
       end if

       allocate( file%hyai(file%nlev+1),  file%hybi(file%nlev+1), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'trcdata_init: file%hyai,file%hybi allocation error = ',astat
          call endrun('trcdata_init: failed to allocate file%hyai and file%hybi arrays')
       end if

       call pio_seterrorhandling(File%curr_fileid, PIO_BCAST_ERROR)
       ierr = pio_inq_varid( file%curr_fileid, 'P0', varid)
       call pio_seterrorhandling(File%curr_fileid, PIO_INTERNAL_ERROR)

       if ( ierr == PIO_NOERR ) then
          ierr = pio_get_var( file%curr_fileid, varid, file%p0 )
       else
          file%p0 = 100000._r8
       endif
       ierr = pio_inq_varid( file%curr_fileid, 'hyam', varid )
       ierr = pio_get_var( file%curr_fileid, varid, file%hyam )
       ierr = pio_inq_varid( file%curr_fileid, 'hybm', varid )
       ierr = pio_get_var( file%curr_fileid, varid, file%hybm )
       if (file%conserve_column) then
          ierr = pio_inq_varid( file%curr_fileid, 'hyai', varid )
          ierr = pio_get_var( file%curr_fileid, varid, file%hyai )
          ierr = pio_inq_varid( file%curr_fileid, 'hybi', varid )
          ierr = pio_get_var( file%curr_fileid, varid, file%hybi )
       endif

       allocate( file         %ps  (pcols,begchunk:endchunk), stat=astat   )
       if( astat/= 0 ) then
          write(iulog,*) 'trcdata_init: failed to allocate file%ps array; error = ',astat
          call endrun
       end if
       allocate( file%ps_in(1)%data(pcols,begchunk:endchunk), stat=astat   )
       if( astat/= 0 ) then
          write(iulog,*) 'trcdata_init: failed to allocate file%ps_in(1)%data array; error = ',astat
          call endrun
       end if
       allocate( file%ps_in(2)%data(pcols,begchunk:endchunk), stat=astat   )
       if( astat/= 0 ) then
          write(iulog,*) 'trcdata_init: failed to allocate file%ps_in(2)%data array; error = ',astat
          call endrun
       end if
       if( file%fill_in_months ) then
          allocate( file%ps_in(3)%data(pcols,begchunk:endchunk), stat=astat   )
          if( astat/= 0 ) then
             write(iulog,*) 'trcdata_init: failed to allocate file%ps_in(3)%data array; error = ',astat
             call endrun
          end if
          allocate( file%ps_in(4)%data(pcols,begchunk:endchunk), stat=astat   )
          if( astat/= 0 ) then
             write(iulog,*) 'trcdata_init: failed to allocate file%ps_in(4)%data array; error = ',astat
             call endrun
          end if
       end if
    endif

    flds_loop: do f = 1,mxnflds

       ! get netcdf variable id for the field
       ierr = pio_inq_varid( file%curr_fileid, flds(f)%srcnam, flds(f)%var_id )

       ! determine if the field has a vertical dimension

       if (lev_dimid>0) then
          ierr = pio_inquire_variable(  file%curr_fileid, flds(f)%var_id, ndims=nvardims )
          ierr = pio_inquire_variable(  file%curr_fileid, flds(f)%var_id, dimids=vardimids(:nvardims) )
          flds(f)%srf_fld = .not.any(vardimids(:nvardims)==lev_dimid)
       else
          flds(f)%srf_fld = .true.
       endif

       ! allocate memory only if not already in pbuf2d

       if ( .not. file%in_pbuf(f) ) then 
          if ( flds(f)%srf_fld .or. file%top_bndry ) then
             allocate( flds(f)         %data(pcols,1,begchunk:endchunk), stat=astat   )
          else
             allocate( flds(f)         %data(pcols,pver,begchunk:endchunk), stat=astat   )
          endif
          if( astat/= 0 ) then
             write(iulog,*) 'trcdata_init: failed to allocate flds(f)%data array; error = ',astat
             call endrun
          end if
       else
          flds(f)%pbuf_ndx = pbuf_get_index(flds(f)%fldnam,errcode)
       endif
   
       if (flds(f)%srf_fld) then
          allocate( flds(f)%input(1)%data(pcols,1,begchunk:endchunk), stat=astat   )
       else
          allocate( flds(f)%input(1)%data(pcols,file%nlev,begchunk:endchunk), stat=astat   )
       endif
       if( astat/= 0 ) then
          write(iulog,*) 'trcdata_init: failed to allocate flds(f)%input(1)%data array; error = ',astat
          call endrun
       end if
       if (flds(f)%srf_fld) then
          allocate( flds(f)%input(2)%data(pcols,1,begchunk:endchunk), stat=astat   )
       else
          allocate( flds(f)%input(2)%data(pcols,file%nlev,begchunk:endchunk), stat=astat   )
       endif
       if( astat/= 0 ) then
          write(iulog,*) 'trcdata_init: failed to allocate flds(f)%input(2)%data array; error = ',astat
          call endrun
       end if

       if( file%fill_in_months ) then
          if (flds(f)%srf_fld) then
             allocate( flds(f)%input(3)%data(pcols,1,begchunk:endchunk), stat=astat   )
          else
             allocate( flds(f)%input(3)%data(pcols,file%nlev,begchunk:endchunk), stat=astat   )
          endif
          if( astat/= 0 ) then
             write(iulog,*) 'trcdata_init: failed to allocate flds(f)%input(3)%data array; error = ',astat
             call endrun
          end if
          if (flds(f)%srf_fld) then
             allocate( flds(f)%input(4)%data(pcols,1,begchunk:endchunk), stat=astat   )
          else
             allocate( flds(f)%input(4)%data(pcols,file%nlev,begchunk:endchunk), stat=astat   )
          endif
          if( astat/= 0 ) then
             write(iulog,*) 'trcdata_init: failed to allocate flds(f)%input(4)%data array; error = ',astat
             call endrun
          end if
       endif

       if ( file%zonal_ave ) then
          ierr = pio_inq_vardimid (file%curr_fileid, flds(f)%var_id, dimids(1:3))
          do did = 1,3
             if      ( dimids(did) == lat_dimid ) then
                flds(f)%coords(ZA_LATDIM) = did
                flds(f)%order(did) = ZA_LATDIM
             else if ( dimids(did) == lev_dimid ) then
                flds(f)%coords(ZA_LEVDIM) = did
                flds(f)%order(did) = ZA_LEVDIM
             else if ( dimids(did) == tim_dimid ) then
                flds(f)%coords(ZA_TIMDIM) = did
                flds(f)%order(did) = ZA_TIMDIM
             endif
          enddo
       else if ( flds(f)%srf_fld ) then
          ierr = pio_inq_vardimid (file%curr_fileid, flds(f)%var_id, dimids(1:3))
          do did = 1,3
             if      ( dimids(did) == lon_dimid ) then
                flds(f)%coords(LONDIM) = did
                flds(f)%order(did) = LONDIM
             else if ( dimids(did) == lat_dimid ) then
                flds(f)%coords(LATDIM) = did
                flds(f)%order(did) = LATDIM
             else if ( dimids(did) == tim_dimid ) then
                flds(f)%coords(PS_TIMDIM) = did
                flds(f)%order(did) = PS_TIMDIM
             endif
          enddo
       else
          ierr = pio_inq_vardimid (file%curr_fileid, flds(f)%var_id, dimids)
          do did = 1,4
             if      ( dimids(did) == lon_dimid ) then
                flds(f)%coords(LONDIM) = did
                flds(f)%order(did) = LONDIM
             else if ( dimids(did) == lat_dimid ) then
                flds(f)%coords(LATDIM) = did
                flds(f)%order(did) = LATDIM
             else if ( dimids(did) == lev_dimid ) then
                flds(f)%coords(LEVDIM) = did
                flds(f)%order(did) = LEVDIM
             else if ( dimids(did) == tim_dimid ) then
                flds(f)%coords(TIMDIM) = did
                flds(f)%order(did) = TIMDIM
             endif
          enddo
       endif

       ierr = pio_get_att( file%curr_fileid, flds(f)%var_id, 'units', data_units)
       data_units = trim(data_units)
       flds(f)%units = data_units(1:32)

    enddo flds_loop

! if weighting by latitude, compute weighting for horizontal interpolation
    if( file%weight_by_lat ) then
! get dimensions of CAM resolution
        plon = get_dyn_grid_parm('plon')
        plat = get_dyn_grid_parm('plat')
        
! weight_x & weight_y are weighting function for x & y interpolation
        allocate(file%weight_x(plon,file%nlon))
        allocate(file%weight_y(plat,file%nlat))
        allocate(file%count_x(plon))
        allocate(file%count_y(plat))
        allocate(file%index_x(plon,file%nlon))
        allocate(file%index_y(plat,file%nlat))
        file%weight_x(:,:) = 0.0_r8
        file%weight_y(:,:) = 0.0_r8
        file%count_x(:) = 0
        file%count_y(:) = 0
        file%index_x(:,:) = 0
        file%index_y(:,:) = 0

        if(masterproc) then
! compute weighting 
            call xy_interp_init(file%nlon,file%nlat,file%lons,file%lats,plon,plat,file%weight_x,file%weight_y)

            do i2=1,plon
               file%count_x(i2) = 0
               do i1=1,file%nlon
                  if(file%weight_x(i2,i1).gt.0.0_r8 ) then
                     file%count_x(i2) = file%count_x(i2) + 1
                     file%index_x(i2,file%count_x(i2)) = i1
                  endif
               enddo
            enddo

            do j2=1,plat
               file%count_y(j2) = 0
               do j1=1,file%nlat
                  if(file%weight_y(j2,j1).gt.0.0_r8 ) then
                     file%count_y(j2) = file%count_y(j2) + 1
                     file%index_y(j2,file%count_y(j2)) = j1
                  endif
               enddo
            enddo
        endif

#if ( defined SPMD)
        call mpibcast(file%weight_x, plon*file%nlon, mpir8 , 0, mpicom)
        call mpibcast(file%weight_y, plat*file%nlat, mpir8 , 0, mpicom)
        call mpibcast(file%count_x, plon, mpiint , 0, mpicom)
        call mpibcast(file%count_y, plat, mpiint , 0, mpicom)
        call mpibcast(file%index_x, plon*file%nlon, mpiint , 0, mpicom)
        call mpibcast(file%index_y, plat*file%nlat, mpiint , 0, mpicom)
#endif
    endif

  end subroutine trcdata_init

!-----------------------------------------------------------------------
! Reads more data if needed and interpolates data to current model time 
!-----------------------------------------------------------------------
  subroutine advance_trcdata( flds, file, state, pbuf2d )
    use physics_types,only : physics_state
    use physics_buffer, only : physics_buffer_desc
    use ppgrid, only : pver,pcols

    implicit none

    type(trfile),        intent(inout) :: file
    type(trfld),         intent(inout) :: flds(:)
    type(physics_state), intent(in)    :: state(begchunk:endchunk)
    
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    integer :: ncol
    real(r8) :: data_time
    real(r8) :: t(pcols,pver)          ! input temperature (K)
    real(r8) :: rho(pcols,pver)          ! input temperature (K)
    real(r8) :: pmid(pcols,pver)       ! pressure at layer midpoints (pa)
!--------------------------------------BLH-----------------------------------    

    call t_startf('advance_trcdata')
    if ( .not.( file%fixed .and. file%initialized ) ) then

       call get_model_time(file)

       data_time = file%datatimep

       if ( file%cyclical .or. file%cyclical_list ) then
          ! wrap around
          if ( (file%datatimep<file%datatimem) .and. (file%curr_mod_time>file%datatimem) ) then
             data_time = data_time + file%one_yr 
          endif
       endif

    ! For stepTime need to advance if the times are equal
    ! Should not impact other runs?
       if ( file%curr_mod_time >= data_time ) then
          call t_startf('read_next_trcdata')
          call read_next_trcdata(state, flds, file )
          call t_stopf('read_next_trcdata')
          if(masterproc) write(iulog,*) 'READ_NEXT_TRCDATA ', flds%fldnam
       end if

    endif
    
    ! need to interpolate the data, regardless
    ! each mpi task needs to interpolate
    call t_startf('interpolate_trcdata')
    call interpolate_trcdata( state, flds, file, pbuf2d )
    call t_stopf('interpolate_trcdata')

    file%initialized = .true.

    call t_stopf('advance_trcdata')

  end subroutine advance_trcdata

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine get_fld_data( flds, field_name, data, ncol, lchnk, pbuf )

    use physics_buffer, only : physics_buffer_desc, pbuf_get_field

    implicit none

    type(trfld), intent(inout) :: flds(:)
    character(len=*), intent(in) :: field_name
    real(r8), intent(out) :: data(:,:)
    integer, intent(in) :: lchnk
    integer, intent(in) :: ncol
    type(physics_buffer_desc), pointer :: pbuf(:)
    

    integer :: f, nflds
    real(r8),pointer  :: tmpptr(:,:)

    data(:,:) = 0._r8
    nflds = size(flds)

    do f = 1, nflds
       if ( trim(flds(f)%fldnam) == trim(field_name) ) then
          if ( flds(f)%pbuf_ndx>0 ) then
             call pbuf_get_field(pbuf, flds(f)%pbuf_ndx, tmpptr)
             data(:ncol,:) = tmpptr(:ncol,:)
          else
             data(:ncol,:) = flds(f)%data(:ncol,:,lchnk)
          endif
       endif
    enddo

 end subroutine get_fld_data

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine put_fld_data( flds, field_name, data, ncol, lchnk, pbuf )

    use physics_buffer, only : physics_buffer_desc, pbuf_get_field

    implicit none

    type(trfld), intent(inout) :: flds(:)
    character(len=*), intent(in) :: field_name
    real(r8), intent(in) :: data(:,:)
    integer, intent(in) :: lchnk
    integer, intent(in) :: ncol
    type(physics_buffer_desc), pointer :: pbuf(:)
    

    integer :: f, nflds
    real(r8),pointer  :: tmpptr(:,:)

    nflds = size(flds)

    do f = 1, nflds
       if ( trim(flds(f)%fldnam) == trim(field_name) ) then
          if ( flds(f)%pbuf_ndx>0 ) then
             call pbuf_get_field(pbuf, flds(f)%pbuf_ndx, tmpptr)
             tmpptr(:ncol,:) = data(:ncol,:)
          else
             flds(f)%data(:ncol,:,lchnk) = data(:ncol,:)
          endif
       endif
    enddo

 end subroutine put_fld_data

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine get_fld_ndx( flds, field_name, idx  )

    implicit none

    type(trfld), intent(in) :: flds(:)
    character(len=*), intent(in) :: field_name
    integer, intent(out) :: idx    
    integer :: f, nflds

    idx = -1
    nflds = size(flds)

    do f = 1, nflds
       if ( trim(flds(f)%fldnam) == trim(field_name) ) then
          idx = f
          return
       endif
    enddo

  end subroutine get_fld_ndx

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  subroutine get_model_time(file)
    implicit none
    type(trfile), intent(inout) :: file

    integer yr, mon, day, ncsec  ! components of a date

    call get_curr_date(yr, mon, day, ncsec)

    if ( file%cyclical .or. file%cyclical_list) yr = file%cyc_yr
    call set_time_float_from_date( file%curr_mod_time, yr, mon, day, ncsec )
    file%next_mod_time = file%curr_mod_time + get_step_size()/86400._r8

  end subroutine get_model_time

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  subroutine check_files( file, fids, itms, times_found)

    implicit none

    type(trfile),      intent(inout) :: file
    type(file_desc_t), intent(out)   :: fids(2) ! ids of files that contains these recs
    integer, optional, intent(out)   :: itms(2)
    logical, optional, intent(inout) :: times_found

    !-----------------------------------------------------------------------
    ! 	... local variables
    !-----------------------------------------------------------------------
    logical            :: list_cycled

    list_cycled = .false.

   !-----------------------------------------------------------------------
   !        If next time beyond the end of the time list,
   !        then increment the filename and move on to the next file
   !-----------------------------------------------------------------------
    if ((file%next_mod_time > file%curr_data_times(size(file%curr_data_times))).or.file%cyclical_list) then
       if (file%cyclical_list) then
          if ( associated(file%next_data_times) ) then
              if ((file%curr_mod_time > file%datatimep)) then

               call advance_file(file)
     
            endif
         endif

       endif
       if ( .not. associated(file%next_data_times) ) then
          ! open next file if not already opened...
          if (file%cyclical_list) then
              file%next_filename = incr_filename( file%curr_filename, filenames_list=file%filenames_list, datapath=file%pathname ,&
                                                  cyclical_list=file%cyclical_list, list_cycled=list_cycled)
          else
              file%next_filename = incr_filename( file%curr_filename, filenames_list=file%filenames_list, datapath=file%pathname)
          endif
          call open_trc_datafile( file%next_filename, file%pathname, file%next_fileid, file%next_data_times )
          file%next_data_times = file%next_data_times - file%offset_time
       endif
    endif
    
    !-----------------------------------------------------------------------
    !        If using next_data_times and the current is greater than or equal to the next, then
    !        close the current file, and set up for next file.
    !-----------------------------------------------------------------------
    if ( associated(file%next_data_times) ) then
       if (file%cyclical_list .and. list_cycled) then    ! special case - list cycled 

          file%datatimem = file%curr_data_times(size(file%curr_data_times)) 
          itms(1)=size(file%curr_data_times)
          fids(1)=file%curr_fileid

          file%datatimep = file%next_data_times(1)
          itms(2)=1
          fids(2) = file%next_fileid

          times_found = .true.

       else if (file%curr_mod_time >= file%next_data_times(1)) then

          call advance_file(file)

       endif
    endif

  end subroutine check_files

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  function incr_filename( filename, filenames_list, datapath, cyclical_list, list_cycled )

    !-----------------------------------------------------------------------
    ! 	... Increment or decrement a date string withing a filename
    !           the filename date section is assumed to be of the form
    !           yyyy-dd-mm
    !-----------------------------------------------------------------------

    use string_utils, only : incstr
    use shr_file_mod, only : shr_file_getunit, shr_file_freeunit

    implicit none


    character(len=*),           intent(in)    :: filename ! present dynamical dataset filename
    character(len=*), optional, intent(in)    :: filenames_list 
    character(len=*), optional, intent(in)    :: datapath
    logical         , optional, intent(in)    :: cyclical_list  ! If true, allow list to cycle
    logical         , optional, intent(out)   :: list_cycled
    character(len=shr_kind_cl)                :: incr_filename         ! next filename in the sequence


    ! set new next_filename ...

    !-----------------------------------------------------------------------
    !	... local variables
    !-----------------------------------------------------------------------
    integer :: pos, pos1, istat
    character(len=shr_kind_cl) :: fn_new, line, filepath
    character(len=6)   :: seconds
    character(len=5)   :: num
    integer :: ios,unitnumber

    if (present(list_cycled)) list_cycled = .false.

    if (( .not. present(filenames_list)) .or.(len_trim(filenames_list) == 0)) then
       !-----------------------------------------------------------------------
       !	... ccm type filename
       !-----------------------------------------------------------------------
       pos = len_trim( filename )
       fn_new = filename(:pos)
       if ( masterproc ) write(iulog,*) 'incr_flnm: old filename = ',trim(fn_new)
       if( fn_new(pos-2:) == '.nc' ) then
          pos = pos - 3
       end if
       istat = incstr( fn_new(:pos), 1 )
       if( istat /= 0 ) then
          write(iulog,*) 'incr_flnm: incstr returned ', istat
          write(iulog,*) '           while trying to decrement ',trim( fn_new )
          call endrun
       end if

    else

       !-------------------------------------------------------------------
       !  ... open filenames_list
       !-------------------------------------------------------------------
       if ( masterproc ) write(iulog,*) 'incr_flnm: old filename = ',trim(filename)
       if ( masterproc ) write(iulog,*) 'incr_flnm: open filenames_list : ',trim(filenames_list)
       unitnumber = shr_file_getUnit()
       if ( present(datapath) ) then
         filepath = trim(datapath) //'/'// trim(filenames_list)
       else
         filepath = trim(datapath)
       endif

       open( unit=unitnumber, file=filepath, iostat=ios, status="OLD")
       if (ios /= 0) then
          call endrun('not able to open filenames_list file: '//trim(filepath))
       endif

       !-------------------------------------------------------------------
       !  ...  read file names
       !-------------------------------------------------------------------
       read( unit=unitnumber, fmt='(A)', iostat=ios ) line 
       if (ios /= 0) then
          call endrun('not able to increment file name from filenames_list file: '//trim(filenames_list))
       endif

       !-------------------------------------------------------------------
       !      If current filename is '', then initialize with the first filename read in
       !      and skip this section.
       !-------------------------------------------------------------------
       if (filename /= '') then 

          !-------------------------------------------------------------------
          !       otherwise read until find current filename
          !-------------------------------------------------------------------
          do while( trim(line) /= trim(filename) )
             read( unit=unitnumber, fmt='(A)', iostat=ios ) line 
             if (ios /= 0) then
                call endrun('not able to increment file name from filenames_list file: '//trim(filenames_list))
             endif
          enddo
   
          !-------------------------------------------------------------------
          !      Read next filename
          !-------------------------------------------------------------------
          read( unit=unitnumber, fmt='(A)', iostat=ios ) line 

          !---------------------------------------------------------------------------------
          !       If cyclical_list, then an end of file is not an error, but rather 
          !       a signal to rewind and start over
          !---------------------------------------------------------------------------------

          if (ios /= 0) then
             if (present(cyclical_list)) then
                if (cyclical_list) then
                   list_cycled=.true.
                   rewind(unitnumber)
                   read( unit=unitnumber, fmt='(A)', iostat=ios ) line 
                     ! Error here should never happen, but check just in case
                   if (ios /= 0) then
                      call endrun('not able to increment file name from filenames_list file: '//trim(filenames_list))
                   endif
                else
                   call endrun('not able to increment file name from filenames_list file: '//trim(filenames_list))
                endif
             else
                call endrun('not able to increment file name from filenames_list file: '//trim(filenames_list))
             endif
          endif

       endif

       !---------------------------------------------------------------------------------
       !     Assign the current filename and close the filelist
       !---------------------------------------------------------------------------------
       fn_new = trim(line)

       close(unit=unitnumber)
       call shr_file_freeUnit(unitnumber)
    endif

    !---------------------------------------------------------------------------------
    !      return the current filename 
    !---------------------------------------------------------------------------------
    incr_filename = trim(fn_new)
    if ( masterproc ) write(iulog,*) 'incr_flnm: new filename = ',trim(incr_filename)

  end function incr_filename

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  subroutine find_times( itms, fids, time, file, datatimem, datatimep, times_found )

    implicit none

    type(trfile), intent(in) :: file
    real(r8), intent(out) :: datatimem, datatimep

    integer, intent(out) :: itms(2) ! record numbers that bracket time
    type(file_desc_t), intent(out) :: fids(2) ! ids of files that contains these recs

    real(r8), intent(in) :: time    ! time of interest
    logical, intent(inout)  :: times_found

    integer :: np1        ! current forward time index of dataset
    integer :: n,i      ! 
    integer :: curr_tsize, next_tsize, all_tsize
    integer :: astat
    integer :: cyc_tsize

    real(r8), allocatable, dimension(:):: all_data_times

    curr_tsize = size(file%curr_data_times)
    next_tsize = 0
    if ( associated(file%next_data_times)) next_tsize = size(file%next_data_times)

    all_tsize = curr_tsize + next_tsize

    allocate( all_data_times( all_tsize ), stat=astat )
    if( astat/= 0 ) then
       write(iulog,*) 'find_times: failed to allocate all_data_times array; error = ',astat
       call endrun
    end if

    all_data_times(:curr_tsize) = file%curr_data_times(:)
    if (next_tsize > 0) all_data_times(curr_tsize+1:all_tsize) = file%next_data_times(:)

    if ( .not. file%cyclical ) then
       if ( all( all_data_times(:) > time ) ) then
          write(iulog,*) 'FIND_TIMES: ALL data times are after ', time
          write(iulog,*) 'FIND_TIMES: data times: ',all_data_times(:)
          write(iulog,*) 'FIND_TIMES: time: ',time
          call endrun('find_times: all(all_data_times(:) > time) '// trim(file%curr_filename) )
       endif

       ! find bracketing times 
       find_times_loop : do n=1, all_tsize-1
          np1 = n + 1
          datatimem = all_data_times(n)   !+ file%offset_time
          datatimep = all_data_times(np1) !+ file%offset_time
       ! When stepTime, datatimep may not equal the time (as only datatimem is used)
       ! Should not break other runs?
          if ( (time .ge. datatimem) .and. (time .lt. datatimep) ) then
             times_found = .true.
             exit find_times_loop
          endif
       enddo find_times_loop

    else  ! file%cyclical

       cyc_tsize = file%cyc_ndx_end - file%cyc_ndx_beg + 1

       if ( cyc_tsize > 1 ) then

          call findplb(all_data_times(file%cyc_ndx_beg:file%cyc_ndx_end),cyc_tsize, time, n )

          if (n == cyc_tsize) then
             np1 = 1
          else
             np1 = n+1
          endif

          datatimem = all_data_times(n  +file%cyc_ndx_beg-1)   
          datatimep = all_data_times(np1+file%cyc_ndx_beg-1) 
          times_found = .true.

       endif
    endif

    if ( .not. times_found ) then
       if (masterproc) then
          write(iulog,*)'FIND_TIMES: Failed to find dates bracketing desired time =', time
          write(iulog,*)' datatimem = ',file%datatimem
          write(iulog,*)' datatimep = ',file%datatimep
          write(iulog,*)' all_data_times = ',all_data_times
          !call endrun()
          return
       endif
    endif

    deallocate( all_data_times, stat=astat )
    if( astat/= 0 ) then
       write(iulog,*) 'find_times: failed to deallocate all_data_times array; error = ',astat
       call endrun
    end if
  
    if ( .not. file%cyclical ) then
      itms(1) = n
      itms(2) = np1
    else
      itms(1) = n   +file%cyc_ndx_beg-1
      itms(2) = np1 +file%cyc_ndx_beg-1
    endif

    fids(:) = file%curr_fileid

    do i=1,2
       if ( itms(i) > curr_tsize ) then 
          itms(i) = itms(i) - curr_tsize 
          fids(i) = file%next_fileid
       endif
    enddo

  end subroutine find_times

!------------------------------------------------------------------------
!------------------------------------------------------------------------
  subroutine read_next_trcdata(state, flds, file )
    
    use shr_const_mod, only:pi => shr_const_pi
    use physics_types,only : physics_state
    use ppgrid,           only: pcols, pver, pverp,begchunk,endchunk
    use physconst,        only: rair
    use scamMod
    use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_aer_props, &
                            rad_cnst_get_mode_props, rad_cnst_get_mode_num
    
    implicit none

    type (trfile), intent(inout) :: file
    type (trfld),intent(inout) :: flds(:)
    type(physics_state), intent(in)    :: state(begchunk:endchunk)
    
    integer :: recnos(4),i,f,nflds      ! 
    integer :: cnt4(4)            ! array of counts for each dimension
    integer :: strt4(4)           ! array of starting indices
    integer :: cnt3(3)            ! array of counts for each dimension
    integer :: strt3(3)           ! array of starting indices
    type(file_desc_t) :: fids(4)
    logical :: times_found 

    integer :: cur_yr, cur_mon, cur_day, cur_sec, yr1, yr2, mon, date, sec
    real(r8) :: series1_time, series2_time
    type(file_desc_t) :: fid1, fid2
   
    integer :: nspec,nmodes,n,ii,kk,l1,k,lchnk
    real(r8) :: profile_p(pver),pp(pver),meanP,meanO,sumii,volfrac,specdens,aero_den(6)
    real(r8) :: q_a(3,6)
    real(r8) :: rho(pcols,pver)        ! air density (kg m-3)
    character(len=20) :: aername
    character(len=3) :: arnam(7) = (/'so4','pom','soa','bc ','dst','ncl','num'/)

    nflds = size(flds)
    times_found = .false.

    if(single_column .and. scm_observed_aero) then

      call ver_profile_aero(pp, state(begchunk)%pmid)
      ! The following do loop gets the species properties and calculates the aerosol
      ! mass mixing ratio of each species from observed total number and size distribution
      ! properties in the unit kg/m^3 for the 3 modes. Data is read from the forcing file.
      ! For mode 1 (accumulation mode) q_a(1,1)=q(so4),q_a(1,2)=q(pom),q_a(1,3)=q(soa),
      ! q_a(1,4)=q(bc),q_a(1,5)=q(dst),q_a(1,6)=q(ncl).
      ! For mode 2 (aitken mode) q_a(2,1)=q(so4),q_a(2,2)=q(soa),q_a(2,3)=q(ncl).
      ! For mode 3 (coarse mode) q_a(3,1)=q(dst),q_a(3,2)=q(ncl),q_a(3,3)=q(so4).

      call rad_cnst_get_info(0, nmodes=nmodes)
      do n=1, nmodes
        call rad_cnst_get_info(0, n, nspec=nspec)
        do l1 = 1, nspec
          call rad_cnst_get_aer_props(0, n,l1,density_aer=specdens)
          call rad_cnst_get_aer_props(0, n, l1, aername=aername)
          aero_den(l1)=specdens
          q_a(n,l1) = specdens*scm_div(n,l1)*scm_num(n)*((pi/6.0_r8)*( &
                      scm_dgnum(n)**3)*exp(4.5_r8*(log(scm_std(n))**2) )) 
        enddo
      enddo
      
      do k = 1, pver
        do ii = 1, state(begchunk)%ncol
          rho(ii,k) = state(begchunk)%pmid(ii,k)/(rair*state(begchunk)%t(ii,k))
        enddo
      enddo
    end if

    do while( .not. times_found )
      call find_times( recnos, fids, file%curr_mod_time,file,file%datatimem,file%datatimep, times_found )
      if ( .not. times_found ) then
        call check_files( file, fids, recnos, times_found )
      endif
    enddo

    ! If single column do not interpolate aerosol data, just use the step function
    if(single_column) then
      file%stepTime = .true.
    endif

    if (file%stepTime) then
       file%interp_recs = 1
    else
       file%interp_recs = 2
    end if

    if ( file%fill_in_months ) then

       if( file%datatimep-file%datatimem > file%one_yr ) then

          call get_curr_date(cur_yr, cur_mon, cur_day, cur_sec)

          call set_date_from_time_float(file%datatimem, yr1, mon, date, sec )
          call set_date_from_time_float(file%datatimep, yr2, mon, date, sec )

          call set_time_float_from_date( series1_time, yr1, cur_mon, cur_day, cur_sec )
          call set_time_float_from_date( series2_time, yr2, cur_mon, cur_day, cur_sec )

          fid1 = fids(1)
          fid2 = fids(2)
          file%cyclical = .true.
          call set_cycle_indices( fid1, file%cyc_ndx_beg, file%cyc_ndx_end, yr1)
          call find_times( recnos(1:2), fids(1:2), series1_time, file, file%datatimes(1), file%datatimes(2), times_found )

          if ( .not. times_found ) then
              call endrun('read_next_trcdata: time not found for series1_time')
          endif
          call set_cycle_indices( fid2, file%cyc_ndx_beg, file%cyc_ndx_end, yr2)

          if ( fid1%fh /= fid2%fh ) then
            file%cyc_ndx_beg = file%cyc_ndx_beg + size(file%curr_data_times)
            file%cyc_ndx_end = file%cyc_ndx_end + size(file%curr_data_times)
          endif
          call find_times( recnos(3:4), fids(3:4), series2_time, file, file%datatimes(3), file%datatimes(4), times_found )
          if ( .not. times_found ) then
              call endrun('read_next_trcdata: time not found for series2_time')
          endif
          file%cyclical = .false.
          file%interp_recs = 4

          call set_date_from_time_float( file%datatimes(1), yr1, mon, date, sec )
          call set_time_float_from_date( file%datatimem, cur_yr,  mon, date, sec )
          if (file%datatimes(1) > file%datatimes(2) ) then ! wrap around
            if ( cur_mon == 1 ) then 
               call set_time_float_from_date( file%datatimem, cur_yr-1,  mon, date, sec )
            endif 
          endif

          call set_date_from_time_float( file%datatimes(2), yr1, mon, date, sec )
          call set_time_float_from_date( file%datatimep, cur_yr,  mon, date, sec )
          if (file%datatimes(1) > file%datatimes(2) ) then ! wrap around
            if ( cur_mon == 12 ) then
              call set_time_float_from_date( file%datatimep, cur_yr+1,  mon, date, sec )
            endif 
          endif

       endif

    endif

    !
    ! Set up hyperslab corners
    !
    strt4(:) = 1
    strt3(:) = 1

    do i=1,file%interp_recs

       do f = 1,nflds
          if ( file%zonal_ave ) then
             cnt3(flds(f)%coords(ZA_LATDIM)) = file%nlat
             if (flds(f)%srf_fld) then
                cnt3(flds(f)%coords(ZA_LEVDIM)) = 1
             else
                cnt3(flds(f)%coords(ZA_LEVDIM)) = file%nlev
             endif
             cnt3(flds(f)%coords(ZA_TIMDIM)) = 1
             strt3(flds(f)%coords(ZA_TIMDIM)) = recnos(i)
             call read_za_trc( fids(i), flds(f)%var_id, flds(f)%input(i)%data, strt3, cnt3, file, &
                  (/ flds(f)%order(ZA_LATDIM),flds(f)%order(ZA_LEVDIM) /) )
          else if ( flds(f)%srf_fld ) then
             cnt3( flds(f)%coords(LONDIM)) = file%nlon
             cnt3( flds(f)%coords(LATDIM)) = file%nlat
             cnt3( flds(f)%coords(PS_TIMDIM)) = 1
             strt3(flds(f)%coords(PS_TIMDIM)) = recnos(i)
             call read_2d_trc( fids(i), flds(f)%var_id, flds(f)%input(i)%data(:,1,:), strt3, cnt3, file, &
                 (/ flds(f)%order(LONDIM),flds(f)%order(LATDIM) /) )
          else
             cnt4(flds(f)%coords(LONDIM)) = file%nlon
             cnt4(flds(f)%coords(LATDIM)) = file%nlat
             cnt4(flds(f)%coords(LEVDIM)) = file%nlev
             cnt4(flds(f)%coords(TIMDIM)) = 1
             strt4(flds(f)%coords(TIMDIM)) = recnos(i)
             call read_3d_trc( fids(i), flds(f)%var_id, flds(f)%input(i)%data, strt4, cnt4, file, &
                  (/ flds(f)%order(LONDIM),flds(f)%order(LATDIM),flds(f)%order(LEVDIM) /))

             !
             ! This section sets the observed aersol mass and number mixing ratios in the
             ! appropriate variables. The observed aerosol inforamtion is read from the
             ! forcing file as total number and size distribution parameters. Then the total 
             ! volume is calculated.Using the desnsity of each species, the mass is calculated. 
             ! Finally the number is partition among the each species using the species fraction 
             ! data read from the forcing file.  
             !
             if(single_column .and. scm_observed_aero) then
                kk=index(trim(flds(f)%fldnam),'_')-1
                if(index(trim(flds(f)%fldnam),'1') > 0 .and.index(trim(flds(f)%fldnam),'log') < 1) then
                     if(flds(f)%fldnam(1:kk).eq.arnam(1).and.index(trim(flds(f)%fldnam),'log') < 1) then
                           call replace_aero_data(flds(f)%fldnam,arnam(1),flds(f)%input(i)%data, &
                                rho,pp,q_a(1,1),state(begchunk)%ncol)
                     elseif(flds(f)%fldnam(1:kk).eq.arnam(2).and.index(trim(flds(f)%fldnam),'log') < 1) then
                           call replace_aero_data(flds(f)%fldnam,arnam(2),flds(f)%input(i)%data, &
                                rho,pp,q_a(1,2),state(begchunk)%ncol)
                     elseif(flds(f)%fldnam(1:kk).eq.arnam(3).and.index(trim(flds(f)%fldnam),'log') < 1) then
                           call replace_aero_data(flds(f)%fldnam,arnam(3),flds(f)%input(i)%data, &
                               rho,pp,q_a(1,3),state(begchunk)%ncol)
                    elseif(flds(f)%fldnam(1:kk).eq.arnam(4).and.index(trim(flds(f)%fldnam),'log') < 1) then
                           call replace_aero_data(flds(f)%fldnam,arnam(4),flds(f)%input(i)%data, &
                                rho,pp,q_a(1,4),state(begchunk)%ncol)
                     elseif(flds(f)%fldnam(1:kk).eq.arnam(5).and.index(trim(flds(f)%fldnam),'log') < 1) then
                           call replace_aero_data(flds(f)%fldnam,arnam(5),flds(f)%input(i)%data, &
                                rho,pp,q_a(1,5),state(begchunk)%ncol)
                     elseif(flds(f)%fldnam(1:kk).eq.arnam(6).and.index(trim(flds(f)%fldnam),'log') < 1) then
                           call replace_aero_data(flds(f)%fldnam,arnam(6),flds(f)%input(i)%data, &
                                rho,pp,q_a(1,6),state(begchunk)%ncol)
                     elseif(flds(f)%fldnam(1:kk).eq.arnam(7).and.index(trim(flds(f)%fldnam),'log') < 1) then
                           call replace_aero_data(flds(f)%fldnam,arnam(7),flds(f)%input(i)%data, &
                                rho,pp,scm_num(1),state(begchunk)%ncol)
                     endif
                elseif(index(trim(flds(f)%fldnam),'2') > 0 .and.index(trim(flds(f)%fldnam),'log') < 1) then
                     if(flds(f)%fldnam(1:kk).eq.arnam(1).and.index(trim(flds(f)%fldnam),'log') < 1) then
                           call replace_aero_data(flds(f)%fldnam,arnam(1),flds(f)%input(i)%data, &
                                rho,pp,q_a(2,1),state(begchunk)%ncol)
                     elseif(flds(f)%fldnam(1:kk).eq.arnam(3).and.index(trim(flds(f)%fldnam),'log') < 1) then
                           call replace_aero_data(flds(f)%fldnam,arnam(3),flds(f)%input(i)%data, &
                               rho,pp,q_a(2,2),state(begchunk)%ncol)
                     elseif(flds(f)%fldnam(1:kk).eq.arnam(6).and.index(trim(flds(f)%fldnam),'log') < 1) then
                           call replace_aero_data(flds(f)%fldnam,arnam(6),flds(f)%input(i)%data, &
                                rho,pp,q_a(2,3),state(begchunk)%ncol)
                     elseif(flds(f)%fldnam(1:kk).eq.arnam(7).and.index(trim(flds(f)%fldnam),'log') < 1) then
                           call replace_aero_data(flds(f)%fldnam,arnam(7),flds(f)%input(i)%data, &
                                rho,pp,scm_num(2),state(begchunk)%ncol)
                     endif
                elseif(index(trim(flds(f)%fldnam),'3') > 0 .and.index(trim(flds(f)%fldnam),'log') < 1) then
                     if(flds(f)%fldnam(1:kk).eq.arnam(1).and.index(trim(flds(f)%fldnam),'log') < 1) then
                           call replace_aero_data(flds(f)%fldnam,arnam(1),flds(f)%input(i)%data, &
                                rho,pp,q_a(3,3),state(begchunk)%ncol)
                     elseif(flds(f)%fldnam(1:kk).eq.arnam(5).and.index(trim(flds(f)%fldnam),'log') < 1) then
                           call replace_aero_data(flds(f)%fldnam,arnam(5),flds(f)%input(i)%data, &
                                rho,pp,q_a(3,1),state(begchunk)%ncol)
                     elseif(flds(f)%fldnam(1:kk).eq.arnam(6).and.index(trim(flds(f)%fldnam),'log') < 1) then
                           call replace_aero_data(flds(f)%fldnam,arnam(6),flds(f)%input(i)%data, &
                                rho,pp,q_a(3,2),state(begchunk)%ncol)
                     elseif(flds(f)%fldnam(1:kk).eq.arnam(7).and.index(trim(flds(f)%fldnam),'log') < 1) then
                           call replace_aero_data(flds(f)%fldnam,arnam(7),flds(f)%input(i)%data, &
                                rho,pp,scm_num(3),state(begchunk)%ncol)
                     endif
                endif 
             endif !scm_observed_aero
          endif

       enddo

       if ( file%has_ps ) then
          cnt3(file%ps_coords(LONDIM)) = file%nlon
          cnt3(file%ps_coords(LATDIM)) = file%nlat
          cnt3(file%ps_coords(PS_TIMDIM)) = 1
          strt3(file%ps_coords(PS_TIMDIM)) = recnos(i)
          call read_2d_trc( fids(i), file%ps_id, file%ps_in(i)%data, strt3, cnt3, file, &
               (/ file%ps_order(LONDIM),file%ps_order(LATDIM) /) )
       endif

    enddo

  end subroutine read_next_trcdata

!--------------------------------------------------------------------------------
!This subroutine replaces the climatological aerosol information by the observed
!once after they are read
!
   subroutine   replace_aero_data(aerofulnam,spnam,aero_q_data,rho,pp,q_mix,ncoli)
         use ppgrid,           only:  pcols,pver,begchunk,endchunk
         
         implicit none
         real(r8), intent(inout) :: aero_q_data(pcols,pver,begchunk:endchunk)
         real(r8), intent(in) :: rho(pcols,pver),pp(pver) 
         real(r8) :: sumii,meanO
         real(r8), intent(in) :: q_mix
         character(len=32),intent(in) ::aerofulnam
         character(len=3),intent(in) :: spnam
         character(len=32) ::aerosubnam
         integer, intent(in) :: ncoli 
         integer :: ii,k,countj
               
                if(trim(aerofulnam(1:2)).eq.'bc') then
                  aerosubnam=aerofulnam(1:4)
                else
                  aerosubnam=aerofulnam(1:5)
                endif
           
           if((trim(aerosubnam).eq.(trim(spnam)//'_a').or.trim(aerosubnam).eq.(trim(spnam)//'_c')) &
               .and.index(trim(aerofulnam),'log') < 1) then

                     aero_q_data=q_mix
                   sumii=0._r8
                   countj =0
                  do ii = 1, ncoli
                    do k = 1, pver
                      if(pp(k).gt.0._r8) then
                        countj=countj+1
                      endif
                      if(trim(spnam).ne.'num') then
                          aero_q_data(ii,k,2)=(aero_q_data(ii,k,2)*pp(k)) /rho(ii,k)
                      endif
                         sumii=sumii + aero_q_data(ii,k,2)
                    enddo
                  enddo
                   meanO=sumii/countj
                  do k = 1, pver
                    if(meanO.ne.0.) then
                          aero_q_data(1,k,2)=pp(k)*meanO
                     else                 
                        aero_q_data(1,k,2)=0._r8
                    endif
                  enddo
              endif

   end subroutine replace_aero_data
   
!---------------------------------------------------------------------------
!This subroutine generates a heavyside type profiles for the observed aerosol.
!This setting is constant profile upto 500mb and then exponentially decreasing to zero
!at the top of the atmosphere. The level of initial decay is controlled by
!"initial_val". Larger then -3.5 pushes the decay point up and smaller brings it
!closer to the surface.
!   
   subroutine ver_profile_aero(vertprof_aero, pmid_aero)
              use mo_constants, only : pi
              use ppgrid,       only:  pcols,pver
         
         implicit none
         real(r8), intent(inout) :: vertprof_aero(pver)
         real(r8), intent(in) :: pmid_aero(pcols,pver)
         real(r8) :: initial_val = -3.5_r8
         integer :: k,counti
         
         counti=0   
         do k=1,pver
           if(pmid_aero(1,k)/100..le.500._r8) then
             counti=counti+1
           endif
          if(k==1) then
           vertprof_aero(k)=0._r8
          elseif(k==2) then
            vertprof_aero(k)=1._r8/(1._r8 + exp(-2._r8*initial_val * pi)) 
          else
            vertprof_aero(k)=1._r8/(1._r8 + exp(-2._r8*(initial_val * pi + pi/4._r8*(k-2)))) 
          endif
        enddo 
   end subroutine ver_profile_aero

!------------------------------------------------------------------------


  subroutine read_2d_trc( fid, vid, loc_arr, strt, cnt, file, order )
    use interpolate_data,  only : lininterp_init, lininterp, interp_type, lininterp_finish

    use phys_grid,    only : pcols, begchunk, endchunk, get_ncols_p, get_rlat_all_p, get_rlon_all_p, get_lon_all_p, get_lat_all_p 
    use mo_constants, only : pi
    use dycore,       only: dycore_is		
    use polar_avg,    only: polar_average
    use horizontal_interpolate, only : xy_interp

    implicit none
    type(file_desc_t), intent(in) :: fid
    type(var_desc_t), intent(in) :: vid
    integer, intent(in) :: strt(:), cnt(:), order(2)
    real(r8),intent(out)  :: loc_arr(:,:)
    type (trfile), intent(in) :: file

    real(r8) :: to_lats(pcols), to_lons(pcols), wrk(pcols)
    real(r8), allocatable, target :: wrk2d(:,:)
    real(r8), pointer :: wrk2d_in(:,:)

    integer :: tsize, c, i, j, ierr, ncols
    real(r8), parameter :: zero=0_r8, twopi=2_r8*pi
    type(interp_type) :: lon_wgts, lat_wgts    
    integer :: lons(pcols), lats(pcols)

     nullify(wrk2d_in)
     allocate( wrk2d(cnt(1),cnt(2)), stat=ierr )
     if( ierr /= 0 ) then
        write(iulog,*) 'read_2d_trc: wrk2d allocation error = ',ierr
        call endrun
     end if

     if(order(1)/=1 .or. order(2)/=2 .or. cnt(1)/=file%nlon .or. cnt(2)/=file%nlat) then
        allocate( wrk2d_in(file%nlon, file%nlat), stat=ierr )
        if( ierr /= 0 ) then
           write(iulog,*) 'read_2d_trc: wrk2d_in allocation error = ',ierr
           call endrun
        end if
     end if


    ierr = pio_get_var( fid, vid, strt, cnt, wrk2d )
    if(associated(wrk2d_in)) then
       wrk2d_in = reshape( wrk2d(:,:),(/file%nlon,file%nlat/), order=order )
       deallocate(wrk2d)
    else
       wrk2d_in => wrk2d
    end if

    j=1

! if weighting by latitude, the perform horizontal interpolation by using weight_x, weight_y

   if(file%weight_by_lat) then

      call t_startf('xy_interp')

      do c = begchunk,endchunk
        ncols = get_ncols_p(c)
        call get_lon_all_p(c,ncols,lons)
        call get_lat_all_p(c,ncols,lats)

        call xy_interp(file%nlon,file%nlat,1,plon,plat,pcols,ncols,file%weight_x,file%weight_y,wrk2d_in,loc_arr(:,c-begchunk+1),  &
                            lons,lats,file%count_x,file%count_y,file%index_x,file%index_y) 
      enddo

      call t_stopf('xy_interp')

    else
      do c=begchunk,endchunk
        ncols = get_ncols_p(c)
        call get_rlat_all_p(c, pcols, to_lats)
        call get_rlon_all_p(c, pcols, to_lons)

        call lininterp_init(file%lons, file%nlon, to_lons, ncols, 2, lon_wgts, zero, twopi)
        call lininterp_init(file%lats, file%nlat, to_lats, ncols, 1, lat_wgts)

        call lininterp(wrk2d_in, file%nlon, file%nlat, loc_arr(1:ncols,c-begchunk+1), ncols, lon_wgts, lat_wgts)    
       
        call lininterp_finish(lon_wgts)
        call lininterp_finish(lat_wgts)
      end do
    endif

    if(allocated(wrk2d)) then
       deallocate(wrk2d)
    else
       deallocate(wrk2d_in)
    end if
    if(dycore_is('LR')) call polar_average(loc_arr)
  end subroutine read_2d_trc

!------------------------------------------------------------------------

  subroutine read_za_trc( fid, vid, loc_arr, strt, cnt, file, order )
    use interpolate_data, only : lininterp_init, lininterp, interp_type, lininterp_finish
    use phys_grid,        only : pcols, begchunk, endchunk, get_ncols_p, get_rlat_all_p, get_rlon_all_p	
    use mo_constants,     only : pi
    use dycore,           only : dycore_is		
    use polar_avg,        only : polar_average

    implicit none
    type(file_desc_t), intent(in) :: fid
    type(var_desc_t),  intent(in) :: vid
    integer,           intent(in) :: strt(:), cnt(:)
    integer,           intent(in) :: order(2)
    real(r8),          intent(out):: loc_arr(:,:,:)
    type (trfile),     intent(in) :: file

    type(interp_type) :: lat_wgts
    real(r8) :: to_lats(pcols), to_lons(pcols), wrk(pcols)
    real(r8), allocatable, target :: wrk2d(:,:)
    real(r8), pointer :: wrk2d_in(:,:)
    integer :: c, k, ierr, ncols

     nullify(wrk2d_in)
     allocate( wrk2d(cnt(1),cnt(2)), stat=ierr )
     if( ierr /= 0 ) then
        write(iulog,*) 'read_2d_trc: wrk2d allocation error = ',ierr
        call endrun
     end if

     if(order(1)/=1 .or. order(2)/=2 .or. cnt(1)/=file%nlat .or. cnt(2)/=file%nlev) then
        allocate( wrk2d_in(file%nlat, file%nlev), stat=ierr )
        if( ierr /= 0 ) then
           write(iulog,*) 'read_2d_trc: wrk2d_in allocation error = ',ierr
           call endrun
        end if
     end if


    ierr = pio_get_var( fid, vid, strt, cnt, wrk2d )
    if(associated(wrk2d_in)) then
       wrk2d_in = reshape( wrk2d(:,:),(/file%nlat,file%nlev/), order=order )
       deallocate(wrk2d)
    else
       wrk2d_in => wrk2d
    end if

    do c=begchunk,endchunk
       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, pcols, to_lats)

       call lininterp_init(file%lats, file%nlat, to_lats, ncols, 1, lat_wgts)
       do k=1,file%nlev
          call lininterp(wrk2d_in(:,k), file%nlat, wrk(1:ncols), ncols, lat_wgts)    
          loc_arr(1:ncols,k,c-begchunk+1) = wrk(1:ncols)
       end do
       call lininterp_finish(lat_wgts)
    end do

    if(allocated(wrk2d)) then
       deallocate(wrk2d)
    else
       deallocate(wrk2d_in)
    end if
!    if(dycore_is('LR')) call polar_average(loc_arr)
  end subroutine read_za_trc

!------------------------------------------------------------------------

  subroutine read_3d_trc( fid, vid, loc_arr, strt, cnt, file, order)
    use interpolate_data, only : lininterp_init, lininterp, interp_type, lininterp_finish
    use phys_grid,        only : pcols, begchunk, endchunk, get_ncols_p, get_rlat_all_p, get_rlon_all_p, get_lon_all_p,&
                                 get_lat_all_p 
    use mo_constants,     only : pi
    use dycore,           only : dycore_is		
    use polar_avg,        only : polar_average
    use dycore,           only  : dycore_is
    use horizontal_interpolate, only : xy_interp

    implicit none

    type(file_desc_t), intent(in) :: fid
    type(var_desc_t), intent(in) :: vid
    integer, intent(in) :: strt(:), cnt(:), order(3)
    real(r8),intent(out)  :: loc_arr(:,:,:)
    
    type (trfile), intent(in) :: file

    integer :: i,j,k, astat, c, ncols
    integer :: lons(pcols), lats(pcols)

    integer                     :: jlim(2), jl, ju, ierr
    integer                     :: gndx

    real(r8), allocatable, target :: wrk3d(:,:,:)
    real(r8), pointer :: wrk3d_in(:,:,:)
    real(r8) :: to_lons(pcols), to_lats(pcols)
    real(r8), parameter :: zero=0_r8, twopi=2_r8*pi
    type(interp_type) :: lon_wgts, lat_wgts    

    loc_arr(:,:,:) = 0._r8
    nullify(wrk3d_in)
    allocate(wrk3d(cnt(1),cnt(2),cnt(3)), stat=ierr)
    if( ierr /= 0 ) then
       write(iulog,*) 'read_3d_trc: wrk3d allocation error = ',ierr
       call endrun
    end if

    ierr = pio_get_var( fid, vid, strt, cnt, wrk3d )

    if(order(1)/=1 .or. order(2)/=2 .or. order(3)/=3 .or. &
         cnt(1)/=file%nlon.or.cnt(2)/=file%nlat.or.cnt(3)/=file%nlev) then
       allocate(wrk3d_in(file%nlon,file%nlat,file%nlev),stat=ierr)
       if( ierr /= 0 ) then
          write(iulog,*) 'read_3d_trc: wrk3d allocation error = ',ierr
          call endrun
       end if
       wrk3d_in = reshape( wrk3d(:,:,:),(/file%nlon,file%nlat,file%nlev/), order=order )
       deallocate(wrk3d)
    else
       wrk3d_in => wrk3d
    end if

    j=1

! If weighting by latitude, then perform horizontal interpolation by using weight_x, weight_y

   if(file%weight_by_lat) then

     call t_startf('xy_interp')

     do c = begchunk,endchunk
        ncols = get_ncols_p(c)
        call get_lon_all_p(c,ncols,lons)
        call get_lat_all_p(c,ncols,lats)

        call xy_interp(file%nlon,file%nlat,file%nlev,plon,plat,pcols,ncols,file%weight_x,file%weight_y,wrk3d_in, &
             loc_arr(:,:,c-begchunk+1), lons,lats,file%count_x,file%count_y,file%index_x,file%index_y) 
     enddo

     call t_stopf('xy_interp')

   else
    do c=begchunk,endchunk
       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, pcols, to_lats)
       call get_rlon_all_p(c, pcols, to_lons)

       call lininterp_init(file%lons, file%nlon, to_lons(1:ncols), ncols, 2, lon_wgts, zero, twopi)
       call lininterp_init(file%lats, file%nlat, to_lats(1:ncols), ncols, 1, lat_wgts)


       call lininterp(wrk3d_in, file%nlon, file%nlat, file%nlev, loc_arr(:,:,c-begchunk+1), ncols, pcols, lon_wgts, lat_wgts)    	


       call lininterp_finish(lon_wgts)
       call lininterp_finish(lat_wgts)
    end do
   endif

    if(allocated(wrk3d)) then
       deallocate( wrk3d, stat=astat )
    else
       deallocate( wrk3d_in, stat=astat )
    end if
    if( astat/= 0 ) then
       write(iulog,*) 'read_3d_trc: failed to deallocate wrk3d array; error = ',astat
       call endrun
    endif
    if(dycore_is('LR')) call polar_average(file%nlev, loc_arr)
  end subroutine read_3d_trc

!------------------------------------------------------------------------------

  subroutine interpolate_trcdata( state, flds, file, pbuf2d )
    use mo_util,      only : rebin
    use physics_types,only : physics_state
    use physconst,    only : cday
    use physics_buffer, only : physics_buffer_desc, pbuf_get_field

    implicit none

    type(physics_state), intent(in) :: state(begchunk:endchunk)                 
    type (trfld),        intent(inout) :: flds(:)
    type (trfile),       intent(inout) :: file
    
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)


    real(r8) :: fact1, fact2
    real(r8) :: deltat
    integer :: f,nflds,c,ncol, i,k
    real(r8) :: ps(pcols)
    real(r8) :: datain(pcols,file%nlev)
    real(r8) :: pin(pcols,file%nlev)
    real(r8)            :: model_z(pverp)
    real(r8), parameter :: m2km  = 1.e-3_r8
    real(r8), pointer :: data_out3d(:,:,:)
    real(r8), pointer :: data_out(:,:)
    integer :: chnk_offset

    nflds = size(flds)

    if ( file%interp_recs == 4 ) then
       deltat = file%datatimes(3) - file%datatimes(1)
       fact1 = (file%datatimes(3) - file%datatimem)/deltat
       fact2 = 1._r8-fact1
!$OMP PARALLEL DO PRIVATE (C, NCOL, F)
       do c = begchunk,endchunk
          ncol = state(c)%ncol
          if ( file%has_ps ) then
             file%ps_in(1)%data(:ncol,c) = fact1*file%ps_in(1)%data(:ncol,c) + fact2*file%ps_in(3)%data(:ncol,c) 
          endif
          do f = 1,nflds
             flds(f)%input(1)%data(:ncol,:,c) = fact1*flds(f)%input(1)%data(:ncol,:,c) + fact2*flds(f)%input(3)%data(:ncol,:,c) 
          enddo
       enddo

       deltat = file%datatimes(4) - file%datatimes(2)
       fact1 = (file%datatimes(4) - file%datatimep)/deltat
       fact2 = 1._r8-fact1

!$OMP PARALLEL DO PRIVATE (C, NCOL, F)
       do c = begchunk,endchunk
          ncol = state(c)%ncol
          if ( file%has_ps ) then
             file%ps_in(2)%data(:ncol,c) = fact1*file%ps_in(2)%data(:ncol,c) + fact2*file%ps_in(4)%data(:ncol,c) 
          endif
          do f = 1,nflds
             flds(f)%input(2)%data(:ncol,:,c) = fact1*flds(f)%input(2)%data(:ncol,:,c) + fact2*flds(f)%input(4)%data(:ncol,:,c) 
          enddo
       enddo

    endif
    !-------------------------------------------------------------------------
    !       If file%interp_recs=1 then no time interpolation -- set
    !       fact1=1 and fact2=0 and will just use first value unmodified
    !-------------------------------------------------------------------------

    if (file%interp_recs == 1) then
       fact1=1._r8
       fact2=0._r8
    else
       file%interp_recs = 2

       deltat = file%datatimep - file%datatimem

       if ( file%cyclical .and. (deltat < 0._r8) ) then
          deltat = deltat+file%one_yr
          if ( file%datatimep >= file%curr_mod_time ) then
             fact1 = (file%datatimep - file%curr_mod_time)/deltat
          else
             fact1 = (file%datatimep+file%one_yr - file%curr_mod_time)/deltat
          endif
       else
             fact1 = (file%datatimep - file%curr_mod_time)/deltat
       endif

       ! this assures that FIXED data are b4b on restarts
       if ( file%fixed ) then
          fact1 = dble(int(fact1*cday+.5_r8))/dble(cday)
       endif
       fact2 = 1._r8-fact1
    endif

    chnk_offset=-begchunk+1

    fld_loop: do f = 1,nflds

       if (flds(f)%pbuf_ndx<=0) then
          data_out3d => flds(f)%data(:,:,:)
       endif

!$OMP PARALLEL DO PRIVATE (C, NCOL, PS, I, K, PIN, DATAIN, MODEL_Z, DATA_OUT)
       do c = begchunk,endchunk
          if (flds(f)%pbuf_ndx>0) then
             call pbuf_get_field(pbuf2d, c, flds(f)%pbuf_ndx, data_out)
          else
             data_out => data_out3d(:,:,c+chnk_offset)
          endif
          ncol = state(c)%ncol
          if (file%alt_data) then

             if (fact2 == 0) then  ! This needed as %data is not set if fact2=0 (and lahey compiler core dumps)
                datain(:ncol,:) = fact1*flds(f)%input(nm)%data(:ncol,:,c)
             else
                datain(:ncol,:) = fact1*flds(f)%input(nm)%data(:ncol,:,c) + fact2*flds(f)%input(np)%data(:ncol,:,c) 
             end if
             do i = 1,ncol
                model_z(1:pverp) = m2km * state(c)%zi(i,pverp:1:-1)
                call rebin( file%nlev, pver, file%ilevs, model_z, datain(i,:), data_out(i,:) )
             enddo

          else

             if ( file%nlev>1 ) then
                if ( file%has_ps ) then
                   if (fact2 == 0) then  ! This needed as %data is not set if fact2=0 (and lahey compiler core dumps)
                      ps(:ncol) = fact1*file%ps_in(nm)%data(:ncol,c)
                   else
                      ps(:ncol) = fact1*file%ps_in(nm)%data(:ncol,c) + fact2*file%ps_in(np)%data(:ncol,c) 
                   end if
                   do i = 1,ncol
                      do k = 1,file%nlev
                         pin(i,k) = file%p0*file%hyam(k) + ps(i)*file%hybm(k)
                      enddo
                   enddo
                else
                   do k = 1,file%nlev
                      pin(:,k) = file%levs(k)
                   enddo
                endif
             endif

             if (flds(f)%srf_fld) then
                do i = 1,ncol
                   if (fact2 == 0) then  ! This needed as %data is not set if fact2=0 (and lahey compiler core dumps)
                      data_out(i,1) = &
                           fact1*flds(f)%input(nm)%data(i,1,c)
                   else
                      data_out(i,1) = &
                           fact1*flds(f)%input(nm)%data(i,1,c) + fact2*flds(f)%input(np)%data(i,1,c) 
                   endif
                enddo
             else
                if (fact2 == 0) then  ! This needed as %data is not set if fact2=0 (and lahey compiler core dumps)
                   datain(:ncol,:) = fact1*flds(f)%input(nm)%data(:ncol,:,c)
                else
                   datain(:ncol,:) = fact1*flds(f)%input(nm)%data(:ncol,:,c) + fact2*flds(f)%input(np)%data(:ncol,:,c)
                end if
                if ( file%top_bndry ) then
                   call vert_interp_ub(ncol, file%nlev, file%levs,  datain(:ncol,:), data_out(:ncol,:) )
                else if(file%conserve_column) then
                   call vert_interp_mixrat(ncol,file%nlev,pver,state(c)%pint, &
                        datain, data_out(:,:), &
                        file%p0,ps,file%hyai,file%hybi)
                else
                   call vert_interp(ncol, file%nlev, pin, state(c)%pmid, datain, data_out(:,:) )
                endif
             endif

          endif
       enddo

    enddo fld_loop

  end subroutine interpolate_trcdata

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine get_dimension( fid, dname, dsize, dimid, data )
    implicit none
    type(file_desc_t), intent(inout) :: fid
    character(*), intent(in) :: dname
    integer, intent(out) :: dsize

    integer, optional, intent(out) :: dimid
    real(r8), optional, pointer, dimension(:) :: data

    integer :: vid, ierr, id

    call pio_seterrorhandling( fid, PIO_BCAST_ERROR)
    ierr = pio_inq_dimid( fid, dname, id )
    call pio_seterrorhandling( fid, PIO_INTERNAL_ERROR)

    if ( ierr==PIO_NOERR ) then

       ierr = pio_inq_dimlen( fid, id, dsize )

       if ( present(dimid) ) then
          dimid = id
       endif

       if ( present(data) ) then
          if ( associated(data) ) then
             deallocate(data, stat=ierr)
             if( ierr /= 0 ) then
                write(iulog,*) 'get_dimension: data deallocation error = ',ierr
                call endrun('get_dimension: failed to deallocate data array')
             end if
          endif
          allocate( data(dsize), stat=ierr )
          if( ierr /= 0 ) then
             write(iulog,*) 'get_dimension: data allocation error = ',ierr
             call endrun('get_dimension: failed to allocate data array')
          end if

          ierr =  pio_inq_varid( fid, dname, vid )
          ierr =  pio_get_var( fid, vid, data )
       endif
    else
       dsize = 1
       if ( present(dimid) ) then
          dimid = -1
       endif
    endif

  end subroutine get_dimension

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine set_cycle_indices( fileid, cyc_ndx_beg, cyc_ndx_end, cyc_yr )

    implicit none

    type(file_desc_t), intent(inout)  :: fileid
    integer, intent(out) :: cyc_ndx_beg
    integer, intent(out) :: cyc_ndx_end
    integer, intent(in)  :: cyc_yr

    integer, allocatable , dimension(:) :: dates, datesecs
    integer :: timesize, i, astat, year, ierr
    type(var_desc_T) :: dateid
    call get_dimension( fileid, 'time', timesize )
    cyc_ndx_beg=-1

    allocate( dates(timesize), stat=astat  )
    if( astat/= 0 ) then
       write(*,*) 'set_cycle_indices: failed to allocate dates array; error = ',astat
       call endrun
    end if

    ierr = pio_inq_varid(   fileid, 'date',  dateid  )
    ierr = pio_get_var( fileid, dateid, dates )

    do i=1,timesize
       year = dates(i) / 10000
       if ( year == cyc_yr ) then
          if  (cyc_ndx_beg < 0)  then
             cyc_ndx_beg = i
          endif
          cyc_ndx_end = i
       endif
    enddo
    deallocate( dates, stat=astat  )
    if( astat/= 0 ) then
       write(*,*) 'set_cycle_indices: failed to deallocate dates array; error = ',astat
       call endrun
    end if
    if (cyc_ndx_beg < 0) then
       write(*,*) 'set_cycle_indices: cycle year not found : ' , cyc_yr
       call endrun('set_cycle_indices: cycle year not found')
    endif

  end subroutine set_cycle_indices
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  subroutine open_trc_datafile( fname, path, piofile, times, cyc_ndx_beg, cyc_ndx_end, cyc_yr )

    use ioFileMod,     only: getfil
    use cam_pio_utils, only: cam_pio_openfile

    implicit none

    character(*), intent(in) :: fname
    character(*), intent(in) :: path
    type(file_desc_t), intent(inout) :: piofile
    real(r8), pointer :: times(:)

    integer, optional, intent(out) :: cyc_ndx_beg
    integer, optional, intent(out) :: cyc_ndx_end
    integer, optional, intent(in) :: cyc_yr

    character(len=shr_kind_cl) :: filen, filepath
    integer :: year, month, day, dsize, i, timesize
    integer :: dateid,secid
    integer, allocatable , dimension(:) :: dates, datesecs
    integer :: astat, ierr
    logical :: need_first_ndx

    if (len_trim(path) == 0) then
       filepath = trim(fname)
    else
       filepath = trim(path) // '/' // trim(fname)
    end if
    !
    ! open file and get fileid
    !
    call getfil( filepath, filen, 0 )
    call cam_pio_openfile( piofile, filen, PIO_NOWRITE)
    if(masterproc) write(iulog,*)'open_trc_datafile: ',trim(filen)

    call get_dimension(piofile, 'time', timesize)
    
    if ( associated(times) ) then
       deallocate(times, stat=ierr)
       if( ierr /= 0 ) then
          write(iulog,*) 'open_trc_datafile: data deallocation error = ',ierr
          call endrun('open_trc_datafile: failed to deallocate data array')
       end if
    endif
    allocate( times(timesize), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'open_trc_datafile: data allocation error = ',ierr
       call endrun('open_trc_datafile: failed to allocate data array')
    end if

    allocate( dates(timesize), stat=astat  )
    if( astat/= 0 ) then
       if(masterproc) write(iulog,*) 'open_trc_datafile: failed to allocate dates array; error = ',astat
       call endrun
    end if
    allocate( datesecs(timesize), stat=astat  )
    if( astat/= 0 ) then
       if(masterproc) write(iulog,*) 'open_trc_datafile: failed to allocate datesec array; error = ',astat
       call endrun
    end if

    ierr =  pio_inq_varid( piofile, 'date',    dateid  )
    call pio_seterrorhandling( piofile, PIO_BCAST_ERROR)
    ierr = pio_inq_varid( piofile, 'datesec', secid  )
    call pio_seterrorhandling( piofile, PIO_INTERNAL_ERROR)
    
    if(ierr==PIO_NOERR) then
       ierr = pio_get_var( piofile, secid,  datesecs  )
    else
       datesecs=0
    end if

    ierr =  pio_get_var( piofile, dateid, dates )
    need_first_ndx=.true.

    do i=1,timesize
       year = dates(i) / 10000
       month = mod(dates(i),10000)/100
       day = mod(dates(i),100)
       call set_time_float_from_date( times(i), year, month, day, datesecs(i) )
       if ( present(cyc_yr) ) then
          if ( year == cyc_yr ) then
             if ( present(cyc_ndx_beg) .and. need_first_ndx ) then
                cyc_ndx_beg = i
                need_first_ndx = .false.
             endif
             if ( present(cyc_ndx_end) ) then
                cyc_ndx_end = i
             endif
          endif
       endif
    enddo

    deallocate( dates, stat=astat  )
    if( astat/= 0 ) then
       if(masterproc) write(iulog,*) 'open_trc_datafile: failed to deallocate dates array; error = ',astat
       call endrun
    end if
    deallocate( datesecs, stat=astat  )       
    if( astat/= 0 ) then
       if(masterproc) write(iulog,*) 'open_trc_datafile: failed to deallocate datesec array; error = ',astat
       call endrun
    end if
       
    if ( present(cyc_yr) .and. present(cyc_ndx_beg) ) then
       if (cyc_ndx_beg < 0) then
          write(iulog,*) 'open_trc_datafile: cycle year not found : ' , cyc_yr
          call endrun('open_trc_datafile: cycle year not found')
       endif
    endif

  end subroutine open_trc_datafile

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
  subroutine specify_fields( specifier, fields )

    implicit none

    character(len=*), intent(in) :: specifier(:)
    type(trfld), pointer, dimension(:) :: fields

    integer :: fld_cnt, astat
    integer :: i,j
    character(len=shr_kind_cl) :: str1, str2
    character(len=32), allocatable, dimension(:) :: fld_name,  src_name
    integer :: nflds

    nflds = size(specifier)

    allocate(fld_name(nflds),  src_name(nflds), stat=astat )
    if( astat/= 0 ) then
       write(iulog,*) 'specify_fields: failed to allocate fld_name, src_name arrays; error = ',astat
       call endrun
    end if

    fld_cnt = 0

    count_cnst: do i = 1, nflds

       if ( len_trim( specifier(i) ) == 0 ) then
          exit count_cnst
       endif

       j = scan( specifier(i),':')

       if (j > 0) then
          str1 = trim(adjustl( specifier(i)(:j-1) ))
          str2 = trim(adjustl( specifier(i)(j+1:) ))
          fld_name(i) = trim(adjustl( str1 ))
          src_name(i) = trim(adjustl( str2 ))
       else
          fld_name(i) = trim(adjustl( specifier(i) ))
          src_name(i) = trim(adjustl( specifier(i) ))
       endif

       fld_cnt = fld_cnt + 1

    enddo count_cnst

    if( fld_cnt < 1 ) then
       nullify(fields)
       return
    end if

    !-----------------------------------------------------------------------
    ! 	... allocate field type array
    !-----------------------------------------------------------------------
    allocate( fields(fld_cnt), stat=astat )
    if( astat/= 0 ) then
       write(iulog,*) 'specify_fields: failed to allocate fields array; error = ',astat
       call endrun
    end if

    do i = 1,fld_cnt
       fields(i)%fldnam = fld_name(i)
       fields(i)%srcnam = src_name(i)
    enddo

    deallocate(fld_name, src_name)

  end subroutine specify_fields

!------------------------------------------------------------------------------

  subroutine init_trc_restart( whence, piofile, tr_file )

    implicit none
    character(len=*), intent(in) :: whence
    type(file_desc_t), intent(inout) :: piofile
    type(trfile), intent(inout) :: tr_file
 
    character(len=32) :: name
    integer :: ioerr, mcdimid, maxlen
 

    ! Dimension should already be defined in restart file
    call pio_seterrorhandling(pioFile, PIO_BCAST_ERROR)
    ioerr = pio_inq_dimid(pioFile,'max_chars', mcdimid)
    call pio_seterrorhandling(pioFile, PIO_INTERNAL_ERROR)
    ! but define it if nessasary
    if(ioerr/= PIO_NOERR) then
       ioerr = pio_def_dim(pioFile, 'max_chars', SHR_KIND_CL, mcdimid)
    end if

    if(len_trim(tr_file%curr_filename)>1) then
       allocate(tr_file%currfnameid)
       name = trim(whence)//'_curr_fname'
       ioerr = pio_def_var(pioFile, name,pio_char,  (/mcdimid/), tr_file%currfnameid)
       ioerr = pio_put_att(pioFile, tr_file%currfnameid, 'offset_time', tr_file%offset_time)
       maxlen = len_trim(tr_file%curr_filename)
       ioerr = pio_put_att(pioFile, tr_file%currfnameid, 'actual_len', maxlen)	
    else
       nullify(tr_file%currfnameid)
    end if

    if(len_trim(tr_file%next_filename)>1) then
       allocate(tr_file%nextfnameid)
       name = trim(whence)//'_next_fname'
       ioerr = pio_def_var(pioFile, name,pio_char,  (/mcdimid/), tr_file%nextfnameid)
       maxlen = len_trim(tr_file%next_filename)
       ioerr = pio_put_att(pioFile, tr_file%nextfnameid, 'actual_len', maxlen)
    else
       nullify(tr_file%nextfnameid)
    end if
  end subroutine init_trc_restart
!-------------------------------------------------------------------------
! writes file names to restart file
!-------------------------------------------------------------------------
  subroutine write_trc_restart( piofile, tr_file )

    implicit none

    type(file_desc_t), intent(inout) :: piofile
    type(trfile), intent(inout) :: tr_file

    integer :: ioerr, slen   ! error status
    if(associated(tr_file%currfnameid)) then
       ioerr = pio_put_var(pioFile, tr_file%currfnameid, tr_file%curr_filename)
       deallocate(tr_file%currfnameid)
       nullify(tr_file%currfnameid)
    end if
    if(associated(tr_file%nextfnameid)) then
       ioerr = pio_put_var(pioFile, tr_file%nextfnameid, tr_file%next_filename)
       deallocate(tr_file%nextfnameid)
       nullify(tr_file%nextfnameid)
    end if
  end subroutine write_trc_restart

!-------------------------------------------------------------------------
! reads file names from restart file
!-------------------------------------------------------------------------
  subroutine read_trc_restart( whence, piofile, tr_file )

    implicit none

    character(len=*), intent(in) :: whence
    type(file_desc_t), intent(inout) :: piofile
    type(trfile), intent(inout) :: tr_file
    type(var_desc_t) :: vdesc
    character(len=64) :: name
    integer :: ioerr   ! error status
    integer :: slen

    call PIO_SetErrorHandling(piofile, PIO_BCAST_ERROR)
    name = trim(whence)//'_curr_fname'
    ioerr = pio_inq_varid(piofile, name, vdesc)
    if(ioerr==PIO_NOERR) then
       tr_file%curr_filename=' '
       ioerr = pio_get_att(piofile, vdesc, 'offset_time', tr_file%offset_time)
       ioerr = pio_get_att(piofile, vdesc, 'actual_len', slen)
       ioerr = pio_get_var(piofile, vdesc, tr_file%curr_filename)
       if(slen<SHR_KIND_CL) tr_file%curr_filename(slen+1:)=' '
    end if

    name = trim(whence)//'_next_fname'
    ioerr = pio_inq_varid(piofile, name, vdesc)
    if(ioerr==PIO_NOERR) then
       tr_file%next_filename=' '
       ioerr = pio_get_att(piofile, vdesc, 'actual_len', slen)
       ioerr = pio_get_var(piofile, vdesc, tr_file%next_filename)
       if(slen<SHR_KIND_CL) tr_file%next_filename(slen+1:)=' '
    end if
    call PIO_SetErrorHandling(piofile, PIO_INTERNAL_ERROR)



  end subroutine read_trc_restart
!------------------------------------------------------------------------------
   subroutine vert_interp_mixrat( ncol, nsrc, ntrg, trg_x, src, trg, p0, ps, hyai, hybi)
  
    implicit none

    integer, intent(in)   :: ncol 
    integer, intent(in)   :: nsrc                  ! dimension source array
    integer, intent(in)   :: ntrg                  ! dimension target array
    real(r8)              :: src_x(nsrc+1)         ! source coordinates
    real(r8), intent(in)      :: trg_x(pcols,ntrg+1)         ! target coordinates
    real(r8), intent(in)      :: src(pcols,nsrc)             ! source array
    real(r8), intent(out)     :: trg(pcols,ntrg)             ! target array

    real(r8) :: ps(pcols), p0, hyai(nsrc+1), hybi(nsrc+1)
    !---------------------------------------------------------------
    !   ... local variables
    !---------------------------------------------------------------
    integer  :: i, j, n
    integer  :: sil
    real(r8)     :: tl, y
    real(r8)     :: bot, top
   
  
    
    do n = 1,ncol
    
    do i=1,nsrc+1
     src_x(i) = p0*hyai(i)+ps(n)*hybi(i)
    enddo

    do i = 1, ntrg
       tl = trg_x(n,i+1)
       if( (tl.gt.src_x(1)).and.(trg_x(n,i).lt.src_x(nsrc+1)) ) then
          do sil = 1,nsrc
             if( (tl-src_x(sil))*(tl-src_x(sil+1)).le.0.0_r8 ) then
                exit
             end if
          end do

          if( tl.gt.src_x(nsrc+1)) sil = nsrc

          y = 0.0_r8
          bot = min(tl,src_x(nsrc+1))   
          top = trg_x(n,i)
          do j = sil,1,-1
           if( top.lt.src_x(j) ) then
             y = y+(bot-src_x(j))*src(n,j)
            bot = src_x(j)
           else
            y = y+(bot-top)*src(n,j)
            exit
           endif
          enddo
          trg(n,i) = y
       else
        trg(n,i) = 0.0_r8
       end if
    end do

    if( trg_x(n,ntrg+1).lt.src_x(nsrc+1) ) then
     top = trg_x(n,ntrg+1)
     bot = src_x(nsrc+1)
     y = 0.0_r8
     do j=nsrc,1,-1
      if( top.lt.src_x(j) ) then
       y = y+(bot-src_x(j))*src(n,j)
       bot = src_x(j)
      else
       y = y+(bot-top)*src(n,j)
       exit
      endif
     enddo
     trg(n,ntrg) = trg(n,ntrg)+y
    endif

! turn mass into mixing ratio 
    do i=1,ntrg
     trg(n,i) = trg(n,i)/(trg_x(n,i+1)-trg_x(n,i))
    enddo
    
    enddo

   end subroutine vert_interp_mixrat
!------------------------------------------------------------------------------
  subroutine vert_interp( ncol, levsiz, pin, pmid, datain, dataout )
    !-------------------------------------------------------------------------- 
    ! 
    ! Interpolate data from current time-interpolated values to model levels
    !--------------------------------------------------------------------------
    implicit none
    ! Arguments
    !
    integer,  intent(in)  :: ncol                ! number of atmospheric columns
    integer,  intent(in)  :: levsiz
    real(r8), intent(in)  :: pin(pcols,levsiz)
    real(r8), intent(in)  :: pmid(pcols,pver)          ! level pressures 
    real(r8), intent(in)  :: datain(pcols,levsiz)
    real(r8), intent(out) :: dataout(pcols,pver)     

    !
    ! local storage
    !

    integer ::  i                   ! longitude index
    integer ::  k, kk, kkstart      ! level indices
    integer ::  kupper(pcols)       ! Level indices for interpolation
    real(r8) :: dpu                ! upper level pressure difference
    real(r8) :: dpl                ! lower level pressure difference



    !--------------------------------------------------------------------------
    !
    ! Initialize index array
    !
    do i=1,ncol
       kupper(i) = 1
    end do

    do k=1,pver
       !
       ! Top level we need to start looking is the top level for the previous k
       ! for all column points
       !
       kkstart = levsiz
       do i=1,ncol
          kkstart = min0(kkstart,kupper(i))
       end do
       !
       ! Store level indices for interpolation
       !
       do kk=kkstart,levsiz-1
          do i=1,ncol
             if (pin(i,kk).lt.pmid(i,k) .and. pmid(i,k).le.pin(i,kk+1)) then
                kupper(i) = kk
             end if
          end do
       end do
       ! interpolate or extrapolate...
       do i=1,ncol
          if (pmid(i,k) .lt. pin(i,1)) then
             dataout(i,k) = datain(i,1)*pmid(i,k)/pin(i,1)
          else if (pmid(i,k) .gt. pin(i,levsiz)) then
             dataout(i,k) = datain(i,levsiz)
          else
             dpu = pmid(i,k) - pin(i,kupper(i))
             dpl = pin(i,kupper(i)+1) - pmid(i,k)
             dataout(i,k) = (datain(i,kupper(i) )*dpl + &
                  datain(i,kupper(i)+1)*dpu)/(dpl + dpu)
          end if
       end do
    end do


  end subroutine vert_interp

!------------------------------------------------------------------------------
  subroutine vert_interp_ub( ncol, nlevs, plevs,  datain, dataout )
    use ref_pres, only : ptop_ref


    !----------------------------------------------------------------------- 
    ! 
    ! Interpolate data from current time-interpolated values to top interface pressure
    !  -- from mo_tgcm_ubc.F90
    !--------------------------------------------------------------------------
    implicit none
    ! Arguments
    !
    integer,  intent(in)  :: ncol
    integer,  intent(in)  :: nlevs
    real(r8), intent(in)  :: plevs(nlevs)
    real(r8), intent(in)  :: datain(ncol,nlevs)
    real(r8), intent(out) :: dataout(ncol)   

    !
    ! local variables
    !
    integer  :: i,ku,kl,kk
    real(r8) :: pinterp, delp
    
    pinterp = ptop_ref

    if( pinterp <= plevs(1) ) then
       kl = 1
       ku = 1
       delp = 0._r8
    else if( pinterp >= plevs(nlevs) ) then
       kl = nlevs
       ku = nlevs
       delp = 0._r8
    else

       do kk = 2,nlevs
          if( pinterp <= plevs(kk) ) then
             ku = kk
             kl = kk - 1
             delp = log( pinterp/plevs(kk) ) / log( plevs(kk-1)/plevs(kk) )
             exit
          end if
       end do

    end if

    do i = 1,ncol
       dataout(i) = datain(i,kl) + delp * (datain(i,ku) - datain(i,kl))
    end do

  end subroutine vert_interp_ub
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  subroutine advance_file(file)

    !------------------------------------------------------------------------------
    !   This routine advances to the next file
    !------------------------------------------------------------------------------

    use shr_sys_mod, only: shr_sys_system
    use ioFileMod, only: getfil

    implicit none

    type(trfile), intent(inout) :: file

    !-----------------------------------------------------------------------
    !   local variables
    !-----------------------------------------------------------------------
    character(len=shr_kind_cl) :: ctmp
    character(len=shr_kind_cl) :: loc_fname   
    integer            :: istat, astat

    !-----------------------------------------------------------------------
    !   close current file ...
    !-----------------------------------------------------------------------
    call pio_closefile( file%curr_fileid )

    !-----------------------------------------------------------------------
    !   remove if requested
    !-----------------------------------------------------------------------
    if( file%remove_trc_file ) then
       call getfil( file%curr_filename, loc_fname, 0 )
       write(iulog,*) 'advance_file: removing file = ',trim(loc_fname) 
       ctmp = 'rm -f ' // trim(loc_fname) 
       write(iulog,*) 'advance_file: fsystem issuing command - '
       write(iulog,*) trim(ctmp)
       call shr_sys_system( ctmp, istat )
    end if
   
    !-----------------------------------------------------------------------
    !   Advance the filename and file id
    !-----------------------------------------------------------------------
    file%curr_filename = file%next_filename
    file%curr_fileid = file%next_fileid
   
    !-----------------------------------------------------------------------
    !   Advance the curr_data_times
    !-----------------------------------------------------------------------
    deallocate( file%curr_data_times, stat=astat )
    if( astat/= 0 ) then
       write(iulog,*) 'advance_file: failed to deallocate file%curr_data_times array; error = ',astat
       call endrun
    end if
    allocate( file%curr_data_times( size( file%next_data_times ) ), stat=astat )
    if( astat/= 0 ) then
       write(iulog,*) 'advance_file: failed to allocate file%curr_data_times array; error = ',astat
       call endrun
    end if
    file%curr_data_times(:) = file%next_data_times(:)
    
    !-----------------------------------------------------------------------
    !   delete information about next file (as was just assigned to current)
    !-----------------------------------------------------------------------
    file%next_filename = ''
    
    deallocate( file%next_data_times, stat=astat )
    if( astat/= 0 ) then
       write(iulog,*) 'advance_file: failed to deallocate file%next_data_times array; error = ',astat
       call endrun
    end if
    nullify( file%next_data_times )

  end subroutine advance_file

end module tracer_data

!-----------------------------------------------------------------------
! $Id: output_netcdf.F90 7169 2014-08-05 21:42:25Z dschanen@uwm.edu $
!===============================================================================
module output_netcdf
#ifdef NETCDF

! Description:
!   Functions and subroutines for writing NetCDF files

! References:
!   <http://www.unidata.ucar.edu/software/netcdf/docs/>
!-------------------------------------------------------------------------------

  implicit none

  public :: open_netcdf, write_netcdf, close_netcdf

  private :: define_netcdf, write_grid, first_write, format_date

  ! Constant parameters
  ! This will truncate all timesteps smaller than 1 mn to a minute for 
  ! the purposes of viewing the data in grads
  logical, parameter, private :: &
    l_grads_netcdf_boost_ts = .false. 

  private ! Default scope

  contains
!-------------------------------------------------------------------------------
  subroutine open_netcdf( nlat, nlon, fdir, fname, ia, iz, zgrid,  & 
                          day, month, year, rlat, rlon, & 
                          time, dtwrite, nvar, ncf )

! Description:
!   Defines the structure used to reference the file `ncf'

! References:
!   None
!-------------------------------------------------------------------------------
    use netcdf, only: & 
      NF90_CLOBBER, & ! Variable(s)
      NF90_NOERR,   & 
      nf90_create,  & ! Procedure
      nf90_strerror

    use stat_file_module, only: & 
      stat_file ! Type

    use clubb_precision, only:  & 
      time_precision, & ! Variable(s)
      core_rknd

    use constants_clubb, only:  & 
      fstderr, & ! Variable(s)
      sec_per_min

    use stats_variables, only: &
      l_allow_small_stats_tout

    implicit none

    ! Input Variables
    character(len=*), intent(in) ::  & 
      fdir,   & ! Directory name of file
      fname     ! File name

    integer, intent(in) ::  & 
      nlat, nlon,       & ! Number of points in the X and Y
      day, month, year, & ! Time
      ia, iz,           & ! First and last grid point
      nvar                ! Number of variables

    real( kind = core_rknd ), dimension(nlat), intent(in) ::  & 
      rlat ! Latitudes   [degrees_E]

    real( kind = core_rknd ), dimension(nlon), intent(in) ::  & 
      rlon ! Longitudes  [degrees_N]

    real( kind = core_rknd ), intent(in) :: & 
      dtwrite ! Time between write intervals   [s]

    real( kind = time_precision ), intent(in) ::  & 
     time   ! Current time                    [s]

    real( kind = core_rknd ), dimension(:), intent(in) ::  & 
      zgrid  ! The model grid                  [m]

    ! Input/output Variables
    type (stat_file), intent(inout) :: ncf

    ! Local Variables
    integer :: stat  ! Error status
    integer :: k     ! Array index

    ! ---- Begin Code ----

    ncf%nvar    = nvar

    ! If there is no data to write, then return
    if ( ncf%nvar == 0 ) then
      return
    end if

    ! Initialization for NetCDF
    ncf%l_defined = .false.

    ! Define file (compatability with GrADS writing)
    ncf%fdir   = fdir
    ncf%fname  = fname
    ncf%ia     = ia
    ncf%iz     = iz
    ncf%day    = day
    ncf%month  = month
    ncf%year   = year
    ncf%nlat   = nlat
    ncf%nlon   = nlon
    ncf%time   = time

    ncf%dtwrite = dtwrite

    ! Check to make sure the timestep is appropriate. The GrADS program does not support an
    ! output timestep less than 1 minute.  Other programs can read netCDF files like this
    if ( dtwrite < sec_per_min ) then
      write(fstderr,*) "Warning: GrADS program requires an output timestep of at least &
                       &one minute, but the requested output timestep &
                       &(stats_tout) is less than one minute."
      if ( .not. l_allow_small_stats_tout ) then
        write(fstderr,*) "To override this warning, set l_allow_small_stats_tout = &
                         &.true. in the stats_setting namelist in the &
                         &appropriate *_model.in file."
        stop "Fatal error in open_netcdf"
      end if
    end if ! dtwrite < sec_per_min

    ! From open_grads.
    ! This probably for the case of a reversed grid as in COAMPS
    if ( ia <= iz ) then
      do k=1,iz-ia+1
        ncf%z(k) = zgrid(ia+k-1)
      end do
    else ! Always this for CLUBB
      do k=1,ia-iz+1
        ncf%z(k) = zgrid(ia-k+1)
      end do
    end if

    allocate( ncf%rlat(1:nlat), ncf%rlon(1:nlon) )

    ncf%rlat = rlat
    ncf%rlon = rlon

    ! Create NetCDF dataset: enter define mode
    stat = nf90_create( path = trim( fdir )//trim( fname )//'.nc',  & 
                        cmode = NF90_CLOBBER,  & ! overwrite existing file
                        ncid = ncf%iounit )
    if ( stat /= NF90_NOERR ) then
      write(unit=fstderr,fmt=*) "Error opening file: ",  & 
        trim( fdir )//trim( fname )//'.nc', & 
        trim( nf90_strerror( stat ) )
      stop "Fatal Error"
    end if

    call define_netcdf( ncf%iounit, ncf%nlat, ncf%nlon, ncf%iz, & ! In
                        ncf%day, ncf%month, ncf%year, ncf%time, & ! In 
                        ncf%LatDimId, ncf%LongDimId, ncf%AltDimId, ncf%TimeDimId, &  ! Out
                        ncf%LatVarId, ncf%LongVarId, ncf%AltVarId, ncf%TimeVarId ) ! Out

    return
  end subroutine open_netcdf

!-------------------------------------------------------------------------------

  subroutine write_netcdf( ncf )

! Description:
!   Writes some data to the NetCDF dataset, but doesn't close it.
!
! References:
!   None   
!-------------------------------------------------------------------------------

    use netcdf, only: & 
        NF90_NOERR,  & ! Variable(s)
        nf90_put_var,  & ! Procedure
        nf90_strerror

    use stat_file_module, only: & 
        stat_file ! Variable

    use constants_clubb, only:  & 
        fstderr, & ! Variable
        sec_per_min

    use clubb_precision, only: &
      time_precision ! Constant(s)

    implicit none

    ! Input
    type (stat_file), intent(inout) :: ncf    ! The file

    ! Local Variables
    integer, dimension(:), allocatable :: stat ! Error status
    real(kind=8), dimension(1) :: time         ! Time          [s]

    integer :: i ! Array index

    ! ---- Begin Code ----

    ! If there is no data to write, then return
    if ( ncf%nvar == 0 ) then
      return
    end if

    ncf%ntimes = ncf%ntimes + 1

    if ( .not. ncf%l_defined ) then
      call first_write( ncf ) ! finalize the variable definitions
      call write_grid( ncf )  ! define lat., long., and grid
      ncf%l_defined = .true.
    end if

    allocate( stat( ncf%nvar ) )
    if ( l_grads_netcdf_boost_ts ) then
      time = real( nint( real(ncf%ntimes, kind=time_precision) &
                            * real(ncf%dtwrite / sec_per_min, time_precision) ), &
                              kind=time_precision ) ! minutes(rounded)
    else
      time = real( ncf%ntimes, kind=time_precision ) &
           * real( ncf%dtwrite, kind=time_precision )  ! seconds
    end if

    stat(1) = nf90_put_var( ncid=ncf%iounit, varid=ncf%TimeVarId,  & 
                            values=time(1), start=(/ncf%ntimes/) )
    if ( stat(1) /= NF90_NOERR ) then
      stop "time variable nf90_put_var failed"
    end if

    do i = 1, ncf%nvar, 1
      stat(i)  & 
      = nf90_put_var( ncid=ncf%iounit, varid=ncf%var(i)%indx,  & 
                      values=ncf%var(i)%ptr(:,:,ncf%ia:ncf%iz),  & 
                      start=(/1,1,1,ncf%ntimes/), & 
                      count=(/ncf%nlon,ncf%nlat,ncf%iz,1/) )

    end do ! i=1..nvar

    if ( any (stat /= NF90_NOERR ) ) then
      do i=1,ncf%nvar,1
        if( stat(i) /= NF90_NOERR ) then
          write(unit=fstderr,fmt=*) ncf%var(i)%name,  & 
            trim( nf90_strerror( stat(i) ) )
        end if
      end do
      stop "nf90_put_var error"
    end if


    deallocate( stat )

    return
  end subroutine write_netcdf

!-------------------------------------------------------------------------------
  subroutine define_netcdf( ncid, nlat, nlon, iz, &
                            day, month, year, time, & 
                            LatDimId, LongDimId, AltDimId, TimeDimId, & 
                            LatVarId, LongVarId, AltVarId, TimeVarId )

! Description:
!   Used internally to create a definition for the NetCDF dataset
!
! References:
!   None
!-------------------------------------------------------------------------------
    use netcdf, only: & 
      NF90_NOERR,   & ! Constants
      NF90_FLOAT, & 
      NF90_DOUBLE, & 
      NF90_UNLIMITED

    use netcdf, only: & 
      nf90_def_dim,  & ! Functions
      nf90_strerror, & 
      nf90_def_var, & 
      nf90_put_att

    use clubb_precision, only:  & 
      time_precision ! Variable(s)

    use constants_clubb, only:  & 
      fstderr ! Variable(s)

    implicit none

    integer, intent(in) ::  & 
      nlat,   & ! Number of points in the N/S direction
      nlon      ! Number of points in the E/W direction

    ! Input Variables
    integer, intent(in) ::  & 
      day, month, year,  & ! Time of year
      ncid,              & ! Number used by NetCDF for ref. the file
      iz                   ! Dimension in z

    real(kind=time_precision), intent(in) ::  & 
      time    ! Current model time [s]

    ! Output Variables
    integer, intent(out) ::  & 
      LatDimId, LongDimId, AltDimId, TimeDimId  ! NetCDF id's for dimensions

    ! NetCDF id's for data (e.g. longitude) associated with each dimension
    integer, intent(out) ::  & 
      LatVarId, LongVarId, AltVarId, TimeVarId

    ! Local variables
    integer :: stat
    character(len=35) :: TimeUnits

    ! ---- Begin Code ----

    ! Define the dimensions for the variables
    stat = nf90_def_dim( ncid, "longitude", nlon, LongDimId )

    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining longitude: ", & 
        trim( nf90_strerror( stat ) )
      stop
    end if

    stat =  nf90_def_dim( ncid, "latitude", nlat, LatDimId )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining latitude: ", & 
        trim( nf90_strerror( stat ) )
      stop
    end if

    stat = nf90_def_dim( ncid, "altitude", iz, AltDimId )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining altitude: ", & 
      trim( nf90_strerror( stat ) )
      stop
    end if

    stat =  nf90_def_dim( ncid, "time", NF90_UNLIMITED, TimeDimId )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining time: ", & 
        trim( nf90_strerror( stat ) )
      stop
    end if

    ! Define the initial variables for the dimensions
    ! Longitude = deg_E = X
    stat = nf90_def_var( ncid, "longitude", NF90_FLOAT, & 
                         (/LongDimId/), LongVarId )

    ! Latitude = deg_N = Y
    stat = nf90_def_var( ncid, "latitude", NF90_FLOAT, & 
                         (/LatDimId/), LatVarId )

    ! Altitude = meters above the surface = Z
    stat = nf90_def_var( ncid, "altitude", NF90_FLOAT, & 
                        (/AltDimId/), AltVarId )

    ! grads2nc stores time as a double prec. value, so we follow that
    stat = nf90_def_var( ncid, "time", NF90_DOUBLE, & 
                         (/TimeDimId/), TimeVarId )

    ! Assign attribute values

    ! Time attribute
    stat = nf90_put_att( ncid, TimeVarId, "cartesian_axis", "T" )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining time: ", trim( nf90_strerror( stat ) )
      stop
    end if

    call format_date( day, month, year, time, TimeUnits )

    stat = nf90_put_att( ncid, TimeVarId, "units", TimeUnits )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining time: ", trim( nf90_strerror( stat ) )
      stop
    end if

    stat = nf90_put_att( ncid, TimeVarId, "ipositive", 1 )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining time: ", trim( nf90_strerror( stat ) )
      stop
    end if

    stat = nf90_put_att( ncid, TimeVarId, "calendar_type", "Gregorian" )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining time", trim( nf90_strerror( stat ) )
      stop
    end if

    ! Define Location
    ! X & Y coordinates
    stat = nf90_put_att( ncid, LongVarId, "cartesian_axis", "X" )

    stat = nf90_put_att( ncid, LongVarId, "units",  "degrees_E" )

    stat = nf90_put_att( ncid, LongVarId, "ipositive",  1 )

    stat = nf90_put_att( ncid, LatVarId, "cartesian_axis",  "Y" )

    stat = nf90_put_att( ncid, LatVarId, "units", "degrees_N" )

    stat = nf90_put_att( ncid, LatVarId, "ipositive", 1 )

    ! Altitude, Z coordinate
    stat = nf90_put_att( ncid, AltVarId, "cartesian_axis",  "Z" )

    stat = nf90_put_att( ncid, AltVarId, "units", "meters" )

    stat = nf90_put_att( ncid, AltVarId, "positive",  "up" )

    stat = nf90_put_att( ncid, AltVarId, "ipositive", 1 )

    return
  end subroutine define_netcdf

!-------------------------------------------------------------------------------
  subroutine close_netcdf( ncf )

! Description:
!   Close a previously opened stats file.

! Notes:
!   I assume nf90_close() exists so that the NetCDF libraries can do a
!   form of buffered I/O, but I don't know the implementation
!   details. -dschanen
!-------------------------------------------------------------------------------

    use stat_file_module, only: & 
        stat_file ! Type

    use netcdf, only: & 
        NF90_NOERR,  & ! Variable
        nf90_close,  & ! Procedure(s)
        nf90_strerror

    use constants_clubb, only:  & 
        fstderr  ! Variable

    implicit none

    ! Input/Output Variables
    type (stat_file), intent(inout) :: ncf

    ! Local Variables
    integer :: stat

    ! ---- Begin Code ----

    ! If there is no data to write, then return
    if ( ncf%nvar == 0 ) then
      return
    end if

    stat = nf90_close( ncf%iounit )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error closing file "//  & 
        trim( ncf%fname )//": ", trim( nf90_strerror( stat ) )
      stop "Fatal error"
    end if

    return
  end subroutine close_netcdf

!-------------------------------------------------------------------------------
  subroutine first_write( ncf )

! Description:
!   Used on the first call to write_nc to finalize definitions
!   for the dataset, including the attributes for variable records.
! References:
!   None
!-------------------------------------------------------------------------------

    use netcdf, only: & 
      NF90_NOERR,  & ! Constants
      NF90_FLOAT,  &
      NF90_DOUBLE, & 
      NF90_GLOBAL, &
      nf90_def_var,  & ! Procedure(s)
      nf90_strerror, & 
      nf90_put_att, & 
      nf90_enddef

    use stat_file_module, only: &
      stat_file ! Derived type

    use constants_clubb, only:  &
      fstderr ! Variable

    use parameters_model, only: &
      T0, &       ! Real variables
      ts_nudge, &
      sclr_tol    ! Real array variable

    use parameters_tunable, only: &
      params_list ! Variable names (characters)

    use parameters_tunable, only: &
      get_parameters ! Subroutine

    use parameter_indices, only: &
      nparams ! Integer

    use model_flags, only: &
      l_pos_def, &
      l_hole_fill, &
      l_clip_semi_implicit, &
      l_standard_term_ta, &
      l_single_C2_Skw, &
      l_gamma_Skw, &
      l_uv_nudge, &
      l_tke_aniso

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: date_and_time, huge, selected_real_kind, size, any, trim

    ! Enabling l_output_file_run_date allows the date and time that the netCDF
    ! output file is created to be included in the netCDF output file.
    ! Disabling l_output_file_run_date means that this information will not be
    ! included in the netCDF output file.  The advantage of disabling this
    ! output is that it allows for a check for binary differences between two
    ! netCDF output files.
    logical, parameter :: &
      l_output_file_run_date = .false.

    ! Input/Output Variables
    type (stat_file), intent(inout) :: ncf

    ! Local Variables
    integer, dimension(:), allocatable :: stat
    
    integer :: netcdf_precision ! Level of precision for netCDF output

    real( kind = core_rknd ), dimension(nparams) :: params ! Tunable parameters

    integer :: i     ! Array index
    logical :: l_error ! Error stat

    character(len=10) :: current_time
    character(len=8)  :: current_date
    ! Range for NetCDF variables
    real(kind=4), dimension(2) :: var_range

    ! Dimensions for variables
    integer, dimension(4) :: var_dim


!-------------------------------------------------------------------------------
!      Typical valid ranges (IEEE 754)

!      real(kind=4): +/- 3.4028235E+38
!      real(kind=8): +/- 1.797693134862316E+308
!      real(kind=16):+/- 1.189731495357231765085759326628007E+4932

!      We use a 4 byte data model for NetCDF and GrADS to save disk space
!-------------------------------------------------------------------------------

    ! ---- Begin Code ----

    var_range(1) = -huge( var_range(1) )
    var_range(2) =  huge( var_range(2) )

! var_range = (/ -1.e31, 1.e31 /)

! Explanation:  The NetCDF documentation claims the NF90_UNLIMITED
!   variable should be the first dimension, but def_var is somehow
!   inverted and requires the opposite.  After writing, these
!   dimensions are all in the opposite order of this in the file.
!   -dschanen

    var_dim(1) = ncf%LongDimId ! X
    var_dim(2) = ncf%LatDimId  ! Y
    var_dim(3) = ncf%AltDimId  ! Z
    var_dim(4) = ncf%TimeDimId ! The NF90_UNLIMITED dimension

    allocate( stat( ncf%nvar ) )

    l_error = .false.


    select case (core_rknd)
      case ( selected_real_kind( p=5 ) )
        netcdf_precision = NF90_FLOAT
      case ( selected_real_kind( p=12 ) )
        netcdf_precision = NF90_DOUBLE
      case default
        netcdf_precision = NF90_DOUBLE
    end select

    do i = 1, ncf%nvar, 1
!     stat(i) = nf90_def_var( ncf%iounit, trim( ncf%var(i)%name ), &
!                  NF90_FLOAT, (/ncf%TimeDimId, ncf%AltDimId, &
!                  ncf%LatDimId, ncf%LongDimId/), ncf%var(i)%indx )
      stat(i) = nf90_def_var( ncf%iounit, trim( ncf%var(i)%name ), & 
                netcdf_precision, var_dim(:), ncf%var(i)%indx )
      if ( stat(i) /= NF90_NOERR ) then
        write(fstderr,*) "Error defining variable ",  & 
          ncf%var(i)%name //": ", trim( nf90_strerror( stat(i) ) )
        l_error = .true.
      end if

      stat(i) = nf90_put_att( ncf%iounit, ncf%var(i)%indx, & 
                "valid_range", var_range(1:2) )
      if ( stat(i) /= NF90_NOERR ) then
        write(fstderr,*) "Error defining valid range", & 
          trim( nf90_strerror( stat(i) ) )
        l_error = .true.
      end if

      stat(i) = nf90_put_att( ncf%iounit, ncf%var(i)%indx, "long_name",  & 
                trim( ncf%var(i)%description ) )
      if ( stat(i) /= NF90_NOERR ) then
        write(fstderr,*) "Error in description", & 
          trim( nf90_strerror( stat(i) ) )
        l_error = .true.
      end if

      stat(i) = nf90_put_att( ncf%iounit, ncf%var(i)%indx, "units",  & 
                trim( ncf%var(i)%units ) )
      if ( stat(i) /= NF90_NOERR ) then
        write(fstderr,*) "Error in units", & 
          trim( nf90_strerror( stat(i) ) )
        l_error = .true.
      end if
    end do

    if ( l_error ) stop "Error in netCDF file definition."

    deallocate( stat )

    if ( l_output_file_run_date ) then
      allocate( stat(3) )
    else
      allocate( stat(2) )
    end if

    ! Define global attributes of the file, for reproducing the results and
    ! determining how a run was configured
    stat(1) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "Conventions", "COARDS" )
    stat(2) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "model", "CLUBB" )

    if ( l_output_file_run_date ) then

      ! Enabling l_output_file_run_date allows the date and time that the
      ! netCDF output file is created to be included in the netCDF output file.
      ! Disabling l_output_file_run_date means that this information will not
      ! be included in the netCDF output file.  The advantage of disabling this
      ! output is that it allows for a check for binary differences between two
      ! netCDF output files.

      ! Figure out when the model is producing this file
      call date_and_time( current_date, current_time )

      stat(3) = nf90_put_att(ncf%iounit, NF90_GLOBAL, "created_on", &
                             current_date(1:4)//'-'//current_date(5:6)//'-'// &
                             current_date(7:8)//' '// &
                             current_time(1:2)//':'//current_time(3:4) )

    end if ! l_output_file_run_date

    if ( any( stat /= NF90_NOERR ) ) then
      write(fstderr,*) "Error writing model information"
      do i = 1, size( stat ), 1
        write(fstderr,*) trim( nf90_strerror( stat(i) ) )
      end do
      stop
    end if

    ! Write the model flags to the file
    deallocate( stat )
    allocate( stat(8) ) ! # of model flags

    stat(1) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_pos_def", lchar( l_pos_def ) )
    stat(2) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_hole_fill", lchar( l_hole_fill ) )
    stat(3) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_clip_semi_implicit", &
      lchar( l_clip_semi_implicit ) )
    stat(4) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_standard_term_ta", &
      lchar( l_standard_term_ta ) )
    stat(5) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_single_C2_Skw", &
      lchar( l_single_C2_Skw ) )
    stat(6) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_gamma_Skw", lchar( l_gamma_Skw ) )
    stat(7) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_uv_nudge", lchar( l_uv_nudge ) )
    stat(8) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_tke_aniso", lchar( l_tke_aniso ) )

    if ( any( stat /= NF90_NOERR ) ) then
      write(fstderr,*) "Error writing model flags"
      do i = 1, size( stat ), 1
        write(fstderr,*) i, trim( nf90_strerror( stat(i) ) )
      end do
      stop
    end if

    ! Write model parameter values to the file
    deallocate( stat )
    allocate( stat(nparams) )

    stat(1) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "T0", T0 )
    stat(2) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "ts_nudge", ts_nudge )
    stat(3) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "sclr_tol", sclr_tol )

    call get_parameters( params )

    do i = 1, nparams, 1
      stat(i) = nf90_put_att( ncf%iounit, NF90_GLOBAL, params_list(i), params(i) )
    end do

    if ( any( stat /= NF90_NOERR ) ) then
      write(fstderr,*) "Error writing parameters"
      do i = 1, nparams, 1
        write(fstderr,*) i, trim( nf90_strerror( stat(i) ) )
      end do
      stop
    end if

    stat(1) = nf90_enddef( ncf%iounit ) ! end definitions
    if ( stat(1) /= NF90_NOERR ) then
      write(fstderr,*) "Error finalizing definitions", & 
        trim( nf90_strerror( stat(1) ) )
      stop
    end if

    deallocate( stat )

    return
  end subroutine first_write

!-------------------------------------------------------------------------------
  subroutine write_grid( ncf )

! Description:
!   Writes inforation about latitude, longitude and the grid
! References:
!   None
!-------------------------------------------------------------------------------

    use netcdf, only: & 
        NF90_NOERR,   & ! Variable(s)
        nf90_put_var,  & ! Procedure(s)
        nf90_strerror
    use stat_file_module, only: & 
        stat_file ! Type
    use constants_clubb, only:  & 
        fstderr ! Variable

    implicit none

    ! Input Variable(s)
    type (stat_file), intent(inout) :: ncf

    integer :: stat

    ! ---- Begin Code ----

    stat = nf90_put_var( ncid=ncf%iounit, varid=ncf%AltVarId,  & 
                         values=ncf%z(ncf%ia:ncf%iz) )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error entering grid: ",  & 
        trim( nf90_strerror( stat ) )
      stop
    end if

    stat = nf90_put_var( ncid=ncf%iounit, varid=ncf%LongVarId,  & 
                         values=ncf%rlon )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error entering longitude: ",  & 
        trim( nf90_strerror( stat ) )
      stop
    end if

    stat = nf90_put_var( ncid=ncf%iounit, varid=ncf%LatVarId,  & 
                         values=ncf%rlat )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error entering latitude: ",  & 
        trim( nf90_strerror( stat ) )
      stop
    end if

    return
  end subroutine write_grid

!-------------------------------------------------------------------------------

  subroutine format_date & 
             ( day_in, month_in, year_in, time_in, date )

! Description:
!   Put the model date in a format that udunits and NetCDF can easily
!   handle.  GrADSnc is dumb and apparently cannot handle time
!   intervals < 1 minute.

! Notes:
!   Adapted from the original GrADS version written by Chris Golaz.
!   Uses Fortran `internal' files to write the string output.
!-------------------------------------------------------------------------------

    use calendar, only:  &
      compute_current_date ! Procedure(s)

    use clubb_precision, only:  & 
        time_precision ! Variable(s)

    implicit none

    ! External
    intrinsic :: floor, int, mod, nint

    ! Input Variables
    integer, intent(in) ::  & 
      day_in,           & ! Day of Month at Model Start   [dd]
      month_in,         & ! Month of Year at Model Start  [mm]
      year_in             ! Year at Model Start         [yyyy]

    real(kind=time_precision), intent(in) :: time_in ! Start time [s]

    ! Output Variables
    character(len=35), intent(out) :: date

    integer::  & 
      iday, imonth, iyear  ! Integer for day, month and year.

    real(kind=time_precision) :: st_time ! Start time [s]

    call compute_current_date( day_in, month_in,  & 
                               year_in, & 
                               time_in, & 
                               iday, imonth, & 
                               iyear, & 
                               st_time )

    if ( .not. l_grads_netcdf_boost_ts ) then
      date = "seconds since YYYY-MM-DD HH:MM:00.0"
    else
      date = "minutes since YYYY-MM-DD HH:MM:00.0"
    end if
    write(date(15:18),'(i4.4)') iyear
    write(date(20:21),'(i2.2)') imonth
    write(date(23:24),'(i2.2)') iday
    write(date(26:27),'(i2.2)') floor( st_time / 3600._time_precision )
    write(date(29:30),'(i2.2)') int( mod( nint( st_time ),3600 ) / 60 )

    return
  end subroutine format_date

!===============================================================================
  character function lchar( l_input )
! Description:
!   Cast a logical to a character data type.
!
! References:
!   None
!-------------------------------------------------------------------------------

    implicit none

    ! Input Variable
    logical, intent(in) :: l_input

    ! ---- Begin Code ----

    if ( l_input ) then
      lchar = 'T'
    else
      lchar = 'F'
    end if

    return
  end function lchar

#endif /*NETCDF*/
end module output_netcdf

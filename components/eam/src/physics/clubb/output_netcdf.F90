!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module output_netcdf
#ifdef NETCDF

! Description:
!   Functions and subroutines for writing NetCDF files

! References:
!   <http://www.unidata.ucar.edu/software/netcdf/docs/>
!-------------------------------------------------------------------------------

  implicit none

  public :: open_netcdf_for_writing, write_netcdf, close_netcdf

  private :: define_netcdf, write_grid, first_write, format_date

  ! Constant parameters
  ! This will truncate all timesteps smaller than 1 mn to a minute for 
  ! the purposes of viewing the data in grads
  logical, parameter, private :: &
    l_grads_netcdf_boost_ts = .false. 

  private ! Default scope

  contains
!-------------------------------------------------------------------------------
  subroutine open_netcdf_for_writing( nlat, nlon, fdir, fname, ia, iz, zgrid,  & 
                          day, month, year, lat_vals, lon_vals, & 
                          time, dtwrite, nvar, &
                          ncf, &
                          nsamp) ! optional

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
        fstderr ! Variable(s)

    use error_code, only: &
        err_code, &           ! Error Indicator
        clubb_fatal_error     ! Constant

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
      lat_vals ! Latitudes   [degrees_E]

    real( kind = core_rknd ), dimension(nlon), intent(in) ::  & 
      lon_vals ! Longitudes  [degrees_N]

    real( kind = core_rknd ), intent(in) :: & 
      dtwrite ! Time between write intervals   [s]

    real( kind = time_precision ), intent(in) ::  & 
     time   ! Current time                    [s]

    real( kind = core_rknd ), dimension(:), intent(in) ::  & 
      zgrid  ! The model grid                  [m]

    ! Input/output Variables
    type (stat_file), intent(inout) :: ncf

    ! Number of SILHS samples, used only for SILHS sample outputting
    integer, optional, intent(in) :: nsamp

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
    ncf%ntimes = 0

! According to Chris Vogl, netcdf can handle time steps < 1 min.  So this check is unneeded.
!    ! Check to make sure the timestep is appropriate. The GrADS program does not support an
!    ! output timestep less than 1 minute.  Other programs can read netCDF files like this
!    if ( dtwrite < sec_per_min ) then
!      write(fstderr,*) "Warning: GrADS program requires an output timestep of at least &
!                       &one minute, but the requested output timestep &
!                       &(stats_tout) is less than one minute."
!      if ( .not. l_allow_small_stats_tout ) then
!        write(fstderr,*) "To override this warning, set l_allow_small_stats_tout = &
!                         &.true. in the stats_setting namelist in the &
!                         &appropriate *_model.in file."
!        write(fstderr,*) "Fatal error in open_netcdf_for_writing"
!        err_code = clubb_fatal_error
!        return
!      end if
!    end if ! dtwrite < sec_per_min

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

    allocate( ncf%lat_vals(1:nlat), ncf%lon_vals(1:nlon) )

    ncf%lat_vals = lat_vals
    ncf%lon_vals = lon_vals

    ! If nsamp is present, SILHS samples are being handled.  Therefore set
    ! ncf%nsamp and ncf%samp_idx.  samp_idx holds the SILHS sample indices.
    if ( present(nsamp) ) then
      ncf%nsamp = nsamp
      allocate( ncf%samp_idx(1:nsamp) )
      forall( k=1:nsamp )
        ncf%samp_idx(k) = real( k, kind = core_rknd )
      end forall
    endif

    ! Create NetCDF dataset: enter define mode
    stat = nf90_create( path = trim( fdir )//trim( fname )//'.nc',  & 
                        cmode = NF90_CLOBBER,  & ! overwrite existing file
                        ncid = ncf%iounit )
    if ( stat /= NF90_NOERR ) then
      write(unit=fstderr,fmt=*) "Error opening file: ",  & 
        trim( fdir )//trim( fname )//'.nc', & 
        trim( nf90_strerror( stat ) )
      err_code = clubb_fatal_error
      return
    end if

    call define_netcdf( ncf%iounit, ncf%nlat, ncf%nlon, ncf%iz, ncf%nsamp, & ! In
                  ncf%day, ncf%month, ncf%year, ncf%time, & ! In
                  ncf%SampDimId, ncf%LatDimId, ncf%LongDimId, ncf%AltDimId, ncf%TimeDimId, &  ! Out
                  ncf%SampVarId, ncf%LatVarId, ncf%LongVarId, ncf%AltVarId, ncf%TimeVarId ) ! Out

    return
  end subroutine open_netcdf_for_writing

!-------------------------------------------------------------------------------

  subroutine write_netcdf( clubb_params, &
                           l_uv_nudge, &
                           l_tke_aniso, &
                           l_standard_term_ta, &
                           ncf )

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

    use error_code, only: &
        err_code, &           ! Error Indicator
        clubb_fatal_error     ! Constant

    use parameter_indices, only: &
        nparams    ! Variable(s)

    use clubb_precision, only: &
        time_precision, & ! Constant(s)
        core_rknd

    implicit none

    ! Input
    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    logical, intent(in) :: &
      l_uv_nudge,         & ! For wind speed nudging
      l_tke_aniso,        & ! For anisotropic turbulent kinetic energy, i.e. TKE = 1/2
                            ! (u'^2 + v'^2 + w'^2)
      l_standard_term_ta    ! Use the standard discretization for the turbulent advection terms.
                            ! Setting to .false. means that a_1 and a_3 are pulled outside of the
                            ! derivative in advance_wp2_wp3_module.F90 and in
                            ! advance_xp2_xpyp_module.F90.

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
      call first_write( clubb_params, & ! intent(in)
                        l_uv_nudge, & ! intent(in)
                        l_tke_aniso, & ! intent(in)
                        l_standard_term_ta, & ! intent(in)
                        ncf ) ! finalize the variable definitions intent(inout)
      call write_grid( ncf )  ! define lat., long., and grid intent(inout)
      ncf%l_defined = .true.
      if ( err_code == clubb_fatal_error ) return
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
      write(fstderr,*) "time variable nf90_put_var failed"
      err_code = clubb_fatal_error
      return
    end if

    ! If the grid_avg_var are allocated, then print to 4d netcdf.
    ! Otherwise, if the samples_of_var are allocated, print to 5d
    do i = 1, ncf%nvar, 1
      if ( allocated(ncf%grid_avg_var) ) then
         stat(i)  &
         = nf90_put_var( ncid=ncf%iounit, varid=ncf%grid_avg_var(i)%indx,  &
                         values=ncf%grid_avg_var(i)%ptr(:,:,ncf%ia:ncf%iz),  &
                         start=(/1,1,1,ncf%ntimes/), &
                         count=(/ncf%nlon,ncf%nlat,ncf%iz,1/) )
      elseif ( allocated(ncf%samples_of_var) ) then
        stat(i)  &
        = nf90_put_var( ncid=ncf%iounit, varid=ncf%samples_of_var(i)%indx,  &
                        values=ncf%samples_of_var(i)%ptr(:,:,:,ncf%ia:ncf%iz),  &
                        start=(/1,1,1,1,ncf%ntimes/), &
                        count=(/ncf%nsamp,ncf%nlon,ncf%nlat,ncf%iz,1/) )
      endif
    enddo ! i=1..nvar

    if ( any (stat /= NF90_NOERR ) ) then
      do i=1,ncf%nvar,1
        if( stat(i) /= NF90_NOERR ) then
          if ( allocated(ncf%grid_avg_var) ) then
            write(unit=fstderr,fmt=*) ncf%grid_avg_var(i)%name,  &
              trim( nf90_strerror( stat(i) ) )
          elseif ( allocated(ncf%samples_of_var) ) then
            write(unit=fstderr,fmt=*) ncf%samples_of_var(i)%name,  &
              trim( nf90_strerror( stat(i) ) )
          endif
        end if
      end do
      write(fstderr,*) "nf90_put_var error"
      err_code = clubb_fatal_error
      return
    end if

    deallocate( stat )

    return
  end subroutine write_netcdf

!-------------------------------------------------------------------------------
  subroutine define_netcdf( ncid, nlat, nlon, iz, nsamp, &
                            day, month, year, time, & 
                            SampDimId, LatDimId, LongDimId, AltDimId, TimeDimId, &
                            SampVarId, LatVarId, LongVarId, AltVarId, TimeVarId )

! Description:
!   Used internally to create a definition for the NetCDF dataset
!
! References:
!   None
!-------------------------------------------------------------------------------
    use netcdf, only: & 
        NF90_NOERR,   & ! Constants
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

    use error_code, only: &
        err_code, &           ! Error Indicator
        clubb_fatal_error     ! Constant

    implicit none

    integer, intent(in) ::  & 
      nlat,   & ! Number of points in the N/S direction
      nlon,   & ! Number of points in the E/W direction
      nsamp     ! Number of SILHS samples

    ! Input Variables
    integer, intent(in) ::  & 
      day, month, year,  & ! Time of year
      ncid,              & ! Number used by NetCDF for ref. the file
      iz                   ! Dimension in z

    real(kind=time_precision), intent(in) ::  & 
      time    ! Current model time [s]

    ! Output Variables
    integer, intent(out) ::  &
    ! NetCDF id's for dimensions, including for SILHS samples if needed
      SampDimId, LatDimId, LongDimId, AltDimId, TimeDimId

    ! NetCDF id's for data (e.g. longitude) associated with each dimension,
    ! including for SILHS samples if needed
    integer, intent(out) ::  & 
      SampVarId, LatVarId, LongVarId, AltVarId, TimeVarId

    ! Local variables
    integer :: stat
    character(len=35) :: TimeUnits

    ! ---- Begin Code ----

    ! Define the dimensions for the variables.
    ! Start with SILHS samples so this dimension is listed first in the netCDF
    ! file.  Since ncf not present to test allocation, test using nsamp. nsamp
    ! is initialized to zero, will only be nonzero if printing SILHS samples.
    if ( nsamp > 0 ) then
      ! Define SILHS sample dimension
      stat =  nf90_def_dim( ncid, "lh_sample_number", nsamp, SampDimId )
      if ( stat /= NF90_NOERR ) then
        write(fstderr,*) "Error defining lh_sample_number: ", &
          trim( nf90_strerror( stat ) )
        err_code = clubb_fatal_error
        return
      end if
      ! Define SILHS sample number variable
      stat = nf90_def_var( ncid, "lh_sample_number", NF90_DOUBLE, &
                          (/SampDimId/), SampVarId )
      ! Attributes for SILHS sample number variable
      stat = nf90_put_att( ncid, SampVarId, "description", "SILHS sample (i.e. subcolumn) index" )
      stat = nf90_put_att( ncid, SampVarId, "units", "number" )
    endif !if nsamp>0

    stat = nf90_def_dim( ncid, "longitude", nlon, LongDimId )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining longitude: ", & 
        trim( nf90_strerror( stat ) )
      err_code = clubb_fatal_error
      return
    end if

    stat =  nf90_def_dim( ncid, "latitude", nlat, LatDimId )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining latitude: ", & 
        trim( nf90_strerror( stat ) )
      err_code = clubb_fatal_error
      return
    end if

    stat = nf90_def_dim( ncid, "altitude", iz, AltDimId )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining altitude: ", & 
      trim( nf90_strerror( stat ) )
      err_code = clubb_fatal_error
      return
    end if

    stat =  nf90_def_dim( ncid, "time", NF90_UNLIMITED, TimeDimId )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining time: ", & 
        trim( nf90_strerror( stat ) )
      err_code = clubb_fatal_error
      return
    end if

    ! Define the initial variables for the dimensions
    ! Longitude = deg_E = X
    stat = nf90_def_var( ncid, "longitude", NF90_DOUBLE, & 
                         (/LongDimId/), LongVarId )

    ! Latitude = deg_N = Y
    stat = nf90_def_var( ncid, "latitude", NF90_DOUBLE, & 
                         (/LatDimId/), LatVarId )

    ! Altitude = meters above the surface = Z
    stat = nf90_def_var( ncid, "altitude", NF90_DOUBLE, & 
                        (/AltDimId/), AltVarId )

    ! grads2nc stores time as a double prec. value, so we follow that
    stat = nf90_def_var( ncid, "time", NF90_DOUBLE, & 
                         (/TimeDimId/), TimeVarId )

    ! Assign attribute values

    ! Time attribute
    stat = nf90_put_att( ncid, TimeVarId, "cartesian_axis", "T" )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining time: ", trim( nf90_strerror( stat ) )
      err_code = clubb_fatal_error
      return
    end if

    call format_date( day, month, year, time, & ! intent(in)
                      TimeUnits ) ! intent(out)

    stat = nf90_put_att( ncid, TimeVarId, "units", TimeUnits )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining time: ", trim( nf90_strerror( stat ) )
      err_code = clubb_fatal_error
      return
    end if

    stat = nf90_put_att( ncid, TimeVarId, "ipositive", 1 )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining time: ", trim( nf90_strerror( stat ) )
      err_code = clubb_fatal_error
      return
    end if

    stat = nf90_put_att( ncid, TimeVarId, "calendar_type", "Gregorian" )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error defining time", trim( nf90_strerror( stat ) )
      err_code = clubb_fatal_error
      return
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
      error stop "Fatal error"
    end if

    return
  end subroutine close_netcdf

!-------------------------------------------------------------------------------
  subroutine first_write( clubb_params, &
                          l_uv_nudge, &
                          l_tke_aniso, &
                          l_standard_term_ta, &
                          ncf )

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

    use parameter_indices, only: &
        nparams ! Integer

    use model_flags, only: &
        l_pos_def, &
        l_hole_fill, &
        l_gamma_Skw

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use error_code, only: &
        err_code, &           ! Error Indicator
        clubb_fatal_error     ! Constant

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

    ! Input Variables
    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    logical, intent(in) :: &
      l_uv_nudge,         & ! For wind speed nudging
      l_tke_aniso,        & ! For anisotropic turbulent kinetic energy, i.e. TKE = 1/2
                            ! (u'^2 + v'^2 + w'^2)
      l_standard_term_ta    ! Use the standard discretization for the turbulent advection terms.
                            ! Setting to .false. means that a_1 and a_3 are pulled outside of the
                            ! derivative in advance_wp2_wp3_module.F90 and in
                            ! advance_xp2_xpyp_module.F90.

    ! Input/Output Variables
    type (stat_file), intent(inout) :: ncf

    ! Local Variables
    integer, dimension(:), allocatable :: stat
    
    integer :: netcdf_precision ! Level of precision for netCDF output

    integer :: i     ! Array index

    character(len=10) :: current_time
    character(len=8)  :: current_date
    ! Range for NetCDF variables
    real( kind = core_rknd ), dimension(2) :: var_range

    ! Dimensions for variables
    integer, allocatable, dimension(:) :: var_dim

!-------------------------------------------------------------------------------
!      Typical valid ranges (IEEE 754)

!      real(kind=4): +/- 3.4028235E+38
!      real(kind=8): +/- 1.797693134862316E+308
!      real(kind=16):+/- 1.189731495357231765085759326628007E+4932

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

    ! If samples_of_var is allocated, print to 5d netcdf, otherwise 4d.
    if ( allocated(ncf%samples_of_var) ) then
      allocate( var_dim(1:5) )
      var_dim(1) = ncf%SampDimId
      i = 1
    else
      allocate( var_dim(1:4) )
      i = 0
    endif

    var_dim(i+1) = ncf%LongDimId ! X
    var_dim(i+2) = ncf%LatDimId  ! Y
    var_dim(i+3) = ncf%AltDimId  ! Z
    var_dim(i+4) = ncf%TimeDimId ! The NF90_UNLIMITED dimension

    allocate( stat( ncf%nvar ) )

    select case (core_rknd)
      case ( selected_real_kind( p=5 ) )
        netcdf_precision = NF90_FLOAT
      case ( selected_real_kind( p=12 ) )
        netcdf_precision = NF90_DOUBLE
      case default
        netcdf_precision = NF90_DOUBLE
    end select

    ! Specify whether "grid_avg_var" or "samples_of_var"
    do i = 1, ncf%nvar, 1
      if ( allocated(ncf%grid_avg_var) ) then
        stat(i) = nf90_def_var( ncf%iounit, trim( ncf%grid_avg_var(i)%name ), &
                    netcdf_precision, var_dim(:), ncf%grid_avg_var(i)%indx )
      elseif ( allocated(ncf%samples_of_var) ) then
        stat(i) = nf90_def_var( ncf%iounit, trim( ncf%samples_of_var(i)%name ), &
                    netcdf_precision, var_dim(:), ncf%samples_of_var(i)%indx )
      endif
      if ( stat(i) /= NF90_NOERR ) then
        write(fstderr,*) "Error defining variable ",  & 
          ncf%grid_avg_var(i)%name //": ", trim( nf90_strerror( stat(i) ) )
        err_code = clubb_fatal_error
        return
      end if

      if ( allocated(ncf%grid_avg_var) ) then
        stat(i) = nf90_put_att( ncf%iounit, ncf%grid_avg_var(i)%indx, &
                    "valid_range", var_range(1:2) )
      elseif ( allocated(ncf%samples_of_var) ) then
        stat(i) = nf90_put_att( ncf%iounit, ncf%samples_of_var(i)%indx, &
                    "valid_range", var_range(1:2) )
      endif
      if ( stat(i) /= NF90_NOERR ) then
        write(fstderr,*) "Error defining valid range", & 
          trim( nf90_strerror( stat(i) ) )
        err_code = clubb_fatal_error
        return
      end if

      if ( allocated(ncf%grid_avg_var) ) then
        stat(i) = nf90_put_att( ncf%iounit, ncf%grid_avg_var(i)%indx, "long_name",  &
                  trim( ncf%grid_avg_var(i)%description ) )
      elseif ( allocated(ncf%samples_of_var) ) then
        stat(i) = nf90_put_att( ncf%iounit, ncf%samples_of_var(i)%indx, "long_name",  &
                  trim( ncf%samples_of_var(i)%description ) )
      endif
      if ( stat(i) /= NF90_NOERR ) then
        write(fstderr,*) "Error in description", & 
          trim( nf90_strerror( stat(i) ) )
        err_code = clubb_fatal_error
        return
      end if

      if ( allocated(ncf%grid_avg_var) ) then
        stat(i) = nf90_put_att( ncf%iounit, ncf%grid_avg_var(i)%indx, "units",  &
                  trim( ncf%grid_avg_var(i)%units ) )
      elseif ( allocated(ncf%samples_of_var) ) then
        stat(i) = nf90_put_att( ncf%iounit, ncf%samples_of_var(i)%indx, "units",  &
                  trim( ncf%samples_of_var(i)%units ) )
      endif
      if ( stat(i) /= NF90_NOERR ) then
        write(fstderr,*) "Error in units", & 
          trim( nf90_strerror( stat(i) ) )
        err_code = clubb_fatal_error
        return
      end if
    end do

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
      err_code = clubb_fatal_error
      return
    end if

    ! Write the model flags to the file
    deallocate( stat )
    allocate( stat(6) ) ! # of model flags

    stat(1) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_pos_def", lchar( l_pos_def ) )
    stat(2) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_hole_fill", lchar( l_hole_fill ) )
    stat(3) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_standard_term_ta", &
      lchar( l_standard_term_ta ) )
    stat(4) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_gamma_Skw", lchar( l_gamma_Skw ) )
    stat(5) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_uv_nudge", lchar( l_uv_nudge ) )
    stat(6) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "l_tke_aniso", lchar( l_tke_aniso ) )

    if ( any( stat /= NF90_NOERR ) ) then
      write(fstderr,*) "Error writing model flags"
      do i = 1, size( stat ), 1
        write(fstderr,*) i, trim( nf90_strerror( stat(i) ) )
      end do
      err_code = clubb_fatal_error
      return
    end if

    ! Write model parameter values to the file
    deallocate( stat )
    allocate( stat(nparams) )

    stat(1) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "T0", T0 )
    stat(2) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "ts_nudge", ts_nudge )
    stat(3) = nf90_put_att( ncf%iounit, NF90_GLOBAL, "sclr_tol", sclr_tol )

    do i = 1, nparams, 1
      stat(i) = nf90_put_att( ncf%iounit, NF90_GLOBAL, params_list(i), clubb_params(i) )
    end do

    if ( any( stat /= NF90_NOERR ) ) then
      write(fstderr,*) "Error writing parameters"
      do i = 1, nparams, 1
        write(fstderr,*) i, trim( nf90_strerror( stat(i) ) )
      end do
      err_code = clubb_fatal_error
      return
    end if

    stat(1) = nf90_enddef( ncf%iounit ) ! end definitions
    if ( stat(1) /= NF90_NOERR ) then
      write(fstderr,*) "Error finalizing definitions", & 
        trim( nf90_strerror( stat(1) ) )
      err_code = clubb_fatal_error
      return
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

    use error_code, only: &
        err_code, &         ! Error Indicator
        clubb_fatal_error   ! Constant

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
      err_code = clubb_fatal_error
      return
    end if

    stat = nf90_put_var( ncid=ncf%iounit, varid=ncf%LongVarId,  & 
                         values=ncf%lon_vals )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error entering longitude: ",  & 
        trim( nf90_strerror( stat ) )
      err_code = clubb_fatal_error
      return
    end if

    stat = nf90_put_var( ncid=ncf%iounit, varid=ncf%LatVarId,  & 
                         values=ncf%lat_vals )
    if ( stat /= NF90_NOERR ) then
      write(fstderr,*) "Error entering latitude: ",  & 
        trim( nf90_strerror( stat ) )
      err_code = clubb_fatal_error
      return
    end if

    ! Write the SILHS sample indices if samples_of_var allocated
    if ( allocated(ncf%samples_of_var) ) then
      stat = nf90_put_var( ncid=ncf%iounit, varid=ncf%SampVarId,  &
                           values=ncf%samp_idx )
      if ( stat /= NF90_NOERR ) then
        write(fstderr,*) "Error entering grid: ",  &
          trim( nf90_strerror( stat ) )
        err_code = clubb_fatal_error
        return
      end if
    endif

    return
  end subroutine write_grid

!-------------------------------------------------------------------------------

  subroutine format_date & 
             ( day_in, month_in, year_in, time_in, &
               date )

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

    call compute_current_date( day_in, month_in,  & ! intent(in)
                               year_in, &  ! intent(in)
                               time_in, & ! intent(in)
                               iday, imonth, & ! intent(out)
                               iyear, &  ! intent(out)
                               st_time ) ! intent(out)

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
    
    if ( .not. l_grads_netcdf_boost_ts ) then
      write(date(32:33),'(i2.2)') nint(((real(mod( nint( st_time ),3600),kind=time_precision) / &
                     60._time_precision) - (real(int(mod( nint( st_time ),3600 ) / 60 ), & 
                                               kind=time_precision) ) )*60._time_precision)
    end if

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

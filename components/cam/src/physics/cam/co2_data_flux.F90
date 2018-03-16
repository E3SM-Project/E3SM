
module co2_data_flux

!------------------------------------------------------------------------------------------------
! for data reading and interpolation                                           
!------------------------------------------------------------------------------------------------

  use shr_kind_mod,   only : r8 => shr_kind_r8
  use spmd_utils,     only : masterproc
  use ppgrid,         only : begchunk, endchunk, pcols
  use phys_grid,      only : scatter_field_to_chunk, get_ncols_p, get_rlat_all_p, get_rlon_all_p
  use error_messages, only : alloc_err, handle_ncerr, handle_err
  use cam_abortutils,     only : endrun
  use netcdf
  use error_messages, only : handle_ncerr  
  use cam_logfile,    only : iulog
  use interpolate_data,  only : lininterp_init, lininterp, interp_type, lininterp_finish
  use mo_constants, only : pi, d2r
#if ( defined SPMD )
  use mpishorthand, only: mpicom, mpiint, mpir8
#endif

  implicit none

! public data

! public type

  public read_interp

! public interface

  public read_data_flux
  public interp_time_flux

! private data

  private
!  integer,  parameter :: totsz=2000           ! number greater than data time sample
  real(r8), parameter :: zero=0_r8, twopi=2_r8*pi !parameters required for interpolation
  real(r8), parameter :: daysperyear = 365.0_r8  ! Number of days in a year         
 
!--------------------------------------------------------------------------------------------------
TYPE :: read_interp          
      
  real(r8), pointer, dimension(:,:)   :: co2flx
                       ! Interpolated output (pcols,begchunk:endchunk)
  real(r8), pointer, dimension(:,:,:) :: co2bdy
                       ! bracketing data     (pcols,begchunk:endchunk,2)
  real(r8) :: cdayfm   ! Calendar day for prv. month read in
  real(r8) :: cdayfp   ! Calendar day for nxt. month read in
  integer :: nm_f      ! Array indices for prv. month data
  integer :: np_f      ! Array indices for nxt. month data
  integer :: np1_f     ! current forward time index of dataset
  integer :: timesz    ! size of time dimension on dataset
  integer :: lonsiz,latsiz           !size of lats and lons
  integer, pointer :: date_f(:)      ! Date on dataset (YYYYMMDD)
  integer, pointer :: sec_f(:)       ! seconds of date on dataset (0-86399) 
  real(r8), pointer :: lon_f(:)      !lon values
  real(r8), pointer :: lat_f(:)      !lat values
  character(len=256) :: locfn   ! dataset name
  integer :: ncid_f             ! netcdf id for dataset       
  integer :: fluxid             ! netcdf id for dataset flux      
  real(r8), pointer :: xvar(:,:,:) ! work space for dataset  
 
END TYPE read_interp
!-------------------------------------------------------------------------------------------------
                       
contains
 
!===============================================================================

subroutine read_data_flux (input_file, xin)
 
!-------------------------------------------------------------------------------              
! Do initial read of time-varying 2d(lat,lon NetCDF dataset, 
! reading two data bracketing current timestep
!-------------------------------------------------------------------------------
 
  use time_manager, only : get_curr_date, get_curr_calday, &
                           is_perpetual, get_perp_date, get_step_size, is_first_step
  use ioFileMod,    only : getfil

  implicit none
 
!---------------------------Common blocks-------------------------------
! Dummy arguments
  character(len=*),   intent(in)    :: input_file
  TYPE(read_interp),  intent(inout) :: xin   
 
  type(interp_type) :: lon_wgts, lat_wgts !weights -output from horizontal interpolation
  integer lonid                 ! netcdf id for longitude variable
  integer latid                 ! netcdf id for latitude variable
  integer timeid                ! netcdf id for time variable
  integer dateid                ! netcdf id for date variable
  integer secid                 ! netcdf id for seconds variable
  integer londimid              ! netcdf id for longitude variable
  integer latdimid              ! netcdf id for latitude variable
 
  integer dtime                 ! timestep size [seconds]
  integer cnt3(3)               ! array of counts for each dimension
  integer strt3(3)              ! array of starting indices
  integer n                     ! indices
  integer j                     ! latitude index
  integer istat                 ! error return
  integer  :: yr, mon, day      ! components of a date
  integer  :: ncdate            ! current date in integer format [yyyymmdd]
  integer  :: ncsec             ! current time of day [seconds]
  integer :: ncol, lchnk
  real(r8) ::  to_lats(pcols), to_lons(pcols)
  real(r8) calday               ! calendar day (includes yr if no cycling)
  real(r8) caldayloc            ! calendar day (includes yr if no cycling)
 
  xin%nm_f = 1
  xin%np_f = 2
 
! Allocate space for data.
 
  allocate( xin%co2flx(pcols,begchunk:endchunk), stat=istat )
  call alloc_err( istat, 'CO2FLUX_READ', 'co2flx', &
       pcols*(endchunk-begchunk+1) )
 
  allocate( xin%co2bdy(pcols,begchunk:endchunk,2), stat=istat )
  call alloc_err( istat, 'CO2FLUX_READ', 'co2bdy', &
       pcols*(endchunk-begchunk+1)*2 )
 
! SPMD: Master does all the work.
 
  if (masterproc) then
!
! Use year information only if not cycling sst dataset
!
     if (is_first_step()) then
        dtime = get_step_size()
        dtime = -dtime
        calday = get_curr_calday(offset=dtime)
     else
        calday = get_curr_calday()
     endif
     if ( is_perpetual() ) then
        call get_perp_date(yr, mon, day, ncsec)
     else
        if (is_first_step()) then
           call get_curr_date(yr, mon, day, ncsec,offset=dtime)
        else
           call get_curr_date(yr, mon, day, ncsec)
        endif
     end if
 
     ncdate = yr*10000 + mon*100 + day
 
     caldayloc = calday + yr*daysperyear
 
! Open NetCDF File
 
     call getfil(input_file, xin%locfn)
     call handle_ncerr( nf90_open(xin%locfn, 0, xin%ncid_f),&
       'co2_data_flux.F90:154')
     write(iulog,*)'CO2FLUX_READ: NCOPN returns id ',xin%ncid_f,' for file ',trim(xin%locfn)
 
! Get and check dimension info
 
     call handle_ncerr( nf90_inq_dimid( xin%ncid_f, 'lon', londimid ),&
       'co2_data_flux.F90:160')
     call handle_ncerr( nf90_inq_dimid( xin%ncid_f, 'lat', latdimid ),&
       'co2_data_flux.F90:162')
     call handle_ncerr( nf90_inq_dimid( xin%ncid_f, 'time',  timeid ),&
       'co2_data_flux.F90:164')

     call handle_ncerr( nf90_inquire_dimension( xin%ncid_f, londimid, len=xin%lonsiz     ),&
       'co2_data_flux.F90:167')
     call handle_ncerr( nf90_inquire_dimension( xin%ncid_f, latdimid, len=xin%latsiz     ),&
       'co2_data_flux.F90:169')
     call handle_ncerr( nf90_inquire_dimension( xin%ncid_f, timeid,   len=xin%timesz ),&
       'co2_data_flux.F90:171')

     allocate(xin%date_f(xin%timesz), xin%sec_f(xin%timesz))
! Get data id
        
     call handle_ncerr( nf90_inq_varid( xin%ncid_f, 'date',     dateid     ),&
       'co2_data_flux.F90:192')
     call handle_ncerr( nf90_inq_varid( xin%ncid_f, 'datesec',  secid      ),&
       'co2_data_flux.F90:194')
     call handle_ncerr( nf90_inq_varid( xin%ncid_f, 'CO2_flux', xin%fluxid ),&
       'co2_data_flux.F90:196')
 
! Retrieve entire date and sec variables.
 
     call handle_ncerr( nf90_get_var ( xin%ncid_f, dateid, xin%date_f ),&
       'co2_data_flux.F90:201')
     call handle_ncerr( nf90_get_var ( xin%ncid_f, secid,  xin%sec_f  ),&
       'co2_data_flux.F90:203')

!Retrieve lats and lons for horizontal interpolatio

     allocate(xin%lon_f(xin%lonsiz), xin%lat_f(xin%latsiz)) !handle allocation error!!

     call handle_ncerr( nf90_inq_varid( xin%ncid_f, 'lon',     lonid     ),&
          'co2_data_flux.F90:192')
     call handle_ncerr( nf90_inq_varid( xin%ncid_f, 'lat',     latid     ),&
          'co2_data_flux.F90:192')

     call handle_ncerr( nf90_get_var ( xin%ncid_f, lonid,  xin%lon_f  ),&
          'co2_data_flux.F90:200_2')

     call handle_ncerr( nf90_get_var ( xin%ncid_f, latid,  xin%lat_f  ),&
          'co2_data_flux.F90:200_2')

     !convert to radians
     xin%lat_f = xin%lat_f * d2r !degree to radians
     xin%lon_f = xin%lon_f * d2r !degree to radians

! initialize
 
     strt3(1) = 1
     strt3(2) = 1
     strt3(3) = 1
     cnt3(1)  = xin%lonsiz
     cnt3(2)  = xin%latsiz
     cnt3(3)  = 1

  endif

#ifdef SPMD
  call mpibcast( cnt3, 2,     mpiint, 0, mpicom )
#endif
  allocate(xin%xvar(cnt3(1),cnt3(2),2))
  if (masterproc) then
! Normal interpolation between consecutive time slices.

     do n=1,xin%timesz-1

        xin%np1_f = n + 1
 
        call bnddyi(xin%date_f(n  ),       xin%sec_f(n  ),       xin%cdayfm)
        call bnddyi(xin%date_f(xin%np1_f), xin%sec_f(xin%np1_f), xin%cdayfp)
 
        yr         = xin%date_f(n)/10000
        xin%cdayfm = xin%cdayfm + yr*daysperyear

        yr         = xin%date_f(xin%np1_f)/10000
        xin%cdayfp = xin%cdayfp + yr*daysperyear
 
!       read 2 time sample bracketing ncdate

        if ( caldayloc > xin%cdayfm .and. caldayloc <= xin%cdayfp ) then
 
           strt3(3) = n
           call handle_ncerr( nf90_get_var ( xin%ncid_f, xin%fluxid, xin%xvar(:,:,xin%nm_f), strt3, cnt3), &
             'co2_data_flux.F90:235')
           strt3(3) = xin%np1_f                                      
           call handle_ncerr( nf90_get_var ( xin%ncid_f, xin%fluxid, xin%xvar(:,:,xin%np_f), strt3, cnt3),&
             'co2_data_flux.F90:238')

           goto 10

        end if
 
     end do

     write(iulog,*)'CO2FLUX_READ: Failed to find dates bracketing ncdate, ncsec=', ncdate, ncsec
     call endrun
 
10   continue
     write(iulog,*)'CO2FLUX_READ: Read ', trim(xin%locfn), ' for dates ', xin%date_f(n), xin%sec_f(n), &
          ' and ', xin%date_f(xin%np1_f), xin%sec_f(xin%np1_f)
 
#if (defined SPMD )
     call mpibcast( xin%timesz, 1,     mpiint, 0, mpicom )
     call mpibcast( xin%lonsiz, 1,     mpiint, 0, mpicom )
     call mpibcast( xin%latsiz, 1,     mpiint, 0, mpicom )
     call mpibcast( xin%date_f, xin%timesz, mpiint, 0, mpicom )
     call mpibcast( xin%sec_f,  xin%timesz, mpiint, 0, mpicom )
     call mpibcast( xin%lon_f,  xin%lonsiz, mpir8, 0, mpicom )
     call mpibcast( xin%lat_f,  xin%latsiz, mpir8, 0, mpicom )
     call mpibcast( xin%cdayfm, 1,     mpir8 , 0, mpicom )
     call mpibcast( xin%cdayfp, 1,     mpir8,  0, mpicom )
     call mpibcast( xin%np1_f,  1,     mpiint, 0, mpicom )
  else
     call mpibcast( xin%timesz, 1,     mpiint, 0, mpicom )
     call mpibcast( xin%lonsiz, 1,     mpiint, 0, mpicom )
     call mpibcast( xin%latsiz, 1,     mpiint, 0, mpicom )
     allocate(xin%date_f(xin%timesz), xin%sec_f(xin%timesz))
     call mpibcast( xin%date_f, xin%timesz, mpiint, 0, mpicom )
     call mpibcast( xin%sec_f,  xin%timesz, mpiint, 0, mpicom )
     allocate(xin%lon_f(xin%lonsiz), xin%lat_f(xin%latsiz))
     call mpibcast( xin%lon_f,  xin%lonsiz, mpir8, 0, mpicom )
     call mpibcast( xin%lat_f,  xin%latsiz, mpir8, 0, mpicom )
     call mpibcast( xin%cdayfm, 1,     mpir8 , 0, mpicom )
     call mpibcast( xin%cdayfp, 1,     mpir8,  0, mpicom )
     call mpibcast( xin%np1_f,  1,     mpiint, 0, mpicom )
#endif
  end if

  !comment out the following "scatter_field_to_chunk" as we need to do horizontal interpolation first
  !Hosrizontal interpolation will take care of distributing the data to chunks and columns
  !call scatter_field_to_chunk ( 1,1,2,cnt3(1), xin%xvar, xin%co2bdy )

  !Following code is taken from tracer_data.F90
  do lchnk=begchunk,endchunk
     ncol = get_ncols_p(lchnk)
     call get_rlat_all_p(lchnk, pcols, to_lats)
     call get_rlon_all_p(lchnk, pcols, to_lons)

     call lininterp_init(xin%lon_f,xin%lonsiz, to_lons, ncol, 2, lon_wgts, zero, twopi) !2 is for cyclic with limits 0 to 2pi
     call lininterp_init(xin%lat_f,xin%latsiz, to_lats, ncol, 1, lat_wgts)
     
     call lininterp(xin%xvar(:,:,xin%nm_f), xin%lonsiz, xin%latsiz, xin%co2bdy(1:ncol,lchnk,xin%nm_f), ncol, lon_wgts, lat_wgts)
     call lininterp(xin%xvar(:,:,xin%np_f), xin%lonsiz, xin%latsiz, xin%co2bdy(1:ncol,lchnk,xin%np_f), ncol, lon_wgts, lat_wgts)
     
     call lininterp_finish(lon_wgts)
     call lininterp_finish(lat_wgts)
  end do

  return
end subroutine read_data_flux
 
!===============================================================================

subroutine interp_time_flux (xin, prev_timestep)
 
!-----------------------------------------------------------------------
! Time interpolate data to current time.
! Reading in new monthly data if necessary.
!
!-----------------------------------------------------------------------
 
  use time_manager, only : get_curr_date, get_curr_calday, &
                          is_perpetual, get_perp_date, get_step_size, is_first_step
  use interpolate_data,   only : get_timeinterp_factors
 
  logical, intent(in), optional     :: prev_timestep ! If using previous timestep, set to true
  TYPE(read_interp),  intent(inout) :: xin

!---------------------------Local variables-----------------------------
  type(interp_type) :: lon_wgts, lat_wgts     
  integer dtime          ! timestep size [seconds]
  integer cnt3(3)        ! array of counts for each dimension
  integer strt3(3)       ! array of starting indices
  integer i,j,lchnk      ! indices
  integer ncol           ! number of columns in current chunk
  integer ntmp           ! temporary
  real(r8) fact1, fact2  ! time interpolation factors
  integer :: yr, mon, day! components of a date
  integer :: ncdate      ! current date in integer format [yyyymmdd]
  integer :: ncsec       ! current time of day [seconds]
  real(r8) :: calday     ! current calendar day
  real(r8) caldayloc     ! calendar day (includes yr if no cycling)
  real(r8) deltat        ! time (days) between interpolating data
  logical :: previous
  logical :: co2cyc=.false.

  real(r8) ::  to_lats(pcols), to_lons(pcols)


 
!-----------------------------------------------------------------------
 
! SPMD: Master does all the work.  Sends needed info to slaves
 
! Use year information only if a multiyear dataset
 
     if ( .not. present(prev_timestep) ) then
        previous = .false.
     else
        previous = prev_timestep
     end if
 
     if (previous .and. is_first_step()) then
        dtime = get_step_size()
        dtime = -dtime
        calday = get_curr_calday(offset=dtime)
     else
        calday = get_curr_calday()
     endif
 
     if ( is_perpetual() ) then
        call get_perp_date(yr, mon, day, ncsec)
     else
        if (previous .and. is_first_step()) then
           call get_curr_date(yr, mon, day, ncsec,offset=dtime)
        else
           call get_curr_date(yr, mon, day, ncsec)
        endif
     end if
 
     ncdate = yr*10000 + mon*100 + day
 
     caldayloc = calday + yr*daysperyear

     if (masterproc) then

        strt3(1) = 1
        strt3(2) = 1
        strt3(3) = 1
        cnt3(1)  = xin%lonsiz
        cnt3(2)  = xin%latsiz
        cnt3(3)  = 1

     endif

#ifdef SPMD
     call mpibcast(cnt3, 2, mpiint, 0, mpicom)
#endif
 
! If model time is past current forward data timeslice, read in the next
! timeslice for time interpolation. 
 
     if ( caldayloc > xin%cdayfp .and. .not. (xin%np1_f==1 .and. caldayloc > xin%cdayfm) ) then
 
        xin%np1_f = xin%np1_f + 1
 
        if ( xin%np1_f > xin%timesz ) then
           call endrun ('CO2FLUX_INTERP: Attempt to read past end of dataset')
        end if
 
        xin%cdayfm = xin%cdayfp
 
        call bnddyi( xin%date_f(xin%np1_f), xin%sec_f(xin%np1_f), xin%cdayfp )

        yr = xin%date_f(xin%np1_f)/10000
        xin%cdayfp = xin%cdayfp + yr*daysperyear

        if ( .not. (xin%np1_f == 1 .or. caldayloc <= xin%cdayfp) ) then
         
           if (masterproc) then
              write(iulog,*)'CO2FLUX_INTERP: Input data for date', xin%date_f(xin%np1_f), ' sec ', xin%sec_f(xin%np1_f), &
                        ' does not exceed model date', ncdate, ' sec ', ncsec, ' Stopping.'
           end if
           call endrun ()
        end if

        ntmp     = xin%nm_f
        xin%nm_f = xin%np_f
        xin%np_f = ntmp
 
        if (masterproc) then
           strt3(3) = xin%np1_f
           call handle_ncerr( nf90_get_var ( xin%ncid_f, xin%fluxid, xin%xvar(:,:,xin%np_f), strt3, cnt3 ),&
             'co2_data_flux.F90:391')
           write(iulog,*)'CO2FLUX_INTERP: Read ', trim(xin%locfn),' for date (yyyymmdd) ', xin%date_f(xin%np1_f), &
                     ' sec ', xin%sec_f(xin%np1_f)
        endif
 
       !comment out the following "scatter_field_to_chunk" as we need to do horizontal interpolation first
       !Hosrizontal interpolation will take care of distributing the data to chunks and columns 
       !call scatter_field_to_chunk ( 1,1,2,cnt3(1), xin%xvar, xin%co2bdy )

       !Following code is taken from tracer_data.F90
        do lchnk=begchunk,endchunk
           ncol = get_ncols_p(lchnk)
           call get_rlat_all_p(lchnk, pcols, to_lats)
           call get_rlon_all_p(lchnk, pcols, to_lons)

           call lininterp_init(xin%lon_f,xin%lonsiz, to_lons, ncol, 2, lon_wgts, zero, twopi) !2 is for cyclic with limits 0 to 2pi
           call lininterp_init(xin%lat_f,xin%latsiz, to_lats, ncol, 1, lat_wgts)
           
           call lininterp(xin%xvar(:,:,xin%nm_f), xin%lonsiz, xin%latsiz, xin%co2bdy(1:ncol,lchnk,xin%nm_f), ncol, lon_wgts, lat_wgts)
           call lininterp(xin%xvar(:,:,xin%np_f), xin%lonsiz, xin%latsiz, xin%co2bdy(1:ncol,lchnk,xin%np_f), ncol, lon_wgts, lat_wgts)
           
           call lininterp_finish(lon_wgts)
           call lininterp_finish(lat_wgts)
        end do
     end if
 
! Determine time interpolation factors.
 
     call get_timeinterp_factors ( co2cyc, xin%np1_f, xin%cdayfm, xin%cdayfp, caldayloc, fact1, fact2, 'CO2FLUX_INTERP:' )

     do lchnk=begchunk,endchunk
        ncol = get_ncols_p(lchnk)
        do i=1,ncol
           xin%co2flx(i,lchnk) = xin%co2bdy(i,lchnk,xin%nm_f)*fact1 + xin%co2bdy(i,lchnk,xin%np_f)*fact2
        end do
     end do

  return
end subroutine interp_time_flux

!============================================================================================================

end module co2_data_flux


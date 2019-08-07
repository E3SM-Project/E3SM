!BALLI: take care of the units of CO2!!

module co2_data_flux

!------------------------------------------------------------------------------------------------
! for data reading and interpolation                                           
!------------------------------------------------------------------------------------------------

  use shr_kind_mod,   only : r8 => shr_kind_r8,shr_kind_cx, cl => shr_kind_cl
  use spmd_utils,     only : masterproc
  use ppgrid,         only : begchunk, endchunk, pcols
  use phys_grid,      only : scatter_field_to_chunk, get_ncols_p
  use error_messages, only : alloc_err, handle_ncerr, handle_err
  use cam_abortutils,     only : endrun
  use netcdf,         only : nf90_inq_dimid, nf90_inquire_dimension, nf90_inq_varid, nf90_get_var, nf90_open, &
       nf90_noerr
  use cam_logfile,    only : iulog
  use physics_types,  only: physics_state, physics_state_copy
  use dycore,         only: dycore_is
  use shr_log_mod ,   only: errMsg => shr_log_errMsg
  use input_data_utils, only: time_coordinate
#ifdef CO2_BILIN_REGRID
  use tracer_data,    only : trfld, trfile, trcdata_init, advance_trcdata
#endif
  
#if ( defined SPMD )
  use mpishorthand, only: mpicom, mpiint, mpir8
#endif

  implicit none

! public data

! public type

  public co2_data_flux_type

! public interface

  public read_data_flux
  public co2_data_flux_init
  public interp_time_flux
  public co2_data_flux_advance

! private data

  private
!  integer,  parameter :: totsz=2000           ! number greater than data time sample
  real(r8), parameter :: daysperyear = 365.0_r8  ! Number of days in a year         
  integer :: lonsiz  ! size of longitude dimension, if dataset is 2d(lat,lon), in CAM grid
  integer :: latsiz  ! size of latitude dimension

  integer :: ncolsiz ! size of number of columns (if SE grid)

  !Following state data type is declared so that we can send state as an argument 
  !in advance_trcdata call. Only "pint" (or pmid) is used from state variable
  !to facilitate vertical interpolation. Therefore a state computed at the model 
  !initialization should suffice
  type(physics_state),pointer :: state_at_init(:)
 
!--------------------------------------------------------------------------------------------------
type :: co2_data_flux_type          

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
  integer, pointer :: date_f(:)      ! Date on dataset (YYYYMMDD)
  integer, pointer :: sec_f(:)       ! seconds of date on dataset (0-86399) 
  character(len=cl) :: filename   ! dataset name
  integer :: ncid_f             ! netcdf id for dataset       
  integer :: fluxid             ! netcdf id for dataset flux      
  real(r8), pointer :: xvar(:,:,:) ! work space for dataset  
  type(time_coordinate) :: time_coord
  character(len=cl)     :: varname
  logical               :: initialized
#ifdef CO2_BILIN_REGRID
  type(trfld),  pointer, dimension(:)   :: fields
  type(trfile)                          :: file
#endif
 
end type co2_data_flux_type

! dimension names for physics grid (physgrid)
logical           :: dimnames_set = .false.
character(len=8)  :: dim1name, dim2name

!===============================================================================
contains
!===============================================================================

subroutine co2_data_flux_init (input_file, varname, xin)

!-------------------------------------------------------------------------------
! Initialize co2_data_flux_type instance
!   including initial read of input and interpolation to the current timestep
!-------------------------------------------------------------------------------

   use ioFileMod,        only: getfil
   use ppgrid,           only: begchunk, endchunk, pcols
   use cam_grid_support, only: cam_grid_id, cam_grid_check
   use cam_grid_support, only: cam_grid_get_dim_names
   use shr_log_mod,      only: errMsg => shr_log_errMsg

   ! Arguments
   character(len=*),          intent(in)    :: input_file
   character(len=*),          intent(in)    :: varname
   type(co2_data_flux_type),  intent(inout) :: xin

   ! Local variables
   integer  :: grid_id
   real(r8) :: dtime
   !----------------------------------------------------------------------------

   if (.not. dimnames_set) then
      grid_id = cam_grid_id('physgrid')
      if (.not. cam_grid_check(grid_id)) then
         call endrun('ERROR: no "physgrid" grid:'//errmsg(__FILE__,__LINE__))
      endif
      call cam_grid_get_dim_names(grid_id, dim1name, dim2name)
      dimnames_set = .true.
   end if

   call getfil(input_file, xin%filename)
   xin%varname = varname
   xin%initialized = .false.

   dtime = 1.0_r8 - 200.0_r8 / 86400.0_r8
   call xin%time_coord%initialize(input_file, delta_days=dtime)

   allocate( xin%co2bdy(pcols,begchunk:endchunk,2), &
             xin%co2flx(pcols,begchunk:endchunk)    )
   call co2_data_flux_advance(xin)

   xin%initialized = .true.

end subroutine co2_data_flux_init





subroutine read_data_flux (input_file, xin, state, pbuf2d)
 
!-------------------------------------------------------------------------------              
! Do initial read of time-varying 2d(lat,lon NetCDF dataset, 
! reading two data bracketing current timestep
!-------------------------------------------------------------------------------
 
  use time_manager, only : get_curr_date, get_curr_calday, &
                           is_perpetual, get_perp_date, get_step_size, is_first_step
  use ioFileMod,    only : getfil
  
  use physics_buffer, only: physics_buffer_desc
  use dyn_grid,     only: get_horiz_grid_dim_d
  use cam_pio_utils,    only: cam_pio_openfile
  use pio,              only: file_desc_t, pio_nowrite, pio_closefile
  use cam_grid_support, only: cam_grid_id, cam_grid_check, cam_grid_get_dim_names
  use ncdio_atm,        only: infld




  implicit none
 
!---------------------------Common blocks-------------------------------
! Dummy arguments
  character(len=*),   intent(in)    :: input_file
  TYPE(co2_data_flux_type),  intent(inout) :: xin   
  type(physics_state),       pointer :: state(:)
  type(physics_buffer_desc), pointer :: pbuf2d(:,:)
 
  character(len = shr_kind_cx) :: msg
  integer lonid                 ! netcdf id for longitude variable
  integer latid                 ! netcdf id for latitude variable
  integer timeid                ! netcdf id for time variable
  integer dateid                ! netcdf id for date variable
  integer secid                 ! netcdf id for seconds variable
  integer londimid              ! netcdf id for longitude variable
  integer latdimid              ! netcdf id for latitude variable
  integer ncoldimid             ! netcdf id for columns variable
  integer :: hdim1_d,hdim2_d    ! model grid size
 
  integer dtime                 ! timestep size [seconds]
  integer cnt3(3)               ! array of counts for each dimension
  integer strt3(3)              ! array of starting indices
  integer n, c                  ! indices
  integer j                     ! latitude index
  integer istat                 ! error return
  integer  :: yr, mon, day      ! components of a date
  integer  :: ncdate            ! current date in integer format [yyyymmdd]
  integer  :: ncsec             ! current time of day [seconds]
  real(r8) calday               ! calendar day (includes yr if no cycling)
  real(r8) caldayloc            ! calendar day (includes yr if no cycling)
#ifdef CO2_BILIN_REGRID
  character(len=32)  :: specifier(1) = ''
#endif
  type(file_desc_t) :: fh_co2_data_flux
  character(len=8)  :: dim1name, dim2name
  integer  :: grid_id
  logical           :: found
  real(r8) :: dtime1



! Allocate space for data.
 
  allocate( xin%co2flx(pcols,begchunk:endchunk), stat=istat )
  call alloc_err( istat, 'CO2FLUX_READ', 'co2flx', &
       pcols*(endchunk-begchunk+1) )


#ifdef CO2_BILIN_REGRID
  allocate (xin%file%in_pbuf(1))
  xin%file%in_pbuf(1) = .false.

  specifier(1)    = 'CO2_flux' !name of variable to read from file

! Open file and initialize "fields" and "file" derived types
! Some of the arguments passed here are hardwired which can be replaced with variables
! if need be.
  call trcdata_init(specifier, input_file , '', '', xin%fields, xin%file, .false., 0, 0, 0, 'SERIAL')
  
  !allocate state and copy the state at init 
  allocate(state_at_init(begchunk:endchunk))
  do c = begchunk, endchunk
     call physics_state_copy(state(c),state_at_init(c))
  enddo
#else
 
  xin%nm_f = 1
  xin%np_f = 2
 
 
  allocate( xin%co2bdy(pcols,begchunk:endchunk,2), stat=istat )
  call alloc_err( istat, 'CO2FLUX_READ', 'co2bdy', &
       pcols*(endchunk-begchunk+1)*2 )
 
! SPMD: Master does all the work.
 
!  if (masterproc) then
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
     dtime1 = 1.0_r8 - 200.0_r8 / 86400.0_r8
     call xin%time_coord%initialize(input_file, delta_days=dtime1)


     grid_id = cam_grid_id('physgrid')
     call cam_grid_get_dim_names(grid_id, dim1name, dim2name)
     call getfil(input_file, xin%filename)

     call cam_pio_openfile(fh_co2_data_flux, trim(xin%filename), PIO_NOWRITE)

      ! read time-level 1
      ! skip the read if the needed vals are present in time-level 2
     call infld('CO2_flux', fh_co2_data_flux, dim1name, dim2name, &
              1, pcols, begchunk, endchunk, xin%co2bdy(:,:,1), found, &
              gridname='physgrid', timelevel=xin%time_coord%indxs(1))
     
     !varname, ncid, dimname1,                        &
     !  dim1b, dim1e, dim2b, dim2e, field, readvar, gridname, timelevel)
     



     call handle_ncerr( nf90_open(xin%filename, 0, xin%ncid_f),__FILE__,__LINE__)
     write(iulog,*)'CO2FLUX_READ: NCOPN returns id ',xin%ncid_f,' for file ',trim(xin%filename)
 
! Get and check dimension info
! initialize
 
     strt3(1) = 1
     strt3(2) = 1
     strt3(3) = 1

     cnt3(2)  = 1
     cnt3(3)  = 1

     call get_horiz_grid_dim_d(hdim1_d,hdim2_d)
     if( dycore_is('SE') )  then
        istat = nf90_inq_dimid( xin%ncid_f, 'ncol', ncoldimid )
        if(istat /= nf90_noerr) then
           write(msg,*)'Input file:', trim(adjustl(xin%filename)), '  should be on the same grid as the model (SE):',__FILE__,__LINE__
           call endrun(msg)
        endif
        call handle_ncerr( nf90_inquire_dimension( xin%ncid_f, ncoldimid, len=ncolsiz),__FILE__,__LINE__)
        if(hdim1_d .ne. ncolsiz) then
           write(msg,*)'Input file grid size(',ncolsiz,') should be same as model grid size (',hdim1_d,'):',__FILE__,__LINE__
           call endrun(msg)
        endif
        cnt3(1)  = ncolsiz
     elseif( dycore_is('LR')) then
        istat = nf90_inq_dimid( xin%ncid_f, 'lon', londimid )
        if(istat /= nf90_noerr) then
           write(msg,*)'Input file:', trim(adjustl(xin%filename)), '  should be on the same grid as the model (LR or FV):'&
                ,__FILE__,__LINE__
           call endrun(msg)
        endif
        call handle_ncerr( nf90_inq_dimid( xin%ncid_f, 'lat', latdimid ),__FILE__,__LINE__)
        call handle_ncerr( nf90_inquire_dimension( xin%ncid_f, londimid, len=lonsiz),__FILE__,__LINE__)
        call handle_ncerr( nf90_inquire_dimension( xin%ncid_f, latdimid, len=latsiz),__FILE__,__LINE__)
        !if(hdim1_d .ne. lonsiz .or. hdim2_d .ne. latsiz) then
        !   write(msg,*)'Input file grid size(',lonsiz,'x',latsiz,') should be same as model grid size (',hdim1_d,'x',hdim2_d,'):'&
        !        ,__FILE__,__LINE__
        !   call endrun(msg)
        !endif

        cnt3(1)  = lonsiz
        cnt3(2)  = latsiz
     else
        call endrun('Only SE or LR grids are supported currently:'//errmsg(__FILE__,__LINE__))
     endif

     call handle_ncerr( nf90_inq_dimid( xin%ncid_f, 'time',  timeid ),__FILE__,__LINE__)
     call handle_ncerr( nf90_inquire_dimension( xin%ncid_f, timeid,   len=xin%timesz ),__FILE__,__LINE__)

     allocate(xin%date_f(xin%timesz), xin%sec_f(xin%timesz))

! Get data id
     call handle_ncerr( nf90_inq_varid( xin%ncid_f, 'date',     dateid     ),__FILE__,__LINE__)
     call handle_ncerr( nf90_inq_varid( xin%ncid_f, 'datesec',  secid      ),__FILE__,__LINE__)
     call handle_ncerr( nf90_inq_varid( xin%ncid_f, 'CO2_flux', xin%fluxid ),__FILE__,__LINE__)
 
! Retrieve entire date and sec variables.
 
     call handle_ncerr( nf90_get_var ( xin%ncid_f, dateid, xin%date_f ),__FILE__,__LINE__)
     call handle_ncerr( nf90_get_var ( xin%ncid_f, secid,  xin%sec_f  ),__FILE__,__LINE__)

!  endif

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
           call handle_ncerr( nf90_get_var ( xin%ncid_f, xin%fluxid, xin%xvar(:,:,xin%nm_f), strt3, cnt3),__FILE__,__LINE__)
           strt3(3) = xin%np1_f                                      
           call handle_ncerr( nf90_get_var ( xin%ncid_f, xin%fluxid, xin%xvar(:,:,xin%np_f), strt3, cnt3),__FILE__,__LINE__)

           goto 10

        end if
 
     end do

     write(iulog,*)'CO2FLUX_READ: Failed to find dates bracketing ncdate, ncsec=', ncdate, ncsec
     call endrun
 
10   continue
     write(iulog,*)'CO2FLUX_READ: Read ', trim(xin%filename), ' for dates ', xin%date_f(n), xin%sec_f(n), &
          ' and ', xin%date_f(xin%np1_f), xin%sec_f(xin%np1_f)
 
#if (defined SPMD )
     call mpibcast( xin%timesz, 1,     mpiint, 0, mpicom )
     call mpibcast( xin%date_f, xin%timesz, mpiint, 0, mpicom )
     call mpibcast( xin%sec_f,  xin%timesz, mpiint, 0, mpicom )
     call mpibcast( xin%cdayfm, 1,     mpir8 , 0, mpicom )
     call mpibcast( xin%cdayfp, 1,     mpir8,  0, mpicom )
     call mpibcast( xin%np1_f,  1,     mpiint, 0, mpicom )
  else
     call mpibcast( xin%timesz, 1,     mpiint, 0, mpicom )
     allocate(xin%date_f(xin%timesz), xin%sec_f(xin%timesz))
     call mpibcast( xin%date_f, xin%timesz, mpiint, 0, mpicom )
     call mpibcast( xin%sec_f,  xin%timesz, mpiint, 0, mpicom )
     call mpibcast( xin%cdayfm, 1,     mpir8 , 0, mpicom )
     call mpibcast( xin%cdayfp, 1,     mpir8,  0, mpicom )
     call mpibcast( xin%np1_f,  1,     mpiint, 0, mpicom )
#endif
  end if

  call scatter_field_to_chunk ( 1,1,2,cnt3(1), xin%xvar, xin%co2bdy )

#endif
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
#ifdef CO2_BILIN_REGRID
  use ppgrid,           only: pver, pverp
#endif
 
  logical, intent(in), optional     :: prev_timestep ! If using previous timestep, set to true
  TYPE(co2_data_flux_type),  intent(inout) :: xin

!---------------------------Local variables-----------------------------
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
 
#ifdef CO2_BILIN_REGRID
  integer :: icol      ! indices

  !Read next data if needed and interpolate (space and time)
  call advance_trcdata( xin%fields, xin%file, state_at_init)

  !Assign 
  do lchnk   = begchunk, endchunk
     ncol    = get_ncols_p(lchnk)
     do icol = 1, ncol
        xin%co2flx(icol,lchnk) = xin%fields(1)%data(icol, 1,lchnk)
     end do
  end do

#else
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
        
        if( dycore_is('SE') )  then
           cnt3(1)  = ncolsiz
        elseif( dycore_is('LR')) then
           cnt3(1)  = lonsiz
           cnt3(2)  = latsiz
        endif
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
           call handle_ncerr( nf90_get_var ( xin%ncid_f, xin%fluxid, xin%xvar(:,:,xin%np_f), strt3, cnt3 ),__FILE__,__LINE__)
           write(iulog,*)'CO2FLUX_INTERP: Read ', trim(xin%filename),' for date (yyyymmdd) ', xin%date_f(xin%np1_f), &
                     ' sec ', xin%sec_f(xin%np1_f)
        endif
 
        call scatter_field_to_chunk ( 1,1,2,cnt3(1), xin%xvar, xin%co2bdy )
     end if
 
! Determine time interpolation factors.
 
     call get_timeinterp_factors ( co2cyc, xin%np1_f, xin%cdayfm, xin%cdayfp, caldayloc, fact1, fact2, 'CO2FLUX_INTERP:' )

     do lchnk=begchunk,endchunk
        ncol = get_ncols_p(lchnk)
        do i=1,ncol
           xin%co2flx(i,lchnk) = -999.0 !xin%co2bdy(i,lchnk,xin%nm_f)*fact1 + xin%co2bdy(i,lchnk,xin%np_f)*fact2
        end do
     end do
     

#endif
  return
end subroutine interp_time_flux

!============================================================================================================

subroutine co2_data_flux_advance (xin)

!-------------------------------------------------------------------------------
! Advance the contents of a co2_data_flux_type instance
!   including reading new data, if necessary
!-------------------------------------------------------------------------------

   use cam_pio_utils,    only: cam_pio_openfile
   use ncdio_atm,        only: infld
   use pio,              only: file_desc_t, pio_nowrite, pio_closefile
   use ppgrid,           only: begchunk, endchunk, pcols

   ! Arguments
   type(co2_data_flux_type),  intent(inout) :: xin

   ! Local variables
   character(len=*), parameter :: subname = 'co2_data_flux_advance'
   logical           :: read_data
   integer           :: indx2_pre_adv
   type(file_desc_t) :: fh_co2_data_flux
   logical           :: found

   !----------------------------------------------------------------------------

   read_data = xin%time_coord%read_more() .or. .not. xin%initialized

   indx2_pre_adv = xin%time_coord%indxs(2)

   call xin%time_coord%advance()

   if ( read_data ) then

      call cam_pio_openfile(fh_co2_data_flux, trim(xin%filename), PIO_NOWRITE)

      ! read time-level 1
      ! skip the read if the needed vals are present in time-level 2
      if (xin%initialized .and. xin%time_coord%indxs(1) == indx2_pre_adv) then
         xin%co2bdy(:,:,1) = xin%co2bdy(:,:,2)
      else
         call infld(trim(xin%varname), fh_co2_data_flux, dim1name, dim2name, &
              1, pcols, begchunk, endchunk, xin%co2bdy(:,:,1), found, &
              gridname='physgrid', timelevel=xin%time_coord%indxs(1))
         if (.not. found) then
            call endrun(subname // ': ERROR: ' // trim(xin%varname) // ' not found')
         endif
      endif

      ! read time-level 2
      call infld(trim(xin%varname), fh_co2_data_flux, dim1name, dim2name, &
           1, pcols, begchunk, endchunk, xin%co2bdy(:,:,2), found, &
           gridname='physgrid', timelevel=xin%time_coord%indxs(2))
      if (.not. found) then
         call endrun(subname // ': ERROR: ' // trim(xin%varname) // ' not found')
      endif

      call pio_closefile(fh_co2_data_flux)
   endif

   ! interpolate between time-levels
   ! If time:bounds is in the dataset, and the dataset calendar is compatible with CAM's,
   ! then the time_coordinate class will produce time_coord%wghts(2) == 0.0,
   ! generating fluxes that are piecewise constant in time.

   if (xin%time_coord%wghts(2) == 0.0_r8) then
      xin%co2flx(:,:) = xin%co2bdy(:,:,1)
   else
      xin%co2flx(:,:) = xin%co2bdy(:,:,1) + &
           xin%time_coord%wghts(2) * (xin%co2bdy(:,:,2) - xin%co2bdy(:,:,1))
   endif

end subroutine co2_data_flux_advance

end module co2_data_flux


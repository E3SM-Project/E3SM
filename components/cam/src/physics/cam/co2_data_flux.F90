
module co2_data_flux

!------------------------------------------------------------------------------------------------
! for data reading and interpolation                                           
!------------------------------------------------------------------------------------------------

  use shr_kind_mod,   only : r8 => shr_kind_r8
  use spmd_utils,     only : masterproc
  use ppgrid,         only : begchunk, endchunk, pcols
  use phys_grid,      only : scatter_field_to_chunk, get_ncols_p
  use error_messages, only : alloc_err, handle_ncerr, handle_err
  use cam_abortutils,     only : endrun
  use netcdf
  use error_messages, only : handle_ncerr  
  use cam_logfile,    only : iulog
  use tracer_data,  only : trfld, trfile, trcdata_init, advance_trcdata
  
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
  type(trfld), pointer :: fields(:)
  type(trfile)         :: file

!  integer,  parameter :: totsz=2000           ! number greater than data time sample
  real(r8), parameter :: daysperyear = 365.0_r8  ! Number of days in a year         
  integer :: lonsiz  ! size of longitude dimension, dataset is 2d(lat,lon), in CAM grid
  integer :: latsiz  ! size of latitude dimension
 
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
  integer, pointer :: date_f(:)      ! Date on dataset (YYYYMMDD)
  integer, pointer :: sec_f(:)       ! seconds of date on dataset (0-86399) 
  character(len=256) :: locfn   ! dataset name
  integer :: ncid_f             ! netcdf id for dataset       
  integer :: fluxid             ! netcdf id for dataset flux      
  real(r8), pointer :: xvar(:,:,:) ! work space for dataset  
 
END TYPE read_interp
!-------------------------------------------------------------------------------------------------
                       
contains
 
!===============================================================================

subroutine read_data_flux (input_file, xin, state, pbuf2d)
 
!-------------------------------------------------------------------------------              
! Do initial read of time-varying 2d(lat,lon NetCDF dataset, 
! reading two data bracketing current timestep
!-------------------------------------------------------------------------------
 
  use time_manager, only : get_curr_date, get_curr_calday, &
                           is_perpetual, get_perp_date, get_step_size, is_first_step
  use ioFileMod,    only : getfil
  use physics_types,  only: physics_state
  use physics_buffer, only: physics_buffer_desc

  implicit none
 
!---------------------------Common blocks-------------------------------
! Dummy arguments
  character(len=*),   intent(in)    :: input_file
  TYPE(read_interp),  intent(inout) :: xin   
  type(physics_state), pointer       :: state(:)
  type(physics_buffer_desc), pointer :: pbuf2d(:,:)
 
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
  real(r8) calday               ! calendar day (includes yr if no cycling)
  real(r8) caldayloc            ! calendar day (includes yr if no cycling)
 
  character(len=32)  :: specifier(1) = ''

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


  !call trcdata_init( specifier, filename, filelist, datapath, fields, file, &
  !                     rmv_file, cycle_yr, fixed_ymd, fixed_tod, datatype)
  allocate (file%in_pbuf(1))
  file%in_pbuf(1) = .false.
  specifier(1) = 'CO2_flux:CO2_flux'
  call trcdata_init( specifier ,input_file , '', '', fields, file, &
       .false., 0, 0, 0, 'SERIAL' )

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
 
!-----------------------------------------------------------------------


  call advance_trcdata( fields, file)

  do lchnk=begchunk,endchunk
     ncol = get_ncols_p(lchnk)
     do i=1,ncol
        xin%co2flx(i,lchnk) = fields(1)%data(i,1,lchnk)
     end do
  end do

  return
end subroutine interp_time_flux

!============================================================================================================

end module co2_data_flux


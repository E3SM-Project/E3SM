
module co2_data_flux

!------------------------------------------------------------------------------------------------
! for data reading and interpolation                                           
!------------------------------------------------------------------------------------------------

  use shr_kind_mod,   only : r8 => shr_kind_r8
  use ppgrid,         only : begchunk, endchunk, pcols
  use error_messages, only : alloc_err
  use tracer_data,    only : trfld, trfile, trcdata_init, advance_trcdata

  implicit none

  public read_interp
  public read_data_flux
  public interp_time_flux

  private
!--------------------------------------------------------------------------------------------------
TYPE :: read_interp          
  real(r8),     pointer, dimension(:,:) :: co2flx
  type(trfld),  pointer, dimension(:)   :: fields
  type(trfile)                          :: file
END TYPE read_interp
!-------------------------------------------------------------------------------------------------
 
contains
 
!===============================================================================

subroutine read_data_flux (input_file, xin, state, pbuf2d)
 
!-------------------------------------------------------------------------------              
! Do initial read of time-varying 2d(lat,lon NetCDF dataset)
!-------------------------------------------------------------------------------
 
  use physics_types,  only: physics_state
  use physics_buffer, only: physics_buffer_desc

  implicit none
 
!---------------------------Common blocks-------------------------------
! Dummy arguments
  character(len=*),   intent(in)     :: input_file
  TYPE(read_interp),  intent(inout)  :: xin   
  type(physics_state),       pointer :: state(:)
  type(physics_buffer_desc), pointer :: pbuf2d(:,:)
 
  integer istat                 ! error return
  character(len=32)  :: specifier(1) = ''

 
! Allocate space for data.
 
  allocate( xin%co2flx(pcols,begchunk:endchunk), stat=istat )
  call alloc_err( istat, 'CO2FLUX_READ', 'co2flx', pcols*(endchunk-begchunk+1) )

  allocate (xin%file%in_pbuf(1))
  xin%file%in_pbuf(1) = .false.

  specifier(1)    = 'CO2_flux' !name of variable to read from file

! Open file and initialize "fields" and "file" derived types
! Some of the arguments passed here are hardwired which can be replaced with variables
! if need be.
  call trcdata_init(specifier, input_file , '', '', xin%fields, xin%file, .false., 0, 0, 0, 'SERIAL')

  return
end subroutine read_data_flux
 
!===============================================================================

subroutine interp_time_flux (xin)
 
!-----------------------------------------------------------------------
! Time interpolate data to current time.
! Reading in new monthly data if necessary.
!
!-----------------------------------------------------------------------
  use phys_grid,      only : get_ncols_p

  TYPE(read_interp),  intent(inout) :: xin

!---------------------------Local variables-----------------------------
  integer :: icol, lchnk      ! indices
  integer :: ncol             ! number of columns in current chunk
!-----------------------------------------------------------------------

  !Read next data if needed and interpolate (space and time)
  call advance_trcdata( xin%fields, xin%file)

  !Assign 
  do lchnk   = begchunk, endchunk
     ncol    = get_ncols_p(lchnk)
     do icol = 1, ncol
        xin%co2flx(icol,lchnk) = xin%fields(1)%data(icol, 1,lchnk)
     end do
  end do

  return
end subroutine interp_time_flux

!============================================================================================================
end module co2_data_flux

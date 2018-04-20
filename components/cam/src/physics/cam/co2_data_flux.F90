
module co2_data_flux

!------------------------------------------------------------------------------------------------
! for data reading and interpolation                                           
!------------------------------------------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
  use tracer_data,  only : trfld, trfile, trcdata_init, advance_trcdata

  implicit none

! public 
  public read_data_flux
  public interp_time_flux

! private
  private
 
contains
 
!===============================================================================

subroutine read_data_flux (fields, file, input_file, state, pbuf2d)
 
!-------------------------------------------------------------------------------              
! Do initial read of time-varying 2d(lat,lon NetCDF dataset)
!-------------------------------------------------------------------------------
 
  use physics_types,  only: physics_state
  use physics_buffer, only: physics_buffer_desc

  implicit none
 
!---------------------------Common blocks-------------------------------
! Dummy arguments
  type(trfld),               pointer :: fields(:)
  type(trfile),    intent(out)       :: file
  character(len=*),intent(in)        :: input_file
  type(physics_state),       pointer :: state(:)
  type(physics_buffer_desc), pointer :: pbuf2d(:,:)
 
  character(len=32)  :: specifier(1) = ''

 
! Initialize variables needed for tracer_data routines
  allocate (file%in_pbuf(1))
  file%in_pbuf(1) = .false.
  specifier(1)    = 'CO2_flux' !name of variable to read from file

! Open file and initialize "fields" and "file" derived types
! Some of the arguments passed here are hardwired which can be replaced with variables
! if need be.
  call trcdata_init(specifier, input_file , '', '', fields, file, .false., 0, 0, 0, 'SERIAL')

  return
end subroutine read_data_flux
 
!===============================================================================

subroutine interp_time_flux (fields, file, data_flux)
 
!-----------------------------------------------------------------------
! Time interpolate data to current time.
! Reading in new monthly data if necessary.
!
!-----------------------------------------------------------------------
  use phys_grid,      only : get_ncols_p
  use ppgrid,         only : begchunk, endchunk

  type(trfld),               pointer :: fields(:)
  type(trfile), intent(inout)        :: file
  real(r8),     intent(inout)        :: data_flux(:,:)

!---------------------------Local variables-----------------------------
  integer :: icol, lchnk      ! indices
  integer :: ncol             ! number of columns in current chunk
!-----------------------------------------------------------------------

  !Read next data if needed and interpolate (space and time)
  call advance_trcdata( fields, file)

  !Assign 
  do lchnk   = begchunk, endchunk
     ncol    = get_ncols_p(lchnk)
     do icol = 1, ncol
        data_flux(icol, lchnk) = fields(1)%data(icol, 1, lchnk)
     end do
  end do

  return
end subroutine interp_time_flux

!============================================================================================================
end module co2_data_flux

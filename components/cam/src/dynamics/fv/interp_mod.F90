module interp_mod
  use shr_kind_mod,   only : r8=>shr_kind_r8
  use cam_abortutils, only : endrun

  implicit none
  private
  save

  public :: setup_history_interpolation
  public :: set_interp_hfile
  public :: write_interpolated

  interface write_interpolated
     module procedure write_interpolated_scalar
     module procedure write_interpolated_vector
  end interface
  integer, parameter :: nlat=0, nlon=0
contains

  subroutine setup_history_interpolation(interp_ok, mtapes, interp_output,    &
       interp_info)
    use cam_history_support, only: interp_info_t

    ! Dummy arguments
    logical,             intent(inout) :: interp_ok
    integer,             intent(in)    :: mtapes
    logical,             intent(in)    :: interp_output(:)
    type(interp_info_t), intent(inout) :: interp_info(:)

    interp_ok = .false.

  end subroutine setup_history_interpolation

  subroutine set_interp_hfile(hfilenum, interp_info)
    use cam_history_support, only: interp_info_t

    ! Dummy arguments
    integer,             intent(in)    :: hfilenum
    type(interp_info_t), intent(inout) :: interp_info(:)
  end subroutine set_interp_hfile

  subroutine write_interpolated_scalar(File, varid, fld, numlev, data_type, decomp_type) 
    use pio, only : file_desc_t, var_desc_t
    use shr_kind_mod, only : r8=>shr_kind_r8
    implicit none
    type(file_desc_t), intent(inout) :: File
    type(var_desc_t), intent(inout) :: varid
    real(r8), intent(in) :: fld(:,:,:)
    integer, intent(in) :: numlev, data_type, decomp_type
    call endrun('This routine is a stub, you shouldnt get here')

  end subroutine write_interpolated_scalar

  subroutine write_interpolated_vector(File, varidu, varidv, fldu, fldv, numlev, data_type, decomp_type) 
    use pio, only : file_desc_t, var_desc_t
    implicit none
    type(file_desc_t), intent(inout) :: File
    type(var_desc_t), intent(inout) :: varidu, varidv
    real(r8), intent(in) :: fldu(:,:,:), fldv(:,:,:)
    integer, intent(in) :: numlev, data_type, decomp_type
    call endrun('This routine is a stub, you shouldnt get here')

  end subroutine write_interpolated_vector

end module interp_mod

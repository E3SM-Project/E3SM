module interp_mod
  use shr_kind_mod, only : r8=>shr_kind_r8
  use abortutils, only : endrun
  implicit none
  public get_interp_lat, get_interp_lon, setup_history_interpolation, write_interpolated
  public var_is_vector_uvar, var_is_vector_vvar, latlon_interpolation, add_interp_attributes

  interface write_interpolated
     module procedure write_interpolated_scalar
     module procedure write_interpolated_vector
  end interface
  integer, parameter :: nlat=0, nlon=0
contains

  subroutine add_interp_attributes(file)
    use pio, only : file_desc_t
    type(file_desc_t) :: file

    call endrun('This routine is a stub, you shouldnt get here')
    
  end subroutine add_interp_attributes


  subroutine setup_history_interpolation(mtapes)
    integer, intent(in) :: mtapes
    call endrun('This routine is a stub, you shouldnt get here')

  end subroutine setup_history_interpolation

  function latlon_interpolation(t)
    integer, intent(in) :: t
    logical :: latlon_interpolation

    latlon_interpolation = .false.
  end function latlon_interpolation



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

  function var_is_vector_uvar(name)
    character(len=*), intent(in) :: name
    integer ::var_is_vector_uvar
    var_is_vector_uvar=0
  end function var_is_vector_uvar
  function var_is_vector_vvar(name)
    character(len=*), intent(in) :: name
    integer ::var_is_vector_vvar
    var_is_vector_vvar=0
  end function var_is_vector_vvar

  function get_interp_lat() result(thislat)
    real(kind=r8) :: thislat(nlat)
    call endrun('This routine is a stub, you shouldnt get here')

    return
  end function get_interp_lat
  function get_interp_lon() result(thislon)
    real(kind=r8) :: thislon(nlon)
    call endrun('This routine is a stub, you shouldnt get here')

    return
  end function get_interp_lon


end module interp_mod

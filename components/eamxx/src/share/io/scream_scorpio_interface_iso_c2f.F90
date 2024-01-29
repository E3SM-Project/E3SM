module scream_scorpio_interface_iso_c2f
  use iso_c_binding, only: c_int, c_double, c_float, c_bool, c_ptr
  implicit none

!
! This file contains bridges from scream c++ to shoc fortran.
!
contains
!=====================================================================!
  subroutine eam_init_pio_subsystem_c2f(mpicom,compid) bind(c)
    use scream_scorpio_interface, only : eam_init_pio_subsystem
    integer(kind=c_int), value, intent(in) :: mpicom,compid

    call eam_init_pio_subsystem(mpicom,compid)
  end subroutine eam_init_pio_subsystem_c2f
!=====================================================================!
  subroutine eam_pio_finalize_c2f() bind(c)
    use scream_scorpio_interface, only : eam_pio_finalize

    call eam_pio_finalize()
  end subroutine eam_pio_finalize_c2f
!=====================================================================!
  function get_file_ncid_c2f(filename_in) result(ncid) bind(c)
    use scream_scorpio_interface, only : lookup_pio_atm_file, pio_atm_file_t
    type(c_ptr), intent(in)         :: filename_in

    type(pio_atm_file_t), pointer :: atm_file
    character(len=256)      :: filename
    integer(kind=c_int) :: ncid
    logical :: found

    call convert_c_string(filename_in,filename)
    call lookup_pio_atm_file(filename,atm_file,found)
    if (found) then
      ncid = int(atm_file%pioFileDesc%fh,kind=c_int)
    else
      ncid = -1
    endif
  end function get_file_ncid_c2f
!=====================================================================!
  function get_file_mode_c2f(filename_in) result(mode) bind(c)
    use scream_scorpio_interface, only : lookup_pio_atm_file, pio_atm_file_t

    type(c_ptr), intent(in)         :: filename_in

    type(pio_atm_file_t), pointer :: atm_file
    character(len=256)      :: filename
    integer(kind=c_int)     :: mode
    logical :: found

    call convert_c_string(filename_in,filename)
    call lookup_pio_atm_file(filename,atm_file,found)
    if (found) then
      mode = atm_file%purpose
    else
      mode = 0
    endif
  end function get_file_mode_c2f
  function is_file_open_c2f(filename_in,purpose) result(res) bind(c)
    use scream_scorpio_interface, only : lookup_pio_atm_file, pio_atm_file_t

    type(c_ptr), intent(in)         :: filename_in
    integer(kind=c_int), intent(in) :: purpose

    type(pio_atm_file_t), pointer :: atm_file
    character(len=256)      :: filename
    logical (kind=c_bool)   :: res
    logical :: found

    call convert_c_string(filename_in,filename)
    call lookup_pio_atm_file(filename,atm_file,found)
    if (found) then
      res = LOGICAL(purpose .lt. 0 .or. atm_file%purpose .eq. purpose,kind=c_bool)
    else
      res = .false.
    endif
  end function is_file_open_c2f
!=====================================================================!
  subroutine register_file_c2f(filename_in,purpose) bind(c)
    use scream_scorpio_interface, only : register_file
    type(c_ptr), intent(in)         :: filename_in
    integer(kind=c_int), intent(in) :: purpose

    character(len=256)      :: filename

    call convert_c_string(filename_in,filename)
    call register_file(trim(filename),purpose)

  end subroutine register_file_c2f
!=====================================================================!
  subroutine set_decomp_c2f(filename_in) bind(c)
    use scream_scorpio_interface, only : set_decomp
    type(c_ptr), intent(in) :: filename_in

    character(len=256)       :: filename

    call convert_c_string(filename_in,filename)
    call set_decomp(trim(filename))
  end subroutine set_decomp_c2f
!=====================================================================!
  subroutine set_dof_c2f(filename_in,varname_in,dof_len,dof_vec) bind(c)
    use scream_scorpio_interface, only : set_dof, pio_offset_kind
    use iso_c_binding, only: c_int64_t
    type(c_ptr), intent(in)                             :: filename_in
    type(c_ptr), intent(in)                             :: varname_in
    integer(kind=c_int), value, intent(in)              :: dof_len
    integer(kind=c_int64_t), intent(in), dimension(dof_len) :: dof_vec

    character(len=256)            :: filename
    character(len=256)            :: varname
    integer                       :: ii
    integer(kind=pio_offset_kind), allocatable :: dof_vec_f90(:)

    call convert_c_string(filename_in,filename)
    call convert_c_string(varname_in,varname)
    ! Need to add 1 to the dof_vec because C++ starts indices at 0 not 1:
    allocate(dof_vec_f90(dof_len))
    do ii = 1,dof_len
      dof_vec_f90(ii) = dof_vec(ii) + 1
    end do
    call set_dof(trim(filename),trim(varname),dof_len,dof_vec_f90)
    deallocate(dof_vec_f90)
  end subroutine set_dof_c2f
!=====================================================================!
  subroutine eam_pio_closefile_c2f(filename_in) bind(c)
    use scream_scorpio_interface, only : eam_pio_closefile
    type(c_ptr), intent(in) :: filename_in
    character(len=256)      :: filename

    call convert_c_string(filename_in,filename)
    call eam_pio_closefile(trim(filename))

  end subroutine eam_pio_closefile_c2f
!=====================================================================!
  subroutine eam_pio_flush_file_c2f(filename_in) bind(c)
    use scream_scorpio_interface, only : eam_pio_flush_file
    type(c_ptr), intent(in) :: filename_in
    character(len=256)      :: filename

    call convert_c_string(filename_in,filename)
    call eam_pio_flush_file(trim(filename))

  end subroutine eam_pio_flush_file_c2f
!=====================================================================!
  subroutine pio_update_time_c2f(filename_in,time) bind(c)
    use scream_scorpio_interface, only : eam_update_time
    type(c_ptr), intent(in) :: filename_in
    real(kind=c_double), value, intent(in) :: time

    character(len=256)       :: filename

    call convert_c_string(filename_in,filename)
    call eam_update_time(trim(filename),time)

  end subroutine pio_update_time_c2f
!=====================================================================!
  subroutine register_variable_c2f(filename_in, shortname_in, longname_in, &
                                   units_in, numdims, var_dimensions_in,   &
                                   dtype, nc_dtype, pio_decomp_tag_in) bind(c)
    use scream_scorpio_interface, only : register_variable
    type(c_ptr), intent(in)                :: filename_in
    type(c_ptr), intent(in)                :: shortname_in
    type(c_ptr), intent(in)                :: longname_in
    type(c_ptr), intent(in)                :: units_in
    integer(kind=c_int), value, intent(in) :: numdims
    type(c_ptr), intent(in)                :: var_dimensions_in(numdims)
    integer(kind=c_int), value, intent(in) :: dtype, nc_dtype
    type(c_ptr), intent(in)                :: pio_decomp_tag_in

    character(len=256) :: filename
    character(len=256) :: shortname
    character(len=256) :: longname
    character(len=256) :: units
    character(len=256) :: var_dimensions(numdims)
    character(len=256) :: pio_decomp_tag
    integer            :: ii

    call convert_c_string(filename_in,filename)
    call convert_c_string(shortname_in,shortname)
    call convert_c_string(longname_in,longname)
    call convert_c_string(units_in,units)
    call convert_c_string(pio_decomp_tag_in,pio_decomp_tag)
    do ii = 1,numdims
      call convert_c_string(var_dimensions_in(ii), var_dimensions(ii))
    end do

    call register_variable(filename,shortname,longname,units,numdims,var_dimensions,dtype,nc_dtype,pio_decomp_tag)

  end subroutine register_variable_c2f
!=====================================================================!
  subroutine set_variable_metadata_char_c2f(filename_in, varname_in, metaname_in, metaval_in) bind(c)
    use scream_scorpio_interface, only : set_variable_metadata_char
    type(c_ptr), intent(in)                :: filename_in
    type(c_ptr), intent(in)                :: varname_in
    type(c_ptr), intent(in)                :: metaname_in
    type(c_ptr), intent(in)                :: metaval_in

    character(len=256) :: filename
    character(len=256) :: varname
    character(len=256) :: metaname
    character(len=256) :: metaval

    call convert_c_string(filename_in,filename)
    call convert_c_string(varname_in,varname)
    call convert_c_string(metaname_in,metaname)
    call convert_c_string(metaval_in,metaval)

    call set_variable_metadata_char(filename,varname,metaname,metaval)

  end subroutine set_variable_metadata_char_c2f
!=====================================================================!
  subroutine set_variable_metadata_float_c2f(filename_in, varname_in, metaname_in, metaval_in) bind(c)
    use scream_scorpio_interface, only : set_variable_metadata_float
    type(c_ptr), intent(in)                :: filename_in
    type(c_ptr), intent(in)                :: varname_in
    type(c_ptr), intent(in)                :: metaname_in
    real(kind=c_float), value, intent(in)  :: metaval_in

    character(len=256) :: filename
    character(len=256) :: varname
    character(len=256) :: metaname

    call convert_c_string(filename_in,filename)
    call convert_c_string(varname_in,varname)
    call convert_c_string(metaname_in,metaname)

    call set_variable_metadata_float(filename,varname,metaname,metaval_in)

  end subroutine set_variable_metadata_float_c2f
!=====================================================================!
  subroutine set_variable_metadata_double_c2f(filename_in, varname_in, metaname_in, metaval_in) bind(c)
    use scream_scorpio_interface, only : set_variable_metadata_double
    type(c_ptr), intent(in)                :: filename_in
    type(c_ptr), intent(in)                :: varname_in
    type(c_ptr), intent(in)                :: metaname_in
    real(kind=c_double), value, intent(in) :: metaval_in

    character(len=256) :: filename
    character(len=256) :: varname
    character(len=256) :: metaname

    call convert_c_string(filename_in,filename)
    call convert_c_string(varname_in,varname)
    call convert_c_string(metaname_in,metaname)

    call set_variable_metadata_double(filename,varname,metaname,metaval_in)

  end subroutine set_variable_metadata_double_c2f
!=====================================================================!
  function get_variable_metadata_float_c2f(filename_in, varname_in, metaname_in) result(metaval_out) bind(c)
    use scream_scorpio_interface, only : get_variable_metadata_float
    type(c_ptr), intent(in)                :: filename_in
    type(c_ptr), intent(in)                :: varname_in
    type(c_ptr), intent(in)                :: metaname_in
    real(kind=c_float)                     :: metaval_out

    character(len=256) :: filename
    character(len=256) :: varname
    character(len=256) :: metaname

    call convert_c_string(filename_in,filename)
    call convert_c_string(varname_in,varname)
    call convert_c_string(metaname_in,metaname)

    metaval_out = get_variable_metadata_float(filename,varname,metaname)

  end function get_variable_metadata_float_c2f
!=====================================================================!
  function get_variable_metadata_double_c2f(filename_in, varname_in, metaname_in) result(metaval_out) bind(c)
    use scream_scorpio_interface, only : get_variable_metadata_double
    type(c_ptr), intent(in)                :: filename_in
    type(c_ptr), intent(in)                :: varname_in
    type(c_ptr), intent(in)                :: metaname_in
    real(kind=c_double)                    :: metaval_out

    character(len=256) :: filename
    character(len=256) :: varname
    character(len=256) :: metaname

    call convert_c_string(filename_in,filename)
    call convert_c_string(varname_in,varname)
    call convert_c_string(metaname_in,metaname)

    metaval_out = get_variable_metadata_double(filename,varname,metaname)

  end function get_variable_metadata_double_c2f
!=====================================================================!
  subroutine get_variable_metadata_char_c2f(filename_in, varname_in, metaname_in, metaval_out) bind(c)
    use scream_scorpio_interface, only : get_variable_metadata_char
    use iso_c_binding, only: C_NULL_CHAR, c_char, c_f_pointer
    type(c_ptr), intent(in) :: filename_in
    type(c_ptr), intent(in) :: varname_in
    type(c_ptr), intent(in) :: metaname_in
    type(c_ptr), intent(in) :: metaval_out

    character(len=256) :: filename
    character(len=256) :: varname
    character(len=256) :: metaname
    character(len=256) :: metaval
    character(len=256, kind=c_char), pointer :: temp_string
    integer :: slen

    call c_f_pointer(metaval_out,temp_string)

    call convert_c_string(filename_in,filename)
    call convert_c_string(varname_in,varname)
    call convert_c_string(metaname_in,metaname)

    metaval = get_variable_metadata_char(filename,varname,metaname)

    slen = len(trim(metaval))
    ! If string is 255 or less, add terminating char. If not, it's still
    ! ok (the C++ string will have length=max_length=256)
    if (slen .le. 255) then
      temp_string = trim(metaval) // C_NULL_CHAR
    else
      temp_string = metaval
    endif
  end subroutine get_variable_metadata_char_c2f
!=====================================================================!
  subroutine register_dimension_c2f(filename_in, shortname_in, longname_in, length, partitioned) bind(c)
    use scream_scorpio_interface, only : register_dimension

    type(c_ptr), intent(in)                 :: filename_in
    type(c_ptr), intent(in)                 :: shortname_in
    type(c_ptr), intent(in)                 :: longname_in
    integer(kind=c_int), value, intent(in)  :: length
    logical(kind=c_bool), value, intent(in) :: partitioned

    character(len=256) :: filename
    character(len=256) :: shortname
    character(len=256) :: longname

    call convert_c_string(filename_in,filename)
    call convert_c_string(shortname_in,shortname)
    call convert_c_string(longname_in,longname)
    call register_dimension(filename,shortname,longname,length,LOGICAL(partitioned))

  end subroutine register_dimension_c2f
!=====================================================================!
  function read_curr_time_c2f(filename_in) result(val) bind(c)
    use scream_scorpio_interface, only : read_time_at_index
    type(c_ptr), intent(in)                :: filename_in
    real(kind=c_double)                    :: val

    character(len=256) :: filename

    call convert_c_string(filename_in,filename)
    val        = read_time_at_index(filename)

  end function read_curr_time_c2f
!=====================================================================!
  function read_time_at_index_c2f(filename_in,time_index) result(val) bind(c)
    use scream_scorpio_interface, only : read_time_at_index
    type(c_ptr), intent(in)                :: filename_in
    integer(kind=c_int), intent(in)        :: time_index ! zero-based
    real(kind=c_double)                    :: val

    character(len=256) :: filename

    call convert_c_string(filename_in,filename)
    val = read_time_at_index(filename,time_index)

  end function read_time_at_index_c2f
!=====================================================================!
  function is_enddef_c2f(filename_in) bind(c) result(enddef)
    use scream_scorpio_interface, only : lookup_pio_atm_file, pio_atm_file_t
    type(c_ptr), intent(in) :: filename_in

    type(pio_atm_file_t), pointer :: atm_file
    character(len=256)    :: filename
    logical (kind=c_bool) :: enddef
    logical :: found

    call convert_c_string(filename_in,filename)
    call lookup_pio_atm_file(filename,atm_file,found)
    if (found) then
      enddef = LOGICAL(atm_file%is_enddef, kind=c_bool)
    endif
  end function is_enddef_c2f
!=====================================================================!
  subroutine eam_pio_enddef_c2f(filename_in) bind(c)
    use scream_scorpio_interface, only : eam_pio_enddef
    type(c_ptr), intent(in) :: filename_in

    character(len=256)      :: filename

    call convert_c_string(filename_in,filename)
    call eam_pio_enddef(filename)
  end subroutine eam_pio_enddef_c2f
!=====================================================================!
  subroutine eam_pio_redef_c2f(filename_in) bind(c)
    use scream_scorpio_interface, only : eam_pio_redef
    type(c_ptr), intent(in) :: filename_in

    character(len=256)      :: filename

    call convert_c_string(filename_in,filename)
    call eam_pio_redef(filename)
  end subroutine eam_pio_redef_c2f
!=====================================================================!
  subroutine convert_c_string(c_string_ptr,f_string)
    use iso_c_binding, only: c_f_pointer, C_NULL_CHAR
  ! Purpose: To convert a c_string pointer to the proper fortran string format.
    type(c_ptr), intent(in) :: c_string_ptr
    character(len=256), intent(out) :: f_string
    character(len=256), pointer :: temp_string
    integer :: str_len

    call c_f_pointer(c_string_ptr,temp_string)
    str_len = index(temp_string, C_NULL_CHAR) - 1
    f_string = trim(temp_string(1:str_len))

  end subroutine convert_c_string
!=====================================================================!
  subroutine grid_write_data_array_c2f_int(filename_in,varname_in,buf,buf_size) bind(c)
    use scream_scorpio_interface, only: grid_write_data_array

    type(c_ptr), intent(in) :: filename_in
    type(c_ptr), intent(in) :: varname_in
    integer(kind=c_int), intent(in), value :: buf_size
    integer(kind=c_int), intent(in) :: buf(buf_size)

    character(len=256) :: filename
    character(len=256) :: varname

    call convert_c_string(filename_in,filename)
    call convert_c_string(varname_in,varname)
    call grid_write_data_array(filename,varname,buf,buf_size)

  end subroutine grid_write_data_array_c2f_int
  subroutine grid_write_data_array_c2f_float(filename_in,varname_in,buf,buf_size) bind(c)
    use scream_scorpio_interface, only: grid_write_data_array

    type(c_ptr), intent(in) :: filename_in
    type(c_ptr), intent(in) :: varname_in
    integer(kind=c_int), intent(in), value :: buf_size
    real(kind=c_float), intent(in) :: buf(buf_size)

    character(len=256) :: filename
    character(len=256) :: varname

    call convert_c_string(filename_in,filename)
    call convert_c_string(varname_in,varname)
    call grid_write_data_array(filename,varname,buf,buf_size)

  end subroutine grid_write_data_array_c2f_float
  subroutine grid_write_data_array_c2f_double(filename_in,varname_in,buf,buf_size) bind(c)
    use scream_scorpio_interface, only: grid_write_data_array

    type(c_ptr), intent(in) :: filename_in
    type(c_ptr), intent(in) :: varname_in
    integer(kind=c_int), intent(in), value :: buf_size
    real(kind=c_double), intent(in) :: buf(buf_size)

    character(len=256) :: filename
    character(len=256) :: varname

    call convert_c_string(filename_in,filename)
    call convert_c_string(varname_in,varname)
    call grid_write_data_array(filename,varname,buf,buf_size)

  end subroutine grid_write_data_array_c2f_double
!=====================================================================!
  subroutine grid_read_data_array_c2f_int(filename_in,varname_in,time_index,buf,buf_size) bind(c)
    use scream_scorpio_interface, only: grid_read_data_array

    type(c_ptr), intent(in) :: filename_in
    type(c_ptr), intent(in) :: varname_in
    integer(kind=c_int), value, intent(in) :: time_index ! zero-based
    integer(kind=c_int), intent(in), value :: buf_size
    integer(kind=c_int), intent(out) :: buf(buf_size)

    character(len=256) :: filename
    character(len=256) :: varname

    call convert_c_string(filename_in,filename)
    call convert_c_string(varname_in,varname)
    call grid_read_data_array(filename,varname,buf,buf_size,time_index+1)

  end subroutine grid_read_data_array_c2f_int
!=====================================================================!
  subroutine grid_read_data_array_c2f_float(filename_in,varname_in,time_index,buf,buf_size) bind(c)
    use scream_scorpio_interface, only: grid_read_data_array

    type(c_ptr), intent(in) :: filename_in
    type(c_ptr), intent(in) :: varname_in
    integer(kind=c_int), value, intent(in) :: time_index ! zero-based
    integer(kind=c_int), intent(in), value :: buf_size
    real(kind=c_float), intent(out) :: buf(buf_size)

    character(len=256) :: filename
    character(len=256) :: varname

    call convert_c_string(filename_in,filename)
    call convert_c_string(varname_in,varname)
    call grid_read_data_array(filename,varname,buf,buf_size,time_index+1)

  end subroutine grid_read_data_array_c2f_float
!=====================================================================!
  subroutine grid_read_data_array_c2f_double(filename_in,varname_in,time_index,buf,buf_size) bind(c)
    use scream_scorpio_interface, only: grid_read_data_array

    type(c_ptr), intent(in) :: filename_in
    type(c_ptr), intent(in) :: varname_in
    integer(kind=c_int), value, intent(in) :: time_index ! zero-based
    integer(kind=c_int), intent(in), value :: buf_size
    real(kind=c_double), intent(out) :: buf(buf_size)

    character(len=256) :: filename
    character(len=256) :: varname

    call convert_c_string(filename_in,filename)
    call convert_c_string(varname_in,varname)
    call grid_read_data_array(filename,varname,buf,buf_size,time_index+1)

  end subroutine grid_read_data_array_c2f_double
!=====================================================================!
end module scream_scorpio_interface_iso_c2f

module scream_scorpio_interface_iso_c
  use iso_c_binding
  implicit none
     
#include "scream_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
#else
# define c_real c_float
#endif

!
! This file contains bridges from scream c++ to shoc fortran. 
!

contains
!!=====================================================================!
!  subroutine eam_init_pio_1_c(mpicom,compid) bind(c)
!    use scream_scorpio_interface, only : eam_init_pio_subsystem, register_outfile
!    integer(kind=c_int), value, intent(in) :: mpicom
!    integer(kind=c_int), value, intent(in) :: compid
!    write(*,*) "ASD - Got this far: ", mpicom, compid
!!    call eam_init_pio_subsystem(mpicom,compid)
!    call register_outfile("example_pio_structured.nc")
!    call register_outfile("example_pio_structured_v2.nc")
!  end subroutine eam_init_pio_1_c
!=====================================================================!
  subroutine eam_init_pio_subsystem_c(mpicom,compid) bind(c)
    use scream_scorpio_interface, only : eam_init_pio_subsystem
    integer(kind=c_int), value, intent(in) :: mpicom
    integer(kind=c_int), value, intent(in) :: compid

    call eam_init_pio_subsystem(mpicom,compid)
  end subroutine eam_init_pio_subsystem_c
!=====================================================================!
  subroutine register_outfile_c(filename_in) bind(c)
    use scream_scorpio_interface, only : register_outfile
    !character(len=*,kind=c_char), intent(in) :: filename
    type(c_ptr), intent(in) :: filename_in

    character(len=256)       :: filename

    call convert_c_string(filename_in,filename)
    call register_outfile(trim(filename))

  end subroutine register_outfile_c
!=====================================================================!
  subroutine register_dimension_c(filename_in,length) bind(c) !, shortname_in, longname_in, length) bind(c)
    use scream_scorpio_interface, only : register_dimension
    type(c_ptr), intent(in) :: filename_in
!    type(c_ptr), intent(in) :: shortname_in
!    type(c_ptr), intent(in) :: longname_in
    integer(kind=c_int), intent(in) :: length

    character(len=256)       :: filename
    character(len=256)       :: shortname = "x"
    character(len=256)       :: longname  = "horizontal distance"

    write(*,*) "ASD - in : ", length
    call convert_c_string(filename_in,filename)
    write(*,*) "ASD f90 : ", trim(filename)
!    call convert_c_string(shortname_in,shortname)
!    call convert_c_string(longname_in,longname)
    write(*,*) "ASD - co : ", trim(filename), trim(shortname), trim(longname), length
    call register_dimension(filename,shortname,longname,length)
    
  end subroutine register_dimension_c
!=====================================================================!
  subroutine eam_init_pio_2_c() bind(c)
    use scream_scorpio_interface, only : eam_init_pio_2
    call eam_init_pio_2()
  end subroutine eam_init_pio_2_c
!=====================================================================!
  subroutine eam_history_write_c() bind(c)
    use scream_scorpio_interface, only : eam_history_write
    call eam_history_write()
  end subroutine eam_history_write_c
!======================================================================
  subroutine convert_c_string(c_string_ptr,f_string)
  ! Purpose: To convert a c_string pointer to the proper fortran string format.
    type(c_ptr), intent(in) :: c_string_ptr
    character(len=256), intent(out) :: f_string
    character(len=256), pointer :: temp_string
    integer :: str_len

    call c_f_pointer(c_string_ptr,temp_string)
    str_len = index(temp_string, C_NULL_CHAR) - 1
    f_string = trim(temp_string(1:str_len))
    
    return
  end subroutine convert_c_string
!=====================================================================!

!  subroutine grid_write_data_array_c(filename,hbuf,varname)
!    use scream_scorpio_interface, only: grid_write_data_array
!
!    character(kind=c_char,len=*), intent(in) :: filename
!    real(kind=c_real), intent(in),
!    character(kind=c_char,len=*), intent(in) :: varname
!
!    call grid_write_data_array()
!
!  end subroutine grid_write_data_array_c

end module scream_scorpio_interface_iso_c

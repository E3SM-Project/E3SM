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
  subroutine register_outfile_c(filename) bind(c)
    use scream_scorpio_interface, only : register_outfile
    character(len=*,kind=c_char), intent(in) :: filename

    call register_outfile(filename)
  end subroutine register_outfile_c
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

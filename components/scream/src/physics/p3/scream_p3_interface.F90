#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module scream_p3_interface_mod

  use iso_c_binding, only: c_ptr, c_f_pointer, c_int, c_double, c_bool,C_NULL_CHAR

  public :: p3_init_f90
  public :: p3_main_f90
  public :: p3_finalize_f90

  integer :: test

contains

  !====================================================================!
  subroutine p3_init_f90 () bind(c)

    test = 0
    print *, 'P3 test = ', test

  end subroutine p3_init_f90
  !====================================================================!
  subroutine p3_main_f90 () bind(c)

    test = test + 1
    print *, 'P3 test = ', test

  end subroutine p3_main_f90
  !====================================================================!
  subroutine p3_finalize_f90 () bind(c)

    test = -999
    print *, 'P3 test = ', test

  end subroutine p3_finalize_f90
  !====================================================================!

end module scream_p3_interface_mod

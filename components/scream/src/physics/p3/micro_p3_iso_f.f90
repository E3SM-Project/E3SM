module micro_p3_iso_f
  use iso_c_binding
  implicit none

#include "scream_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
#else
# define c_real c_float
#endif

!
! This file contains bridges from micro_p3 fortran to scream c++.
!

interface

  subroutine find_lookuptable_indices_1a_f(dumi,dumjj,dumii,dumzz,dum1,dum4,dum5,dum6,      &
       qitot,nitot,qirim,rhop) bind(C)
    use iso_c_binding

    ! arguments:
    integer(kind=c_int), intent(out) :: dumi,dumjj,dumii,dumzz
    real(kind=c_real),   intent(out) :: dum1,dum4,dum5,dum6
    real(kind=c_real),   value, intent(in)  :: qitot,nitot,qirim,rhop
  end subroutine find_lookuptable_indices_1a_f

  subroutine find_lookuptable_indices_1b_f(dumj,dum3,qr,nr) bind(C)
    use iso_c_binding

    integer(kind=c_int), intent(out) :: dumj
    real(kind=c_real),   intent(out) :: dum3
    real(kind=c_real),   value, intent(in) :: qr, nr
  end subroutine find_lookuptable_indices_1b_f

  subroutine access_lookup_table_f(dumjj,dumii,dumi,index,dum1,dum4,dum5,proc) bind(C)
    use iso_c_binding

    integer(kind=c_int), value, intent(in) :: dumjj, dumii, dumi, index
    real(kind=c_real),   value, intent(in) :: dum1, dum4, dum5
    real(kind=c_real),   intent(out) :: proc
  end subroutine access_lookup_table_f

  subroutine access_lookup_table_coll_f(dumjj,dumii,dumj,dumi,index,dum1,dum3,dum4,dum5,proc) bind(C)
    use iso_c_binding

    integer(kind=c_int), value, intent(in) :: dumjj,dumii,dumj,dumi,index
    real(kind=c_real),   value, intent(in) :: dum1,dum3,dum4,dum5
    real(kind=c_real),   intent(out) :: proc
  end subroutine access_lookup_table_coll_f

end interface

end module micro_p3_iso_f

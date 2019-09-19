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

  subroutine get_cloud_dsd2_f(qc,nc,mu_c,rho,nu,lamc,cdist,cdist1,lcldm) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in)         :: qc,rho,lcldm
    real(kind=c_real), intent(inout)             :: nc
    real(kind=c_real), intent(out)               :: mu_c,nu,lamc,cdist,cdist1
  end subroutine get_cloud_dsd2_f

  subroutine get_rain_dsd2_f(qr,nr,mu_r,lamr,cdistr,logn0r,rcldm) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: qr,rcldm
    real(kind=c_real), intent(inout)     :: nr
    real(kind=c_real), intent(out)       :: lamr,mu_r,cdistr,logn0r
  end subroutine get_rain_dsd2_f

  !
  ! These are some routine math operations that are not BFB between
  ! fortran and C++ on all platforms, so fortran will need to use
  ! the C++ versions in order to stay BFB.
  !

  function cxx_pow(base, exp) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in)  :: base
    real(kind=c_real), value, intent(in)  :: exp

    ! return
    real(kind=c_real)               :: cxx_pow
  end function cxx_pow

  function cxx_cbrt(base) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in)  :: base

    ! return
    real(kind=c_real)               :: cxx_cbrt
  end function cxx_cbrt

  function cxx_gamma(input) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: input

    ! return
    real(kind=c_real)            :: cxx_gamma
  end function cxx_gamma

  function cxx_log(input) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: input

    ! return
    real(kind=c_real)            :: cxx_log
  end function cxx_log

  function cxx_log10(input) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: input

    ! return
    real(kind=c_real)            :: cxx_log10
  end function cxx_log10

  function cxx_exp(input) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: input

    ! return
    real(kind=c_real)            :: cxx_exp
  end function cxx_exp

end interface

end module micro_p3_iso_f

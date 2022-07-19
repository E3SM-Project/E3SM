
module scream_trcmix_c2f

  !----------------------------------------------------------------------------------------------
  ! Purpose: bridge to trcmix
  !----------------------------------------------------------------------------------------------

  use iso_c_binding
  implicit none

#include "scream_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
#else
# define c_real c_float
#endif

contains

  subroutine trcmix_c2f(name, ncol, pcols, pver, clat, pmid, q, mwdry, mwco2, mwn2o, mwch4, mwf11, mwf12) bind(c)

    use scream_trcmix, only : trcmix

    ! Arguments
    type(c_ptr),                intent(in)                          :: name  ! constituent name
    integer(kind=c_int), value, intent(in)                          :: ncol  ! number of columns
    integer(kind=c_int), value, intent(in)                          :: pcols ! Maximum number of columns in a chunk (physics data structure)
    integer(kind=c_int), value, intent(in)                          :: pver  ! PLEV
    real(kind=c_real),          intent(in),  dimension(pcols)       :: clat  ! latitude in radians for columns
    real(kind=c_real),          intent(in),  dimension(pcols,pver)  :: pmid  ! model pressures
    real(kind=c_real),          intent(out), dimension(pcols,pver)  :: q     ! constituent mass mixing ratio

    real(kind=c_real), value, intent(in) :: mwdry
    real(kind=c_real), value, intent(in) :: mwco2
    real(kind=c_real), value, intent(in) :: mwn2o
    real(kind=c_real), value, intent(in) :: mwch4
    real(kind=c_real), value, intent(in) :: mwf11
    real(kind=c_real), value, intent(in) :: mwf12

    character(len=256), pointer :: name_f
    integer :: len

    call c_f_pointer(name, name_f)
    len = index(name_f, C_NULL_CHAR) - 1

    call trcmix(name_f(1:len), ncol, pcols, pver, clat, pmid, q, mwdry, mwco2, mwn2o, mwch4, mwf11, mwf12)

  end subroutine trcmix_c2f

end module scream_trcmix_c2f

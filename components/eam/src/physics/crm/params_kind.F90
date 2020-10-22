
module params_kind
#ifndef MMF_STANDALONE
  use shr_kind_mod,  only: shr_kind_r8
#endif
  implicit none

#ifdef CRM_SINGLE_PRECISION
  integer, parameter :: crm_rknd = selected_real_kind( 6) ! 4 byte real
#else
  ! default precision of real - kind(1.d0)
  integer, parameter :: crm_rknd = selected_real_kind(12) ! 8 byte real
#endif

#ifdef MMF_STANDALONE
  integer, parameter :: r8 = selected_real_kind(12)
#else
  integer, parameter :: r8 = shr_kind_r8
#endif

end module params_kind

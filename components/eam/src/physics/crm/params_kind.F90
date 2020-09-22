
module params_kind
  implicit none
#ifdef CRM_SINGLE_PRECISION
  integer, parameter :: crm_rknd = selected_real_kind( 6) ! 4 byte real
#else
  ! default precision of real - kind(1.d0)
  integer, parameter :: crm_rknd = selected_real_kind(12) ! 8 byte real
#endif
end module params_kind

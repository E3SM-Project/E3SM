module zm_iso_c
  use iso_c_binding
  implicit none

#include "eamxx_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
#else
# define c_real c_float
#endif

!===================================================================================================
contains
!===================================================================================================

subroutine zm_find_mse_max_c( pcols, ncol, pver, num_msg, msemax_top_k, pergro_active, temperature, zmid, sp_humidity, msemax_klev, mse_max_val ) bind(C)
  use zm_conv_cape,   only: find_mse_max
  use zm_conv_types,  only: zm_const_t, zm_param_t
  use zm_conv_types,  only: zm_param_set_for_testing, zm_const_set_for_testing
  !-----------------------------------------------------------------------------
  ! Interface Arguments
  integer(kind=c_int), value,                intent(in) :: pcols           ! number of atmospheric columns (max)
  integer(kind=c_int), value,                intent(in) :: ncol            ! number of atmospheric columns (actual)
  integer(kind=c_int), value,                intent(in) :: pver            ! number of mid-point vertical levels
  integer(kind=c_int), value,                intent(in) :: num_msg         ! number of missing moisture levels at the top of model
  integer(kind=c_int), dimension(ncol),      intent(in) :: msemax_top_k    ! upper limit index of max MSE search
  logical(kind=c_bool),value,                intent(in) :: pergro_active   ! flag for perturbation growth test (pergro)
  real(kind=c_real),   dimension(ncol,pver), intent(in) :: temperature     ! environement temperature
  real(kind=c_real),   dimension(ncol,pver), intent(in) :: zmid            ! height/altitude at mid-levels
  real(kind=c_real),   dimension(ncol,pver), intent(in) :: sp_humidity     ! specific humidity
  integer(kind=c_int), dimension(ncol),      intent(out):: msemax_klev     ! index of max MSE at parcel launch level
  real(kind=c_real),   dimension(ncol),      intent(out):: mse_max_val     ! value of max MSE at parcel launch level
  !-----------------------------------------------------------------------------
  ! Local Variables
  type(zm_const_t) :: zm_const ! derived type to hold ZM constants
  type(zm_param_t) :: zm_param ! derived type to hold ZM tunable parameters
  logical :: pergro_active_f
  !-----------------------------------------------------------------------------
  call zm_param_set_for_testing(zm_param)
  call zm_const_set_for_testing(zm_const)
  !-----------------------------------------------------------------------------
  pergro_active_f = pergro_active
  call find_mse_max( pcols, ncol, pver, num_msg, msemax_top_k, pergro_active_f, temperature, zmid, sp_humidity, zm_const, zm_param, msemax_klev, mse_max_val )
  !-----------------------------------------------------------------------------
end subroutine zm_find_mse_max_c

!===================================================================================================

end module zm_iso_c

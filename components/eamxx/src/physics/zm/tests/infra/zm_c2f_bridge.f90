module zm_c2f_bridge
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

subroutine zm_common_init_bridge_f() bind(C)
  use zm_eamxx_bridge_wv_saturation, only: wv_sat_init

  call wv_sat_init()
end subroutine zm_common_init_bridge_f

subroutine zm_common_finalize_bridge_f() bind(C)
  use zm_eamxx_bridge_wv_saturation, only: wv_sat_final

  call wv_sat_final()
end subroutine zm_common_finalize_bridge_f

subroutine ientropy_bridge_f(s, p, qt, t, qst, tfg) bind(C)
  use zm_conv_util, only : ientropy
  use zm_conv_types,  only: zm_const_t, zm_const_set_for_testing

  real(kind=c_real) , value, intent(in) :: s, p, qt, tfg
  real(kind=c_real) , intent(out) :: t, qst

  type(zm_const_t) :: zm_const

  call zm_const_set_for_testing(zm_const)
  call ientropy(1, s, p, qt, t, qst, tfg, zm_const)
end subroutine ientropy_bridge_f

subroutine entropy_bridge_f(tk, p, qtot, entropy_rv) bind(C)
  use zm_conv_util, only : entropy
  use zm_conv_types, only: zm_const_t, zm_const_set_for_testing

  real(kind=c_real) , value, intent(in) :: tk, p, qtot
  real(kind=c_real) , intent(out) :: entropy_rv

  type(zm_const_t) :: zm_const

  call zm_const_set_for_testing(zm_const)
  entropy_rv = entropy(tk, p, qtot, zm_const)
end subroutine entropy_bridge_f

subroutine zm_transport_tracer_bridge_f(pcols, pver, doconvtran, q, ncnst, mu, md, du, eu, ed, dp, jt, mx, ideep, il1g, il2g, fracis, dqdt, dpdry, dt) bind(C)
  use zm_transport, only : zm_transport_tracer

  integer(kind=c_int) , value, intent(in) :: pcols, pver, ncnst, il1g, il2g
  logical(kind=c_bool) , intent(in), dimension(ncnst) :: doconvtran
  real(kind=c_real) , intent(in), dimension(pcols, pver, ncnst) :: q, fracis
  real(kind=c_real) , intent(in), dimension(pcols, pver) :: mu, md, du, eu, ed, dp, dpdry
  integer(kind=c_int) , intent(in), dimension(pcols) :: jt, mx, ideep
  real(kind=c_real) , intent(out), dimension(pcols, pver, ncnst) :: dqdt
  real(kind=c_real) , value, intent(in) :: dt

  call zm_transport_tracer(pcols, pver, doconvtran, q, ncnst, mu, md, du, eu, ed, dp, jt, mx, ideep, il1g, il2g, fracis, dqdt, dpdry, dt)
end subroutine zm_transport_tracer_bridge_f

subroutine zm_transport_momentum_bridge_f(pcols, ncol, pver, pverp, wind_in, nwind, mu, md, du, eu, ed, dp, jt, mx, ideep, il1g, il2g, wind_tend, pguall, pgdall, icwu, icwd, dt, seten) bind(C)
  use zm_transport, only : zm_transport_momentum

  integer(kind=c_int) , value, intent(in) :: pcols, ncol, pver, pverp, nwind, il1g, il2g
  real(kind=c_real) , intent(in), dimension(pcols, pver, nwind) :: wind_in
  real(kind=c_real) , intent(in), dimension(pcols, pver) :: mu, md, du, eu, ed, dp
  integer(kind=c_int) , intent(in), dimension(pcols) :: jt, mx, ideep
  real(kind=c_real) , intent(out), dimension(pcols, pver, nwind) :: wind_tend, pguall, pgdall, icwu, icwd
  real(kind=c_real) , value, intent(in) :: dt
  real(kind=c_real) , intent(out), dimension(pcols, pver) :: seten

  call zm_transport_momentum(pcols, ncol, pver, pverp, wind_in, nwind, mu, md, du, eu, ed, dp, jt, mx, ideep, il1g, il2g, wind_tend, pguall, pgdall, icwu, icwd, dt, seten)
end subroutine zm_transport_momentum_bridge_f

subroutine compute_dilute_cape_bridge_f(pcols, ncol, pver, pverp, num_cin, num_msg, sp_humidity_in, temperature_in, zmid, pmid, pint, pblt, tpert, parcel_temp, parcel_qsat, msemax_klev, lcl_temperature, lcl_klev, eql_klev, cape, calc_msemax_klev, prev_msemax_klev, use_input_tq_mx, q_mx, t_mx) bind(C)
  use zm_conv_cape, only : compute_dilute_cape
  use zm_conv_types,  only: zm_const_t, zm_param_t
  use zm_conv_types,  only: zm_param_set_for_testing, zm_const_set_for_testing

  integer(kind=c_int) , value, intent(in) :: pcols, ncol, pver, pverp, num_cin, num_msg
  real(kind=c_real) , intent(in), dimension(pcols, pver) :: sp_humidity_in, temperature_in, zmid, pmid
  real(kind=c_real) , intent(in), dimension(pcols, pverp) :: pint
  integer(kind=c_int) , intent(in), dimension(pcols) :: pblt, prev_msemax_klev
  real(kind=c_real) , intent(in), dimension(pcols) :: tpert
  real(kind=c_real) , intent(out), dimension(pcols, pver) :: parcel_temp
  real(kind=c_real) , intent(inout), dimension(pcols, pver) :: parcel_qsat
  integer(kind=c_int) , intent(inout), dimension(pcols) :: msemax_klev, lcl_klev, eql_klev
  real(kind=c_real) , intent(out), dimension(pcols) :: lcl_temperature
  real(kind=c_real) , intent(inout), dimension(pcols) :: cape, q_mx, t_mx
  logical(kind=c_bool) , value, intent(in) :: calc_msemax_klev, use_input_tq_mx

  type(zm_const_t) :: zm_const ! derived type to hold ZM constants
  type(zm_param_t) :: zm_param ! derived type to hold ZM tunable parameters
  !-----------------------------------------------------------------------------
  call zm_param_set_for_testing(zm_param)
  call zm_const_set_for_testing(zm_const)

  call compute_dilute_cape(pcols, ncol, pver, pverp, num_cin, num_msg, sp_humidity_in, temperature_in, zmid, pmid, pint, pblt, tpert, parcel_temp, parcel_qsat, msemax_klev, lcl_temperature, lcl_klev, eql_klev, cape, zm_const, zm_param, calc_msemax_klev, prev_msemax_klev, use_input_tq_mx, q_mx, t_mx)
end subroutine compute_dilute_cape_bridge_f

subroutine find_mse_max_bridge_f(pcols, ncol, pver, num_msg, msemax_top_k, pergro_active, temperature, zmid, sp_humidity, msemax_klev, mse_max_val) bind(C)
  use zm_conv_cape, only : find_mse_max
  use zm_conv_types,  only: zm_const_t, zm_param_t
  use zm_conv_types,  only: zm_param_set_for_testing, zm_const_set_for_testing

  integer(kind=c_int) , value, intent(in) :: pcols, ncol, pver, num_msg
  integer(kind=c_int) , intent(in), dimension(pcols) :: msemax_top_k
  logical(kind=c_bool) , value, intent(in) :: pergro_active
  real(kind=c_real) , intent(in), dimension(pcols, pver) :: temperature, zmid, sp_humidity
  integer(kind=c_int) , intent(inout), dimension(pcols) :: msemax_klev
  real(kind=c_real) , intent(inout), dimension(pcols) :: mse_max_val

  type(zm_const_t) :: zm_const ! derived type to hold ZM constants
  type(zm_param_t) :: zm_param ! derived type to hold ZM tunable parameters
  !-----------------------------------------------------------------------------
  call zm_param_set_for_testing(zm_param)
  call zm_const_set_for_testing(zm_const)

  call find_mse_max(pcols, ncol, pver, num_msg, msemax_top_k, pergro_active, temperature, zmid, sp_humidity, zm_const, zm_param, msemax_klev, mse_max_val)
end subroutine find_mse_max_bridge_f

subroutine compute_dilute_parcel_bridge_f(pcols, ncol, pver, num_msg, klaunch, pmid, temperature, sp_humidity, tpert, pblt, parcel_temp, parcel_vtemp, parcel_qsat, lcl_pmid, lcl_temperature, lcl_klev) bind(C)
  use zm_conv_cape, only : compute_dilute_parcel
  use zm_conv_types,  only: zm_const_t, zm_param_t
  use zm_conv_types,  only: zm_param_set_for_testing, zm_const_set_for_testing

  integer(kind=c_int) , value, intent(in) :: pcols, ncol, pver, num_msg
  integer(kind=c_int) , intent(in), dimension(pcols) :: klaunch, pblt
  real(kind=c_real) , intent(in), dimension(pcols, pver) :: pmid, temperature, sp_humidity
  real(kind=c_real) , intent(in), dimension(pcols) :: tpert
  real(kind=c_real) , intent(inout), dimension(pcols, pver) :: parcel_temp, parcel_vtemp, parcel_qsat
  real(kind=c_real) , intent(inout), dimension(pcols) :: lcl_pmid, lcl_temperature
  integer(kind=c_int) , intent(inout), dimension(pcols) :: lcl_klev

  type(zm_const_t) :: zm_const ! derived type to hold ZM constants
  type(zm_param_t) :: zm_param ! derived type to hold ZM tunable parameters
  !-----------------------------------------------------------------------------
  call zm_param_set_for_testing(zm_param)
  call zm_const_set_for_testing(zm_const)

  call compute_dilute_parcel(pcols, ncol, pver, num_msg, klaunch, pmid, temperature, sp_humidity, tpert, pblt, zm_const, zm_param, parcel_temp, parcel_vtemp, parcel_qsat, lcl_pmid, lcl_temperature, lcl_klev)
end subroutine compute_dilute_parcel_bridge_f

subroutine compute_cape_from_parcel_bridge_f(pcols, ncol, pver, pverp, num_cin, num_msg, temperature, tv, zmid, sp_humidity, pint, msemax_klev, lcl_pmid, lcl_klev, parcel_qsat, parcel_temp, parcel_vtemp, eql_klev, cape) bind(C)
  use zm_conv_cape, only : compute_cape_from_parcel
  use zm_conv_types,  only: zm_const_t, zm_param_t
  use zm_conv_types,  only: zm_param_set_for_testing, zm_const_set_for_testing

  integer(kind=c_int) , value, intent(in) :: pcols, ncol, pver, pverp, num_cin, num_msg
  real(kind=c_real) , intent(in), dimension(pcols, pver) :: temperature, tv, zmid, sp_humidity
  real(kind=c_real) , intent(in), dimension(pcols, pverp) :: pint
  integer(kind=c_int) , intent(in), dimension(pcols) :: msemax_klev, lcl_klev
  real(kind=c_real) , intent(in), dimension(pcols) :: lcl_pmid
  real(kind=c_real) , intent(inout), dimension(pcols, pver) :: parcel_qsat, parcel_temp, parcel_vtemp
  integer(kind=c_int) , intent(inout), dimension(pcols) :: eql_klev
  real(kind=c_real) , intent(inout), dimension(pcols) :: cape

  type(zm_const_t) :: zm_const ! derived type to hold ZM constants
  type(zm_param_t) :: zm_param ! derived type to hold ZM tunable parameters
  !-----------------------------------------------------------------------------
  call zm_param_set_for_testing(zm_param)
  call zm_const_set_for_testing(zm_const)

  call compute_cape_from_parcel(pcols, ncol, pver, pverp, num_cin, num_msg, temperature, tv, zmid, sp_humidity, pint, msemax_klev, lcl_pmid, lcl_klev, zm_const, zm_param, parcel_qsat, parcel_temp, parcel_vtemp, eql_klev, cape)
end subroutine compute_cape_from_parcel_bridge_f
end module zm_c2f_bridge

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

subroutine compute_cape_from_parcel_bridge_f(pcols, ncol, pver, pverp, num_cin, num_msg, temperature, tv, sp_humidity, pint, msemax_klev, lcl_pmid, lcl_klev, parcel_qsat, parcel_temp, parcel_vtemp, eql_klev, cape) bind(C)
  use zm_conv_cape, only : compute_cape_from_parcel
  use zm_conv_types,  only: zm_const_t, zm_param_t
  use zm_conv_types,  only: zm_param_set_for_testing, zm_const_set_for_testing

  integer(kind=c_int) , value, intent(in) :: pcols, ncol, pver, pverp, num_cin, num_msg
  real(kind=c_real) , intent(in), dimension(pcols, pver) :: temperature, tv, sp_humidity
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

  call compute_cape_from_parcel(pcols, ncol, pver, pverp, num_cin, num_msg, temperature, tv, sp_humidity, pint, msemax_klev, lcl_pmid, lcl_klev, zm_const, zm_param, parcel_qsat, parcel_temp, parcel_vtemp, eql_klev, cape)
end subroutine compute_cape_from_parcel_bridge_f

subroutine zm_conv_mcsp_calculate_shear_bridge_f(pcols, ncol, pver, state_pmid, state_u, state_v, mcsp_shear) bind(C)
  use zm_conv_mcsp, only : zm_conv_mcsp_calculate_shear

  integer(kind=c_int) , value, intent(in) :: pcols, ncol, pver
  real(kind=c_real) , intent(in), dimension(pcols, pver) :: state_pmid, state_u, state_v
  real(kind=c_real) , intent(out), dimension(pcols) :: mcsp_shear

  call zm_conv_mcsp_calculate_shear(pcols, ncol, pver, state_pmid, state_u, state_v, mcsp_shear)
end subroutine zm_conv_mcsp_calculate_shear_bridge_f

subroutine zm_conv_mcsp_tend_bridge_f(pcols, ncol, pver, pverp, ztodt, jctop, state_pmid, state_pint, state_pdel, state_s, state_q, state_u, state_v, ptend_zm_s, ptend_zm_q, ptend_s, ptend_q, ptend_u, ptend_v, mcsp_dt_out, mcsp_dq_out, mcsp_du_out, mcsp_dv_out, mcsp_freq, mcsp_shear, zm_depth) bind(C)
  use zm_conv_mcsp, only : zm_conv_mcsp_tend
  use zm_conv_types,  only: zm_const_t, zm_param_t
  use zm_conv_types,  only: zm_param_set_for_testing, zm_const_set_for_testing

  integer(kind=c_int) , value, intent(in) :: pcols, ncol, pver, pverp
  real(kind=c_real) , value, intent(in) :: ztodt
  integer(kind=c_int) , intent(in), dimension(pcols) :: jctop
  real(kind=c_real) , intent(in), dimension(pcols, pver) :: state_pmid, state_pdel, state_s, state_q, state_u, state_v, ptend_zm_s, ptend_zm_q
  real(kind=c_real) , intent(in), dimension(pcols, pverp) :: state_pint
  real(kind=c_real) , intent(inout), dimension(pcols, pver) :: ptend_s, ptend_q, ptend_u, ptend_v
  real(kind=c_real) , intent(out), dimension(pcols, pver) :: mcsp_dt_out, mcsp_dq_out, mcsp_du_out, mcsp_dv_out
  real(kind=c_real) , intent(out), dimension(pcols) :: mcsp_freq, mcsp_shear, zm_depth

  type(zm_const_t) :: zm_const ! derived type to hold ZM constants
  type(zm_param_t) :: zm_param ! derived type to hold ZM tunable parameters
  !-----------------------------------------------------------------------------
  call zm_param_set_for_testing(zm_param)
  call zm_const_set_for_testing(zm_const)

  call zm_conv_mcsp_tend(pcols, ncol, pver, pverp, ztodt, jctop, zm_const, zm_param, state_pmid, state_pint, state_pdel, state_s, state_q, state_u, state_v, ptend_zm_s, ptend_zm_q, ptend_s, ptend_q, ptend_u, ptend_v, mcsp_dt_out, mcsp_dq_out, mcsp_du_out, mcsp_dv_out, mcsp_freq, mcsp_shear, zm_depth)
end subroutine zm_conv_mcsp_tend_bridge_f

subroutine zm_conv_main_bridge_f(pcols, ncol, pver, pverp, is_first_step, time_step, t_mid, q_mid_in, omega, p_mid_in, p_int_in, p_del_in, geos, z_mid_in, z_int_in, pbl_hgt, tpert, landfrac, t_star, q_star, lengath, gather_index, msemax_klev_g, jctop, jcbot, jt, prec, heat, qtnd, cape, dcape, mcon, pflx, zdu, mflx_up, entr_up, detr_up, mflx_dn, entr_dn, p_del, dsubcld, ql, rliq, rprd, dlf) bind(C)
  use zm_conv, only : zm_conv_main, zm_const, zm_param
  use zm_conv_types,  only: zm_param_set_for_testing, zm_const_set_for_testing
  use zm_aero_type,           only: zm_aero_t
  use zm_microphysics_state,  only: zm_microp_st

  integer(kind=c_int) , value, intent(in) :: pcols, ncol, pver, pverp
  logical(kind=c_bool) , value, intent(in) :: is_first_step
  real(kind=c_real) , value, intent(in) :: time_step
  real(kind=c_real) , intent(in), dimension(pcols, pver) :: t_mid, q_mid_in, omega, p_mid_in, p_del_in, z_mid_in, t_star, q_star
  real(kind=c_real) , intent(in), dimension(pcols, pverp) :: p_int_in, z_int_in
  real(kind=c_real) , intent(in), dimension(pcols) :: geos, pbl_hgt, tpert, landfrac
  integer(kind=c_int) , intent(out) :: lengath
  integer(kind=c_int) , intent(out), dimension(pcols) :: gather_index, msemax_klev_g, jctop, jcbot, jt
  real(kind=c_real) , intent(out), dimension(pcols) :: prec, cape, dcape, dsubcld, rliq
  real(kind=c_real) , intent(out), dimension(pcols, pver) :: heat, qtnd, zdu, mflx_up, entr_up, detr_up, mflx_dn, entr_dn, p_del, ql, rprd, dlf
  real(kind=c_real) , intent(out), dimension(pcols, pverp) :: mcon, pflx

  type(zm_aero_t)    :: aero            ! aerosol object
  type(zm_microp_st) :: microp_st

  call zm_param_set_for_testing(zm_param)
  call zm_const_set_for_testing(zm_const)

  call zm_conv_main(pcols, ncol, pver, pverp, is_first_step, time_step, t_mid, q_mid_in, omega, p_mid_in, p_int_in, p_del_in, geos, z_mid_in, z_int_in, pbl_hgt, tpert, landfrac, t_star, q_star, lengath, gather_index, msemax_klev_g, jctop, jcbot, jt, prec, heat, qtnd, cape, dcape, mcon, pflx, zdu, mflx_up, entr_up, detr_up, mflx_dn, entr_dn, p_del, dsubcld, ql, rliq, rprd, dlf) ! aero, microp_st)
end subroutine zm_conv_main_bridge_f

subroutine zm_conv_evap_bridge_f(pcols, ncol, pver, pverp, time_step, p_mid, p_del, t_mid, q_mid, prdprec, cldfrc, tend_s, tend_q, tend_s_snwprd, tend_s_snwevmlt, prec, snow, ntprprd, ntsnprd, flxprec, flxsnow) bind(C)
  use zm_conv, only : zm_conv_evap, zm_const, zm_param
  use zm_conv_types,  only: zm_param_set_for_testing, zm_const_set_for_testing

  integer(kind=c_int) , value, intent(in) :: pcols, ncol, pver, pverp
  real(kind=c_real) , value, intent(in) :: time_step
  real(kind=c_real) , intent(in), dimension(pcols, pver) :: p_mid, p_del, t_mid, q_mid, prdprec, cldfrc
  real(kind=c_real) , intent(inout), dimension(pcols, pver) :: tend_s, tend_q
  real(kind=c_real) , intent(out), dimension(pcols, pver) :: tend_s_snwprd, tend_s_snwevmlt, ntprprd, ntsnprd
  real(kind=c_real) , intent(inout), dimension(pcols) :: prec
  real(kind=c_real) , intent(out), dimension(pcols) :: snow
  real(kind=c_real) , intent(out), dimension(pcols, pverp) :: flxprec, flxsnow

  call zm_param_set_for_testing(zm_param)
  call zm_const_set_for_testing(zm_const)

  call zm_conv_evap(pcols, ncol, pver, pverp, time_step, p_mid, p_del, t_mid, q_mid, prdprec, cldfrc, tend_s, tend_q, tend_s_snwprd, tend_s_snwevmlt, prec, snow, ntprprd, ntsnprd, flxprec, flxsnow)
end subroutine zm_conv_evap_bridge_f

subroutine zm_calc_fractional_entrainment_bridge_f(pcols, ncol, pver, pverp, msg, jb, jt, j0, z_mid, z_int, dz, h_env, h_env_sat, h_env_min, lambda, lambda_max) bind(C)
  use zm_conv, only : zm_calc_fractional_entrainment, zm_const, zm_param
  use zm_conv_types,  only: zm_param_set_for_testing, zm_const_set_for_testing

  integer(kind=c_int) , value, intent(in) :: pcols, ncol, pver, pverp, msg
  integer(kind=c_int) , intent(in), dimension(pcols) :: jb, jt
  integer(kind=c_int) , intent(inout), dimension(pcols) :: j0
  real(kind=c_real) , intent(in), dimension(pcols, pver) :: z_mid, dz, h_env, h_env_sat
  real(kind=c_real) , intent(in), dimension(pcols, pverp) :: z_int
  real(kind=c_real) , intent(inout), dimension(pcols) :: h_env_min
  real(kind=c_real) , intent(out), dimension(pcols, pver) :: lambda
  real(kind=c_real) , intent(out), dimension(pcols) :: lambda_max

  call zm_param_set_for_testing(zm_param)
  call zm_const_set_for_testing(zm_const)

  call zm_calc_fractional_entrainment(pcols, ncol, pver, pverp, msg, jb, jt, j0, z_mid, z_int, dz, h_env, h_env_sat, h_env_min, lambda, lambda_max)
end subroutine zm_calc_fractional_entrainment_bridge_f

subroutine zm_downdraft_properties_bridge_f(pcols, ncol, pver, pverp, msg, jb, jt, j0, jd, z_int, dz, s_mid, q_mid, h_env, lambda, lambda_max, qsthat, hsthat, gamhat, rprd, mflx_up, mflx_dn, entr_dn, s_dnd, q_dnd, h_dnd, q_dnd_sat, evp, totevp) bind(C)
  use zm_conv, only : zm_downdraft_properties, zm_const, zm_param
  use zm_conv_types,  only: zm_param_set_for_testing, zm_const_set_for_testing

  integer(kind=c_int) , value, intent(in) :: pcols, ncol, pver, pverp, msg
  integer(kind=c_int) , intent(in), dimension(pcols) :: jb, j0
  integer(kind=c_int) , intent(inout), dimension(pcols) :: jt, jd
  real(kind=c_real) , intent(in), dimension(pcols, pverp) :: z_int
  real(kind=c_real) , intent(in), dimension(pcols, pver) :: dz, s_mid, q_mid, h_env, lambda, qsthat, hsthat, gamhat, rprd, mflx_up
  real(kind=c_real) , intent(in), dimension(pcols) :: lambda_max
  real(kind=c_real) , intent(inout), dimension(pcols, pver) :: mflx_dn, entr_dn, s_dnd, q_dnd, h_dnd, q_dnd_sat, evp
  real(kind=c_real) , intent(inout), dimension(pcols) :: totevp

  call zm_param_set_for_testing(zm_param)
  call zm_const_set_for_testing(zm_const)

  call zm_downdraft_properties(pcols, ncol, pver, pverp, msg, jb, jt, j0, jd, z_int, dz, s_mid, q_mid, h_env, lambda, lambda_max, qsthat, hsthat, gamhat, rprd, mflx_up, mflx_dn, entr_dn, s_dnd, q_dnd, h_dnd, q_dnd_sat, evp, totevp)
end subroutine zm_downdraft_properties_bridge_f

subroutine zm_cloud_properties_bridge_f(pcols, ncol, pver, pverp, msg, limcnv, p_mid, z_mid, z_int, t_mid, s_mid, s_int, q_mid, landfrac, tpert_g, jb, lel, jt, jlcl, j0, jd, mflx_up, entr_up, detr_up, mflx_dn, entr_dn, mflx_net, s_upd, q_upd, ql, s_dnd, q_dnd, qst, cu, evp, pflx, rprd) bind(C)
  use zm_conv, only : zm_cloud_properties, zm_const, zm_param
  use zm_conv_types,  only: zm_param_set_for_testing, zm_const_set_for_testing

  integer(kind=c_int) , value, intent(in) :: pcols, ncol, pver, pverp, msg, limcnv
  real(kind=c_real) , intent(in), dimension(pcols, pver) :: p_mid, z_mid, t_mid, s_mid, s_int, q_mid
  real(kind=c_real) , intent(in), dimension(pcols, pverp) :: z_int
  real(kind=c_real) , intent(in), dimension(pcols) :: landfrac, tpert_g
  integer(kind=c_int) , intent(in), dimension(pcols) :: jb, lel
  integer(kind=c_int) , intent(out), dimension(pcols) :: jt, jlcl, j0, jd
  real(kind=c_real) , intent(out), dimension(pcols, pver) :: mflx_up, entr_up, detr_up, mflx_dn, entr_dn, mflx_net, s_upd, q_upd, ql, s_dnd, q_dnd, qst, cu, evp, rprd
  real(kind=c_real) , intent(out), dimension(pcols, pverp) :: pflx

  call zm_param_set_for_testing(zm_param)
  call zm_const_set_for_testing(zm_const)

  call zm_cloud_properties(pcols, ncol, pver, pverp, msg, limcnv, p_mid, z_mid, z_int, t_mid, s_mid, s_int, q_mid, landfrac, tpert_g, jb, lel, jt, jlcl, j0, jd, mflx_up, entr_up, detr_up, mflx_dn, entr_dn, mflx_net, s_upd, q_upd, ql, s_dnd, q_dnd, qst, cu, evp, pflx, rprd)
end subroutine zm_cloud_properties_bridge_f

subroutine zm_closure_bridge_f(pcols, ncol, pver, pverp, msg, cape_threshold_in, lcl, lel, jt, mx, dsubcld, z_mid, z_int, p_mid, p_del, t_mid, s_mid, q_mid, qs, ql, s_int, q_int, t_pcl_lcl, t_pcl, q_pcl_sat, s_upd, q_upd, mflx_net, detr_up, mflx_up, mflx_dn, q_dnd, s_dnd, cape, cld_base_mass_flux) bind(C)
  use zm_conv, only : zm_closure, zm_const, zm_param
  use zm_conv_types,  only: zm_param_set_for_testing, zm_const_set_for_testing

  integer(kind=c_int) , value, intent(in) :: pcols, ncol, pver, pverp, msg
  real(kind=c_real) , value, intent(in) :: cape_threshold_in
  integer(kind=c_int) , intent(in), dimension(pcols) :: lcl, lel, jt, mx
  real(kind=c_real) , intent(in), dimension(pcols) :: dsubcld, t_pcl_lcl, cape
  real(kind=c_real) , intent(in), dimension(pcols, pver) :: z_mid, p_mid, p_del, t_mid, s_mid, q_mid, qs, ql, s_int, q_int, t_pcl, q_pcl_sat, s_upd, q_upd, mflx_net, detr_up, mflx_up, mflx_dn, q_dnd, s_dnd
  real(kind=c_real) , intent(in), dimension(pcols, pverp) :: z_int
  real(kind=c_real) , intent(out), dimension(pcols) :: cld_base_mass_flux

  call zm_const_set_for_testing(zm_const)
  call zm_param_set_for_testing(zm_param)

  call zm_closure(pcols, ncol, pver, pverp, msg, cape_threshold_in, lcl, lel, jt, mx, dsubcld, z_mid, z_int, p_mid, p_del, t_mid, s_mid, q_mid, qs, ql, s_int, q_int, t_pcl_lcl, t_pcl, q_pcl_sat, s_upd, q_upd, mflx_net, detr_up, mflx_up, mflx_dn, q_dnd, s_dnd, cape, cld_base_mass_flux)
end subroutine zm_closure_bridge_f

subroutine zm_calc_output_tend_bridge_f(pcols, ncol, pver, pverp, msg, jt, mx, dsubcld, p_del, s_int, q_int, s_upd, q_upd, mflx_up, detr_up, mflx_dn, s_dnd, q_dnd, ql, evp, cu, dsdt, dqdt, dl) bind(C)
  use zm_conv, only : zm_calc_output_tend, zm_const, zm_param
  use zm_conv_types,  only: zm_param_set_for_testing, zm_const_set_for_testing

  integer(kind=c_int) , value, intent(in) :: pcols, ncol, pver, pverp, msg
  integer(kind=c_int) , intent(in), dimension(pcols) :: jt, mx
  real(kind=c_real) , intent(in), dimension(pcols) :: dsubcld
  real(kind=c_real) , intent(in), dimension(pcols, pver) :: p_del, s_int, q_int, s_upd, q_upd, mflx_up, detr_up, mflx_dn, s_dnd, q_dnd, ql, evp, cu
  real(kind=c_real) , intent(out), dimension(pcols, pver) :: dsdt, dqdt, dl

  call zm_param_set_for_testing(zm_param)
  call zm_const_set_for_testing(zm_const)

  call zm_calc_output_tend(pcols, ncol, pver, pverp, msg, jt, mx, dsubcld, p_del, s_int, q_int, s_upd, q_upd, mflx_up, detr_up, mflx_dn, s_dnd, q_dnd, ql, evp, cu, dsdt, dqdt, dl)
end subroutine zm_calc_output_tend_bridge_f

end module zm_c2f_bridge

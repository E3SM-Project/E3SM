module dp_iso_c
  use iso_c_binding
  implicit none

#include "scream_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
#else
# define c_real c_float
#endif

!
! This file contains bridges from scream c++ to DP fortran.
!

contains
  subroutine advance_iop_forcing_c(plev, pcnst, scm_dt, ps_in, u_in, v_in, t_in, q_in, t_phys_frc, u_update, v_update, t_update, q_update) bind(C)
    !use dp, only : advance_iop_forcing

    integer(kind=c_int) , value, intent(in) :: plev, pcnst
    real(kind=c_real) , value, intent(in) :: scm_dt, ps_in
    real(kind=c_real) , intent(in), dimension(plev) :: u_in, v_in, t_in, t_phys_frc
    real(kind=c_real) , intent(in), dimension(plev, pcnst) :: q_in
    real(kind=c_real) , intent(out), dimension(plev) :: u_update, v_update, t_update
    real(kind=c_real) , intent(out), dimension(plev, pcnst) :: q_update

    !call advance_iop_forcing(plev, pcnst, scm_dt, ps_in, u_in, v_in, t_in, q_in, t_phys_frc, u_update, v_update, t_update, q_update)
  end subroutine advance_iop_forcing_c
  subroutine advance_iop_nudging_c(plev, scm_dt, ps_in, t_in, q_in, t_update, q_update, relaxt, relaxq) bind(C)
    !use dp, only : advance_iop_nudging

    integer(kind=c_int) , value, intent(in) :: plev
    real(kind=c_real) , value, intent(in) :: scm_dt, ps_in
    real(kind=c_real) , intent(in), dimension(plev) :: t_in, q_in
    real(kind=c_real) , intent(out), dimension(plev) :: t_update, q_update, relaxt, relaxq

    !call advance_iop_nudging(plev, scm_dt, ps_in, t_in, q_in, t_update, q_update, relaxt, relaxq)
  end subroutine advance_iop_nudging_c
  subroutine advance_iop_subsidence_c(plev, pcnst, scm_dt, ps_in, u_in, v_in, t_in, q_in, u_update, v_update, t_update, q_update) bind(C)
    !use dp, only : advance_iop_subsidence

    integer(kind=c_int) , value, intent(in) :: plev, pcnst
    real(kind=c_real) , value, intent(in) :: scm_dt, ps_in
    real(kind=c_real) , intent(in), dimension(plev) :: u_in, v_in, t_in
    real(kind=c_real) , intent(in), dimension(plev, pcnst) :: q_in
    real(kind=c_real) , intent(out), dimension(plev) :: u_update, v_update, t_update
    real(kind=c_real) , intent(out), dimension(plev, pcnst) :: q_update

    !call advance_iop_subsidence(plev, pcnst, scm_dt, ps_in, u_in, v_in, t_in, q_in, u_update, v_update, t_update, q_update)
  end subroutine advance_iop_subsidence_c
  subroutine iop_setinitial_c(nelemd, elem) bind(C)
    !use dp, only : iop_setinitial

    integer(kind=c_int) , value, intent(in) :: nelemd
    type(c_ptr) , intent(inout), dimension(nelemd) :: elem

    !call iop_setinitial(nelemd, elem)
  end subroutine iop_setinitial_c
  subroutine iop_broadcast_c() bind(C)
    !use dp, only : iop_broadcast

    !call iop_broadcast()
  end subroutine iop_broadcast_c
  subroutine apply_iop_forcing_c(nelemd, elem, hvcoord, hybrid, tl, n, t_before_advance, nets, nete) bind(C)
    !use dp, only : apply_iop_forcing

    integer(kind=c_int) , value, intent(in) :: nelemd, n, nets, nete
    type(c_ptr) , intent(inout), dimension(nelemd) :: elem
    type(c_ptr) , intent(inout) :: hvcoord
    type(c_ptr) , intent(in) :: hybrid
    type(c_ptr) , intent(in) :: tl
    logical(kind=c_bool) , value, intent(in) :: t_before_advance

    !call apply_iop_forcing(nelemd, elem, hvcoord, hybrid, tl, n, t_before_advance, nets, nete)
  end subroutine apply_iop_forcing_c
  subroutine iop_domain_relaxation_c(nelemd, np, nlev, elem, hvcoord, hybrid, t1, dp, nelemd_todo, np_todo, dt) bind(C)
    !use dp, only : iop_domain_relaxation

    integer(kind=c_int) , value, intent(in) :: nelemd, np, nlev, t1, nelemd_todo, np_todo
    type(c_ptr) , intent(inout), dimension(nelemd) :: elem
    type(c_ptr) , intent(in) :: hvcoord
    type(c_ptr) , intent(in) :: hybrid
    real(kind=c_real) , intent(inout), dimension(np, np, nlev) :: dp
    real(kind=c_real) , value, intent(in) :: dt

    !call iop_domain_relaxation(nelemd, np, nlev, elem, hvcoord, hybrid, t1, dp, nelemd_todo, np_todo, dt)
  end subroutine iop_domain_relaxation_c
  subroutine crm_resolved_turb_c(nelemd, elem, hvcoord, hybrid, t1, nelemd_todo, np_todo) bind(C)
    !use dp, only : crm_resolved_turb

    integer(kind=c_int) , value, intent(in) :: nelemd, t1, nelemd_todo, np_todo
    type(c_ptr) , intent(inout), dimension(nelemd) :: elem
    type(c_ptr) , intent(in) :: hvcoord
    type(c_ptr) , intent(in) :: hybrid

    !call crm_resolved_turb(nelemd, elem, hvcoord, hybrid, t1, nelemd_todo, np_todo)
  end subroutine crm_resolved_turb_c
  subroutine iop_default_opts_c(scmlat_out, scmlon_out, iopfile_out,       &
       single_column_out, scm_iop_srf_prop_out, iop_nudge_tq_out,          &
       iop_nudge_uv_out, iop_nudge_tq_low_out, iop_nudge_tq_high_out,      &
       iop_nudge_tscale_out, scm_observed_aero_out, iop_dosubsidence_out,  &
       scm_multcols_out, dp_crm_out, iop_perturb_high_out, precip_off_out, &
       scm_zero_non_iop_tracers_out) bind(C)
    !use dp, only : iop_default_opts

    real(kind=c_real) , intent(out) :: scmlat_out, scmlon_out, iop_nudge_tq_low_out, iop_nudge_tq_high_out, iop_nudge_tscale_out, iop_perturb_high_out
    type(c_ptr) , intent(out) :: iopfile_out
    logical(kind=c_bool) , intent(out) :: single_column_out, &
         scm_iop_srf_prop_out, iop_nudge_tq_out, iop_nudge_uv_out, &
         scm_observed_aero_out, iop_dosubsidence_out, scm_multcols_out, &
         dp_crm_out, precip_off_out, scm_zero_non_iop_tracers_out

    !call iop_default_opts(scmlat_out, scmlon_out, iopfile_out, single_column_out, scm_iop_srf_prop_out, iop_nudge_tq_out, iop_nudge_uv_out, iop_nudge_tq_low_out, iop_nudge_tq_high_out, iop_nudge_tscale_out, scm_observed_aero_out, iop_dosubsidence_out, scm_multcols_out, dp_crm_out, iop_perturb_high_out, precip_off_out, scm_zero_non_iop_tracers_out)
  end subroutine iop_default_opts_c
  subroutine iop_setopts_c(scmlat_in, scmlon_in, iopfile_in, &
       single_column_in, scm_iop_srf_prop_in, iop_nudge_tq_in, &
       iop_nudge_uv_in, iop_nudge_tq_low_in, iop_nudge_tq_high_in, &
       iop_nudge_tscale_in, scm_observed_aero_in, iop_dosubsidence_in, &
       scm_multcols_in, dp_crm_in, iop_perturb_high_in, precip_off_in, &
       scm_zero_non_iop_tracers_in) bind(C)
    !use dp, only : iop_setopts

    real(kind=c_real) , value, intent(in) :: scmlat_in, scmlon_in, iop_nudge_tq_low_in, iop_nudge_tq_high_in, iop_nudge_tscale_in, iop_perturb_high_in
    type(c_ptr) , intent(in) :: iopfile_in
    logical(kind=c_bool) , value, intent(in) :: single_column_in, &
         scm_iop_srf_prop_in, iop_nudge_tq_in, iop_nudge_uv_in, &
         scm_observed_aero_in, iop_dosubsidence_in, scm_multcols_in, &
         dp_crm_in, precip_off_in, scm_zero_non_iop_tracers_in

    !call iop_setopts(scmlat_in, scmlon_in, iopfile_in, single_column_in, scm_iop_srf_prop_in, iop_nudge_tq_in, iop_nudge_uv_in, iop_nudge_tq_low_in, iop_nudge_tq_high_in, iop_nudge_tscale_in, scm_observed_aero_in, iop_dosubsidence_in, scm_multcols_in, dp_crm_in, iop_perturb_high_in, precip_off_in, scm_zero_non_iop_tracers_in)
  end subroutine iop_setopts_c
  subroutine setiopupdate_init_c() bind(C)
    !use dp, only : setiopupdate_init

    ! call setiopupdate_init()
  end subroutine setiopupdate_init_c
  subroutine setiopupdate_c() bind(C)
    !use dp, only : setiopupdate

    !call setiopupdate()
  end subroutine setiopupdate_c
  subroutine readiopdata_c(plev, iop_update_phase1, hyam, hybm) bind(C)
    !use dp, only : readiopdata

    integer(kind=c_int) , value, intent(in) :: plev
    logical(kind=c_bool) , value, intent(in) :: iop_update_phase1
    real(kind=c_real) , intent(in), dimension(plev) :: hyam, hybm

    !call readiopdata(plev, iop_update_phase1, hyam, hybm)
  end subroutine readiopdata_c
  subroutine iop_intht_c() bind(C)
    !use dp, only : iop_intht

    !call iop_intht()
  end subroutine iop_intht_c
end module dp_iso_c

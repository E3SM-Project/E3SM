module dp_iso_f
  use iso_c_binding
  implicit none

#include "scream_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
#else
# define c_real c_float
#endif

!
! This file contains bridges from DP fortran to scream c++.
!

interface

  subroutine advance_iop_forcing_f(plev, pcnst, scm_dt, ps_in, u_in, v_in, t_in, q_in, t_phys_frc, u_update, v_update, t_update, q_update) bind(C)
    use iso_c_binding

    integer(kind=c_int) , value, intent(in) :: plev, pcnst
    real(kind=c_real) , value, intent(in) :: scm_dt, ps_in
    real(kind=c_real) , intent(in), dimension(plev) :: u_in, v_in, t_in, t_phys_frc
    real(kind=c_real) , intent(in), dimension(plev, pcnst) :: q_in
    real(kind=c_real) , intent(out), dimension(plev) :: u_update, v_update, t_update
    real(kind=c_real) , intent(out), dimension(plev, pcnst) :: q_update
  end subroutine advance_iop_forcing_f
  subroutine advance_iop_nudging_f(plev, scm_dt, ps_in, t_in, q_in, t_update, q_update, relaxt, relaxq) bind(C)
    use iso_c_binding

    integer(kind=c_int) , value, intent(in) :: plev
    real(kind=c_real) , value, intent(in) :: scm_dt, ps_in
    real(kind=c_real) , intent(in), dimension(plev) :: t_in, q_in
    real(kind=c_real) , intent(out), dimension(plev) :: t_update, q_update, relaxt, relaxq
  end subroutine advance_iop_nudging_f
  subroutine advance_iop_subsidence_f(plev, pcnst, scm_dt, ps_in, u_in, v_in, t_in, q_in, u_update, v_update, t_update, q_update) bind(C)
    use iso_c_binding

    integer(kind=c_int) , value, intent(in) :: plev, pcnst
    real(kind=c_real) , value, intent(in) :: scm_dt, ps_in
    real(kind=c_real) , intent(in), dimension(plev) :: u_in, v_in, t_in
    real(kind=c_real) , intent(in), dimension(plev, pcnst) :: q_in
    real(kind=c_real) , intent(out), dimension(plev) :: u_update, v_update, t_update
    real(kind=c_real) , intent(out), dimension(plev, pcnst) :: q_update
  end subroutine advance_iop_subsidence_f
  subroutine iop_setinitial_f(nelemd, elem) bind(C)
    use iso_c_binding

    integer(kind=c_int) , value, intent(in) :: nelemd
    type(c_ptr) , intent(inout), dimension(nelemd) :: elem
  end subroutine iop_setinitial_f
  subroutine iop_broadcast_f() bind(C)
    use iso_c_binding

    
  end subroutine iop_broadcast_f
  subroutine apply_iop_forcing_f(nelemd, elem, hvcoord, hybrid, tl, n, t_before_advance, nets, nete) bind(C)
    use iso_c_binding

    integer(kind=c_int) , value, intent(in) :: nelemd, n, nets, nete
    type(c_ptr) , intent(inout), dimension(nelemd) :: elem
    type(c_ptr) , intent(inout) :: hvcoord
    type(c_ptr) , intent(in) :: hybrid
    type(c_ptr) , intent(in) :: tl
    logical(kind=c_bool) , value, intent(in) :: t_before_advance
  end subroutine apply_iop_forcing_f
  subroutine iop_domain_relaxation_f(nelemd, np, nlev, elem, hvcoord, hybrid, t1, dp, nelemd_todo, np_todo, dt) bind(C)
    use iso_c_binding

    integer(kind=c_int) , value, intent(in) :: nelemd, np, nlev, t1, nelemd_todo, np_todo
    type(c_ptr) , intent(inout), dimension(nelemd) :: elem
    type(c_ptr) , intent(in) :: hvcoord
    type(c_ptr) , intent(in) :: hybrid
    real(kind=c_real) , intent(inout), dimension(np, np, nlev) :: dp
    real(kind=c_real) , value, intent(in) :: dt
  end subroutine iop_domain_relaxation_f
  subroutine crm_resolved_turb_f(nelemd, elem, hvcoord, hybrid, t1, nelemd_todo, np_todo) bind(C)
    use iso_c_binding

    integer(kind=c_int) , value, intent(in) :: nelemd, t1, nelemd_todo, np_todo
    type(c_ptr) , intent(inout), dimension(nelemd) :: elem
    type(c_ptr) , intent(in) :: hvcoord
    type(c_ptr) , intent(in) :: hybrid
  end subroutine crm_resolved_turb_f
end interface

end module dp_iso_f

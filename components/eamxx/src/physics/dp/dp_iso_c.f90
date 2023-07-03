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
end module dp_iso_c

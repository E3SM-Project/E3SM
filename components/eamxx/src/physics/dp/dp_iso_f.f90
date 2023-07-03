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
end interface

end module dp_iso_f

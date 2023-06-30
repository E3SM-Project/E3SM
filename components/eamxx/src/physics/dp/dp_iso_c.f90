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
  subroutine advance_iop_forcing_c(scm_dt, ps_in, u_in, v_in, t_in, q_in, t_phys_frc, u_update, v_update, t_update, q_update) bind(C)
    !use apply_iop_forcing, only : advance_iop_forcing

    real(kind=c_real) , value, intent(in) :: scm_dt, ps_in
    real(kind=c_real) , intent(in), dimension(plev) :: u_in, v_in, t_in, t_phys_frc
    real(kind=c_real) , intent(in), dimension(plev, pcnst) :: q_in
    real(kind=c_real) , intent(out), dimension(plev) :: u_update, v_update, t_update
    real(kind=c_real) , intent(out), dimension(plev, pcnst) :: q_update

    !call advance_iop_forcing(scm_dt, ps_in, u_in, v_in, t_in, q_in, t_phys_frc, u_update, v_update, t_update, q_update)
  end subroutine advance_iop_forcing_c
end module dp_iso_c


module tms_iso_c
  use iso_c_binding
  implicit none

#include "eamxx_config.f"

!
! This file contains bridges from scream c++ to tms fortran.
!

contains
  subroutine init_tms_c(orocnst, z0fac, karman, gravit, rair) bind(c)
    use trb_mtn_stress, only: init_tms

    real(kind=c_double), value, intent(in)  :: orocnst, z0fac, karman, gravit, rair
    character(len=128) :: errstring

    integer,  parameter :: r8 = selected_real_kind(12) ! 8 byte real

    call init_tms(r8, orocnst, z0fac, karman, gravit, rair, errstring)
  end subroutine init_tms_c

  subroutine compute_tms_c(ncols, nlevs, u_wind, v_wind, t_mid, p_mid, exner, &
                           zm, sgh, landfrac, ksrf, taux, tauy) bind(c)
    use trb_mtn_stress, only: compute_tms

    integer(kind=c_int), value, intent(in) :: ncols, nlevs
    real(kind=c_double) , intent(in), dimension(ncols, nlevs) :: u_wind,v_wind,t_mid,p_mid,exner,zm
    real(kind=c_double) , intent(in), dimension(ncols) :: sgh,landfrac
    real(kind=c_double) , intent(out), dimension(ncols) :: ksrf, taux, tauy

    call compute_tms(ncols, nlevs, ncols, u_wind, v_wind, t_mid, p_mid, exner, zm, sgh, ksrf, taux, tauy, landfrac)
  end subroutine compute_tms_c

end module tms_iso_c

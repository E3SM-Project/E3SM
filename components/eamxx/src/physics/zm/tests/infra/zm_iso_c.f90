module zm_iso_c
  use iso_c_binding
  implicit none

#include "eamxx_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
#else
# define c_real c_float
#endif

! This file contains bridges from scream c++ to gw fortran.

!===================================================================================================
contains
!===================================================================================================

! subroutine zm_init_c( pver_in ) bind(C)
!   use zm_common, only : zm_common_init
!   !-----------------------------------------------------------------------------
!   ! Interface Arguments
!   integer(kind=c_int), intent(in), value                 :: ncol_in
!   integer(kind=c_int), intent(in), value                 :: pver_in
!   integer(kind=c_int), intent(in), dimension(ncol)       :: ?
!   real(kind=c_real),   intent(in), dimension(ncol, pver) :: ?
!   !-----------------------------------------------------------------------------
!   ! Local Variables
!   character(len=128) :: errstring
!   !-----------------------------------------------------------------------------
  
!   !-----------------------------------------------------------------------------
! end subroutine zm_init_c

!===================================================================================================

! subroutine zm_compute_dilute_cape_c( ncol, pver ) bind(C)
!   use zm_conv_cape, only : compute_dilute_cape
!   !-----------------------------------------------------------------------------
!   ! Interface Arguments
!   integer(kind=c_int), intent(in), value :: ncol
!   integer(kind=c_int), intent(in), value :: pver
!   ! logical(kind=c_bool),intent(in), value :: do_taper
!   ! real(kind=c_real),   intent(in), value :: dt
!   ! integer(kind=c_int), intent(in), dimension(ncol) :: tend_level
!   ! real(kind=c_real),   intent(in), dimension(ncol) :: lat
!   ! real(kind=c_real),   intent(in), dimension(ncol, pver) :: dpm
!   !-----------------------------------------------------------------------------
!   ! Local Variables
!   ! ???
!   !-----------------------------------------------------------------------------
!   call compute_dilute_cape( pcols, ncol, pver, pverp, &
!                             num_cin, num_msg, &
!                             sp_humidity_in, temperature_in, &
!                             zmid, pmid, pint, pblt, tpert, &
!                             parcel_temp, parcel_qsat, msemax_klev, &
!                             lcl_temperature, lcl_klev, &
!                             eql_klev, cape, &
!                             zm_const, zm_param, &
!                             iclosure, dcapemx, &
!                             use_input_tq_mx, q_mx, t_mx )
!   !-----------------------------------------------------------------------------
! end subroutine zm_compute_dilute_cape_c

!===================================================================================================

! subroutine zm_find_mse_max_c( ncol, pver ) bind(C)
!   use zm_conv_cape,   only: find_mse_max
!   use zm_conv_types,  only: zm_param_set_for_testing, zm_const_set_for_testing
!   !-----------------------------------------------------------------------------
!   ! Interface Arguments
!   integer(kind=c_int), intent(in), value                 :: ncol
!   integer(kind=c_int), intent(in), value                 :: pver
!   real(kind=c_real),   intent(in), dimension(ncol, pver) :: zmid_in
!   real(kind=c_real),   intent(in), dimension(ncol, pver) :: temperature_in
!   real(kind=c_real),   intent(in), dimension(ncol, pver) :: sp_humidity_in
!   !-----------------------------------------------------------------------------
!   ! Local Variables
!   integer(kind=c_int)                       :: pcols           ! number of atmospheric columns (max)
!   real(kind=c_real)                         :: num_msg         ! number of missing moisture levels at the top of model
!   integer(kind=c_int), dimension(ncol)      :: msemax_top_k    ! upper limit index of max MSE search
!   logical(kind=c_bool)                      :: pergro_active   ! flag for perturbation growth test (pergro)
!   real(kind=c_real),   dimension(ncol,pver) :: zmid            ! height/altitude at mid-levels
!   real(kind=c_real),   dimension(ncol,pver) :: temperature     ! environement temperature
!   real(kind=c_real),   dimension(ncol,pver) :: sp_humidity     ! specific humidity
!   type(zm_const_t)                          :: zm_const        ! derived type to hold ZM constants
!   type(zm_param_t)                          :: zm_param        ! derived type to hold ZM tunable parameters
!   integer(kind=c_int), dimension(ncol)      :: msemax_klev     ! index of max MSE at parcel launch level
!   real(kind=c_real),   dimension(ncol)      :: mse_max_val     ! value of max MSE at parcel launch level
!   !-----------------------------------------------------------------------------
!   ! initialize input data
!   pcols           = ncol
!   num_msg         = 1
!   pergro_active   = .false.
!   msemax_top_k(:) = 1 ! allow all levels to be searched

!   do i = 1,ncol
!     do k = 1,pver
!       zmid       (i,k) = 
!       temperature(i,k) = 
!       sp_humidity(i,k) = 
!     end do
!     msemax_klev(i) = -999
!     mse_max_val(i) = -999
!   end do
!   !-----------------------------------------------------------------------------
!   call zm_param_set_for_testing(zm_param)
!   call zm_const_set_for_testing(zm_const)
!   !-----------------------------------------------------------------------------
!   call find_mse_max( pcols, ncol, pver, num_msg, &
!                      msemax_top_k, pergro_active, &
!                      temperature, zmid, sp_humidity, &
!                      zm_const, zm_param, &
!                      msemax_klev, mse_max_val )
!   !-----------------------------------------------------------------------------
! end subroutine zm_find_mse_max_c


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

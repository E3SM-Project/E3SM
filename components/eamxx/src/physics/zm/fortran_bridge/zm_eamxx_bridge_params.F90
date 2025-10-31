module zm_eamxx_bridge_params
  !-----------------------------------------------------------------------------
  ! Purpose: 
  !-----------------------------------------------------------------------------
  use iso_c_binding
  !-----------------------------------------------------------------------------
#include "eamxx_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
#else
  ZM BRIDGE DOES NOT SUPPORT SINGLE PRECISION
# define c_real c_float
#endif
  !-----------------------------------------------------------------------------
  ! public variables
  integer, public, parameter :: r8 = c_real
  integer, public, parameter :: btype = c_bool
  logical, public :: masterproc

  integer, public :: pcols
  integer, public :: pver
  integer, public :: pverp
  integer, public :: top_lev

end module zm_eamxx_bridge_params

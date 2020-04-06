! physics_utils: common elements for SCREAM physics models.
module physics_utils

#ifdef SCREAM_CONFIG_IS_CMAKE
  use iso_c_binding, only: c_double, c_float, c_bool
#else
  use shr_kind_mod,   only: rtype=>shr_kind_r8, itype=>shr_kind_i8
#endif

  implicit none
  private
  save

#ifdef SCREAM_CONFIG_IS_CMAKE
#include "scream_config.f"

  integer,parameter,public :: rtype8 = c_double ! 8 byte real, compatible with c type double
  integer,parameter,public :: btype  = c_bool ! boolean type, compatible with c
  integer,parameter,public :: itype = selected_int_kind (13) ! 8 byte integer

#  ifdef SCREAM_DOUBLE_PRECISION
  integer,parameter,public :: rtype = c_double ! 8 byte real, compatible with c type double
#  else
  integer,parameter,public :: rtype = c_float ! 4 byte real, compatible with c type float
#  endif

#else
  integer,parameter,public :: btype = kind(.true.) ! native logical
  public :: rtype
  integer,parameter,public :: rtype8 = selected_real_kind(15, 307) ! 8 byte real, compatible with c type double

#endif

    contains

end module physics_utils

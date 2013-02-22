#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module kinds
#ifdef CAM
      use shr_kind_mod, only : SHR_KIND_I4, SHR_KIND_R8, SHR_KIND_I8, SHR_KIND_CL
      use cam_logfile, only : iulog ! _EXTERNAL
#endif
implicit none
private
!
!  most floating point variables should be of type real_kind = real*8
!  For higher precision, we also have quad_kind = real*16, but this
!  is only supported on IBM systems
! 
#if defined(CAM)  
  integer (kind=4), public, parameter::  &
  int_kind     = SHR_KIND_I4,            &
  log_kind     = kind(.true.),           &
  long_kind    = SHR_KIND_I8,            &
  real_kind    = SHR_KIND_R8,            &
  longdouble_kind    = 8
  public :: shr_kind_cl, iulog

#else 
  integer (kind=4), public, parameter::  &
  int_kind     = 4,                      &
  long_kind    = 8,                      &
  log_kind     = 4,                      &
  real_kind    = 8,                      &
  iulog        = 6,                      & ! stderr file handle
#if HOMME_NO_QUAD_PREC
  longdouble_kind    = 8
#else 
  longdouble_kind    = 16
#endif
#endif

end module kinds


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module kinds

! EAM uses internal iulog. SCREAM/HOMME use one defined here (defaulting to stdout)
#ifdef CAM
  use cam_logfile, only : iulog ! _EXTERNAL
#else
  use iso_fortran_env, only: output_unit
#endif

! EAM/SCREAM builds can use kinds from shr_kind_mod
#if defined(CAM) || defined(SCREAM)
  use shr_kind_mod, only : SHR_KIND_I4, SHR_KIND_R8, SHR_KIND_I8, SHR_KIND_CL
#endif

implicit none
private
!
!  most floating point variables should be of type real_kind = real*8
!  For higher precision, we also have quad_kind = real*16, but this
!  is only supported on IBM systems
! 
#if defined(CAM) || defined(SCREAM)
  integer (kind=4), public, parameter::  &
  int_kind     = SHR_KIND_I4,            &
  log_kind     = kind(.true.),           &
  long_kind    = SHR_KIND_I8,            &
  real_kind    = SHR_KIND_R8
#else 
  ! STANDALONE HOMME
  integer (kind=4), public, parameter::  &
  int_kind     = 4,                      &
  long_kind    = 8,                      &
  log_kind     = 4,                      &
  real_kind    = 8
#endif

! EAM uses iulog from cam_logfile, SCREAM/Homme will declare iulog in this file
#if defined(CAM)
  public :: iulog
#else
  ! default to output unit to stdout, but can be changed at runtime
  integer (kind=4), public :: iulog = output_unit
#endif

  ! QUAD precision can be enabled via CMake var (or compiler macro for EAM).
  ! Currently, EAM does not use quad precision, while HommeStandalone/SCREAM can.
  integer (kind=4), public, parameter::  &
#ifdef HOMME_QUAD_PREC
  longdouble_kind    = 16
#else 
  longdouble_kind    = 8
#endif

end module kinds


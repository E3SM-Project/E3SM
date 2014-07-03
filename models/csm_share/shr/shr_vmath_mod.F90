!===============================================================================
! SVN $Id: shr_vmath_mod.F90 6752 2007-10-04 21:02:15Z jwolfe $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_140509/shr/shr_vmath_mod.F90 $
!===============================================================================
! PURPOSE: 
!   provides a uniform, platform-independent API for vector math functions
!===============================================================================

module shr_vmath_mod

   !----------------------------------------------------------------------------
   ! routines that evaluate various math functions for vector arguments
   ! intended to provide platform independent access to vendor optimized code
   !----------------------------------------------------------------------------

   use shr_kind_mod
   use shr_log_mod, only: s_loglev  => shr_log_Level
   use shr_log_mod, only: s_logunit => shr_log_Unit

   implicit none

   private
   public :: shr_vmath_sqrt, &
      shr_vmath_exp, shr_vmath_log, &
      shr_vmath_sin, shr_vmath_cos, &
      shr_vmath_rsqrt, shr_vmath_div

   contains

!===============================================================================

subroutine shr_vmath_sqrt(X, Y, n)

   !----- arguments ---
   integer(SHR_KIND_IN),intent(in)  ::   n  ! vector length
   real   (SHR_KIND_R8),intent(in)  :: X(n) ! input vector argument
   real   (SHR_KIND_R8),intent(out) :: Y(n) ! output vector argument

!-------------------------------------------------------------------------------
! PURPOSE: sqrt for vector arguments, optimized on different platforms
!-------------------------------------------------------------------------------

#if (defined NO_SHR_VMATH)
   Y = sqrt(X)
#else

#if (defined AIX)
   call vsqrt(Y, X, n)
#endif

#if (defined IRIX64)
   call shr_vmath_fwrap_vsqrt(X, Y, n)
#endif

#if (defined OSF1)
   call vsqrt(X, 1, Y, 1, n)
#endif

#if (!defined AIX && !defined IRIX64 && !defined OSF1)
   Y = sqrt(X)
#endif
#endif

end subroutine shr_vmath_sqrt

!===============================================================================

subroutine shr_vmath_rsqrt(X, Y, n)

   !----- arguments ---
   integer(SHR_KIND_IN),intent(in)  ::   n  ! vector length
   real   (SHR_KIND_R8),intent(in)  :: X(n) ! input vector argument
   real   (SHR_KIND_R8),intent(out) :: Y(n) ! output vector argument

!-------------------------------------------------------------------------------
! PURPOSE: sqrt for vector arguments, optimized on different platforms
!-------------------------------------------------------------------------------

#if (defined NO_SHR_VMATH)
   Y = 1.0_SHR_KIND_R8/sqrt(X)
#else

#if (defined AIX)
   call vrsqrt(Y, X, n)
#endif

#if (!defined AIX)
   Y = 1.0_SHR_KIND_R8/sqrt(X)
#endif
#endif

end subroutine shr_vmath_rsqrt

!===============================================================================

subroutine shr_vmath_exp(X, Y, n)

   !----- arguments ---
   integer(SHR_KIND_IN),intent(in)  ::   n  ! vector length
   real   (SHR_KIND_R8),intent(in)  :: X(n) ! input vector argument
   real   (SHR_KIND_R8),intent(out) :: Y(n) ! output vector argument

!-------------------------------------------------------------------------------
! PURPOSE: exp for vector arguments, optimized on different platforms
!-------------------------------------------------------------------------------

#if (defined NO_SHR_VMATH)
   Y = exp(X)
#else

#if (defined AIX)
   call vexp(Y, X, n)
#endif

#if (defined IRIX64)
   call shr_vmath_fwrap_vexp(X, Y, n)
#endif

#if (defined OSF1)
   call vexp(X, 1, Y, 1, n)
#endif

#if (!defined AIX && !defined IRIX64 && !defined OSF1)
   Y = exp(X)
#endif
#endif

end subroutine shr_vmath_exp

!===============================================================================

subroutine shr_vmath_div(X, Y, Z, n)
   !----- arguments ---
   integer(SHR_KIND_IN),intent(in)  ::   n  ! vector length
   real   (SHR_KIND_R8),intent(in)  :: X(n) ! input vector argument
   real   (SHR_KIND_R8),intent(in)  :: Y(n) ! input vector argument
   real   (SHR_KIND_R8),intent(out) :: Z(n) ! output vector argument

#if (defined NO_SHR_VMATH)
   integer :: i
   do i=1,n
      Z(i) = X(i)/Y(i)
   enddo
#else
#if (defined AIX)
   call vdiv(Z,X,Y,n)
#else
   integer :: i
   do i=1,n
      Z(i) = X(i)/Y(i)
   enddo
#endif
#endif
   return
 end subroutine shr_vmath_div

!===============================================================================

subroutine shr_vmath_log(X, Y, n)

   !----- arguments ---
   integer(SHR_KIND_IN),intent(in)  ::   n  ! vector length
   real   (SHR_KIND_R8),intent(in)  :: X(n) ! input vector argument
   real   (SHR_KIND_R8),intent(out) :: Y(n) ! output vector argument

!-------------------------------------------------------------------------------
! PURPOSE: log for vector arguments, optimized on different platforms
!-------------------------------------------------------------------------------

#if (defined NO_SHR_VMATH)
   Y = log(X)
#else

#if (defined AIX)
   call vlog(Y, X, n)
#endif

#if (defined IRIX64)
   call shr_vmath_fwrap_vlog(X, Y, n)
#endif

#if (defined OSF1)
   call vlog(X, 1, Y, 1, n)
#endif

#if (!defined AIX && !defined IRIX64 && !defined OSF1)
   Y = log(X)
#endif
#endif

end subroutine shr_vmath_log

!===============================================================================

subroutine shr_vmath_sin(X, Y, n)

   !----- arguments ---
   integer(SHR_KIND_IN),intent(in)  ::   n  ! vector length
   real   (SHR_KIND_R8),intent(in)  :: X(n) ! input vector argument
   real   (SHR_KIND_R8),intent(out) :: Y(n) ! output vector argument

!-------------------------------------------------------------------------------
! PURPOSE: sin for vector arguments, optimized on different platforms
!-------------------------------------------------------------------------------

#if (defined NO_SHR_VMATH)
   Y = sin(X)
#else

#if (defined AIX)
   call vsin(Y, X, n)
#endif

#if (defined IRIX64)
   call shr_vmath_fwrap_vsin(X, Y, n)
#endif

#if (defined OSF1)
   call vsin(X, 1, Y, 1, n)
#endif

#if (!defined AIX && !defined IRIX64 && !defined OSF1)
   Y = sin(X)
#endif
#endif

end subroutine shr_vmath_sin

!===============================================================================

subroutine shr_vmath_cos(X, Y, n)

   !----- arguments ---
   integer(SHR_KIND_IN),intent(in)  ::   n  ! vector length
   real   (SHR_KIND_R8),intent(in)  :: X(n) ! input vector argument
   real   (SHR_KIND_R8),intent(out) :: Y(n) ! output vector argument

!-------------------------------------------------------------------------------
! PURPOSE: cos for vector arguments, optimized on different platforms
!-------------------------------------------------------------------------------

#if (defined NO_SHR_VMATH)
   Y = cos(X)
#else

#if (defined AIX)
   call vcos(Y, X, n)
#endif

#if (defined IRIX64)
   call shr_vmath_fwrap_vcos(X, Y, n)
#endif

#if (defined OSF1)
   call vcos(X, 1, Y, 1, n)
#endif

#if (!defined AIX && !defined IRIX64 && !defined OSF1)
   Y = cos(X)
#endif
#endif

end subroutine shr_vmath_cos

!===============================================================================

end module shr_vmath_mod

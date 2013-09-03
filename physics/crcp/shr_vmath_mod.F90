!===============================================================================
! PURPOSE: 
!   provides a uniform, platform-independent API for vector math functions
!===============================================================================

module shr_vmath_mod

   !----------------------------------------------------------------------------
   ! routines that evaluate various math functions for vector arguments
   ! intended to provide platform independent access to vendor optimized code
   !----------------------------------------------------------------------------

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
   integer,intent(in)  ::   n  ! vector length
   real   ,intent(in)  :: X(n) ! input vector argument
   real   ,intent(out) :: Y(n) ! output vector argument

!-------------------------------------------------------------------------------
! PURPOSE: sqrt for vector arguments, optimized on different platforms
!-------------------------------------------------------------------------------

#if (defined NO_SHR_VMATH)
   Y = sqrt(X)
#else

#if (defined _AIX)
   call vsqrt(Y, X, n)
#endif

#if (defined IRIX64)
   call shr_vmath_fwrap_vsqrt(X, Y, n)
#endif

#if (defined OSF1)
   call vsqrt(X, 1, Y, 1, n)
#endif

#if (!defined _AIX && !defined IRIX64 && !defined OSF1)
   Y = sqrt(X)
#endif
#endif

end subroutine shr_vmath_sqrt

!===============================================================================

subroutine shr_vmath_rsqrt(X, Y, n)

   !----- arguments ---
   integer,intent(in)  ::   n  ! vector length
   real   ,intent(in)  :: X(n) ! input vector argument
   real   ,intent(out) :: Y(n) ! output vector argument

!-------------------------------------------------------------------------------
! PURPOSE: sqrt for vector arguments, optimized on different platforms
!-------------------------------------------------------------------------------

#if (defined NO_SHR_VMATH)
   Y = 1.0/sqrt(X)
#else

#if (defined _AIX)
   call vrsqrt(Y, X, n)
#endif

#if (!defined _AIX)
   Y = 1.0/sqrt(X)
#endif
#endif

end subroutine shr_vmath_rsqrt

!===============================================================================

subroutine shr_vmath_exp(X, Y, n)

   !----- arguments ---
   integer,intent(in)  ::   n  ! vector length
   real   ,intent(in)  :: X(n) ! input vector argument
   real   ,intent(out) :: Y(n) ! output vector argument

!-------------------------------------------------------------------------------
! PURPOSE: exp for vector arguments, optimized on different platforms
!-------------------------------------------------------------------------------

#if (defined NO_SHR_VMATH)
   Y = exp(X)
#else

#if (defined _AIX)
   call vexp(Y, X, n)
#endif

#if (defined IRIX64)
   call shr_vmath_fwrap_vexp(X, Y, n)
#endif

#if (defined OSF1)
   call vexp(X, 1, Y, 1, n)
#endif

#if (!defined _AIX && !defined IRIX64 && !defined OSF1)
   Y = exp(X)
#endif
#endif

end subroutine shr_vmath_exp

!===============================================================================

subroutine shr_vmath_div(X, Y, Z, n)
   !----- arguments ---
   integer,intent(in)  ::   n  ! vector length
   real   ,intent(in)  :: X(n) ! input vector argument
   real   ,intent(in)  :: Y(n) ! input vector argument
   real   ,intent(out) :: Z(n) ! output vector argument

#if (defined NO_SHR_VMATH)
   integer :: i
   do i=1,n
      Z(i) = X(i)/Y(i)
   enddo
#else
#if (defined _AIX)
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
   integer,intent(in)  ::   n  ! vector length
   real   ,intent(in)  :: X(n) ! input vector argument
   real   ,intent(out) :: Y(n) ! output vector argument

!-------------------------------------------------------------------------------
! PURPOSE: log for vector arguments, optimized on different platforms
!-------------------------------------------------------------------------------

#if (defined NO_SHR_VMATH)
   Y = log(X)
#else

#if (defined _AIX)
   call vlog(Y, X, n)
#endif

#if (defined IRIX64)
   call shr_vmath_fwrap_vlog(X, Y, n)
#endif

#if (defined OSF1)
   call vlog(X, 1, Y, 1, n)
#endif

#if (!defined _AIX && !defined IRIX64 && !defined OSF1)
   Y = log(X)
#endif
#endif

end subroutine shr_vmath_log

!===============================================================================

subroutine shr_vmath_sin(X, Y, n)

   !----- arguments ---
   integer,intent(in)  ::   n  ! vector length
   real   ,intent(in)  :: X(n) ! input vector argument
   real   ,intent(out) :: Y(n) ! output vector argument

!-------------------------------------------------------------------------------
! PURPOSE: sin for vector arguments, optimized on different platforms
!-------------------------------------------------------------------------------

#if (defined NO_SHR_VMATH)
   Y = sin(X)
#else

#if (defined _AIX)
   call vsin(Y, X, n)
#endif

#if (defined IRIX64)
   call shr_vmath_fwrap_vsin(X, Y, n)
#endif

#if (defined OSF1)
   call vsin(X, 1, Y, 1, n)
#endif

#if (!defined _AIX && !defined IRIX64 && !defined OSF1)
   Y = sin(X)
#endif
#endif

end subroutine shr_vmath_sin

!===============================================================================

subroutine shr_vmath_cos(X, Y, n)

   !----- arguments ---
   integer,intent(in)  ::   n  ! vector length
   real   ,intent(in)  :: X(n) ! input vector argument
   real   ,intent(out) :: Y(n) ! output vector argument

!-------------------------------------------------------------------------------
! PURPOSE: cos for vector arguments, optimized on different platforms
!-------------------------------------------------------------------------------

#if (defined NO_SHR_VMATH)
   Y = cos(X)
#else

#if (defined _AIX)
   call vcos(Y, X, n)
#endif

#if (defined IRIX64)
   call shr_vmath_fwrap_vcos(X, Y, n)
#endif

#if (defined OSF1)
   call vcos(X, 1, Y, 1, n)
#endif

#if (!defined _AIX && !defined IRIX64 && !defined OSF1)
   Y = cos(X)
#endif
#endif

end subroutine shr_vmath_cos

!===============================================================================

end module shr_vmath_mod

!===============================================================================
! SVN $Id: shr_vmath_mod.F90 66411 2014-12-19 22:40:08Z santos@ucar.edu $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_150116/shr/shr_vmath_mod.F90 $
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
#ifndef NO_SHR_VMATH
#if (defined CPRINTEL)
   call vdsqrt(n, X, Y)
   return
#endif

#if (defined AIX)
   call vsqrt(Y, X, n)
   return
#endif

#endif
   Y = sqrt(X)
   return

end subroutine shr_vmath_sqrt

!===============================================================================

subroutine shr_vmath_rsqrt(X, Y, n)

   !----- arguments ---
   integer(SHR_KIND_IN),intent(in)  ::   n  ! vector length
   real   (SHR_KIND_R8),intent(in)  :: X(n) ! input vector argument
   real   (SHR_KIND_R8),intent(out) :: Y(n) ! output vector argument

!-------------------------------------------------------------------------------
! PURPOSE: reciprical sqrt for vector arguments, optimized on different platforms
!-------------------------------------------------------------------------------

#ifndef NO_SHR_VMATH
#if (defined AIX)
   call vrsqrt(Y, X, n)
   return
#endif
!#ifdef CPRINTEL
! Does not pass unit tests
!   real   (SHR_KIND_R8)  :: RX(n) !
!   call vdsqrt(n, X, RX)
!   call vddiv(n, 1.0_SHR_KIND_R8,RX, Y)
!   return
!#endif
#endif
   Y = 1.0_SHR_KIND_R8/sqrt(X)


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

#ifndef NO_SHR_VMATH
#if (defined CPRINTEL)
   call vdexp(n, X, Y)
   return
#endif
#if (defined AIX)
   call vexp(Y, X, n)
   return
#endif
#endif

   Y = exp(X)
   return

end subroutine shr_vmath_exp

!===============================================================================

subroutine shr_vmath_div(X, Y, Z, n)
   !----- arguments ---
   integer(SHR_KIND_IN),intent(in)  ::   n  ! vector length
   real   (SHR_KIND_R8),intent(in)  :: X(n) ! input vector argument
   real   (SHR_KIND_R8),intent(in)  :: Y(n) ! input vector argument
   real   (SHR_KIND_R8),intent(out) :: Z(n) ! output vector argument
   integer :: i
#ifndef NO_SHR_VMATH
#if (defined CPRINTEL)
   call vddiv(n, X, Y, Z)
   return
#endif

#if (defined AIX)
   call vdiv(Z,X,Y,n)
   return
#endif
#endif

   do i=1,n
      Z(i) = X(i)/Y(i)
   enddo
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
#ifndef NO_SHR_VMATH
#if (defined AIX)
   call vlog(Y, X, n)
   return
#endif
#if (defined CPRINTEL)
   call vdln(n, X, Y)
   return
#endif
#endif
   Y = log(X)
   return


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

#ifndef NO_SHR_VMATH
#if (defined AIX)
   call vsin(Y, X, n)
   return
#endif

#if (defined CPRINTEL)
   call vdsin(n, X, Y)
   return
#endif
#endif
   Y = sin(X)
   return

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

#ifndef NO_SHR_VMATH
#if (defined AIX)
   call vcos(Y, X, n)
   return
#endif
#if (defined CPRINTEL)
   call vdcos(n, X, Y)
   return
#endif
#endif
   Y = cos(X)
   return

end subroutine shr_vmath_cos

!===============================================================================

end module shr_vmath_mod

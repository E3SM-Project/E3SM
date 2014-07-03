!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module POP_KindsMod

!BOP
! !MODULE: POP_KindsMod
!
! !DESCRIPTION:
!  This module defines default numerical data types for all common data
!  types like integer, character, logical, real4 and real8.
!
! !USERDOC:
!  Users should not need to adjust anything in this module.  If various
!  character strings like long paths to files exceed the default
!  character length, the default value may be increased.
!
! !REFDOC:
!  This module is supplied to provide consistent data representation
!  across machine architectures.  It is meant to replace the old
!  Fortran double precision and real *X declarations that were
!  implementation-specific.
!  Users should not need to adjust anything in this module.  If various
!  character strings like long paths to files exceed the default
!  character length, the default value may be increased.
!
! !REVISION HISTORY:
!  SVN:$Id: POP_KindsMod.F90 22881 2010-05-11 04:23:39Z njn01 $

! !USES:
!  uses no other modules

   implicit none
   private
   save

! !DEFINED PARAMETERS:

   integer, parameter, public ::                &
      POP_CharLength  = 256                    ,&
      POP_Logical     = kind(.true.)           ,&
      POP_i4          = selected_int_kind(6)   ,&
      POP_i8          = selected_int_kind(13)  ,&
      POP_r4          = selected_real_kind(6)  ,&
      POP_r8          = selected_real_kind(13) ,&
      POP_r16         = selected_real_kind(26)

   integer, parameter, public ::               &
#ifdef TAVG_R8
      POP_rtavg          = POP_r8 ! nonstandard r8 for debugging purposes only
#else
      POP_rtavg          = POP_r4 ! standard, single-precision
#endif

!EOP
!BOC
!EOC
!***********************************************************************

 end module POP_KindsMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

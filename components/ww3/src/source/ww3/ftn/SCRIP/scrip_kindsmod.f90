!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module SCRIP_KindsMod

!BOP
! !MODULE: SCRIP_KindsMod
!
! !DESCRIPTION:
! This module defines default numerical data types for all common data
! types like integer, character, logical, real4 and real8.
!
! !USERDOC:
! Users should not need to adjust anything in this module. If various
! character strings like long paths to files exceed the default
! character length, the default value may be increased.
!
! !REFDOC:
! This module is supplied to provide consistent data representation
! across machine architectures. It is meant to replace the old
! Fortran double precision and real *X declarations that were
! implementation-specific.
! Users should not need to adjust anything in this module. If various
! character strings like long paths to files exceed the default
! character length, the default value may be increased.
!
! !REVISION HISTORY:
! SVN:$Id: SCRIP_KindsMod.F90 82 2008-02-14 19:36:07Z pwjones $

! !USES:
! uses no other modules

   implicit none
   private
   save

! !DEFINED PARAMETERS:

   integer, parameter, public :: &
      SCRIP_CharLength = 100 ,&
      SCRIP_Logical = kind(.true.) ,&
      SCRIP_i4 = selected_int_kind(6) ,&
      SCRIP_i8 = selected_int_kind(13) ,&
      SCRIP_r4 = selected_real_kind(6) ,&
      SCRIP_r8 = selected_real_kind(13) ,&
      SCRIP_r16 = selected_real_kind(26)

!EOP
!BOC
!EOC
!***********************************************************************

 end module SCRIP_KindsMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module kinds_mod

!BOP
! !MODULE: kinds_mod
!
! !DESCRIPTION:
!  This module defines default numerical data types for all common data
!  types like integer, character, logical, real4 and real8.
!
! !REVISION HISTORY:
!  SVN:$Id: kinds_mod.F90 12674 2008-10-31 22:21:32Z njn01 $

! !USES:
!  uses no other modules

   implicit none
   private
   save

! !DEFINED PARAMETERS:

   integer, parameter, public ::               &
      char_len       = 256                    ,&
      char_len_long  = 512                    ,&
      log_kind       = kind(.true.)           ,&
      int_kind       = kind(1)                ,&
      i4             = selected_int_kind(6)   ,&
      i8             = selected_int_kind(13)  ,&
      r4             = selected_real_kind(6)  ,&
      r8             = selected_real_kind(13)

   integer, parameter, public ::               &
#ifdef TAVG_R8
      rtavg          = r8 ! nonstandard r8 for debugging purposes only
#else
      rtavg          = r4 ! standard, single-precision
#endif


!EOP
!BOC
!EOC
!***********************************************************************

 end module kinds_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

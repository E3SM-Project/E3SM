!>
!! @file pio_kinds.F90
!! @brief basic data types
!!
!<
 module pio_kinds

!BOP
! !MODULE: pio_kinds
!
! !DESCRIPTION:
!  This module defines default numerical data types for all common data
!  types like integer, character, logical, real4 and real8.
!
! !REVISION HISTORY:
!  CVS:$Id: pio_kinds.F90,v 1.1.1.1 2006/07/31 16:15:30 dennis Exp $
!  CVS:$Name:  $

! !USES:
!  uses mpi if available
#ifndef NO_MPIMOD
   use mpi, only : MPI_OFFSET_KIND ! _EXTERNAL
#endif

   implicit none
   private
#ifdef NO_MPIMOD
#include <mpif.h>
#endif
! !DEFINED PARAMETERS:

   integer, parameter, public ::          &
      char_len  = 360                    ,&
      log_kind  = kind(.true.)           ,&
      int_kind  = kind(1)                ,&
      i4        = selected_int_kind(6)   ,&
      i8        = selected_int_kind(13)  ,&
      r4        = selected_real_kind(6)  ,&
      r8        = selected_real_kind(13)
!
!  MPI defines MPI_OFFSET_KIND as the byte size of the
!  type, which is not nessasarily the type kind
!

   integer, parameter, public :: PIO_OFFSET_KIND=MPI_OFFSET_KIND

!EOP
!BOC
!EOC
!***********************************************************************

 end module pio_kinds

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

!>
!! @file pio_kinds.F90
!! @brief basic data types
!!
!! $Revision$
!! $LastChangedDate$
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
   include 'mpif.h'   ! _EXTERNAL
#endif
! !DEFINED PARAMETERS:

   integer, parameter, public ::          &
      log_kind  = kind(.true.)           ,&
      int_kind  = kind(1)                ,&
      i4        = selected_int_kind(6)   ,&
      i8        = selected_int_kind(13)  ,&
      r4        = selected_real_kind(6)  ,&
      r8        = selected_real_kind(13)


   integer, parameter, public ::          &
      PIO_OFFSET = MPI_OFFSET_KIND

!EOP
!BOC
!EOC
!***********************************************************************

 end module pio_kinds

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

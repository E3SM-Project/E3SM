!-------------------------------------------------------------------------------
! NASA/GSFC Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: Params
!
!> @brief
!! <BriefDescription>
!!
!! @author
!! Tom Clune, NASA/GSFC 
!!
!! @date
!! 07 Nov 2013
!! 
!! @note <A note here.>
!! <Or starting here...>
!
! REVISION HISTORY:
!
! 07 Nov 2013 - Added the prologue for the compliance with Doxygen. 
!
!-------------------------------------------------------------------------------

module Params_mod
  use ISO_FORTRAN_ENV
  implicit none

  integer, parameter, public :: MAX_LENGTH_NAME = 128

  integer, parameter :: R32 = selected_real_kind(p=6)
  integer, parameter :: R64 = selected_real_kind(p=14)
  integer, parameter :: C32 = selected_real_kind(p=6)
  integer, parameter :: C64 = selected_real_kind(p=14)

  integer, parameter :: I32 = INT32
  integer, parameter :: I64 = INT64

!  integer, parameter :: I32 = selected_int_kind()
!  integer, parameter :: I64 = selected_int_kind()

  integer, parameter :: NEQP=0, EQP=1, GTP=2, GEP=3, LTP=4, LEP=5, &
  &  RELEQP=6

contains

end module Params_mod


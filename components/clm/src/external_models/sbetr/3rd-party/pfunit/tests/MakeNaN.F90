!-------------------------------------------------------------------------------
! NASA/GSFC Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: MakeNaN
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

module MakeNaN_mod
  use Params_mod, only: r32, r64
  implicit none
  private

  public :: makeNaN_32
  public :: makeNaN_64

contains

     function makeNaN_32() result(NaN_32)
        real(r32) :: NaN_32
        integer, parameter :: i32 = selected_int_kind(8)
        integer(i32), parameter :: nan_bits_32 = int(Z'7FA00000',i32)

        NaN_32 = transfer(nan_bits_32, NaN_32)

      end function makeNaN_32

      function makeNaN_64() result(NaN_64)
        real(r64) :: NaN_64
        integer, parameter :: i64 = selected_int_kind(18)
        integer(i64), parameter :: nan_bits_64 = int(Z'7FF4000000000000',i64)

        NaN_64 = transfer(nan_bits_64, NaN_64)

      end function makeNaN_64

end module MakeNaN_mod

!-------------------------------------------------------------------------------
! NASA/GSFC Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: MakeInfinity
!
!> @brief
!! <BriefDescription>
!!
!! @author
!! Tom Clune, NASA/GSFC SIVO
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

module MakeInfinity_mod
   implicit none
   private

   public :: makeInf_32
   public :: makeInf_64

contains

   function makeInf_32() result(Inf_32)
      use Params_mod, only: r32
      real(r32) :: Inf_32
      integer, parameter :: i32 = selected_int_kind(8)
      integer(i32), parameter :: inf_bits_32 = int(Z'7F800000',i32)
      
      Inf_32 = transfer(inf_bits_32, Inf_32)
      
   end function makeInf_32
   
   function makeInf_64() result(Inf_64)
      use Params_mod, only: r64
      real(r64) :: Inf_64
      integer, parameter :: i64 = selected_int_kind(18)
      integer(i64), parameter :: inf_bits_64 = int(Z'7FF0000000000000',i64)
      
      Inf_64 = transfer(inf_bits_64, Inf_64)
      
   end function makeInf_64

end module MakeInfinity_mod


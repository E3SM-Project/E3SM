
module units

   use abortutils, only: endrun
   use shr_file_mod, only: shr_file_getUnit, shr_file_freeUnit

implicit none

PRIVATE

   public :: getunit, freeunit

CONTAINS

   integer function getunit (iu)
!
! Arguments
!
   integer, intent(in), optional :: iu   ! desired unit number
!
! Local workspace
!

     getunit = shr_file_getUnit( iu )

   end function getunit

!#######################################################################

   subroutine freeunit (iu)
!
! Arguments
!
   integer, intent(in) :: iu       ! unit number to be freed

   call shr_file_freeUnit( iu )

   return
   end subroutine freeunit

end module units

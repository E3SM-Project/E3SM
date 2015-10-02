
module srchutil

!----------------------------------------------------------------------- 
! 
! Purpose: Module containing Fortran equivalents to Cray library functions
!          NOTE: some aspects of this code may not meet the CCM coding standard
! 
! Author: Lifted from Cray manuals
! 
!-----------------------------------------------------------------------
#if (! defined UNICOSMP )

CONTAINS

   integer function ismax (n, sx, incx)
!----------------------------------------------------------------------- 
! 
! Purpose: determine index of array maximum
! 
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8

   implicit none
!
! Arguments
!   
   real(r8), intent(in) :: sx(*)        ! array to be searched

   integer, intent(in) :: n             ! number of values to be searched
   integer, intent(in) :: incx          ! increment through array        
!
! Local workspace
!
   integer imax          ! tmp index of max value
   integer i             ! loop index

   real(r8) max          ! maximum value

   max = sx(1)
   imax = 1
   do 100 i=1+incx,n,incx
      if (max .lt. sx(i)) then
         max = sx(i)
         imax = i
      endif
100   continue
      ismax = imax
      return
   end function ismax
 
!===============================================================================

   integer function ismin(n,sx,incx)

   use shr_kind_mod, only: r8 => shr_kind_r8
   
   implicit none
!
! Arguments
!   
   real(r8), intent(in) :: sx(*)  ! array to be searched           

   integer n                      ! number of values to be searched
   integer incx                   ! increment through array        
!
! Local workspace
!
   integer imin                   ! tmp index of min value         
   integer i                      ! loop index                     

   real(r8) min                   ! minimum value                  
   
   min = sx(1)
   imin = 1
   do 100 i=1+incx,n,incx
      if (min .gt. sx(i)) then
         min = sx(i)
         imin = i
      endif
100   continue
      ismin = imin
      return
   end function ismin
 
!===============================================================================

   subroutine whenieq (n, array, inc, target, index, nval)
!----------------------------------------------------------------------- 
! 
! Purpose: Determine indices of "array" which equal "target"
! 
!-----------------------------------------------------------------------

      implicit none
!
! Arguments
!
      integer, intent(in) :: array(*)    ! array to be searched
      integer, intent(in) :: target      ! value to compare against
      integer, intent(in) :: inc         ! increment to move through array

      integer, intent(out) :: nval       ! number of values meeting criteria
      integer, intent(out) :: index(*)   ! output array of indices
!
! Local workspace
!
      integer :: i
      integer :: n
      integer :: ina

      ina=1
      nval=0
      if (inc .lt. 0) ina=(-inc)*(n-1)+1
      do i=1,n
         if(array(ina) .eq. target) then
           nval=nval+1
           index(nval)=i
         end if
         ina=ina+inc
      enddo
      return
   end subroutine whenieq

#endif

end module srchutil


subroutine chktime( time, nrec )

   !----------------------------------------------------------------------- 
   ! Purpose: 
   ! Make sure the time coordinate looks like calander day, and is increasing.
   ! Calendar day can either start with 1 Jan 0Z = day 1.0  or
   !  1 Jan 0Z = day 0.0
   !
   ! Author: B. Eaton
   !----------------------------------------------------------------------- 

   use shr_kind_mod, only: r8 => shr_kind_r8
   use abortutils,   only: endrun
   use cam_logfile,  only: iulog

   implicit none

   integer, intent(in) :: nrec                 ! size of time array
   real(r8), intent(in) :: time(nrec)           ! time coordinate expressed as calendar day.
   
   ! Local varibles:
   integer :: i
   !-----------------------------------------------------------------------

   if ( time(1) .lt. 0._r8 .or. time(1) .ge. 366._r8 ) then
      write(iulog,*)'chktime: illegal time coordinate ',time(1)
      call endrun
   end if
   do i = 2, nrec
      if ( time(i) .lt. 0._r8 .or. time(i) .ge. 366._r8 ) then
         write(iulog,*)'chktime: illegal time coordinate ', time(i)
         call endrun
      end if
      if ( time(i) .le. time(i-1) ) then
         write(iulog,*)'chktime: time not increasing ', time(i-1), time(i)
         call endrun
      end if
   end do

   return

end subroutine chktime

!#######################################################################

subroutine findplb( x, nx, xval, index )

   !----------------------------------------------------------------------- 
   ! Purpose: 
   ! "find periodic lower bound"
   ! Search the input array for the lower bound of the interval that
   ! contains the input value.  The returned index satifies:
   ! x(index) .le. xval .lt. x(index+1)
   ! Assume the array represents values in one cycle of a periodic coordinate.
   ! So, if xval .lt. x(1), or xval .ge. x(nx), then the index returned is nx.
   !
   ! Author: B. Eaton
   !----------------------------------------------------------------------- 

   use shr_kind_mod, only: r8 => shr_kind_r8

   implicit none

   integer, intent(in) ::   nx         ! size of x
   real(r8), intent(in) ::  x(nx)      ! strictly increasing array
   real(r8), intent(in) ::  xval       ! value to be searched for in x
   
   integer, intent(out) ::  index

   ! Local variables:
   integer i
   !-----------------------------------------------------------------------

   if ( xval .lt. x(1) .or. xval .ge. x(nx) ) then
      index = nx
      return
   end if

   do i = 2, nx
      if ( xval .lt. x(i) ) then
         index = i-1
         return
      end if
   end do

end subroutine findplb

!#######################################################################

subroutine linintp( npts, t1, t2, tint, f1, f2, fint )


   !-----------------------------------------------------------------------
   ! Purpose:
   ! Linearly interpolate between f1(t1) and f2(t2) to fint(tint),
   ! where f1, f2, and f3 are chunked data structures
   !
   ! Author: B. Eaton
   ! Chunked by P. Worley
   ! Un-Chunked by P. Rasch (it wasnt right)
   !-----------------------------------------------------------------------

! Linearly interpolate between f1(t1) and f2(t2) to fint(tint).


   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid
   use phys_grid, only: get_ncols_p

      implicit none

! Input arguments:
      integer,   intent(in)::   npts
      real(r8),  intent(in)::   t1               ! time level of f1
      real(r8),  intent(in)::   t2               ! time level of f2
      real(r8),  intent(in)::   tint             ! interpolant time
      real(r8),  intent(in)::   f1(npts)         ! field at time t1
      real(r8),  intent(in)::   f2(npts)         ! field at time t2

! Output arguments:
      real(r8), intent(out)::  fint(npts)       ! field at time tint

! Local variables:
      integer i
      real(r8) factor
!------------------------------------------------------------------------------

!      call t_startf ('linintp')

      factor = ( tint - t1 )/( t2 - t1)

      do i = 1, npts
         fint(i) = f1(i) + ( f2(i) - f1(i) )*factor
      end do

!      call t_stopf ('linintp')

      return
end subroutine linintp

!#######################################################################

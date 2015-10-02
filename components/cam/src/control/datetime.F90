   subroutine datetime(cdate, ctime) 
!-----------------------------------------------------------------------
!
! Purpose:
!
!  A generic Date and Time routine
!
! Author: CCM Core group
!
!-----------------------------------------------------------------------
!
! $Id$
!
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!
!-----------------------------Arguments---------------------------------
   character , intent(out) :: cdate*8 
   character , intent(out) :: ctime*8 
!-----------------------------------------------------------------------
!
!---------------------------Local Variables------------------------------
   integer, dimension(8) :: values 
   character :: date*8, time*10, zone*5 
!-----------------------------------------------------------------------
 
   call date_and_time (date, time, zone, values) 
   cdate(1:2) = date(5:6) 
   cdate(3:3) = '/' 
   cdate(4:5) = date(7:8) 
   cdate(6:6) = '/' 
   cdate(7:8) = date(3:4) 
   ctime(1:2) = time(1:2) 
   ctime(3:3) = ':' 
   ctime(4:5) = time(3:4) 
   ctime(6:6) = ':' 
   ctime(7:8) = time(5:6) 
 
   return  
   end subroutine datetime 
 

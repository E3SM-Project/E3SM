!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module io_tools
 
!BOP
!
! !MODULE:  io_tools
!
! !DESCRIPTION:
!   This module contains routines intended to facilitate io
!   Presently, only routines used to document output are included
!
! !REVISION HISTORY:
!   SVN:$Id: io_tools.F90 923 2006-05-10 22:25:10Z njn01 $
!
! !USES 
   use kinds_mod
   use grid
   use io    
   use shr_sys_mod
      
   implicit none
   save
!EOP
!BOC
 
!-----------------------------------------------------------------------
!     interfaces
!-----------------------------------------------------------------------
 
   interface document
       module procedure document_char, &
                        document_int,  &
                        document_log,  &
                        document_dbl,  &
                        document_real 
   end interface
 
!EOC
!***********************************************************************

 contains
 
!***********************************************************************
!BOP
! !IROUTINE: document
! !INTERFACE:

 subroutine document_char (sub_name,message,char_val)
 
! !DESCRIPTION:
!  This routine writes out the calling subroutine name and two
!  associated messages
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic document
!  routine corresponding to a character string

! !INPUT PARAMETERS:
   character (*) :: sub_name, message 
   character (*), optional :: char_val

!EOP
!BOC
 
   character(*),parameter :: fmt1 = "(  5x, '(',a,') ',a)"
   character(*),parameter :: fmt2 = "(  5x, '(',a,') ',a,' = ',a)" 
 
   if (my_task == master_task) then
     if (present(char_val)) then
       write(stdout,fmt2)  sub_name, message, trim(char_val)
     else
       write(stdout,fmt1)  sub_name, trim(message)
     endif
     call shr_sys_flush (stdout)
   endif
 
!EOC
 
 end subroutine document_char
 
!***********************************************************************
!BOP
! !IROUTINE: document
! !INTERFACE:
 subroutine document_int (sub_name,message,ival)
 
! !DESCRIPTION:
!  This routine writes out the calling subroutine name and an
!  associated message
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic document
!  routine corresponding to an integer

! !INPUT PARAMETERS:
   character (*)       :: sub_name, message
   integer (int_kind)  :: ival

!EOP
!BOC

   character(*),parameter :: fmt = "(  5x, '(',a,') ',a,' = ',i10)" 
 
   if (my_task == master_task) then
      write(stdout,fmt)  sub_name, message, ival
      call shr_sys_flush (stdout)
   endif
 
!EOC
 
 end subroutine document_int
 
!***********************************************************************
!BOP
! !IROUTINE: document
! !INTERFACE:
 subroutine document_log (sub_name,message,lval)
 

! !DESCRIPTION:
!  This routine writes out the calling subroutine name and an
!  associated message
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic document
!  routine corresponding to a logical variable

! !INPUT PARAMETERS:
   character (*)      :: sub_name, message
   logical (log_kind) :: lval

!EOP
!BOC

   character(*),parameter :: fmt = "(  5x, '(',a,') ',a,' = ',L3)" 
 
   if (my_task == master_task) then
      write(stdout,fmt)  sub_name, message, lval
      call shr_sys_flush (stdout)
   endif
 
!EOC
 
 end subroutine document_log
 
!***********************************************************************
!BOP
! !IROUTINE: document
! !INTERFACE:
 subroutine document_dbl (sub_name,message,dval)

! !DESCRIPTION:
!  This routine writes out the calling subroutine name and an
!  associated message
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic document
!  routine corresponding to a r8 variable

! !INPUT PARAMETERS:
   character (*) :: sub_name, message
   real (r8)     :: dval
 
!EOP
!BOC

   character(*),parameter :: fmt = "(  5x, '(',a,') ',a,' = ',1pe23.16)" 
 
   if (my_task == master_task) then
      write(stdout,fmt)  sub_name, message, dval
      call shr_sys_flush (stdout)
   endif
 
!EOC
 
 end subroutine document_dbl
 
!***********************************************************************
!BOP
! !IROUTINE: document
! !INTERFACE:
 subroutine document_real (sub_name,message,rval)
 

! !DESCRIPTION:
!  This routine writes out the calling subroutine name and an
!  associated message
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic document
!  routine corresponding to a r8 variable

! !INPUT PARAMETERS:
   character (*) :: sub_name, message
   real(r4)    :: rval
 
!EOP
!BOC

   character(*),parameter :: fmt = "(  5x, '(',a,') ',a,' = ',1pe23.16)" 
 
   if (my_task == master_task) then
      write(stdout,fmt)  sub_name, message, rval
      call shr_sys_flush (stdout)
   endif
 
!EOC
 
 end subroutine document_real
 
 
 end module io_tools

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

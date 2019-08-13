module shr_nl_mod

! Utilities for namelist reading
! Adapted Fall 2012 from CAM's namelist_utils.

implicit none
private

save

public :: &
     shr_nl_find_group_name ! seek through a file to find a specified namelist

contains

! This routine probably discards more error code information than it needs to.

subroutine shr_nl_find_group_name(unit, group, status)

  use shr_string_mod, only: shr_string_toLower

!---------------------------------------------------------------------------------------
! Purpose:
! Search a file that contains namelist input for the specified namelist group name.
! Leave the file positioned so that the current record is the first record of the
! input for the specified group.
!
! Method:
! Read the file line by line.  Each line is searched for an '&' which may only
! be preceded by blanks, immediately followed by the group name which is case
! insensitive.  If found then backspace the file so the current record is the
! one containing the group name and return success.  Otherwise return -1.
!
! Author:  B. Eaton, August 2007
!---------------------------------------------------------------------------------------

   integer,          intent(in)  :: unit     ! fortran unit attached to file
   character(len=*), intent(in)  :: group    ! namelist group name
   integer,          intent(out) :: status   ! 0 for success, -1 if group name not found

   ! Local variables

   integer           :: len_grp
   integer           :: ios    ! io status
   character(len=80) :: inrec  ! first 80 characters of input record
   character(len=80) :: inrec2 ! left adjusted input record
   character(len=len(group)) :: lc_group

   !---------------------------------------------------------------------------

   len_grp = len_trim(group)
   lc_group = shr_string_toLower(group)

   ios = 0
   do while (ios <= 0)

      read(unit, '(a)', iostat=ios, end=100) inrec

      if (ios <= 0) then  ! ios < 0  indicates an end of record condition

         ! look for group name in this record

         ! remove leading blanks
         inrec2 = adjustl(inrec)

         ! check for leading '&'
         if (inrec2(1:1) == '&') then

            ! check for case insensitive group name
            if (trim(lc_group) == shr_string_toLower(inrec2(2:len_grp+1))) then

               ! found group name.  backspace to leave file position at this record
               backspace(unit)
               status = 0
               return

            end if
         end if
      end if

   end do

   100 continue  ! end of file processing
   status = -1

end subroutine shr_nl_find_group_name

end module shr_nl_mod

module elm_nlUtilsMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_nltUtilsMod
!
! !DESCRIPTION:
! Utilities to handle namelists.
!
! !USES:

! !PUBLIC TYPES:
  implicit none
  save

  private    ! By default everything is private

! !PUBLIC MEMBER FUNCTIONS:
  public :: find_nlgroup_name   ! find a specified namelist group in a file
!
! !REVISION HISTORY:
! Created by B. Eaton
! Move to CLM by E. Kluzek
!
! !PRIVATE MEMBER FUNCTIONS: None
!-----------------------------------------------------------------------
! !PRIVATE DATA MEMBERS: None

!EOP
!-----------------------------------------------------------------------
contains

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: find_nlgroup_name
!
! !INTERFACE:
  subroutine find_nlgroup_name(unit, group, status)
!
! !DESCRIPTION:
! Search a file that contains namelist input for the specified namelist group name.
! Leave the file positioned so that the current record is the first record of the
! input for the specified group.
! 
! METHOD: 
! Read the file line by line.  Each line is searched for an '&' which may only
! be preceded by blanks, immediately followed by the group name which is case
! insensitive.  If found then backspace the file so the current record is the
! one containing the group name and return success.  Otherwise return -1.
!
! !USES:
   use shr_kind_mod  , only : CS => shr_kind_cs
   use shr_string_mod, only : shr_string_toLower
!
! !ARGUMENTS:
   integer,          intent(in)  :: unit     ! fortran unit attached to file
   character(len=*), intent(in)  :: group    ! namelist group name
   integer,          intent(out) :: status   ! 0 for success, -1 if group name not found
!
! !REVISION HISTORY:
! Created by B. Eaton, August 2007
! Move to CLM E. Kluzek, August 2012
!
!
! !LOCAL VARIABLES:
!EOP
   integer           :: len_grp  ! length of the groupname
   integer           :: ios      ! io status
   character(len=CS) :: inrec    ! first shr_kind_CS characters of input record
   character(len=CS) :: inrec2   ! left adjusted input record
   character(len=len(group)) :: lc_group ! lower-case group name
   character(len=32) :: subname = 'find_nlgroup_name' ! subroutine name
!-----------------------------------------------------------------------
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

end subroutine find_nlgroup_name

!-----------------------------------------------------------------------

end module elm_nlUtilsMod

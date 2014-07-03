!===============================================================================
! SVN $Id: shr_msg_mod.F90 6752 2007-10-04 21:02:15Z jwolfe $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_140509/shr/shr_msg_mod.F90 $
!===============================================================================

!BOP ===========================================================================
!
! !MODULE: shr_msg_mod.F90 --- Module to handle various file utilily functions.
!
! !DESCRIPTION:
!
! Miscilaneous methods to handle file and directory utilities needed for 
! message passing CCSM.
!
! NOTE:
!
! This module can be replaced by the direct use of shr_file_mod.F90 as all
! real functionality was moved to shr_file_mod.
!
! !REVISION HISTORY:
!   2006-05-08 E. Kluzek, move functionality to shr_file_mod.F90.
!   2000-??-?? B. Kauffman, original version circa 2000
!
! !INTERFACE: ------------------------------------------------------------------

MODULE shr_msg_mod

! ! USES:
   use shr_file_mod   ! The real guts of everything here
   use shr_log_mod, only: s_loglev  => shr_log_Level
   use shr_log_mod, only: s_logunit => shr_log_Unit

   IMPLICIT none
   PRIVATE           ! By default everything is private to this module

! !PUBLIC TYPES:               

   ! no public types

! !PUBLIC MEMBER FUNCTIONS:
   public :: shr_msg_chDir    ! change current working directory
   public :: shr_msg_chStdIn  ! change stdin  (attach to a file)
   public :: shr_msg_chStdOut ! change stdout (attach to a file)
   public :: shr_msg_stdio    ! change dir and stdin and stdout
   public :: shr_msg_dirio    ! change stdin and stdout

! !PUBLIC DATA MEMBERS:

   ! no public data members

!EOP

!===============================================================================
CONTAINS
!===============================================================================

!BOP ===========================================================================
!
! !IROUTINE: shr_msg_stdio -- Change working directory, and redirect stdin/stdout
!
! !DESCRIPTION:
!   1) change the cwd (current working directory) and 
!   2) redirect stdin & stdout (units 5 & 6) to named files,
!   where the desired cwd & files are specified by namelist file.
!
! !INTERFACE: ------------------------------------------------------------------  
SUBROUTINE shr_msg_stdio(model)

   implicit none

! !INPUT/OUTPUT PARAMETERS:
   !--- arguments ---
   character(len=*),intent(in) :: model ! used to construct env varible name
!EOP

!-------------------------------------------------------------------------------
! Notes:
!
!-------------------------------------------------------------------------------

     call shr_file_stdio(model)
 
END SUBROUTINE shr_msg_stdio

!===============================================================================

!BOP ===========================================================================
!
! !IROUTINE: shr_msg_chdir -- Change working directory.
!
! !DESCRIPTION:
!   change the cwd (current working directory), see shr_msg_stdio for notes
!
! !INTERFACE: ------------------------------------------------------------------  

SUBROUTINE shr_msg_chdir(model)

! !USES:
   use shr_sys_mod, only: shr_sys_chdir

   implicit none

! !INPUT/OUTPUT PARAMETERS:
   !--- arguments ---
   character(len=*),intent(in) :: model ! used to construct env varible name
!EOP

   !--- local ---

!-------------------------------------------------------------------------------
! Notes:
!
!-------------------------------------------------------------------------------

   call shr_file_chdir( model )
 
END SUBROUTINE shr_msg_chdir

!===============================================================================

!BOP ===========================================================================
!
! !IROUTINE: shr_msg_dirio --- Change stdin and stdout.
!
! !DESCRIPTION:
!   change the stdin & stdout (units 5 & 6), see shr_msg_stdio for notes
!
! !INTERFACE: ------------------------------------------------------------------  
SUBROUTINE shr_msg_dirio(model)

   implicit none

! !INPUT/OUTPUT PARAMETERS:
   !--- arguments ---
   character(len=*),intent(in) :: model ! used to construct env varible name
!EOP

   !--- local ---

!-------------------------------------------------------------------------------
! Notes:
!
!-------------------------------------------------------------------------------

   call shr_file_dirio(model)
 
END SUBROUTINE shr_msg_dirio

!===============================================================================

!BOP ===========================================================================
!
! !IROUTINE: shr_msg_chStdIn -- Change stdin
!
! !DESCRIPTION:
!   change the stdin (unit 5), see shr_msg_stdio for notes
!
! !INTERFACE: ------------------------------------------------------------------  
SUBROUTINE shr_msg_chStdIn(model)

   implicit none

! !INPUT/OUTPUT PARAMETERS:
   !--- arguments ---
   character(*),intent(in) :: model ! used to construct env varible name
!EOP
   !--- local ---

!-------------------------------------------------------------------------------
! Notes:
!
!-------------------------------------------------------------------------------

    call shr_file_chStdIn(model)
 
END SUBROUTINE shr_msg_chStdIn

!===============================================================================

!BOP ===========================================================================
!
! !IROUTINE: shr_msg_stdout -- Change stdout
!
! !DESCRIPTION:
!   change the stdout (unit 6), see shr_msg_stdio for notes
!
! !INTERFACE: ------------------------------------------------------------------  

SUBROUTINE shr_msg_chStdOut(model)

   implicit none

! !INPUT/OUTPUT PARAMETERS:
   !--- arguments ---
   character(*),intent(in) :: model ! used to construct env varible name
!EOP

   !--- local ---

!-------------------------------------------------------------------------------
! Notes:
!
!-------------------------------------------------------------------------------

    call shr_file_chStdOut(model)

END SUBROUTINE shr_msg_chStdOut

!===============================================================================

END MODULE shr_msg_mod

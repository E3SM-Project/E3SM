!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module SCRIP_IOUnitsMod

!BOP
!
! !MODULE: SCRIP_IOUnitsMod
!
! !DESCRIPTION:
! This module contains an I/O unit manager for tracking, assigning
! and reserving I/O unit numbers.
!
! There are three reserved I/O units set as parameters in this
! module. The default units for standard input (stdin), standard
! output (stdout) and standard error (stderr). These are currently
! set as units 5,6,6, respectively as that is the most commonly
! used among vendors. However, the user may change these if those
! default units are conflicting with other models or if the
! vendor is using different values.
!
! The maximum number of I/O units per node is currently set by
! the parameter SCRIP\_IOMaxUnits.
!
! !REFDOC:
!
! !REVISION HISTORY:
! SVN:$Id: SCRIP_IOUnitsMod.F90 83 2008-02-22 17:26:54Z pwjones $

! !USES:

   use SCRIP_KindsMod

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: SCRIP_IOUnitsGet, &
             SCRIP_IOUnitsRelease, &
             SCRIP_IOUnitsReserve, &
             SCRIP_IOUnitsRedirect, &
             SCRIP_IOUnitsFlush

! !PUBLIC DATA MEMBERS:

   integer (SCRIP_i4), parameter, public :: &
      SCRIP_stdin = 5, &! reserved unit for standard input
      SCRIP_stdout = 6, &! reserved unit for standard output
      SCRIP_stderr = 6 ! reserved unit for standard error

   ! common formats for writing to stdout, stderr

   character (9), parameter, public :: &
      SCRIP_delimFormat = "(72('-'))"

   character (5), parameter, public :: &
      SCRIP_blankFormat = "(' ')"

!EOP
!BOC
!-----------------------------------------------------------------------
!
! private io unit manager variables
!
!-----------------------------------------------------------------------

   integer (SCRIP_i4), parameter :: &
      SCRIP_IOUnitsMinUnits = 11, & ! do not use unit numbers below this
      SCRIP_IOUnitsMaxUnits = 99 ! maximum number of open units

   logical (SCRIP_Logical) :: &
      SCRIP_IOUnitsInitialized = .false.

   logical (SCRIP_Logical), dimension(SCRIP_IOUnitsMaxUnits) :: &
      SCRIP_IOUnitsInUse ! flag=.true. if unit currently open

!EOC
!***********************************************************************

contains

!***********************************************************************
!BOP
! !IROUTINE: SCRIP_IOUnitsGet
! !INTERFACE:

 subroutine SCRIP_IOUnitsGet(iunit)

! !DESCRIPTION:
! This routine returns the next available i/o unit and marks it as
! in use to prevent any later use.
! Note that {\em all} processors must call this routine even if only
! the master task is doing the i/o. This is necessary insure that
! the units remain synchronized for other parallel I/O functions.
!
! !REVISION HISTORY:
! same as module

! !OUTPUT PARAMETERS:

   integer (SCRIP_i4), intent(out) :: &
      iunit ! next free i/o unit

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer (SCRIP_i4) :: n ! dummy loop index

   logical (SCRIP_Logical) :: alreadyInUse

!-----------------------------------------------------------------------
!
! check to see if units initialized and initialize if necessary
!
!-----------------------------------------------------------------------

   if (.not. SCRIP_IOUnitsInitialized) then
      SCRIP_IOUnitsInUse = .false.
      SCRIP_IOUnitsInUse(SCRIP_stdin) = .true.
      SCRIP_IOUnitsInUse(SCRIP_stdout) = .true.
      SCRIP_IOUnitsInUse(SCRIP_stderr) = .true.

      SCRIP_IOUnitsInitialized = .true.
   endif

!-----------------------------------------------------------------------
!
! find next free unit
!
!-----------------------------------------------------------------------

   srch_units: do n=SCRIP_IOUnitsMinUnits, SCRIP_IOUnitsMaxUnits
      if (.not. SCRIP_IOUnitsInUse(n)) then ! I found one, I found one

         !*** make sure not in use by library or calling routines
         INQUIRE (unit=n,OPENED=alreadyInUse)

         if (.not. alreadyInUse) then
            iunit = n ! return the free unit number
            SCRIP_IOUnitsInUse(iunit) = .true. ! mark iunit as being in use
            exit srch_units
         else
            !*** if inquire shows this unit in use, mark it as
            !*** in use to prevent further queries
            SCRIP_IOUnitsInUse(n) = .true.
         endif
      endif
   end do srch_units

   if (iunit > SCRIP_IOUnitsMaxUnits) &
      stop 'SCRIP_IOUnitsGet: No free units'

!-----------------------------------------------------------------------
!EOC

 end subroutine SCRIP_IOUnitsGet

!***********************************************************************
!BOP
! !IROUTINE: SCRIP_IOUnitsRelease
! !INTERFACE:

 subroutine SCRIP_IOUnitsRelease(iunit)

! !DESCRIPTION:
! This routine releases an i/o unit (marks it as available).
! Note that {\em all} processors must call this routine even if only
! the master task is doing the i/o. This is necessary insure that
! the units remain synchronized for other parallel I/O functions.
!
! !REVISION HISTORY:
! same as module

! !INPUT PARAMETER:

   integer (SCRIP_i4), intent(in) :: &
      iunit ! i/o unit to be released

!EOP
!BOC
!-----------------------------------------------------------------------
!
! check for proper unit number
!
!-----------------------------------------------------------------------

   if (iunit < 1 .or. iunit > SCRIP_IOUnitsMaxUnits) then
      stop 'SCRIP_IOUnitsRelease: bad unit'
   endif

!-----------------------------------------------------------------------
!
! mark the unit as not in use
!
!-----------------------------------------------------------------------

   SCRIP_IOUnitsInUse(iunit) = .false. ! that was easy...

!-----------------------------------------------------------------------
!EOC

 end subroutine SCRIP_IOUnitsRelease

!***********************************************************************
!BOP
! !IROUTINE: SCRIP_IOUnitsReserve
! !INTERFACE:

 subroutine SCRIP_IOUnitsReserve(iunit)

! !DESCRIPTION:
! This routine marks an IO unit as in use to reserve its use
! for purposes outside of SCRIP IO. This is necessary for
! cases where you might be importing code developed elsewhere
! that performs its own I/O and open/closes units.
! Note that {\em all} processors must call this routine even if only
! the master task is doing the i/o. This is necessary insure that
! the units remains synchronized for other parallel I/O functions.
!
! !REVISION HISTORY:
! same as module

! !INPUT PARAMETER:

   integer (SCRIP_i4), intent(in) :: &
      iunit ! i/o unit to be reserved

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   logical (SCRIP_Logical) :: alreadyInUse

!-----------------------------------------------------------------------
!
! check for proper unit number
!
!-----------------------------------------------------------------------

   if (iunit < SCRIP_IOUnitsMinUnits .or. &
       iunit > SCRIP_IOUnitsMaxUnits) then
      stop 'SCRIP_IOUnitsReserve: invalid unit'
   endif

!-----------------------------------------------------------------------
!
! check to see if SCRIP already using this unit
!
!-----------------------------------------------------------------------

   if (SCRIP_IOUnitsInUse(iunit)) then
      stop 'SCRIP_IOUnitsReserve: unit already in use by SCRIP'
   endif

!-----------------------------------------------------------------------
!
! check to see if others already using this unit
!
!-----------------------------------------------------------------------

   INQUIRE (unit=iunit, OPENED=alreadyInUse)
   if (alreadyInUse) then
      stop 'SCRIP_IOUnitsReserve: unit already in use by others'
   endif

!-----------------------------------------------------------------------
!
! mark the unit as in use
!
!-----------------------------------------------------------------------

   SCRIP_IOUnitsInUse(iunit) = .true.

!-----------------------------------------------------------------------
!EOC

 end subroutine SCRIP_IOUnitsReserve

!***********************************************************************
!BOP
! !IROUTINE: SCRIP_IOUnitsRedirect
! !INTERFACE:

 subroutine SCRIP_IOUnitsRedirect(iunit, filename)

! !DESCRIPTION:
! This routine enables a user to redirect stdin, stdout, stderr to
! a file instead of to the terminal. It is only permitted for these
! special units. The SCRIP IO file operators should be used for
! normal I/O.
! Note that {\em all} processors must call this routine even if only
! the master task is doing the i/o. This is necessary insure that
! the units remains synchronized for other parallel I/O functions.
!
! !REVISION HISTORY:
! same as module

! !INPUT PARAMETER:

   integer (SCRIP_i4), intent(in) :: &
      iunit ! i/o unit to be redirected to file

   character (*), intent(in) :: &
      filename ! filename, including path, to which
                               ! i/o should be directed

!EOP
!BOC
!-----------------------------------------------------------------------
!
! check for proper unit number and open file
!
!-----------------------------------------------------------------------

   if (iunit == SCRIP_stdin) then ! open input file for stdin
      open(unit=iunit, file=filename, status='old', form='formatted')

   else if (iunit == SCRIP_stdout) then ! open output file for stdout
      open(unit=iunit, file=filename, status='unknown', form='formatted')

   else if (iunit == SCRIP_stderr .and. SCRIP_stderr /= SCRIP_stdout) then
      ! open output file for stderr
      open(unit=iunit, file=filename, status='unknown', form='formatted')

   else
      stop 'SCRIP_IOUnitsRedirect: invalid unit'

   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine SCRIP_IOUnitsRedirect

!***********************************************************************
!BOP
! !IROUTINE: SCRIP_IOUnitsFlush
! !INTERFACE:

 subroutine SCRIP_IOUnitsFlush(iunit)

! !DESCRIPTION:
! This routine enables a user to flush the output from an IO unit
! (typically stdout) to force output when the system is buffering
! such output. Because this system function is system dependent,
! we only support this wrapper and users are welcome to insert the
! code relevant to their local machine.
!
! !REVISION HISTORY:
! same as module

! !INPUT PARAMETER:

   integer (SCRIP_i4), intent(in) :: &
      iunit ! i/o unit to be flushed

!EOP
!BOC
!-----------------------------------------------------------------------
!
! insert your system code here
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!EOC

 end subroutine SCRIP_IOUnitsFlush

!***********************************************************************

 end module SCRIP_IOUnitsMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module POP_IOUnitsMod

!BOP
!
! !MODULE:  POP_IOUnitsMod
!
! !DESCRIPTION:
!  This module contains an I/O unit manager for tracking, assigning
!  and reserving I/O unit numbers.
!
!  There are three reserved I/O units set as parameters in this
!  module.  The default units for standard input (stdin), standard
!  output (stdout) and standard error (stderr).  These are currently
!  set as units 5,6,6, respectively as that is the most commonly
!  used among vendors. However, the user may change these if those
!  default units are conflicting with other models or if the
!  vendor is using different values.
!
!  The maximum number of I/O units per node is currently set by
!  the parameter POP\_IOMaxUnits.
!
! !REFDOC:
!
! !REVISION HISTORY:
!  SVN:$Id: POP_IOUnitsMod.F90 29587 2011-08-03 22:59:39Z mvertens $
!  2006-08-21: Phil Jones
!     added wrapper for system flush routine
!  2006-08-15: Phil Jones
!     fixed problem in case construct for duplicate unit numbers
!     stripped DOS line feed-CR stuff
!  2006-07-05: Phil Jones
!     added new IO unit manager to follow new naming conventions
!         and check to see if unit assigned by another component
!         if POP is being called as subroutine in a coupled context

! !USES:

   use POP_KindsMod
#ifdef CCSMCOUPLED
   use shr_sys_mod
   use shr_file_mod
#endif

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: POP_IOUnitsGet,                &
             POP_IOUnitsRelease,            &
             POP_IOUnitsReserve,            &
             POP_IOUnitsRedirect,           &
             POP_IOUnitsFlush

! !PUBLIC DATA MEMBERS:

#ifdef CCSMCOUPLED
   integer (POP_i4),            public :: &
#else
   integer (POP_i4), parameter, public :: &
#endif
      POP_stdin  =  5,  &! reserved unit for standard input
      POP_stdout =  6,  &! reserved unit for standard output
      POP_stderr =  6    ! reserved unit for standard error

   ! common formats for writing to stdout, stderr

   character (9), parameter, public ::   &
      POP_delimFormat    = "(72('-'))",  &
      POP_delimFormatNew = "(72('='))"

   character (5), parameter, public :: &
      POP_blankFormat = "(' ')" 

#ifdef CCSMCOUPLED
      ! instance control
      integer (POP_i4)  , public :: inst_index
      character(len=16) , public :: inst_name
      character(len=16) , public :: inst_suffix
#endif

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  private io unit manager variables
!
!-----------------------------------------------------------------------

   integer (POP_i4), parameter :: &
      POP_IOUnitsMinUnits = 11,   & ! do not use unit numbers below this
      POP_IOUnitsMaxUnits = 99      ! maximum number of open units

   logical (POP_Logical) :: &
      POP_IOUnitsInitialized = .false.

   logical (POP_Logical), dimension(POP_IOUnitsMaxUnits) :: &
      POP_IOUnitsInUse       ! flag=.true. if unit currently open

!EOC
!***********************************************************************

contains

!***********************************************************************
!BOP
! !IROUTINE: POP_IOUnitsGet
! !INTERFACE:

 subroutine POP_IOUnitsGet(iunit)

! !DESCRIPTION:
!  This routine returns the next available i/o unit and marks it as
!  in use to prevent any later use.
!  Note that {\em all} processors must call this routine even if only
!  the master task is doing the i/o.  This is necessary insure that
!  the units remain synchronized for other parallel I/O functions.
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      iunit                     ! next free i/o unit

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: n  ! dummy loop index

   logical (POP_Logical) :: alreadyInUse

!-----------------------------------------------------------------------
!
!  check to see if units initialized and initialize if necessary
!
!-----------------------------------------------------------------------

   if (.not. POP_IOUnitsInitialized) then
      POP_IOUnitsInUse = .false.
      POP_IOUnitsInUse(POP_stdin) = .true.
      POP_IOUnitsInUse(POP_stdout) = .true.
      POP_IOUnitsInUse(POP_stderr) = .true.

      POP_IOUnitsInitialized = .true.
   endif

!-----------------------------------------------------------------------
!
!  find next free unit
!
!-----------------------------------------------------------------------

#ifdef CCSMCOUPLED
   iunit = shr_file_getUnit()
#else
   srch_units: do n=POP_IOUnitsMinUnits, POP_IOUnitsMaxUnits
      if (.not. POP_IOUnitsInUse(n)) then   ! I found one, I found one

         !*** make sure not in use by library or calling routines
         INQUIRE (unit=n,OPENED=alreadyInUse)

         if (.not. alreadyInUse) then
            iunit = n        ! return the free unit number
            POP_IOUnitsInUse(iunit) = .true.  ! mark iunit as being in use
            exit srch_units
         else
            !*** if inquire shows this unit in use, mark it as
            !***    in use to prevent further queries
            POP_IOUnitsInUse(n) = .true.
         endif
      endif
   end do srch_units

   if (iunit > POP_IOUnitsMaxUnits) stop 'POP_IOUnitsGet: No free units'
#endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_IOUnitsGet

!***********************************************************************
!BOP
! !IROUTINE: POP_IOUnitsRelease
! !INTERFACE:

 subroutine POP_IOUnitsRelease(iunit)

! !DESCRIPTION:
!  This routine releases an i/o unit (marks it as available).
!  Note that {\em all} processors must call this routine even if only
!  the master task is doing the i/o.  This is necessary insure that
!  the units remain synchronized for other parallel I/O functions.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETER:

   integer (POP_i4), intent(in) :: &
      iunit                    ! i/o unit to be released

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  check for proper unit number
!
!-----------------------------------------------------------------------

   if (iunit < 1 .or. iunit > POP_IOUnitsMaxUnits) then
      stop 'POP_IOUnitsRelease: bad unit'
   endif

!-----------------------------------------------------------------------
!
!  mark the unit as not in use
!
!-----------------------------------------------------------------------

#ifdef CCSMCOUPLED
   call shr_file_freeUnit(iunit)
#else
   POP_IOUnitsInUse(iunit) = .false.  !  that was easy...
#endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_IOUnitsRelease

!***********************************************************************
!BOP
! !IROUTINE: POP_IOUnitsReserve
! !INTERFACE:

 subroutine POP_IOUnitsReserve(iunit)

! !DESCRIPTION:
!  This routine marks an IO unit as in use to reserve its use
!  for purposes outside of POP IO.  This is necessary for
!  cases where you might be importing code developed elsewhere
!  that performs its own I/O and open/closes units.
!  Note that {\em all} processors must call this routine even if only
!  the master task is doing the i/o.  This is necessary insure that
!  the units remains synchronized for other parallel I/O functions.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETER:

   integer (POP_i4), intent(in) :: &
      iunit                    ! i/o unit to be reserved

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   logical (POP_Logical) :: alreadyInUse

!-----------------------------------------------------------------------
!
!  check for proper unit number
!
!-----------------------------------------------------------------------

   if (iunit < POP_IOUnitsMinUnits .or. iunit > POP_IOUnitsMaxUnits) then
      stop 'POP_IOUnitsReserve: invalid unit'
   endif

!-----------------------------------------------------------------------
!
!  check to see if POP already using this unit
!
!-----------------------------------------------------------------------

   if (POP_IOUnitsInUse(iunit)) then
      stop 'POP_IOUnitsReserve: unit already in use by POP'
   endif

!-----------------------------------------------------------------------
!
!  check to see if others already using this unit
!
!-----------------------------------------------------------------------

   INQUIRE (unit=iunit, OPENED=alreadyInUse)
   if (alreadyInUse) then
      stop 'POP_IOUnitsReserve: unit already in use by others'
   endif

!-----------------------------------------------------------------------
!
!  mark the unit as in use
!
!-----------------------------------------------------------------------

   POP_IOUnitsInUse(iunit) = .true.  !  that was easy...

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_IOUnitsReserve

!***********************************************************************
!BOP
! !IROUTINE: POP_IOUnitsRedirect
! !INTERFACE:

 subroutine POP_IOUnitsRedirect(iunit, filename)

! !DESCRIPTION:
!  This routine enables a user to redirect stdin, stdout, stderr to
!  a file instead of to the terminal.  It is only permitted for these
!  special units.  The POP IO file operators should be used for
!  normal I/O.
!  Note that {\em all} processors must call this routine even if only
!  the master task is doing the i/o.  This is necessary insure that
!  the units remains synchronized for other parallel I/O functions.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETER:

   integer (POP_i4), intent(in) :: &
      iunit                    ! i/o unit to be redirected to file

   character (*), intent(in) :: &
      filename                 ! filename, including path, to which
                               !   i/o should be directed

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  check for proper unit number and open file
!
!-----------------------------------------------------------------------

   if (iunit == POP_stdin) then ! open input file for stdin
      open(unit=iunit, file=filename, status='old', form='formatted')

   else if (iunit == POP_stdout) then ! open output file for stdout
      open(unit=iunit, file=filename, status='unknown', form='formatted')

   else if (iunit == POP_stderr .and. POP_stderr /= POP_stdout)  then
      ! open output file for stderr
      open(unit=iunit, file=filename, status='unknown', form='formatted')

   else
      stop 'POP_IOUnitsRedirect: invalid unit'

   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_IOUnitsRedirect

!***********************************************************************
!BOP
! !IROUTINE: POP_IOUnitsFlush
! !INTERFACE:

 subroutine POP_IOUnitsFlush(iunit)

! !DESCRIPTION:
!  This routine enables a user to flush the output from an IO unit
!  (typically stdout) to force output when the system is buffering
!  such output.  Because this system function is system dependent,
!  we only support this wrapper and users are welcome to insert the
!  code relevant to their local machine.  In the case where the CCSM
!  libraries are available, the shared routine for sys flush can be
!  used (and is provided here under a preprocessor option).
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETER:

   integer (POP_i4), intent(in) :: &
      iunit                    ! i/o unit to be flushed

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  insert your system code here
!
!-----------------------------------------------------------------------

#ifdef CCSMCOUPLED

   call shr_sys_flush(iunit)

#endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_IOUnitsFlush

!***********************************************************************

 end module POP_IOUnitsMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

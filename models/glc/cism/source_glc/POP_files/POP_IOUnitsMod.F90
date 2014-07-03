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
! !USERDOC:
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
!  SVN:$Id: POP_IOUnitsMod.F90 808 2006-04-28 17:06:38Z njn01 $

! !USES:

   use POP_KindsMod

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: POP_IOUnitsGet,                &
             POP_IOUnitsRelease,            &
             POP_IOUnitsReserve

! !PUBLIC DATA MEMBERS:

   integer (POP_i4), parameter, public :: &
      POP_stdin  =  5,  &! reserved unit for standard input
      POP_stdout =  6,  &! reserved unit for standard output
      POP_stderr =  6    ! reserved unit for standard error

   ! common formats for writing to stdout, stderr

   character (9), parameter, public :: &
      POP_delimFormat = "(72('-'))"

   character (5), parameter, public :: &
      POP_blankFormat = "(' ')" 

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  io unit manager variables
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
!  the units remains synchronized for other parallel I/O functions.
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

   srch_units: do n=POP_IOUnitsMinUnits, POP_IOUnitsMaxUnits
      if (.not. POP_IOUnitsInUse(n)) then   ! I found one, I found one

         !*** make sure not in use by library routines
         INQUIRE (unit=n,OPENED=alreadyInUse)

         if (.not. alreadyInUse) then         
            iunit = n
            POP_IOUnitsInUse(iunit) = .true.  ! mark iunit as being in use
            exit srch_units
         endif
      endif
   end do srch_units

   if (iunit > POP_IOUnitsMaxUnits) stop 'POP_IOUnitsGet: No free units'

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

   POP_IOUnitsInUse(iunit) = .false.  !  that was easy...

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_IOUnitsRelease

!***********************************************************************
!BOP
! !IROUTINE: POP_IOUnitsReserve
! !INTERFACE:

 subroutine POP_IOUnitsReserve(iunit)

! !DESCRIPTION:
!  This routine releases an i/o unit (marks it as available).
!  Note that {\em all} processors must call this routine even if only
!  the master task is doing the i/o.  This is necessary insure that
!  the units remains synchronized for other parallel I/O functions.
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

 end module POP_IOUnitsMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
! CVS $Id: m_stdio.F90 8708 2008-01-29 03:51:08Z robj $
! CVS $Name: MCT_2_5_0 $  
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: m_stdio - a F90 module defines std. I/O parameters
!
! !DESCRIPTION:
!	Define system dependent I/O parameters.
!
! !INTERFACE:

	module m_stdio
	implicit none
	private

	public	:: stdin	! a unit linked to UNIX stdin
	public	:: stdout	! a unit linked to UNIX stdout
	public	:: stderr	! a unit linked to UNIX stderr

	public	:: LEN_FILENAME

! !REVISION HISTORY:
!	10oct96 - Jing G.	- Defined
!       25Jul02 - J. Larson     - Changed cpp define token HP-UX to
!                                 HP_UX for compatibility with Fujitsu
!                                 cpp.
!EOP
!_______________________________________________________________________

!    Defines standar i/o units.

	integer, parameter :: stdin  = 5
	integer, parameter :: stdout = 6

#ifdef	sysHP_UX
		! Special setting for HP-UX

	integer, parameter :: stderr = 6
#else
		! Generic setting for UNIX other than HP-UX

	integer, parameter :: stderr = 6
#endif

	integer, parameter :: LEN_FILENAME = 128

!-----------------------------------------------------------------------
end module m_stdio
!.

!! This module defines the prescision used for real variables within the CARMA code.
!! It uses defintion for a real that is used within CAM.
!!
!! NOTE: This file is adapted for use within CAM and is different than the file of the
!! same name that is part of the standard CARMA distribution.
module carma_precision_mod

   use shr_kind_mod, only: shr_kind_r8

	implicit none
	
	integer, parameter      :: f = shr_kind_r8
	real(kind=f), parameter :: powmax = 308._f
	
	! Precision control strategy
	! JAS CU-Boulder June 8, 2006
	!
	! I imagine using these statements bracketed with some CPP statements
	! to control the overall precision of a model.  All variables would be
	! declared as real(f).  All physical constants would have a
	! a suffix of _f, e.g. 2._f, to force them into the proper precision.
	!
	! I do wonder if it would be more accurate to declare variables as
	! real( kind=f ), but real(f) is how Chivers and Sleightholme
	! declare in their F90 text.
	!
	! Both real(f) and real( kind=f ) seem to work, but I'm more comfortable
	! with real( kind=f ), so I'm using that in all declarations.
	
	!--
	! Numerical constants
	!! Define 1 in the specified precision.
	real(kind=f), parameter :: ONE = 1._f
	
	!!  Define smallest possible number such that ONE + ALMOST_ZERO > ONE
	real(kind=f), parameter :: ALMOST_ZERO = epsilon( ONE )
	real(kind=f), parameter :: ALMOST_ONE  = ONE - ALMOST_ZERO
end module

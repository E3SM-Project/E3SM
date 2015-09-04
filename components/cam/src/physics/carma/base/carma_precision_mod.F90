module carma_precision_mod

	implicit none
	
#ifdef SINGLE
	
	! For floats commonly referred to as 'real'
	! -at least 6 places of precision past the decimal
	! -must span at least 10**(-37) to 10**(37) 
	integer, parameter      :: f = selected_real_kind(6,37)
	real(kind=f), parameter :: powmax = 85._f
	
#else
	
	! For floats commonly referred to as 'double precision'
	! -at least 15 places of precision past the decimal
	! -must span at least 10**(-307) to 10**(307) 
	
	integer, parameter      :: f = selected_real_kind(15,307)
	real(kind=f), parameter :: powmax = 706._f
	
#endif
	
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

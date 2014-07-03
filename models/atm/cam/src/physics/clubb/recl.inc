!-------------------------------------------------------------------------------
! $Id: recl.inc 5623 2012-01-17 17:55:26Z connork@uwm.edu $
! Description:
!   Preprocessing rules for determining how large an unformatted 
!   data record is when using Fortran write.  This does not affect
!   netCDF output at all.

! Notes:
!   New directives will need to be added to port CLUBB GrADS output
!   to new compilers that do not use byte size record lengths.

!   Early Alpha processors lacked the ability to work with anything
!   smaller than a 32 bit word, so DEC Fortran and its successors
!   (Compaq Visual Fortran, newer Intel Fortran, etc.) all use 4 
!   byte records.  Note that specifying byterecl on Alpha still 
!   results in a performance hit, even on newer chips.
!-------------------------------------------------------------------------------
#if defined GFDL  /* F_RECL should be 4 for the GFDL SCM-CLUBB */
#  define F_RECL 4
#elif defined __INTEL_COMPILER && __INTEL_COMPILER >= 800 /* Versions of Intel fortran > 8.0_core_rknd */
#  define F_RECL 1
#elif defined(__alpha) /* Assume 4 byte word on Alpha processors */
#  define F_RECL 1
#else
#  define F_RECL 4 /* Most compilers and computers */
#endif

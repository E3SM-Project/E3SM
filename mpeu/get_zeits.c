/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
! CVS $Id$
! CVS $Name$
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: get_zeits - a C interface to times for Fortran calls
!
! !DESCRIPTION:
!
! !INTERFACE: */
 /*
  System times() dependencies:
 */


#include <sys/types.h>
#ifndef SYSCATAMOUNT /* times is not implemented and will always fail on Catamount */ 
#include <sys/times.h>
#endif

#include <time.h> /* POSIX standard says CLOCKS_PER_SEC is here */

/*
 *  CLK_TCK is obsolete - replace with CLOCKS_PER_SEC
 */

#define ZCLK_TCK ((double)CLOCKS_PER_SEC)



 /*
  The default is FORTRAN_UNDERSCORE_, but not explicitly used.
 */

#ifdef FORTRAN_CAPS_
#  define	get_zeits_		GET_ZEITS
#  define	get_ztick_		GET_ZTICK
#endif

#ifdef FORTRAN_SAME
#  define	get_zeits_		get_zeits
#  define	get_ztick_		get_ztick
#endif

#ifdef FORTRAN_GNUF2C
#  define	get_zeits_		get_zeits__
#  define	get_ztick_		get_ztick__
#endif

 /*  Prototype: */

   void get_zeits_(double *zts);
   void get_ztick_(double *tic);

/*!REVISION HISTORY:
! 	12Mar98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! 	06Jul99 - J.W. Larson <jlarson@dao> - support for AIX platform
!EOP */

/*  Implementations: */

void get_zeits_(zts)
  double *zts;
{

#ifndef SYSCATAMOUNT
  struct tms tm;
  double secs;
  secs=1./ZCLK_TCK;

  zts[0]=times(&tm)*secs;
  zts[1]=tm.tms_utime*secs;
  zts[2]=tm.tms_stime*secs;
  zts[3]=tm.tms_cutime*secs;
  zts[4]=tm.tms_cstime*secs;
#else
  zts[0]=0.;
  zts[1]=0.;
  zts[2]=0.;
  zts[3]=0.;
  zts[4]=0.;
#endif

}

void get_ztick_(tic)
  double *tic;
{
  tic[0]=1./ZCLK_TCK;
}


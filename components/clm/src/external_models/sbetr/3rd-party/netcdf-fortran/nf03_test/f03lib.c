/*
 * $Id: fort03.c,v 1.9 2007/07/28 13:14:00 rick Exp $
 *
 * This file contains support functions for FORTRAN code.  For example,
 * under HP-UX A.09.05, the U77 library doesn't contain the exit()
 * routine -- so we create one here.
 */

/*
   Modified fortlib.c - We remove all cfortran.h stuff to make 
   compiling easier. Also functions are modified to be extern
   and not static so FORTRAN can see them

   Version 1. July  2007 first vesion 
   Version 2. April 2009 - modified for netCDF 4.0.1
      
   Modified by: Richard Weed, Ph.D
   Center for Advanced Vehicular Systems
   Misssissippi State University
   rweed@.cavs.msstate.edu
*/

#include <stdlib.h>
#include <limits.h>
#include <float.h>


extern double
myrand(int iflag)
{
    if (iflag != 0)
	srand(iflag);

    /*
     * Return a pseudo-random value between 0.0 and 1.0.
     *
     * We don't use RAND_MAX here because not all compilation
     * environments define it (e.g. gcc(1) under SunOS 4.1.3).
     */
    return (rand() % 32768) / 32767.0;
}


extern int
myshift(int value, int amount)
{
    if (amount < 0)
	value >>= -amount;
    else
    if (amount > 0)
	value <<= amount;
    return value;
}

#include <signal.h>
extern void
nc_ignorefpe(int doit)
{
	if(doit)
		(void) signal(SIGFPE, SIG_IGN);
}

extern double cmax_uchar()
{
    return UCHAR_MAX;
}

extern double cmin_schar()
{
    return SCHAR_MIN;
}

extern double cmax_schar()
{
    return SCHAR_MAX;
}

extern double cmin_short()
{
    return SHRT_MIN;
}

extern double cmax_short()
{
    return SHRT_MAX;
}

extern double cmin_int()
{
    return INT_MIN;
}

extern double cmax_int()
{
    return INT_MAX;
}

extern double cmin_long()
{
    return LONG_MIN;
}

extern double cmax_long()
{
    return LONG_MAX;
}

extern double cmax_float()
{
    return FLT_MAX;
}

extern double cmax_double()
{
    return DBL_MAX;
}

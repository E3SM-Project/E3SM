/* =============================================================================
** SVN $Id: shr_vmath_fwrap.c 239 2006-02-08 19:02:33Z kauff $
** SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070423/shr/shr_vmath_fwrap.c $
** =============================================================================
** Fortran wrappers for shr_vmath calls for systems that only
** provide a C interface to the system vector math routines.
** =============================================================================
*/

#if (defined IRIX64)

void shr_vmath_fwrap_vsqrt_(double *X, double *Y, int *n)
{
   vsqrt(X, Y, *n, 1, 1);
}

void shr_vmath_fwrap_vexp_(double *X, double *Y, int *n)
{
   vexp(X, Y, *n, 1, 1);
}

void shr_vmath_fwrap_vlog_(double *X, double *Y, int *n)
{
   vlog(X, Y, *n, 1, 1);
}

void shr_vmath_fwrap_vsin_(double *X, double *Y, int *n)
{
   vsin(X, Y, *n, 1, 1);
}

void shr_vmath_fwrap_vcos_(double *X, double *Y, int *n)
{
   vcos(X, Y, *n, 1, 1);
}

#endif

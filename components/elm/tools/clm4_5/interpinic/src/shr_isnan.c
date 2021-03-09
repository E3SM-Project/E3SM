#include <math.h>
#include <shr_isnan.h>

/*
  Only define the "C" function for isnan if a FORTRAN version is NOT available
*/
#ifdef NOFTN_INTRINSIC

/*
** Use the MCT logic to figure out FORTRAN name mangling
*/
#if   defined(FORTRAN_UNDERSCORE_) || defined(FORTRANUNDERSCORE)
#define FORT_NAME(lower,upper) lower##_
#elif defined(FORTRAN_GNUF2C)
#define FORT_NAME(lower,upper) lower##__
#elif defined(FORTRAN_SAME)
#define FORT_NAME(lower,upper) lower
#elif defined(FORTRAN_CAPS_)
#define FORT_NAME(lower,upper) upper
#else
#error "Unrecognized Fortran-mangle type"
/* set to something reasonable to avoid cascade of cc errors */
#define FORT_NAME(lower,upper) lower##_
#endif

/*
** This uses the "C" isnan function and makes it callable from FORTRAN
*/
int FORT_NAME( shr_sisnan , SHR_SISNAN )(float *real )
/*
* Check if the input single precisioun float value is a NaN or not
*/
{
  return( isnanf( *real ) );
}

int FORT_NAME( shr_disnan , SHR_DISNAN )(double *real )
/*
* Check if the input double precisioun float value is a NaN or not
*/
{
  return( isnan( *real ) );
}

#endif

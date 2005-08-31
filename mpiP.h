
/*
 * Private .h file for MPI
 */


#include <malloc.h>
#include <stdio.h>

#include "listops.h"
#include "mpi.h"


/*
 * Fortran name mangling
 *
 * the configure.ac specifies these
 *
 * cpp does not have the ability to change the case
 * of the argument, so the invocation of the macro
 * has to be give both e.g. FORT_NAME(hello,HELLO)
 * and maps to "hello_", "hello", and "HELLO" repectively.
 */


#if   defined(FORTRAN_UNDERSCORE_)
#define FORT_NAME(lower,upper) lower##_
#elif defined(FORTRAN_SAME)
#define FORT_NAME(lower,upper) lower
#elif defined(FORTRAN_CAPS)
#define FORT_NAME(lower,upper) upper
#else
#define FORT_NAME(lower,upper) FORTRAN_MANGLE_ERROR
#error "Unrecognized Fortran-mangle type"
#endif



/****************************************************************************/


typedef struct
{
  pList sendlist;
  pList recvlist;

  int num;

} Comm;



typedef struct
{
  pListitem listitem;        /* to allow Req to be removed from list */

  int *buf;
  int tag;
  int complete;

} Req;






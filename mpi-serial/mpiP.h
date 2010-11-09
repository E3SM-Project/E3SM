
/*
 * Private .h file for MPI
 */


#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

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
 *
 * IMPORTANT NOTE:
 * In the case of FORTRAN_GNUF2C (e.g. g95), the rule is this:
 *    name does not contain an underscore -> append *one* underscore
 *    name contains an underscore -> append *two* underscore
 * Since all the mpi-serial names exported to fortran start with "mpi_",
 * we only support the latter.
 *
 * Note: FORTRANUNDERSCORE is needed by ccsm
 *
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
 * MPI_GROUP_ONE must not conflict with MPI_GROUP_NULL or
 * MPI_GROUP_EMPTY
 */

#define MPI_GROUP_ONE  (1)


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
  int source;
  int tag;
  int complete;

} Req;




/****************************************************************************/


extern void *mpi_malloc(int size);
extern void mpi_free(void *ptr);

extern MPI_Comm mpi_comm_new(void);

extern void mpi_destroy_handles(void);
extern void mpi_alloc_handle(int *handle, void **data);
extern void *mpi_handle_to_ptr(int handle);
extern void mpi_free_handle(int handle);


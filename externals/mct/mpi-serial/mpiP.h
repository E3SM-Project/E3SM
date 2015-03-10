
/*
 * Private .h file for MPI
 */


#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "listops.h"
#include "mpi.h"
#include "config.h"

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


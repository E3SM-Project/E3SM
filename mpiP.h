
/*
 * Private .h file for MPI
 */


#include <malloc.h>
#include <stdio.h>

#include "listops.h"
#include "mpi.h"


#define FORT_NAME(x) x##_


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






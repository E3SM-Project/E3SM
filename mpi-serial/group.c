
#include "mpiP.h"


/*********/


FORT_NAME( mpi_group_incl, MPI_GROUP_INCL )
     (int *group, int *n, int *ranks, int *newgroup, int *ierror)
{
  *ierror= MPI_Group_incl(*group, *n, ranks, newgroup);
}


int MPI_Group_incl(MPI_Group group, int n, int *ranks, MPI_Group *newgroup)
{

  if (group==MPI_GROUP_NULL || group==MPI_GROUP_EMPTY)
    *newgroup=group;
  else
    if (n==1 && ranks[0]==0)
      *newgroup=MPI_GROUP_ONE;
    else
      *newgroup=MPI_GROUP_NULL;


  return(MPI_SUCCESS);
}


/*********/


FORT_NAME( mpi_group_free, MPI_GROUP_FREE )(int *group, int *ierror)
{
  *ierror= MPI_Group_free(group);
}


int MPI_Group_free(MPI_Group *group)
{
  *group= MPI_GROUP_NULL;

  return(MPI_SUCCESS);
}

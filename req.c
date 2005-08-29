#include "mpiP.h"


/*
 * COMPLETION
 */



FORT_NAME( mpi_test )(int *request, int *flag, int *status, int *ierror)
{
  *ierror=MPI_Test( (void *)request ,flag,(MPI_Status *)status);
}



int MPI_Test(MPI_Request *request, int *flag, MPI_Status *status)
{
  Req *req;

  if (*request==MPI_REQUEST_NULL)
    {
      status->MPI_TAG= MPI_ANY_TAG;
      status->MPI_SOURCE= MPI_ANY_SOURCE;
      *flag=1;
      return(MPI_SUCCESS);
    }


  req=mpi_handle_to_ptr(*request);

  *flag=req->complete;

  if (*flag)
    {
      status->MPI_SOURCE= 0;
      status->MPI_TAG= req->tag;

      mpi_free_handle(*request);
      *request=MPI_REQUEST_NULL;
    }

  return(MPI_SUCCESS);
}



FORT_NAME( mpi_wait )(int *request, int *status, int *ierror)
{
  *ierror=MPI_Wait( (void *)request, (MPI_Status *)status );
}



int MPI_Wait(MPI_Request *request, MPI_Status *status)
{
  int flag;

  MPI_Test(request,&flag,status);

  if (!flag)
    {
      fprintf(stderr,"MPI_Wait: request not complete, deadlock\n");
      abort();
    }

  return(MPI_SUCCESS);
}


/*********/


FORT_NAME( mpi_waitany )(int *count, int *requests,
			 int *index, int *status, int *ierror)
{

  *ierror=MPI_Waitany(*count, (void *)requests,index,(MPI_Status *)status);
}



int MPI_Waitany(int count, MPI_Request *array_of_requests,
		int *index, MPI_Status *status)
{
  int i;
  int flag;

  for (i=0; i<count; i++)
    {
      MPI_Test(&array_of_requests[i],&flag,status);
      
      if (flag)
	return(MPI_SUCCESS);
    }

  /* none are completed */

  fprintf(stderr,"MPI_Waitany: no requests complete, deadlock\n");
  abort();
}


/*********/


FORT_NAME( mpi_waitall )(int *count, int *array_of_requests,
			 int *array_of_statuses, int *ierror)
{
  *ierror=MPI_Waitall(*count, (void *)array_of_requests,
		      (MPI_Status *)array_of_statuses);

}



int MPI_Waitall(int count, MPI_Request *array_of_requests,
		MPI_Status *array_of_statuses)
{
  int i;
  int flag;

  for (i=0; i<count; i++)
    {
      MPI_Test(&array_of_requests[i],&flag,&array_of_statuses[i]);

      if (!flag)
	{
	  fprintf(stderr,"MPI_Waitall: request not complete, deadlock\n");
	  abort();
	}
    }

  return(MPI_SUCCESS);
}




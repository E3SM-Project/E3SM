

#include "mpiP.h"


/****************************************************************************/

static int mpi_initialized=0;


/****************************************************************************/


void *mpi_malloc(int size)
{
  void *ret;

  ret=malloc(size);

  if (!ret)
    {
      fprintf(stderr,"mpi_malloc: failed to allocate %d bytes\n",size);
      abort();
    }
      
  return(ret);
}


void mpi_free(void *ptr)
{
  free(ptr);
}


/*
 * Communicators
 *
 */



static MPI_Comm mpi_comm_new(void)
{
  MPI_Comm chandle;
  Comm *cptr;
  static int num=0;

  mpi_alloc_handle(&chandle,(void **) &cptr);

  cptr->sendlist=AP_list_new();
  cptr->recvlist=AP_list_new();

  cptr->num=num++;

  return(chandle);
}


/*********/


FORT_NAME( mpi_comm_free )(int *comm, int *ierror)
{
  *ierror=MPI_Comm_free(comm);
}


/*
 * MPI_Comm_free()
 *
 * Note: will NOT free any pending MPI_Request handles
 * that are allocated...  correct user code should have
 * already done a Wait or Test to free them.
 *
 */


int MPI_Comm_free(MPI_Comm *comm)
{
  pList sendlist, recvlist;
  int size;
  Comm *mycomm;

  mycomm=mpi_handle_to_ptr(*comm);   /* (Comm *)(*comm) */

  sendlist=mycomm->sendlist;
  recvlist=mycomm->recvlist;

  size=AP_list_size(sendlist);
  if (size!=0)
    fprintf(stderr,"MPI_Comm_free: warning: %d pending send reqs\n",
	    size);
  AP_list_free(sendlist);


  size=AP_list_size(recvlist);
  if (size!=0)
    fprintf(stderr,"MPI_Comm_free: warning: %d pending receive reqs\n",
	    size);
  AP_list_free(recvlist);

  mpi_free_handle(*comm);            /* free(mycomm); */
  *comm=MPI_COMM_NULL;

  return(MPI_SUCCESS);
}


/*********/


FORT_NAME( mpi_comm_dup )(int *comm, int *newcomm, int *ierror)
{

  *ierror=MPI_Comm_dup( *comm, newcomm);

}


int MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm)
{
  *newcomm= mpi_comm_new();

#ifdef INFO
  fflush(stdout);
  fprintf(stderr,"MPI_Comm_dup: new comm handle=%d\n",*newcomm);
#endif

  return(MPI_SUCCESS);
}


/*********/



FORT_NAME( mpi_comm_size )(int *comm, int *size, int *ierror)
{
  *ierror=MPI_Comm_size(*comm, size);
}



int MPI_Comm_size(MPI_Comm comm, int *size)
{
  *size=1;

  return(MPI_SUCCESS);
}


/*********/


FORT_NAME( mpi_comm_rank )(int *comm, int *rank, int *ierror)
{
  *ierror=MPI_Comm_rank( *comm, rank);
}


int MPI_Comm_rank(MPI_Comm comm, int *rank)
{
  *rank=0;

  return(MPI_SUCCESS);
}



/**************************************************************************/

/*
 * INIT/FINALIZE
 *
 */



FORT_NAME( mpi_init_fort )(int *f_MPI_COMM_WORLD,
                           int *f_MPI_ANY_SOURCE, int *f_MPI_ANY_TAG,
                           int *f_MPI_COMM_NULL, int *f_MPI_REQUEST_NULL,
                           int *f_MPI_MAX_ERROR_STRING, 
                           int *f_MPI_STATUS_SIZE, 
                           int *f_MPI_SOURCE, int *f_MPI_TAG, int *f_MPI_ERROR,
			   int *f_status,
			   int *fsource, int *ftag, int *ferror,
                           int *f_MPI_INTEGER, void *fint1, void *fint2,
                           int *f_MPI_LOGICAL, void *flog1, void *flog2,
                           int *f_MPI_REAL, void *freal1, void *freal2,
                           int *f_MPI_DOUBLE_PRECISION,
			   void *fdub1, void *fdub2,
			   int *f_MPI_COMPLEX, void *fcomp1, void *fcomp2,
                           int *ierror)
{
  int err;
  int size;
  int offset;

  *ierror=MPI_Init(NULL,NULL);

  err=0;

#define verify_eq(name)  \
  if (*f_##name != name) \
    { fprintf(stderr,"mpi_init_fort: mpif.h %s (%d) " \
                     "is not consistant with mpi.h value (%d)\n", \
                     #name,*f_##name,name); \
      err=1; }

#define verify_size(name,p1,p2) \
  if ( (size=((char *)(p2) - (char *)(p1))) != *f_##name ) \
    { fprintf(stderr,"mpi_init_fort: mpif.h %s (%d) " \
                     "does not match fortran size %d\n", \
                     #name,*f_##name,size); \
      err=1; }

#define verify_field(name) \
  { offset= (char *)&((MPI_Status *)f_status)->name - (char *)f_status; \
    if ( offset != (*f_##name-1)*sizeof(int) ) \
    { fprintf(stderr,"mpi_init_fort: mpif.h %s (%d) (%d bytes) " \
                     "is inconsistant w/offset in MPI_Status (%d bytes)\n", \
                    #name,*f_##name,(*f_##name-1)*sizeof(int),offset); \
      err=1; }}


  verify_eq(MPI_COMM_WORLD);
  verify_eq(MPI_ANY_SOURCE);
  verify_eq(MPI_ANY_TAG);
  verify_eq(MPI_COMM_NULL);
  verify_eq(MPI_REQUEST_NULL);
  verify_eq(MPI_MAX_ERROR_STRING);

  verify_eq(MPI_STATUS_SIZE);
  verify_field(MPI_SOURCE);
  verify_field(MPI_TAG);
  verify_field(MPI_ERROR);

  verify_eq(MPI_INTEGER);
  verify_size(MPI_INTEGER,fint1,fint2);

  verify_size(MPI_LOGICAL,flog1,flog2);

  verify_eq(MPI_REAL);
  verify_size(MPI_REAL,freal1,freal2);

  verify_eq(MPI_DOUBLE_PRECISION);
  verify_size(MPI_DOUBLE_PRECISION,fdub1,fdub2);

  verify_size(MPI_COMPLEX,fcomp1,fcomp2);

  if (err)
    abort();
}



int MPI_Init(int *argc, char **argv[]) 
{
  MPI_Comm my_comm_world;

  my_comm_world=mpi_comm_new();

  if (my_comm_world != MPI_COMM_WORLD)
    {
      fprintf(stderr,"MPI_Init: conflicting MPI_COMM_WORLD\n");
      abort();
    }

  mpi_initialized=1;
  return(MPI_SUCCESS);
}


/*********/


FORT_NAME(mpi_finalize)(int *ierror)
{
  *ierror=MPI_Finalize();
}


/*
 * MPI_Finalize()
 *
 * this library doesn't support re-initializing MPI, so
 * the finalize will just leave everythign as it is...
 *
 */


int MPI_Finalize(void)
{
  mpi_initialized=0;

  mpi_destroy_handles();

  return(MPI_SUCCESS);
}


/*********/


FORT_NAME( mpi_abort )(int *comm, int *errorcode, int *ierror)
{
  *ierror=MPI_Abort( *comm, *errorcode);
}



int MPI_Abort(MPI_Comm comm, int errorcode)
{
  fprintf(stderr,"MPI_Abort: error code = %d\n",errorcode);
  exit(errorcode);
}


/*********/



FORT_NAME( mpi_error_string )(int *errorcode, char *string,
			      int *resultlen, int *ierror)
{
  *ierror=MPI_Error_string(*errorcode, string, resultlen);
}


int MPI_Error_string(int errorcode, char *string, int *resultlen)
{
  sprintf(string,"MPI Error: code %d\n",errorcode);
  *resultlen=strlen(string);

  return(MPI_SUCCESS);
}


/*********/


FORT_NAME( mpi_initialized )(int *flag, int *ierror)
{
  *ierror=MPI_Initialized(flag);
}


int MPI_Initialized(int *flag)
{
  *flag= mpi_initialized;

  return(MPI_SUCCESS);
}


/***********************************************************************/


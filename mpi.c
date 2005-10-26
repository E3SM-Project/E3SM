

#include "mpiP.h"


/****************************************************************************/

static int initialized=0;


/****************************************************************************/


/*
 * INIT/FINALIZE
 *
 */



FORT_NAME( mpi_init_fort , MPI_INIT_FORT)
                          (int *f_MPI_COMM_WORLD,
                           int *f_MPI_ANY_SOURCE, int *f_MPI_ANY_TAG,
                           int *f_MPI_COMM_NULL, int *f_MPI_REQUEST_NULL,
			   int *f_MPI_GROUP_NULL, int *f_MPI_GROUP_EMPTY,
			   int *f_MPI_UNDEFINED,
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

  /*
   * These 3 macros compare things from mpif.h (as passed in by the f_
   * arguments) to the values in C (from #including mpi.h).
   *
   * Unfortunately, this kind of thing is done most easily in a nasty
   * looking macto.
   *
   */


  /*
   * verify_eq
   *   compare value of constants in C and fortran
   *   i.e. compare *f_<name> to <name>
   */

#define verify_eq(name)  \
  if (*f_##name != name) \
    { fprintf(stderr,"mpi_init_fort: mpif.h %s (%d) " \
                     "is not consistant with mpi.h value (%d)\n", \
                     #name,*f_##name,name); \
      err=1; }

  /*
   * verify_size
   *   verify that the type name in fortran has the correct
   *   value (i.e. the size of that data type).
   *   Determine size by subtracting the pointer values of two
   *   consecutive array locations.
   */

#define verify_size(name,p1,p2) \
  if ( (size=((char *)(p2) - (char *)(p1))) != *f_##name ) \
    { fprintf(stderr,"mpi_init_fort: mpif.h %s (%d) " \
                     "does not match fortran size (%d)\n", \
                     #name,*f_##name,size); \
      err=1; }

  /*
   * verify_field
   *   check the struct member offsets for MPI_Status vs. the
   *   fortan integer array offsets.  E.g. the location of
   *   status->MPI_SOURCE should be the same as STATUS(MPI_SOURCE)
   */

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
  verify_eq(MPI_GROUP_NULL);
  verify_eq(MPI_GROUP_EMPTY);
  verify_eq(MPI_UNDEFINED);
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

  initialized=1;
  return(MPI_SUCCESS);
}


/*********/


FORT_NAME( mpi_finalize, MPI_FINALIZE )(int *ierror)
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
  initialized=0;

  mpi_destroy_handles();

  return(MPI_SUCCESS);
}


/*********/


FORT_NAME( mpi_abort , MPI_ABORT )(int *comm, int *errorcode, int *ierror)
{
  *ierror=MPI_Abort( *comm, *errorcode);
}



int MPI_Abort(MPI_Comm comm, int errorcode)
{
  fprintf(stderr,"MPI_Abort: error code = %d\n",errorcode);
  exit(errorcode);
}


/*********/



FORT_NAME( mpi_error_string , MPI_ERROR_STRING)
                             (int *errorcode, char *string,
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


FORT_NAME( mpi_initialized , MPI_INITIALIZED )(int *flag, int *ierror)
{
  *ierror=MPI_Initialized(flag);
}


int MPI_Initialized(int *flag)
{
  *flag= initialized;

  return(MPI_SUCCESS);
}




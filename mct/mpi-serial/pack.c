
#include "mpiP.h"

/*
 *
 */


FORT_NAME( mpi_pack , MPI_PACK )
     ( void *inbuf, int *incount, int *datatype,
       void *outbuf, int *outsize, int *position, int *comm, int *ierror)
{
  *ierror=MPI_Pack(inbuf, *incount,* datatype,
		   outbuf, *outsize, position, *comm);
}



int MPI_Pack( void *inbuf, int incount, MPI_Datatype datatype,
              void *outbuf, int outsize, int *position, MPI_Comm comm)
{
  int size;

  size=datatype*incount;

  if ( (*position)+size > outsize)
    {
      fprintf(stderr,"MPI_Pack: ran out of space in outbuf\n");
      abort();
    }

  memcpy( (char *)outbuf+(*position), inbuf, size);
  (*position)+=size;

  return(MPI_SUCCESS);
}



/*
 *
 */


FORT_NAME( mpi_unpack , MPI_UNPACK )
     ( void *inbuf, int *insize, int *position,
       void *outbuf, int *outcount, int *datatype,
       int *comm, int *ierror )
{
  *ierror=MPI_Unpack( inbuf, *insize, position,
                      outbuf, *outcount, *datatype, *comm);
}



int MPI_Unpack( void *inbuf, int insize, int *position,
                void *outbuf, int outcount, MPI_Datatype datatype,
                MPI_Comm comm )
{
  int size;

  size=datatype*outcount;

  if ( (*position)+size > insize )
    {
      fprintf(stderr,"MPI_Unpack: ran out of data in inbuf\n");
      abort();
    }


  memcpy(outbuf, (char *)inbuf+(*position) , size);
  (*position)+=size;

  return(MPI_SUCCESS);

}

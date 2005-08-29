
#include "mpiP.h"


/*
 * SENDING
 *
 */



static int mpi_match_recv(void *r, void *tag)
{
  return( ((Req *)r)->tag == MPI_ANY_TAG ||
	  ((Req *)r)->tag == *((int *)tag) );
}


/*
 *
 */



FORT_NAME( mpi_isend )(void *buf, int *count, int *datatype,
		       int *dest, int *tag, int *comm, int *req, int *ierror)
{

  *ierror=MPI_Isend(buf,*count,*datatype,*dest,*tag,
		    *comm, (void *)req);

}



int MPI_Isend(void *buf, int count, MPI_Datatype datatype,
	      int dest, int tag, MPI_Comm comm, MPI_Request *request) 
{
  pListitem match;
  Comm *mycomm;
  Req *rreq, *sreq;

  mycomm=mpi_handle_to_ptr(comm);   /* (Comm *)comm; */

#ifdef INFO
  fflush(stdout);
  fprintf(stderr,"MPI_Isend: Comm=%d  tag=%d  count=%d type=%d\n",
	 mycomm->num,tag,count,datatype);
#endif

  if (dest!=0)
    {
      fprintf(stderr,"MPI_Isend: send to %d\n",dest);
      abort();
    }

  mpi_alloc_handle(request,(void **) &sreq);


  if ( match=AP_list_search_func(mycomm->recvlist,mpi_match_recv,&tag) )
    {
      rreq=(Req *)AP_listitem_data(match);
      AP_list_delete_item(mycomm->recvlist,match);

      memcpy(rreq->buf,buf,count * datatype);
      rreq->complete=1;
      rreq->tag=tag;                    /* in case rreq->tag was MPI_ANY_TAG */

      sreq->complete=1;

#ifdef DEBUG
      printf("Completion(send) value=%d tag=%d\n",
	     *((int *)buf),rreq->tag);
#endif

      return(MPI_SUCCESS);
    }

  sreq->buf=buf;
  sreq->tag=tag;
  sreq->complete=0;
  sreq->listitem=AP_list_append(mycomm->sendlist,sreq);

#ifdef INFO
  print_list(mycomm->sendlist,"sendlist for comm ",mycomm->num);
#endif

  return(MPI_SUCCESS);
}


/*********/


FORT_NAME(mpi_send) (void *buf, int *count, int *datatype,
		     int *dest, int *tag, int *comm, int *ierror)
{
  *ierror=MPI_Send(buf, *count, *datatype, *dest, *tag, *comm);
}



int MPI_Send(void* buf, int count, MPI_Datatype datatype,
	     int dest, int tag, MPI_Comm comm)
{
  MPI_Request request;
  MPI_Status status;
  int flag;

#ifdef INFO
  fflush(stdout);
  fprintf(stderr,"MPI_Send: ");
#endif

  MPI_Isend(buf,count,datatype,dest,tag,comm,&request);
  MPI_Wait(&request,&status);


  return(MPI_SUCCESS);
}




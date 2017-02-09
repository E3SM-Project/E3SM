/**
 * @file pio_spmd.c
 * @author Algorithms modeled after spmd_utils in the Community Atmosphere Model; C translation Jim Edwards
 * @date 2014
 * @brief MPI_Gather, MPI_Gatherv, and MPI_Alltoallw with flow control options
 */
 
#ifdef TESTSWAPM
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#define PIO_NOERR 0
#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })
#define MAX_GATHER_BLOCK_SIZE 32
#else
#include <pio.h>
#include <pio_internal.h>
#endif

/** 
 ** @brief Wrapper for MPI calls to print the Error string on error
 */
void CheckMPIReturn(const int ierr,const char file[],const int line)
{
  
  if(ierr != MPI_SUCCESS){
    char errstring[MPI_MAX_ERROR_STRING];
    int errstrlen;
    int mpierr = MPI_Error_string( ierr, errstring, &errstrlen);

    fprintf(stderr, "MPI ERROR: %s in file %s at line %d\n",errstring, file, line);
    
  }
}


/**
 **  @brief Provides the functionality of MPI_Gather with flow control options
 */

int pio_fc_gather( void *sendbuf, const int sendcnt, const MPI_Datatype sendtype,
		   void *recvbuf, const int recvcnt, const MPI_Datatype recvtype, const int root, 
		   const MPI_Comm comm, const int flow_cntl)
{
  bool fc_gather;
  int gather_block_size;
  int mytask, nprocs;
  int mtag;
  MPI_Status status;
  int ierr;
  int hs;
  int displs;
  int dsize;
  


  if(flow_cntl > 0){
    fc_gather = true;
    gather_block_size = min(flow_cntl,MAX_GATHER_BLOCK_SIZE);
  }else{
    fc_gather = false;
  }

  if(fc_gather){
    CheckMPIReturn(MPI_Comm_rank (comm, &mytask), __FILE__,__LINE__);
    CheckMPIReturn(MPI_Comm_size (comm, &nprocs), __FILE__,__LINE__);

    mtag = 2*nprocs;
    hs = 1;

    if(mytask == root){
      int preposts = min(nprocs-1, gather_block_size);
      int head=0;
      int count=0;
      int tail = 0;
      MPI_Request rcvid[gather_block_size];

      CheckMPIReturn(MPI_Type_size(recvtype, &dsize), __FILE__,__LINE__);

      for(int p=0;p<nprocs;p++){
	if(p != root){
	  if(recvcnt > 0){
	    count++;
	    if(count > preposts){
	      CheckMPIReturn(MPI_Wait(rcvid+tail, &status), __FILE__,__LINE__);
	      tail = (tail+1) % preposts;
	    }
	    displs = p*recvcnt*dsize;

	    char *ptr = (char *) recvbuf + displs;

	    CheckMPIReturn(MPI_Irecv(  ptr, recvcnt, recvtype, p, mtag, comm, rcvid+head), __FILE__,__LINE__);
	    head= (head+1) % preposts;
	    CheckMPIReturn(MPI_Send( &hs, 1, MPI_INT, p, mtag, comm), __FILE__,__LINE__);
	  }
	}
      }

      // copy local data
      CheckMPIReturn(MPI_Type_size(sendtype, &dsize), __FILE__,__LINE__);      
      memcpy(recvbuf, sendbuf, sendcnt*dsize );

      count = min(count, preposts);

      if(count>0){
	CheckMPIReturn(MPI_Waitall( count, rcvid, MPI_STATUSES_IGNORE),__FILE__,__LINE__);
      }
    }else{
      if(sendcnt > 0){
	CheckMPIReturn(MPI_Recv( &hs, 1, MPI_INT, root, mtag, comm, &status), __FILE__,__LINE__);
	CheckMPIReturn(MPI_Send( sendbuf, sendcnt, sendtype, root, mtag, comm), __FILE__,__LINE__);
      }
    }
  }else{
    //    printf("%s %d %ld %d %ld %d %d\n",__FILE__,__LINE__,sendbuf,sendcnt,recvbuf,recvcnt,root);
    CheckMPIReturn(MPI_Gather ( sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm), __FILE__,__LINE__);
  }

  return PIO_NOERR;
}
  


/**
 **  @brief Provides the functionality of MPI_Gatherv with flow control options
 */

int pio_fc_gatherv( void *sendbuf, const int sendcnt, const MPI_Datatype sendtype,
		    void *recvbuf, const int recvcnts[], const int displs[],
		    const MPI_Datatype recvtype, const int root, 
		    const MPI_Comm comm, const int flow_cntl)
{
  bool fc_gather;
  int gather_block_size;
  int mytask, nprocs;
  int mtag;
  MPI_Status status;
  int ierr;
  int hs;
  int dsize;
  

  if(flow_cntl > 0){
    fc_gather = true;
    gather_block_size = min(flow_cntl,MAX_GATHER_BLOCK_SIZE);
  }else{
    fc_gather = false;
  }

  if(fc_gather){
    CheckMPIReturn(MPI_Comm_rank (comm, &mytask), __FILE__,__LINE__);
    CheckMPIReturn(MPI_Comm_size (comm, &nprocs), __FILE__,__LINE__);

    mtag = 2*nprocs;
    hs = 1;

    if(mytask == root){
      int preposts = min(nprocs-1, gather_block_size);
      int head=0;
      int count=0;
      int tail = 0;
      MPI_Request rcvid[gather_block_size];
      // printf("%s %d %d\n",__FILE__,__LINE__,(int) recvtype);
      CheckMPIReturn(MPI_Type_size(recvtype, &dsize), __FILE__,__LINE__);

      for(int p=0;p<nprocs;p++){
	if(p != root){
	  if(recvcnts[p] > 0){
	    count++;
	    if(count > preposts){
	      CheckMPIReturn(MPI_Wait(rcvid+tail, &status), __FILE__,__LINE__);
	      tail = (tail+1) % preposts;
	    }

	    void *ptr = (void *)((char *) recvbuf + dsize*displs[p]);

	    //	  printf("%s %d %d %d\n",__FILE__,__LINE__,p,(int) recvtype);
	    CheckMPIReturn(MPI_Irecv(  ptr, recvcnts[p], recvtype, p, mtag, comm, rcvid+head), __FILE__,__LINE__);
	    head= (head+1) % preposts;
	    CheckMPIReturn(MPI_Send( &hs, 1, MPI_INT, p, mtag, comm), __FILE__,__LINE__);
	  }
	}
      }

      // copy local data
      CheckMPIReturn(MPI_Type_size(sendtype, &dsize), __FILE__,__LINE__);      
      CheckMPIReturn(MPI_Sendrecv(sendbuf, sendcnt, sendtype,
				  mytask, 102, recvbuf, recvcnts[mytask], recvtype,
				  mytask, 102, comm, &status),__FILE__,__LINE__);

      count = min(count, preposts);
      if(count>0)
	CheckMPIReturn(MPI_Waitall( count, rcvid, MPI_STATUSES_IGNORE),__FILE__,__LINE__);
    }else{
      if(sendcnt > 0){
	CheckMPIReturn(MPI_Recv( &hs, 1, MPI_INT, root, mtag, comm, &status), __FILE__,__LINE__);
	CheckMPIReturn(MPI_Send( sendbuf, sendcnt, sendtype, root, mtag, comm), __FILE__,__LINE__);
      }
    }
  }else{
    CheckMPIReturn(MPI_Gatherv ( sendbuf, sendcnt, sendtype, recvbuf, recvcnts, displs, recvtype, root, comm), __FILE__,__LINE__);
  }

  return PIO_NOERR;
}

///
///  @brief Returns the smallest power of 2 greater than i
///  
int ceil2(const int i)
{
  int p=1;
  while(p<i){
    p*=2;
  }
  return(p);
}

///
///  @brief Given integers p and k between 0 and np-1,  
///  if (p+1)^k <= np-1 then return (p+1)^k else -1
int pair(const int np, const int p, const int k)
{
  int q = (p+1) ^ k ;
  int pair = (q > np-1)? -1: q;
  return pair;
}

/**
 **  @brief Provides the functionality of MPI_Alltoallw with flow control options
 */
int pio_swapm(void *sndbuf,   int sndlths[], int sdispls[],  MPI_Datatype stypes[], 
	      void *rcvbuf,  int rcvlths[],  int rdispls[],  MPI_Datatype rtypes[], 
	       MPI_Comm comm,const  bool handshake, bool isend,const  int max_requests)
{

  int nprocs;
  int mytask;

  int maxsend=0;
  int maxrecv=0;

  CheckMPIReturn(MPI_Comm_size(comm, &nprocs),__FILE__,__LINE__);
  CheckMPIReturn(MPI_Comm_rank(comm, &mytask),__FILE__,__LINE__);


  if(max_requests == 0) {
#ifdef DEBUG
    int totalrecv=0;
    int totalsend=0;
    for(int i=0;i<nprocs;i++){
      //      printf("%d sndlths %d %d %d %d\n",i,sndlths[i],sdispls[i],rcvlths[i],rdispls[i]);
      totalsend+=sndlths[i];
      totalrecv+=rcvlths[i];
    }
    printf("%s %d totalsend %d totalrecv %d \n",__FILE__,__LINE__,totalsend,totalrecv);

#endif
#ifdef OPEN_MPI
    for(int i=0;i<nprocs;i++){
      if(stypes[i]==MPI_DATATYPE_NULL){
	stypes[i]=MPI_CHAR;
      }
      if(rtypes[i]==MPI_DATATYPE_NULL){
	rtypes[i]=MPI_CHAR;
      }
    }
#endif
    CheckMPIReturn(MPI_Alltoallw( sndbuf, sndlths, sdispls, stypes, rcvbuf, rcvlths, rdispls, rtypes, comm),__FILE__,__LINE__);
    
#ifdef OPEN_MPI
    for(int i=0;i<nprocs;i++){
      if(stypes[i]==MPI_CHAR){
        stypes[i]=MPI_DATATYPE_NULL;
      }
      if(rtypes[i]==MPI_CHAR){
	rtypes[i]=MPI_DATATYPE_NULL;
      }
    }
#endif
    return PIO_NOERR;
  }


  int tag;
  int offset_t;
  int ierr;
  MPI_Status status;
  int steps;
  int istep;
  int rstep;
  int p;
  int maxreq;
  int maxreqh;
  int hs;
  int cnt;
  void *ptr;
  MPI_Request rcvids[nprocs];



  offset_t = nprocs;
  // send to self
  if(sndlths[mytask] > 0){
    void *sptr, *rptr;
    int extent, lb;
    tag = mytask + offset_t;
    sptr = (void *)((char *) sndbuf + sdispls[mytask]);
    rptr = (void *)((char *) rcvbuf  + rdispls[mytask]);

    /*
      MPI_Type_get_extent(stypes[mytask], &lb, &extent);
      printf("%s %d %d %d\n",__FILE__,__LINE__,extent, lb);
      MPI_Type_get_extent(rtypes[mytask], &lb, &extent);
      printf("%s %d %d %d\n",__FILE__,__LINE__,extent, lb);
    */
#ifdef ONEWAY
    CheckMPIReturn(MPI_Sendrecv(sptr, sndlths[mytask],stypes[mytask],
				mytask, tag, rptr, rcvlths[mytask], rtypes[mytask],
				mytask, tag, comm, &status),__FILE__,__LINE__);
#else
    //   printf("%s %d \n",__FILE__,__LINE__);
    CheckMPIReturn(MPI_Irecv(rptr, rcvlths[mytask], rtypes[mytask],
			     mytask, tag, comm, rcvids),__FILE__,__LINE__);
    //printf("%s %d \n",__FILE__,__LINE__);
    CheckMPIReturn(MPI_Send(sptr, sndlths[mytask], stypes[mytask],
			     mytask, tag, comm),__FILE__,__LINE__);

    //printf("%s %d %d\n",__FILE__,__LINE__,rcvids[0]);
    CheckMPIReturn(MPI_Wait(rcvids, &status),__FILE__,__LINE__);

#endif


  }
  if(nprocs==1)
    return PIO_NOERR;

  int swapids[nprocs];
  MPI_Request sndids[nprocs];
  MPI_Request hs_rcvids[nprocs];

  for(int i=0;i<nprocs;i++){
    rcvids[i] = MPI_REQUEST_NULL;
    swapids[i]=0;
  }
  if(isend)
    for(int i=0;i<nprocs;i++)
      sndids[i]=MPI_REQUEST_NULL;
  if(handshake)
    for(int i=0;i<nprocs;i++)
      hs_rcvids[i]=MPI_REQUEST_NULL;
  
  steps = 0;
  for(istep=0;istep<ceil2(nprocs)-1;istep++){
    p = pair(nprocs, istep, mytask) ;
    if( p >= 0 && (sndlths[p] > 0 || rcvlths[p] > 0)){
      swapids[steps++] = p;
    }
  }

  if(steps == 1){
    maxreq = 1;
    maxreqh = 1;
  }else{
    if(max_requests > 1 && max_requests<steps){
      maxreq = max_requests;
      maxreqh = maxreq/2;
    }else if(max_requests>=steps){
      maxreq = steps;
      maxreqh = steps;
    }else{
      maxreq = 2;
      maxreqh = 1;
    }
  } 
  if(handshake){
    hs = 1;
    for(istep=0; istep<maxreq; istep++){
      p = swapids[istep];
      if( sndlths[p] > 0){
	tag = mytask+offset_t;
	CheckMPIReturn(MPI_Irecv( &hs, 1, MPI_INT, p, tag, comm, hs_rcvids+istep), __FILE__,__LINE__);
      }
    }
  }
  for(istep=0;istep < maxreq; istep++){
    p = swapids[istep];
    if(rcvlths[p] > 0){
      tag = p + offset_t;
      ptr = (void *)((char *) rcvbuf + rdispls[p]);

      //	  printf("%s %d %d %d\n",__FILE__,__LINE__,p,(int) rtypes[p]);
      CheckMPIReturn(MPI_Irecv( ptr, rcvlths[p], rtypes[p], p, tag, comm, rcvids+istep), __FILE__,__LINE__);

      if(handshake)
	CheckMPIReturn(MPI_Send( &hs, 1, MPI_INT, p, tag, comm), __FILE__,__LINE__);
    }
  }

  rstep = maxreq;
  for(istep = 0; istep < steps; istep++){
    p = swapids[istep];
    if(sndlths[p] > 0){
      tag = mytask + offset_t;
      if(handshake){
	CheckMPIReturn(MPI_Wait ( hs_rcvids+istep, &status), __FILE__,__LINE__);
	hs_rcvids[istep] = MPI_REQUEST_NULL;
      }
      ptr = (void *)((char *) sndbuf + sdispls[p]);

      if(isend){
	CheckMPIReturn(MPI_Irsend(ptr, sndlths[p], stypes[p], p, tag, comm,sndids+istep), __FILE__,__LINE__);
      }else{
	CheckMPIReturn(MPI_Send(ptr, sndlths[p], stypes[p], p, tag, comm), __FILE__,__LINE__);
      }

    }
    if(istep > maxreqh){
      p = istep - maxreqh;
      if(rcvids[p] != MPI_REQUEST_NULL){
	CheckMPIReturn(MPI_Wait(rcvids+p, &status), __FILE__,__LINE__);
	rcvids[p] = MPI_REQUEST_NULL;
      }
      if(rstep < steps){
	p = swapids[rstep];
	if(handshake && sndlths[p] > 0){
	  tag = mytask + offset_t;
	  CheckMPIReturn(MPI_Irecv( &hs, 1, MPI_INT, p, tag, comm, hs_rcvids+rstep), __FILE__,__LINE__);
	}
	if(rcvlths[p] > 0){
	  tag = p + offset_t;
	  
	  ptr = (void *)((char *) rcvbuf + rdispls[p]);
	  CheckMPIReturn(MPI_Irecv( ptr, rcvlths[p], rtypes[p], p, tag, comm, rcvids+rstep), __FILE__,__LINE__);
	  if(handshake)
	    CheckMPIReturn(MPI_Send( &hs, 1, MPI_INT, p, tag, comm), __FILE__,__LINE__);
	}
	rstep++;
      }
    }
  }
  //     printf("%s %d %d \n",__FILE__,__LINE__,nprocs);
  if(steps>0){
    CheckMPIReturn(MPI_Waitall(steps, rcvids, MPI_STATUSES_IGNORE), __FILE__,__LINE__);
    if(isend)
      CheckMPIReturn(MPI_Waitall(steps, sndids, MPI_STATUSES_IGNORE), __FILE__,__LINE__);
  }
  //      printf("%s %d %d \n",__FILE__,__LINE__,nprocs);
  
  return PIO_NOERR;
}

#ifdef TESTSWAPM_OLD
#include <sys/time.h>
/*
This program tests MPI_Alltoallw by having processor i send different
amounts of data to each processor.
The first test sends i items to processor i from all processors.
*/
int main( int argc, char **argv )
{
    MPI_Comm comm;
    int *sbuf, *rbuf;
    int rank, size;
    int *sendcounts, *recvcounts, *rdispls, *sdispls;
    int i, j, *p, err;
    MPI_Datatype *sendtypes, *recvtypes;
    struct timeval t1, t2;
    int msg_cnt;

    MPI_Init( &argc, &argv );
    err = 0;
    comm = MPI_COMM_WORLD;
    /* Create the buffer */
    MPI_Comm_size( comm, &size );
    MPI_Comm_rank( comm, &rank );
    sbuf = (int *)malloc( size * size * sizeof(int) );
    rbuf = (int *)malloc( size * size * sizeof(int) );
    if (!sbuf || !rbuf) {
        fprintf( stderr, "Could not allocated buffers!\n" );
        fflush(stderr);
        MPI_Abort( comm, 1 );
    }
    /* Test pio_fc_gather */
    
    msg_cnt=0;
    while(msg_cnt<= MAX_GATHER_BLOCK_SIZE){
      /* Load up the buffers */
      for (i=0; i<size*size; i++) {
        sbuf[i] = i + 100*rank;
        rbuf[i] = -i;
      }

      MPI_Barrier(comm);
      if(rank==0) printf("Start gather test %d\n",msg_cnt);
      if(rank == 0) gettimeofday(&t1, NULL);
      
      err = pio_fc_gather( sbuf, size, MPI_INT, rbuf, size, MPI_INT, 0, comm, msg_cnt);


      if(rank == 0){
	gettimeofday(&t2, NULL);
	printf("time = %f\n",t2.tv_sec - t1.tv_sec + 1.e-6*( t2.tv_usec - t1.tv_usec));
      }
      if(rank==0){
	for(j=0;j<size;j++){
	  for(i=0;i<size;i++){
	    if(rbuf[i + j*size] != i + 100*j){ 
	      printf("got %d expected %d\n",rbuf[i+j*size] , i+100*j);
	    }
	  }
	}
      }

      MPI_Barrier(comm);


      if(msg_cnt==0)
	msg_cnt=1;
      else
	msg_cnt*=2;




    }

    

    /* Test pio_swapm Create and load the arguments to alltoallv */
    sendcounts = (int *)malloc( size * sizeof(int) );
    recvcounts = (int *)malloc( size * sizeof(int) );
    rdispls = (int *)malloc( size * sizeof(int) );
    sdispls = (int *)malloc( size * sizeof(int) );
    sendtypes = (MPI_Datatype *)malloc( size * sizeof(MPI_Datatype) );
    recvtypes = (MPI_Datatype *)malloc( size * sizeof(MPI_Datatype) );
    if (!sendcounts || !recvcounts || !rdispls || !sdispls || !sendtypes || !recvtypes) {
        fprintf( stderr, "Could not allocate arg items!\n" );
        fflush(stderr);
        MPI_Abort( comm, 1 );
    }

    for (i=0; i<size; i++) {
        sendcounts[i] = i + 1;
        recvcounts[i] = rank +1;
        rdispls[i] = i * (rank+1) * sizeof(int) ;
        sdispls[i] = (((i+1) * (i))/2) * sizeof(int) ;
        sendtypes[i] = recvtypes[i] = MPI_INT;
    }
    
    //    for(int msg_cnt=4; msg_cnt<size; msg_cnt*=2){
    //   if(rank==0) printf("message count %d\n",msg_cnt);
    msg_cnt = 0;
    for(int itest=0;itest<5; itest++){
      bool hs=false; 
      bool isend=false;
      /* Load up the buffers */
      for (i=0; i<size*size; i++) {
        sbuf[i] = i + 100*rank;
        rbuf[i] = -i;
      }
      MPI_Barrier(comm);

      if(rank==0) printf("Start itest %d\n",itest);
      if(rank == 0) gettimeofday(&t1, NULL);

      if(itest==0){
	err = pio_swapm( size, rank, sbuf,  sendcounts, sdispls, sendtypes, 
			 rbuf,  recvcounts, rdispls, recvtypes, comm, hs, isend, 0);
      }else if(itest==1){
	hs = true;
	isend = true;
	err = pio_swapm( size, rank, sbuf,  sendcounts, sdispls, sendtypes, 
			 rbuf,  recvcounts, rdispls, recvtypes, comm, hs, isend, msg_cnt);
      }else if(itest==2){
	hs = false;
	isend = true;
	err = pio_swapm( size, rank, sbuf, sendcounts, sdispls, sendtypes, 
			 rbuf, recvcounts, rdispls, recvtypes, comm, hs, isend, msg_cnt);

      }else if(itest==3){
	hs = false;
	isend = false;
	err = pio_swapm( size, rank, sbuf, sendcounts, sdispls, sendtypes, 
			 rbuf, recvcounts, rdispls, recvtypes, comm, hs, isend, msg_cnt);

      }else if(itest==4){
	hs = true;
	isend = false;
	err = pio_swapm( size, rank, sbuf,  sendcounts, sdispls, sendtypes, 
			 rbuf,  recvcounts, rdispls, recvtypes, comm, hs, isend, msg_cnt);

      }

      if(rank == 0){
	gettimeofday(&t2, NULL);
	printf("itest = %d time = %f\n",itest,t2.tv_sec - t1.tv_sec + 1.e-6*( t2.tv_usec - t1.tv_usec));
      }
      /*
      printf("scnt: %d %d %d %d\n",sendcounts[0],sendcounts[1],sendcounts[2],sendcounts[3]);
      printf("sdispls: %d %d %d %d\n",sdispls[0],sdispls[1],sdispls[2],sdispls[3]);
      printf("rcnt: %d %d %d %d\n",recvcounts[0],recvcounts[1],recvcounts[2],recvcounts[3]);
      printf("rdispls: %d %d %d %d\n",rdispls[0],rdispls[1],rdispls[2],rdispls[3]);

      printf("send: ");
      for(i=0;i<size*size;i++)
	printf("%d ",sbuf[i]);
      printf("\n");
      printf("recv: ");
      for(i=0;i<size*size;i++)
	printf("%d ",rbuf[i]);
      printf("\n");
      */
      MPI_Barrier(comm);
      /* Check rbuf */
      for (i=0; i<size; i++) {
	p = rbuf + rdispls[i]/sizeof(int);
	for (j=0; j<rank+1; j++) {
	  if (p[j] != i * 100 + (rank*(rank+1))/2 + j) {
	    fprintf( stderr, "[%d] got %d expected %d for %d %dth in itest=%d\n",
		     rank, p[j],(i*100 + (rank*(rank+1))/2+j), i, j, itest);
	    fflush(stderr);
	    err++;
	  }
	}
      }
    }

    //    }

    free( sendtypes );
    free( recvtypes );
    free( sdispls );
    free( rdispls );
    free( recvcounts );
    free( sendcounts );
    free( rbuf );
    free( sbuf );
    MPI_Finalize();
    return 0;
}


#endif

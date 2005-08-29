
#include <stdio.h>
#include "mpi.h"





main(int argc, char *argv[])
{
  MPI_Request sreq[10], sreq2[10], rreq[10], rreq2[10];
  int sbuf[10],sbuf2[10],rbuf[10],rbuf2[10];
  int tag;
  MPI_Status status[10];
  int i;
  MPI_Comm comm2;
  int flag;


  MPI_Initialized(&flag);
  printf("MPI is initialized = %d\n",flag);

  MPI_Init(NULL,NULL);

  MPI_Comm_dup(MPI_COMM_WORLD,&comm2);

  MPI_Initialized(&flag);
  printf("MPI is initialized = %d\n",flag);

  for (i=0; i<5; i++)
    {
      tag=100+i;
      printf("Post receive tag %d\n",tag);

      MPI_Irecv(&rbuf[2*i],1,MPI_2INT,
		0,tag,MPI_COMM_WORLD,&rreq[i]);
    }


  for (i=0; i<5; i++)
    {
      sbuf2[i]=1000+10*i;
      tag=100+i;
      printf("Send %d tag %d\n",sbuf2[i],tag);
      MPI_Isend(&sbuf2[i],1,MPI_INT,0,tag,comm2,&sreq2[i]);
    }

  for (i=0; i<5; i++)
    {
      sbuf[2*i]=10*i;
      sbuf[2*i+1]=10*i+1;
      tag=100+(4-i);
      printf("Send %d tag %d\n",sbuf[i],tag);
      MPI_Isend(&sbuf[2*i],1,MPI_2INT,0,tag,MPI_COMM_WORLD,&sreq[i]);
    }


  MPI_Waitall(5,sreq,status);
  MPI_Waitall(5,rreq,status);

  for (i=0; i<5; i++)
    printf("tag %d rbuf= %d %d\n",status[i].MPI_TAG,rbuf[2*i],rbuf[2*i+1]);


  for (i=0; i<5; i++)
    {
      tag=100+i;
      printf("Post receive tag %d\n",tag);

      MPI_Irecv(&rbuf2[i],1,MPI_INT,
		0,tag,comm2,&rreq2[i]);
    }


  MPI_Waitall(5,sreq2,status);
  MPI_Waitall(5,rreq2,status);

  MPI_Finalize();

}




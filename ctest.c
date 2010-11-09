
#include <stdio.h>
#include "mpi.h"





main(int argc, char *argv[])
{
  MPI_Request sreq[10], sreq2[10], rreq[10], rreq2[10];
  int sbuf[10],sbuf2[10],rbuf[10],rbuf2[10];
  int tag;
  MPI_Status status[10], sr_status;
  int i,j;
  MPI_Comm comm2;
  int flag;
  MPI_Group mygroup;
  char pname[MPI_MAX_PROCESSOR_NAME];
  int pnamelen;

  int position, temp;


  printf("Time: %f\n",MPI_Wtime());

  MPI_Initialized(&flag);
  printf("MPI is initialized = %d\n",flag);

  MPI_Init(NULL,NULL);

  MPI_Get_processor_name(pname,&pnamelen);
  printf("Processor name: %s (len=%d)\n",pname,pnamelen);

#if 0
  MPI_Comm_dup(MPI_COMM_WORLD,&comm2);
#endif

#if 0
  MPI_Comm_split(MPI_COMM_WORLD,42,99,&comm2);
#endif

#if 1
  MPI_Comm_group(MPI_COMM_WORLD,&mygroup);
  MPI_Comm_create(MPI_COMM_WORLD,mygroup,&comm2);
#endif

  MPI_Initialized(&flag);
  printf("MPI is initialized = %d\n",flag);

  for (i=0; i<5; i++)
    {
      tag=100+i;
      printf("COMWORLD Post ireceive tag %d\n",tag);

      MPI_Irecv(&rbuf[2*i],1,MPI_2INT,
		0,tag,MPI_COMM_WORLD,&rreq[i]);
    }

  for (i=0; i<5; i++)
    {
      sbuf2[i]=1000+10*i;
      tag=100+i;
      printf("COM2 Post isend %d tag %d\n",sbuf2[i],tag);
      MPI_Isend(&sbuf2[i],1,MPI_INT,0,tag,comm2,&sreq2[i]);
    }

  for (i=0; i<5; i++)
    {
      sbuf[2*i]=10*i;
      sbuf[2*i+1]=10*i+1;
      tag=100+(4-i);
      printf("COMWORLD Post irsend %d tag %d\n",sbuf[i],tag);
      MPI_Irsend(&sbuf[2*i],1,MPI_2INT,0,tag,MPI_COMM_WORLD,&sreq[i]);
    }


  printf("Time: %f\n",MPI_Wtime());
  MPI_Waitall(5,sreq,status);
  MPI_Waitall(5,rreq,status);

  printf("Waiting for COMWORLD send/receives\n");

  for (i=0; i<5; i++)
    printf("tag %d rbuf= %d %d\n",status[i].MPI_TAG,rbuf[2*i],rbuf[2*i+1]);


  for (i=0; i<5; i++)
    {
      tag=100+i;
      printf("COM2 Post receive tag %d\n",tag);

      MPI_Irecv(&rbuf2[i],1,MPI_INT,
		0,tag,comm2,&rreq2[i]);
    }


  MPI_Waitall(5,sreq2,status);
  MPI_Waitall(5,rreq2,status);

  printf("Waiting for COM2 send/receive\n");

  for (i=0; i<5; i++)
    printf("tag %d rbuf= %d\n",status[i].MPI_TAG,rbuf2[i]);


  /*
   * pack/unpack
   */

  position=0;
  for (i=0; i<5; i++)
    {
      temp=100+i;
      MPI_Pack(&temp, 1, MPI_INT, sbuf, 20, &position, MPI_COMM_WORLD); 
    }

  MPI_Isend( sbuf, position, MPI_PACKED, 0, 0, MPI_COMM_WORLD,&sreq[0]); 

  MPI_Irecv( rbuf, position, MPI_PACKED, 0, 0, MPI_COMM_WORLD, &rreq[0] );
  MPI_Waitall(1,rreq,status);

  printf("Pack/send/unpack: \n");

  position=0;
  for (i=0; i<5; i++)
    {
      MPI_Unpack(rbuf,20,&position,&temp,1,MPI_INT,MPI_COMM_WORLD);
      printf("%d\n",temp);
    }


  /*
   * sendrecv
   */


  sbuf[0]=42;
  rbuf[0]=0;
  sr_status.MPI_SOURCE= -1;
  sr_status.MPI_TAG= -1;

  MPI_Sendrecv(sbuf,1,MPI_INT,0,127,
	       rbuf,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,
	       MPI_COMM_WORLD,&sr_status);

  printf("Done with MPI_Sendrecv, rbuf=%d  source=%d  tag=%d\n",
	 rbuf[0],sr_status.MPI_SOURCE,sr_status.MPI_TAG);


  /*
   * Send to nowhere
   */


  printf("Send to MPI_PROC_NULL\n"); 
  MPI_Send(sbuf, 1, MPI_INT, MPI_PROC_NULL, 77, MPI_COMM_WORLD);
  


  /*
   * Finish up
   */


  MPI_Finalize();


  for (i=0; i<5; i++)
    {
      printf("Time: %f\n",MPI_Wtime());
      sleep(1);
    }
}




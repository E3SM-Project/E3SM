#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "../gptl.h"

static const MPI_Comm comm = MPI_COMM_WORLD;
static int iam;

int main (int argc, char **argv)
{
  int i, ret;
  int commsize;
  int val;
  const int count = 1024;
  const int tag = 98;
  int sendbuf[count];
  int recvbuf[count];
  int *gsbuf;
  int *atoabufsend, *atoabufrecv;
  int sum;
  MPI_Status status;
  MPI_Request sendreq, recvreq;
  int dest;
  int source;
  int resultlen;                      /* returned length of string from MPI routine */
  int provided;                       /* level of threading support in this MPI lib */
  char string[MPI_MAX_ERROR_STRING];  /* character string returned from MPI routine */
  const char *mpiroutine[] = {"MPI_Ssend", "MPI_Send", "MPI_Recv", "MPI_Sendrecv", "MPI_Irecv",
			      "MPI_Isend", "MPI_Waitall", "MPI_Barrier", "MPI_Bcast", "MPI_Allreduce",
			      "MPI_Gather", "MPI_Scatter", "MPI_Alltoall", "MPI_Reduce", "MPI_Issend"};
  const int nroutines = sizeof (mpiroutine) / sizeof (char *);
  double wallclock;

  void chkbuf (const char *, int *, const int, const int);

  /*
  int DebugWait = 1;
  while (DebugWait) {
  }
  */

  ret = GPTLsetoption (GPTLoverhead, 0);       /* Don't print overhead stats */
  ret = GPTLsetoption (GPTLpercent, 0);        /* Don't print percentage stats */
  ret = GPTLsetoption (GPTLabort_on_error, 1); /* Abort on any GPTL error */

  /* 
  ** Only initialize GPTL if ENABLE_PMPI is false.
  ** If it is true, the library will be initialized in MPI_Init()
  */
#ifndef ENABLE_PMPI
  ret = GPTLinitialize ();                     /* Initialize GPTL */
  ret = GPTLstart ("total");                   /* Time the whole program */
#endif

  /* Initialize MPI by using MPI_Init_thread: report back level of MPI support */
  if ((ret = MPI_Init_thread (&argc, &argv, MPI_THREAD_SINGLE, &provided)) != 0) {
    MPI_Error_string (ret, string, &resultlen);
    printf ("%s: error from MPI_Init_thread: %s\n", argv[0], string);
    MPI_Abort (comm, -1);
  }
  
  ret = MPI_Comm_rank (comm, &iam);            /* Get my rank */
  ret = MPI_Comm_size (comm, &commsize);       /* Get communicator size */

  if (iam == 0) {
    printf ("%s: testing suite of MPI routines for auto-instrumentation via GPTL PMPI layer\n", argv[0]);
    switch (provided) {
    case MPI_THREAD_SINGLE:
      printf ("MPI support level is MPI_THREAD_SINGLE\n");
      break;
    case MPI_THREAD_SERIALIZED:
      printf ("MPI support level is MPI_THREAD_SERIALIZED\n");
      break;
    case MPI_THREAD_MULTIPLE:
      printf ("MPI support level is MPI_THREAD_MULTIPLE\n");
      break;
    default:
      printf ("MPI support level is not known\n");
      MPI_Abort (comm, -1);
    }
  }
       
  for (i = 0; i < count; ++i)
    sendbuf[i] = iam;

  dest = (iam + 1)%commsize;
  source = iam - 1;
  if (source < 0)
    source = commsize - 1;

  /* Send, Ssend */
  if (commsize % 2 == 0) {
    if (iam % 2 == 0) {
      ret = MPI_Send (sendbuf, count, MPI_INT, dest, tag, comm);
      ret = MPI_Recv (recvbuf, count, MPI_INT, source, tag, comm, &status);
    } else {
      ret = MPI_Recv (recvbuf, count, MPI_INT, source, tag, comm, &status);
      ret = MPI_Send (sendbuf, count, MPI_INT, dest, tag, comm);
    }
    chkbuf ("MPI_Send + MPI_Recv", recvbuf, count, source);

    if (iam % 2 == 0) {
      ret = MPI_Ssend (sendbuf, count, MPI_INT, dest, tag, comm);
      ret = MPI_Recv (recvbuf, count, MPI_INT, source, tag, comm, &status);
    } else {
      ret = MPI_Recv (recvbuf, count, MPI_INT, source, tag, comm, &status);
      ret = MPI_Ssend (sendbuf, count, MPI_INT, dest, tag, comm);
    }
    chkbuf ("MPI_Ssend + MPI_Recv", recvbuf, count, source);

    if (iam % 2 == 0) {
      ret = MPI_Issend (sendbuf, count, MPI_INT, dest, tag, comm, &sendreq);
      ret = MPI_Recv (recvbuf, count, MPI_INT, source, tag, comm, &status);
    } else {
      ret = MPI_Recv (recvbuf, count, MPI_INT, source, tag, comm, &status);
      ret = MPI_Issend (sendbuf, count, MPI_INT, dest, tag, comm, &sendreq);
    }
    chkbuf ("MPI_Issend + MPI_Recv", recvbuf, count, source);
  }

  ret = MPI_Sendrecv (sendbuf, count, MPI_INT, dest, tag, 
		      recvbuf, count, MPI_INT, source, tag, 
		      comm, &status);
  chkbuf ("MPI_Sendrecv", recvbuf, count, source);

  ret = MPI_Irecv (recvbuf, count, MPI_INT, source, tag, 
		   comm, &recvreq);
  ret = MPI_Isend (sendbuf, count, MPI_INT, dest, tag, 
		   comm, &sendreq);
  ret = MPI_Wait (&recvreq, &status);
  ret = MPI_Wait (&sendreq, &status);
  chkbuf ("MPI_Isend + MPI_Irecv + MPI_Wait", recvbuf, count, source);

  ret = MPI_Irecv (recvbuf, count, MPI_INT, source, tag, 
		   comm, &recvreq);
  ret = MPI_Isend (sendbuf, count, MPI_INT, dest, tag, 
		   comm, &sendreq);
  ret = MPI_Waitall (1, &recvreq, &status);
  ret = MPI_Waitall (1, &sendreq, &status);
  chkbuf ("MPI_Waitall", recvbuf, count, source);

  ret = MPI_Barrier (comm);

  ret = MPI_Bcast (sendbuf, count, MPI_INT, 0, comm);
  chkbuf ("MPI_Bcast", sendbuf, count, 0);

  for (i = 0; i < count; ++i)
    sendbuf[i] = iam;

  ret = MPI_Allreduce (sendbuf, recvbuf, count, MPI_INT, MPI_SUM, comm);
  sum = 0.;
  for (i = 0; i < commsize; ++i) 
    sum += i;
  chkbuf ("MPI_Allreduce", recvbuf, count, sum);

  gsbuf = (int *) malloc (commsize * count * sizeof (int));
  ret = MPI_Gather (sendbuf, count, MPI_INT,
		    gsbuf, count, MPI_INT,
		    0, comm);
  if (iam == 0) {
    val = 0;
    for (i = 0; i < commsize*count; ++i) {
      if (gsbuf[i] != val) {
	printf ("iam=%d MPI_Gather: bad gsbuf[%d]=%d != %d\n", iam, i, gsbuf[i], val);
	MPI_Abort (comm, -1);
      }
      if ((i+1) % count == 0)
	++val;
    }
  }

  ret = MPI_Scatter (gsbuf, count, MPI_INT,
		     recvbuf, count, MPI_INT,
		     0, comm);
  chkbuf ("MPI_Scatter", recvbuf, count, iam);

  atoabufsend = (int *) malloc (commsize * sizeof (int));
  atoabufrecv = (int *) malloc (commsize * sizeof (int));
  for (i = 0; i < commsize; ++i)
    atoabufsend[i] = i;

  ret = MPI_Alltoall (atoabufsend, 1, MPI_INT,
		      atoabufrecv, 1, MPI_INT,
		      comm);

  for (i = 0; i < commsize; ++i)
    if (atoabufrecv[i] != iam) {
      printf ("iam=%d MPI_Alltoall: bad atoabufrecv[%d]=%d != %d\n", 
	      iam, i, atoabufrecv[i], i);
      MPI_Abort (comm, -1);
    }

  ret = MPI_Reduce (sendbuf, recvbuf, count, MPI_INT, MPI_SUM, 0, comm);
  if (iam == 0) {
    sum = 0.;
    for (i = 0; i < commsize; ++i) 
      sum += i;
    chkbuf ("MPI_Reduce", recvbuf, count, sum);
  }

  ret = MPI_Finalize ();          /* Clean up MPI */

#ifndef ENABLE_PMPI
  ret = GPTLstop ("total");
  ret = GPTLpr (iam);             /* Print the results */
#endif

  /* Check that PMPI entries were generated for all expected routines */
  if (iam == 0) {
    for (i = 0; i < nroutines; ++i) {
      printf ("%s: checking that there is a GPTL entry for MPI routine %s...\n", argv[0], mpiroutine[i]);
      ret = GPTLget_wallclock (mpiroutine[i], 0, &wallclock);
      if (ret < 0) {
	printf ("Failure\n");
	return -1;
      }
      printf("Success\n");
    }
  }
  return 0;
}

void chkbuf (const char *msg, int *recvbuf, const int count, const int source)
{
  int i;
  for (i = 0; i < count; ++i)
    if (recvbuf[i] != source) {
      printf ("iam=%d %s:bad recvbuf[%d]=%d != %d\n", iam, msg, i, recvbuf[i], source);
      MPI_Abort (comm, -1);
    }
}

#include <stdio.h>
#include <unistd.h>
#ifdef HAVE_MPI
#include <mpi.h>
#endif
#ifdef THREADED_OMP
#include <omp.h>
#endif
#include "../gptl.h"

int main (int argc, char **argv)
{
  int nthreads = 1;  /* Value is 1 if no threading */
  int iam = 0;       /* Value is 0 if no MPI */
  int commsize = 1;  /* Value is 1 if no MPI */
  int provided = -1; /* level of threading support in this MPI lib */
  int n;
  int ret;

#ifdef HAVE_MPI
  int resultlen;                      /* returned length of string from MPI routine */
  char string[MPI_MAX_ERROR_STRING];  /* character string returned from MPI routine */

  /* Initialize MPI by using MPI_Init_thread: report back level of MPI support */
  if ((ret = MPI_Init_thread (&argc, &argv, MPI_THREAD_SINGLE, &provided)) != 0) {
    MPI_Error_string (ret, string, &resultlen);
    printf ("%s: error from MPI_Init_thread: %s\n", argv[0], string);
    MPI_Abort (MPI_COMM_WORLD, -1);
  }
  
  ret = MPI_Comm_rank (MPI_COMM_WORLD, &iam);            /* Get my rank */
  ret = MPI_Comm_size (MPI_COMM_WORLD, &commsize);       /* Get communicator size */
#endif

  if (iam == 0) {
    printf ("%s: testing GPTLpr() and GPTLpr_summary() with some timers ON\n", argv[0]);
    printf ("Check timing.* files: 1st and last ranks, 1st and last threads should print error\n");
#ifdef HAVE_MPI
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
      MPI_Abort (MPI_COMM_WORLD, -1);
    }
#endif
  }

  ret = GPTLsetoption (GPTLoverhead, 0);       /* Don't print overhead stats */
  ret = GPTLsetoption (GPTLpercent, 0);        /* Don't print percentage stats */
  ret = GPTLinitialize ();                     /* Initialize GPTL */

  ret = GPTLstart ("total");
  /* Everyone starts "sub", but 1st and last ranks erroneously start it twice */
  ret = GPTLstart ("sub");
  if (iam == 0 || iam == commsize-1)
    ret = GPTLstart ("sub");

#ifdef THREADED_OMP
  nthreads = omp_get_max_threads ();
#endif

  if (iam == 0)
    printf ("nthreads=%d ntasks=%d\n", nthreads, commsize);

#pragma omp parallel for private (ret)
  for (n = 0; n < nthreads; ++n) {
    ret = GPTLstart ("threaded_region");
    ret = GPTLstart ("threaded_region_sub");

    /* sleep a short time so timings are meaningful */
    ret = sleep (iam+n);

    /* Everyone starts "threaded_region_sub", but 1st and last threads erroneously start it twice */
    if (n == 0 || n == nthreads-1)
      ret = GPTLstart ("threaded_region_sub");

    ret = GPTLstop ("threaded_region_sub");
    ret = GPTLstop ("threaded_region");
  }

  ret = GPTLstop ("sub");
  ret = GPTLstop ("total");
  ret = GPTLpr (iam);
#ifdef HAVE_MPI
  ret = GPTLpr_summary (MPI_COMM_WORLD);
  ret = MPI_Finalize ();
#else
  ret = GPTLpr_summary ();
#endif
  return 0;
}
